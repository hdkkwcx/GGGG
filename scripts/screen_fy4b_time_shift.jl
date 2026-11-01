#!/usr/bin/env julia

# Screen a whole-hour time shift for FY4B: does the gauge at hour t match FY4B at hour t+k better
# than FY4B at t? It reads the rebuilt FY4B series (w3, before any shift or quantile map), the
# gauges and GPM_uncal on the qmfusion run's grid, and fits nothing but the choice of k.
#
#   julia --project=. scripts/screen_fy4b_time_shift.jl
#
# Hours are hour-ending Beijing labels on both sides (`FY4BPreprocessing.aligned_hour_end`), so a
# lag found here is in the data, not in the labels. k > 0 means FY4B's later hour; k runs over
# `LAGS`. The strata are season × regime, the regime read from fy4b_hour_meta.csv: R1 is 133°E
# software V1.0, R2 133°E V1.0.0, R3 105°E V1.0.0.
#
# Writes output/fy4b_time_shift_screen/:
# - lag_scores.csv: per product (FY4B, and GPM_uncal as the control), stratum and k, the scores at
#   k and at k = 0 on the cells finite at both, so each k is compared on its own shared cells.
# - diurnal_phase.csv: per stratum, the first-harmonic phase of the mean hourly amount for the
#   gauges, FY4B and GPM_uncal on the cells all three have, and each product's lag behind the gauges.
# - picked_lags.csv: for each balanced_spatial fold and stratum, the k that most raises r on the
#   training stations. A stratum with fewer than `MIN_HOURS` hours keeps k = 0.
# - heldout.csv: the held-out stations at their fold's picked k against k = 0, pooled over the
#   folds, with a paired day-block interval (`DayBlockBootstrap`) on RMSE, r, CSI_0.1 and CSI_2,
#   for every stratum with at least `MIN_HOURS` hours. MAM and DJF get one only pooled over the
#   regimes (`with_interval`): within one regime they hold so few rain days that some resamples
#   have none, and a score with no events has no value.
# - provenance.csv.
#
# Rule, fixed before the run: a per-fold shift is worth building into FY4B_QM only if the pooled
# held-out r rises with an interval above 0, CSI at 0.1 and 2 mm/h is not worse beyond its interval,
# and GPM_uncal's best k over all hours is 0.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("RunProvenance", "TableIO", "DayBlockBootstrap", "SatelliteCalibration",
    "SatelliteTemporalEvaluation")

using CSV, Printf, Statistics
using Main.DayBlockBootstrap: group_sums, sum_stats, day_groups, compare_days, threshold_label
using Main.SatelliteCalibration: product_regimes
using Main.SatelliteTemporalEvaluation: SEASONS, season_of, harmonic_phase,
    circular_hour_difference
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const REBUILD_DIR = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024")
const FY4B_FILE = joinpath(REBUILD_DIR, "fy4b_w3_h0.csv")
const HOUR_META = joinpath(REBUILD_DIR, "fy4b_hour_meta.csv")
const OUTDIR = joinpath(ROOT, "output", "fy4b_time_shift_screen")
const THRESHOLDS = [0.1, 2.0]
const LAGS = sort(collect(-3:3); by=k -> (abs(k), k))
const MIN_HOURS = 500
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261101
const CI_STATS = (:RMSE, :r, Symbol("CSI_0.1"), :CSI_2)
const PRODUCTS = ["FY4B", "GPM_uncal"]
"""`product_regimes` labels, renamed in time order."""
const REGIMES = Dict("133.0E_V1.0" => "R1", "133.0E_V1.0.0" => "R2", "105.0E_V1.0.0" => "R3")
const REGIME_NAMES = ["R1", "R2", "R3"]

"""
The qmfusion run's configuration (its grid, files and folds), and a copy of its `MGERConfig` that
writes to `OUTDIR`: loading the grid writes `global_common_time_qc.csv`, which must not land in
the run.
"""
function run_config()
    cfg = benchmark_config(:full; nested_covariates=true, satellite_wet_blend=true,
        blend_axes=[:constant, :agreement_envelope], fused_anchor=true, gpm_uncal=true,
        fy4b=true, fy4b_file=FY4B_FILE, fy4b_qm=true, fy4b_qm_adaptive=true, fy4b_qm_fusion=true)
    fields = fieldnames(MGERConfig)
    mger = MGERConfig(; merge(NamedTuple{fields}(Tuple(getfield(cfg.mger, f) for f in fields)),
        (; outdir=OUTDIR))...)
    return cfg, mger
end

"""Each cell's value at hour t+k, or NaN when hour t+k is not on the grid."""
function shift_hours(matrix::AbstractMatrix{Float64}, times::AbstractVector{DateTime}, k::Int)
    column = Dict(time => j for (j, time) in enumerate(times))
    out = fill(NaN, size(matrix))
    for (j, time) in enumerate(times)
        other = get(column, time + Hour(k), 0)
        other == 0 || (out[:, j] .= view(matrix, :, other))
    end
    return out
end

"""
Each grid hour's index into `REGIME_NAMES`, read as the benchmark reads it
(`SatelliteCalibration.product_regimes`, which errors on an hour with no regime).
"""
function hour_regimes(times)
    labels, regime_of_time = product_regimes(HOUR_META, times)
    names = [REGIMES[label] for label in labels]
    return [findfirst(==(names[regime]), REGIME_NAMES) for regime in regime_of_time]
end

"""Every stratum as (season, regime, member groups), the pooled ones first."""
function strata()
    group(s, r) = (s - 1) * length(REGIME_NAMES) + r
    out = Tuple{String,String,Vector{Int}}[]
    for (season, seasons) in vcat([("all", eachindex(SEASONS))],
            [(SEASONS[s], [s]) for s in eachindex(SEASONS)])
        for (regime, regimes) in vcat([("all", eachindex(REGIME_NAMES))],
                [(REGIME_NAMES[r], [r]) for r in eachindex(REGIME_NAMES)])
            push!(out, (season, regime, [group(s, r) for s in seasons for r in regimes]))
        end
    end
    return out
end

_stat(stats, name) = getproperty(stats, Symbol(name))

"""Whether a stratum gets a bootstrap interval; see the header."""
with_interval(season, regime) = regime == "all" || season in ("all", "JJA", "SON")

function lag_scores(y, sat, times, column_group, n_group, hours_in)
    groups = repeat(permutedims(column_group), size(y, 1))
    rows = NamedTuple[]
    for product in PRODUCTS
        base = sat[product]
        for k in sort(LAGS)
            shifted = shift_hours(base, times, k)
            mask = .!isnan.(y) .& .!isnan.(base) .& .!isnan.(shifted)
            sums_k = group_sums(shifted, y, mask, groups, n_group, THRESHOLDS)
            sums_0 = group_sums(base, y, mask, groups, n_group, THRESHOLDS)
            for (season, regime, members) in strata()
                total_k = vec(sum(sums_k[members, :]; dims=1))
                total_k[1] == 0 && continue
                at_k = sum_stats(total_k, THRESHOLDS)
                at_0 = sum_stats(vec(sum(sums_0[members, :]; dims=1)), THRESHOLDS)
                fields = Pair{Symbol,Any}[:product => product, :season => season,
                    :regime => regime, :k => k, :hours => sum(hours_in[members]), :n => at_k.n,
                    :r => at_k.r, :r_k0 => at_0.r, :delta_r => at_k.r - at_0.r,
                    :RMSE => at_k.RMSE, :RMSE_k0 => at_0.RMSE, :bias => at_k.bias]
                for threshold in THRESHOLDS, score in ("POD", "FAR", "CSI")
                    name = string(score, "_", threshold_label(threshold))
                    push!(fields, Symbol(name) => _stat(at_k, name))
                    score == "CSI" && push!(fields,
                        Symbol("delta_", name) => _stat(at_k, name) - _stat(at_0, name))
                end
                push!(rows, (; fields...))
            end
        end
    end
    return DataFrame(rows)
end

"""The lag with the largest gain in r per product and stratum (ties go to the smaller |k|)."""
function best_lags(scores)
    rows = NamedTuple[]
    for group in groupby(scores, [:product, :season, :regime]; sort=false)
        ordered = sort(group, [order(:delta_r; rev=true), order(:k; by=abs)])
        push!(rows, (; product=ordered.product[1], season=ordered.season[1],
            regime=ordered.regime[1], hours=ordered.hours[1], best_k=ordered.k[1],
            delta_r=ordered.delta_r[1]))
    end
    return DataFrame(rows)
end

function diurnal_phase(y, sat, times, column_group)
    hour_of_day = hour.(times)
    sources = ["gauge" => y, "FY4B" => sat["FY4B"], "GPM_uncal" => sat["GPM_uncal"]]
    common = .!isnan.(y) .& .!isnan.(sat["FY4B"]) .& .!isnan.(sat["GPM_uncal"])
    rows = NamedTuple[]
    for (season, regime, members) in strata()
        columns = findall(in(members), column_group)
        phases = Dict{String,Any}()
        for (name, values) in sources
            total, count_ = zeros(24), zeros(Int, 24)
            for j in columns, i in axes(values, 1)
                common[i, j] || continue
                total[hour_of_day[j] + 1] += values[i, j]
                count_[hour_of_day[j] + 1] += 1
            end
            present = findall(>(0), count_)
            length(present) >= 3 || continue
            means = total[present] ./ count_[present]
            fit = harmonic_phase(present .- 1, means)
            phases[name] = (; fit..., peak_hour=present[argmax(means)] - 1, n=sum(count_))
        end
        haskey(phases, "gauge") || continue
        for (name, _) in sources
            haskey(phases, name) || continue
            p = phases[name]
            push!(rows, (; season, regime, source=name, n=p.n, mean=p.mean, amplitude=p.amplitude,
                phase_hour=p.phase_hour, peak_hour=p.peak_hour,
                lag_behind_gauge_h=circular_hour_difference(p.phase_hour,
                    phases["gauge"].phase_hour)))
        end
    end
    return DataFrame(rows)
end

function heldout_test(y, fy4b, times, column_group, n_group, hours_in, station_fold)
    shifted = Dict(k => (k == 0 ? fy4b : shift_hours(fy4b, times, k)) for k in LAGS)
    candidate = copy(fy4b)
    picked = NamedTuple[]
    for fold in sort(unique(station_fold))
        train = station_fold .!= fold
        for g in 1:n_group
            columns = findall(==(g), column_group)
            season = SEASONS[div(g - 1, length(REGIME_NAMES)) + 1]
            regime = REGIME_NAMES[mod(g - 1, length(REGIME_NAMES)) + 1]
            best_k, best_gain = 0, 0.0
            if hours_in[g] >= MIN_HOURS
                yy = y[train, columns]
                base = fy4b[train, columns]
                for k in LAGS
                    k == 0 && continue
                    x = shifted[k][train, columns]
                    mask = .!isnan.(yy) .& .!isnan.(base) .& .!isnan.(x)
                    gain = cor(x[mask], yy[mask]) - cor(base[mask], yy[mask])
                    gain > best_gain && ((best_k, best_gain) = (k, gain))
                end
            end
            candidate[.!train, columns] .= shifted[best_k][.!train, columns]
            push!(picked, (; fold, season, regime, hours=hours_in[g], k=best_k,
                training_gain_r=best_gain))
        end
    end
    mask = .!isnan.(y) .& .!isnan.(fy4b) .& .!isnan.(candidate)
    day_group, n_day = day_groups(times, size(y, 1))
    groups = repeat(permutedims(column_group), size(y, 1))
    rows = NamedTuple[]
    for (season, regime, members) in strata()
        sum(hours_in[members]) >= MIN_HOURS && with_interval(season, regime) || continue
        in_stratum = mask .& in.(groups, Ref(Set(members)))
        any(in_stratum) || continue
        sums(x) = group_sums(x, y, in_stratum, day_group, n_day, THRESHOLDS)
        row = compare_days(sums(candidate), sums(fy4b), THRESHOLDS; ci_stats=CI_STATS,
            reps=BOOTSTRAP_REPS, seed=BOOTSTRAP_SEED)
        push!(rows, merge((; season, regime, hours=sum(hours_in[members])), row))
    end
    return DataFrame(picked), DataFrame(rows)
end

function screen_main()
    mkpath(OUTDIR)
    cfg, mger = run_config()
    run_dir = cfg.mger.outdir
    split_path = joinpath(run_dir, "balanced_spatial", "split_common.csv")
    isfile(split_path) || error("the qmfusion run is missing: $run_dir")
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path,
        collect(values(mger.sat_paths)), HOUR_META, split_path)
    provenance = start_provenance(ROOT, OUTDIR; inputs)

    _, ids, data = load_global_common_product_data(mger)
    times = data["FY4B"].times
    y = Matrix{Float64}(data["FY4B"].Y_obs)
    sat = Dict(product => Matrix{Float64}(data[product].Y_sat) for product in PRODUCTS)
    regime = hour_regimes(times)
    season = [findfirst(==(season_of(Date(time))), SEASONS) for time in times]
    column_group = (season .- 1) .* length(REGIME_NAMES) .+ regime
    n_group = length(SEASONS) * length(REGIME_NAMES)
    hours_in = [count(==(g), column_group) for g in 1:n_group]
    splits = CSV.read(split_path, DataFrame)
    fold_of = Dict(string(row.station_id) => row.fold for row in eachrow(splits)
        if row.scheme == "balanced_spatial")
    Set(keys(fold_of)) == Set(ids) || error("the split's stations differ from the grid's")
    station_fold = [fold_of[id] for id in ids]
    @printf("%d stations, %d hours; hours per regime: %s\n", length(ids), length(times),
        join(["$(REGIME_NAMES[r]) $(count(==(r), regime))" for r in eachindex(REGIME_NAMES)], ", "))

    scores = lag_scores(y, sat, times, column_group, n_group, hours_in)
    best = best_lags(scores)
    phases = diurnal_phase(y, sat, times, column_group)
    picked, heldout = heldout_test(y, sat["FY4B"], times, column_group, n_group,
        hours_in, station_fold)

    write_csv_atomic(joinpath(OUTDIR, "lag_scores.csv"), scores)
    write_csv_atomic(joinpath(OUTDIR, "diurnal_phase.csv"), phases)
    write_csv_atomic(joinpath(OUTDIR, "picked_lags.csv"), picked)
    write_csv_atomic(joinpath(OUTDIR, "heldout.csv"), heldout)
    CSV.write(joinpath(OUTDIR, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println("\nBest k by the gain in r (all stations; strata with at least $MIN_HOURS hours):")
    for row in eachrow(best)
        row.hours >= MIN_HOURS || continue
        @printf("  %-9s %-4s %-4s hours %5d  best k %+d  Δr %+.4f\n", row.product, row.season,
            row.regime, row.hours, row.best_k, row.delta_r)
    end
    println("\nDiurnal phase, lag behind the gauges in hours (first harmonic of the mean amount):")
    for row in eachrow(phases)
        row.source == "gauge" && continue
        @printf("  %-4s %-4s %-9s phase %5.2f h  peak %2d  lag %+.2f h\n", row.season,
            row.regime, row.source, row.phase_hour, row.peak_hour, row.lag_behind_gauge_h)
    end
    println("\nPicked k per fold (non-zero only):")
    for row in eachrow(picked)
        row.k == 0 || @printf("  fold %d %-4s %-3s k %+d  training gain r %+.4f\n", row.fold,
            row.season, row.regime, row.k, row.training_gain_r)
    end
    println("\nHeld-out stations, picked k against k = 0 (95% day-block interval):")
    for row in eachrow(heldout)
        @printf("  %-4s %-4s Δr %+.4f (%+.4f..%+.4f)  RMSE %+.2f%% (%+.2f..%+.2f)  CSI0.1 %+.4f (%+.4f..%+.4f)  CSI2 %+.4f (%+.4f..%+.4f)\n",
            row.season, row.regime, row.delta_r, row.delta_r_lo, row.delta_r_hi,
            100 * row.rel_RMSE, 100 * row.rel_RMSE_lo, 100 * row.rel_RMSE_hi,
            row[Symbol("delta_CSI_0.1")], row[Symbol("delta_CSI_0.1_lo")],
            row[Symbol("delta_CSI_0.1_hi")], row.delta_CSI_2, row.delta_CSI_2_lo,
            row.delta_CSI_2_hi)
    end
    pooled = only(eachrow(heldout[(heldout.season .== "all") .& (heldout.regime .== "all"), :]))
    gpm_best = only(eachrow(best[(best.product .== "GPM_uncal") .& (best.season .== "all") .&
        (best.regime .== "all"), :])).best_k
    passes = pooled.delta_r_lo > 0 && pooled[Symbol("delta_CSI_0.1_hi")] >= 0 &&
        pooled.delta_CSI_2_hi >= 0 && gpm_best == 0
    println("\nRule (fixed before the run): GPM_uncal's best k over all hours is $gpm_best; ",
        passes ? "a per-fold shift clears it" : "a per-fold shift does not clear it")
    println("Wrote $(relpath(OUTDIR, ROOT))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    screen_main()
end

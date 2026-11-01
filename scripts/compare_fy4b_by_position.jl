#!/usr/bin/env julia

# Is FY4B_QM better from 105°E? Score the qmfusion run's FY4B fields per regime, against
# GPM_uncal's same fields as the yardstick for each year's weather. Fits nothing: it reads the
# run's out-of-fold files.
#
#   julia --project=. scripts/compare_fy4b_by_position.jl
#
# The run is the 2026-10-31 qmfusion run:
#
#   julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates \
#     --satellite-wet-blend --blend-axis agreement_envelope --fused-anchor --with-fy4b \
#     --fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv --fy4b-qm --fy4b-qm-adaptive \
#     --fy4b-qm-fusion
#
# Regimes, read as the benchmark reads them (`SatelliteCalibration.product_regimes`): R1 is 133°E
# software V1.0 (Jun 2022 - May 2023), R2 133°E V1.0.0 (May 2023 - Feb 2024), R3 105°E V1.0.0
# (Mar - Dec 2024). R2 to R3 changes only the position, R1 to R2 only the software. Each regime
# holds one JJA and one SON, so the comparison between regimes uses those two seasons
# (`PRIMARY`) and the season mix cannot decide it.
#
# Cells: both FY4B_QM's and GPM_uncal's evaluation masks, a gauge value, and every field in
# `FIELDS` finite. `blend_mgwr` is the own-anchor blend, so FY4B_QM's score there is its own and
# owes nothing to the other products' bands.
#
# Writes output/fy4b_by_position/:
# - scores.csv: per scheme, season set and regime, every field's scores (RMSE, bias, r, and POD,
#   FAR, CSI and FBI at 0.1, 2 and 8 mm/h; FAR_0.1 is the false-rain share).
# - gaps.csv: within each of those strata, the `PAIRS` with a paired day-block interval
#   (`DayBlockBootstrap.compare_days`) on RMSE, r, FAR_0.1 and CSI_0.1.
# - differences.csv: on `PRIMARY`, how each pair's gap changes from R2 to R3 (position) and from
#   R1 to R2 (software). The gaps are the RMSE ratio − 1, r difference and FAR_0.1 difference;
#   days are resampled independently within each regime.
# - provenance.csv.
#
# Reading, fixed before the run: "FY4B_QM is better from 105°E" holds only if, on balanced_spatial,
# its gap to GPM_uncal in RMSE and in r shrinks from R2 to R3 with an interval excluding 0, for
# both `raw` and `blend_mgwr`. Each regime is one year, so even then it is "consistent with".
#
# Checks; the script errors if one fails: the run exists, every grid hour has a regime, the
# out-of-fold tables' hours equal the grid's, and FY4B_QM's `raw` scored on its own mask reproduces
# the run's metrics_pooled.csv RMSE.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("RunProvenance", "TableIO", "DayBlockBootstrap", "SatelliteCalibration",
    "SatelliteTemporalEvaluation")

using CSV, Printf, Random, Statistics
using Main.DayBlockBootstrap: read_wide_values, read_wide_mask, group_sums, sum_stats,
    day_groups, compare_days
using Main.SatelliteCalibration: product_regimes
using Main.SatelliteTemporalEvaluation: SEASONS, season_of
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const REBUILD_DIR = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024")
const FY4B_FILE = joinpath(REBUILD_DIR, "fy4b_w3_h0.csv")
const HOUR_META = joinpath(REBUILD_DIR, "fy4b_hour_meta.csv")
const OUTDIR = joinpath(ROOT, "output", "fy4b_by_position")
const THRESHOLDS = [0.1, 2.0, 8.0]
const CI_STATS = (:RMSE, :r, Symbol("FAR_0.1"), Symbol("CSI_0.1"))
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261101
"""`product_regimes` labels, renamed in time order."""
const REGIMES = Dict("133.0E_V1.0" => "R1", "133.0E_V1.0.0" => "R2", "105.0E_V1.0.0" => "R3")
const REGIME_NAMES = ["R1", "R2", "R3"]
const QM = FY4B_QM_PRODUCT
const GPM = "GPM_uncal"
const FIELDS = [("FY4B", "raw"), (QM, "raw"), (QM, "mgwr"), (QM, "blend_mgwr"), (GPM, "raw"),
    (GPM, "mgwr"), (GPM, "blend_mgwr"), (QM, "adw")]
"""(candidate, baseline) pairs: FY4B_QM against GPM_uncal field for field, against adw, and plain FY4B."""
const PAIRS = [((QM, "raw"), (GPM, "raw")), ((QM, "mgwr"), (GPM, "mgwr")),
    ((QM, "blend_mgwr"), (GPM, "blend_mgwr")), ((QM, "blend_mgwr"), (QM, "adw")),
    (("FY4B", "raw"), (GPM, "raw"))]
const PRIMARY = "JJA+SON"
const SEASON_SETS = vcat([("all", SEASONS), (PRIMARY, ["JJA", "SON"])], [(s, [s]) for s in SEASONS])
const REGIME_SETS = vcat([("all", REGIME_NAMES)], [(r, [r]) for r in REGIME_NAMES])
const DIFFERENCES = [("position", "R2", "R3"), ("software", "R1", "R2")]

field_label((product, method)) = "$product:$method"

"""
The qmfusion run's configuration, and a copy of its `MGERConfig` that writes to `OUTDIR`: loading
the grid writes `global_common_time_qc.csv`, which must not land in the run.
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

"""Each grid hour's regime name, read as the benchmark reads it."""
function hour_regimes(times)
    labels, regime_of_time = product_regimes(HOUR_META, times)
    return [REGIMES[labels[regime]] for regime in regime_of_time]
end

"""Candidate's gap to the baseline from two rows of summed `group_sums`."""
function gap(candidate, baseline)
    c, b = sum_stats(candidate, THRESHOLDS), sum_stats(baseline, THRESHOLDS)
    far = Symbol("FAR_0.1")
    return [c.RMSE / b.RMSE - 1, c.r - b.r, getproperty(c, far) - getproperty(b, far)]
end

const GAP_NAMES = ["rmse_ratio", "r_gap", "far_gap"]

"""
How the candidate's gap to the baseline changes from regime `a` to regime `b`, each given as
per-day sums. Days are resampled independently within each regime.
"""
function regime_difference(candidate_a, baseline_a, candidate_b, baseline_b)
    keep_a, keep_b = findall(>(0), candidate_a[:, 1]), findall(>(0), candidate_b[:, 1])
    total(m, rows) = vec(sum(view(m, rows, :); dims=1))
    from = gap(total(candidate_a, keep_a), total(baseline_a, keep_a))
    to = gap(total(candidate_b, keep_b), total(baseline_b, keep_b))
    rng = MersenneTwister(BOOTSTRAP_SEED)
    draws = zeros(BOOTSTRAP_REPS, length(GAP_NAMES))
    for rep in 1:BOOTSTRAP_REPS
        sample_a = rand(rng, keep_a, length(keep_a))
        sample_b = rand(rng, keep_b, length(keep_b))
        draws[rep, :] = gap(total(candidate_b, sample_b), total(baseline_b, sample_b)) .-
            gap(total(candidate_a, sample_a), total(baseline_a, sample_a))
    end
    fields = Pair{Symbol,Float64}[]
    for (j, name) in enumerate(GAP_NAMES)
        append!(fields, [Symbol(name, "_from") => from[j], Symbol(name, "_to") => to[j],
            Symbol("change_", name) => to[j] - from[j],
            Symbol("change_", name, "_lo") => quantile(draws[:, j], 0.025),
            Symbol("change_", name, "_hi") => quantile(draws[:, j], 0.975)])
    end
    return (; fields...)
end

function compare_main()
    mkpath(OUTDIR)
    cfg, mger = run_config()
    run_dir = cfg.mger.outdir
    pooled_path = joinpath(run_dir, "metrics_pooled.csv")
    isfile(pooled_path) || error("the qmfusion run is missing: $run_dir")
    schemes = string.(cfg.cv_schemes)
    path(scheme, product, file) = joinpath(run_dir, scheme, lowercase(product), file)
    mask_path(scheme, product) = path(scheme, product, "common_evaluation_mask.csv")
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path,
        collect(values(mger.sat_paths)), HOUR_META, pooled_path,
        [path(s, p, "oof_$m.csv") for s in schemes for (p, m) in FIELDS],
        [mask_path(s, p) for s in schemes for p in (QM, GPM)])
    provenance = start_provenance(ROOT, OUTDIR; inputs)

    _, ids, data = load_global_common_product_data(mger)
    times = data["FY4B"].times
    y = Matrix{Float64}(data["FY4B"].Y_obs)
    regime = hour_regimes(times)
    season = [season_of(Date(time)) for time in times]
    println("Hours per regime: ", join(["$r $(count(==(r), regime))" for r in REGIME_NAMES], ", "))
    day_group, n_day = day_groups(times, length(ids))
    pooled = CSV.read(pooled_path, DataFrame)

    score_rows, gap_rows, difference_rows = NamedTuple[], NamedTuple[], NamedTuple[]
    for scheme in schemes
        series = Dict(field_label(field) =>
            read_wide_values(path(scheme, field[1], "oof_$(field[2]).csv"), ids, times)
            for field in FIELDS)
        qm_mask = read_wide_mask(mask_path(scheme, QM), ids, times)
        own = qm_mask .& .!isnan.(y) .& .!isnan.(series["$QM:raw"])
        rmse = sqrt(mean((series["$QM:raw"][own] .- y[own]) .^ 2))
        expected = only(pooled[(pooled.scheme .== scheme) .& (pooled.product .== QM) .&
            (pooled.method .== "raw") .& (pooled.group .== "overall"), :RMSE])
        isapprox(rmse, expected; rtol=1e-6) ||
            error("$scheme: FY4B_QM raw scores $rmse here, $expected in metrics_pooled.csv")
        mask = qm_mask .& read_wide_mask(mask_path(scheme, GPM), ids, times) .& .!isnan.(y)
        for field in FIELDS
            mask .&= .!isnan.(series[field_label(field)])
        end
        @printf("%s: FY4B_QM raw reproduces metrics_pooled.csv (RMSE %.6f); %d shared cells\n",
            scheme, rmse, count(mask))
        primary = Dict{String,Dict{String,Matrix{Float64}}}()
        for (season_set, seasons) in SEASON_SETS, (regime_set, regimes) in REGIME_SETS
            columns = in.(season, Ref(seasons)) .& in.(regime, Ref(regimes))
            cells = mask .& permutedims(columns)
            any(cells) || continue
            sums = Dict(field_label(field) =>
                group_sums(series[field_label(field)], y, cells, day_group, n_day, THRESHOLDS)
                for field in FIELDS)
            for field in FIELDS
                stats = sum_stats(vec(sum(sums[field_label(field)]; dims=1)), THRESHOLDS)
                push!(score_rows, merge((; scheme, seasons=season_set, regime=regime_set,
                    field=field_label(field), hours=count(columns)), stats))
            end
            for (candidate, baseline) in PAIRS
                row = compare_days(sums[field_label(candidate)], sums[field_label(baseline)],
                    THRESHOLDS; ci_stats=CI_STATS, reps=BOOTSTRAP_REPS, seed=BOOTSTRAP_SEED)
                push!(gap_rows, merge((; scheme, seasons=season_set, regime=regime_set,
                    candidate=field_label(candidate), baseline=field_label(baseline)), row))
            end
            season_set == PRIMARY && regime_set in REGIME_NAMES && (primary[regime_set] = sums)
        end
        for (comparison, from, to) in DIFFERENCES, (candidate, baseline) in PAIRS
            c, b = field_label(candidate), field_label(baseline)
            push!(difference_rows, merge((; scheme, comparison, from, to, seasons=PRIMARY,
                candidate=c, baseline=b),
                regime_difference(primary[from][c], primary[from][b], primary[to][c],
                    primary[to][b])))
        end
    end

    write_csv_atomic(joinpath(OUTDIR, "scores.csv"), DataFrame(score_rows))
    write_csv_atomic(joinpath(OUTDIR, "gaps.csv"), DataFrame(gap_rows))
    write_csv_atomic(joinpath(OUTDIR, "differences.csv"), DataFrame(difference_rows))
    CSV.write(joinpath(OUTDIR, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println("\nScores on $PRIMARY by regime (RMSE, r, FAR_0.1):")
    for row in score_rows
        row.seasons == PRIMARY || continue
        @printf("  %-16s %-4s %-20s n %8d  RMSE %.4f  r %.3f  FAR0.1 %.3f  bias %+.4f\n",
            row.scheme, row.regime, row.field, row.n, row.RMSE, row.r,
            getproperty(row, Symbol("FAR_0.1")), row.bias)
    end
    println("\nGaps on $PRIMARY by regime (95% day-block interval):")
    for row in gap_rows
        row.seasons == PRIMARY || continue
        @printf("  %-16s %-4s %-20s vs %-20s RMSE %+.2f%% (%+.2f..%+.2f)  Δr %+.3f (%+.3f..%+.3f)  ΔFAR0.1 %+.3f\n",
            row.scheme, row.regime, row.candidate, row.baseline, 100 * row.rel_RMSE,
            100 * row.rel_RMSE_lo, 100 * row.rel_RMSE_hi, row.delta_r, row.delta_r_lo,
            row.delta_r_hi, getproperty(row, Symbol("delta_FAR_0.1")))
    end
    println("\nChange in the gap between regimes on $PRIMARY (95% interval, days resampled per regime):")
    for row in difference_rows
        @printf("  %-16s %-8s %s->%s %-20s vs %-20s RMSE ratio %+.3f->%+.3f change %+.3f (%+.3f..%+.3f)  r gap %+.3f->%+.3f change %+.3f (%+.3f..%+.3f)  FAR gap change %+.3f (%+.3f..%+.3f)\n",
            row.scheme, row.comparison, row.from, row.to, row.candidate, row.baseline,
            row.rmse_ratio_from, row.rmse_ratio_to, row.change_rmse_ratio,
            row.change_rmse_ratio_lo, row.change_rmse_ratio_hi, row.r_gap_from, row.r_gap_to,
            row.change_r_gap, row.change_r_gap_lo, row.change_r_gap_hi, row.change_far_gap,
            row.change_far_gap_lo, row.change_far_gap_hi)
    end
    judged = [row for row in difference_rows if row.scheme == "balanced_spatial" &&
        row.comparison == "position" && row.baseline in ("$GPM:raw", "$GPM:blend_mgwr") &&
        row.candidate in ("$QM:raw", "$QM:blend_mgwr")]
    length(judged) == 2 || error("expected two judged rows, found $(length(judged))")
    better = all(row -> row.change_rmse_ratio_hi < 0 && row.change_r_gap_lo > 0, judged)
    println("\nReading (fixed before the run): ", better ?
        "FY4B_QM's gap to GPM_uncal shrinks from R2 to R3 in both raw and blend_mgwr" :
        "not shown that FY4B_QM is better from 105°E")
    println("Wrote $(relpath(OUTDIR, ROOT))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_main()
end

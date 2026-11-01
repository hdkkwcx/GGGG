#!/usr/bin/env julia

# Screen the two-stage FY4B correction (`SatelliteTwoStage`) on held-out raw FY4B: per fold,
# counting tables of P(gauge wet) and the wet amount per context class, fitted on the training
# gauges. A screen, not a result: no GWR is fitted and nothing downstream (MGWR, the blend) is
# scored, so it runs in minutes.
#
#   julia --project=. scripts/screen_fy4b_two_stage.jl
#
# It reuses the benchmark's own grid, folds and per-fold FY4B_QM: `benchmark_config` with the flags
# of the 2026-10-31 run (`--with-fy4b --fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv
# --fy4b-qm --fy4b-qm-adaptive --fy4b-qm-fusion`, IMERG Uncal). As a guard, the held-out FY4B_QM
# built here must equal that run's `oof_raw.csv` for FY4B_QM, cell by cell, in both schemes, which
# shows the grid, the folds and the calibration are the benchmark's.
#
# Scored on the cells where the gauge and all six series are finite. That is not the benchmark's
# common mask (which also needs every GWR method), so raw numbers differ slightly from its
# metrics_pooled.csv. Writes output/fy4b_two_stage_screen/:
# - metrics.csv: r, RMSE, bias, and POD/FAR/CSI/FBI at 0.1, 2 and 8 mm/h, overall and per
#   regime x season;
# - paired.csv: the change in r and in RMSE of each two-stage anchor against FY4B_QM, FY4B,
#   GPM_uncal and GSMaP, with a day-block bootstrap interval;
# - context_tables.csv: every fold's counts per cell and class;
# - provenance.csv.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("SatelliteTwoStage", "RunProvenance", "TableIO")

using CSV, Printf, Random, Statistics
using Main.SatelliteTwoStage
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const FY4B_FILE = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024", "fy4b_w3_h0.csv")
const OUTDIR = joinpath(ROOT, "output", "fy4b_two_stage_screen")
const THRESHOLDS = [0.1, 2.0, 8.0]
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261101
const SERIES = ["FY4B", "FY4B_QM", "FY4B_TS_E", "FY4B_TS_O", "GPM_uncal", "GSMaP"]
const CANDIDATES = ["FY4B_TS_E", "FY4B_TS_O"]
const BASELINES = ["FY4B_QM", "FY4B", "GPM_uncal", "GSMaP"]

"""
The configuration of the run being reproduced, and a copy of its `MGERConfig` that writes to
`OUTDIR`: loading the grid writes `global_common_time_qc.csv`, which must not land in that run.
"""
function screen_config()
    reference = benchmark_config(:full; nested_covariates=true, satellite_wet_blend=true,
        blend_axes=[:constant, :agreement_envelope], fused_anchor=true, gpm_uncal=true,
        fy4b=true, fy4b_file=FY4B_FILE, fy4b_qm=true, fy4b_qm_adaptive=true,
        fy4b_qm_fusion=true)
    fields = fieldnames(MGERConfig)
    values = NamedTuple{fields}(Tuple(getfield(reference.mger, field) for field in fields))
    return reference, MGERConfig(; merge(values, (; outdir=OUTDIR))...)
end

"""Error unless `held` equals the stored out-of-fold table at `path` in every cell."""
function check_reproduction(path, held, times, ids)
    table = CSV.read(path, DataFrame)
    DateTime.(string.(table.time)) == times || error("$path: its hours differ from the grid's")
    stored = Matrix{Float64}(undef, length(ids), length(times))
    for (i, id) in enumerate(ids)
        stored[i, :] = [ismissing(v) ? NaN : Float64(v) for v in table[!, id]]
    end
    mismatched = count(.!isequal.(stored, held))
    mismatched == 0 || error("held-out FY4B_QM differs from $path in $mismatched cells")
    println("  FY4B_QM reproduces $(relpath(path, ROOT)) in all $(length(held)) cells")
end

_label(threshold) = isinteger(threshold) ? string(Int(threshold)) : string(threshold)

function series_metrics(x, y, mask)
    xs = x[mask]
    ys = y[mask]
    values = Pair{Symbol,Float64}[:r => cor(xs, ys), :RMSE => sqrt(mean(abs2, xs .- ys)),
        :bias => mean(xs .- ys), :mean_obs => mean(ys), :mean_pred => mean(xs)]
    for threshold in THRESHOLDS
        observed = ys .>= threshold
        estimated = xs .>= threshold
        hits = count(observed .& estimated)
        misses = count(observed .& .!estimated)
        false_alarms = count(.!observed .& estimated)
        tag = _label(threshold)
        append!(values, [Symbol("POD_", tag) => hits / (hits + misses),
            Symbol("FAR_", tag) => false_alarms / (hits + false_alarms),
            Symbol("CSI_", tag) => hits / (hits + misses + false_alarms),
            Symbol("FBI_", tag) => (hits + false_alarms) / (hits + misses)])
    end
    return (; n=length(xs), values...)
end

"""Per day: n, Σx, Σy, Σx², Σy², Σxy and Σ(x-y)² over the masked cells."""
function daily_sums(x, y, mask, day_of, n_day)
    sums = zeros(n_day, 7)
    @inbounds for t in axes(y, 2), i in axes(y, 1)
        mask[i, t] || continue
        d = day_of[t]
        xv, yv = x[i, t], y[i, t]
        sums[d, 1] += 1
        sums[d, 2] += xv
        sums[d, 3] += yv
        sums[d, 4] += xv * xv
        sums[d, 5] += yv * yv
        sums[d, 6] += xv * yv
        sums[d, 7] += (xv - yv)^2
    end
    return sums
end

function _r_rmse(total)
    n = total[1]
    mx, my = total[2] / n, total[3] / n
    r = (total[6] / n - mx * my) / sqrt((total[4] / n - mx^2) * (total[5] / n - my^2))
    return r, sqrt(total[7] / n)
end

"""Each candidate against each baseline: the change in r and RMSE, with a paired day-block
bootstrap (whole days resampled for both series at once)."""
function paired_rows(scheme, series, y, mask, times)
    days = Date.(times)
    unique_days = sort(unique(days))
    day_of = [searchsortedfirst(unique_days, day) for day in days]
    sums = Dict(name => daily_sums(series[name], y, mask, day_of, length(unique_days))
        for name in vcat(CANDIDATES, BASELINES))
    keep = findall(>(0), sums[first(CANDIDATES)][:, 1])
    rows = NamedTuple[]
    for candidate in CANDIDATES, baseline in BASELINES
        rng = MersenneTwister(BOOTSTRAP_SEED)
        r1, rmse1 = _r_rmse(vec(sum(sums[candidate]; dims=1)))
        r0, rmse0 = _r_rmse(vec(sum(sums[baseline]; dims=1)))
        delta_r = Float64[]
        delta_rmse = Float64[]
        for _ in 1:BOOTSTRAP_REPS
            sampled = rand(rng, keep, length(keep))
            b1 = _r_rmse(vec(sum(view(sums[candidate], sampled, :); dims=1)))
            b0 = _r_rmse(vec(sum(view(sums[baseline], sampled, :); dims=1)))
            push!(delta_r, b1[1] - b0[1])
            push!(delta_rmse, b1[2] - b0[2])
        end
        push!(rows, (; scheme, candidate, baseline, n=Int(sum(sums[candidate][:, 1])),
            r_baseline=r0, r_candidate=r1, delta_r=r1 - r0,
            delta_r_lo=quantile(delta_r, 0.025), delta_r_hi=quantile(delta_r, 0.975),
            rmse_baseline=rmse0, rmse_candidate=rmse1, delta_rmse=rmse1 - rmse0,
            rel_rmse=(rmse1 - rmse0) / rmse0, delta_rmse_lo=quantile(delta_rmse, 0.025),
            delta_rmse_hi=quantile(delta_rmse, 0.975)))
    end
    return rows
end

function screen_main()
    mkpath(OUTDIR)
    reference, mger = screen_config()
    oof_qm(scheme) = joinpath(reference.mger.outdir, scheme, lowercase(FY4B_QM_PRODUCT),
        "oof_raw.csv")
    grid_path = something(reference.fy4b_qm_grid)
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path,
        [mger.sat_paths[name] for name in sort(collect(keys(mger.sat_paths)))], grid_path,
        joinpath(dirname(grid_path), "fy4b_hour_meta.csv"),
        [oof_qm(string(scheme)) for scheme in reference.cv_schemes])
    provenance = start_provenance(ROOT, OUTDIR; inputs)
    started = time()

    _, ids, product_data = load_global_common_product_data(mger)
    times = product_data["FY4B"].times
    y = product_data["FY4B"].Y_obs
    station_meta = load_station_meta(mger.station_meta_path; station_id_col=mger.station_id_col,
        lon_col=mger.lon_col, lat_col=mger.lat_col)
    lonlat = build_X_lonlat(station_meta, ids)
    println("Grid: $(length(ids)) stations x $(length(times)) hours; reading the shift sums")
    calibration = _prepare_fy4b_calibration!(reference, ids, product_data)
    nearby_cache = Dict{Int,Tuple{Matrix{UInt8},Int}}()
    id_map = Dict(id => index for (index, id) in enumerate(ids))
    groups = vcat(["overall" => trues(length(times))],
        ["$(calibration.regimes[(cell - 1) ÷ 4 + 1]) $(Main.SatelliteCalibration.SEASONS[(cell - 1) % 4 + 1])" =>
            calibration.cells .== cell for cell in sort(unique(calibration.cells))])

    metric_rows = NamedTuple[]
    comparison_rows = NamedTuple[]
    table_rows = NamedTuple[]
    for scheme_symbol in reference.cv_schemes
        scheme = string(scheme_symbol)
        folds = benchmark_folds(scheme_symbol, ids, lonlat; k=reference.k,
            seed=first(benchmark_seeds(reference)), center_init=reference.fold_center_init,
            rotation=0)
        held_out = Dict(name => fill(NaN, size(y)) for name in ("FY4B_QM", CANDIDATES...))
        for fold in 1:reference.k
            _, _, train, val = _fold_station_indices(folds, id_map, reference.k, fold)
            # Exactly the benchmark's call, so the guard below can hold it to that run.
            qm = first(calibrate_fold(
                calibration.grid, calibration.sums, y, train, calibration.regimes,
                calibration.regime_of_time, calibration.cells;
                cache=calibration.cache, wet=mger.rain_threshold,
                adaptive=reference.fy4b_qm_adaptive, choice_of_time=calibration.month_of_time,
                choice_labels=calibration.month_labels,
            ))
            expected, occurrence, rows = two_stage_fold(
                calibration.grid, calibration.sums, y, train, calibration.regimes,
                calibration.regime_of_time, calibration.cells, times, qm;
                cache=calibration.cache, nearby_cache, wet=mger.rain_threshold,
            )
            for (name, anchor) in zip(("FY4B_QM", CANDIDATES...), (qm, expected, occurrence))
                held_out[name][val, :] = anchor[val, :]
            end
            append!(table_rows, [merge((; scheme, fold, n_train_station=length(train)), row)
                for row in rows])
            println("  $scheme fold $fold done ($(round((time() - started) / 60; digits=1)) min)")
        end
        check_reproduction(oof_qm(scheme), held_out["FY4B_QM"], times, ids)
        series = Dict{String,Matrix{Float64}}(held_out)
        for name in ("FY4B", "GPM_uncal", "GSMaP")
            series[name] = product_data[name].Y_sat
        end
        mask = isfinite.(y)
        for name in SERIES
            mask .&= isfinite.(series[name])
        end
        for (group, columns) in groups, name in SERIES
            push!(metric_rows, merge((; scheme, group, series=name),
                series_metrics(series[name], y, mask .& permutedims(columns))))
        end
        append!(comparison_rows, paired_rows(scheme, series, y, mask, times))
    end

    write_csv_atomic(joinpath(OUTDIR, "metrics.csv"), DataFrame(metric_rows))
    write_csv_atomic(joinpath(OUTDIR, "paired.csv"), DataFrame(comparison_rows))
    write_csv_atomic(joinpath(OUTDIR, "context_tables.csv"), DataFrame(table_rows))
    CSV.write(joinpath(OUTDIR, "provenance.csv"), finish_provenance!(provenance, ROOT))

    for row in metric_rows
        row.group == "overall" || continue
        event(name) = join([@sprintf("%.3f", getproperty(row, Symbol(name, "_", _label(threshold))))
            for threshold in THRESHOLDS], "/")
        @printf("%-17s %-10s n=%8d  r %.3f  RMSE %.3f  bias %+.3f  CSI %s  FBI %s\n",
            row.scheme, row.series, row.n, row.r, row.RMSE, row.bias, event("CSI"), event("FBI"))
    end
    for row in comparison_rows
        @printf("%-17s %-10s vs %-10s r %+.4f (%+.4f..%+.4f)  RMSE %+.2f%% (%+.4f..%+.4f)\n",
            row.scheme, row.candidate, row.baseline, row.delta_r, row.delta_r_lo, row.delta_r_hi,
            100 * row.rel_rmse, row.delta_rmse_lo, row.delta_rmse_hi)
    end
    println("Wrote $(relpath(OUTDIR, ROOT)) in $(round((time() - started) / 60; digits=1)) min")
end

if abspath(PROGRAM_FILE) == @__FILE__
    screen_main()
end

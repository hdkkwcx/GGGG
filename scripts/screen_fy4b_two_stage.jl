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
# Three anchors: FY4B_TS_E (P x amount), FY4B_TS_O (the occurrence rule keeping FY4B_QM's amounts)
# and FY4B_TS_I, the amount-only ablation (P x amount on the intensity class alone, no rain nearby
# and no adjacent hours). TS_I against FY4B_QM is what a conditional mean gains without context;
# TS_E against TS_I is what the context adds.
#
# Scored on the cells where the gauge and all seven series are finite. That is not the benchmark's
# common mask (which also needs every GWR method), so raw numbers differ slightly from its
# metrics_pooled.csv. Writes output/fy4b_two_stage_screen/:
# - metrics.csv: r, RMSE, bias, and POD/FAR/CSI/FBI at 0.1, 2 and 8 mm/h, overall and per
#   regime x season;
# - paired.csv: the change in r and in RMSE for each pair in `PAIRS`, with a day-block bootstrap
#   interval;
# - context_tables.csv: every fold's counts per cell and class, for TS_E and TS_O;
# - provenance.csv.
#
#   julia --project=. scripts/screen_fy4b_two_stage.jl --crossed
#
# `--crossed` holds out weeks as well as gauges, as `check_fy4b_qm_time_holdout.jl` does for FY4B_QM,
# and writes to output/fy4b_two_stage_screen_crossed/. Why: that check found FY4B_QM's raw r gain
# (+0.023) gone once weeks are held out, because the map learns a few false-rain episodes, and the
# tables here are counted on the same weeks too. A diagnosis of the screen, not a new fix.
# - Gauges only: the default screen, unchanged. Its point estimates must equal the default run's
#   paired.csv (only the intervals change, to whole weeks).
# - Crossed (suffix `_x`): the same time folds as that check (`TimeBlockFolds`, 7-day blocks, 5 folds,
#   balanced within regime x season, seed 20261107). Each of the 25 pairs of gauge fold and time fold
#   sets every gauge to NaN outside its training hours (held-out blocks and 1 day either side), calls
#   `calibrate_fold` and `two_stage_fold` unchanged, and keeps its held-out gauges in its held-out
#   blocks. The pixel shift is not held out in time (it comes from the run's own shift sums).
# - Intervals resample whole 7-day blocks, in both versions.
# The reading, fixed before the run (balanced spatial judged, random split reported only):
# 1. MATERIAL if TS_E's gain over FY4B_QM moves, between the two versions, by more than 0.01 in r or
#    2 points in RMSE (the limits of the FY4B_QM check).
# 2. SURVIVES if, crossed, TS_E's change in r against FY4B_QM and against plain FY4B both have the
#    lower end of their 95% interval above 0; otherwise it does not survive.
# What the context adds crossed (TS_E_x against TS_I_x) is reported. Also writes reading.csv.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("SatelliteTwoStage", "RunProvenance", "TableIO", "SatelliteCalibration",
    "SatelliteTemporalEvaluation", "TimeBlockFolds")

using CSV, Printf, Random, Statistics
using Main.SatelliteTwoStage
using Main.SatelliteCalibration: product_regimes
using Main.SatelliteTemporalEvaluation: season_of
using Main.TimeBlockFolds: label_day, time_block_folds, crossed_splits
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const FY4B_FILE = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024", "fy4b_w3_h0.csv")
const OUTDIR = joinpath(ROOT, "output", "fy4b_two_stage_screen")
const REGIMES = Dict("133.0E_V1.0" => "R1", "133.0E_V1.0.0" => "R2", "105.0E_V1.0.0" => "R3")
const TIME_FOLD_SEED = 20261107
const MATERIAL_RMSE_POINTS = 2.0
const MATERIAL_R = 0.01
const JUDGED = "balanced_spatial"
const THRESHOLDS = [0.1, 2.0, 8.0]
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261101
const SERIES = ["FY4B", "FY4B_QM", "FY4B_TS_I", "FY4B_TS_E", "FY4B_TS_O", "GPM_uncal", "GSMaP"]
const ANCHORS = ["FY4B_QM", "FY4B_TS_I", "FY4B_TS_E", "FY4B_TS_O"]
"""(candidate, baseline) pairs compared in `paired.csv`."""
const PAIRS = [
    ("FY4B_TS_I", "FY4B_QM"), ("FY4B_TS_E", "FY4B_TS_I"), ("FY4B_TS_E", "FY4B_QM"),
    ("FY4B_TS_O", "FY4B_QM"), ("FY4B_TS_E", "FY4B"), ("FY4B_TS_E", "GPM_uncal"),
    ("FY4B_TS_E", "GSMaP"),
]
"""`--crossed`: the crossed anchors' names, and the pairs compared, gauges only then crossed."""
crossed_name(name) = name * "_x"
const CROSSED_PAIRS = vcat(
    [(candidate, baseline) for (candidate, baseline) in PAIRS[1:5]],
    [(crossed_name(candidate), baseline == "FY4B" ? baseline : crossed_name(baseline))
        for (candidate, baseline) in PAIRS[1:5]],
    [(crossed_name("FY4B_TS_E"), "FY4B_TS_E")],
)

"""
The configuration of the run being reproduced, and a copy of its `MGERConfig` that writes to
`OUTDIR`: loading the grid writes `global_common_time_qc.csv`, which must not land in that run.
"""
function screen_config(outdir=OUTDIR)
    reference = benchmark_config(:full; nested_covariates=true, satellite_wet_blend=true,
        blend_axes=[:constant, :agreement_envelope], fused_anchor=true, gpm_uncal=true,
        fy4b=true, fy4b_file=FY4B_FILE, fy4b_qm=true, fy4b_qm_adaptive=true,
        fy4b_qm_fusion=true)
    fields = fieldnames(MGERConfig)
    values = NamedTuple{fields}(Tuple(getfield(reference.mger, field) for field in fields))
    return reference, MGERConfig(; merge(values, (; outdir))...)
end

"""Error unless the gauges-only rows' point estimates equal the default screen's `paired.csv`."""
function check_default_points(path, rows)
    stored = CSV.read(path, DataFrame)
    for row in rows
        expected = only(filter(r -> r.scheme == row.scheme && r.candidate == row.candidate &&
            r.baseline == row.baseline, stored))
        ok = expected.n == row.n && isapprox(expected.delta_r, row.delta_r; rtol=1e-9, atol=1e-12) &&
            isapprox(expected.delta_rmse, row.delta_rmse; rtol=1e-9, atol=1e-12)
        ok || error("$(row.scheme) $(row.candidate) vs $(row.baseline) differs from $path")
    end
    println("  gauges-only point estimates reproduce $(relpath(path, ROOT)) ($(length(rows)) pairs)")
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

"""Each pair in `pairs`: the change in r and RMSE, with a paired day-block bootstrap (whole days
resampled for both series at once). `group_of_time` (each hour's block, from 1) resamples those
blocks instead of days."""
function paired_rows(scheme, series, y, mask, times; pairs=PAIRS, group_of_time=nothing)
    if group_of_time === nothing
        days = Date.(times)
        unique_days = sort(unique(days))
        day_of = [searchsortedfirst(unique_days, day) for day in days]
        n_group = length(unique_days)
    else
        day_of = group_of_time
        n_group = maximum(group_of_time)
    end
    sums = Dict(name => daily_sums(series[name], y, mask, day_of, n_group)
        for name in unique(vcat(first.(pairs), last.(pairs))))
    keep = findall(>(0), sums[first(first(pairs))][:, 1])
    rows = NamedTuple[]
    for (candidate, baseline) in pairs
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

function screen_main(args=ARGS)
    crossed = "--crossed" in args
    outdir = crossed ? OUTDIR * "_crossed" : OUTDIR
    default_paired = joinpath(OUTDIR, "paired.csv")
    mkpath(outdir)
    reference, mger = screen_config(outdir)
    oof_qm(scheme) = joinpath(reference.mger.outdir, scheme, lowercase(FY4B_QM_PRODUCT),
        "oof_raw.csv")
    grid_path = something(reference.fy4b_qm_grid)
    hour_meta = joinpath(dirname(grid_path), "fy4b_hour_meta.csv")
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path,
        [mger.sat_paths[name] for name in sort(collect(keys(mger.sat_paths)))], grid_path,
        hour_meta, [oof_qm(string(scheme)) for scheme in reference.cv_schemes],
        crossed ? [default_paired] : String[])
    provenance = start_provenance(ROOT, outdir; inputs)
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
    if crossed
        # The time folds of check_fy4b_qm_time_holdout.jl: the same strata, block length and seed.
        labels, regime_index = product_regimes(hour_meta, times)
        regime = [REGIMES[labels[index]] for index in regime_index]
        blocks = time_block_folds(times, string.(regime, "_", season_of.(label_day.(times)));
            seed=TIME_FOLD_SEED)
        write_csv_atomic(joinpath(outdir, "time_folds.csv"), DataFrame(; time=times,
            block=blocks.block_of_time, time_fold=blocks.fold_of_time))
        println("  $(length(blocks.fold_of_block)) week blocks, hours per time fold: ",
            join([count(==(f), blocks.fold_of_time) for f in 1:reference.k], ", "))
    end
    two_stage(observed, train, qm; context=true) = two_stage_fold(
        calibration.grid, calibration.sums, observed, train, calibration.regimes,
        calibration.regime_of_time, calibration.cells, times, qm;
        cache=calibration.cache, nearby_cache, wet=mger.rain_threshold, context,
    )

    metric_rows = NamedTuple[]
    comparison_rows = NamedTuple[]
    table_rows = NamedTuple[]
    for scheme_symbol in reference.cv_schemes
        scheme = string(scheme_symbol)
        folds = benchmark_folds(scheme_symbol, ids, lonlat; k=reference.k,
            seed=first(benchmark_seeds(reference)), center_init=reference.fold_center_init,
            rotation=0)
        held_out = Dict(name => fill(NaN, size(y)) for name in ANCHORS)
        fold_of_station = zeros(Int, length(ids))
        for fold in 1:reference.k
            _, _, train, val = _fold_station_indices(folds, id_map, reference.k, fold)
            fold_of_station[val] .= fold
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
            amount_only = first(two_stage_fold(
                calibration.grid, calibration.sums, y, train, calibration.regimes,
                calibration.regime_of_time, calibration.cells, times, qm;
                cache=calibration.cache, nearby_cache, wet=mger.rain_threshold, context=false,
            ))
            for (name, anchor) in zip(ANCHORS, (qm, amount_only, expected, occurrence))
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
        scored = SERIES
        if crossed
            for name in ANCHORS
                series[crossed_name(name)] = fill(NaN, size(y))
            end
            for split in crossed_splits(fold_of_station, blocks.fold_of_time, times)
                observed = fill(NaN, size(y))
                observed[:, split.train_hours] = y[:, split.train_hours]
                qm = first(calibrate_fold(
                    calibration.grid, calibration.sums, observed, split.train_stations,
                    calibration.regimes, calibration.regime_of_time, calibration.cells;
                    cache=calibration.cache, wet=mger.rain_threshold,
                    adaptive=reference.fy4b_qm_adaptive, choice_of_time=calibration.month_of_time,
                    choice_labels=calibration.month_labels,
                ))
                expected, occurrence, _ = two_stage(observed, split.train_stations, qm)
                amount_only = first(two_stage(observed, split.train_stations, qm; context=false))
                for (name, anchor) in zip(ANCHORS, (qm, amount_only, expected, occurrence))
                    series[crossed_name(name)][split.test_stations, split.test_hours] =
                        anchor[split.test_stations, split.test_hours]
                end
            end
            for name in ANCHORS
                count(isnan.(series[crossed_name(name)]) .!= isnan.(series[name])) == 0 ||
                    error("$scheme: crossed $name is missing cells the gauges-only one has")
            end
            println("  $scheme: 25 crossed fits done ($(round((time() - started) / 60; digits=1)) min)")
            scored = vcat(SERIES, crossed_name.(ANCHORS))
        end
        mask = isfinite.(y)
        for name in scored
            mask .&= isfinite.(series[name])
        end
        for (group, columns) in groups, name in scored
            push!(metric_rows, merge((; scheme, group, series=name),
                series_metrics(series[name], y, mask .& permutedims(columns))))
        end
        if crossed
            rows = paired_rows(scheme, series, y, mask, times; pairs=CROSSED_PAIRS,
                group_of_time=blocks.block_of_time)
            check_default_points(default_paired, rows[1:5])
            append!(comparison_rows, rows)
        else
            append!(comparison_rows, paired_rows(scheme, series, y, mask, times))
        end
    end

    reading = crossed ? crossed_reading(comparison_rows, string.(reference.cv_schemes)) : nothing
    write_csv_atomic(joinpath(outdir, "metrics.csv"), DataFrame(metric_rows))
    write_csv_atomic(joinpath(outdir, "paired.csv"), DataFrame(comparison_rows))
    write_csv_atomic(joinpath(outdir, "context_tables.csv"), DataFrame(table_rows))
    reading === nothing || write_csv_atomic(joinpath(outdir, "reading.csv"), DataFrame(reading))
    CSV.write(joinpath(outdir, "provenance.csv"), finish_provenance!(provenance, ROOT))

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
    if reading !== nothing
        println("\nReading:")
        for row in reading
            @printf("  %-17s %-44s gauges only %+.4f  crossed %+.4f (%+.4f..%+.4f)  %s%s\n", row.scheme,
                row.check, row.gauges_only, row.crossed, row.crossed_lo, row.crossed_hi, row.verdict,
                row.judged ? "" : "  (reported only)")
        end
    end
    println("Wrote $(relpath(outdir, ROOT)) in $(round((time() - started) / 60; digits=1)) min")
end

"""
The reading fixed in the header, per scheme: whether TS_E's gain over FY4B_QM moves materially
when weeks are held out, and whether it survives against FY4B_QM and plain FY4B; the context's
share is reported with no verdict.
"""
function crossed_reading(rows, schemes)
    pick(scheme, c, b) = only(filter(r -> r.scheme == scheme && r.candidate == c && r.baseline == b, rows))
    x = crossed_name
    reading = NamedTuple[]
    for scheme in schemes
        judged = scheme == JUDGED
        g, c = pick(scheme, "FY4B_TS_E", "FY4B_QM"), pick(scheme, x("FY4B_TS_E"), x("FY4B_QM"))
        plain = pick(scheme, x("FY4B_TS_E"), "FY4B")
        context = pick(scheme, x("FY4B_TS_E"), x("FY4B_TS_I"))
        r_move = c.delta_r - g.delta_r
        rmse_move = 100 * (c.rel_rmse - g.rel_rmse)
        append!(reading, [
            (; scheme, judged, check="TS_E vs FY4B_QM: change in r", gauges_only=g.delta_r,
                crossed=c.delta_r, crossed_lo=c.delta_r_lo, crossed_hi=c.delta_r_hi, move=r_move,
                verdict=abs(r_move) > MATERIAL_R ? "MATERIAL" : "not material"),
            (; scheme, judged, check="TS_E vs FY4B_QM: change in RMSE (%)", gauges_only=100 * g.rel_rmse,
                crossed=100 * c.rel_rmse, crossed_lo=100 * c.delta_rmse_lo / c.rmse_baseline,
                crossed_hi=100 * c.delta_rmse_hi / c.rmse_baseline, move=rmse_move,
                verdict=abs(rmse_move) > MATERIAL_RMSE_POINTS ? "MATERIAL" : "not material"),
            (; scheme, judged, check="crossed TS_E vs plain FY4B: change in r",
                gauges_only=pick(scheme, "FY4B_TS_E", "FY4B").delta_r, crossed=plain.delta_r,
                crossed_lo=plain.delta_r_lo, crossed_hi=plain.delta_r_hi, move=NaN,
                verdict=plain.delta_r_lo > 0 ? "above 0" : "not above 0"),
            (; scheme, judged, check="the gain survives (both r intervals above 0)",
                gauges_only=NaN, crossed=NaN, crossed_lo=NaN, crossed_hi=NaN, move=NaN,
                verdict=c.delta_r_lo > 0 && plain.delta_r_lo > 0 ? "SURVIVES" : "DOES NOT SURVIVE"),
            (; scheme, judged, check="context: TS_E vs TS_I, change in r (reported)",
                gauges_only=pick(scheme, "FY4B_TS_E", "FY4B_TS_I").delta_r, crossed=context.delta_r,
                crossed_lo=context.delta_r_lo, crossed_hi=context.delta_r_hi, move=NaN, verdict="-"),
        ])
    end
    return reading
end

if abspath(PROGRAM_FILE) == @__FILE__
    screen_main()
end

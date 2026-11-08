#!/usr/bin/env julia

# Does FY4B_QM's raw gain survive when time is held out as well as gauges? A check, not a result:
# no paper text or number changes, and nothing downstream (MGWR, the blend) is refitted.
#
#   julia --project=. scripts/check_fy4b_qm_time_holdout.jl
#   julia --project=. scripts/check_fy4b_qm_time_holdout.jl <run dir> --matched-products
#
# A run folder reads the run without `--fy4b-qm-fusion` instead (`..._fy4bw3h0_fy4bqm_fy4bqmadapt_gpmuncal`,
# whose name must be what `benchmark_config` gives), and writes to `fy4b_qm_time_holdout` plus that run's
# options from `_fusedanchor` on. FY4B_QM's raw field does not depend on the fusion, so the guard below
# holds it to that run's `oof_raw.csv` as well. `--matched-products` scores on the cells every product of
# the run is scored on (`BenchmarkDiagnostics.matched_evaluation_mask`), adds `_matched`, and checks
# against the quality comparison written with the same run and flag (`compare_fy4b_quality.jl <run dir>
# --matched-products --block-days 7`), whose n, RMSE and r do not depend on the block length.
#
# Why: FY4B_QM is fitted per fold on training gauges, so a held-out gauge never shapes the anchor it
# is scored against. But the two gauges share the hour, and false rain often covers the network at
# once. The map strength is also chosen per month on that same month's training gauges. So §3.4's raw
# gain could be flattered by weather both sides of the split saw.
#
# Two versions of FY4B_QM, each scored at held-out gauges on the same cells:
# - gauges only: exactly the benchmark's `calibrate_fold`. It must equal the 2026-10-31 run's
#   `oof_raw.csv` in every cell, and against plain FY4B it must reproduce output/fy4b_quality/paired.csv.
# - crossed (`TimeBlockFolds`): 7-day blocks in 5 time folds, balanced within regime x season. Each of
#   the 25 pairs of gauge fold and time fold sets every gauge to NaN in the held-out blocks and 1 day
#   either side, calls `calibrate_fold` unchanged, and keeps its held-out gauges in its held-out blocks.
#
# The pixel shift is not held out in time: it comes from the run's own `shift_sums`. It is one choice
# per regime and stable across folds, and recomputing it means rereading the whole grid per time fold.
# So this measures the map and its monthly strength, which is where the freedom is.
#
# The reading, fixed before the run (balanced spatial judged, random split reported only). The gain
# is MATERIAL if, between the two versions against plain FY4B,
# 1. the change in raw RMSE moves by more than 2 points (paper: -16.9%), or
# 2. the change in raw r moves by more than 0.01 (paper: +0.023).
# Either way the paper is unchanged here; if material, the user decides what §3.4 says.
#
# Intervals resample whole 7-day blocks (`compare_days` on week groups), not days.
#
# Writes output/fy4b_qm_time_holdout/: paired.csv, strata.csv (season and regime), monthly.csv,
# strengths.csv (every fit's monthly choice), time_folds.csv, reading.csv and provenance.csv.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("DayBlockBootstrap", "RunProvenance", "TableIO", "SatelliteCalibration",
    "SatelliteTemporalEvaluation", "TimeBlockFolds", "BenchmarkDiagnostics")

using CSV, Printf, Statistics
using Main.BenchmarkDiagnostics: matched_evaluation_mask, matched_fold_failures
using Main.DayBlockBootstrap: read_wide_values, read_wide_mask, group_sums, sum_stats, compare_days
using Main.SatelliteCalibration: product_regimes
using Main.SatelliteTemporalEvaluation: SEASONS, season_of
using Main.TimeBlockFolds: label_day, time_block_folds, crossed_splits
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const REBUILD_DIR = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024")
const FY4B_FILE = joinpath(REBUILD_DIR, "fy4b_w3_h0.csv")
const HOUR_META = joinpath(REBUILD_DIR, "fy4b_hour_meta.csv")
const QUALITY_FILE = joinpath(ROOT, "output", "fy4b_quality", "paired.csv")
const OUTDIR = joinpath(ROOT, "output", "fy4b_qm_time_holdout")
const REGIMES = Dict("133.0E_V1.0" => "R1", "133.0E_V1.0.0" => "R2", "105.0E_V1.0.0" => "R3")
const TIME_FOLD_SEED = 20261107
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261101
const THRESHOLDS = [0.1, 2.0, 8.0]
const CI_STATS = [:RMSE, :r, :abs_bias, Symbol("CSI_0.1"), :CSI_2, :CSI_8, :POD_8]
const PLAIN, GAUGES, CROSSED = "FY4B", "FY4B_QM_gauges", "FY4B_QM_crossed"
const PAIRS = [(GAUGES, PLAIN), (CROSSED, PLAIN), (CROSSED, GAUGES)]
const MATERIAL_RMSE_POINTS = 2.0
const MATERIAL_R = 0.01

"""
The configuration of the run being reproduced, and a copy of its `MGERConfig` that writes to
`outdir`: loading the grid writes `global_common_time_qc.csv`, which must not land in that run.
`run_dir` picks the run with or without `--fy4b-qm-fusion` (by its `_qmfusion` suffix) and must be
the folder that configuration names; without it, the `--fy4b-qm-fusion` run.
"""
function check_config(outdir=OUTDIR; run_dir=nothing)
    qm_fusion = run_dir === nothing || occursin("_qmfusion", basename(rstrip(run_dir, '/')))
    reference = benchmark_config(:full; nested_covariates=true, satellite_wet_blend=true,
        blend_axes=[:constant, :agreement_envelope], fused_anchor=true, gpm_uncal=true,
        fy4b=true, fy4b_file=FY4B_FILE, fy4b_qm=true, fy4b_qm_adaptive=true,
        fy4b_qm_fusion=qm_fusion)
    run_dir === nothing || rstrip(abspath(run_dir), '/') == rstrip(reference.mger.outdir, '/') ||
        error("$run_dir is not a run this check reads; the matching one is $(reference.mger.outdir)")
    fields = fieldnames(MGERConfig)
    values = NamedTuple{fields}(Tuple(getfield(reference.mger, field) for field in fields))
    return reference, MGERConfig(; merge(values, (; outdir))...)
end

"""The run's options from `_fusedanchor` on, which an output folder for that run carries."""
run_tag(run_dir) = last(split(basename(rstrip(run_dir, '/')), "_fusedanchor"; limit=2))

"""Error unless `held` equals the stored out-of-fold table at `path` in every cell."""
function check_reproduction(path, held, times, ids)
    stored = read_wide_values(path, ids, times)
    mismatched = count(.!isequal.(stored, held))
    mismatched == 0 || error("held-out FY4B_QM differs from $path in $mismatched cells")
    println("  gauges-only FY4B_QM reproduces $(relpath(path, ROOT)) in all $(length(held)) cells")
end

"""Error unless gauges-only FY4B_QM against plain FY4B reproduces `quality_file`'s raw row."""
function check_quality(scheme, row, quality_file=QUALITY_FILE)
    quality = CSV.read(quality_file, DataFrame)
    expected = only(filter(r -> r.scheme == scheme && r.product == "FY4B_QM" && r.baseline == "FY4B" &&
        r.method == "raw", quality))
    ok = row.candidate_n == expected.candidate_n &&
        isapprox(row.rel_RMSE, expected.rel_RMSE; rtol=1e-9) &&
        isapprox(row.candidate_r, expected.candidate_r; rtol=1e-9) &&
        isapprox(row.baseline_r, expected.baseline_r; rtol=1e-9)
    ok || error("$scheme: gauges-only FY4B_QM does not reproduce $quality_file")
    @printf("  %s reproduces fy4b_quality: n %d, RMSE %+.2f%%, r %.3f vs %.3f\n", scheme,
        row.candidate_n, 100 * row.rel_RMSE, row.candidate_r, row.baseline_r)
end

"""One row per month and fit: what the fit chose, with its design and folds."""
strength_rows(scheme, design, gauge_fold, time_fold, rows) =
    [(; scheme, design, gauge_fold, time_fold, row.regime, month=row.cell, row.n_pairs,
        row.threshold_applied, row.strength, row.train_rmse_none, row.train_rmse_full,
        row.train_rmse_chosen) for row in rows]

"""Per month, the share of crossed fits that chose what the same gauge fold chose gauges-only."""
function agreement_by_month(rows)
    choice(row) = (row.threshold_applied, row.strength)
    gauges = Dict((row.scheme, row.gauge_fold, row.regime, row.month) => choice(row)
        for row in rows if row.design == "gauges_only")
    crossed = filter(row -> row.design == "crossed", rows)
    keys_seen = unique((row.scheme, row.regime, row.month) for row in crossed)
    return [begin
        fits = filter(row -> (row.scheme, row.regime, row.month) == key, crossed)
        (; scheme=key[1], regime=key[2], month=key[3], n_fits=length(fits),
            share_same=mean(choice(row) == gauges[(row.scheme, row.gauge_fold, row.regime, row.month)]
                for row in fits))
    end for key in keys_seen]
end

function check_main(args=ARGS)
    matched = "--matched-products" in args
    positional = filter(!startswith("--"), args)
    chosen = isempty(positional) ? nothing : abspath(only(positional))
    tag = (chosen === nothing ? "" : run_tag(chosen)) * (matched ? "_matched" : "")
    outdir = OUTDIR * tag
    quality_file = isempty(tag) ? QUALITY_FILE :
        joinpath(ROOT, "output", "fy4b_quality$(tag)_block7d", "paired.csv")
    mkpath(outdir)
    reference, mger = check_config(outdir; run_dir=chosen)
    run_dir = reference.mger.outdir
    path(scheme, product, file) = joinpath(run_dir, scheme, lowercase(product), file)
    schemes = string.(reference.cv_schemes)
    grid_path = something(reference.fy4b_qm_grid)
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path,
        [mger.sat_paths[name] for name in sort(collect(keys(mger.sat_paths)))], grid_path,
        HOUR_META, quality_file,
        [path(s, p, file) for s in schemes for p in (PLAIN, "FY4B_QM")
            for file in ("oof_raw.csv", "common_evaluation_mask.csv")],
        matched ? [joinpath(run_dir, "run_status.csv")] : String[])
    provenance = start_provenance(ROOT, outdir; inputs)
    # Matched: every product the run scored defines the sample, as in run_claim_reassessment.jl.
    run_products = sort(unique(String.(CSV.read(joinpath(run_dir, "metrics_pooled.csv"), DataFrame;
        select=[:product]).product)))
    status = matched ? CSV.read(joinpath(run_dir, "run_status.csv"), DataFrame) : nothing
    started = time()

    _, ids, product_data = load_global_common_product_data(mger)
    times = product_data[PLAIN].times
    y = Matrix{Float64}(product_data[PLAIN].Y_obs)
    station_meta = load_station_meta(mger.station_meta_path; station_id_col=mger.station_id_col,
        lon_col=mger.lon_col, lat_col=mger.lat_col)
    lonlat = build_X_lonlat(station_meta, ids)
    calibration = _prepare_fy4b_calibration!(reference, ids, product_data)
    id_map = Dict(id => index for (index, id) in enumerate(ids))
    println("Grid: $(length(ids)) stations x $(length(times)) hours")

    labels, regime_index = product_regimes(HOUR_META, times)
    regime = [REGIMES[labels[index]] for index in regime_index]
    day = label_day.(times)
    season = season_of.(day)
    blocks = time_block_folds(times, string.(regime, "_", season); seed=TIME_FOLD_SEED)
    write_csv_atomic(joinpath(outdir, "time_folds.csv"), DataFrame(; time=times,
        block=blocks.block_of_time, time_fold=blocks.fold_of_time))
    n_block = length(blocks.fold_of_block)
    block_group = repeat(permutedims(blocks.block_of_time), length(ids))
    println("  $n_block blocks, hours per time fold: ",
        join([count(==(f), blocks.fold_of_time) for f in 1:reference.k], ", "))
    calibrate(observed, train) = calibrate_fold(
        calibration.grid, calibration.sums, observed, train, calibration.regimes,
        calibration.regime_of_time, calibration.cells;
        cache=calibration.cache, wet=mger.rain_threshold,
        adaptive=reference.fy4b_qm_adaptive, choice_of_time=calibration.month_of_time,
        choice_labels=calibration.month_labels,
    )

    paired, strata_rows, monthly, strengths = (NamedTuple[] for _ in 1:4)
    for scheme_symbol in reference.cv_schemes
        scheme = string(scheme_symbol)
        folds = benchmark_folds(scheme_symbol, ids, lonlat; k=reference.k,
            seed=first(benchmark_seeds(reference)), center_init=reference.fold_center_init,
            rotation=0)
        gauges = fill(NaN, size(y))
        fold_of_station = zeros(Int, length(ids))
        for fold in 1:reference.k
            _, _, train, val = _fold_station_indices(folds, id_map, reference.k, fold)
            # Exactly the benchmark's call, so the guard below can hold it to that run.
            qm, _, _, rows = calibrate(y, train)
            gauges[val, :] = qm[val, :]
            fold_of_station[val] .= fold
            append!(strengths, strength_rows(scheme, "gauges_only", fold, 0, rows))
        end
        check_reproduction(path(scheme, "FY4B_QM", "oof_raw.csv"), gauges, times, ids)

        crossed = fill(NaN, size(y))
        for split in crossed_splits(fold_of_station, blocks.fold_of_time, times)
            observed = fill(NaN, size(y))
            observed[:, split.train_hours] = y[:, split.train_hours]
            qm, _, _, rows = calibrate(observed, split.train_stations)
            crossed[split.test_stations, split.test_hours] = qm[split.test_stations, split.test_hours]
            append!(strengths, strength_rows(scheme, "crossed", split.gauge_fold, split.time_fold, rows))
        end
        count(isnan.(crossed) .!= isnan.(gauges)) == 0 ||
            error("$scheme: crossed FY4B_QM is missing cells the gauges-only one has")
        @printf("  %s: 25 crossed fits done (%.1f min)\n", scheme, (time() - started) / 60)

        plain = read_wide_values(path(scheme, PLAIN, "oof_raw.csv"), ids, times)
        series = Dict(PLAIN => plain, GAUGES => gauges, CROSSED => crossed)
        mask = read_wide_mask(path(scheme, "FY4B_QM", "common_evaluation_mask.csv"), ids, times) .&
            read_wide_mask(path(scheme, PLAIN, "common_evaluation_mask.csv"), ids, times) .&
            isfinite.(y) .& isfinite.(plain) .& isfinite.(gauges) .& isfinite.(crossed)
        if matched
            failures = matched_fold_failures(status; scheme, products=run_products)
            isempty(failures) || error("--matched-products: $scheme folds failed:\n" *
                sprint(show, select(failures, :product, :fold, :method, :status)))
            mask .&= matched_evaluation_mask(joinpath(run_dir, scheme), run_products, ids, times).matched
        end
        sums(x, cells) = group_sums(x, y, cells, block_group, n_block, THRESHOLDS)
        compare(candidate, baseline, cells) = compare_days(sums(series[candidate], cells),
            sums(series[baseline], cells), THRESHOLDS; ci_stats=CI_STATS, reps=BOOTSTRAP_REPS,
            seed=BOOTSTRAP_SEED)
        for (candidate, baseline) in PAIRS
            push!(paired, merge((; scheme, candidate, baseline, n_block), compare(candidate, baseline, mask)))
        end
        check_quality(scheme, only(filter(r -> r.scheme == scheme && r.candidate == GAUGES &&
            r.baseline == PLAIN, paired)), quality_file)
        for (unit, label, hours) in vcat([("season", s, season .== s) for s in SEASONS],
                [("regime", r, regime .== r) for r in ("R1", "R2", "R3")])
            cells = mask .& permutedims(hours)
            for (candidate, baseline) in PAIRS
                push!(strata_rows, merge((; scheme, unit, label, candidate, baseline),
                    compare(candidate, baseline, cells)))
            end
        end
        months = sort(unique(firstdayofmonth.(day)))
        month_of_time = [searchsortedfirst(months, firstdayofmonth(d)) for d in day]
        month_group = repeat(permutedims(month_of_time), length(ids))
        for (name, x) in series
            totals = group_sums(x, y, mask, month_group, length(months), THRESHOLDS)
            for m in eachindex(months)
                totals[m, 1] > 0 || continue
                push!(monthly, merge((; scheme, series=name, month=months[m]), sum_stats(totals[m, :], THRESHOLDS)))
            end
        end
    end

    judged = "balanced_spatial"
    pick(c, b, scheme=judged) = only(filter(r -> r.scheme == scheme && r.candidate == c && r.baseline == b, paired))
    reading = NamedTuple[]
    for scheme in schemes
        g, c = pick(GAUGES, PLAIN, scheme), pick(CROSSED, PLAIN, scheme)
        rmse_move = 100 * (c.rel_RMSE - g.rel_RMSE)
        r_move = c.delta_r - g.delta_r
        append!(reading, [
            (; scheme, judged=scheme == judged, check="change in raw RMSE (%)", gauges_only=100 * g.rel_RMSE,
                crossed=100 * c.rel_RMSE, move=rmse_move, limit=MATERIAL_RMSE_POINTS,
                material=abs(rmse_move) > MATERIAL_RMSE_POINTS),
            (; scheme, judged=scheme == judged, check="change in raw r", gauges_only=g.delta_r,
                crossed=c.delta_r, move=r_move, limit=MATERIAL_R, material=abs(r_move) > MATERIAL_R),
        ])
    end

    write_csv_atomic(joinpath(outdir, "paired.csv"), DataFrame(paired))
    write_csv_atomic(joinpath(outdir, "strata.csv"), DataFrame(strata_rows))
    write_csv_atomic(joinpath(outdir, "monthly.csv"), DataFrame(monthly))
    write_csv_atomic(joinpath(outdir, "strengths.csv"), DataFrame(strengths))
    write_csv_atomic(joinpath(outdir, "reading.csv"), DataFrame(reading))
    CSV.write(joinpath(outdir, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println("\nPaired (whole 7-day blocks, 95%):")
    for row in paired
        @printf("  %-17s %-16s vs %-15s n=%7d RMSE %.4f vs %.4f %+.2f%% (%+.2f..%+.2f)  r %.4f vs %.4f %+.4f (%+.4f..%+.4f)  CSI_0.1 %+.4f  CSI_2 %+.4f  CSI_8 %+.4f  POD_8 %+.4f\n",
            row.scheme, row.candidate, row.baseline, row.candidate_n, row.candidate_RMSE, row.baseline_RMSE,
            100 * row.rel_RMSE, 100 * row.rel_RMSE_lo, 100 * row.rel_RMSE_hi, row.candidate_r, row.baseline_r,
            row.delta_r, row.delta_r_lo, row.delta_r_hi, getproperty(row, Symbol("delta_CSI_0.1")),
            row.delta_CSI_2, row.delta_CSI_8, row.delta_POD_8)
    end
    println("\nBy season and regime ($judged):")
    for row in filter(r -> r.scheme == judged, strata_rows)
        @printf("  %-6s %-4s %-16s vs %-15s RMSE %+.2f%% (%+.2f..%+.2f)  r %+.4f (%+.4f..%+.4f)\n", row.unit,
            row.label, row.candidate, row.baseline, 100 * row.rel_RMSE, 100 * row.rel_RMSE_lo,
            100 * row.rel_RMSE_hi, row.delta_r, row.delta_r_lo, row.delta_r_hi)
    end
    println("\nRaw RMSE by month ($judged): plain, gauges only, crossed; share of crossed fits choosing the gauges-only strength:")
    agreement = Dict((row.regime, row.month) => row.share_same
        for row in agreement_by_month(strengths) if row.scheme == judged)
    label_of(month) = Dates.format(month, "yyyy-mm")
    wide(name, month) = only(filter(r -> r.scheme == judged && r.series == name && r.month == month, monthly))
    for month in sort(unique(r.month for r in monthly))
        shares = [share for ((_, cell), share) in agreement if cell == label_of(month)]
        @printf("  %s  %.3f  %.3f  %.3f  %s\n", label_of(month), wide(PLAIN, month).RMSE, wide(GAUGES, month).RMSE,
            wide(CROSSED, month).RMSE, isempty(shares) ? "-" : @sprintf("%.2f", mean(shares)))
    end
    println("\nReading:")
    for row in reading
        @printf("  %-17s %-22s gauges only %+.4f  crossed %+.4f  moves %+.4f (limit %.2f)  %s%s\n", row.scheme,
            row.check, row.gauges_only, row.crossed, row.move, row.limit, row.material ? "MATERIAL" : "not material",
            row.judged ? "" : "  (reported only)")
    end
    any(row.material for row in reading if row.judged) ?
        println("Verdict: MATERIAL - the gain §3.4 cites depends on the shared weather; the user decides") :
        println("Verdict: not material - §3.4's raw gain survives holding out time")
    println("Wrote $(relpath(outdir, ROOT)) in $(round((time() - started) / 60; digits=1)) min")
end

if abspath(PROGRAM_FILE) == @__FILE__
    check_main()
end

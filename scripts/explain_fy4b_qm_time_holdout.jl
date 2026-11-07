#!/usr/bin/env julia

# Why does FY4B_QM's raw gain shrink when time is held out? `check_fy4b_qm_time_holdout.jl` found the
# change material at 456b9cf. This separates the possible causes, after that result and with no verdict
# of its own. No paper text or number changes.
#
#   julia --project=. scripts/explain_fy4b_qm_time_holdout.jl
#
# Balanced spatial only. Every version is scored at held-out gauges on the check's cells, against plain
# FY4B, with whole 7-day blocks resampled:
# - gauges: the benchmark's FY4B_QM, held to the run's `oof_raw.csv` in every cell;
# - crossed: the check's version, held to its `paired.csv`;
# - scattered: as crossed, with as many training hours, but the held-out hours are scattered at random
#   rather than in blocks. It sees as little data as crossed but keeps every weather system on both
#   sides. If it stays near gauges, less data is not why crossed falls;
# - map_held: the map fitted without the held-out blocks, with the strengths gauges chose;
# - strength_held: the gauges' map, with strengths chosen without the held-out blocks.
#   Together these show whether the map or the monthly strength loses the gain.
#
# Also compares the versions on the error budget's gross false-rain hours (FY4B's and FY4B_QM's) and
# on all other hours, and tabulates what each map does to large FY4B values, with and without each
# time fold, for one gauge fold.
#
# Writes output/fy4b_qm_time_holdout_explained/: paired.csv, difference_by_day.csv, map_values.csv and
# provenance.csv.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("DayBlockBootstrap", "RunProvenance", "TableIO", "SatelliteCalibration",
    "SatelliteTemporalEvaluation", "TimeBlockFolds")

using CSV, Printf, Random, Statistics
using Main.DayBlockBootstrap: read_wide_values, read_wide_mask, group_sums, compare_days
using Main.SatelliteCalibration: product_regimes, shift_fold, fit_quantile_map, choose_strength
using Main.SatelliteTemporalEvaluation: season_of
using Main.TimeBlockFolds: label_day, time_block_folds, training_hours
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const SC = Main.SatelliteCalibration
const REBUILD_DIR = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024")
const FY4B_FILE = joinpath(REBUILD_DIR, "fy4b_w3_h0.csv")
const HOUR_META = joinpath(REBUILD_DIR, "fy4b_hour_meta.csv")
const CHECK_FILE = joinpath(ROOT, "output", "fy4b_qm_time_holdout", "paired.csv")
const GROSS_FILE = joinpath(ROOT, "output", "fy4b_error_budget", "gross_hours.csv")
const OUTDIR = joinpath(ROOT, "output", "fy4b_qm_time_holdout_explained")
const SCHEME = :balanced_spatial
const REGIMES = Dict("133.0E_V1.0" => "R1", "133.0E_V1.0.0" => "R2", "105.0E_V1.0.0" => "R3")
const TIME_FOLD_SEED = 20261107
const SCATTER_SEED = 20261108
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261101
const THRESHOLDS = [0.1, 2.0, 8.0]
const CI_STATS = [:RMSE, :r, Symbol("CSI_0.1"), :POD_8]
const PLAIN, GAUGES, CROSSED, SCATTERED = "FY4B", "gauges", "crossed", "scattered"
const MAP_HELD, STRENGTH_HELD = "map_held", "strength_held"
const VERSIONS = [GAUGES, CROSSED, SCATTERED, MAP_HELD, STRENGTH_HELD]
const MAP_INPUTS = [2.0, 5.0, 10.0, 16.0, 20.0]

function explain_config()
    reference = benchmark_config(:full; nested_covariates=true, satellite_wet_blend=true,
        blend_axes=[:constant, :agreement_envelope], fused_anchor=true, gpm_uncal=true,
        fy4b=true, fy4b_file=FY4B_FILE, fy4b_qm=true, fy4b_qm_adaptive=true,
        fy4b_qm_fusion=true)
    fields = fieldnames(MGERConfig)
    values = NamedTuple{fields}(Tuple(getfield(reference.mger, field) for field in fields))
    return reference, MGERConfig(; merge(values, (; outdir=OUTDIR))...)
end

"""Each hour's scattered fold, and the hours each fold trains on: as many as crossed's, none its own."""
function scattered_folds(n_time, k, crossed_training_counts)
    rng = MersenneTwister(SCATTER_SEED)
    fold_of_hour = shuffle(rng, mod1.(1:n_time, k))
    train = map(1:k) do fold
        others = shuffle(rng, findall(!=(fold), fold_of_hour))
        sort(others[1:min(length(others), crossed_training_counts[fold])])
    end
    return fold_of_hour, train
end

"""`y` with every gauge set to NaN outside `train_hours`."""
function masked(y, train_hours)
    observed = fill(NaN, size(y))
    observed[:, train_hours] = y[:, train_hours]
    return observed
end

function explain_main()
    mkpath(OUTDIR)
    reference, mger = explain_config()
    scheme = string(SCHEME)
    run_dir = reference.mger.outdir
    path(product, file) = joinpath(run_dir, scheme, lowercase(product), file)
    grid_path = something(reference.fy4b_qm_grid)
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path,
        [mger.sat_paths[name] for name in sort(collect(keys(mger.sat_paths)))], grid_path,
        HOUR_META, CHECK_FILE, GROSS_FILE,
        [path(p, file) for p in (PLAIN, "FY4B_QM") for file in ("oof_raw.csv", "common_evaluation_mask.csv")])
    provenance = start_provenance(ROOT, OUTDIR; inputs)
    started = time()

    _, ids, product_data = load_global_common_product_data(mger)
    times = product_data[PLAIN].times
    y = Matrix{Float64}(product_data[PLAIN].Y_obs)
    station_meta = load_station_meta(mger.station_meta_path; station_id_col=mger.station_id_col,
        lon_col=mger.lon_col, lat_col=mger.lat_col)
    lonlat = build_X_lonlat(station_meta, ids)
    calibration = _prepare_fy4b_calibration!(reference, ids, product_data)
    id_map = Dict(id => index for (index, id) in enumerate(ids))
    k = reference.k
    labels, regime_index = product_regimes(HOUR_META, times)
    regime = [REGIMES[labels[index]] for index in regime_index]
    day = label_day.(times)
    blocks = time_block_folds(times, string.(regime, "_", season_of.(day)); seed=TIME_FOLD_SEED)
    block_group = repeat(permutedims(blocks.block_of_time), length(ids))
    n_block = length(blocks.fold_of_block)
    crossed_train = [findall(training_hours(blocks.fold_of_time, times, f)) for f in 1:k]
    fold_of_hour, scattered_train = scattered_folds(length(times), k, length.(crossed_train))
    @printf("Training hours per fold, of %d: %s (crossed and scattered alike: %s)\n", length(times),
        join(length.(crossed_train), ", "), length.(crossed_train) == length.(scattered_train))

    folds = benchmark_folds(SCHEME, ids, lonlat; k, seed=first(benchmark_seeds(reference)),
        center_init=reference.fold_center_init, rotation=0)
    series = Dict(name => fill(NaN, size(y)) for name in VERSIONS)
    cells, months = calibration.cells, calibration.month_of_time
    map_rows = NamedTuple[]
    fit_map(observed, shifted, train) = fit_quantile_map(observed, shifted, train, cells; wet=mger.rain_threshold)
    strengths(observed, shifted, maps, train) = choose_strength(observed, shifted, maps, train, cells, months)
    observed_crossed = [masked(y, hours) for hours in crossed_train]
    observed_scattered = [masked(y, hours) for hours in scattered_train]
    for gauge_fold in 1:k
        _, _, train, val = _fold_station_indices(folds, id_map, k, gauge_fold)
        shifted = first(shift_fold(calibration.grid, calibration.sums, train, calibration.regimes,
            calibration.regime_of_time; cache=calibration.cache))
        maps_all = fit_map(y, shifted, train)
        strengths_all = strengths(y, shifted, maps_all, train)
        value(maps, chosen, station, t) = begin
            m = maps[cells[t]]
            x = shifted[station, t]
            SC._adjusted_value(x, SC._map_value(m, x), m.tau, chosen[months[t]].candidate)
        end
        fill_in!(name, maps, chosen, hours) = for t in hours, i in val
            series[name][i, t] = value(maps, chosen, i, t)
        end
        fill_in!(GAUGES, maps_all, strengths_all, axes(y, 2))
        for f in 1:k
            test_hours = findall(==(f), blocks.fold_of_time)
            observed = observed_crossed[f]
            maps_held = fit_map(observed, shifted, train)
            fill_in!(CROSSED, maps_held, strengths(observed, shifted, maps_held, train), test_hours)
            fill_in!(MAP_HELD, maps_held, strengths_all, test_hours)
            fill_in!(STRENGTH_HELD, maps_all, strengths(observed, shifted, maps_all, train), test_hours)
            scattered = observed_scattered[f]
            maps_scattered = fit_map(scattered, shifted, train)
            fill_in!(SCATTERED, maps_scattered, strengths(scattered, shifted, maps_scattered, train),
                findall(==(f), fold_of_hour))
            if gauge_fold == 1
                for (cell, m) in sort(collect(maps_held); by=first)
                    all_map = maps_all[cell]
                    for x in MAP_INPUTS
                        push!(map_rows, (; gauge_fold, time_fold=f,
                            regime=REGIMES[calibration.regimes[(cell - 1) ÷ 4 + 1]],
                            season=SC.SEASONS[(cell - 1) % 4 + 1], fy4b=x, mapped_held=SC._map_value(m, x),
                            mapped_all=SC._map_value(all_map, x), top_knot_held=m.qx[end],
                            top_knot_all=all_map.qx[end], gauge_top_held=m.qy[end], gauge_top_all=all_map.qy[end]))
                    end
                end
            end
        end
        @printf("  gauge fold %d done (%.1f min)\n", gauge_fold, (time() - started) / 60)
    end
    stored = read_wide_values(path("FY4B_QM", "oof_raw.csv"), ids, times)
    mismatched = count(.!isequal.(stored, series[GAUGES]))
    mismatched == 0 || error("gauges differs from the run's oof_raw.csv in $mismatched cells")
    println("  gauges reproduces the run's oof_raw.csv in all $(length(stored)) cells")

    plain = read_wide_values(path(PLAIN, "oof_raw.csv"), ids, times)
    series[PLAIN] = plain
    mask = read_wide_mask(path("FY4B_QM", "common_evaluation_mask.csv"), ids, times) .&
        read_wide_mask(path(PLAIN, "common_evaluation_mask.csv"), ids, times) .& isfinite.(y)
    for name in vcat(PLAIN, VERSIONS)
        mask .&= isfinite.(series[name])
    end
    gross_times = Set(DateTime(first(string(t), 19)) for t in CSV.read(GROSS_FILE, DataFrame).time)
    gross_hour = [t in gross_times for t in times]
    @printf("  %d scored cells; %d gross hours on the grid\n", count(mask), count(gross_hour))

    sums(x, cells) = group_sums(x, y, cells, block_group, n_block, THRESHOLDS)
    compare(candidate, baseline, cells) = compare_days(sums(series[candidate], cells),
        sums(series[baseline], cells), THRESHOLDS; ci_stats=CI_STATS, reps=BOOTSTRAP_REPS, seed=BOOTSTRAP_SEED)
    paired = NamedTuple[]
    for (subset, cells) in (("all", mask), ("not_gross", mask .& permutedims(.!gross_hour)),
            ("gross", mask .& permutedims(gross_hour)))
        for (candidate, baseline) in vcat([(name, PLAIN) for name in VERSIONS], [(CROSSED, GAUGES), (SCATTERED, GAUGES)])
            push!(paired, merge((; subset, candidate, baseline), compare(candidate, baseline, cells)))
        end
    end
    expected = only(filter(r -> r.scheme == scheme && r.candidate == "FY4B_QM_crossed" && r.baseline == "FY4B",
        CSV.read(CHECK_FILE, DataFrame)))
    row = only(filter(r -> r.subset == "all" && r.candidate == CROSSED && r.baseline == PLAIN, paired))
    isapprox(row.candidate_RMSE, expected.candidate_RMSE; rtol=1e-9) && row.candidate_n == expected.candidate_n ||
        error("crossed does not reproduce $CHECK_FILE")
    @printf("  crossed reproduces the check: n %d, RMSE %.5f\n", row.candidate_n, row.candidate_RMSE)

    unique_days = sort(unique(day))
    position = Dict(d => i for (i, d) in enumerate(unique_days))
    sse_gauges, sse_crossed, sse_plain = (zeros(length(unique_days)) for _ in 1:3)
    for t in axes(y, 2), i in axes(y, 1)
        mask[i, t] || continue
        d = position[day[t]]
        sse_gauges[d] += (series[GAUGES][i, t] - y[i, t])^2
        sse_crossed[d] += (series[CROSSED][i, t] - y[i, t])^2
        sse_plain[d] += (plain[i, t] - y[i, t])^2
    end
    difference = sse_crossed .- sse_gauges
    order = sortperm(difference; rev=true)
    write_csv_atomic(joinpath(OUTDIR, "paired.csv"), DataFrame(paired))
    write_csv_atomic(joinpath(OUTDIR, "difference_by_day.csv"), DataFrame(; day=unique_days[order],
        sse_plain=sse_plain[order], sse_gauges=sse_gauges[order], sse_crossed=sse_crossed[order],
        difference=difference[order]))
    write_csv_atomic(joinpath(OUTDIR, "map_values.csv"), DataFrame(map_rows))
    CSV.write(joinpath(OUTDIR, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println("\nPaired (whole 7-day blocks, 95%):")
    for row in paired
        @printf("  %-9s %-13s vs %-7s n=%7d RMSE %.4f vs %.4f %+.2f%% (%+.2f..%+.2f)  r %.4f vs %.4f %+.4f (%+.4f..%+.4f)\n",
            row.subset, row.candidate, row.baseline, row.candidate_n, row.candidate_RMSE, row.baseline_RMSE,
            100 * row.rel_RMSE, 100 * row.rel_RMSE_lo, 100 * row.rel_RMSE_hi, row.candidate_r, row.baseline_r,
            row.delta_r, row.delta_r_lo, row.delta_r_hi)
    end
    total = sum(difference)
    @printf("\nCrossed minus gauges, squared error: net %.0f. Largest days (share of the net):\n", total)
    for index in order[1:12]
        @printf("  %s  plain %8.0f  gauges %8.0f  crossed %8.0f  difference %+8.0f (%5.1f%%)\n", unique_days[index],
            sse_plain[index], sse_gauges[index], sse_crossed[index], difference[index], 100 * difference[index] / total)
    end
    println("\nWhat the map does to large FY4B values, gauge fold 1 (all hours -> without each time fold):")
    for (regime_name, season_name) in unique((r.regime, r.season) for r in map_rows)
        for x in (10.0, 16.0)
            rows = filter(r -> r.regime == regime_name && r.season == season_name && r.fy4b == x, map_rows)
            @printf("  %s %s FY4B %2.0f -> %5.2f all hours; without a fold %s  (top knot %.1f; without %s)\n", regime_name,
                season_name, x, rows[1].mapped_all, join((@sprintf("%.2f", r.mapped_held) for r in rows), " "),
                rows[1].top_knot_all, join((@sprintf("%.1f", r.top_knot_held) for r in rows), " "))
        end
    end
    println("Wrote $(relpath(OUTDIR, ROOT)) in $(round((time() - started) / 60; digits=1)) min")
end

if abspath(PROGRAM_FILE) == @__FILE__
    explain_main()
end

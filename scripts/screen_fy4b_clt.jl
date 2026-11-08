#!/usr/bin/env julia

# Pilot screen: does FY4B's Cloud Type product (CLT) tell FY4B_QM's false rain from its real rain
# better than the rain amount does? One month, July 2023 (all regime 133E V1.0.0, all JJA), scored
# on held-out gauges. A screen, not a result: nothing downstream (MGWR, the blend) is refitted.
#
#   julia --project=. scripts/screen_fy4b_clt.jl
#   julia --project=. scripts/screen_fy4b_clt.jl <run dir> --matched-products
#
# A run folder reads the run without `--fy4b-qm-fusion` instead (`..._fy4bw3h0_fy4bqm_fy4bqmadapt_gpmuncal`,
# whose name must be what `benchmark_config` gives); FY4B_QM's raw field does not depend on the fusion, so
# the guard below holds it to that run as well. The output folder adds the run's options from
# `_fusedanchor` on. `--matched-products` also keeps only the cells every product of the run is scored on
# (`BenchmarkDiagnostics.matched_evaluation_mask`) and adds `_matched`.
#
# It reuses the benchmark's grid, folds and per-fold FY4B_QM, exactly as
# `screen_fy4b_two_stage.jl` does, with the same guard: the held-out FY4B_QM built here must
# equal the 2026-10-31 run's `oof_raw.csv` in every cell of both schemes.
#
# The CLT files are read from data/FY4B_CLT (a link to the external drive, kept apart from the
# QPE tree because `find_nc_files` keys on start time only). Each gauge-hour gets the dominant
# CLT class over the 3x3 box at that fold's FY4B_QM pixel shift and the hour's four quarters
# (`FY4BCloudType`). The CLT tree is not hashed; clt_hourly_qc.csv lists every file used.
#
# Three arms, fitted per fold on the training gauges' July cells and scored on held-out gauges:
# - A `FY4B_QM`;
# - B `FY4B_QM_CLT`: FY4B_QM set to 0 in the chosen classes. Greedy: add the class that raises
#   training CSI_0.1 the most, and stop when none does. A class needs `MIN_WET_PAIRS` training
#   FY4B_QM-wet pairs to enter.
# - C `FY4B_QM_AMT`, the amount-matched control: FY4B_QM set to 0 at or below a cut chosen so it
#   removes as many training wet hours as B does, the lightest first.
#
# The rule, fixed before the run (balanced spatial judged, random split reported only). The full
# CLT download is worth it only if
# 1. B against C: delta CSI_0.1 >= +0.01, with its 95% interval above 0, and
# 2. B against A: delta POD_8 not significantly negative (upper bound >= 0).
# If B picks no class in 3 or more of the 5 balanced folds, that is a stop as well.
#
# Writes output/fy4b_clt_screen/: clt_hourly_qc.csv, class_table.csv, class_by_amount.csv,
# diurnal.csv, filters.csv, metrics.csv, paired.csv, rule.csv and provenance.csv.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("FY4BCloudType", "DayBlockBootstrap", "RunProvenance", "TableIO",
    "BenchmarkDiagnostics")

using CSV, NCDatasets, Printf, Statistics
using Main.FY4BCloudType
using Main.FY4BPreprocessing: find_nc_files, parse_nc_filename
using Main.BenchmarkDiagnostics: matched_evaluation_mask, matched_fold_failures
using Main.DayBlockBootstrap: group_sums, sum_stats, day_groups, compare_days
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const FY4B_FILE = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024", "fy4b_w3_h0.csv")
const CLT_DIR = joinpath(ROOT, "data", "FY4B_CLT")
const QPE_DIR = joinpath(ROOT, "data", "FY4B")
const OUTDIR = joinpath(ROOT, "output", "fy4b_clt_screen")
const CLT_YEARS, CLT_MONTHS = 2023:2023, 7:7
const REGIME = "133.0E_V1.0.0"
const WINDOW = 3
const MIN_VALID = cld(WINDOW^2 * 4, 2)    # half of the hour's 36 pixel-quarters
const MIN_WET_PAIRS = 200
const WET = 0.1
const THRESHOLDS = [0.1, 2.0, 8.0]
const AMOUNT_EDGES = [0.1, 0.5, 2.0, 8.0, Inf]
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261101
const CI_STATS = [:RMSE, :r, :abs_bias, Symbol("CSI_0.1"), Symbol("POD_0.1"), Symbol("FAR_0.1"),
    Symbol("FBI_0.1"), :CSI_2, :POD_8, :CSI_8]
const SERIES = ["FY4B_QM", "FY4B_QM_CLT", "FY4B_QM_AMT", "GPM_uncal"]
"""(candidate, baseline) pairs compared in `paired.csv`."""
const PAIRS = [("FY4B_QM_CLT", "FY4B_QM"), ("FY4B_QM_CLT", "FY4B_QM_AMT"),
    ("FY4B_QM_AMT", "FY4B_QM")]
const RULE_CSI_GAIN = 0.01
const RULE_MAX_EMPTY_FOLDS = 2

"""
The configuration of the run being reproduced, and a copy of its `MGERConfig` that writes to
`outdir`: loading the grid writes `global_common_time_qc.csv`, which must not land in that run.
`run_dir` picks the run with or without `--fy4b-qm-fusion` (by its `_qmfusion` suffix) and must be
the folder that configuration names; without it, the `--fy4b-qm-fusion` run.
"""
function screen_config(outdir=OUTDIR; run_dir=nothing)
    qm_fusion = run_dir === nothing || occursin("_qmfusion", basename(rstrip(run_dir, '/')))
    reference = benchmark_config(:full; nested_covariates=true, satellite_wet_blend=true,
        blend_axes=[:constant, :agreement_envelope], fused_anchor=true, gpm_uncal=true,
        fy4b=true, fy4b_file=FY4B_FILE, fy4b_qm=true, fy4b_qm_adaptive=true,
        fy4b_qm_fusion=qm_fusion)
    run_dir === nothing || rstrip(abspath(run_dir), '/') == rstrip(reference.mger.outdir, '/') ||
        error("$run_dir is not a run this screen reads; the matching one is $(reference.mger.outdir)")
    fields = fieldnames(MGERConfig)
    values = NamedTuple{fields}(Tuple(getfield(reference.mger, field) for field in fields))
    return reference, MGERConfig(; merge(values, (; outdir))...)
end

"""The run's options from `_fusedanchor` on, which an output folder for that run carries."""
run_tag(run_dir) = last(split(basename(rstrip(run_dir, '/')), "_fusedanchor"; limit=2))

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

"""
Guard on the pixel alignment: on the first CLT quarter, the share of QPE-wet pixels (>= 1 mm/h,
accepted) that CLT calls clear, as read and with CLT transposed. Read as `[x, y]`, the two
products must agree far better than with the axes swapped.
"""
function check_alignment(groups)
    clt_file = first(groups[minimum(keys(groups))])
    start = first(parse_nc_filename(basename(clt_file)))
    qpe_groups = find_nc_files(QPE_DIR; years=year(start):year(start), months=month(start):month(start))
    qpe_file = only(filter(file -> first(parse_nc_filename(basename(file))) == start,
        qpe_groups[floor(start, Hour)]))
    rain, quality = NCDataset(qpe_file) do ds
        ds["Precipitation"].var[:, :], ds["DQF"].var[:, :]
    end
    clt = NCDataset(ds -> ds["CLT"].var[:, :], clt_file)
    wet = (rain .>= 1) .& (rain .<= 30) .& (quality .< 2)
    as_read = count(wet .& (clt .== 0)) / count(wet)
    transposed = count(wet .& (permutedims(clt) .== 0)) / count(wet)
    @printf("  Alignment: %d QPE-wet pixels at %s; CLT calls %.4f of them clear (%.4f transposed)\n",
        count(wet), start, as_read, transposed)
    as_read < transposed / 2 || error("CLT and QPE do not line up as [x, y]")
    return (; start, qpe_wet_pixels=count(wet), clear_share=as_read, clear_share_transposed=transposed)
end

"""Per class: training wet pairs, hits and false alarms at `WET`, and the misses overall."""
function training_counts(q, y, class, cells)
    wet_pairs = Dict{Int8,Int}(); hits = Dict{Int8,Int}(); false_alarms = Dict{Int8,Int}()
    misses = 0
    for index in cells
        c = class[index]
        observed, estimated = y[index] >= WET, q[index] >= WET
        if estimated
            wet_pairs[c] = get(wet_pairs, c, 0) + 1
            observed ? (hits[c] = get(hits, c, 0) + 1) :
                (false_alarms[c] = get(false_alarms, c, 0) + 1)
        elseif observed
            misses += 1
        end
    end
    return wet_pairs, hits, false_alarms, misses
end

"""
The classes the greedy filter zeroes, in the order they were added, and the training CSI before
and after. Zeroing class c turns its hits into misses and its false alarms into dry hours.
"""
function choose_classes(wet_pairs, hits, false_alarms, misses)
    H, F = sum(values(hits); init=0), sum(values(false_alarms); init=0)
    csi(set) = begin
        h = H - sum(get(hits, c, 0) for c in set; init=0)
        f = F - sum(get(false_alarms, c, 0) for c in set; init=0)
        h / (H + misses + f)
    end
    chosen = Int8[]
    current = csi(chosen)
    start = current
    candidates = sort([c for (c, n) in wet_pairs if n >= MIN_WET_PAIRS])
    while true
        remaining = setdiff(candidates, chosen)
        isempty(remaining) && break
        scores = [csi(vcat(chosen, c)) for c in remaining]
        best = argmax(scores)
        scores[best] > current || break
        push!(chosen, remaining[best])
        current = scores[best]
    end
    return chosen, start, current
end

"""
The amount cut that removes closest to `target` of the training wet values `wet_values` (zeroing
every value at or below it), the smaller cut on a tie; `-Inf` when `target` is 0.
"""
function matched_cut(wet_values, target)
    target == 0 && return -Inf, 0
    sorted = sort(wet_values)
    best_cut, best_removed = -Inf, 0
    for cut in unique(sorted)
        removed = searchsortedlast(sorted, cut)
        abs(removed - target) < abs(best_removed - target) && ((best_cut, best_removed) = (cut, removed))
        removed >= target && break
    end
    return best_cut, best_removed
end

_class_name(code) = code < 0 ? "none" : CLT_CLASS_NAMES[findfirst(==(code), CLT_CLASSES)]

"""Event counts at `WET` for one series over `cells`: wet hours, hits, false alarms."""
function wet_counts(x, y, cells)
    wet = hits = false_alarms = 0
    for index in cells
        x[index] >= WET || continue
        wet += 1
        y[index] >= WET ? (hits += 1) : (false_alarms += 1)
    end
    return wet, hits, false_alarms
end

function class_rows(scheme, q, gpm, y, class, mask)
    rows = NamedTuple[]
    all_cells = findall(mask)
    _, total_hits, total_fa = wet_counts(q, y, all_cells)
    for code in Int8.(CLT_CLASSES)
        cells = findall(mask .& (class .== code))
        isempty(cells) && continue
        wet, hits, fa = wet_counts(q, y, cells)
        gpm_wet, gpm_hits, gpm_fa = wet_counts(gpm, y, cells)
        push!(rows, (; scheme, class=Int(code), name=_class_name(code), n_cells=length(cells),
            share_cells=length(cells) / length(all_cells),
            gauge_wet_fraction=count(>=(WET), y[cells]) / length(cells),
            gauge_mean=mean(y[cells]), qm_wet=wet, qm_hits=hits, qm_false_alarms=fa,
            qm_FAR=fa / wet, share_of_false_alarms=fa / total_fa, share_of_hits=hits / total_hits,
            gpm_wet, gpm_hits, gpm_false_alarms=gpm_fa, gpm_FAR=gpm_fa / gpm_wet))
    end
    return rows
end

function amount_rows(scheme, q, y, class, mask)
    rows = NamedTuple[]
    for code in Int8.(CLT_CLASSES), b in 1:(length(AMOUNT_EDGES) - 1)
        lo, hi = AMOUNT_EDGES[b], AMOUNT_EDGES[b + 1]
        cells = findall(mask .& (class .== code) .& (q .>= lo) .& (q .< hi))
        isempty(cells) && continue
        false_alarms = count(<(WET), y[cells])
        push!(rows, (; scheme, class=Int(code), name=_class_name(code), qm_lo=lo, qm_hi=hi,
            n_wet=length(cells), false_alarms, FAR=false_alarms / length(cells),
            gauge_mean=mean(y[cells])))
    end
    return rows
end

function diurnal_rows(scheme, q, y, class, mask, times)
    rows = NamedTuple[]
    for h in 0:23, code in Int8.(CLT_CLASSES)
        columns = hour.(times) .== h
        cells = findall(mask .& (class .== code) .& permutedims(columns))
        wet, hits, fa = wet_counts(q, y, cells)
        wet == 0 && continue
        push!(rows, (; scheme, hour_end_bjt=h, class=Int(code), name=_class_name(code), qm_wet=wet,
            qm_hits=hits, qm_false_alarms=fa))
    end
    return rows
end

_stat(row, name) = getproperty(row, Symbol(name))

function screen_main(args=ARGS)
    matched = "--matched-products" in args
    positional = filter(!startswith("--"), args)
    chosen = isempty(positional) ? nothing : abspath(only(positional))
    outdir = OUTDIR * (chosen === nothing ? "" : run_tag(chosen)) * (matched ? "_matched" : "")
    mkpath(outdir)
    reference, mger = screen_config(outdir; run_dir=chosen)
    oof_qm(scheme) = joinpath(reference.mger.outdir, scheme, lowercase(FY4B_QM_PRODUCT),
        "oof_raw.csv")
    grid_path = something(reference.fy4b_qm_grid)
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path,
        [mger.sat_paths[name] for name in sort(collect(keys(mger.sat_paths)))], grid_path,
        joinpath(dirname(grid_path), "fy4b_hour_meta.csv"),
        [oof_qm(string(scheme)) for scheme in reference.cv_schemes],
        matched ? [joinpath(reference.mger.outdir, "run_status.csv")] : String[])
    provenance = start_provenance(ROOT, outdir; inputs)
    started = time()

    isdir(CLT_DIR) || error("CLT tree not found at $CLT_DIR - is the external drive mounted?")
    groups, clt_qc = find_clt_files(CLT_DIR; years=CLT_YEARS, months=CLT_MONTHS)
    write_csv_atomic(joinpath(outdir, "clt_hourly_qc.csv"), clt_qc)
    println("CLT: $(nrow(clt_qc)) hours, $(count(clt_qc.complete)) complete, " *
        "$(sum(clt_qc.file_count)) files")
    alignment = check_alignment(groups)

    _, ids, product_data = load_global_common_product_data(mger)
    times = product_data["FY4B"].times
    y = product_data["FY4B"].Y_obs
    station_meta = load_station_meta(mger.station_meta_path; station_id_col=mger.station_id_col,
        lon_col=mger.lon_col, lat_col=mger.lat_col)
    lonlat = build_X_lonlat(station_meta, ids)
    calibration = _prepare_fy4b_calibration!(reference, ids, product_data)
    id_map = Dict(id => index for (index, id) in enumerate(ids))

    # The CLT-complete hours on the FY4B grid, all in one regime.
    complete_times = Set(clt_qc.time[clt_qc.complete])
    columns = findall(in(complete_times), times)
    clt_month = [year(t) in CLT_YEARS && month(t) in CLT_MONTHS for t in times]
    println("Grid: $(length(ids)) stations x $(length(times)) hours; $(length(columns)) of " *
        "$(length(complete_times)) CLT-complete hours are on it; $(count(clt_month)) grid hours " *
        "are labelled $(CLT_YEARS) month $(CLT_MONTHS)")
    all(calibration.regimes[calibration.regime_of_time[columns]] .== REGIME) ||
        error("the CLT hours are not all in regime $REGIME")
    sub_times = times[columns]
    y_sub = y[:, columns]
    gpm_sub = product_data["GPM_uncal"].Y_sat[:, columns]

    # Per scheme and fold: the held-out FY4B_QM (checked against the run) and the fold's shift.
    held_out = Dict{String,Matrix{Float64}}()
    fold_info = Dict{String,Vector{NamedTuple}}()
    for scheme_symbol in reference.cv_schemes
        scheme = string(scheme_symbol)
        folds = benchmark_folds(scheme_symbol, ids, lonlat; k=reference.k,
            seed=first(benchmark_seeds(reference)), center_init=reference.fold_center_init,
            rotation=0)
        held = fill(NaN, size(y))
        fold_info[scheme] = NamedTuple[]
        for fold in 1:reference.k
            _, _, train, val = _fold_station_indices(folds, id_map, reference.k, fold)
            # Exactly the benchmark's call, so the guard below can hold it to that run.
            qm, shift_rows, _, _ = calibrate_fold(
                calibration.grid, calibration.sums, y, train, calibration.regimes,
                calibration.regime_of_time, calibration.cells;
                cache=calibration.cache, wet=mger.rain_threshold,
                adaptive=reference.fy4b_qm_adaptive, choice_of_time=calibration.month_of_time,
                choice_labels=calibration.month_labels,
            )
            held[val, :] = qm[val, :]
            shift = only(filter(row -> row.regime == REGIME, shift_rows))
            push!(fold_info[scheme], (; fold, train, val, qm=qm[:, columns],
                shift=(Int(shift.dx), Int(shift.dy))))
        end
        check_reproduction(oof_qm(scheme), held, times, ids)
        held_out[scheme] = held[:, columns]
    end

    shifts = sort(unique(info.shift for infos in values(fold_info) for info in infos))
    stations = DataFrame(lat=lonlat[:, 2], lon=lonlat[:, 1])
    println("Reading CLT at shifts $(shifts) ($(round((time() - started) / 60; digits=1)) min)")
    clt = clt_class_counts(groups, stations, shifts; window=WINDOW)
    dominant = dominant_class(clt.counts; min_valid=MIN_VALID)
    clt_column = Dict(t => h for (h, t) in enumerate(clt.time))
    to_clt = [clt_column[t] for t in sub_times]
    println("  done ($(round((time() - started) / 60; digits=1)) min)")

    filter_rows = NamedTuple[]
    class_table = NamedTuple[]
    amount_table = NamedTuple[]
    diurnal_table = NamedTuple[]
    metric_rows = NamedTuple[]
    paired = NamedTuple[]
    for scheme_symbol in reference.cv_schemes
        scheme = string(scheme_symbol)
        q = held_out[scheme]
        series = Dict(name => fill(NaN, size(q)) for name in ("FY4B_QM_CLT", "FY4B_QM_AMT"))
        held_class = fill(Int8(-1), size(q))
        for info in fold_info[scheme]
            class = dominant[:, findfirst(==(info.shift), shifts), to_clt]
            usable = isfinite.(y_sub) .& isfinite.(info.qm) .& isfinite.(gpm_sub) .& (class .>= 0)
            train_cells = [CartesianIndex(i, t) for i in info.train, t in axes(q, 2)
                if usable[i, t]]
            wet_pairs, hits, false_alarms, misses = training_counts(info.qm, y_sub, class,
                train_cells)
            chosen, csi_before, csi_after = choose_classes(wet_pairs, hits, false_alarms, misses)
            removed_train = sum(get(wet_pairs, c, 0) for c in chosen; init=0)
            wet_values = [info.qm[index] for index in train_cells if info.qm[index] >= WET]
            cut, removed_amount = matched_cut(wet_values, removed_train)

            for i in info.val, t in axes(q, 2)
                x = q[i, t]
                series["FY4B_QM_CLT"][i, t] = class[i, t] in chosen ? 0.0 : x
                series["FY4B_QM_AMT"][i, t] = x <= cut ? 0.0 : x
                held_class[i, t] = class[i, t]
            end
            val_cells = [CartesianIndex(i, t) for i in info.val, t in axes(q, 2) if usable[i, t]]
            removed_val_clt = count(index -> q[index] >= WET && class[index] in chosen, val_cells)
            removed_val_amt = count(index -> WET <= q[index] <= cut, val_cells)
            push!(filter_rows, (; scheme, fold=info.fold, n_train_station=length(info.train),
                dx=info.shift[1], dy=info.shift[2],
                chosen=join(_class_name.(chosen), ";"), chosen_codes=join(Int.(chosen), ";"),
                train_csi_before=csi_before, train_csi_after=csi_after,
                train_wet=length(wet_values), train_removed_clt=removed_train,
                amount_cut=cut, train_removed_amount=removed_amount,
                val_wet=count(index -> q[index] >= WET, val_cells),
                val_removed_clt=removed_val_clt, val_removed_amount=removed_val_amt))
        end
        series["FY4B_QM"] = q
        series["GPM_uncal"] = gpm_sub
        mask = isfinite.(y_sub) .& (held_class .>= 0)
        for name in SERIES
            mask .&= isfinite.(series[name])
        end
        if matched
            run_dir = reference.mger.outdir
            products = sort(unique(String.(CSV.read(joinpath(run_dir, "metrics_pooled.csv"),
                DataFrame; select=[:product]).product)))
            failures = matched_fold_failures(CSV.read(joinpath(run_dir, "run_status.csv"), DataFrame);
                scheme, products)
            isempty(failures) || error("--matched-products: $scheme folds failed:\n" *
                sprint(show, select(failures, :product, :fold, :method, :status)))
            mask .&= matched_evaluation_mask(joinpath(run_dir, scheme), products, ids,
                times).matched[:, columns]
        end

        append!(class_table, class_rows(scheme, q, gpm_sub, y_sub, held_class, mask))
        append!(amount_table, amount_rows(scheme, q, y_sub, held_class, mask))
        append!(diurnal_table, diurnal_rows(scheme, q, y_sub, held_class, mask, sub_times))
        one_group = ones(Int, size(q))
        for name in SERIES
            total = vec(group_sums(series[name], y_sub, mask, one_group, 1, THRESHOLDS))
            push!(metric_rows, merge((; scheme, series=name), sum_stats(total, THRESHOLDS)))
        end
        day_group, n_day = day_groups(sub_times, size(q, 1))
        sums(x) = group_sums(x, y_sub, mask, day_group, n_day, THRESHOLDS)
        for (candidate, baseline) in PAIRS
            row = compare_days(sums(series[candidate]), sums(series[baseline]), THRESHOLDS;
                ci_stats=CI_STATS, reps=BOOTSTRAP_REPS, seed=BOOTSTRAP_SEED)
            push!(paired, merge((; scheme, candidate, baseline, n_day), row))
        end
    end

    # The rule, on the balanced spatial scheme.
    judged = "balanced_spatial"
    pick(c, b) = only(filter(r -> r.scheme == judged && r.candidate == c && r.baseline == b, paired))
    vs_amount = pick("FY4B_QM_CLT", "FY4B_QM_AMT")
    vs_qm = pick("FY4B_QM_CLT", "FY4B_QM")
    empty_folds = count(r -> r.scheme == judged && isempty(r.chosen), filter_rows)
    rule = [
        (; check="B vs C: delta CSI_0.1 >= $RULE_CSI_GAIN and interval above 0",
            value=_stat(vs_amount, "delta_CSI_0.1"), lo=_stat(vs_amount, "delta_CSI_0.1_lo"),
            hi=_stat(vs_amount, "delta_CSI_0.1_hi"),
            pass=_stat(vs_amount, "delta_CSI_0.1") >= RULE_CSI_GAIN &&
                _stat(vs_amount, "delta_CSI_0.1_lo") > 0),
        (; check="B vs A: delta POD_8 not significantly negative",
            value=vs_qm.delta_POD_8, lo=vs_qm.delta_POD_8_lo, hi=vs_qm.delta_POD_8_hi,
            pass=vs_qm.delta_POD_8_hi >= 0),
        (; check="B picks a class in at least $(reference.k - RULE_MAX_EMPTY_FOLDS) of $(reference.k) folds",
            value=Float64(reference.k - empty_folds), lo=NaN, hi=NaN,
            pass=empty_folds <= RULE_MAX_EMPTY_FOLDS),
    ]

    write_csv_atomic(joinpath(outdir, "filters.csv"), DataFrame(filter_rows))
    write_csv_atomic(joinpath(outdir, "class_table.csv"), DataFrame(class_table))
    write_csv_atomic(joinpath(outdir, "class_by_amount.csv"), DataFrame(amount_table))
    write_csv_atomic(joinpath(outdir, "diurnal.csv"), DataFrame(diurnal_table))
    write_csv_atomic(joinpath(outdir, "metrics.csv"), DataFrame(metric_rows))
    write_csv_atomic(joinpath(outdir, "paired.csv"), DataFrame(paired))
    write_csv_atomic(joinpath(outdir, "rule.csv"), DataFrame(rule))
    write_csv_atomic(joinpath(outdir, "alignment.csv"), DataFrame([alignment]))
    CSV.write(joinpath(outdir, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println("\nFilters (training stations):")
    for row in filter_rows
        @printf("  %-17s fold %d shift (%+d,%+d)  chosen %-24s CSI %.4f -> %.4f  removed %5d / %6d  cut %.3f (%d)\n",
            row.scheme, row.fold, row.dx, row.dy, isempty(row.chosen) ? "-" : row.chosen,
            row.train_csi_before, row.train_csi_after, row.train_removed_clt, row.train_wet,
            row.amount_cut, row.train_removed_amount)
    end
    println("\nHeld-out FY4B_QM by dominant class:")
    for row in class_table
        @printf("  %-17s %-12s cells %6d  wet %5d  FAR %.3f (GPM_uncal %.3f)  FA share %.3f  hit share %.3f\n",
            row.scheme, row.name, row.n_cells, row.qm_wet, row.qm_FAR, row.gpm_FAR,
            row.share_of_false_alarms, row.share_of_hits)
    end
    println("\nHeld-out scores:")
    for row in metric_rows
        @printf("  %-17s %-12s n=%7d  r %.3f  RMSE %.3f  bias %+.4f  CSI %.3f/%.3f/%.3f  POD_8 %.3f  FAR_0.1 %.3f  FBI_0.1 %.2f\n",
            row.scheme, row.series, row.n, row.r, row.RMSE, row.bias, _stat(row, "CSI_0.1"),
            row.CSI_2, row.CSI_8, row.POD_8, _stat(row, "FAR_0.1"), _stat(row, "FBI_0.1"))
    end
    println("\nPaired (day-block bootstrap, 95%):")
    for row in paired
        @printf("  %-17s %-12s vs %-12s CSI_0.1 %+.4f (%+.4f..%+.4f)  POD_8 %+.4f (%+.4f..%+.4f)  RMSE %+.2f%% (%+.2f..%+.2f)  r %+.4f\n",
            row.scheme, row.candidate, row.baseline, _stat(row, "delta_CSI_0.1"),
            _stat(row, "delta_CSI_0.1_lo"), _stat(row, "delta_CSI_0.1_hi"), row.delta_POD_8,
            row.delta_POD_8_lo, row.delta_POD_8_hi, 100 * row.rel_RMSE, 100 * row.rel_RMSE_lo,
            100 * row.rel_RMSE_hi, row.delta_r)
    end
    println("\nRule ($judged):")
    for row in rule
        @printf("  %-4s %s: %+.4f (%+.4f..%+.4f)\n", row.pass ? "pass" : "FAIL", row.check,
            row.value, row.lo, row.hi)
    end
    verdict = all(row.pass for row in rule) ? "GO: the full CLT download is worth it" :
        "STOP: the rule fails (see the checks above); FY4B_QM stays the quality step"
    println("Verdict: $verdict")
    println("Wrote $(relpath(outdir, ROOT)) in $(round((time() - started) / 60; digits=1)) min")
end

if abspath(PROGRAM_FILE) == @__FILE__
    screen_main()
end

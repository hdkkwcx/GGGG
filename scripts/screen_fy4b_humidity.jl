#!/usr/bin/env julia

# Screen: does near-surface humidity tell FY4B_QM's false rain from its real rain better than the
# rain amount does? All grid hours, scored on held-out gauges. A screen, not a result: nothing
# downstream (MGWR, the blend) is refitted. ERA5-Land is reanalysis, not FY4B, so anything kept is
# "FY4B + reanalysis" and could not run in real time.
#
#   julia --project=. scripts/screen_fy4b_humidity.jl
#
# It reuses the benchmark's grid, folds and per-fold FY4B_QM, exactly as `screen_fy4b_clt.jl` does,
# with the same guard: the held-out FY4B_QM built here must equal the 2026-10-31 run's `oof_raw.csv`
# in every cell of both schemes.
#
# Humidity is the relative humidity of ERA5-Land's 2 m temperature and dew point, read as the
# benchmark reads its covariates: `load_era5_panel` at `time_utc + 9 h`, so an hour labelled L by its
# Beijing end takes the value at its start. The script errors unless the alignment holds:
# - mean t2m over the label's hour of day peaks within 13-17 BJT;
# - RH's anomaly from its station x month x hour-of-day mean correlates best with the gauge's wet
#   call at a lag within -2..+4 h of that alignment (alignment.csv has the scan).
#
# Fitted per fold on the training gauges, over all grid hours, and scored on held-out gauges:
# - A `FY4B_QM`;
# - B `FY4B_QM_RH`: FY4B_QM set to 0 where RH < cut. The cut is the one in `CUTS` that raises
#   training CSI_0.1 the most, or no gate if none does;
# - C `FY4B_QM_AMT`, B's amount-matched control: FY4B_QM set to 0 at or below a cut chosen so it
#   removes as many training wet hours as B does, the lightest first;
# - B8 `FY4B_QM_RH8` and C8 `FY4B_QM_AMT8`: as B and C, but the gate leaves values >= 8 mm/h alone;
# - D `GPM_uncal_RH`: the same gate fitted on GPM_uncal. A control for whether humidity is specific
#   to FY4B; in no rule.
#
# The rule, fixed before the run, is the CLT screen's (balanced spatial judged, random split reported
# only). For B against its control C, and separately B8 against C8:
# 1. delta CSI_0.1 >= +0.01, with its 95% interval above 0, and
# 2. against A, delta POD_8 not significantly negative (upper bound >= 0), and
# 3. a gate is chosen in at least 3 of the 5 balanced folds.
# B8 meets 2 by construction. GO if B or B8 meets all three; it would justify a benchmark run, as
# nothing here shows what survives MGWR and the blend. `d2m_c` is already a local MGWR covariate.
#
# Writes output/fy4b_humidity_screen/: alignment.csv, cuts.csv, rh_table.csv, metrics.csv,
# paired.csv, strata.csv, gross_effect.csv (the error-budget's gross false-rain hours), rule.csv and
# provenance.csv.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("DayBlockBootstrap", "RunProvenance", "TableIO", "ERA5VariableSelection",
    "ERA5LandProcessing", "SatelliteCalibration", "SatelliteTemporalEvaluation")

using CSV, Printf, Statistics
using Main.DayBlockBootstrap: group_sums, sum_stats, day_groups, compare_days
using Main.ERA5VariableSelection: load_era5_panel
using Main.ERA5LandProcessing: relative_humidity
using Main.SatelliteCalibration: product_regimes
using Main.SatelliteTemporalEvaluation: SEASONS, season_of
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const REBUILD_DIR = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024")
const FY4B_FILE = joinpath(REBUILD_DIR, "fy4b_w3_h0.csv")
const HOUR_META = joinpath(REBUILD_DIR, "fy4b_hour_meta.csv")
const GROSS_FILE = joinpath(ROOT, "output", "fy4b_error_budget", "gross_hours.csv")
const OUTDIR = joinpath(ROOT, "output", "fy4b_humidity_screen")
const ERA5_PATHS = Dict(year_value => joinpath(ROOT, "data", "processed", "covariates", "era5_land",
    "era5_land_station_hourly_utc_$(year_value).csv") for year_value in 2022:2024)
const OFFSET_HOURS = 9
const T2M_PEAK_HOURS = 13:17
const LAG_SCAN = -12:12
const LAG_ALLOWED = -2:4
const CUTS = 40.0:5.0:90.0
const HEAVY = 8.0
const WET = 0.1
const THRESHOLDS = [0.1, 2.0, 8.0]
const RH_EDGES = [0, 50, 60, 70, 80, 90, 95, 101]
const REGIMES = Dict("133.0E_V1.0" => "R1", "133.0E_V1.0.0" => "R2", "105.0E_V1.0.0" => "R3")
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261106
const CI_STATS = [:RMSE, :r, :abs_bias, Symbol("CSI_0.1"), Symbol("POD_0.1"), Symbol("FAR_0.1"),
    Symbol("FBI_0.1"), :CSI_2, :POD_8, :CSI_8]
const QM, GATE, CONTROL = "FY4B_QM", "FY4B_QM_RH", "FY4B_QM_AMT"
const GATE8, CONTROL8 = "FY4B_QM_RH8", "FY4B_QM_AMT8"
const GPM, GPM_GATE = "GPM_uncal", "GPM_uncal_RH"
const SERIES = [QM, GATE, CONTROL, GATE8, CONTROL8, GPM, GPM_GATE]
"""(candidate, baseline) pairs compared in `paired.csv`."""
const PAIRS = [(GATE, QM), (GATE, CONTROL), (CONTROL, QM), (GATE8, QM), (GATE8, CONTROL8),
    (CONTROL8, QM), (GPM_GATE, GPM)]
const STRATA_PAIRS = [(GATE, QM), (GATE, CONTROL), (GATE8, CONTROL8), (GPM_GATE, GPM)]
const RULE_CSI_GAIN = 0.01
const RULE_MAX_EMPTY_FOLDS = 2

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

"""
RH and t2m over every UTC hour of the ERA5 files, stations x hours, and those hours' labels: a
value's label is `time_utc + OFFSET_HOURS`, the benchmark's alignment.
"""
function load_humidity(ids)
    labels = collect(DateTime(2022, 1, 1, OFFSET_HOURS):Hour(1):(DateTime(2025, 1, 1) + Hour(OFFSET_HOURS - 1)))
    panel = load_era5_panel(ERA5_PATHS, ids, labels; variables=[:t2m_c, :d2m_c],
        offset_hours=OFFSET_HOURS)
    t2m, d2m = panel.values[:t2m_c], panel.values[:d2m_c]
    return labels, relative_humidity.(t2m .+ 273.15, d2m .+ 273.15), t2m
end

"""RH minus its mean for that station, month and hour of day, over the whole record."""
function humidity_anomaly(rh, labels)
    group = [(month(t) - 1) * 24 + hour(t) + 1 for t in labels]
    anomaly = similar(rh)
    for i in axes(rh, 1)
        totals, counts = zeros(288), zeros(Int, 288)
        for t in axes(rh, 2)
            totals[group[t]] += rh[i, t]
            counts[group[t]] += 1
        end
        for t in axes(rh, 2)
            anomaly[i, t] = rh[i, t] - totals[group[t]] / counts[group[t]]
        end
    end
    return anomaly
end

"""
The alignment checks in the header: t2m's mean by the label's hour of day, and the correlation of
the gauge's wet call with RH's anomaly at each lag. Errors if either falls outside its range.
"""
function check_alignment(y, times, labels, rh, t2m)
    position = Dict(label => index for (index, label) in enumerate(labels))
    columns = [position[t] for t in times]
    rows = NamedTuple[]
    hour_mean = [mean(view(t2m, :, columns[hour.(times) .== h])) for h in 0:22]
    for (h, value) in zip(0:22, hour_mean)
        push!(rows, (; check="t2m_by_hour_end_bjt", key=h, value))
    end
    peak_hour = (0:22)[argmax(hour_mean)]
    anomaly = humidity_anomaly(rh, labels)
    wet = Float64.(y .>= WET)
    finite = isfinite.(y)
    scan = Float64[]
    for lag in LAG_SCAN
        shifted = [get(position, t + Hour(lag), 0) for t in times]
        usable = findall(>(0), shifted)
        cells = finite[:, usable]
        value = cor(wet[:, usable][cells], anomaly[:, shifted[usable]][cells])
        push!(scan, value)
        push!(rows, (; check="wet_vs_rh_anomaly_at_lag", key=lag, value))
    end
    peak_lag = LAG_SCAN[argmax(scan)]
    println("  Alignment: t2m peaks at label $peak_hour BJT; the gauge's wet call best matches RH's " *
        "anomaly at lag $peak_lag h (r $(round(maximum(scan); digits=3)); lag 0 r " *
        "$(round(scan[findfirst(==(0), LAG_SCAN)]; digits=3)))")
    peak_hour in T2M_PEAK_HOURS || error("t2m peaks at label $peak_hour BJT: the alignment is wrong")
    peak_lag in LAG_ALLOWED || error("RH's anomaly peaks at lag $peak_lag h: the alignment is wrong")
    return rows, columns
end

"""
The cut in `CUTS` that raises CSI_0.1 the most over `cells` when `x` is zeroed where RH is below it
and `x < spared`. Returns the cut (`-Inf` for no gate), the CSI before and after, and the wet cells
it removes.
"""
function choose_cut(x, y, rh, cells; spared=Inf)
    wet = [index for index in cells if x[index] >= WET]
    hits = count(index -> y[index] >= WET, wet)
    false_alarms = length(wet) - hits
    misses = count(index -> y[index] >= WET && x[index] < WET, cells)
    before = hits / (hits + misses + false_alarms)
    best_cut, best_csi, best_removed = -Inf, before, 0
    for cut in CUTS
        removed_hits = removed_false = 0
        for index in wet
            (rh[index] < cut && x[index] < spared) || continue
            y[index] >= WET ? (removed_hits += 1) : (removed_false += 1)
        end
        csi = (hits - removed_hits) / (hits + misses + false_alarms - removed_false)
        csi > best_csi && ((best_cut, best_csi, best_removed) = (cut, csi, removed_hits + removed_false))
    end
    return best_cut, before, best_csi, best_removed
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

"""Per RH bin, P(gauge wet) by the product's wet call and its shares of false alarms and hits."""
function rh_rows(scheme, product, x, y, rh, mask)
    observed, estimate, humidity = y[mask], x[mask], rh[mask]
    wet = estimate .>= WET
    gauge_wet = observed .>= WET
    total_false, total_hits = count(wet .& .!gauge_wet), count(wet .& gauge_wet)
    return [begin
        lo, hi = RH_EDGES[b], RH_EDGES[b + 1]
        bin = lo .<= humidity .< hi
        (; scheme, product, rh_lo=lo, rh_hi=hi, n_cells=count(bin), gauge_wet=mean(gauge_wet[bin]),
            gauge_wet_when_wet=mean(gauge_wet[bin .& wet]), gauge_wet_when_dry=mean(gauge_wet[bin .& .!wet]),
            share_of_false_alarms=count(bin .& wet .& .!gauge_wet) / total_false,
            share_of_hits=count(bin .& wet .& gauge_wet) / total_hits)
    end for b in 1:(length(RH_EDGES) - 1)]
end

"""Per series, false-rain squared error and wet cells in the gross hours and in all hours."""
function gross_rows(scheme, series, y, mask, gross_hour)
    rows = NamedTuple[]
    baseline = Dict{String,Float64}()
    for name in (QM, GATE, CONTROL, GATE8, CONTROL8)
        x = series[name]
        totals = Dict("gross" => [0.0, 0.0, 0.0], "other" => [0.0, 0.0, 0.0])
        for t in axes(y, 2), i in axes(y, 1)
            mask[i, t] || continue
            bucket = totals[gross_hour[t] ? "gross" : "other"]
            if x[i, t] >= WET
                y[i, t] >= WET ? (bucket[3] += 1) : (bucket[1] += (x[i, t] - y[i, t])^2; bucket[2] += 1)
            end
        end
        for (hours, (sse, false_alarms, hits)) in totals
            name == QM && (baseline[hours] = sse)
            push!(rows, (; scheme, series=name, hours, false_rain_sse=sse,
                share_of_qm_false_rain_sse=sse / baseline[hours], false_alarms=Int(false_alarms),
                hits=Int(hits)))
        end
    end
    return rows
end

_stat(row, name) = getproperty(row, Symbol(name))

"""The three checks of the rule for one gate against its control, on the judged scheme."""
function rule_rows(arm, vs_control, vs_qm, empty_folds, k)
    return [
        (; arm, check="vs control: delta CSI_0.1 >= $RULE_CSI_GAIN and interval above 0",
            value=_stat(vs_control, "delta_CSI_0.1"), lo=_stat(vs_control, "delta_CSI_0.1_lo"),
            hi=_stat(vs_control, "delta_CSI_0.1_hi"),
            pass=_stat(vs_control, "delta_CSI_0.1") >= RULE_CSI_GAIN &&
                _stat(vs_control, "delta_CSI_0.1_lo") > 0),
        (; arm, check="vs FY4B_QM: delta POD_8 not significantly negative",
            value=vs_qm.delta_POD_8, lo=vs_qm.delta_POD_8_lo, hi=vs_qm.delta_POD_8_hi,
            pass=vs_qm.delta_POD_8_hi >= 0),
        (; arm, check="a gate in at least $(k - RULE_MAX_EMPTY_FOLDS) of $k folds",
            value=Float64(k - empty_folds), lo=NaN, hi=NaN, pass=empty_folds <= RULE_MAX_EMPTY_FOLDS),
    ]
end

function screen_main()
    mkpath(OUTDIR)
    reference, mger = screen_config()
    oof_qm(scheme) = joinpath(reference.mger.outdir, scheme, lowercase(QM), "oof_raw.csv")
    grid_path = something(reference.fy4b_qm_grid)
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path,
        [mger.sat_paths[name] for name in sort(collect(keys(mger.sat_paths)))], grid_path,
        HOUR_META, GROSS_FILE, [ERA5_PATHS[year_value] for year_value in sort(collect(keys(ERA5_PATHS)))],
        [oof_qm(string(scheme)) for scheme in reference.cv_schemes])
    provenance = start_provenance(ROOT, OUTDIR; inputs)
    started = time()

    _, ids, product_data = load_global_common_product_data(mger)
    times = product_data["FY4B"].times
    y = Matrix{Float64}(product_data["FY4B"].Y_obs)
    gpm = Matrix{Float64}(product_data[GPM].Y_sat)
    station_meta = load_station_meta(mger.station_meta_path; station_id_col=mger.station_id_col,
        lon_col=mger.lon_col, lat_col=mger.lat_col)
    lonlat = build_X_lonlat(station_meta, ids)
    calibration = _prepare_fy4b_calibration!(reference, ids, product_data)
    id_map = Dict(id => index for (index, id) in enumerate(ids))
    println("Grid: $(length(ids)) stations x $(length(times)) hours")

    labels, rh_all, t2m_all = load_humidity(ids)
    alignment, columns = check_alignment(y, times, labels, rh_all, t2m_all)
    rh = rh_all[:, columns]
    write_csv_atomic(joinpath(OUTDIR, "alignment.csv"), DataFrame(alignment))

    day = Date.(times .- Hour(1))
    labels_regime, regime_index = product_regimes(HOUR_META, times)
    regime = [REGIMES[labels_regime[index]] for index in regime_index]
    season = season_of.(day)
    gross_times = Set(DateTime(first(string(t), 19)) for t in
        filter(row -> row.product == QM, CSV.read(GROSS_FILE, DataFrame)).time)
    gross_hour = [t in gross_times for t in times]
    println("  $(count(gross_hour)) of the budget's $(length(gross_times)) gross FY4B_QM hours are on the grid")

    cut_rows, rh_table, metric_rows, paired, strata_rows_all = (NamedTuple[] for _ in 1:5)
    gross_table = NamedTuple[]
    for scheme_symbol in reference.cv_schemes
        scheme = string(scheme_symbol)
        folds = benchmark_folds(scheme_symbol, ids, lonlat; k=reference.k,
            seed=first(benchmark_seeds(reference)), center_init=reference.fold_center_init,
            rotation=0)
        series = Dict(name => fill(NaN, size(y)) for name in SERIES)
        series[GPM] = gpm
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
            usable = isfinite.(y) .& isfinite.(qm) .& isfinite.(gpm) .& isfinite.(rh)
            train_cells = [CartesianIndex(i, t) for i in train, t in axes(y, 2) if usable[i, t]]
            wet_values = [qm[index] for index in train_cells if qm[index] >= WET]
            cut, csi_before, csi_after, removed = choose_cut(qm, y, rh, train_cells)
            cut8, _, csi8_after, removed8 = choose_cut(qm, y, rh, train_cells; spared=HEAVY)
            gpm_cut, gpm_before, gpm_after, gpm_removed = choose_cut(gpm, y, rh, train_cells)
            amount_cut, amount_removed = matched_cut(wet_values, removed)
            amount_cut8, amount_removed8 = matched_cut(wet_values, removed8)
            series[QM][val, :] = qm[val, :]
            for i in val, t in axes(y, 2)
                x = qm[i, t]
                series[GATE][i, t] = rh[i, t] < cut ? 0.0 : x
                series[GATE8][i, t] = rh[i, t] < cut8 && x < HEAVY ? 0.0 : x
                series[CONTROL][i, t] = x <= amount_cut ? 0.0 : x
                series[CONTROL8][i, t] = x <= amount_cut8 ? 0.0 : x
                series[GPM_GATE][i, t] = rh[i, t] < gpm_cut ? 0.0 : gpm[i, t]
            end
            push!(cut_rows, (; scheme, fold, n_train_station=length(train), train_wet=length(wet_values),
                csi_before, cut, csi_after, removed, amount_cut, amount_removed,
                cut8, csi8_after, removed8, amount_cut8, amount_removed8,
                gpm_before, gpm_cut, gpm_after, gpm_removed))
        end
        check_reproduction(oof_qm(scheme), series[QM], times, ids)

        mask = isfinite.(y) .& isfinite.(rh)
        for name in SERIES
            mask .&= isfinite.(series[name])
        end
        append!(rh_table, rh_rows(scheme, QM, series[QM], y, rh, mask))
        append!(rh_table, rh_rows(scheme, GPM, gpm, y, rh, mask))
        append!(gross_table, gross_rows(scheme, series, y, mask, gross_hour))
        one_group = ones(Int, size(y))
        for name in SERIES
            total = vec(group_sums(series[name], y, mask, one_group, 1, THRESHOLDS))
            push!(metric_rows, merge((; scheme, series=name), sum_stats(total, THRESHOLDS)))
        end
        day_group, n_day = day_groups(times .- Hour(1), size(y, 1))
        sums(x, cells) = group_sums(x, y, cells, day_group, n_day, THRESHOLDS)
        compare(candidate, baseline, cells) = compare_days(sums(series[candidate], cells),
            sums(series[baseline], cells), THRESHOLDS; ci_stats=CI_STATS, reps=BOOTSTRAP_REPS,
            seed=BOOTSTRAP_SEED)
        for (candidate, baseline) in PAIRS
            push!(paired, merge((; scheme, candidate, baseline, n_day), compare(candidate, baseline, mask)))
        end
        for (unit, label, hours) in vcat([("season", s, season .== s) for s in SEASONS],
                [("regime", r, regime .== r) for r in ("R1", "R2", "R3")])
            cells = mask .& permutedims(hours)
            for (candidate, baseline) in STRATA_PAIRS
                push!(strata_rows_all, merge((; scheme, unit, label, candidate, baseline),
                    compare(candidate, baseline, cells)))
            end
        end
    end

    # The rule, on the balanced spatial scheme.
    judged = "balanced_spatial"
    pick(c, b) = only(filter(r -> r.scheme == judged && r.candidate == c && r.baseline == b, paired))
    judged_cuts = filter(r -> r.scheme == judged, cut_rows)
    empty_gate = count(r -> r.cut == -Inf, judged_cuts)
    empty_gate8 = count(r -> r.cut8 == -Inf, judged_cuts)
    rule = vcat(rule_rows(GATE, pick(GATE, CONTROL), pick(GATE, QM), empty_gate, reference.k),
        rule_rows(GATE8, pick(GATE8, CONTROL8), pick(GATE8, QM), empty_gate8, reference.k))

    write_csv_atomic(joinpath(OUTDIR, "cuts.csv"), DataFrame(cut_rows))
    write_csv_atomic(joinpath(OUTDIR, "rh_table.csv"), DataFrame(rh_table))
    write_csv_atomic(joinpath(OUTDIR, "metrics.csv"), DataFrame(metric_rows))
    write_csv_atomic(joinpath(OUTDIR, "paired.csv"), DataFrame(paired))
    write_csv_atomic(joinpath(OUTDIR, "strata.csv"), DataFrame(strata_rows_all))
    write_csv_atomic(joinpath(OUTDIR, "gross_effect.csv"), DataFrame(gross_table))
    write_csv_atomic(joinpath(OUTDIR, "rule.csv"), DataFrame(rule))
    CSV.write(joinpath(OUTDIR, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println("\nCuts (training stations):")
    for row in cut_rows
        @printf("  %-17s fold %d  RH<%-4s CSI %.4f -> %.4f  removed %6d / %6d  amount %.3f | RH8<%-4s CSI %.4f removed %6d | GPM RH<%-4s CSI %.4f -> %.4f\n",
            row.scheme, row.fold, row.cut, row.csi_before, row.csi_after, row.removed, row.train_wet,
            row.amount_cut, row.cut8, row.csi8_after, row.removed8, row.gpm_cut, row.gpm_before, row.gpm_after)
    end
    println("\nHeld-out by RH ($judged):")
    for row in filter(r -> r.scheme == judged, rh_table)
        @printf("  %-10s RH %2d-%3d  cells %7d  gauge wet %.3f | when wet %.3f, when dry %.3f | false alarms %.3f, hits %.3f\n",
            row.product, row.rh_lo, row.rh_hi, row.n_cells, row.gauge_wet, row.gauge_wet_when_wet,
            row.gauge_wet_when_dry, row.share_of_false_alarms, row.share_of_hits)
    end
    println("\nHeld-out scores:")
    for row in metric_rows
        @printf("  %-17s %-13s n=%7d  r %.3f  RMSE %.3f  bias %+.4f  CSI %.3f/%.3f/%.3f  POD_8 %.3f  FAR_0.1 %.3f  FBI_0.1 %.2f\n",
            row.scheme, row.series, row.n, row.r, row.RMSE, row.bias, _stat(row, "CSI_0.1"),
            row.CSI_2, row.CSI_8, row.POD_8, _stat(row, "FAR_0.1"), _stat(row, "FBI_0.1"))
    end
    println("\nPaired (day-block bootstrap, 95%):")
    for row in paired
        @printf("  %-17s %-13s vs %-13s CSI_0.1 %+.4f (%+.4f..%+.4f)  POD_8 %+.4f (%+.4f..%+.4f)  RMSE %+.2f%% (%+.2f..%+.2f)  r %+.4f (%+.4f..%+.4f)\n",
            row.scheme, row.candidate, row.baseline, _stat(row, "delta_CSI_0.1"),
            _stat(row, "delta_CSI_0.1_lo"), _stat(row, "delta_CSI_0.1_hi"), row.delta_POD_8,
            row.delta_POD_8_lo, row.delta_POD_8_hi, 100 * row.rel_RMSE, 100 * row.rel_RMSE_lo,
            100 * row.rel_RMSE_hi, row.delta_r, row.delta_r_lo, row.delta_r_hi)
    end
    println("\nBy season and regime ($judged), CSI_0.1:")
    for row in filter(r -> r.scheme == judged, strata_rows_all)
        @printf("  %-6s %-4s %-13s vs %-13s %+.4f (%+.4f..%+.4f)  RMSE %+.2f%%\n", row.unit, row.label,
            row.candidate, row.baseline, _stat(row, "delta_CSI_0.1"), _stat(row, "delta_CSI_0.1_lo"),
            _stat(row, "delta_CSI_0.1_hi"), 100 * row.rel_RMSE)
    end
    println("\nGross FY4B_QM hours against the rest ($judged), false-rain squared error:")
    for row in filter(r -> r.scheme == judged, gross_table)
        @printf("  %-13s %-5s %.0f (%.1f%% of FY4B_QM's)  false alarms %6d  hits %6d\n", row.series, row.hours,
            row.false_rain_sse, 100 * row.share_of_qm_false_rain_sse, row.false_alarms, row.hits)
    end
    println("\nRule ($judged):")
    for row in rule
        @printf("  %-4s %-12s %s: %+.4f (%+.4f..%+.4f)\n", row.pass ? "pass" : "FAIL", row.arm, row.check,
            row.value, row.lo, row.hi)
    end
    passed(arm) = all(row.pass for row in rule if row.arm == arm)
    verdict = passed(GATE) || passed(GATE8) ?
        "GO: $(join(filter(passed, [GATE, GATE8]), " and ")) meet(s) the rule; only a benchmark run shows what survives MGWR and the blend" :
        "STOP: neither gate meets the rule"
    println("Verdict: $verdict")
    println("Wrote $(relpath(OUTDIR, ROOT)) in $(round((time() - started) / 60; digits=1)) min")
end

if abspath(PROGRAM_FILE) == @__FILE__
    screen_main()
end

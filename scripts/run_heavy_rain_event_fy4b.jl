#!/usr/bin/env julia

# FY4B on the RQ1 heavy-rain events: plain and corrected FY4B against GPM and GSMaP, each event
# scored over the hours FY4B shares with the gauges and the other products.
#
#   julia --project=. scripts/run_heavy_rain_event_fy4b.jl
#   julia --project=. scripts/run_heavy_rain_event_fy4b.jl \
#     --fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv --from 2022-06
#
# With RQ1's `--fy4b-file` (and `--from`), the events are those RQ1 itself selected with FY4B among
# its products (`run_heavy_rain_event_evaluation.jl` with the same flags), FY4B_QM comes from the run
# without `--fy4b-qm-fusion` (plain FY4B in the fusion), and the results go to RQ1's folder plus
# `_fy4bqm`. Each event's window over the four products must then be RQ1's own, and the FY4B, GPM and
# GSMaP rows must reproduce RQ1's `common_hours` rows; both are checked before anything is written.
#
# The events are the paper's six, read from the RQ1 run (`run_heavy_rain_event_evaluation.jl
# --gpm-gsmap`), so FY4B never chooses them and RQ1's GPM-vs-GSMaP results stay as they are. For
# each event:
# - the window is `event_hours` over the gauges and all four products: the hours of the 08-08 day
#   where each has at least 90% of the gauges;
# - the event is kept only when that window passes RQ1's own rule, the `eligible` column of
#   `select_representative_events` (at least 20 hours, at most 20% of the day's gauge rain left out,
#   over the network and at the wettest gauge);
# - a kept event is scored over that window (`fy4b_hours`, all four products), and GPM and GSMaP
#   again over the full day (`full_day`). The full-day rows must reproduce RQ1's own, which is
#   checked before anything is written.
#
# Products:
# - FY4B: the rebuilt 3x3 series of `prepare_fy4b_rebuild.jl`, the FY4B of §2.5;
# - FY4B_QM: corrected FY4B, pixel-shifted and quantile-mapped per fold on training gauges, taken
#   from the FY4B benchmark run's out-of-fold `raw` field. Each gauge's value comes from the
#   balanced spatial fold that held it out, so no gauge is scored against a correction fitted on it;
# - GPM_uncal and GSMaP, as in RQ1.
#
# Writes output/heavy_rain_events_fy4bw3h0_gpmuncal/: event_windows.csv, event_metrics.csv,
# event_metrics_summary.csv, station_event_totals.csv, provenance.csv.

# RQ1's script, for its inputs, analysis window and `summarize_metrics`; its `main` does not run.
include(joinpath(@__DIR__, "run_heavy_rain_event_evaluation.jl"))

# RQ1's events: the GPM/GSMaP run's, or with `--fy4b-file` the run that selected them with FY4B.
const WITH_RQ1_FY4B = FY4B_FILE !== nothing
const RQ1_DIR = WITH_RQ1_FY4B ? OUTDIR : joinpath(ROOT, "output", "heavy_rain_events_gpm_gsmap_gpmuncal")
# The FY4B benchmark run: with `--fy4b-file` the one with plain FY4B in the fusion, otherwise the
# `--fy4b-qm-fusion` run of the earlier `fy4b_benchmark` step.
const FY4B_BENCHMARK_DIR = joinpath(ROOT, "output",
    "interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only_satwetblend_blendagrenv_" *
    "fusedanchor_fy4bw3h0_fy4bqm_fy4bqmadapt_" * (WITH_RQ1_FY4B ? "" : "qmfusion_") * "gpmuncal")
const EVENT_PRODUCT_FILES = [
    "FY4B" => something(FY4B_FILE, joinpath(ROOT, "output", "fy4b_rebuild_2022_2024", "fy4b_w3_h0.csv")),
    "FY4B_QM" => joinpath(FY4B_BENCHMARK_DIR, "balanced_spatial", "fy4b_qm", "oof_raw.csv"),
    [product => joinpath(STUDY_DATA, file) for (product, file) in ALL_PRODUCT_FILES
        if product in FULL_DAY_PRODUCTS]...,
]
const OBS_FILE = joinpath(STUDY_DATA, "hubei_obs_hourly_2022_2025_JunSep.csv")
const STATION_META_FILE = joinpath(STUDY_DATA, "station_meta.csv")
const FY4B_EVENT_DIR = WITH_RQ1_FY4B ? OUTDIR * "_fy4bqm" :
    joinpath(ROOT, "output", "heavy_rain_events_fy4bw3h0_gpmuncal")
# RQ1's metrics a full-day row must reproduce.
const CHECKED_METRICS = [:n, :r, :Bias, :RMSE, :CSI_50]

"""The gauges and the four products on RQ1's time axis, as `load_inputs` reads RQ1's."""
function load_event_inputs()
    obs_times, obs_ids, Y_obs = read_hourly_wide(OBS_FILE)
    obs_times, Y_obs = subset_time_window(obs_times, Y_obs, ANALYSIS_START, ANALYSIS_END)
    products = Dict{String,Matrix{Float64}}()
    for (product, path) in EVENT_PRODUCT_FILES
        times, ids, Y = read_hourly_wide(path)
        times, Y = subset_time_window(times, Y, ANALYSIS_START, ANALYSIS_END)
        products[product] = align_to_reference(obs_times, obs_ids, times, ids, Y)
    end
    lonlat = build_X_lonlat(load_station_meta(STATION_META_FILE), obs_ids)
    return obs_times, obs_ids, Y_obs, products, lonlat
end

"""
Fail unless each event's window is RQ1's own (`selected_events.csv`'s `common_hours`, when RQ1
selected with FY4B) and the FY4B, GPM and GSMaP rows over it equal RQ1's `common_hours` rows.
"""
function check_common_window(windows::DataFrame, metrics::DataFrame, rq1::DataFrame)
    for event in eachrow(windows)
        event.common_hours == event.rq1_common_hours ||
            error("$(event.day): $(event.common_hours) hours with FY4B_QM, $(event.rq1_common_hours) in RQ1")
    end
    for row in eachrow(metrics[(metrics.window .== "fy4b_hours") .& (metrics.product .!= "FY4B_QM"), :])
        match = rq1[(rq1.window .== "common_hours") .& (rq1.event_day .== row.event_day) .&
            (rq1.product .== row.product), :]
        nrow(match) == 1 || error("RQ1 has $(nrow(match)) common-hours rows for $(row.event_day) $(row.product)")
        for column in CHECKED_METRICS
            a, b = row[column], match[1, column]
            isequal(a, b) || isapprox(a, b; rtol=1e-9) ||
                error("$(row.event_day) $(row.product) $column: $a here, $b in RQ1 (FY4B_QM drops gauges?)")
        end
    end
    return nothing
end

"""Fail unless every full-day row equals RQ1's row for the same event and product."""
function check_full_day(metrics::DataFrame, rq1::DataFrame)
    for row in eachrow(metrics[metrics.window .== "full_day", :])
        match = rq1[(rq1.window .== "full_day") .& (rq1.event_day .== row.event_day) .&
            (rq1.product .== row.product), :]
        nrow(match) == 1 || error("RQ1 has $(nrow(match)) full-day rows for $(row.event_day) $(row.product)")
        for column in CHECKED_METRICS
            a, b = row[column], match[1, column]
            isequal(a, b) || isapprox(a, b; rtol=1e-9) ||
                error("$(row.event_day) $(row.product) $column: $a here, $b in RQ1")
        end
    end
    return nothing
end

function fy4b_events_main()
    isempty(ARGS) || WITH_RQ1_FY4B || throw(ArgumentError(
        "run_heavy_rain_event_fy4b.jl takes no arguments, or RQ1's --fy4b-file <path> [--from YYYY-MM]"))
    events_file = joinpath(RQ1_DIR, "selected_events.csv")
    rq1_metrics_file = joinpath(RQ1_DIR, "event_metrics.csv")
    provenance = start_provenance(ROOT, FY4B_EVENT_DIR; inputs=[OBS_FILE, last.(EVENT_PRODUCT_FILES)...,
        STATION_META_FILE, events_file, rq1_metrics_file])
    events = CSV.read(events_file, DataFrame)
    events.day = Date.(string.(events.day))
    # RQ1's own window, which the four-product window must equal when RQ1 selected with FY4B.
    rq1_common_hours = Dict(zip(events.day, events.common_hours))
    rq1 = CSV.read(rq1_metrics_file, DataFrame)
    rq1.event_day = Date.(string.(rq1.event_day))

    println("Reading study-area hourly tables ...")
    obs_times, obs_ids, Y_obs, products, lonlat = load_event_inputs()
    names = first.(EVENT_PRODUCT_FILES)
    println("  $(length(obs_ids)) gauges x $(length(obs_times)) hours, products: $(join(names, ", "))")

    coverage = Dict(day => event_hours(obs_times, Y_obs, products, day) for day in events.day)
    windows = copy(events)
    for product in names
        windows[!, Symbol(lowercase(product), "_hours")] = [coverage[day].product_hours[product] for day in windows.day]
    end
    windows.common_hours = [coverage[day].common_hours for day in windows.day]
    windows.excluded_rain_share = [coverage[day].excluded_rain_share for day in windows.day]
    windows.peak_excluded_rain_share = [coverage[day].peak_excluded_rain_share for day in windows.day]
    # RQ1's eligibility rule on FY4B's window. Nothing is picked (n = 0): the events stay RQ1's.
    windows.eligible = select_representative_events(windows; n_widespread=0, n_localized=0).eligible
    windows.rq1_common_hours = [rq1_common_hours[day] for day in windows.day]
    for event in eachrow(windows)
        println("  $(event.day)  $(rpad(event.event_type, 10))  common hours=$(event.common_hours)  " *
            "excluded rain=$(round(100 * event.excluded_rain_share; digits=1))% " *
            "(peak gauge $(round(100 * event.peak_excluded_rain_share; digits=1))%)  " *
            (event.eligible ? "kept" : "left out"))
    end

    days, obs_daily = daily_totals(obs_times, Y_obs)
    day_column = Dict(day => k for (k, day) in enumerate(days))
    metric_rows = NamedTuple[]
    station_rows = NamedTuple[]
    for event in eachrow(windows[windows.eligible, :])
        hours = coverage[event.day]
        for (window, hour_idx, window_products) in (
            ("fy4b_hours", hours.hour_idx, names),
            ("full_day", hours.day_hour_idx, FULL_DAY_PRODUCTS),
        )
            totals = event_station_totals(Y_obs, Dict(p => products[p] for p in window_products), hour_idx)
            for product in window_products
                scores = spatial_metrics(totals.obs, totals.sat[product], lonlat)
                push!(metric_rows, merge(
                    (; event_day=event.day, event_type=event.event_type, product, window, hours=length(hour_idx)),
                    scores,
                ))
                for i in findall(totals.keep)
                    push!(station_rows, (;
                        event_day=event.day, event_type=event.event_type, window, product,
                        station_id=obs_ids[i], lon=lonlat[i, 1], lat=lonlat[i, 2],
                        obs_24h_mm=obs_daily[i, day_column[event.day]],
                        obs_mm=totals.obs[i], sat_mm=totals.sat[product][i],
                        diff_mm=totals.sat[product][i] - totals.obs[i],
                    ))
                end
            end
        end
    end
    metrics = DataFrame(metric_rows)
    stations = DataFrame(station_rows)
    check_full_day(metrics, rq1)
    WITH_RQ1_FY4B && check_common_window(windows, metrics, rq1)
    lonlat_by_id = Dict(id => (lonlat[i, 1], lonlat[i, 2]) for (i, id) in enumerate(obs_ids))
    summary = summarize_metrics(metrics, stations, lonlat_by_id)

    window_columns = [:day, :event_type, Symbol.(lowercase.(names), "_hours")..., :common_hours,
        :excluded_rain_share, :peak_excluded_rain_share, :eligible]
    write_csv_atomic(joinpath(FY4B_EVENT_DIR, "event_windows.csv"),
        rename(select(windows, window_columns), :day => :event_day))
    write_csv_atomic(joinpath(FY4B_EVENT_DIR, "event_metrics.csv"), metrics)
    write_csv_atomic(joinpath(FY4B_EVENT_DIR, "event_metrics_summary.csv"), summary)
    write_csv_atomic(joinpath(FY4B_EVENT_DIR, "station_event_totals.csv"), stations)

    println("\nSpatial agreement over FY4B's window (gauge-pixel pairs):")
    show(select(metrics[metrics.window .== "fy4b_hours", :],
        :event_day, :event_type, :product, :hours, :n, :r, :RMSE, :Bias, :CSI_50, :centroid_shift_km);
        allrows=true, allcols=true)
    println("\n\nMean over the kept events:")
    show(select(summary[(summary.statistic .== "mean_over_events"), :],
        :product, :window, :n_events, :r, :RMSE, :Bias, :CSI_50); allrows=true, allcols=true)
    write_csv_atomic(joinpath(FY4B_EVENT_DIR, "provenance.csv"), finish_provenance!(provenance, ROOT))
    println("\n\nWrote $(relpath(FY4B_EVENT_DIR, ROOT))")
    return FY4B_EVENT_DIR
end

if abspath(PROGRAM_FILE) == @__FILE__
    fy4b_events_main()
end

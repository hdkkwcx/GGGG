#!/usr/bin/env julia

# Trace FY4B from the raw 15-minute files to the numbers the paper cites, stage by stage, and check
# the stages against each other. Reads existing outputs only: nothing is refitted, and the raw drive
# is not needed (the rebuild's fy4b_hourly_qc.csv lists every file used).
#
#   julia --project=. scripts/trace_fy4b_process.jl          # -> output/fy4b_process_trace/*.csv
#   $MIXEDGWR_PYTHON scripts/plot_fy4b_process.py            # -> output/fy4b_process_trace/fy4b_process.html
#
# It follows the one path the paper uses: the rebuilt 3x3 series (fy4b_w3_h0.csv) as the third
# product of RQ1-RQ3, and corrected FY4B (FY4B_QM) from the `benchmark` step of paper/sources.toml.
# One CSV per stage, and checks.csv with every cross-check as pass or fail; the script exits with an
# error when a check fails, after writing all tables. Nothing in the paper reads these tables, so
# this is not a step of sources.toml.

const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "src", "load_modules.jl"))
load_standalone_modules("RunProvenance", "TableIO", "SatelliteCalibration")

using CSV, DataFrames, Dates, Statistics, TOML
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic
using Main.SatelliteCalibration: product_regimes, season_index, SEASONS

const OUTDIR = joinpath(ROOT, "output", "fy4b_process_trace")
const REBUILD = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024")
const REBUILD_LOG = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024.log")
# scripts/list_fy4b_raw_files.jl's record of the raw drive; optional, since it needs the drive.
const RAW_MANIFEST = joinpath(ROOT, "output", "fy4b_raw_manifest")
const STUDY = joinpath(ROOT, "data", "processed", "study_area")
const SOURCES = joinpath(ROOT, "paper", "sources.toml")
const VARIANT = "fy4b_w3_h0"
const WET = 0.1
const FY4B_CAP = 20.0
# The first and last hour-ending Beijing labels of the 08-08 days 2022-06-01 .. 2024-12-31, the
# record of the paper's matched sample (`--from 2022-06`).
const RECORD = DateTime(2022, 6, 1, 9):Hour(1):DateTime(2025, 1, 1, 8)
# A file's UTC hour start to the hour-ending Beijing label the series carry.
const LABEL_OFFSET = Hour(9)
# The steps of paper/sources.toml that read or make FY4B, in the order the page walks them.
const FY4B_STEPS = ["raw_fy4b", "fy4b_rebuild", "rq1_events", "rq1_cores", "rq2", "benchmark", "claims",
    "fusion", "fy4b_quality_week", "fy4b_time_holdout", "fy4b_clt", "rq1_fy4b", "benchmark_qmfusion",
    "fy4b_fusion_configs"]

parse_label(text) = DateTime(first(string(text), 19), dateformat"yyyy-mm-ddTHH:MM:SS")
met_day(label::DateTime) = Date(label - Hour(9))
relative(path) = replace(relpath(path, ROOT), '\\' => '/')

"""A wide `time, <station>...` table as `(times, ids, values[time, station])`, blanks as NaN."""
function read_wide(path)
    table = CSV.read(path, DataFrame; types=Dict(:time => String))
    ids = [name for name in names(table) if name != "time"]
    values = Matrix{Float64}(undef, nrow(table), length(ids))
    for (j, id) in enumerate(ids)
        values[:, j] = Float64.(coalesce.(table[!, id], NaN))
    end
    return parse_label.(table.time), ids, values
end

"""`values` on the hours `times` and stations `ids`, NaN where the file has no row or column."""
function align(times, ids, file_times, file_ids, values)
    row_of = Dict(t => i for (i, t) in enumerate(file_times))
    col_of = Dict(id => j for (j, id) in enumerate(file_ids))
    out = fill(NaN, length(times), length(ids))
    for (j, id) in enumerate(ids)
        haskey(col_of, id) || continue
        column = view(values, :, col_of[id])
        for (i, t) in enumerate(times)
            row = get(row_of, t, 0)
            row > 0 && (out[i, j] = column[row])
        end
    end
    return out
end

"""
Why an hour has no value, from its row of the rebuild's QC table. The causes do not overlap: an
hour is given the first that applies, in this order.
"""
function hour_cause(status, hour_start, file_count, first_file, gap)
    status == "complete" && return "complete"
    hour_start < first_file && return "before_first_file"
    gap[1] < hour_start < gap[2] && return "relocation_gap"
    hour(hour_start) == 14 && return "quarter_1415_utc"
    file_count == 0 && return "no_files"
    return "quarters_missing"
end

const CAUSE_TEXT = Dict(
    "complete" => "All four 15-minute files exist: the hour is kept",
    "before_first_file" => "Before the first file (the record starts on 1 June 2022)",
    "relocation_gap" => "Between the last hour at 133°E and the first at 105°E",
    "quarter_1415_utc" => "The hour 14-15 UTC (22-23 Beijing time): its 14:15 file is never produced",
    "no_files" => "No file at all in the hour",
    "quarters_missing" => "One to three of the four files missing or unreadable",
)

function qc_tables(qc, meta_times, meta)
    hour_start = parse_label.(qc.hour_start)
    label = hour_start .+ LABEL_OFFSET
    first_file = minimum(hour_start[qc.file_count .> 0])
    lon_of = Dict(parse_label(row.time) => row.sat_lon for row in eachrow(meta))
    last_133 = maximum(t for t in meta_times if lon_of[t] > 120) - LABEL_OFFSET
    first_105 = minimum(t for t in meta_times if lon_of[t] < 120) - LABEL_OFFSET
    cause = [hour_cause(qc.status[i], hour_start[i], qc.file_count[i], first_file, (last_133, first_105))
             for i in 1:nrow(qc)]

    order = ["before_first_file", "relocation_gap", "quarter_1415_utc", "no_files", "quarters_missing", "complete"]
    funnel = DataFrame(cause=String[], what=String[], hours=Int[], files_in_them=Int[],
        first_hour_utc=String[], last_hour_utc=String[])
    for key in order
        rows = findall(==(key), cause)
        push!(funnel, (key, CAUSE_TEXT[key], length(rows), sum(qc.file_count[rows]; init=0),
            isempty(rows) ? "" : string(minimum(hour_start[rows])),
            isempty(rows) ? "" : string(maximum(hour_start[rows]))))
    end

    months = DataFrame(month=String[], hours_possible=Int[], files_expected=Int[], files_found=Int[],
        hours_complete=Int[], unreadable_files=Int[])
    for key in sort(unique(Dates.format.(hour_start, "yyyy-mm")))
        rows = findall(==(key), Dates.format.(hour_start, "yyyy-mm"))
        push!(months, (key, length(rows), 4 * length(rows), sum(qc.file_count[rows]),
            count(==("complete"), qc.status[rows]), sum(qc.unreadable_count[rows])))
    end

    in_record = [c in ("complete", "quarter_1415_utc", "no_files", "quarters_missing") for c in cause]
    by_hour = DataFrame(utc_hour=Int[], beijing_hour_ending=Int[], hours_in_record=Int[], hours_complete=Int[])
    for h in 0:23
        rows = findall(i -> in_record[i] && hour(hour_start[i]) == h, 1:nrow(qc))
        push!(by_hour, (h, mod(h + 9, 24), length(rows), count(==("complete"), qc.status[rows])))
    end

    partial = findall(==("quarters_missing"), cause)
    pattern = string.(coalesce.(qc.missing_minutes[partial], ""))
    quarters = sort(combine(groupby(DataFrame(missing_minutes=pattern), :missing_minutes), nrow => :hours),
        :hours; rev=true)
    return (; hour_start, label, cause, funnel, months, by_hour, quarters, first_file, last_133, first_105)
end

"""Scores of FY4B against the gauges on the cells of `mask`."""
function cell_scores(obs, sat, cap, mask)
    o, s = obs[mask], sat[mask]
    n = length(o)
    n == 0 && return nothing
    wet_o, wet_s = o .>= WET, s .>= WET
    hits, false_alarms, misses = count(wet_o .& wet_s), count(.!wet_o .& wet_s), count(wet_o .& .!wet_s)
    capped = cap[mask]
    return (; cells=n, gauge_wet_share=mean(wet_o), fy4b_wet_share=mean(wet_s),
        gauge_mean=mean(o), fy4b_mean=mean(s), relative_bias_pct=100 * (sum(s) - sum(o)) / sum(o),
        r=cor(o, s), rmse=sqrt(mean(abs2.(s .- o))),
        pod=hits / (hits + misses), far=false_alarms / (hits + false_alarms),
        csi=hits / (hits + misses + false_alarms),
        gauge_max=maximum(o), fy4b_max=maximum(s),
        gauge_cells_ge20=count(>=(FY4B_CAP), o), fy4b_cells_ge20=count(>=(FY4B_CAP), s),
        fy4b_wet_cells_touching_cap=count((capped .> 0) .& wet_s))
end

function provenance_value(table, key)
    rows = table[table.key .== key, :value]
    return isempty(rows) ? "" : string(only(rows))
end

function step_rows(provenance_rows)
    steps = Dict(step["id"] => step for step in TOML.parsefile(SOURCES)["step"])
    is_ancestor(commit) = success(pipeline(`git -C $ROOT merge-base --is-ancestor $commit HEAD`;
        stdout=devnull, stderr=devnull))
    rows = NamedTuple[]
    for id in FY4B_STEPS
        step = steps[id]
        marker = get(step, "marker", "")
        commit = dirty = dirty_end = fy4b_sha = ""
        ancestor = ""
        if !isempty(marker) && isfile(joinpath(ROOT, marker))
            table = CSV.read(joinpath(ROOT, marker), DataFrame; types=String)
            commit = provenance_value(table, "git_commit")
            dirty = provenance_value(table, "git_dirty")
            dirty_end = provenance_value(table, "git_dirty_at_end")
            fy4b_sha = provenance_value(table, "input_sha256:output/fy4b_rebuild_2022_2024/$VARIANT.csv")
            ancestor = isempty(commit) ? "" : string(is_ancestor(commit))
        end
        push!(rows, (; step=id, kind=step["kind"], command=join(get(step, "command", String[]), " "),
            needs=join(get(step, "needs", String[]), " "), marker,
            marker_exists=isempty(marker) ? "" : string(isfile(joinpath(ROOT, marker))),
            git_commit=commit, git_dirty=dirty, git_dirty_at_end=dirty_end, commit_in_history=ancestor,
            fy4b_file_sha256=fy4b_sha, note=replace(get(step, "note", ""), r"\s+" => " ")))
    end
    return DataFrame(rows)
end

function main()
    run_dir = dirname(joinpath(ROOT, TOML.parsefile(SOURCES)["step"][findfirst(
        step -> step["id"] == "benchmark", TOML.parsefile(SOURCES)["step"])]["marker"]))
    rq2_coverage = joinpath(ROOT, "output", "satellite_temporal_evaluation_fy4bw3h0_gpmuncal_from202206_block7d",
        "sample_coverage.csv")
    fusion_sample = joinpath(ROOT, "output",
        "gauge_satellite_fusion_fy4bw3h0_fy4bqm_fy4bqmadapt_gpmuncal_matched_block7d", "matched_sample.csv")
    files = (;
        qc=joinpath(REBUILD, "fy4b_hourly_qc.csv"), meta=joinpath(REBUILD, "fy4b_hour_meta.csv"),
        series=joinpath(REBUILD, "$VARIANT.csv"), capfrac=joinpath(REBUILD, "$(VARIANT)_capfrac.csv"),
        centre=joinpath(REBUILD, "fy4b_w1_h0.csv"), rebuild_provenance=joinpath(REBUILD, "provenance.csv"),
        stations=joinpath(STUDY, "station_meta.csv"),
        obs=joinpath(STUDY, "hubei_obs_hourly_2022_2025_JunSep.csv"),
        gpm=joinpath(STUDY, "hubei_gpm_uncal_hourly_2022_2024_full_aligned.csv"),
        gsmap=joinpath(STUDY, "hubei_gsmap_hourly_2022_2024_full_aligned.csv"),
        scope=joinpath(run_dir, "benchmark_scope.csv"),
        shifts=joinpath(run_dir, "fy4b_qm_shift_selection.csv"), maps=joinpath(run_dir, "fy4b_qm_maps.csv"),
        strengths=joinpath(run_dir, "fy4b_qm_strengths.csv"),
        fused=joinpath(run_dir, "fused_anchor_selection.csv"), blend=joinpath(run_dir, "blend_selection.csv"),
        rq2_coverage, fusion_sample, sources=SOURCES, tries=joinpath(ROOT, "paper", "fy4b_tries.csv"),
    )
    provenance = start_provenance(ROOT, OUTDIR; inputs=collect(values(files)))
    checks = DataFrame(check=String[], expected=String[], found=String[], pass=Bool[])
    check!(name, expected, found) = push!(checks, (name, string(expected), string(found), expected == found))
    write(name, table) = write_csv_atomic(joinpath(OUTDIR, name), table)

    # --- Stage 1-2: raw files and the strict hour -------------------------------------------------
    qc = CSV.read(files.qc, DataFrame; types=Dict(:hour_start => String, :status => String,
        :missing_minutes => String, :present_minutes => String))
    meta = CSV.read(files.meta, DataFrame; types=Dict(:time => String, :software => String))
    meta_times = parse_label.(meta.time)
    q = qc_tables(qc, meta_times, meta)
    log_counts = Dict{String,String}()
    if isfile(REBUILD_LOG)
        text = read(REBUILD_LOG, String)
        for (key, pattern) in ("duplicate_files_dropped" => r"duplicate_file_count = (\d+)",
                               "files_validated" => r"Validated \d+ / (\d+) NetCDF files",
                               "unreadable_segments_warned" => r"Unreadable FY4B segment")
            found = key == "unreadable_segments_warned" ? string(length(collect(eachmatch(pattern, text)))) :
                (m = match(pattern, text); m === nothing ? "" : m.captures[1])
            log_counts[key] = found
        end
    end
    complete = q.cause .== "complete"
    raw = DataFrame(item=String[], value=String[], source=String[])
    push!(raw, ("hours of 2022-2024 (UTC)", string(nrow(qc)), relative(files.qc)))
    push!(raw, ("15-minute files expected", string(4 * nrow(qc)), relative(files.qc)))
    push!(raw, ("15-minute files used, after duplicates are dropped", string(sum(qc.file_count)), relative(files.qc)))
    push!(raw, ("unreadable files, treated as missing", string(sum(qc.unreadable_count)), relative(files.qc)))
    push!(raw, ("first hour with a file (UTC)", string(q.first_file), relative(files.qc)))
    push!(raw, ("last complete hour at 133°E (UTC start)", string(q.last_133), relative(files.meta)))
    push!(raw, ("first complete hour at 105°E (UTC start)", string(q.first_105), relative(files.meta)))
    push!(raw, ("hours with all four files", string(count(complete)), relative(files.qc)))
    for key in sort(collect(keys(log_counts)))
        push!(raw, (replace(key, "_" => " ") * " (run log)", log_counts[key], relative(REBUILD_LOG)))
    end
    manifest = joinpath(RAW_MANIFEST, "summary.csv")
    if isfile(manifest)
        for row in eachrow(CSV.read(manifest, DataFrame; types=String))
            push!(raw, (row.item * " (raw drive)", row.value, relative(manifest)))
        end
        listed = CSV.read(joinpath(RAW_MANIFEST, "provenance.csv"), DataFrame; types=String)
        push!(raw, ("raw drive listed at commit", provenance_value(listed, "git_commit"),
            relative(joinpath(RAW_MANIFEST, "provenance.csv"))))
        for row in eachrow(CSV.read(joinpath(RAW_MANIFEST, "checks.csv"), DataFrame; types=String))
            push!(checks, (row.check, row.expected, row.found, row.pass == "true"))
        end
    else
        push!(raw, ("raw drive", "not listed: run scripts/list_fy4b_raw_files.jl with the drive mounted", ""))
    end
    write("raw_files.csv", raw)
    write("hour_funnel.csv", q.funnel)
    write("hours_by_month.csv", q.months)
    write("hours_by_utc_hour.csv", q.by_hour)
    write("missing_quarters.csv", q.quarters)

    check!("QC table covers every hour of 2022-2024", 26_304, nrow(qc))
    check!("every complete hour has exactly the quarters 00;15;30;45",
        count(complete), count(complete .& (coalesce.(qc.present_minutes, "") .== "00;15;30;45") .& (qc.file_count .== 4)))
    check!("no complete hour starts at 14 UTC (the hour ending 23:00 Beijing time)", 0,
        count(complete .& (hour.(q.hour_start) .== 14)))
    check!("no complete hour before 1 June 2022 00 UTC", 0, count(complete .& (q.hour_start .< DateTime(2022, 6, 1))))
    check!("no complete hour in February 2024", 0,
        count(complete .& (Dates.format.(q.hour_start, "yyyy-mm") .== "2024-02")))
    check!("no file at all in the relocation gap", 0,
        only(q.funnel[q.funnel.cause .== "relocation_gap", :files_in_them]))

    # --- Stage 3-4: the series and its regimes ---------------------------------------------------
    stations = CSV.read(files.stations, DataFrame; types=Dict(:station_id => String))
    ids = stations.station_id
    fy_times, fy_ids, fy_values = read_wide(files.series)
    _, cap_ids, cap_values = read_wide(files.capfrac)
    centre_times, centre_ids, centre_values = read_wide(files.centre)
    regime_labels, regime_of_time = product_regimes(files.meta, fy_times)
    complete_labels = sort(q.label[complete])
    check!("series hours = strict complete hours of the QC table", complete_labels == sort(fy_times), true)
    check!("regime table hours = series hours", sort(meta_times) == sort(fy_times), true)
    check!("series has one column per gauge of the study area", sort(ids) == sort(fy_ids), true)
    check!("every series hour lies in the record June 2022 - December 2024", true,
        all(t -> first(RECORD) <= t <= last(RECORD), fy_times))
    finite_values = filter(isfinite, fy_values)
    check!("no negative value", 0, count(<(0), finite_values))
    check!("no value above the product's 20 mm/h ceiling", 0, count(>(FY4B_CAP + 1e-6), finite_values))

    regimes = DataFrame(regime=String[], sat_lon=String[], software=String[], first_hour=String[],
        last_hour=String[], hours=Int[])
    for g in sortperm([minimum(fy_times[regime_of_time .== k]) for k in eachindex(regime_labels)])
        times = fy_times[regime_of_time .== g]
        lon, software = split(regime_labels[g], "E_")
        push!(regimes, (regime_labels[g], lon, software, string(minimum(times)), string(maximum(times)), length(times)))
    end

    # --- Stage 5: the matched sample ---------------------------------------------------------------
    times = collect(RECORD)
    fy = align(times, ids, fy_times, fy_ids, fy_values)
    cap = align(times, ids, fy_times, cap_ids, cap_values)
    centre = align(times, ids, centre_times, centre_ids, centre_values)
    obs = align(times, ids, read_wide(files.obs)...)
    gpm = align(times, ids, read_wide(files.gpm)...)
    gsmap = align(times, ids, read_wide(files.gsmap)...)
    has_hour = in(Set(fy_times)).(times)
    days = met_day.(times)
    funnel = DataFrame(stage=String[], what=String[], hours=Int[], days=Int[], gauge_hours=Int[])
    function push_stage!(stage, what, mask)
        rows = vec(any(mask; dims=2))
        push!(funnel, (stage, what, count(rows), length(unique(days[rows])), count(mask)))
    end
    everything = trues(size(fy))
    in_strict = everything .& has_hour
    m_fy = isfinite.(fy)
    m_obs = m_fy .& isfinite.(obs)
    m_gpm = m_obs .& isfinite.(gpm)
    matched = m_gpm .& isfinite.(gsmap)
    push_stage!("record", "Every hour of the 08-08 days 1 June 2022 - 31 December 2024, at every gauge", everything)
    push_stage!("strict_hour", "Hours with all four FY4B files", in_strict)
    push_stage!("fy4b_value", "FY4B has a value at the gauge (5 of the 9 pixels accepted in every quarter)", m_fy)
    push_stage!("gauge_value", "and the gauge has a value", m_obs)
    push_stage!("gpm_value", "and GPM (IMERG Uncal) has a value", m_gpm)
    push_stage!("matched_sample", "and GSMaP has a value: the matched sample of RQ1-RQ3", matched)
    sample = CSV.read(files.fusion_sample, DataFrame; types=Dict(:fold => String))
    scored = only(eachrow(sample[(sample.scheme .== "balanced_spatial") .& (sample.product .== "FY4B") .&
        (sample.fold .== "all"), :]))
    push!(funnel, ("scored", "and every anchor's prediction by every method exists (Table 3, balanced spatial CV)",
        scored.hours_matched, 0, scored.n_matched))
    write("sample_funnel.csv", funnel)
    write("sample_by_fold.csv", sample[in.(sample.product, Ref(("FY4B", "FY4B_QM"))), :])

    coverage = CSV.read(files.rq2_coverage, DataFrame)
    paper = only(eachrow(coverage[(coverage.sample .== "all_products") .& (coverage.season .== "all") .&
        (coverage.region .== "all"), :]))
    last_row = funnel[funnel.stage .== "matched_sample", :]
    check!("matched hours = RQ2's all_products sample", paper.n_hours, only(last_row.hours))
    check!("matched days = RQ2's all_products sample", paper.n_days, only(last_row.days))
    check!("matched gauge-hours = RQ2's all_products sample", paper.n_cells, only(last_row.gauge_hours))
    check!("strict hours inside the record = series hours", length(fy_times), only(funnel[funnel.stage .== "strict_hour", :hours]))

    # Complete hours with no value: every file exists, but the values at the gauges are rejected.
    strict_rows = findall(has_hour)
    none(matrix) = vec(.!any(isfinite.(matrix); dims=2))
    some(matrix) = vec(any(isfinite.(matrix); dims=2) .& .!all(isfinite.(matrix); dims=2))
    no_value, partial_value, no_centre = none(fy), some(fy), none(centre)
    regime_at = Dict(t => regime_labels[g] for (t, g) in zip(fy_times, regime_of_time))
    regime_of_row = [get(regime_at, t, "") for t in times]
    valueless = NamedTuple[]
    function add_valueless!(group, level, rows)
        keep = intersect(strict_rows, findall(rows))
        push!(valueless, (; group, level, complete_hours=length(keep),
            no_value_at_any_gauge=count(no_value[keep]), value_at_some_gauges_only=count(partial_value[keep]),
            centre_pixel_no_value_at_any_gauge=count(no_centre[keep]),
            gauge_hours_without_value=count(.!isfinite.(fy[keep, :]))))
    end
    add_valueless!("all", "all", trues(length(times)))
    for label in regimes.regime
        add_valueless!("regime", label, regime_of_row .== label)
    end
    for key in sort(unique(Dates.format.(times, "yyyy-mm")))
        add_valueless!("month", key, Dates.format.(times, "yyyy-mm") .== key)
    end
    for h in 0:23
        add_valueless!("utc_hour", lpad(h, 2, '0'), hour.(times .- LABEL_OFFSET) .== h)
    end
    write("valueless_hours.csv", DataFrame(valueless))
    check!("complete hours - hours with no value at any gauge = matched hours",
        only(funnel[funnel.stage .== "matched_sample", :hours]), length(strict_rows) - count(no_value[strict_rows]))

    regimes.matched_hours = [count(vec(any(matched; dims=2)) .& in(Set(fy_times[regime_of_time .== findfirst(==(label), regime_labels)])).(times))
                             for label in regimes.regime]
    write("regimes.csv", regimes)

    # --- The series against the gauges, on the matched sample ---------------------------------------
    season_of_row = [SEASONS[season_index(t)] for t in times]
    month_of_row = Dates.format.(times, "yyyy-mm")
    summary = NamedTuple[]
    add!(table, group, level, rows) = begin
        scores = cell_scores(obs, fy, cap, matched .& rows)
        scores === nothing || push!(table, merge((; group, level, hours=count(vec(any(matched .& rows; dims=2)))), scores))
    end
    add!(summary, "all", "all", trues(length(times)))
    for label in regimes.regime
        add!(summary, "regime", label, regime_of_row .== label)
    end
    for season in SEASONS
        add!(summary, "season", season, season_of_row .== season)
    end
    for label in regimes.regime, season in SEASONS
        add!(summary, "regime_season", "$label $season", (regime_of_row .== label) .& (season_of_row .== season))
    end
    write("series_summary.csv", DataFrame(summary))
    monthly = NamedTuple[]
    for key in sort(unique(month_of_row))
        add!(monthly, "month", key, month_of_row .== key)
    end
    write("series_by_month.csv", DataFrame(monthly))
    # The 3x3 mean against the centre pixel, raw, on the matched cells where both have a value.
    both = matched .& isfinite.(centre)
    window = DataFrame(variant=["w3_h0 (3x3 mean, used)", "w1_h0 (centre pixel)"],
        cells=[count(both), count(both)],
        r=[cor(obs[both], fy[both]), cor(obs[both], centre[both])],
        rmse=[sqrt(mean(abs2.(fy[both] .- obs[both]))), sqrt(mean(abs2.(centre[both] .- obs[both])))],
        wet_share=[mean(fy[both] .>= WET), mean(centre[both] .>= WET)],
        mean=[mean(fy[both]), mean(centre[both])])
    write("window_choice.csv", window)

    # --- Stage 6-8: the per-fold correction and FY4B's place in the fusion -------------------------
    shifts = CSV.read(files.shifts, DataFrame)
    maps = CSV.read(files.maps, DataFrame)
    strengths = CSV.read(files.strengths, DataFrame; types=Dict(:cell => String))
    write("shifts.csv", shifts)
    write("maps.csv", maps)
    write("strengths.csv", strengths)
    write("fusion_merge_weights.csv", CSV.read(files.fused, DataFrame))
    blend = CSV.read(files.blend, DataFrame)
    write("fusion_blend_weights.csv", blend[(blend.method .== "blend_agrenv_mgwr") .&
        in.(blend.product, Ref(("FY4B", "FY4B_QM"))), :])
    scope = CSV.read(files.scope, DataFrame; types=String)
    write("benchmark_scope.csv", scope)
    for row in eachrow(regimes)
        selected = unique(shifts[shifts.regime .== row.regime, :n_hours])
        check!("regime $(row.regime): hours = n_hours of the shift selection", [row.hours], selected)
    end
    check!("five folds x three regimes of shifts, balanced spatial", 15,
        count(==("balanced_spatial"), shifts.scheme))
    check!("shift never lowers the training correlation", 0, count(shifts.r_train .< shifts.r_train_unshifted))
    check!("every strength is one of 0, 0.25, 0.5, 0.75, 1", true,
        all(in((0.0, 0.25, 0.5, 0.75, 1.0)), strengths.strength))
    check!("the benchmark ran from a clean tree", ("false", "false"),
        (provenance_value(scope, "git_dirty"), provenance_value(scope, "git_dirty_at_end")))

    # --- Steps and provenance ---------------------------------------------------------------------
    steps = step_rows(provenance)
    write("steps.csv", steps)
    series_sha = only(row.value for row in provenance if row.key == "input_sha256:" * relative(files.series))
    for row in eachrow(steps)
        isempty(row.marker) && continue
        check!("step $(row.step): its marker exists", "true", row.marker_exists)
        row.marker_exists == "true" || continue
        check!("step $(row.step): run from a clean tree", "false", row.git_dirty)
        check!("step $(row.step): its commit is in this branch's history", "true", row.commit_in_history)
        isempty(row.fy4b_file_sha256) ||
            check!("step $(row.step): read this fy4b_w3_h0.csv (SHA-256)", series_sha, row.fy4b_file_sha256)
    end

    write("checks.csv", checks)
    CSV.write(joinpath(OUTDIR, "provenance.csv"), finish_provenance!(provenance, ROOT))
    failed = checks[.!checks.pass, :]
    println("$(nrow(checks) - nrow(failed)) of $(nrow(checks)) checks pass; tables in $(relative(OUTDIR))")
    if nrow(failed) > 0
        show(failed; allrows=true, truncate=0)
        println()
        error("$(nrow(failed)) FY4B trace check(s) failed; see checks.csv")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

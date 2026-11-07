#!/usr/bin/env julia

# Where does FY4B's error sit, and what could correcting FY4B's own rain number reach? An error
# budget of the raw fields in one benchmark run's out-of-fold files. Fits nothing; a diagnostic with
# no decision rule. The gross-hour thresholds were chosen after a first look at these data, so
# gross_summary.csv also gives them under other thresholds.
#
#   julia --project=. scripts/diagnose_fy4b_error_budget.jl
#
# The run is the paper's `fy4b_benchmark` step, as in `compare_fy4b_quality.jl`.
#
# Every table scores FY4B, FY4B_QM, GSMaP and GPM_uncal on one mask per scheme: all four evaluation
# masks, a gauge value, and all four raw fields finite. Hour labels are Beijing hour ends, so a
# cell's calendar day is the day of `label - 1 h`. The wet threshold is 0.1 mm/h.
#
# Writes output/fy4b_error_budget/:
# - budget.csv: per product, the scores, and the squared error split by false rain (gauge dry,
#   product wet), misses (gauge wet, product dry), hit amounts (both wet) and both dry.
# - budget_gap.csv: each product's MSE minus GPM_uncal's, per part, with a 95% day-block interval.
# - oracles.csv: r and RMSE with one kind of error removed using the gauges. Bounds, not methods.
# - time_scales.csv: r over 1, 3, 6 and 24 h blocks from Beijing midnight, pooled over stations and
#   for the network mean. A station's block needs 75% of its hours; a network value 100 stations.
# - strata.csv: per season and regime, total against the gauges, hourly r and network-day r.
#   Regimes are `compare_fy4b_by_position.jl`'s R1, R2 and R3. A day takes its first hour's regime.
# - concentration.csv: the share of squared error in each product's worst 10, 30 and 90 days.
# - gross_summary.csv, gross_counts.csv and gross_hours.csv: gross false-rain hours, where the
#   product is wet at >= 50% of gauges and the gauges are wet at <= 5%. An hour needs 100 stations.
#   gross_hours.csv (balanced spatial) lists FY4B's and FY4B_QM's, with the UTC start for raw files.
# - wet_share_bins.csv and gross_vs_real.csv: whether a network-wide wet call is real, per product.
# - overlap.csv: P(gauge wet) by the wet calls of FY4B or FY4B_QM and of GPM_uncal.
# - ceiling.csv: squared error at gauges >= 20 mm/h, against FY4B's 20 mm/h cap.
# - gauge_pairs.csv: pooled hourly r between gauges by distance, a bound on any field's point skill.
# - provenance.csv.
#
# Checks; the script errors if one fails: the run exists, the parts of the budget sum to the squared
# error, and FY4B_QM against FY4B on their shared mask reproduces output/fy4b_quality/paired.csv.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("RunProvenance", "TableIO", "DayBlockBootstrap", "SatelliteCalibration",
    "SatelliteTemporalEvaluation")

using CSV, Printf, Random, Statistics
using Main.DayBlockBootstrap: read_wide_values, read_wide_mask, group_sums, sum_stats
using Main.SatelliteCalibration: product_regimes
using Main.SatelliteTemporalEvaluation: SEASONS, season_of
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const REBUILD_DIR = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024")
const FY4B_FILE = joinpath(REBUILD_DIR, "fy4b_w3_h0.csv")
const HOUR_META = joinpath(REBUILD_DIR, "fy4b_hour_meta.csv")
const QUALITY_FILE = joinpath(ROOT, "output", "fy4b_quality", "paired.csv")
const OUTDIR = joinpath(ROOT, "output", "fy4b_error_budget")
const PRODUCTS = ["FY4B", "FY4B_QM", "GSMaP", "GPM_uncal"]
const REFERENCE = "GPM_uncal"
const WET = 0.1
const THRESHOLDS = [0.1, 2.0, 8.0]
const PARTS = ["false_rain", "misses", "hit_amounts", "both_dry"]
const REGIMES = Dict("133.0E_V1.0" => "R1", "133.0E_V1.0.0" => "R2", "105.0E_V1.0.0" => "R3")
const BLOCK_HOURS = [1, 3, 6, 24]
const MIN_BLOCK_SHARE = 0.75
const MIN_NETWORK = 100
const MIN_DAY_HOURS = 18
const MIN_STRATUM_DAYS = 30
const TOP_DAYS = [10, 30, 90]
const DROP_DAYS = 10
const GROSS_WET, GROSS_DRY = 0.5, 0.05
const GROSS_SENSITIVITY = [(0.25, 0.05), (0.5, 0.05), (0.75, 0.05), (0.5, 0.02), (0.5, 0.1)]
const WET_SHARE_BINS = [(0.1, 0.25), (0.25, 0.5), (0.5, 0.75), (0.75, 0.9), (0.9, 1.01)]
const CAP = 20.0
const PAIR_BINS = [(0, 5), (5, 10), (10, 20), (20, 40), (40, 80)]
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261106

"""
The run's configuration, and a copy of its `MGERConfig` that writes to `OUTDIR`: loading the grid
writes `global_common_time_qc.csv`, which must not land in the run.
"""
function run_config()
    config = benchmark_config(:full; nested_covariates=true, satellite_wet_blend=true,
        fused_anchor=true, gpm_uncal=true, blend_axes=[:constant, :agreement_envelope],
        fy4b=true, fy4b_file=FY4B_FILE, fy4b_qm=true, fy4b_qm_adaptive=true, fy4b_qm_fusion=true)
    fields = fieldnames(MGERConfig)
    mger = MGERConfig(; merge(
        NamedTuple{fields}(Tuple(getfield(config.mger, field) for field in fields)),
        (; outdir=OUTDIR))...)
    return config, mger
end

"""Which of `PARTS` a cell falls in: 1 false rain, 2 misses, 3 hit amounts, 4 both dry."""
part_of(x, y) = y >= WET ? (x >= WET ? 3 : 2) : (x >= WET ? 1 : 4)

"""Per day and part, the squared error of `x` on the masked cells."""
function daily_parts(x, y, mask, day_of_time, n_day)
    parts = zeros(n_day, length(PARTS))
    @inbounds for t in axes(y, 2), i in axes(y, 1)
        mask[i, t] || continue
        parts[day_of_time[t], part_of(x[i, t], y[i, t])] += (x[i, t] - y[i, t])^2
    end
    return parts
end

function budget_rows(scheme, x, y, mask, parts, product)
    total = vec(group_sums(x, y, mask, ones(Int, size(y)), 1, THRESHOLDS))
    stats = sum_stats(total, THRESHOLDS)
    by_part = vec(sum(parts; dims=1))
    isapprox(sum(by_part), total[3]; rtol=1e-9) ||
        error("$product: the budget's parts do not sum to the squared error")
    shares = NamedTuple{Tuple(Symbol("share_", p) for p in PARTS)}(Tuple(by_part ./ sum(by_part)))
    mses = NamedTuple{Tuple(Symbol("mse_", p) for p in PARTS)}(Tuple(by_part ./ stats.n))
    return merge((; scheme, product), stats, (; mse=total[3] / stats.n), shares, mses)
end

"""`product`'s MSE minus `REFERENCE`'s, per part and in total, with whole days resampled."""
function gap_rows(scheme, parts, day_n, product)
    keep = findall(>(0), day_n)
    rng = MersenneTwister(BOOTSTRAP_SEED)
    labels = vcat(PARTS, "total")
    draws = Matrix{Float64}(undef, BOOTSTRAP_REPS, length(labels))
    difference = parts[product] .- parts[REFERENCE]
    gap(days) = begin
        part_gap = vec(sum(view(difference, days, :); dims=1)) ./ sum(view(day_n, days))
        vcat(part_gap, sum(part_gap))
    end
    for rep in 1:BOOTSTRAP_REPS
        draws[rep, :] = gap(rand(rng, keep, length(keep)))
    end
    point = gap(keep)
    return [(; scheme, product, baseline=REFERENCE, part=labels[j], delta_mse=point[j],
        lo=quantile(view(draws, :, j), 0.025), hi=quantile(view(draws, :, j), 0.975))
        for j in eachindex(labels)]
end

function oracle_rows(scheme, x, y, mask, product)
    observed, estimate = y[mask], x[mask]
    gauge_wet, product_wet = observed .>= WET, estimate .>= WET
    no_false = ifelse.(.!gauge_wet .& product_wet, 0.0, estimate)
    no_misses = ifelse.(gauge_wet .& .!product_wet, observed, estimate)
    perfect_hits = ifelse.(gauge_wet .& product_wet, observed, estimate)
    score(z) = (cor(observed, z), sqrt(mean((z .- observed) .^ 2)))
    return [(; scheme, product, oracle=name, r=score(z)[1], RMSE=score(z)[2])
        for (name, z) in (("as_is", estimate), ("no_false_rain", no_false), ("no_misses", no_misses),
            ("perfect_hit_amounts", perfect_hits))]
end

"""Per station and block, the sums of `x` over masked cells and how many cells they hold."""
function block_sums(x, mask, block_of_time, n_block)
    sums, counts = zeros(size(x, 1), n_block), zeros(size(x, 1), n_block)
    @inbounds for t in axes(x, 2), i in axes(x, 1)
        mask[i, t] || continue
        sums[i, block_of_time[t]] += x[i, t]
        counts[i, block_of_time[t]] += 1
    end
    return sums, counts
end

"""Average ranks, with ties sharing theirs."""
function ranks(values)
    order = sortperm(values)
    result = zeros(length(values))
    first_tied = 1
    while first_tied <= length(values)
        last_tied = first_tied
        while last_tied < length(values) && values[order[last_tied + 1]] == values[order[first_tied]]
            last_tied += 1
        end
        result[order[first_tied:last_tied]] .= (first_tied + last_tied) / 2
        first_tied = last_tied + 1
    end
    return result
end

spearman(a, b) = cor(ranks(a), ranks(b))

function time_scale_rows(scheme, series, y, mask, times)
    start_hours = [Dates.value(t - Hour(1) - DateTime(2022)) ÷ 3_600_000 for t in times]
    rows = NamedTuple[]
    for hours in BLOCK_HOURS
        labels = fld.(start_hours, hours)
        unique_labels = sort(unique(labels))
        block_of_time = [searchsortedfirst(unique_labels, label) for label in labels]
        gauge, counts = block_sums(y, mask, block_of_time, length(unique_labels))
        complete = counts .>= MIN_BLOCK_SHARE * hours
        network = findall(>=(MIN_NETWORK), vec(sum(complete; dims=1)))
        network_mean(m) = [mean(m[complete[:, b], b]) for b in network]
        gauge_network = network_mean(gauge)
        for product in PRODUCTS
            field, _ = block_sums(series[product], mask, block_of_time, length(unique_labels))
            field_network = network_mean(field)
            push!(rows, (; scheme, product, block_hours=hours, n_station_blocks=count(complete),
                n_network_blocks=length(network), r_station=cor(field[complete], gauge[complete]),
                r_network=cor(field_network, gauge_network),
                spearman_network=spearman(field_network, gauge_network)))
        end
    end
    return rows
end

"""The network mean of every hour with `MIN_NETWORK` scored stations, `NaN` otherwise."""
function hourly_network(x, mask)
    means = fill(NaN, size(x, 2))
    for t in axes(x, 2)
        cells = view(mask, :, t)
        count(cells) >= MIN_NETWORK && (means[t] = mean(view(x, cells, t)))
    end
    return means
end

"""Per day, the sum of the hourly network means, and how many hours it holds."""
function daily_network(hourly, day_of_time, n_day)
    totals, hours = zeros(n_day), zeros(Int, n_day)
    for t in eachindex(hourly)
        isnan(hourly[t]) && continue
        totals[day_of_time[t]] += hourly[t]
        hours[day_of_time[t]] += 1
    end
    return totals, hours
end

function strata_rows(scheme, series, y, mask, day_of_time, days, regime)
    gauge_daily, day_hours = daily_network(hourly_network(y, mask), day_of_time, length(days))
    first_hour = [findfirst(==(d), day_of_time) for d in eachindex(days)]
    day_regime, day_season = regime[first_hour], season_of.(days)
    hour_season = day_season[day_of_time]
    complete_day = day_hours .>= MIN_DAY_HOURS
    daily = Dict(p => first(daily_network(hourly_network(series[p], mask), day_of_time, length(days)))
        for p in PRODUCTS)
    rows = NamedTuple[]
    for season in vcat("all", SEASONS), stratum in ("all", "R1", "R2", "R3")
        picked(s, r) = (season == "all" || s == season) && (stratum == "all" || r == stratum)
        hour_cells = mask .& permutedims(picked.(hour_season, regime))
        chosen_days = findall(complete_day .& picked.(day_season, day_regime))
        length(chosen_days) >= MIN_STRATUM_DAYS || continue
        for product in PRODUCTS
            x = series[product][hour_cells]
            gauge = y[hour_cells]
            field_day, gauge_day = daily[product][chosen_days], gauge_daily[chosen_days]
            worst = partialsortperm(abs.(field_day .- gauge_day), 1:DROP_DAYS; rev=true)
            kept = setdiff(eachindex(field_day), worst)
            push!(rows, (; scheme, season, regime=stratum, product, n_cells=length(x),
                n_days=length(chosen_days), gauge_mean=mean(gauge), total_ratio=sum(x) / sum(gauge),
                r_hourly=cor(x, gauge), r_network_day=cor(field_day, gauge_day),
                spearman_network_day=spearman(field_day, gauge_day),
                r_network_day_drop_worst=cor(field_day[kept], gauge_day[kept])))
        end
    end
    return rows
end

function concentration_rows(scheme, parts, days, product)
    all_error, false_rain = vec(sum(parts; dims=2)), parts[:, 1]
    by_all, by_false = sortperm(all_error; rev=true), sortperm(false_rain; rev=true)
    share(values, order, top) = sum(values[order[1:top]]) / sum(values)
    return [(; scheme, product, n_days=count(>(0), all_error), top_days=top,
        share_squared_error=share(all_error, by_all, top),
        share_false_rain=share(false_rain, by_false, top),
        worst_false_rain_days=join(string.(days[by_false[1:min(top, 10)]]), ";")) for top in TOP_DAYS]
end

"""Per hour with `MIN_NETWORK` scored stations: wet shares, means and squared errors."""
function hour_profile(x, y, mask, t)
    cells = view(mask, :, t)
    n = count(cells)
    n >= MIN_NETWORK || return nothing
    xs, ys = view(x, cells, t), view(y, cells, t)
    false_alarm = (ys .< WET) .& (xs .>= WET)
    errors = (xs .- ys) .^ 2
    return (; n, product_wet=mean(>=(WET), xs), gauge_wet=mean(>=(WET), ys), product_mean=mean(xs),
        gauge_mean=mean(ys), sse=sum(errors), false_rain_sse=sum(errors[false_alarm]))
end

function gross_tables(scheme, series, y, mask, times, regime)
    summary, counts, hours = NamedTuple[], NamedTuple[], NamedTuple[]
    bins, compared = NamedTuple[], NamedTuple[]
    day = Date.(times .- Hour(1))
    gsmap, gpm = hourly_network(series["GSMaP"], mask), hourly_network(series["GPM_uncal"], mask)
    for product in PRODUCTS
        profile = [hour_profile(series[product], y, mask, t) for t in axes(y, 2)]
        scored = findall(!isnothing, profile)
        total_sse = sum(profile[t].sse for t in scored)
        total_false = sum(profile[t].false_rain_sse for t in scored)
        for (wet_at_least, dry_at_most) in GROSS_SENSITIVITY
            gross = [t for t in scored if profile[t].product_wet >= wet_at_least &&
                profile[t].gauge_wet <= dry_at_most]
            push!(summary, (; scheme, product, product_wet_at_least=wet_at_least,
                gauge_wet_at_most=dry_at_most, primary=(wet_at_least, dry_at_most) == (GROSS_WET, GROSS_DRY),
                hours=length(gross), share_hours=length(gross) / length(times),
                days=length(unique(day[gross])),
                share_squared_error=sum(profile[t].sse for t in gross; init=0.0) / total_sse,
                share_false_rain=sum(profile[t].false_rain_sse for t in gross; init=0.0) / total_false))
        end
        gross = [t for t in scored if profile[t].product_wet >= GROSS_WET &&
            profile[t].gauge_wet <= GROSS_DRY]
        for (unit, keys) in (("month", month.(day[gross])), ("hour_end_bjt", hour.(times[gross])),
                ("year", year.(day[gross])), ("regime", regime[gross]))
            for key in sort(unique(keys))
                push!(counts, (; scheme, product, unit, key=string(key), hours=count(==(key), keys)))
            end
        end
        if product in ("FY4B", "FY4B_QM")
            for t in gross
                p = profile[t]
                push!(hours, (; product, time=times[t], utc_hour_start=times[t] - Hour(9),
                    regime=regime[t], n_station=p.n, product_mean=p.product_mean,
                    gauge_mean=p.gauge_mean, gpm_uncal_mean=gpm[t], gsmap_mean=gsmap[t],
                    product_wet_share=p.product_wet, gauge_wet_share=p.gauge_wet, sse=p.sse,
                    false_rain_sse=p.false_rain_sse))
            end
        end
        for (lo, hi) in WET_SHARE_BINS
            chosen = [t for t in scored if lo <= profile[t].product_wet < hi]
            push!(bins, (; scheme, product, wet_share_lo=lo, wet_share_hi=hi, hours=length(chosen),
                share_gross=mean(profile[t].gauge_wet <= GROSS_DRY for t in chosen),
                median_gauge_wet_share=median(profile[t].gauge_wet for t in chosen)))
        end
        widespread = [t for t in scored if profile[t].product_wet >= GROSS_WET]
        for (kind, chosen) in (("gross", filter(t -> profile[t].gauge_wet <= GROSS_DRY, widespread)),
                ("real", filter(t -> profile[t].gauge_wet > GROSS_DRY, widespread)))
            amounts = [profile[t].product_mean for t in chosen]
            push!(compared, (; scheme, product, kind, hours=length(chosen),
                median_product_mean=median(amounts), p90_product_mean=quantile(amounts, 0.9)))
        end
    end
    return summary, counts, hours, bins, compared
end

function overlap_rows(scheme, series, y, mask)
    observed = y[mask]
    other = series[REFERENCE][mask] .>= WET
    rows = NamedTuple[]
    for product in ("FY4B", "FY4B_QM")
        estimate = series[product][mask]
        wet = estimate .>= WET
        correlation = cor(estimate, series[REFERENCE][mask])
        for (call, cells) in (("both_wet", wet .& other), ("product_only", wet .& .!other),
                ("reference_only", .!wet .& other), ("both_dry", .!wet .& .!other))
            push!(rows, (; scheme, product, reference=REFERENCE, call, n=count(cells),
                gauge_wet=mean(>=(WET), observed[cells]), gauge_mean=mean(observed[cells]),
                r_with_reference=correlation))
        end
    end
    return rows
end

function ceiling_rows(scheme, series, y, mask)
    observed = y[mask]
    heavy = observed .>= CAP
    return [begin
        estimate = series[product][mask]
        errors = (estimate .- observed) .^ 2
        (; scheme, product, gauge_at_cap=count(heavy), share_cells=mean(heavy),
            share_squared_error=sum(errors[heavy]) / sum(errors),
            product_at_cap=count(>=(CAP), estimate), both_at_cap=count(heavy .& (estimate .>= CAP)),
            product_max=maximum(estimate))
    end for product in PRODUCTS]
end

function haversine_km(a, b)
    dlat, dlon = deg2rad(b[2] - a[2]), deg2rad(b[1] - a[1])
    s = sin(dlat / 2)^2 + cos(deg2rad(a[2])) * cos(deg2rad(b[2])) * sin(dlon / 2)^2
    return 2 * 6371 * asin(sqrt(s))
end

"""Pooled hourly r between gauges by distance, from running sums: the pairs hold 10^8 cells."""
function gauge_pair_rows(y, mask, lonlat)
    totals = zeros(length(PAIR_BINS), 6)
    pairs = zeros(Int, length(PAIR_BINS))
    for a in axes(y, 1), b in (a + 1):size(y, 1)
        distance = haversine_km(lonlat[a, :], lonlat[b, :])
        k = findfirst(bin -> bin[1] <= distance < bin[2], PAIR_BINS)
        k === nothing && continue
        pairs[k] += 1
        @inbounds for t in axes(y, 2)
            (mask[a, t] && mask[b, t]) || continue
            u, v = y[a, t], y[b, t]
            totals[k, 1] += 1
            totals[k, 2] += u
            totals[k, 3] += v
            totals[k, 4] += u * u
            totals[k, 5] += v * v
            totals[k, 6] += u * v
        end
    end
    return [begin
        n, su, sv, suu, svv, suv = totals[k, :]
        r = (suv / n - su / n * sv / n) / sqrt((suu / n - (su / n)^2) * (svv / n - (sv / n)^2))
        (; km_lo=PAIR_BINS[k][1], km_hi=PAIR_BINS[k][2], pairs=pairs[k], n_cells=Int(n), r)
    end for k in eachindex(PAIR_BINS)]
end

"""Error unless FY4B_QM and FY4B reproduce `QUALITY_FILE` on their shared mask."""
function check_quality(scheme, series, y, masks)
    quality = CSV.read(QUALITY_FILE, DataFrame)
    row = only(filter(r -> r.scheme == scheme && r.product == "FY4B_QM" && r.baseline == "FY4B" &&
        r.method == "raw", quality))
    shared = masks["FY4B_QM"] .& masks["FY4B"] .& .!isnan.(y) .& .!isnan.(series["FY4B_QM"]) .&
        .!isnan.(series["FY4B"])
    rmse(x) = sqrt(mean((x[shared] .- y[shared]) .^ 2))
    ok = count(shared) == row.candidate_n && isapprox(rmse(series["FY4B_QM"]), row.candidate_RMSE; rtol=1e-9) &&
        isapprox(rmse(series["FY4B"]), row.baseline_RMSE; rtol=1e-9)
    ok || error("$scheme: FY4B_QM and FY4B do not reproduce $QUALITY_FILE")
    @printf("  %s reproduces fy4b_quality: n %d, RMSE %.5f and %.5f\n", scheme, count(shared),
        rmse(series["FY4B_QM"]), rmse(series["FY4B"]))
end

function budget_main()
    mkpath(OUTDIR)
    config, mger = run_config()
    run_dir = config.mger.outdir
    isfile(joinpath(run_dir, "metrics_pooled.csv")) || error("the run is not finished: $run_dir")
    schemes = string.(config.cv_schemes)
    path(scheme, product, file) = joinpath(run_dir, scheme, lowercase(product), file)
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path,
        collect(values(mger.sat_paths)), HOUR_META, QUALITY_FILE,
        [path(s, p, file) for s in schemes for p in PRODUCTS
            for file in ("oof_raw.csv", "common_evaluation_mask.csv")])
    provenance = start_provenance(ROOT, OUTDIR; inputs)

    _, ids, data = load_global_common_product_data(mger)
    times = data[REFERENCE].times
    y = Matrix{Float64}(data[REFERENCE].Y_obs)
    labels, regime_index = product_regimes(HOUR_META, times)
    regime = [REGIMES[labels[index]] for index in regime_index]
    day_of_label = Date.(times .- Hour(1))
    days = sort(unique(day_of_label))
    day_of_time = [searchsortedfirst(days, day) for day in day_of_label]
    station_meta = load_station_meta(mger.station_meta_path; station_id_col=mger.station_id_col,
        lon_col=mger.lon_col, lat_col=mger.lat_col)
    lonlat = build_X_lonlat(station_meta, ids)
    @printf("Grid: %d stations x %d hours, %d Beijing days\n", length(ids), length(times), length(days))

    tables = Dict(name => NamedTuple[] for name in ("budget", "budget_gap", "oracles", "time_scales",
        "strata", "concentration", "gross_summary", "gross_counts", "wet_share_bins",
        "gross_vs_real", "overlap", "ceiling"))
    gross_hours = NamedTuple[]
    pair_rows = NamedTuple[]
    for scheme in schemes
        series = Dict(p => read_wide_values(path(scheme, p, "oof_raw.csv"), ids, times) for p in PRODUCTS)
        masks = Dict(p => read_wide_mask(path(scheme, p, "common_evaluation_mask.csv"), ids, times)
            for p in PRODUCTS)
        check_quality(scheme, series, y, masks)
        mask = .!isnan.(y)
        for p in PRODUCTS
            mask .&= masks[p] .& .!isnan.(series[p])
        end
        @printf("%s: %d cells scored of %d\n", scheme, count(mask), length(mask))
        day_n = zeros(length(days))
        for t in axes(y, 2)
            day_n[day_of_time[t]] += count(view(mask, :, t))
        end
        parts = Dict(p => daily_parts(series[p], y, mask, day_of_time, length(days)) for p in PRODUCTS)
        for p in PRODUCTS
            push!(tables["budget"], budget_rows(scheme, series[p], y, mask, parts[p], p))
            append!(tables["oracles"], oracle_rows(scheme, series[p], y, mask, p))
            append!(tables["concentration"], concentration_rows(scheme, parts[p], days, p))
            p == REFERENCE || append!(tables["budget_gap"], gap_rows(scheme, parts, day_n, p))
        end
        append!(tables["time_scales"], time_scale_rows(scheme, series, y, mask, times))
        append!(tables["strata"], strata_rows(scheme, series, y, mask, day_of_time, days, regime))
        summary, counts, hours, bins, compared = gross_tables(scheme, series, y, mask, times, regime)
        append!(tables["gross_summary"], summary)
        append!(tables["gross_counts"], counts)
        append!(tables["wet_share_bins"], bins)
        append!(tables["gross_vs_real"], compared)
        append!(tables["overlap"], overlap_rows(scheme, series, y, mask))
        append!(tables["ceiling"], ceiling_rows(scheme, series, y, mask))
        if scheme == "balanced_spatial"
            gross_hours = hours
            pair_rows = gauge_pair_rows(y, mask, lonlat)
        end
    end

    for (name, rows) in tables
        write_csv_atomic(joinpath(OUTDIR, "$name.csv"), DataFrame(rows))
    end
    write_csv_atomic(joinpath(OUTDIR, "gross_hours.csv"), DataFrame(gross_hours))
    write_csv_atomic(joinpath(OUTDIR, "gauge_pairs.csv"), DataFrame(pair_rows))
    CSV.write(joinpath(OUTDIR, "provenance.csv"), finish_provenance!(provenance, ROOT))

    judged = "balanced_spatial"
    println("\nBudget ($judged), shares of squared error:")
    for row in filter(r -> r.scheme == judged, tables["budget"])
        @printf("  %-10s r %.3f RMSE %.3f POD %.3f FAR %.3f | false rain %.1f%% misses %.1f%% hits %.1f%%\n",
            row.product, row.r, row.RMSE, getproperty(row, Symbol("POD_0.1")),
            getproperty(row, Symbol("FAR_0.1")), 100row.share_false_rain, 100row.share_misses,
            100row.share_hit_amounts)
    end
    println("\nMSE gap to $REFERENCE ($judged, 95% day-block interval):")
    for row in filter(r -> r.scheme == judged, tables["budget_gap"])
        @printf("  %-10s %-12s %+.3f (%+.3f..%+.3f)\n", row.product, row.part, row.delta_mse, row.lo, row.hi)
    end
    println("\nOracles ($judged):")
    for row in filter(r -> r.scheme == judged, tables["oracles"])
        @printf("  %-10s %-20s r %.3f RMSE %.3f\n", row.product, row.oracle, row.r, row.RMSE)
    end
    println("\nTime scales ($judged), r station / network:")
    for row in filter(r -> r.scheme == judged, tables["time_scales"])
        @printf("  %-10s %2d h  %.3f / %.3f\n", row.product, row.block_hours, row.r_station, row.r_network)
    end
    println("\nStrata ($judged), total ratio, network-day r, Spearman, r without the worst $DROP_DAYS days:")
    for row in filter(r -> r.scheme == judged && (r.season == "all" || r.regime == "all" ||
            r.season == "JJA"), tables["strata"])
        @printf("  %-4s %-4s %-10s days %3d  x%.2f  r %.2f  rho %.2f  drop %.2f\n", row.season, row.regime,
            row.product, row.n_days, row.total_ratio, row.r_network_day, row.spearman_network_day,
            row.r_network_day_drop_worst)
    end
    println("\nConcentration ($judged):")
    for row in filter(r -> r.scheme == judged, tables["concentration"])
        @printf("  %-10s top %2d of %d days: squared error %.1f%%, false rain %.1f%%\n", row.product,
            row.top_days, row.n_days, 100row.share_squared_error, 100row.share_false_rain)
    end
    println("\nGross false-rain hours ($judged):")
    for row in filter(r -> r.scheme == judged, tables["gross_summary"])
        @printf("  %-10s wet>=%.2f gauges<=%.2f: %4d hours (%.2f%%), %3d days, squared error %.1f%%, false rain %.1f%%\n",
            row.product, row.product_wet_at_least, row.gauge_wet_at_most, row.hours, 100row.share_hours,
            row.days, 100row.share_squared_error, 100row.share_false_rain)
    end
    println("\nWrote $(relpath(OUTDIR, ROOT))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    budget_main()
end

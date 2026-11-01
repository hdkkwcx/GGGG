#!/usr/bin/env julia

# Exploratory screen of FY4B_QM's final blend. Two questions: how much of FY4B_QM's edge over the
# gauge-only `adw` is FY4B's own, and could a FY4B-aware blend weight move the blend at all? No GWR
# is refitted. It reads the out-of-fold (OOF) predictions of the 2026-10-31 qmfusion run and runs
# in minutes.
#
#   julia --project=. scripts/screen_fy4b_blend_headroom.jl
#
# The run is `benchmark_config` with `--with-fy4b --fy4b-file output/fy4b_rebuild_2022_2024/
# fy4b_w3_h0.csv --fy4b-qm --fy4b-qm-adaptive --fy4b-qm-fusion` (IMERG Uncal, agreement blend).
# Guard: the blends are rebuilt from the OOF files and the weights in `blend_selection.csv`. They
# must equal `oof_blend_agrenv_mgwr.csv` and `oof_blend_mgwr.csv` in every cell, for FY4B, FY4B_QM
# and GPM_uncal in both schemes. That shows this screen's bands and blends are the benchmark's.
#
# `headroom.csv` is exploratory, not leakage-free. Each held-out fold's weights are chosen on the
# other folds' OOF cells, and those predictions came from models trained on stations that include
# the held-out fold. So it only ranks the keys. A key it picks must be confirmed by a benchmark run
# whose weights are chosen on the inner split. K0 is the control for this protocol only; how far
# it lands from the benchmark's own blend is not a measure of the benchmark's tuning bias.
# `attribution.csv` only scores the run's own nested rows.
#
# Writes output/fy4b_blend_headroom/:
# - attribution.csv: two blends against `adw` and against each other, each on the product's own
#   evaluation mask:
#   - the own-anchor blend `blend_mgwr`, which uses only the target product;
#   - the three-product blend `blend_agrenv_mgwr`, banded on FY4B_QM, GPM_uncal and GSMaP.
#   Also FY4B_QM against FY4B on their shared mask.
# - attenuation.csv: raw, mgwr, adw and both blends for FY4B, FY4B_QM and GPM_uncal on their
#   shared mask, overall, per agreement band and per month, with the skill over `adw`.
# - gap_by_band.csv: where FY4B_QM's blend loses to GPM_uncal's, by agreement band × gauge
#   dry/wet × FY4B_QM dry/wet.
# - headroom.csv: cross-fitted weights under the band keys in `KEYS`, each against K0 and against
#   `adw`. `meets_rule` applies the decision rule: on balanced_spatial against K0, RMSE falls by at
#   least `MIN_REL_GAIN` with an interval below 0, CSI and |FBI − 1| at 8 mm/h are not worse
#   beyond their interval, and |bias| grows by at most `BIAS_TOLERANCE`.
# - headroom_lambdas.csv: the weights each key chose per fold.
# - provenance.csv.
#
# Every comparison carries a paired day-block bootstrap interval: whole days are resampled for both
# series at once. The interval is on the change in RMSE, bias, |bias|, and CSI, FBI and |FBI − 1|
# at 8 mm/h.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("SatelliteCalibration", "RunProvenance", "TableIO", "DayBlockBootstrap")

using CSV, Printf, Random, Statistics
using Main: DayBlockBootstrap
using Main.DayBlockBootstrap: read_wide_values, read_wide_mask, day_groups
using Main.SatelliteCalibration: product_regimes
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const FY4B_FILE = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024", "fy4b_w3_h0.csv")
const OUTDIR = joinpath(ROOT, "output", "fy4b_blend_headroom")
const THRESHOLDS = [0.1, 2.0, 8.0]
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261101
"""Products whose OOF files are read, guarded and scored; GSMaP is read for its raw only."""
const PRODUCTS = ["FY4B", FY4B_QM_PRODUCT, "GPM_uncal"]
const METHODS = ["raw", "mgwr", "adw", "blend_mgwr", "blend_agrenv_mgwr"]
const N_BANDS = blend_band_count(:agreement_envelope)
const N_ENVELOPE = length(BLEND_ENVELOPE_EDGES)
const ENVELOPE_LABELS = ["<2.5", "2.5-8", ">=8"]
"""Statistics that get a bootstrap interval on their change."""
const CI_STATS = (:RMSE, :bias, :abs_bias, :CSI_8, :FBI_8, :FBI_8_dev)
"""Decision rule (plan): RMSE gain needed, and the |bias| growth counted as material (mm/h)."""
const MIN_REL_GAIN = 0.003
const BIAS_TOLERANCE = 0.005

const KEYS = [
    (key="K0", uses_other_products=true, description="the 9 agreement bands, as the benchmark"),
    (key="K1", uses_other_products=true, description="9 bands x FY4B_QM itself dry / wet"),
    (key="K2", uses_other_products=true,
        description="9 bands x the month's QM strength: mapped / none or threshold only"),
    (key="K3", uses_other_products=true, description="9 bands plus band 0 as its own blended ID 10"),
    (key="K4", uses_other_products=false,
        description="FY4B_QM-only bands: its own value on the envelope edges"),
]

"""
The configuration of the run being read, and a copy of its `MGERConfig` that writes to `OUTDIR`.
Loading the grid writes `global_common_time_qc.csv`, which must not land in that run.
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

"""`DayBlockBootstrap`'s functions at this screen's thresholds and bootstrap settings."""
group_sums(x, y, mask, group, n_group) = DayBlockBootstrap.group_sums(x, y, mask, group,
    n_group, THRESHOLDS)
sum_stats(total) = DayBlockBootstrap.sum_stats(total, THRESHOLDS)
compare_days(candidate, baseline) = DayBlockBootstrap.compare_days(candidate, baseline,
    THRESHOLDS; ci_stats=CI_STATS, reps=BOOTSTRAP_REPS, seed=BOOTSTRAP_SEED)

"""Agreement band `b` (1-9) as "<products wet>w <envelope>"."""
band_label(b) = b == 0 ? "band 0 (no product wet)" :
    "$((b - 1) ÷ N_ENVELOPE + 1)w $(ENVELOPE_LABELS[(b - 1) % N_ENVELOPE + 1])"

"""
Band IDs for a band × class key. Band 0 stays 0, which means never blended, and band `b` ≥ 1 of
class `c` becomes `b + N_BANDS·(c − 1)`, so the IDs run over 1..N_BANDS·n_class.

Asserts that the map is one-to-one on the (band, class) pairs that occur, and that every band-0
cell keeps ID 0. Returns the IDs and the (band, class) of each ID used.
"""
function combine_band_class(band::Matrix{Int}, class::Matrix{Int}, n_class::Int)
    size(band) == size(class) || throw(DimensionMismatch("band and class differ in shape"))
    out = zeros(Int, size(band))
    pair_of = Dict{Int,Tuple{Int,Int}}()
    @inbounds for index in eachindex(band)
        b, c = band[index], class[index]
        0 <= b <= N_BANDS || error("band $b is outside 0..$N_BANDS")
        1 <= c <= n_class || error("class $c is outside 1..$n_class")
        b == 0 && continue
        id = b + N_BANDS * (c - 1)
        get!(pair_of, id, (b, c)) == (b, c) ||
            error("band ID $id collides: $(pair_of[id]) and $((b, c))")
        out[index] = id
    end
    all(out[band .== 0] .== 0) || error("a band-0 cell got a blended ID")
    all(out[band .> 0] .> 0) || error("a blended band lost its ID")
    return out, pair_of
end

"""
K3's IDs: band `b` ≥ 1 keeps its ID, and band 0 becomes the blended ID `N_BANDS + 1`. Asserts that
exactly the band-0 cells carry that ID.
"""
function band_zero_blended(band::Matrix{Int})
    out = ifelse.(band .== 0, N_BANDS + 1, band)
    (out .== N_BANDS + 1) == (band .== 0) || error("ID $(N_BANDS + 1) is not exactly band 0")
    pair_of = Dict(id => (id == N_BANDS + 1 ? 0 : id, 1) for id in unique(out))
    return out, pair_of
end

"""
Rebuild every product's two blends from the OOF files and the run's weights, and error unless both
equal the stored blends in every cell.
"""
function guard_blends(scheme, selection, data, band, fold_of, threshold)
    for product in PRODUCTS
        raw, mgwr, adw = (data[(product, method)] for method in ("raw", "mgwr", "adw"))
        rebuilt = Dict(method => fill(NaN, size(mgwr)) for method in ("blend_mgwr", "blend_agrenv_mgwr"))
        for fold in sort(unique(fold_of))
            val = findall(==(fold), fold_of)
            function chosen(method)
                matches = findall(row -> row.scheme == scheme && row.product == product &&
                    row.fold == fold && row.method == method, eachrow(selection))
                length(matches) <= 1 || error("$scheme $product fold $fold: " *
                    "$(length(matches)) $method weight rows, expected at most one")
                return isempty(matches) ? nothing : only(matches)
            end
            # A fold the benchmark skipped has no weights and its blends stay NaN, as stored.
            agreement = chosen("blend_agrenv_mgwr")
            if agreement !== nothing
                lambdas = parse.(Float64, split(selection.lambdas[agreement], "|"))
                rebuilt["blend_agrenv_mgwr"][val, :] =
                    banded_blend_prediction(band[val, :], mgwr[val, :], adw[val, :], lambdas)
            end
            constant = chosen("blend_mgwr")
            if constant !== nothing
                rebuilt["blend_mgwr"][val, :] = satellite_wet_blend_prediction(raw[val, :],
                    mgwr[val, :], adw[val, :], Float64(selection.lambda[constant]), threshold)
            end
        end
        for (method, blend) in rebuilt
            mismatched = count(.!isequal.(blend, data[(product, method)]))
            mismatched == 0 ||
                error("$scheme $product: rebuilt $method differs from the OOF file in $mismatched cells")
        end
        println("  guard: $scheme $product, both blends rebuilt in all $(length(mgwr)) cells")
    end
end

"""
Error unless this screen's scores match `metrics_pooled.csv` for each product's `adw` and blends,
on the product's own mask. The scores come from both the run's own `metric_continuous` /
`metric_event` and this screen's sums.
"""
function check_pooled(scheme, pooled, data, masks, y)
    whole = ones(Int, size(y))
    for product in PRODUCTS, method in ("adw", "blend_mgwr", "blend_agrenv_mgwr")
        x = data[(product, method)]
        scored = masks[product] .& .!isnan.(x)
        continuous = metric_continuous(y, x; mask=scored)
        event = metric_event(y, x; mask=scored, thr=8.0)
        summed = sum_stats(vec(group_sums(x, y, scored, whole, 1)))
        stored(group, level) = only(eachrow(pooled[(pooled.scheme .== scheme) .&
            (pooled.product .== product) .& (pooled.method .== method) .& ismissing.(pooled.fold) .&
            (pooled.group .== group) .& (pooled.level .== level), :]))
        overall = stored("overall", "all")
        at8 = stored("event_threshold", "8.0")
        checks = [(continuous.RMSE, overall.RMSE), (continuous.Bias, overall.Bias),
            (continuous.r, overall.r), (event.CSI, at8.CSI), (event.POD, at8.POD),
            (event.FAR, at8.FAR), (summed.RMSE, overall.RMSE), (summed.bias, overall.Bias),
            (summed.r, overall.r), (summed.CSI_8, at8.CSI)]
        all(isapprox(a, b; rtol=1e-8) for (a, b) in checks) ||
            error("$scheme $product $method: scores differ from metrics_pooled.csv: $checks")
        Int(summed.n) == overall.n || error("$scheme $product $method: n differs")
    end
    println("  check: $scheme scores match metrics_pooled.csv (adw and both blends)")
end

"""
Per fold: weights chosen with `select_blend_lambdas` on the other folds' masked OOF cells, then
applied to the fold with `banded_blend_prediction`. Exploratory, see the header.
"""
function crossfit_blend(y, y_sat, band, mgwr, adw, mask, fold_of, n_bands)
    scored_obs = ifelse.(mask, y, NaN)
    out = fill(NaN, size(y))
    rows = NamedTuple[]
    for fold in sort(unique(fold_of))
        train = findall(!=(fold), fold_of)
        val = findall(==(fold), fold_of)
        choice = select_blend_lambdas(scored_obs[train, :], y_sat[train, :], band[train, :],
            mgwr[train, :], adw[train, :], BLEND_LAMBDAS, n_bands)
        choice === nothing && error("fold $fold: no scorable training cell")
        out[val, :] = banded_blend_prediction(band[val, :], mgwr[val, :], adw[val, :],
            choice.lambdas)
        counts = zeros(Int, n_bands)
        for t in axes(y, 2), i in train
            id = band[i, t]
            id > 0 && mask[i, t] && (counts[id] += 1)
        end
        append!(rows, [(; fold, id, lambda=choice.lambdas[id], n_train_cells=counts[id])
            for id in 1:n_bands])
    end
    return out, rows
end

"""The band IDs, ID count and (band, class) labels for each key in `KEYS`."""
function key_bands(band, raw_qm, strength_class, threshold)
    qm_wet = ifelse.(raw_qm .>= threshold, 2, 1)
    k1, k1_pairs = combine_band_class(band, qm_wet, 2)
    k2, k2_pairs = combine_band_class(band, strength_class, 2)
    k3, k3_pairs = band_zero_blended(band)
    k4 = blend_band_matrix(:agreement_envelope, raw_qm, [raw_qm], threshold)
    maximum(k4) <= N_ENVELOPE || error("FY4B_QM-only bands go past $N_ENVELOPE")
    k0_label(id) = band_label(id)
    return Dict(
        "K0" => (band, N_BANDS, k0_label),
        "K1" => (k1, 2 * N_BANDS, id -> haskey(k1_pairs, id) ?
            "$(band_label(k1_pairs[id][1])), FY4B_QM $(k1_pairs[id][2] == 2 ? "wet" : "dry")" :
            "unused"),
        "K2" => (k2, 2 * N_BANDS, id -> haskey(k2_pairs, id) ?
            "$(band_label(k2_pairs[id][1])), $(k2_pairs[id][2] == 1 ? "mapped" : "unmapped")" :
            "unused"),
        "K3" => (k3, N_BANDS + 1, id -> band_label(first(get(k3_pairs, id, (id, 1))))),
        "K4" => (k4, N_ENVELOPE, id -> "FY4B_QM alone $(ENVELOPE_LABELS[id])"),
    )
end

"""Per fold × hour: class 1 where the fold's training gauges chose a map strength above 0, 2
where they chose none or threshold only (`fy4b_qm_strengths.csv`)."""
function strength_class_matrix(strengths, scheme, fold_of, regimes, regime_of_time, month_label)
    mapped = Dict((row.fold, row.regime, row.cell) => row.strength > 0
        for row in eachrow(strengths) if row.scheme == scheme)
    folds = sort(unique(fold_of))
    per_fold = zeros(Int, maximum(folds), length(month_label))
    for fold in folds, t in eachindex(month_label)
        key = (fold, regimes[regime_of_time[t]], month_label[t])
        haskey(mapped, key) || error("$scheme: no strength for fold $fold, $(key[2]), $(key[3])")
        per_fold[fold, t] = mapped[key] ? 1 : 2
    end
    return per_fold[fold_of, :]
end

function screen_main()
    mkpath(OUTDIR)
    reference, mger = screen_config()
    run_dir = reference.mger.outdir
    isdir(run_dir) || error("the run to screen is missing: $run_dir")
    schemes = string.(reference.cv_schemes)
    oof(scheme, product, file) = joinpath(run_dir, scheme, lowercase(product), file)
    oof_files = vcat(
        [oof(s, p, "oof_$m.csv") for s in schemes for p in PRODUCTS for m in METHODS],
        [oof(s, p, "common_evaluation_mask.csv") for s in schemes for p in PRODUCTS],
        [oof(s, "GSMaP", "oof_raw.csv") for s in schemes],
        [joinpath(run_dir, s, "split_common.csv") for s in schemes])
    run_tables = [joinpath(run_dir, name) for name in
        ("blend_selection.csv", "fy4b_qm_strengths.csv", "metrics_pooled.csv")]
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path,
        [mger.sat_paths[name] for name in sort(collect(keys(mger.sat_paths)))],
        _fy4b_hour_meta_path(reference), run_tables, oof_files)
    provenance = start_provenance(ROOT, OUTDIR; inputs)
    started = time()
    elapsed() = round((time() - started) / 60; digits=1)

    _, ids, product_data = load_global_common_product_data(mger)
    times = product_data["FY4B"].times
    y = product_data["FY4B"].Y_obs
    for name in ("GPM_uncal", "GSMaP")
        product_data[name].times == times && isequal(product_data[name].Y_obs, y) ||
            error("$name's gauges or hours differ from FY4B's")
    end
    regimes, regime_of_time = product_regimes(_fy4b_hour_meta_path(reference), times)
    month_label = [@sprintf("%d-%02d", year(t), month(t)) for t in times]
    selection = CSV.read(run_tables[1], DataFrame; types=Dict(:lambdas => String))
    strengths = CSV.read(run_tables[2], DataFrame;
        types=Dict(:cell => String, :regime => String))
    pooled = CSV.read(run_tables[3], DataFrame; types=Dict(:level => String))
    threshold = mger.rain_threshold
    day_group, n_day = day_groups(times, length(ids))
    day_sums(x, mask) = group_sums(x, y, mask, day_group, n_day)
    println("Grid: $(length(ids)) stations x $(length(times)) hours ($(elapsed()) min)")

    attribution_rows = NamedTuple[]
    attenuation_rows = NamedTuple[]
    gap_rows = NamedTuple[]
    headroom_rows = NamedTuple[]
    lambda_rows = NamedTuple[]
    for scheme in schemes
        split_table = CSV.read(joinpath(run_dir, scheme, "split_common.csv"), DataFrame)
        fold_by_id = Dict(string(row.station_id) => row.fold
            for row in eachrow(split_table) if row.scheme == scheme)
        fold_of = [fold_by_id[id] for id in ids]
        data = Dict{Tuple{String,String},Matrix{Float64}}()
        for product in PRODUCTS, method in METHODS
            data[(product, method)] =
                read_wide_values(oof(scheme, product, "oof_$method.csv"), ids, times)
        end
        data[("GSMaP", "raw")] = read_wide_values(oof(scheme, "GSMaP", "oof_raw.csv"), ids, times)
        masks = Dict(product => read_wide_mask(oof(scheme, product, "common_evaluation_mask.csv"),
            ids, times) .& .!isnan.(y) for product in PRODUCTS)
        println("$scheme: OOF files read ($(elapsed()) min)")

        raw_qm = data[(FY4B_QM_PRODUCT, "raw")]
        band = blend_band_matrix(:agreement_envelope, raw_qm,
            [raw_qm, data[("GPM_uncal", "raw")], data[("GSMaP", "raw")]], threshold)
        guard_blends(scheme, selection, data, band, fold_of, threshold)
        check_pooled(scheme, pooled, data, masks, y)

        # Attribution: the run's own nested rows, each pair on the shared cells of its products.
        uses_other(method) = method == "blend_agrenv_mgwr" ? "yes" :
            method == "blend_mgwr" ? "no" : "gauges only"
        pairs = vcat(
            [((p, "blend_mgwr"), (p, "adw")) for p in PRODUCTS],
            [((p, "blend_agrenv_mgwr"), (p, "adw")) for p in PRODUCTS],
            [((p, "blend_agrenv_mgwr"), (p, "blend_mgwr")) for p in PRODUCTS],
            [((FY4B_QM_PRODUCT, m), ("FY4B", m)) for m in ("blend_mgwr", "blend_agrenv_mgwr")])
        for (candidate, baseline) in pairs
            mask = masks[first(candidate)] .& masks[first(baseline)] .&
                .!isnan.(data[candidate]) .& .!isnan.(data[baseline])
            push!(attribution_rows, merge(
                (; scheme, candidate_product=candidate[1], candidate_method=candidate[2],
                    candidate_uses_other_products=uses_other(candidate[2]),
                    baseline_product=baseline[1], baseline_method=baseline[2],
                    mask=first(candidate) == first(baseline) ? "own" : "shared"),
                compare_days(day_sums(data[candidate], mask), day_sums(data[baseline], mask))))
        end

        # Attenuation: every product and method on one shared mask, by band and by month.
        shared = reduce(.&, [masks[p] .& .!isnan.(data[(p, m)]) for p in PRODUCTS for m in METHODS])
        groupings = [("overall", ones(Int, size(y)), ["all"]),
            ("band", band .+ 1, [band_label(b) for b in 0:N_BANDS]),
            ("month", repeat(permutedims(month.(times)), length(ids)),
                [lpad(m, 2, '0') for m in 1:12])]
        for (group, index, levels) in groupings, product in PRODUCTS
            stats = Dict(m => group_sums(data[(product, m)], y, shared, index, length(levels))
                for m in METHODS)
            for (g, level) in enumerate(levels), method in METHODS
                s = sum_stats(stats[method][g, :])
                adw = sum_stats(stats["adw"][g, :])
                push!(attenuation_rows, (; scheme, product, method, group, level, n=Int(s.n),
                    RMSE=s.RMSE, bias=s.bias, r=s.r, CSI_8=s.CSI_8, FBI_8=s.FBI_8,
                    skill_vs_adw=1 - s.RMSE / adw.RMSE))
            end
        end

        # Where FY4B_QM's blend loses to GPM_uncal's: band × gauge wet × FY4B_QM wet.
        qm_blend = data[(FY4B_QM_PRODUCT, "blend_agrenv_mgwr")]
        gpm_blend = data[("GPM_uncal", "blend_agrenv_mgwr")]
        gap_mask = masks[FY4B_QM_PRODUCT] .& masks["GPM_uncal"] .& .!isnan.(qm_blend) .&
            .!isnan.(gpm_blend)
        gap_group = 1 .+ band .+ (N_BANDS + 1) .* (y .>= threshold) .+
            2 * (N_BANDS + 1) .* (raw_qm .>= threshold)
        n_gap = 4 * (N_BANDS + 1)
        qm_sums = group_sums(qm_blend, y, gap_mask, gap_group, n_gap)
        gpm_sums = group_sums(gpm_blend, y, gap_mask, gap_group, n_gap)
        total_gap = sum(qm_sums[:, 3]) - sum(gpm_sums[:, 3])
        for g in 1:n_gap
            qm_sums[g, 1] == 0 && continue
            b = (g - 1) % (N_BANDS + 1)
            gauge_wet = ((g - 1) ÷ (N_BANDS + 1)) % 2 == 1
            qm_wet = (g - 1) ÷ (2 * (N_BANDS + 1)) == 1
            push!(gap_rows, (; scheme, band=b, band_label=band_label(b),
                gauge=gauge_wet ? "wet" : "dry", fy4b_qm=qm_wet ? "wet" : "dry",
                n=Int(qm_sums[g, 1]), RMSE_fy4b_qm=sqrt(qm_sums[g, 3] / qm_sums[g, 1]),
                RMSE_gpm_uncal=sqrt(gpm_sums[g, 3] / gpm_sums[g, 1]),
                delta_sse=qm_sums[g, 3] - gpm_sums[g, 3],
                share_of_gap=(qm_sums[g, 3] - gpm_sums[g, 3]) / total_gap))
        end

        # Headroom: cross-fitted weights under each key, on FY4B_QM's own mask (exploratory).
        strength_class = strength_class_matrix(strengths, scheme, fold_of, regimes,
            regime_of_time, month_label)
        bands = key_bands(band, raw_qm, strength_class, threshold)
        qm(method) = data[(FY4B_QM_PRODUCT, method)]
        series = Dict("adw" => qm("adw"), "benchmark_agrenv" => qm("blend_agrenv_mgwr"),
            "benchmark_constant" => qm("blend_mgwr"))
        for spec in KEYS
            ids_key, n_bands, label = bands[spec.key]
            series[spec.key], rows = crossfit_blend(y, raw_qm, ids_key, qm("mgwr"), qm("adw"),
                masks[FY4B_QM_PRODUCT], fold_of, n_bands)
            append!(lambda_rows, [merge((; scheme, key=spec.key), row, (; label=label(row.id)))
                for row in rows])
            println("  $scheme $(spec.key) cross-fitted ($(elapsed()) min)")
        end
        head_mask = reduce(.&, [isfinite.(x) for x in values(series)]; init=masks[FY4B_QM_PRODUCT])
        sums = Dict(name => day_sums(x, head_mask) for (name, x) in series)
        described = Dict(spec.key => spec for spec in KEYS)
        for name in vcat([spec.key for spec in KEYS], ["benchmark_agrenv", "benchmark_constant"]),
                baseline in ("K0", "adw")
            name == baseline && continue
            spec = get(described, name, nothing)
            row = merge((; scheme, candidate=name, baseline,
                    uses_other_products=spec === nothing ?
                        (name == "benchmark_agrenv" ? "yes" : "no") :
                        (spec.uses_other_products ? "yes" : "no"),
                    n_bands=spec === nothing ? (name == "benchmark_agrenv" ? N_BANDS : 1) :
                        bands[name][2],
                    description=spec === nothing ? "the run's own $(name == "benchmark_agrenv" ?
                        "blend_agrenv_mgwr" : "blend_mgwr") (nested)" : spec.description),
                compare_days(sums[name], sums[baseline]))
            # The rule is judged on the spatial hold-out only; random rows are context.
            meets_rule = scheme == "balanced_spatial" && baseline == "K0" &&
                name in keys(described) && name != "K0" &&
                row.rel_RMSE <= -MIN_REL_GAIN && row.rel_RMSE_hi < 0 &&
                row.delta_CSI_8_hi >= 0 && row.delta_FBI_8_dev_lo <= 0 &&
                row.delta_abs_bias <= BIAS_TOLERANCE
            push!(headroom_rows, merge(row, (; meets_rule)))
        end
    end

    write_csv_atomic(joinpath(OUTDIR, "attribution.csv"), DataFrame(attribution_rows))
    write_csv_atomic(joinpath(OUTDIR, "attenuation.csv"), DataFrame(attenuation_rows))
    write_csv_atomic(joinpath(OUTDIR, "gap_by_band.csv"), DataFrame(gap_rows))
    write_csv_atomic(joinpath(OUTDIR, "headroom.csv"), DataFrame(headroom_rows))
    write_csv_atomic(joinpath(OUTDIR, "headroom_lambdas.csv"), DataFrame(lambda_rows))
    CSV.write(joinpath(OUTDIR, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println("\nAttribution (RMSE change, 95% day-block interval):")
    for row in attribution_rows
        @printf("%-17s %-8s %-17s vs %-8s %-17s %-5s n=%8d  RMSE %.4f vs %.4f  %+.2f%% (%+.2f..%+.2f)  CSI8 %+.4f  FBI8 %+.3f\n",
            row.scheme, row.candidate_product, row.candidate_method, row.baseline_product,
            row.baseline_method, row.mask, Int(row.candidate_n), row.candidate_RMSE,
            row.baseline_RMSE, 100 * row.rel_RMSE, 100 * row.rel_RMSE_lo, 100 * row.rel_RMSE_hi,
            row.delta_CSI_8, row.delta_FBI_8)
    end
    println("\nHeadroom, exploratory (FY4B_QM, cross-fitted weights):")
    for row in headroom_rows
        @printf("%-17s %-18s vs %-4s other=%-3s RMSE %.4f  %+.2f%% (%+.2f..%+.2f)  bias %+.4f  CSI8 %+.4f (%+.4f..%+.4f)  FBI8 %+.3f  rule=%s\n",
            row.scheme, row.candidate, row.baseline, row.uses_other_products, row.candidate_RMSE,
            100 * row.rel_RMSE, 100 * row.rel_RMSE_lo, 100 * row.rel_RMSE_hi, row.delta_bias,
            row.delta_CSI_8, row.delta_CSI_8_lo, row.delta_CSI_8_hi, row.delta_FBI_8,
            row.meets_rule)
    end
    println("Wrote $(relpath(OUTDIR, ROOT)) in $(elapsed()) min")
end

if abspath(PROGRAM_FILE) == @__FILE__
    screen_main()
end

#!/usr/bin/env julia

# Does FY4B help the fusion, with everything else held equal? Scores the qmfusion run (FY4B_QM in
# the merges and the agreement bands) against the out-of-fusion run (FY4B scored, in neither). Fits
# nothing: it reads the two runs' output files.
#
#   julia --project=. scripts/compare_fy4b_fusion_same_hours.jl          # the full runs
#   julia --project=. scripts/compare_fy4b_fusion_same_hours.jl smoke    # the smoke runs, to test this
#
# Why not `compare_fy4b_in_fusion.jl`, which stays as the record of the first check: that one set the
# qmfusion run against the paper run, and those differ in their hours, and so their tuning samples,
# and in every product's selection seed, which depends on the product's place in the list
# (`run_seed` in `InterpolationBenchmarkFold.jl`). These two runs share the product list, the hours
# and the seeds. The runs:
# - with FY4B in the fusion, the paper's `fy4b_benchmark` step:
#     julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates \
#       --satellite-wet-blend --blend-axis agreement_envelope --fused-anchor --with-fy4b \
#       --fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv --fy4b-qm --fy4b-qm-adaptive \
#       --fy4b-qm-fusion
# - out of it: the same with --fy4b-out-of-fusion in place of --fy4b-qm-fusion.
# Both directories come from `benchmark_config`, so the names cannot drift.
#
# Each row of paired.csv is one field in the qmfusion run (candidate) against the same field in the
# out-of-fusion run (baseline), on the cells both score, with paired day-block intervals
# (`DayBlockBootstrap`) on RMSE, bias, |bias|, and CSI, FBI and |FBI − 1| at 8 mm/h. `role` says
# where the fusion reaches the field:
# - `guard`: `raw`, `adw`, `mgwr` and `blend_mgwr` of GPM_uncal, GSMaP, FY4B and FY4B_QM. The fusion
#   never reaches them, so they must be identical, NaNs and evaluation mask included, or the script
#   errors: the two runs would not differ in the fusion alone.
# - `control`: MERGED_MEAN's and MERGED_OLS's `adw`. A gauge-only field, but of a product the fusion
#   defines, so it is watched rather than required to be identical.
# - `treatment`: GPM_uncal's and GSMaP's agreement blend (FY4B in the bands), and every field of
#   MERGED_MEAN and MERGED_OLS but `adw` (FY4B in the merged anchor, and in the bands).
# - `information`: FY4B's and FY4B_QM's agreement blend. Reported, and in no verdict: neither would
#   be an anchor of the main comparison.
#
# `meets_rule` applies the blend screen's decision rule, unchanged from `compare_fy4b_in_fusion.jl`,
# to treatments on balanced_spatial: RMSE −`MIN_REL_GAIN` or better with an interval below 0, CSI and
# |FBI − 1| at 8 mm/h not worse beyond their interval, and |bias| up by at most `BIAS_TOLERANCE` mm/h.
# The verdict, fixed before the full out-of-fusion run (see `verdict`), is printed at the end.
#
# What the verdict does not isolate: a merged product's covariates are selected per fold on its own
# anchor's residuals, so taking FY4B out of the anchor can change the selection even on equal seeds.
# That is a real path from FY4B to the field, but a discrete one, and the day-block intervals do not
# cover it. selection.csv lists both runs' selection per fold; if a merged product meets the rule,
# read its verdict together with how many of its folds selected differently.
#
# Checks; the script errors if one fails: both runs' time grids (`global_common_time_qc.csv`) equal,
# the same station IDs, the gauges equal, split_common.csv equal per scheme, and the guards.
#
# Writes output/fy4b_same_hours_comparison/ (`_smoke` for the smoke runs): paired.csv,
# fold_status.csv (each run's joint-method folds that did not succeed), selection.csv (mgwr's
# covariates and roles per fold in both runs), provenance.csv.
#
# `--crossed` compares another pair: corrected FY4B fitted out of time and fed to the fusion
# (candidate) against the paper run, plain FY4B in the fusion (baseline). Same product list, hours
# and seeds:
#     julia --project=. scripts/compare_fy4b_fusion_same_hours.jl --crossed          # full runs
#     julia --project=. scripts/compare_fy4b_fusion_same_hours.jl smoke --crossed    # smoke runs
# - the candidate is the paper's `benchmark` step plus --fy4b-qm-crossed --fy4b-qm-fusion. Each
#   7-day time fold's FY4B_QM is fitted on the other time folds' training gauges, so the gauges the
#   fusion fits at an hour never shaped the correction it uses at that hour;
# - the baseline is the paper's `benchmark` step itself.
# Roles change with the pair. FY4B is plain in both runs, so its `raw`, `adw`, `mgwr` and
# `blend_mgwr` are guards beside GPM_uncal's and GSMaP's. Every FY4B_QM field is `information`,
# since its fit is what changes, and its mask check is skipped (scored on shared cells). Treatments
# and controls are as above. Intervals resample whole 7-day blocks (`CROSSED_BLOCK_DAYS`), as the
# paper's do.
#
# Its rule is not the one above. It was chosen by the user and fixed before the full crossed run:
# corrected FY4B replaces plain FY4B in the fusion only if it is no worse. A treatment row is
# `no_worse` on balanced_spatial if all of these hold:
# - the upper end of its RMSE interval is below +`NONINFERIORITY_MARGIN`;
# - CSI at 8 mm/h is not lower beyond its interval;
# - |FBI − 1| at 8 mm/h is not worse beyond its interval;
# - |bias| grows by at most `BIAS_TOLERANCE` mm/h.
# `crossed_verdict` needs that of all four agreement blends. `meets_rule` is still reported. Writes
# output/fy4b_qm_crossed_fusion_comparison/ (`_smoke`), with the same files.
#
# `--confirmed` compares the last FY4B try: confirmed FY4B in the fusion (candidate) against the
# paper run, plain FY4B in it (baseline). Same product list, hours and seeds:
#     julia --project=. scripts/compare_fy4b_fusion_same_hours.jl --confirmed          # full runs
#     julia --project=. scripts/compare_fy4b_fusion_same_hours.jl smoke --confirmed    # smoke runs
# - the candidate is the paper's `benchmark` step plus --fy4b-confirm-fusion. The merges and the
#   agreement bands take FY4B only where GPM_uncal or GSMaP also calls the cell wet, and 0 where
#   FY4B is wet alone (`SatelliteFusion.confirmed_by_others`). It reads no gauge and fits nothing;
# - the baseline is the paper's `benchmark` step itself.
# Roles are the default pair's: neither FY4B's nor FY4B_QM's own fields change, so their `raw`,
# `adw`, `mgwr` and `blend_mgwr` are guards beside GPM_uncal's and GSMaP's, and their agreement
# blends are `information`. Intervals resample whole 7-day blocks, with `CONFIRMED_REPS` draws.
#
# Its rule was chosen by the user and fixed before the full run. About 14 FY4B fixes had already been
# judged on these gauges and hours (`Projects_Pending_Handover/fy4b_tries.md`), so a late pass may be
# luck. The bar is therefore stricter: the usual 5% is shared over `N_TRIES` tries. On
# balanced_spatial, confirmed FY4B replaces plain FY4B in the merges and the bands only if all of
# these hold:
# 1. primary: MERGED_MEAN's agreement blend improves RMSE by `MIN_REL_GAIN` or more, and the upper
#    end of its `STRICT_LEVEL` interval is below 0. FY4B carries a third of MERGED_MEAN's weight;
# 2. no worse: each of the four agreement blends (GPM_uncal, GSMaP, MERGED_MEAN, MERGED_OLS) passes
#    `--crossed`'s checks (`failed_checks`, 95% intervals);
# 3. replication: MERGED_MEAN's agreement blend has a lower RMSE in each of FY4B's periods (R1, R2,
#    R3 in time order, one per regime of `fy4b_hour_meta.csv`; point estimates). In the random split
#    it is not significantly worse: the lower end of its 95% interval is at or below 0. (The user
#    changed this from "lower in random" before the run. A side agent had noted that the random split
#    never moved more than 0.12% in the four earlier tests, so a sign check there is close to a coin
#    toss.)
# A control (MERGED adw) moving RMSE by `MIN_REL_GAIN` or more makes the verdict `confounded`.
# Writes output/fy4b_confirmed_fusion_comparison/ (`_smoke`): the same files, with `rel_RMSE_lo_strict`
# and `rel_RMSE_hi_strict` in paired.csv. It also writes these, information only:
# - periods.csv: every treatment's RMSE per period;
# - confirmed_raw.csv: confirmed against plain FY4B as raw fields, with the share of FY4B's wet
#   cells that confirmation sets to 0;
# - ols_weights.csv: MERGED_OLS's mean coefficients in both runs.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("RunProvenance", "TableIO", "DayBlockBootstrap", "TimeBlockFolds")

using CSV, Printf, Random, Statistics
using Main.DayBlockBootstrap: read_wide_values, read_wide_mask, day_groups, group_sums,
    compare_days, sum_stats
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic
using Main.TimeBlockFolds: block_times

const FY4B_FILE = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024", "fy4b_w3_h0.csv")
const THRESHOLDS = [0.1, 2.0, 8.0]
const CI_STATS = (:RMSE, :bias, :abs_bias, :CSI_8, :FBI_8, :FBI_8_dev)
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261101
const MIN_REL_GAIN = 0.003
const BIAS_TOLERANCE = 0.005
const SINGLE = ["GPM_uncal", "GSMaP"]
const MERGED = ["MERGED_MEAN", "MERGED_OLS"]
const FY4B_PRODUCTS = ["FY4B", "FY4B_QM"]
const AGREEMENT = _blend_method("mgwr", :agreement_envelope)
const METHODS = ["raw", "adw", "mgwr", _blend_method("mgwr"), AGREEMENT]
const JOINT_METHODS = ["residual_gwr", "mixed_gwr", "mgwr"]
"""`--crossed`: the length of the blocks its intervals resample, in days."""
const CROSSED_BLOCK_DAYS = 7
"""`--crossed`'s rule: how much worse RMSE may be at the upper end of its interval."""
const NONINFERIORITY_MARGIN = 0.003
"""
`--confirmed`: the FY4B fixes judged on these gauges and hours, this one included. They are listed
in `Projects_Pending_Handover/fy4b_tries.md`.
"""
const N_TRIES = 15
"""`--confirmed`'s level for its primary check: the usual 5% shared over `N_TRIES` tries."""
const STRICT_LEVEL = 1 - 0.05 / N_TRIES
"""`--confirmed`'s bootstrap draws, enough for `STRICT_LEVEL`'s tails (about 33 draws in each)."""
const CONFIRMED_REPS = 20_000
"""`--confirmed`'s primary product: FY4B carries a third of its weight."""
const CONFIRMED_PRIMARY = "MERGED_MEAN"

"""
Where the fusion reaches `method`'s field for `product`: see the header. Under `--crossed` every
FY4B_QM field is `information`, because its fit is what changes.
"""
function field_role(product, method; crossed::Bool=false)
    if product in MERGED
        return method == "adw" ? "control" : "treatment"
    end
    crossed && product == "FY4B_QM" && return "information"
    method == AGREEMENT || return "guard"
    return product in SINGLE ? "treatment" : "information"
end

"""
The checks of `--crossed`'s rule that a paired row fails, by name. An empty list means the
candidate is no worse.
"""
failed_checks(row) = [name for (name, ok) in (
    "RMSE upper end $(@sprintf("%+.2f", 100 * row.rel_RMSE_hi))% not below " *
        "+$(100 * NONINFERIORITY_MARGIN)%" => row.rel_RMSE_hi < NONINFERIORITY_MARGIN,
    "CSI_8 lower" => row.delta_CSI_8_hi >= 0,
    "|FBI_8 - 1| worse" => row.delta_FBI_8_dev_lo <= 0,
    "|bias| up $(@sprintf("%.4f", row.delta_abs_bias)) mm/h" =>
        row.delta_abs_bias <= BIAS_TOLERANCE,
) if !ok]

"""
`--crossed`'s verdict, chosen by the user and fixed before its full run. It is judged on
balanced_spatial, over the agreement blends of GPM_uncal, GSMaP, MERGED_MEAN and MERGED_OLS:
- `confounded` if a control moves RMSE by `MIN_REL_GAIN` or more either way;
- `use corrected FY4B in the fusion` if every one of them is `no_worse`;
- `keep plain FY4B` otherwise, naming each blend that is not and the checks it fails.
"""
function crossed_verdict(rows)
    balanced = filter(row -> row.scheme == "balanced_spatial", rows)
    any(row -> row.role == "control" && abs(row.rel_RMSE) >= MIN_REL_GAIN, balanced) &&
        return "confounded"
    blends = filter(row -> row.candidate == AGREEMENT && row.role == "treatment", balanced)
    sort([row.product for row in blends]) == sort(vcat(SINGLE, MERGED)) ||
        error("expected one agreement blend for each of $(vcat(SINGLE, MERGED))")
    failed = filter(row -> !row.no_worse, blends)
    isempty(failed) && return "use corrected FY4B in the fusion"
    return "keep plain FY4B: " *
        join(["$(row.product) ($(join(failed_checks(row), ", ")))" for row in failed], "; ")
end

"""
`compare_days`' RMSE interval at `level` instead of 95%. It uses the same seed, draws and
resampling, so it widens the interval `compare_days` reports rather than drawing a new one.
"""
function strict_rmse_interval(candidate, baseline, thresholds; level, reps, seed)
    keep = findall(>(0), candidate[:, 1])
    baseline_rmse = sum_stats(vec(sum(baseline; dims=1)), thresholds).RMSE
    rng = MersenneTwister(seed)
    draws = Vector{Float64}(undef, reps)
    for rep in 1:reps
        sampled = rand(rng, keep, length(keep))
        draws[rep] = sum_stats(vec(sum(view(candidate, sampled, :); dims=1)), thresholds).RMSE -
            sum_stats(vec(sum(view(baseline, sampled, :); dims=1)), thresholds).RMSE
    end
    tail = (1 - level) / 2
    return (; rel_RMSE_lo_strict=quantile(draws, tail) / baseline_rmse,
        rel_RMSE_hi_strict=quantile(draws, 1 - tail) / baseline_rmse)
end

"""
FY4B's periods in time order, R1 first, and each cell's period as a stations × hours matrix. There is
one period per regime of the rebuild's `fy4b_hour_meta.csv` (`product_regimes`), whose labels sort
by name, not by time.
"""
function period_groups(meta_path, times, n_station)
    labels, regime_of_time = product_regimes(meta_path, times)
    order = sortperm([minimum(times[regime_of_time .== regime]) for regime in eachindex(labels)])
    period_of_time = invperm(order)[regime_of_time]
    return labels[order], repeat(permutedims(period_of_time), n_station)
end

"""
The checks of `--confirmed`'s rule that fail, by name; empty means confirmed FY4B passes. `periods`
holds the treatments' RMSE per period, and `rows` the paired rows of both schemes.
"""
function confirmed_failures(rows, periods)
    balanced = filter(row -> row.scheme == "balanced_spatial", rows)
    blends = filter(row -> row.candidate == AGREEMENT && row.role == "treatment", balanced)
    sort([row.product for row in blends]) == sort(vcat(SINGLE, MERGED)) ||
        error("expected one agreement blend for each of $(vcat(SINGLE, MERGED))")
    primary = only(filter(row -> row.product == CONFIRMED_PRIMARY, blends))
    failures = String[]
    primary.rel_RMSE <= -MIN_REL_GAIN || push!(failures, @sprintf(
        "%s gains %+.2f%%, not %.1f%% or more", CONFIRMED_PRIMARY, 100 * primary.rel_RMSE,
        -100 * MIN_REL_GAIN))
    primary.rel_RMSE_hi_strict < 0 || push!(failures, @sprintf(
        "%s's %.2f%% interval reaches %+.2f%%", CONFIRMED_PRIMARY, 100 * STRICT_LEVEL,
        100 * primary.rel_RMSE_hi_strict))
    for row in blends
        failed = failed_checks(row)
        isempty(failed) || push!(failures, "$(row.product) worse ($(join(failed, ", ")))")
    end
    primary_periods = filter(row -> row.scheme == "balanced_spatial" &&
        row.product == CONFIRMED_PRIMARY && row.method == AGREEMENT, periods)
    isempty(primary_periods) && error("no period rows for $CONFIRMED_PRIMARY's agreement blend")
    for row in eachrow(primary_periods)
        row.rel_RMSE < 0 || push!(failures, @sprintf(
            "%s not lower in %s (%+.2f%%)", CONFIRMED_PRIMARY, row.period, 100 * row.rel_RMSE))
    end
    # A run without the random split (the smoke runs) cannot pass: the check is named as missing.
    random = filter(row -> row.scheme == "random" && row.product == CONFIRMED_PRIMARY &&
        row.candidate == AGREEMENT, rows)
    length(random) <= 1 || error("more than one random-split row for $CONFIRMED_PRIMARY's agreement blend")
    if isempty(random)
        push!(failures, "no random split in these runs to check $CONFIRMED_PRIMARY in")
    elseif only(random).rel_RMSE_lo > 0
        push!(failures, @sprintf(
            "%s significantly worse in the random split (%+.2f%%, %+.2f..%+.2f)", CONFIRMED_PRIMARY,
            100 * only(random).rel_RMSE, 100 * only(random).rel_RMSE_lo,
            100 * only(random).rel_RMSE_hi))
    end
    return failures
end

"""
`--confirmed`'s verdict, chosen by the user and fixed before its full run (see the header):
- `confounded` if a control moves RMSE by `MIN_REL_GAIN` or more either way on balanced_spatial;
- `use confirmed FY4B in the merges and bands` if `confirmed_failures` is empty;
- `keep plain FY4B` otherwise, naming each check that fails.
"""
function confirmed_verdict(rows, periods)
    any(row -> row.scheme == "balanced_spatial" && row.role == "control" &&
        abs(row.rel_RMSE) >= MIN_REL_GAIN, rows) && return "confounded"
    failures = confirmed_failures(rows, periods)
    isempty(failures) && return "use confirmed FY4B in the merges and bands"
    return "keep plain FY4B: " * join(failures, "; ")
end

"""MERGED_OLS's mean coefficients per scheme in each run (`fused_anchor_selection.csv`)."""
function ols_weights(dirs)
    rows = NamedTuple[]
    for run in keys(dirs)
        table = CSV.read(joinpath(dirs[run], "fused_anchor_selection.csv"), DataFrame)
        terms = filter(name -> name == "intercept" || startswith(name, "beta_"), names(table))
        for scheme in sort(unique(table.scheme))
            ols = filter(row -> row.scheme == scheme && row.variant == "ols", table)
            push!(rows, merge((; run=String(run), scheme, n_fold=nrow(ols)),
                (; (Symbol(term) => mean(ols[!, term]) for term in terms)...)))
        end
    end
    return DataFrame(rows)
end

"""
Fixed before the full out-of-fusion run. Judged on balanced_spatial, over the agreement blends of
GPM_uncal, GSMaP, MERGED_MEAN and MERGED_OLS:
- `confounded` if a control moves RMSE by `MIN_REL_GAIN` or more either way;
- `add FY4B` if a merged product's agreement blend meets the rule and no agreement blend is worse
  beyond its interval;
- `keep FY4B out` if no agreement blend meets the rule;
- `mixed` otherwise.
"""
function verdict(rows)
    balanced = filter(row -> row.scheme == "balanced_spatial", rows)
    any(row -> row.role == "control" && abs(row.rel_RMSE) >= MIN_REL_GAIN, balanced) &&
        return "confounded"
    blends = filter(row -> row.candidate == AGREEMENT && row.role == "treatment", balanced)
    passes = [row.product for row in blends if row.meets_rule]
    worse = any(row -> row.rel_RMSE_lo > 0, blends)
    any(in(MERGED), passes) && !worse && return "add FY4B"
    isempty(passes) && return "keep FY4B out"
    return "mixed"
end

"""
The candidate run's configuration, the baseline's, and a copy of the former's `MGERConfig` that
writes to `outdir`: loading the data writes `global_common_time_qc.csv`, which must not land in
either run. The candidate is the qmfusion run and the baseline the out-of-fusion run. Under
`crossed`, they are the crossed qmfusion run and the paper run. Under `confirmed`, they are the
confirmed-FY4B run and the paper run.
"""
function run_configs(mode, outdir; crossed::Bool=false, confirmed::Bool=false)
    common = (; nested_covariates=true, satellite_wet_blend=true, fused_anchor=true, gpm_uncal=true,
        blend_axes=[:constant, :agreement_envelope], fy4b=true, fy4b_file=FY4B_FILE, fy4b_qm=true,
        fy4b_qm_adaptive=true)
    with = crossed ? benchmark_config(mode; common..., fy4b_qm_crossed=true, fy4b_qm_fusion=true) :
        confirmed ? benchmark_config(mode; common..., fy4b_confirm_fusion=true) :
        benchmark_config(mode; common..., fy4b_qm_fusion=true)
    without = crossed || confirmed ? benchmark_config(mode; common...) :
        benchmark_config(mode; common..., fy4b_out_of_fusion=true)
    fields = fieldnames(MGERConfig)
    mger = MGERConfig(; merge(
        NamedTuple{fields}(Tuple(getfield(with.mger, field) for field in fields)), (; outdir))...)
    return with, without, mger
end

"""Each run's joint-method folds that did not succeed."""
function fold_status(dirs)
    rows = NamedTuple[]
    for run in keys(dirs), row in eachrow(CSV.read(joinpath(dirs[run], "run_status.csv"), DataFrame))
        row.method in JOINT_METHODS && row.status != "success" && push!(rows, (;
            run=String(run), scheme=row.scheme, product=row.product, fold=row.fold,
            method=row.method, status=row.status))
    end
    return DataFrame(rows)
end

"""
`mgwr`'s selected covariates and roles per scheme, product and fold in both runs, in the two runs'
own words, and whether they are the same set.
"""
function fold_selection(dirs)
    read_roles(run) = Dict((row.scheme, row.product, row.fold) =>
            coalesce(row.covariate_selected_roles, "")
        for row in eachrow(CSV.read(joinpath(dirs[run], "run_status.csv"), DataFrame;
            types=Dict(:covariate_selected_roles => String)))
        if row.method == "mgwr")
    with, without = read_roles(:with), read_roles(:without)
    keys(with) == keys(without) || error("the runs score different folds")
    rows = [(; scheme, product, fold, with=with[(scheme, product, fold)],
            without=without[(scheme, product, fold)],
            same=Set(split(with[(scheme, product, fold)], ";")) ==
                Set(split(without[(scheme, product, fold)], ";")))
        for (scheme, product, fold) in sort(collect(keys(with)))]
    return DataFrame(rows)
end

function compare_main(args=ARGS)
    crossed = "--crossed" in args
    confirmed = "--confirmed" in args
    crossed && confirmed && error("--crossed and --confirmed compare different pairs: give one")
    positional = filter(!startswith("--"), args)
    mode = isempty(positional) ? :full : Symbol(only(positional))
    outdir = joinpath(ROOT, "output",
        (crossed ? "fy4b_qm_crossed_fusion_comparison" :
            confirmed ? "fy4b_confirmed_fusion_comparison" : "fy4b_same_hours_comparison") *
        (mode == :full ? "" : "_$mode"))
    mkpath(outdir)
    with, without, mger = run_configs(mode, outdir; crossed, confirmed)
    dirs = (; with=with.mger.outdir, without=without.mger.outdir)
    println("Candidate: $(relpath(dirs.with, ROOT))\nBaseline:  $(relpath(dirs.without, ROOT))")
    for dir in dirs
        isfile(joinpath(dir, "metrics_pooled.csv")) || error("the run is not finished: $dir")
    end
    schemes = string.(with.cv_schemes)
    string.(without.cv_schemes) == schemes || error("the runs have different schemes")
    products = vcat(SINGLE, MERGED, FY4B_PRODUCTS)
    path(run, scheme, product, file) = joinpath(dirs[run], scheme, lowercase(product), file)
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path, collect(values(mger.sat_paths)),
        [joinpath(dirs[run], "global_common_time_qc.csv") for run in keys(dirs)],
        [joinpath(dirs[run], "run_status.csv") for run in keys(dirs)],
        [joinpath(dirs[run], scheme, "split_common.csv") for run in keys(dirs) for scheme in schemes],
        [path(run, s, p, file) for run in keys(dirs) for s in schemes for p in products
            for file in vcat(["oof_$m.csv" for m in METHODS], "common_evaluation_mask.csv")])
    confirmed && append!(inputs,
        [joinpath(dirs[run], "fused_anchor_selection.csv") for run in keys(dirs)],
        [_fy4b_hour_meta_path(with)])
    provenance = start_provenance(ROOT, outdir; inputs)

    println("Checking that the two runs line up:")
    grid(run) = CSV.read(joinpath(dirs[run], "global_common_time_qc.csv"), DataFrame)
    isequal(grid(:with), grid(:without)) || error("the runs have different time grids")
    _, ids, data = load_global_common_product_data(mger)
    times = data[first(SINGLE)].times
    y = data[first(SINGLE)].Y_obs
    @printf("  %d stations, %d hours on both runs' grid\n", length(ids), length(times))
    for scheme in schemes
        folds(run) = CSV.read(joinpath(dirs[run], scheme, "split_common.csv"), DataFrame)
        isequal(folds(:with), folds(:without)) || error("$scheme: the folds differ")
    end
    println("  split_common.csv equal for ", join(schemes, " and "))

    # `--crossed` and `--confirmed` resample whole 7-day blocks: `block_times` maps each hour to its
    # block's start.
    day_group, n_day = day_groups(crossed || confirmed ?
        block_times(times; block_days=CROSSED_BLOCK_DAYS) : times, length(ids))
    day_sums(x, mask) = group_sums(x, y, mask, day_group, n_day, THRESHOLDS)
    reps = confirmed ? CONFIRMED_REPS : BOOTSTRAP_REPS
    # `--confirmed`'s periods, and what it keeps from the loop below for confirmed_raw.csv.
    period_labels, period_group = confirmed ?
        period_groups(_fy4b_hour_meta_path(with), times, length(ids)) : (String[], nothing)
    period_rows = NamedTuple[]
    raw_fields = Dict{Tuple{String,String},Matrix{Float64}}()
    fy4b_shared = Dict{String,BitMatrix}()
    rows = NamedTuple[]
    for scheme in schemes, product in products
        mask(run) = read_wide_mask(path(run, scheme, product, "common_evaluation_mask.csv"), ids, times)
        mask_with, mask_without = mask(:with), mask(:without)
        series = Dict(method => (;
            with=read_wide_values(path(:with, scheme, product, "oof_$method.csv"), ids, times),
            without=read_wide_values(path(:without, scheme, product, "oof_$method.csv"), ids, times))
            for method in METHODS)
        if product ∉ MERGED && field_role(product, "raw"; crossed) == "guard"
            mask_with == mask_without || error("$scheme $product: the evaluation masks differ")
        end
        for method in METHODS
            field_role(product, method; crossed) == "guard" &&
                !isequal(series[method].with, series[method].without) &&
                error("$scheme $product $method differs between the runs, but the fusion never " *
                    "reaches it")
        end
        shared = mask_with .& mask_without .& .!isnan.(y)
        for method in METHODS, run in (:with, :without)
            shared .&= .!isnan.(series[method][run])
        end
        @printf("  %-16s %-11s shared %d cells (candidate mask %d, baseline mask %d)\n", scheme,
            product, count(shared), count(mask_with), count(mask_without))
        if confirmed && product in ("FY4B", "GPM_uncal", "GSMaP")
            raw_fields[(scheme, product)] = series["raw"].with
            product == "FY4B" && (fy4b_shared[scheme] = shared)
        end
        for method in METHODS
            role = field_role(product, method; crossed)
            candidate_sums = day_sums(series[method].with, shared)
            baseline_sums = day_sums(series[method].without, shared)
            row = merge((; scheme, product, candidate=method, baseline=method, role),
                compare_days(candidate_sums, baseline_sums, THRESHOLDS;
                    ci_stats=CI_STATS, reps, seed=BOOTSTRAP_SEED))
            meets_rule = scheme == "balanced_spatial" && role == "treatment" &&
                row.rel_RMSE <= -MIN_REL_GAIN && row.rel_RMSE_hi < 0 &&
                row.delta_CSI_8_hi >= 0 && row.delta_FBI_8_dev_lo <= 0 &&
                row.delta_abs_bias <= BIAS_TOLERANCE
            no_worse = scheme == "balanced_spatial" && role == "treatment" &&
                isempty(failed_checks(row))
            # `--crossed`'s and `--confirmed`'s rules, in their own columns so the default pair's
            # paired.csv stays as it was.
            if confirmed
                push!(rows, merge(row, strict_rmse_interval(candidate_sums, baseline_sums,
                    THRESHOLDS; level=STRICT_LEVEL, reps, seed=BOOTSTRAP_SEED),
                    (; meets_rule, no_worse)))
                role == "treatment" || continue
                by_period(x) = group_sums(x, y, shared, period_group, length(period_labels),
                    THRESHOLDS)
                with_periods, without_periods = by_period(series[method].with),
                    by_period(series[method].without)
                for (index, label) in enumerate(period_labels)
                    c = sum_stats(with_periods[index, :], THRESHOLDS)
                    b = sum_stats(without_periods[index, :], THRESHOLDS)
                    push!(period_rows, (; scheme, product, method, period="R$index",
                        regime=label, n=Int(c.n), candidate_RMSE=c.RMSE, baseline_RMSE=b.RMSE,
                        rel_RMSE=(c.RMSE - b.RMSE) / b.RMSE))
                end
            else
                push!(rows, crossed ? merge(row, (; meets_rule, no_worse)) :
                    merge(row, (; meets_rule)))
            end
        end
    end
    println("  guards: every field the fusion never reaches is identical in both runs")

    # Information: what confirmation does to FY4B on its own, rebuilt from the raw fields exactly as
    # the run confirms it, and scored on FY4B's shared cells.
    confirmed_raw = NamedTuple[]
    for scheme in (confirmed ? schemes : String[])
        plain = raw_fields[(scheme, "FY4B")]
        others = [raw_fields[(scheme, "GPM_uncal")], raw_fields[(scheme, "GSMaP")]]
        confirmed_fy4b = confirmed_by_others(plain, others, with.mger.rain_threshold)
        shared = fy4b_shared[scheme]
        wet = shared .& (plain .>= with.mger.rain_threshold)
        push!(confirmed_raw, merge(
            (; scheme, fy4b_wet_cells=count(wet),
                zeroed_wet_share=count(wet .& (confirmed_fy4b .== 0)) / count(wet)),
            compare_days(day_sums(confirmed_fy4b, shared), day_sums(plain, shared), THRESHOLDS;
                ci_stats=CI_STATS, reps, seed=BOOTSTRAP_SEED)))
    end
    periods = DataFrame(period_rows)

    statuses = fold_status(dirs)
    selection = fold_selection(dirs)
    write_csv_atomic(joinpath(outdir, "paired.csv"), DataFrame(rows))
    write_csv_atomic(joinpath(outdir, "fold_status.csv"), statuses)
    write_csv_atomic(joinpath(outdir, "selection.csv"), selection)
    if confirmed
        write_csv_atomic(joinpath(outdir, "periods.csv"), periods)
        write_csv_atomic(joinpath(outdir, "confirmed_raw.csv"), DataFrame(confirmed_raw))
        write_csv_atomic(joinpath(outdir, "ols_weights.csv"), ols_weights(dirs))
    end
    CSV.write(joinpath(outdir, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println(crossed ?
        "\nCrossed FY4B_QM in the fusion vs plain FY4B in it (the paper run), same field, shared cells (95% 7-day-block interval):" :
        confirmed ?
        "\nConfirmed FY4B in the fusion vs plain FY4B in it (the paper run), same field, shared cells (95% 7-day-block interval; $(round(100 * STRICT_LEVEL; digits=2))% RMSE interval in brackets):" :
        "\nWith FY4B_QM in the fusion vs out of it, same field, shared cells (95% day-block interval):")
    for row in rows
        @printf("%-17s %-11s %-18s %-11s RMSE %.4f vs %.4f  %+.2f%% (%+.2f..%+.2f)  bias %+.4f vs %+.4f  CSI8 %+.4f (%+.4f..%+.4f)  |FBI8-1| %+.3f (%+.3f..%+.3f)  rule=%s",
            row.scheme, row.product, row.candidate, row.role, row.candidate_RMSE,
            row.baseline_RMSE, 100 * row.rel_RMSE, 100 * row.rel_RMSE_lo, 100 * row.rel_RMSE_hi,
            row.candidate_bias, row.baseline_bias, row.delta_CSI_8, row.delta_CSI_8_lo,
            row.delta_CSI_8_hi, row.delta_FBI_8_dev, row.delta_FBI_8_dev_lo,
            row.delta_FBI_8_dev_hi, row.meets_rule)
        println(crossed ? "  no_worse=$(row.no_worse)" :
            confirmed ? @sprintf("  [%+.2f..%+.2f]  no_worse=%s", 100 * row.rel_RMSE_lo_strict,
                100 * row.rel_RMSE_hi_strict, row.no_worse) : "")
    end
    println("\nJoint-method folds that did not succeed (mgwr; fold_status.csv has the rest):")
    for row in eachrow(statuses)
        row.method == "mgwr" && @printf("  %-8s %-17s %-11s fold %d %s\n", row.run, row.scheme,
            row.product, row.fold, row.status)
    end
    println("\nmgwr's selection per fold, candidate | baseline, where the runs differ (selection.csv has all):")
    for row in eachrow(selection)
        row.same || @printf("  %-17s %-11s fold %d  %s | %s\n", row.scheme, row.product, row.fold,
            row.with, row.without)
    end
    for scheme in schemes, product in MERGED
        cells = filter(row -> row.scheme == scheme && row.product == product, selection)
        @printf("  %-17s %-11s %d of %d folds select the same\n", scheme, product,
            count(cells.same), nrow(cells))
    end
    if crossed
        passing = [row.product for row in rows if row.meets_rule && row.candidate == AGREEMENT]
        println("\nUsual must-improve rule (information only): ",
            isempty(passing) ? "no agreement blend meets it" : "met by " * join(passing, ", "))
        println("Verdict (balanced_spatial, no-worse rule fixed before the run): ",
            crossed_verdict(rows))
    elseif confirmed
        println("\nAgreement blends per FY4B period, candidate vs baseline RMSE (periods.csv has every treatment):")
        for row in eachrow(periods)
            row.method == AGREEMENT && @printf("  %-17s %-11s %s %-15s n %8d  %.4f vs %.4f  %+.2f%%\n",
                row.scheme, row.product, row.period, row.regime, row.n, row.candidate_RMSE,
                row.baseline_RMSE, 100 * row.rel_RMSE)
        end
        println("\nConfirmed vs plain FY4B as raw fields, FY4B's shared cells (information):")
        for row in confirmed_raw
            @printf("  %-17s zeroes %.1f%% of FY4B's %d wet cells; RMSE %.4f vs %.4f %+.2f%% (%+.2f..%+.2f); CSI0.1 %.3f vs %.3f; POD0.1 %.3f vs %.3f; FAR0.1 %.3f vs %.3f; FBI8 %.3f vs %.3f\n",
                row.scheme, 100 * row.zeroed_wet_share, row.fy4b_wet_cells, row.candidate_RMSE,
                row.baseline_RMSE, 100 * row.rel_RMSE, 100 * row.rel_RMSE_lo, 100 * row.rel_RMSE_hi,
                row.var"candidate_CSI_0.1", row.var"baseline_CSI_0.1", row.var"candidate_POD_0.1",
                row.var"baseline_POD_0.1", row.var"candidate_FAR_0.1", row.var"baseline_FAR_0.1",
                row.candidate_FBI_8, row.baseline_FBI_8)
        end
        println("\nMERGED_OLS mean coefficients (ols_weights.csv):")
        show(stdout, ols_weights(dirs); allrows=true, allcols=true, summary=false, eltypes=false)
        passing = [row.product for row in rows if row.meets_rule && row.candidate == AGREEMENT]
        println("\n\nUsual must-improve rule (information only): ",
            isempty(passing) ? "no agreement blend meets it" : "met by " * join(passing, ", "))
        println("Verdict (balanced_spatial, the stricter rule fixed before the run; ",
            "$(round(100 * STRICT_LEVEL; digits=2))% for the primary, $N_TRIES tries): ",
            confirmed_verdict(rows, periods))
    else
        println("\nVerdict (balanced_spatial, fixed before the run): ", verdict(rows))
    end
    println("Wrote $(relpath(outdir, ROOT))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_main()
end

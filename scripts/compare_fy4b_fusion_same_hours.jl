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

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("RunProvenance", "TableIO", "DayBlockBootstrap", "TimeBlockFolds")

using CSV, Printf, Statistics
using Main.DayBlockBootstrap: read_wide_values, read_wide_mask, day_groups, group_sums,
    compare_days
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
either run. The candidate is the qmfusion run and the baseline the out-of-fusion run, or, under
`crossed`, the crossed qmfusion run and the paper run.
"""
function run_configs(mode, outdir; crossed::Bool=false)
    common = (; nested_covariates=true, satellite_wet_blend=true, fused_anchor=true, gpm_uncal=true,
        blend_axes=[:constant, :agreement_envelope], fy4b=true, fy4b_file=FY4B_FILE, fy4b_qm=true,
        fy4b_qm_adaptive=true)
    with = crossed ? benchmark_config(mode; common..., fy4b_qm_crossed=true, fy4b_qm_fusion=true) :
        benchmark_config(mode; common..., fy4b_qm_fusion=true)
    without = crossed ? benchmark_config(mode; common...) :
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
    positional = filter(!startswith("--"), args)
    mode = isempty(positional) ? :full : Symbol(only(positional))
    outdir = joinpath(ROOT, "output",
        (crossed ? "fy4b_qm_crossed_fusion_comparison" : "fy4b_same_hours_comparison") *
        (mode == :full ? "" : "_$mode"))
    mkpath(outdir)
    with, without, mger = run_configs(mode, outdir; crossed)
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

    # `--crossed` resamples whole 7-day blocks: `block_times` maps each hour to its block's start.
    day_group, n_day = day_groups(
        crossed ? block_times(times; block_days=CROSSED_BLOCK_DAYS) : times, length(ids))
    day_sums(x, mask) = group_sums(x, y, mask, day_group, n_day, THRESHOLDS)
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
        for method in METHODS
            role = field_role(product, method; crossed)
            row = merge((; scheme, product, candidate=method, baseline=method, role),
                compare_days(day_sums(series[method].with, shared),
                    day_sums(series[method].without, shared), THRESHOLDS;
                    ci_stats=CI_STATS, reps=BOOTSTRAP_REPS, seed=BOOTSTRAP_SEED))
            meets_rule = scheme == "balanced_spatial" && role == "treatment" &&
                row.rel_RMSE <= -MIN_REL_GAIN && row.rel_RMSE_hi < 0 &&
                row.delta_CSI_8_hi >= 0 && row.delta_FBI_8_dev_lo <= 0 &&
                row.delta_abs_bias <= BIAS_TOLERANCE
            # `--crossed`'s rule, in its own column so the default pair's paired.csv stays as it was.
            push!(rows, crossed ? merge(row, (; meets_rule, no_worse=scheme == "balanced_spatial" &&
                role == "treatment" && isempty(failed_checks(row)))) : merge(row, (; meets_rule)))
        end
    end
    println("  guards: every field the fusion never reaches is identical in both runs")

    statuses = fold_status(dirs)
    selection = fold_selection(dirs)
    write_csv_atomic(joinpath(outdir, "paired.csv"), DataFrame(rows))
    write_csv_atomic(joinpath(outdir, "fold_status.csv"), statuses)
    write_csv_atomic(joinpath(outdir, "selection.csv"), selection)
    CSV.write(joinpath(outdir, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println(crossed ?
        "\nCrossed FY4B_QM in the fusion vs plain FY4B in it (the paper run), same field, shared cells (95% 7-day-block interval):" :
        "\nWith FY4B_QM in the fusion vs out of it, same field, shared cells (95% day-block interval):")
    for row in rows
        @printf("%-17s %-11s %-18s %-11s RMSE %.4f vs %.4f  %+.2f%% (%+.2f..%+.2f)  bias %+.4f vs %+.4f  CSI8 %+.4f (%+.4f..%+.4f)  |FBI8-1| %+.3f (%+.3f..%+.3f)  rule=%s",
            row.scheme, row.product, row.candidate, row.role, row.candidate_RMSE,
            row.baseline_RMSE, 100 * row.rel_RMSE, 100 * row.rel_RMSE_lo, 100 * row.rel_RMSE_hi,
            row.candidate_bias, row.baseline_bias, row.delta_CSI_8, row.delta_CSI_8_lo,
            row.delta_CSI_8_hi, row.delta_FBI_8_dev, row.delta_FBI_8_dev_lo,
            row.delta_FBI_8_dev_hi, row.meets_rule)
        println(crossed ? "  no_worse=$(row.no_worse)" : "")
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
    else
        println("\nVerdict (balanced_spatial, fixed before the run): ", verdict(rows))
    end
    println("Wrote $(relpath(outdir, ROOT))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_main()
end

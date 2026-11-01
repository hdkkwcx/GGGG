#!/usr/bin/env julia

# Does FY4B help the fusion? Score the paper run's fields (no FY4B) against the qmfusion run's
# (FY4B_QM in the merges and the agreement bands) on the cells both runs score. Fits nothing: it
# reads the two runs' output files.
#
#   julia --project=. scripts/compare_fy4b_in_fusion.jl
#
# The two runs:
# - without FY4B, the paper's `benchmark` step:
#     julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates \
#       --satellite-wet-blend --blend-axis agreement_envelope --fused-anchor
# - with FY4B_QM fed to the fusion (2026-10-31): the same plus
#     --with-fy4b --fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv --fy4b-qm \
#     --fy4b-qm-adaptive --fy4b-qm-fusion
# Both directories come from `benchmark_config`, so the names cannot drift.
#
# Shared cells, per scheme and product: FY4B's hours (the qmfusion run's grid), both runs'
# evaluation masks, a gauge value, and every method in `METHODS` finite in both runs.
#
# Each row of paired.csv is one field in the FY4B run (candidate) against the same field in the
# paper run (baseline), with paired day-block intervals (`DayBlockBootstrap`) on RMSE, bias,
# |bias|, and CSI, FBI and |FBI − 1| at 8 mm/h. `role` says where FY4B reaches the field:
# - `guard`: GPM_uncal's and GSMaP's `raw`, which never see FY4B. They must be identical.
# - `control`: fields FY4B never reaches: `adw` everywhere, and GPM_uncal's and GSMaP's `mgwr` and
#   `blend_mgwr`. They can still move, because each run tunes ADW's power, the MGWR bandwidth and
#   the blend weights on hours drawn from its own grid. So they measure that confound.
# - `treatment`: GPM_uncal's and GSMaP's agreement blend (FY4B in the bands), and every field of
#   MERGED_MEAN and MERGED_OLS but `adw` (FY4B in the merged anchor, and in the bands).
#
# `meets_rule` applies the blend screen's decision rule to treatments on balanced_spatial: RMSE
# −`MIN_REL_GAIN` or better with an interval below 0, CSI and |FBI − 1| at 8 mm/h not worse beyond
# their interval, and |bias| up by at most `BIAS_TOLERANCE` mm/h. The verdict, fixed before the
# run (see `verdict`), is printed at the end.
#
# Checks; the script errors if one fails: the same station IDs, every FY4B-run hour in the paper
# run's grid, the gauges equal on those hours, split_common.csv equal per scheme, and the guards.
#
# Writes output/fy4b_in_fusion_comparison/: paired.csv, provenance.csv.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("RunProvenance", "TableIO", "DayBlockBootstrap")

using CSV, Printf, Statistics
using Main.DayBlockBootstrap: read_wide_values, read_wide_mask, day_groups, group_sums,
    compare_days
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const FY4B_FILE = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024", "fy4b_w3_h0.csv")
const OUTDIR = joinpath(ROOT, "output", "fy4b_in_fusion_comparison")
const THRESHOLDS = [0.1, 2.0, 8.0]
const CI_STATS = (:RMSE, :bias, :abs_bias, :CSI_8, :FBI_8, :FBI_8_dev)
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261101
const MIN_REL_GAIN = 0.003
const BIAS_TOLERANCE = 0.005
const SINGLE = ["GPM_uncal", "GSMaP"]
const MERGED = ["MERGED_MEAN", "MERGED_OLS"]
const AGREEMENT = _blend_method("mgwr", :agreement_envelope)
const METHODS = ["raw", "adw", "mgwr", _blend_method("mgwr"), AGREEMENT]

"""Where FY4B reaches `method`'s field for `product`: `guard`, `control` or `treatment`."""
function field_role(product, method)
    method == "adw" && return "control"
    product in MERGED && return "treatment"
    method == "raw" && return "guard"
    return method == AGREEMENT ? "treatment" : "control"
end

"""
The paper run's configuration, the FY4B run's, and copies of their `MGERConfig`s that write to
`OUTDIR`: loading a grid writes `global_common_time_qc.csv`, which must land in neither run.
"""
function run_configs()
    common = (; nested_covariates=true, satellite_wet_blend=true, fused_anchor=true,
        gpm_uncal=true, blend_axes=[:constant, :agreement_envelope])
    without = benchmark_config(:full; common..., fy4b=false)
    with = benchmark_config(:full; common..., fy4b=true, fy4b_file=FY4B_FILE, fy4b_qm=true,
        fy4b_qm_adaptive=true, fy4b_qm_fusion=true)
    fields = fieldnames(MGERConfig)
    redirect(mger) = MGERConfig(; merge(
        NamedTuple{fields}(Tuple(getfield(mger, field) for field in fields)), (; outdir=OUTDIR))...)
    return without, with, redirect(without.mger), redirect(with.mger)
end

"""
Fixed before the run. Judged on balanced_spatial:
- `confounded` if any control moves RMSE by `MIN_REL_GAIN` or more either way, because the
  tuning hours alone would then be as large as the bar;
- `add FY4B` if a merged product's agreement blend meets the rule and no product's agreement blend
  is worse beyond its interval;
- `keep FY4B out` if no product's agreement blend meets the rule;
- `mixed` otherwise.
"""
function verdict(rows)
    balanced = filter(row -> row.scheme == "balanced_spatial", rows)
    any(row -> row.role == "control" && abs(row.rel_RMSE) >= MIN_REL_GAIN, balanced) &&
        return "confounded"
    blends = filter(row -> row.candidate == AGREEMENT, balanced)
    passes = [row.product for row in blends if row.meets_rule]
    worse = any(row -> row.rel_RMSE_lo > 0, blends)
    any(in(MERGED), passes) && !worse && return "add FY4B"
    isempty(passes) && return "keep FY4B out"
    return "mixed"
end

function compare_main()
    mkpath(OUTDIR)
    without, with, mger_without, mger_with = run_configs()
    dirs = (; without=without.mger.outdir, with=with.mger.outdir)
    for dir in dirs
        isfile(joinpath(dir, "metrics_pooled.csv")) || error("the run is not finished: $dir")
    end
    schemes = string.(with.cv_schemes)
    string.(without.cv_schemes) == schemes || error("the runs have different schemes")
    products = vcat(SINGLE, MERGED)
    path(run, scheme, product, file) = joinpath(dirs[run], scheme, lowercase(product), file)
    inputs = vcat(mger_with.station_meta_path, mger_with.obs_hourly_wide_path,
        unique(vcat(collect(values(mger_without.sat_paths)), collect(values(mger_with.sat_paths)))),
        [joinpath(dirs[run], scheme, "split_common.csv") for run in keys(dirs) for scheme in schemes],
        [path(run, s, p, file) for run in keys(dirs) for s in schemes for p in products
            for file in vcat(["oof_$m.csv" for m in METHODS], "common_evaluation_mask.csv")])
    provenance = start_provenance(ROOT, OUTDIR; inputs)

    println("Checking that the two runs line up:")
    _, ids, data_without = load_global_common_product_data(mger_without)
    ids_with, data_with = load_global_common_product_data(mger_with)[2:3]
    ids_with == ids || error("the runs have different stations")
    times_without = data_without[first(SINGLE)].times
    times = data_with[first(SINGLE)].times
    hours = indexin(times, times_without)
    any(isnothing, hours) && error("some FY4B-run hours are not in the paper run's grid")
    hours = Int.(hours)
    y = data_with[first(SINGLE)].Y_obs
    isequal(data_without[first(SINGLE)].Y_obs[:, hours], y) ||
        error("the gauges differ on the shared hours")
    @printf("  %d stations; %d of the paper run's %d hours; gauges equal on them\n",
        length(ids), length(times), length(times_without))
    for scheme in schemes
        folds(run) = CSV.read(joinpath(dirs[run], scheme, "split_common.csv"), DataFrame)
        isequal(folds(:without), folds(:with)) || error("$scheme: the folds differ")
    end
    println("  split_common.csv equal for ", join(schemes, " and "))

    day_group, n_day = day_groups(times, length(ids))
    day_sums(x, mask) = group_sums(x, y, mask, day_group, n_day, THRESHOLDS)
    rows = NamedTuple[]
    for scheme in schemes, product in products
        read_without(file) = read_wide_values(path(:without, scheme, product, file), ids,
            times_without)[:, hours]
        mask_without = read_wide_mask(path(:without, scheme, product, "common_evaluation_mask.csv"),
            ids, times_without)[:, hours]
        mask_with = read_wide_mask(path(:with, scheme, product, "common_evaluation_mask.csv"),
            ids, times)
        series = Dict(method => (; without=read_without("oof_$method.csv"),
            with=read_wide_values(path(:with, scheme, product, "oof_$method.csv"), ids, times))
            for method in METHODS)
        shared = mask_without .& mask_with .& .!isnan.(y)
        for method in METHODS, run in (:without, :with)
            shared .&= .!isnan.(series[method][run])
        end
        @printf("  %-16s %-11s shared %d cells (paper mask %d, FY4B-run mask %d)\n", scheme,
            product, count(shared), count(mask_without), count(mask_with))
        for method in METHODS
            role = field_role(product, method)
            candidate, baseline = series[method].with, series[method].without
            role == "guard" && candidate[shared] != baseline[shared] &&
                error("$scheme $product $method differs between the runs, but FY4B never reaches it")
            row = merge((; scheme, product, candidate=method, baseline=method, role),
                compare_days(day_sums(candidate, shared), day_sums(baseline, shared), THRESHOLDS;
                    ci_stats=CI_STATS, reps=BOOTSTRAP_REPS, seed=BOOTSTRAP_SEED))
            meets_rule = scheme == "balanced_spatial" && role == "treatment" &&
                row.rel_RMSE <= -MIN_REL_GAIN && row.rel_RMSE_hi < 0 &&
                row.delta_CSI_8_hi >= 0 && row.delta_FBI_8_dev_lo <= 0 &&
                row.delta_abs_bias <= BIAS_TOLERANCE
            push!(rows, merge(row, (; meets_rule)))
        end
    end
    println("  guards: GPM_uncal's and GSMaP's raw are identical in both runs")

    write_csv_atomic(joinpath(OUTDIR, "paired.csv"), DataFrame(rows))
    CSV.write(joinpath(OUTDIR, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println("\nWith FY4B_QM vs without, same field, shared cells (95% day-block interval):")
    for row in rows
        @printf("%-17s %-11s %-18s %-9s RMSE %.4f vs %.4f  %+.2f%% (%+.2f..%+.2f)  bias %+.4f vs %+.4f  CSI8 %+.4f (%+.4f..%+.4f)  |FBI8-1| %+.3f (%+.3f..%+.3f)  rule=%s\n",
            row.scheme, row.product, row.candidate, row.role, row.candidate_RMSE,
            row.baseline_RMSE, 100 * row.rel_RMSE, 100 * row.rel_RMSE_lo, 100 * row.rel_RMSE_hi,
            row.candidate_bias, row.baseline_bias, row.delta_CSI_8, row.delta_CSI_8_lo,
            row.delta_CSI_8_hi, row.delta_FBI_8_dev, row.delta_FBI_8_dev_lo,
            row.delta_FBI_8_dev_hi, row.meets_rule)
    end
    println("\nVerdict (balanced_spatial, fixed before the run): ", verdict(rows))
    println("Wrote $(relpath(OUTDIR, ROOT))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_main()
end

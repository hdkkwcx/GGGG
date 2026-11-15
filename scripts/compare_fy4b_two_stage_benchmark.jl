#!/usr/bin/env julia

# Does FY4B's two-stage series help after MGWR and the agreement blend? Scores FY4B_TS_E (the chance
# of rain times the mean wet amount of each hour's context class, `SatelliteTwoStage`) against FY4B_QM
# and plain FY4B in one benchmark run. Fits nothing: it reads the runs' output files.
#
#   julia --project=. scripts/compare_fy4b_two_stage_benchmark.jl          # the full runs
#   julia --project=. scripts/compare_fy4b_two_stage_benchmark.jl smoke    # the smoke runs, to test this
#
# The runs:
# - the paper's `benchmark` step (base):
#     julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates \
#       --satellite-wet-blend --blend-axis agreement_envelope --fused-anchor --with-fy4b \
#       --fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv --fy4b-qm --fy4b-qm-adaptive
# - the same with --fy4b-twostage, which appends FY4B_TS_E and leaves the rest alone.
# Both directories come from `benchmark_config`, so the names cannot drift.
#
# The screen (`scripts/screen_fy4b_two_stage.jl`) scored this series before any gauge correction:
# against FY4B_QM its RMSE is 13% lower, it never reaches 8 mm/h, and its gain in r goes once weeks
# are held out. Neither FY4B_TS_E nor this run's FY4B_QM is held out in time, so the pair is like for
# like, on hours whose weeks both have seen.
#
# Checks; the script errors if one fails:
# - the two-stage run's grid, stations and folds equal the base run's;
# - every product of the base run keeps its `raw`, `adw`, `mgwr` and both blends, and its evaluation
#   mask, cell for cell, NaNs included: the new product moved nothing else. (`residual_gwr` is left
#   out, as threading can move it by a sub-ULP between runs; see CLAUDE.md.) Skipped in smoke mode
#   when there is no base smoke run.
#
# Each row of paired.csv is FY4B_TS_E (candidate) against FY4B_QM or plain FY4B (baseline), same field,
# on the cells both score, with intervals over whole 7-day blocks (`block_times`) on RMSE, r, bias,
# |bias|, CSI at 0.1, 2 and 8 mm/h, FBI_8 and |FBI_8 − 1|. Fields: `raw`, `mgwr`, the ADW blend and the
# agreement blend. Also reported, in no verdict: each FY4B series' agreement blend against GPM_uncal's.
#
# A product's evaluation mask needs every method finite, so a fold in which one method fails for
# FY4B_TS_E leaves the cells both score (in the full run, `mixed_gwr` fails balanced fold 2). The
# rows with role `two_stage_qm_mask` therefore repeat the pair against FY4B_QM on FY4B_QM's own mask,
# wherever the compared field of FY4B_TS_E is finite. They are reported only.
#
# The rule, fixed before the full run, is the decision rule of the dry-band, fusion and humidity-gate
# tests, on the agreement blend in balanced spatial CV, FY4B_TS_E against FY4B_QM. The two-stage
# series helps only if:
# 1. RMSE changes by −`MIN_REL_GAIN` or better, with its interval below 0;
# 2. CSI at 8 mm/h is not significantly lower (upper bound ≥ 0);
# 3. |FBI_8 − 1| is not significantly worse (lower bound ≤ 0);
# 4. |bias| grows by at most `BIAS_TOLERANCE` mm/h.
# Otherwise FY4B_QM stays FY4B's quality step. The random split and the pair against plain FY4B are
# reported only.
#
# Writes output/fy4b_two_stage_benchmark_comparison/ (`_smoke` for the smoke runs): paired.csv,
# rule.csv and provenance.csv.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("RunProvenance", "TableIO", "DayBlockBootstrap", "TimeBlockFolds")

using CSV, Printf, Statistics
using Main.DayBlockBootstrap: read_wide_values, read_wide_mask, day_groups, group_sums,
    compare_days
using Main.TimeBlockFolds: block_times
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const FY4B_FILE = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024", "fy4b_w3_h0.csv")
const THRESHOLDS = [0.1, 2.0, 8.0]
const CI_STATS = (:RMSE, :r, :bias, :abs_bias, Symbol("CSI_0.1"), :CSI_2, :CSI_8, :FBI_8, :FBI_8_dev)
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261101
const BLOCK_DAYS = 7
const MIN_REL_GAIN = 0.003
const BIAS_TOLERANCE = 0.005
const AGREEMENT = _blend_method("mgwr", :agreement_envelope)
const SCORED = ["raw", "mgwr", _blend_method("mgwr"), AGREEMENT]
const GUARDED = ["raw", "adw", "mgwr", _blend_method("mgwr"), AGREEMENT]
const BASE_PRODUCTS = ["FY4B", "GPM_uncal", "GSMaP", "MERGED_MEAN", "MERGED_OLS", FY4B_QM_PRODUCT]

"""
The two runs' configurations, and a copy of the two-stage run's `MGERConfig` that writes to `outdir`:
loading the data writes `global_common_time_qc.csv`, which must not land in either run.
"""
function run_configs(mode, outdir)
    common = (; nested_covariates=true, satellite_wet_blend=true, fused_anchor=true, gpm_uncal=true,
        blend_axes=[:constant, :agreement_envelope], fy4b=true, fy4b_file=FY4B_FILE, fy4b_qm=true,
        fy4b_qm_adaptive=true)
    two_stage = benchmark_config(mode; common..., fy4b_twostage=true)
    base = benchmark_config(mode; common...)
    fields = fieldnames(MGERConfig)
    mger = MGERConfig(; merge(
        NamedTuple{fields}(Tuple(getfield(two_stage.mger, field) for field in fields)), (; outdir))...)
    return two_stage, base, mger
end

"""Error unless every base product keeps its guarded fields and mask in the two-stage run."""
function check_base_unchanged(dirs, schemes, ids, times)
    path(run, scheme, product, file) = joinpath(dirs[run], scheme, lowercase(product), file)
    grid(run) = CSV.read(joinpath(dirs[run], "global_common_time_qc.csv"), DataFrame)
    isequal(grid(:two_stage), grid(:base)) || error("the runs have different time grids")
    for scheme in schemes
        folds(run) = CSV.read(joinpath(dirs[run], scheme, "split_common.csv"), DataFrame)
        isequal(folds(:two_stage), folds(:base)) || error("$scheme: the folds differ")
        for product in BASE_PRODUCTS
            mask(run) = read_wide_mask(path(run, scheme, product, "common_evaluation_mask.csv"), ids, times)
            mask(:two_stage) == mask(:base) || error("$scheme $product: the evaluation mask moved")
            for method in GUARDED
                field(run) = read_wide_values(path(run, scheme, product, "oof_$method.csv"), ids, times)
                isequal(field(:two_stage), field(:base)) ||
                    error("$scheme $product $method moved, but the two-stage product never reaches it")
            end
        end
    end
    println("  every base product keeps $(join(GUARDED, ", ")) and its mask, cell for cell")
end

function compare_main(args=ARGS)
    mode = isempty(args) ? :full : Symbol(only(args))
    outdir = joinpath(ROOT, "output",
        "fy4b_two_stage_benchmark_comparison" * (mode == :full ? "" : "_$mode"))
    mkpath(outdir)
    two_stage, base, mger = run_configs(mode, outdir)
    dirs = (; two_stage=two_stage.mger.outdir, base=base.mger.outdir)
    isfile(joinpath(dirs.two_stage, "metrics_pooled.csv")) ||
        error("the two-stage run is not finished: $(dirs.two_stage)")
    has_base = isfile(joinpath(dirs.base, "metrics_pooled.csv"))
    has_base || mode != :full || error("the base run is missing: $(dirs.base)")
    schemes = string.(two_stage.cv_schemes)
    path(scheme, product, file) = joinpath(dirs.two_stage, scheme, lowercase(product), file)
    products = [FY4B_TS_PRODUCT, FY4B_QM_PRODUCT, "FY4B", "GPM_uncal"]
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path, collect(values(mger.sat_paths)),
        joinpath(dirs.two_stage, "fy4b_twostage_tables.csv"),
        [path(s, p, file) for s in schemes for p in products
            for file in vcat(["oof_$m.csv" for m in SCORED], "common_evaluation_mask.csv")])
    provenance = start_provenance(ROOT, outdir; inputs)

    _, ids, data = load_global_common_product_data(mger)
    times = data["GPM_uncal"].times
    y = data["GPM_uncal"].Y_obs
    @printf("%d stations x %d hours\n", length(ids), length(times))
    println("Checks:")
    has_base ? check_base_unchanged(dirs, schemes, ids, times) :
        println("  no base smoke run, so the base guard is skipped")

    week_group, n_week = day_groups(block_times(times; block_days=BLOCK_DAYS), length(ids))
    week_sums(x, mask) = group_sums(x, y, mask, week_group, n_week, THRESHOLDS)
    compare(candidate, baseline, mask) = compare_days(week_sums(candidate, mask), week_sums(baseline, mask),
        THRESHOLDS; ci_stats=CI_STATS, reps=BOOTSTRAP_REPS, seed=BOOTSTRAP_SEED)
    rows = NamedTuple[]
    for scheme in schemes
        mask(product) = read_wide_mask(path(scheme, product, "common_evaluation_mask.csv"), ids, times)
        field(product, method) = read_wide_values(path(scheme, product, "oof_$method.csv"), ids, times)
        masks = Dict(p => mask(p) for p in products)
        shared_cells(a, b, x, z) = masks[a] .& masks[b] .& .!isnan.(y) .& .!isnan.(x) .& .!isnan.(z)
        for (baseline_product, role) in ((FY4B_QM_PRODUCT, "two_stage"), ("FY4B", "plain")),
            method in SCORED

            candidate, baseline = field(FY4B_TS_PRODUCT, method), field(baseline_product, method)
            row = compare(candidate, baseline,
                shared_cells(FY4B_TS_PRODUCT, baseline_product, candidate, baseline))
            meets_rule = row.rel_RMSE <= -MIN_REL_GAIN && row.rel_RMSE_hi < 0 &&
                row.delta_CSI_8_hi >= 0 && row.delta_FBI_8_dev_lo <= 0 &&
                row.delta_abs_bias <= BIAS_TOLERANCE
            push!(rows, merge((; scheme, candidate=FY4B_TS_PRODUCT, baseline=baseline_product,
                method, role), row, (; meets_rule)))
        end
        for method in SCORED
            candidate, baseline = field(FY4B_TS_PRODUCT, method), field(FY4B_QM_PRODUCT, method)
            push!(rows, merge((; scheme, candidate=FY4B_TS_PRODUCT, baseline=FY4B_QM_PRODUCT, method,
                role="two_stage_qm_mask"), compare(candidate, baseline,
                    shared_cells(FY4B_QM_PRODUCT, FY4B_QM_PRODUCT, candidate, baseline)),
                (; meets_rule=false)))
        end
        gpm = field("GPM_uncal", AGREEMENT)
        for product in (FY4B_TS_PRODUCT, FY4B_QM_PRODUCT, "FY4B")
            candidate = field(product, AGREEMENT)
            push!(rows, merge((; scheme, candidate=product, baseline="GPM_uncal", method=AGREEMENT,
                role="gap"), compare(candidate, gpm, shared_cells(product, "GPM_uncal", candidate, gpm)),
                (; meets_rule=false)))
        end
    end

    judged = only(filter(r -> r.scheme == "balanced_spatial" && r.role == "two_stage" &&
        r.method == AGREEMENT, rows))
    rule = [
        (; check="RMSE -$(100MIN_REL_GAIN)% or better, interval below 0", value=judged.rel_RMSE,
            lo=judged.rel_RMSE_lo, hi=judged.rel_RMSE_hi,
            pass=judged.rel_RMSE <= -MIN_REL_GAIN && judged.rel_RMSE_hi < 0),
        (; check="CSI_8 not significantly lower", value=judged.delta_CSI_8, lo=judged.delta_CSI_8_lo,
            hi=judged.delta_CSI_8_hi, pass=judged.delta_CSI_8_hi >= 0),
        (; check="|FBI_8 - 1| not significantly worse", value=judged.delta_FBI_8_dev,
            lo=judged.delta_FBI_8_dev_lo, hi=judged.delta_FBI_8_dev_hi, pass=judged.delta_FBI_8_dev_lo <= 0),
        (; check="|bias| up by at most $(BIAS_TOLERANCE) mm/h", value=judged.delta_abs_bias,
            lo=judged.delta_abs_bias_lo, hi=judged.delta_abs_bias_hi,
            pass=judged.delta_abs_bias <= BIAS_TOLERANCE),
    ]
    write_csv_atomic(joinpath(outdir, "paired.csv"), DataFrame(rows))
    write_csv_atomic(joinpath(outdir, "rule.csv"), DataFrame(rule))
    CSV.write(joinpath(outdir, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println("\nFY4B_TS_E vs FY4B_QM and plain FY4B, and each against GPM_uncal's agreement blend ",
        "(95% intervals over $(BLOCK_DAYS)-day blocks):")
    for row in rows
        @printf("  %-17s %-17s %-9s vs %-9s %-18s n %d  RMSE %.4f vs %.4f %+.2f%% (%+.2f..%+.2f)  r %+.4f (%+.4f..%+.4f)  bias %+.4f vs %+.4f  CSI0.1 %+.4f  CSI2 %+.4f  CSI8 %+.4f (%+.4f..%+.4f)  |FBI8-1| %+.3f (%+.3f..%+.3f)\n",
            row.scheme, row.role, row.candidate, row.baseline, row.method, row.candidate_n, row.candidate_RMSE,
            row.baseline_RMSE, 100 * row.rel_RMSE, 100 * row.rel_RMSE_lo, 100 * row.rel_RMSE_hi,
            row.delta_r, row.delta_r_lo, row.delta_r_hi, row.candidate_bias, row.baseline_bias,
            getproperty(row, Symbol("delta_CSI_0.1")), row.delta_CSI_2, row.delta_CSI_8,
            row.delta_CSI_8_lo, row.delta_CSI_8_hi, row.delta_FBI_8_dev, row.delta_FBI_8_dev_lo,
            row.delta_FBI_8_dev_hi)
    end
    println("\nRule (balanced spatial, agreement blend, FY4B_TS_E against FY4B_QM):")
    for row in rule
        @printf("  %-4s %s: %+.4f (%+.4f..%+.4f)\n", row.pass ? "pass" : "FAIL", row.check, row.value, row.lo, row.hi)
    end
    println(all(row.pass for row in rule) ?
        "Verdict: the two-stage series HELPS: FY4B_TS_E meets the rule against FY4B_QM" :
        "Verdict: NO: FY4B_TS_E does not meet the rule; FY4B_QM stays the quality step")
    println("Wrote $(relpath(outdir, ROOT))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_main()
end

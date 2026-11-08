#!/usr/bin/env julia

# How far does FY4B's quality rise? Score the corrected FY4B (`FY4B_QM`) against plain FY4B, and
# against GPM_uncal and GSMaP, from one benchmark run's out-of-fold files. Fits nothing.
#
#   julia --project=. scripts/compare_fy4b_quality.jl
#   julia --project=. scripts/compare_fy4b_quality.jl --block-days 7
#
# `--block-days N` resamples whole N-day blocks instead of days and writes to
# output/fy4b_quality_block{N}d/. `--block-days 1` must reproduce the default file; without the flag
# nothing changes (`TimeBlockFolds.block_times`).
#
# The run (2026-10-31, the paper's `fy4b_benchmark` step):
#     julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates \
#       --satellite-wet-blend --blend-axis agreement_envelope --fused-anchor --with-fy4b \
#       --fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv --fy4b-qm --fy4b-qm-adaptive \
#       --fy4b-qm-fusion
# Its directory comes from `benchmark_config`, so the name cannot drift.
#
# Each row of paired.csv is `product`'s field against `baseline`'s, same method, on the cells both
# score: both evaluation masks, a gauge value, and both fields finite. Pairs:
# - FY4B_QM vs FY4B: what the correction gains;
# - FY4B_QM vs GPM_uncal and vs GSMaP: how far the corrected FY4B still is from them.
# Intervals are paired day-block bootstraps (`DayBlockBootstrap`) on RMSE, r, bias and CSI at 2 and
# 8 mm/h.
#
# Writes output/fy4b_quality/: paired.csv, provenance.csv.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("RunProvenance", "TableIO", "DayBlockBootstrap", "TimeBlockFolds")

using CSV, Printf
using Main.DayBlockBootstrap: read_wide_values, read_wide_mask, day_groups, group_sums,
    compare_days
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic
using Main.TimeBlockFolds: block_times, block_days_option

const FY4B_FILE = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024", "fy4b_w3_h0.csv")
const OUTDIR = joinpath(ROOT, "output", "fy4b_quality")
const THRESHOLDS = [0.1, 2.0, 8.0]
const CI_STATS = (:RMSE, :r, :bias, :CSI_2, :CSI_8)
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261101
const CORRECTED = "FY4B_QM"
const BASELINES = ["FY4B", "GPM_uncal", "GSMaP"]
const METHODS = ["raw", "mgwr", _blend_method("mgwr", :agreement_envelope)]

"""
The run's configuration, and a copy of its `MGERConfig` that writes to `outdir`: loading the grid
writes `global_common_time_qc.csv`, which must not land in the run.
"""
function run_config(outdir=OUTDIR)
    config = benchmark_config(:full; nested_covariates=true, satellite_wet_blend=true,
        fused_anchor=true, gpm_uncal=true, blend_axes=[:constant, :agreement_envelope],
        fy4b=true, fy4b_file=FY4B_FILE, fy4b_qm=true, fy4b_qm_adaptive=true, fy4b_qm_fusion=true)
    fields = fieldnames(MGERConfig)
    mger = MGERConfig(; merge(
        NamedTuple{fields}(Tuple(getfield(config.mger, field) for field in fields)),
        (; outdir))...)
    return config, mger
end

function quality_main(args=ARGS)
    block_days, _ = block_days_option(args)
    outdir = block_days === nothing ? OUTDIR : OUTDIR * "_block$(block_days)d"
    mkpath(outdir)
    config, mger = run_config(outdir)
    dir = config.mger.outdir
    isfile(joinpath(dir, "metrics_pooled.csv")) || error("the run is not finished: $dir")
    schemes = string.(config.cv_schemes)
    products = vcat(CORRECTED, BASELINES)
    path(scheme, product, file) = joinpath(dir, scheme, lowercase(product), file)
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path,
        collect(values(mger.sat_paths)),
        [joinpath(dir, scheme, "split_common.csv") for scheme in schemes],
        [path(s, p, file) for s in schemes for p in products
            for file in vcat(["oof_$m.csv" for m in METHODS], "common_evaluation_mask.csv")])
    provenance = start_provenance(ROOT, outdir; inputs)

    _, ids, data = load_global_common_product_data(mger)
    times = data["GPM_uncal"].times
    y = data["GPM_uncal"].Y_obs
    @printf("Run: %s\n  %d stations x %d hours\n", relpath(dir, ROOT), length(ids), length(times))

    day_group, n_day = day_groups(block_days === nothing ? times : block_times(times; block_days),
        length(ids))
    day_sums(x, mask) = group_sums(x, y, mask, day_group, n_day, THRESHOLDS)
    rows = NamedTuple[]
    for scheme in schemes
        masks = Dict(p => read_wide_mask(path(scheme, p, "common_evaluation_mask.csv"), ids, times)
            for p in products)
        for baseline in BASELINES, method in METHODS
            candidate = read_wide_values(path(scheme, CORRECTED, "oof_$method.csv"), ids, times)
            reference = read_wide_values(path(scheme, baseline, "oof_$method.csv"), ids, times)
            shared = masks[CORRECTED] .& masks[baseline] .& .!isnan.(y) .&
                .!isnan.(candidate) .& .!isnan.(reference)
            push!(rows, merge((; scheme, product=CORRECTED, baseline, method),
                compare_days(day_sums(candidate, shared), day_sums(reference, shared), THRESHOLDS;
                    ci_stats=CI_STATS, reps=BOOTSTRAP_REPS, seed=BOOTSTRAP_SEED)))
        end
    end

    write_csv_atomic(joinpath(outdir, "paired.csv"), DataFrame(rows))
    CSV.write(joinpath(outdir, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println("\n$CORRECTED against each baseline, same method, shared cells (95% day-block interval):")
    for row in rows
        @printf("%-17s vs %-10s %-18s n=%8d  RMSE %.4f vs %.4f  %+.2f%% (%+.2f..%+.2f)  r %.3f vs %.3f  CSI2 %+.4f (%+.4f..%+.4f)\n",
            row.scheme, row.baseline, row.method, row.candidate_n, row.candidate_RMSE,
            row.baseline_RMSE, 100 * row.rel_RMSE, 100 * row.rel_RMSE_lo, 100 * row.rel_RMSE_hi,
            row.candidate_r, row.baseline_r, row.delta_CSI_2, row.delta_CSI_2_lo,
            row.delta_CSI_2_hi)
    end
    println("Wrote $(relpath(outdir, ROOT))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    quality_main()
end

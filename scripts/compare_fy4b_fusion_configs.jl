#!/usr/bin/env julia

# The final run's two configurations, scored on one sample. Both feed GPM (IMERG Uncal), GSMaP and a
# FY4B to every merging method (MERGED_MEAN, MERGED_OLS, the agreement bands):
# - `plain`: plain FY4B. The paper's `benchmark` step:
#     julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates \
#       --satellite-wet-blend --blend-axis agreement_envelope --fused-anchor --with-fy4b \
#       --fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv --fy4b-qm --fy4b-qm-adaptive
# - `qm`: each fold's FY4B_QM in its place. The same command plus --fy4b-qm-fusion.
# Fits nothing: it reads the two runs' output files. Both directories come from `benchmark_config`,
# so the names cannot drift.
#
#   julia --project=. scripts/compare_fy4b_fusion_configs.jl          # the full runs
#   julia --project=. scripts/compare_fy4b_fusion_configs.jl smoke    # the smoke runs, to test this
#
# What is equal in the two runs, and across the three products within each: the product list (so
# every product's selection seed, `run_seed` in `InterpolationBenchmarkFold.jl`), the hours (those
# at which the gauges and all three products have a value, all from June 2022), the folds, and the
# tuning hours, which are drawn from the gauges. Everything fitted to gauges (FY4B_QM's shift, map
# and strength, the merge weights, the blend weights, the covariate selection) is fitted inside each
# fold on its training gauges. Only FY4B is mapped to the gauges; GPM and GSMaP stay as delivered.
#
# The sample. A run scores each product on its own mask. Here every configuration, product and
# method is scored on the same held-out gauge-hours: those in the matched sample of both runs
# (`BenchmarkDiagnostics.matched_evaluation_mask`, the cells at which every product's prediction by
# every `MATCHED_METHODS` method exists). Scores are pooled over the whole record, June 2022 –
# December 2024; there are no rows per year, month, season or FY4B period.
#
# Checks; the script errors if one fails: both runs' time grids (`global_common_time_qc.csv`) equal,
# split_common.csv equal per scheme, no method failed a whole fold in either run, and the guards
# (below) identical in the two runs.
#
# Writes output/fy4b_fusion_configs/ (`_smoke` for the smoke runs):
# - scores.csv: one row per scheme, configuration, product and method. `n` is the shared sample, the
#   same in every row of a scheme. Coverage beside the scores: `n_eligible` is the held-out
#   gauge-hours with a gauge value on the run's hours, `n_own` those at which this method's own
#   prediction exists, `own_coverage` their ratio, and `matched_share` is n / n_eligible.
# - paired.csv: `qm` (candidate) against `plain` (baseline), field by field on the shared sample,
#   with 95% paired intervals over whole 7-day blocks. `role` says where the FY4B in the fusion
#   reaches the field:
#   - `guard`: every field of GPM_uncal, GSMaP, FY4B and FY4B_QM but the agreement blends. The
#     fusion never reaches them, so they must be identical, NaNs included, or the script errors;
#   - `control`: the gauge-only fields (idw, adw, tps, gwr) of MERGED_MEAN and MERGED_OLS;
#   - `treatment`: every agreement blend (FY4B in the bands), and every other field of MERGED_MEAN
#     and MERGED_OLS (FY4B in the merged anchor).
#   It is a description. It applies no rule and gives no verdict: this pair was judged as try 13 of
#   `Projects_Pending_Handover/fy4b_tries.md`, and this is not another try.
# - matched_sample.csv: per configuration, product and fold, the cells in the product's own mask and
#   in the shared sample.
# - selection.csv: mgwr's covariates and roles per fold in both runs. A merged product's covariates
#   are selected on its own anchor's residuals, so they can differ between the configurations.
# - fold_status.csv (each run's joint-method folds that did not succeed), provenance.csv.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("RunProvenance", "TableIO", "DayBlockBootstrap", "TimeBlockFolds",
    "BenchmarkDiagnostics")

using CSV, Printf, Statistics
using Main.DayBlockBootstrap: read_wide_values, day_groups, group_sums, compare_days, sum_stats
using Main.BenchmarkDiagnostics: MATCHED_METHODS, matched_evaluation_mask, matched_fold_failures,
    matched_sample_table, read_fold_map
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic
using Main.TimeBlockFolds: block_times

const FY4B_FILE = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024", "fy4b_w3_h0.csv")
const THRESHOLDS = [0.1, 2.0, 8.0]
const CI_STATS = (:RMSE, :bias, :abs_bias, :r, :CSI_2, :CSI_8, :FBI_8, :FBI_8_dev)
"""The intervals resample whole blocks of this many days, 2,000 times, as the paper's do."""
const BLOCK_DAYS = 7
const BOOTSTRAP_REPS = 2000
const BOOTSTRAP_SEED = 20261111
const SINGLE = ["GPM_uncal", "GSMaP", "FY4B", "FY4B_QM"]
const MERGED = ["MERGED_MEAN", "MERGED_OLS"]
const GAUGE_ONLY = ["idw", "adw", "tps", "gwr"]
const AGREEMENT_PREFIX = "blend_agrenv_"
const JOINT_METHODS = ["residual_gwr", "mixed_gwr", "mgwr"]
const CONFIGS = (:plain, :qm)

"""Where the FY4B in the fusion reaches `method`'s field for `product`: see the header."""
function field_role(product, method)
    startswith(method, AGREEMENT_PREFIX) && return "treatment"
    product in MERGED || return "guard"
    return method in GAUGE_ONLY ? "control" : "treatment"
end

"""
The two runs' configurations, and a copy of the first one's `MGERConfig` that writes to `outdir`:
loading the data writes `global_common_time_qc.csv`, which must not land in either run.
"""
function run_configs(mode, outdir)
    common = (; nested_covariates=true, satellite_wet_blend=true, fused_anchor=true, gpm_uncal=true,
        blend_axes=[:constant, :agreement_envelope], fy4b=true, fy4b_file=FY4B_FILE, fy4b_qm=true,
        fy4b_qm_adaptive=true)
    plain = benchmark_config(mode; common...)
    qm = benchmark_config(mode; common..., fy4b_qm_fusion=true)
    fields = fieldnames(MGERConfig)
    mger = MGERConfig(; merge(
        NamedTuple{fields}(Tuple(getfield(plain.mger, field) for field in fields)), (; outdir))...)
    return (; plain, qm), mger
end

"""Each run's joint-method folds that did not succeed."""
function fold_status(dirs)
    rows = NamedTuple[]
    for config in keys(dirs),
            row in eachrow(CSV.read(joinpath(dirs[config], "run_status.csv"), DataFrame))
        row.method in JOINT_METHODS && row.status != "success" && push!(rows, (;
            config=String(config), scheme=row.scheme, product=row.product, fold=row.fold,
            method=row.method, status=row.status))
    end
    return DataFrame(rows)
end

"""`mgwr`'s selected covariates and roles per scheme, product and fold in both runs."""
function fold_selection(dirs)
    read_roles(config) = Dict((row.scheme, row.product, row.fold) =>
            coalesce(row.covariate_selected_roles, "")
        for row in eachrow(CSV.read(joinpath(dirs[config], "run_status.csv"), DataFrame;
            types=Dict(:covariate_selected_roles => String)))
        if row.method == "mgwr")
    plain, qm = read_roles(:plain), read_roles(:qm)
    keys(plain) == keys(qm) || error("the runs score different folds")
    return DataFrame([(; scheme, product, fold, plain=plain[(scheme, product, fold)],
            qm=qm[(scheme, product, fold)],
            same=Set(split(plain[(scheme, product, fold)], ";")) ==
                Set(split(qm[(scheme, product, fold)], ";")))
        for (scheme, product, fold) in sort(collect(keys(plain)))])
end

function configs_main(args=ARGS)
    positional = filter(!startswith("--"), args)
    mode = isempty(positional) ? :full : Symbol(only(positional))
    outdir = joinpath(ROOT, "output", "fy4b_fusion_configs" * (mode == :full ? "" : "_$mode"))
    mkpath(outdir)
    configs, mger = run_configs(mode, outdir)
    dirs = (; plain=configs.plain.mger.outdir, qm=configs.qm.mger.outdir)
    println("plain: $(relpath(dirs.plain, ROOT))\nqm:    $(relpath(dirs.qm, ROOT))")
    for dir in dirs
        isfile(joinpath(dir, "metrics_pooled.csv")) || error("the run is not finished: $dir")
    end
    schemes = string.(configs.plain.cv_schemes)
    string.(configs.qm.cv_schemes) == schemes || error("the runs have different schemes")
    products = vcat(SINGLE, MERGED)
    path(config, scheme, product, file) = joinpath(dirs[config], scheme, lowercase(product), file)
    # A method the run did not write (none in the paper's configuration) is skipped in both.
    methods = filter(method -> isfile(path(:plain, first(schemes), first(products), "oof_$method.csv")),
        MATCHED_METHODS)
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path, collect(values(mger.sat_paths)),
        [joinpath(dirs[config], file) for config in CONFIGS
            for file in ("global_common_time_qc.csv", "run_status.csv")],
        [joinpath(dirs[config], scheme, "split_common.csv") for config in CONFIGS
            for scheme in schemes],
        [path(config, scheme, product, file) for config in CONFIGS for scheme in schemes
            for product in products
            for file in vcat(["oof_$method.csv" for method in methods], "common_evaluation_mask.csv")])
    provenance = start_provenance(ROOT, outdir; inputs)

    println("Checking that the two runs line up:")
    grid(config) = CSV.read(joinpath(dirs[config], "global_common_time_qc.csv"), DataFrame)
    isequal(grid(:plain), grid(:qm)) || error("the runs have different time grids")
    _, ids, data = load_global_common_product_data(mger)
    times = data[first(SINGLE)].times
    y = data[first(SINGLE)].Y_obs
    @printf("  %d stations, %d hours on both runs' grid, %s to %s\n", length(ids), length(times),
        minimum(times), maximum(times))
    statuses = Dict(config => CSV.read(joinpath(dirs[config], "run_status.csv"), DataFrame)
        for config in CONFIGS)

    block_group, n_block = day_groups(block_times(times; block_days=BLOCK_DAYS), length(ids))
    eligible = .!isnan.(y)
    n_eligible = count(eligible)
    score_rows = NamedTuple[]
    paired_rows = NamedTuple[]
    sample_tables = DataFrame[]
    for scheme in schemes
        folds(config) = CSV.read(joinpath(dirs[config], scheme, "split_common.csv"), DataFrame)
        isequal(folds(:plain), folds(:qm)) || error("$scheme: the folds differ")
        fold_map = read_fold_map(joinpath(dirs.plain, scheme, "split_common.csv"))
        fold_of = [fold_map[id] for id in ids]
        shared = copy(eligible)
        samples = Dict{Symbol,Any}()
        for config in CONFIGS
            failures = matched_fold_failures(statuses[config]; scheme, products)
            isempty(failures) || error("$config, $scheme: these folds failed, so the shared " *
                "sample would lose them for every product:\n" *
                sprint(show, select(failures, :product, :fold, :method, :status)))
            samples[config] = matched_evaluation_mask(joinpath(dirs[config], scheme), products, ids,
                times)
            shared .&= samples[config].matched
        end
        for product in SINGLE
            samples[:plain].masks[product] == samples[:qm].masks[product] ||
                error("$scheme $product: the evaluation masks differ between the runs")
        end
        for config in CONFIGS
            table = matched_sample_table(samples[config].masks, shared, fold_of; scheme)
            insertcols!(table, 1, :config => String(config))
            table.n_matched_in_run = [count(view(samples[config].matched,
                    row.fold == "all" ? (1:length(ids)) : findall(==(parse(Int, row.fold)), fold_of), :))
                for row in eachrow(table)]
            push!(sample_tables, table)
        end
        n_shared = count(shared)
        @printf("  %-16s shared sample %d of %d eligible gauge-hours (%.2f%%), on %d hours; plain's own matched sample %d, qm's %d\n",
            scheme, n_shared, n_eligible, 100 * n_shared / n_eligible,
            count(vec(any(shared; dims=1))), count(samples[:plain].matched),
            count(samples[:qm].matched))

        block_sums(x) = group_sums(x, y, shared, block_group, n_block, THRESHOLDS)
        for product in products, method in methods
            fields = (;
                plain=read_wide_values(path(:plain, scheme, product, "oof_$method.csv"), ids, times),
                qm=read_wide_values(path(:qm, scheme, product, "oof_$method.csv"), ids, times))
            role = field_role(product, method)
            role == "guard" && !isequal(fields.plain, fields.qm) && error(
                "$scheme $product $method differs between the runs, but the fusion never reaches it")
            sums = (; plain=block_sums(fields.plain), qm=block_sums(fields.qm))
            for config in CONFIGS
                n_own = count(eligible .& .!isnan.(fields[config]))
                stats = sum_stats(vec(sum(sums[config]; dims=1)), THRESHOLDS)
                Int(stats.n) == n_shared || error("$scheme $product $method: not on the shared sample")
                push!(score_rows, merge((; scheme, config=String(config), product, method, role,
                    n_eligible, n_own, own_coverage=n_own / n_eligible,
                    matched_share=n_shared / n_eligible), stats))
            end
            push!(paired_rows, merge((; scheme, product, method, role),
                compare_days(sums.qm, sums.plain, THRESHOLDS;
                    ci_stats=CI_STATS, reps=BOOTSTRAP_REPS, seed=BOOTSTRAP_SEED)))
        end
    end
    println("  guards: every field the fusion never reaches is identical in both runs")

    fold_statuses = fold_status(dirs)
    selection = fold_selection(dirs)
    write_csv_atomic(joinpath(outdir, "scores.csv"), DataFrame(score_rows))
    write_csv_atomic(joinpath(outdir, "paired.csv"), DataFrame(paired_rows))
    write_csv_atomic(joinpath(outdir, "matched_sample.csv"), reduce(vcat, sample_tables))
    write_csv_atomic(joinpath(outdir, "fold_status.csv"), fold_statuses)
    write_csv_atomic(joinpath(outdir, "selection.csv"), selection)
    CSV.write(joinpath(outdir, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println("\nFY4B_QM in the fusion (qm) vs plain FY4B in it (plain), same field, shared sample ",
        "(95% $(BLOCK_DAYS)-day-block interval). Guards are identical and not listed:")
    for row in paired_rows
        row.role == "guard" && continue
        @printf("%-17s %-11s %-25s %-9s RMSE %.4f vs %.4f  %+.2f%% (%+.2f..%+.2f)  r %.4f vs %.4f  bias %+.4f vs %+.4f  CSI8 %+.4f (%+.4f..%+.4f)  |FBI8-1| %+.3f (%+.3f..%+.3f)\n",
            row.scheme, row.product, row.method, row.role, row.candidate_RMSE, row.baseline_RMSE,
            100 * row.rel_RMSE, 100 * row.rel_RMSE_lo, 100 * row.rel_RMSE_hi, row.candidate_r,
            row.baseline_r, row.candidate_bias, row.baseline_bias, row.delta_CSI_8,
            row.delta_CSI_8_lo, row.delta_CSI_8_hi, row.delta_FBI_8_dev, row.delta_FBI_8_dev_lo,
            row.delta_FBI_8_dev_hi)
    end
    println("\nLowest own prediction coverage per configuration (scores.csv has every row):")
    for config in CONFIGS, scheme in schemes
        rows = filter(row -> row.config == String(config) && row.scheme == scheme, score_rows)
        low = rows[argmin([row.own_coverage for row in rows])]
        @printf("  %-6s %-17s %.4f (%s %s)\n", config, scheme, low.own_coverage, low.product,
            low.method)
    end
    println("\nmgwr's selection per fold, plain | qm, where the runs differ (selection.csv has all):")
    for row in eachrow(selection)
        row.same || @printf("  %-17s %-11s fold %d  %s | %s\n", row.scheme, row.product, row.fold,
            row.plain, row.qm)
    end
    println("Wrote $(relpath(outdir, ROOT))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    configs_main()
end

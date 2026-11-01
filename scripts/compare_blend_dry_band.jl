#!/usr/bin/env julia

# Score the dry-band blends of a `--blend-dry-band` run against the agreement blend they extend,
# and check that the run moved nothing else. Fits nothing: it reads the two runs' output files.
#
#   julia --project=. scripts/compare_blend_dry_band.jl
#
# The run is the 2026-10-31 qmfusion run's command plus --blend-dry-band:
#
#   julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates \
#     --satellite-wet-blend --blend-axis agreement_envelope --blend-dry-band --fused-anchor \
#     --with-fy4b --fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv --fy4b-qm \
#     --fy4b-qm-adaptive --fy4b-qm-fusion
#
# The reference is the same command without --blend-dry-band. Both directories come from
# `benchmark_config`, so the names cannot drift.
#
# Checks; the script errors if one fails:
# - every row of the reference's metrics_pooled.csv, metrics_folds.csv and metrics_stratified.csv
#   is in the run, with the same n and the same scores up to `RTOL`, and the run adds only
#   dry-band rows. `RTOL` allows the sub-ULP drift the threaded prediction path can leave between
#   identical runs (see `InterpolationBenchmarkBootstrap.jl`);
# - blend_selection.csv: the reference's rows are unchanged, and each dry axis's weights for
#   bands 1-9 equal the agreement blend's.
#
# Writes output/blend_dry_band_comparison/:
# - paired.csv: per scheme and product, three comparisons, each on the product's own evaluation
#   mask with paired day-block intervals (`DayBlockBootstrap`) on RMSE, bias, |bias|, and CSI, FBI
#   and |FBI − 1| at 8 mm/h:
#   - `blend_agrenvdry_mgwr` and `blend_agrenvdryraw_mgwr` against `blend_agrenv_mgwr`;
#   - those two against each other;
#   - all three against `adw`.
#   `meets_rule` applies the blend screen's decision rule against `blend_agrenv_mgwr`, on
#   balanced_spatial: RMSE −`MIN_REL_GAIN` or better with an interval below 0, CSI and |FBI − 1|
#   at 8 mm/h not worse beyond their interval, and |bias| up by at most `BIAS_TOLERANCE` mm/h.
# - dry_weights.csv: the dry band's weight per scheme, product, fold and axis, `mgwr` source.
# - provenance.csv.

include(joinpath(@__DIR__, "run_interpolation_benchmark.jl"))
load_standalone_modules("RunProvenance", "TableIO", "DayBlockBootstrap")

using CSV, Printf, Statistics
using Main.DayBlockBootstrap: read_wide_values, read_wide_mask, day_groups, group_sums,
    compare_days
using Main.RunProvenance: start_provenance, finish_provenance!
using Main.TableIO: write_csv_atomic

const FY4B_FILE = joinpath(ROOT, "output", "fy4b_rebuild_2022_2024", "fy4b_w3_h0.csv")
const OUTDIR = joinpath(ROOT, "output", "blend_dry_band_comparison")
const THRESHOLDS = [0.1, 2.0, 8.0]
const CI_STATS = (:RMSE, :bias, :abs_bias, :CSI_8, :FBI_8, :FBI_8_dev)
const BOOTSTRAP_REPS = 1000
const BOOTSTRAP_SEED = 20261101
const MIN_REL_GAIN = 0.003
const BIAS_TOLERANCE = 0.005
const RTOL = 1e-10
const PRODUCTS = ["FY4B", FY4B_QM_PRODUCT, "GPM_uncal", "GSMaP", "MERGED_MEAN", "MERGED_OLS"]
const AGREEMENT = _blend_method("mgwr", :agreement_envelope)
const DRY = [_blend_method("mgwr", axis) for axis in BLEND_DRY_AXES]
"""Every method the dry axes add, for every source."""
const ADDED = Set(_blend_method(source, axis) for axis in BLEND_DRY_AXES
    for source in BLEND_SOURCE_METHODS)
"""(candidate, baseline) pairs in `paired.csv`."""
const PAIRS = vcat([(method, AGREEMENT) for method in DRY], [(DRY[2], DRY[1])],
    [(method, "adw") for method in vcat(AGREEMENT, DRY)])

"""
The reference run's configuration, the dry-band run's, and a copy of the latter's `MGERConfig`
that writes to `OUTDIR`: loading the grid writes `global_common_time_qc.csv`, which must not land
in either run.
"""
function run_configs()
    common = (; nested_covariates=true, satellite_wet_blend=true, fused_anchor=true,
        gpm_uncal=true, fy4b=true, fy4b_file=FY4B_FILE, fy4b_qm=true, fy4b_qm_adaptive=true,
        fy4b_qm_fusion=true)
    reference = benchmark_config(:full; common..., blend_axes=[:constant, :agreement_envelope])
    dry = benchmark_config(:full; common...,
        blend_axes=vcat([:constant, :agreement_envelope], collect(Symbol, BLEND_DRY_AXES)))
    fields = fieldnames(MGERConfig)
    settings = NamedTuple{fields}(Tuple(getfield(dry.mger, field) for field in fields))
    return reference, dry, MGERConfig(; merge(settings, (; outdir=OUTDIR))...)
end

_same(a, b) = isequal(a, b) || (a isa AbstractFloat && b isa AbstractFloat &&
    isapprox(a, b; rtol=RTOL))

"""
Error unless the run's `name` table holds every reference row unchanged, keyed on `keys`, with
`values` the same up to `RTOL`, and adds only rows whose method is in `ADDED`. Returns the
number of rows checked and the largest relative difference seen.
"""
function check_unchanged(name, reference_dir, run_dir, keys, values)
    before = CSV.read(joinpath(reference_dir, name), DataFrame)
    after = CSV.read(joinpath(run_dir, name), DataFrame)
    names(before) == names(after) || error("$name: the columns differ")
    key(row) = join((string(row[column]) for column in keys), "|")
    added = count(in(ADDED), after.method)
    nrow(after) - added == nrow(before) || error(
        "$name: $(nrow(before)) reference rows, but the run has $(nrow(after) - added) others")
    lookup = Dict(key(row) => row for row in eachrow(after) if !(row.method in ADDED))
    length(lookup) == nrow(before) || error("$name: the key is not unique")
    largest = 0.0
    for row in eachrow(before)
        other = get(lookup, key(row), nothing)
        other === nothing && error("$name: the run lacks the row $(key(row))")
        for column in values
            a, b = row[column], other[column]
            _same(a, b) || error("$name $(key(row)) $column: $a in the reference, $b in the run")
            if a isa AbstractFloat && isfinite(a) && a != 0
                largest = max(largest, abs(a - b) / abs(a))
            end
        end
    end
    @printf("  %s: %d reference rows unchanged (largest relative difference %.1e), %d added\n",
        name, nrow(before), largest, added)
end

"""
Error unless blend_selection.csv keeps the reference's rows, and each dry axis's weights for
bands 1-9 equal the agreement blend's. Returns the dry band's weights.
"""
function check_weights(reference_dir, run_dir)
    read_selection(dir) = CSV.read(joinpath(dir, "blend_selection.csv"), DataFrame;
        types=Dict(:lambdas => String))
    before = read_selection(reference_dir)
    after = read_selection(run_dir)
    key(row) = (row.scheme, row.product, row.fold, row.repeat, row.method)
    lookup = Dict(key(row) => row for row in eachrow(after))
    for row in eachrow(before), column in names(before)
        other = get(lookup, key(row), nothing)
        other === nothing && error("blend_selection.csv: the run lacks $(key(row))")
        _same(row[column], other[column]) || error(
            "blend_selection.csv $(key(row)) $column: $(row[column]) vs $(other[column])")
    end
    weights(row) = parse.(Float64, split(row.lambdas, "|"))
    rows = NamedTuple[]
    for row in eachrow(after)
        row.axis in string.(BLEND_DRY_AXES) || continue
        agreement = lookup[(row.scheme, row.product, row.fold, row.repeat,
            _blend_method(row.source, :agreement_envelope))]
        dry, banded = weights(row), weights(agreement)
        length(dry) == length(banded) + 1 && dry[1:end - 1] == banded || error(
            "blend_selection.csv $(key(row)): bands 1-9 differ from the agreement blend's")
        row.source == "mgwr" && push!(rows, (; row.scheme, row.product, row.fold, row.axis,
            row.method, dry_weight=dry[end], n_bands=row.n_bands))
    end
    println("  blend_selection.csv: $(nrow(before)) reference rows unchanged; " *
        "the dry axes' bands 1-9 repeat the agreement blend's weights")
    return rows
end

function compare_main()
    mkpath(OUTDIR)
    reference, dry, mger = run_configs()
    reference_dir, run_dir = reference.mger.outdir, dry.mger.outdir
    isfile(joinpath(run_dir, "metrics_pooled.csv")) || error("the run is not finished: $run_dir")
    schemes = string.(dry.cv_schemes)
    methods = unique(vcat("adw", AGREEMENT, DRY))
    oof(scheme, product, file) = joinpath(run_dir, scheme, lowercase(product), file)
    tables = ("metrics_pooled.csv", "metrics_folds.csv", "metrics_stratified.csv",
        "blend_selection.csv")
    inputs = vcat(mger.station_meta_path, mger.obs_hourly_wide_path,
        [mger.sat_paths[name] for name in sort(collect(keys(mger.sat_paths)))],
        [joinpath(dir, table) for dir in (reference_dir, run_dir) for table in tables],
        [oof(s, p, "oof_$m.csv") for s in schemes for p in PRODUCTS for m in methods],
        [oof(s, p, "common_evaluation_mask.csv") for s in schemes for p in PRODUCTS])
    provenance = start_provenance(ROOT, OUTDIR; inputs)

    println("Checking that the run moved nothing but the dry-band methods:")
    metric_keys = ["scheme", "product", "method", "fold", "repeat", "seed", "group", "level",
        "threshold"]
    metric_values = ["n", "coverage", "RMSE", "MAE", "Bias", "r", "POD", "FAR", "CSI"]
    for table in tables[1:3]
        check_unchanged(table, reference_dir, run_dir, metric_keys, metric_values)
    end
    weight_rows = check_weights(reference_dir, run_dir)

    _, ids, product_data = load_global_common_product_data(mger)
    times = product_data["FY4B"].times
    y = product_data["FY4B"].Y_obs
    day_group, n_day = day_groups(times, length(ids))
    day_sums(x, mask) = group_sums(x, y, mask, day_group, n_day, THRESHOLDS)
    rows = NamedTuple[]
    for scheme in schemes, product in PRODUCTS
        mask = read_wide_mask(oof(scheme, product, "common_evaluation_mask.csv"), ids, times) .&
            .!isnan.(y)
        series = Dict(m => read_wide_values(oof(scheme, product, "oof_$m.csv"), ids, times)
            for m in methods)
        for (candidate, baseline) in PAIRS
            shared = mask .& .!isnan.(series[candidate]) .& .!isnan.(series[baseline])
            row = merge((; scheme, product, candidate, baseline),
                compare_days(day_sums(series[candidate], shared),
                    day_sums(series[baseline], shared), THRESHOLDS;
                    ci_stats=CI_STATS, reps=BOOTSTRAP_REPS, seed=BOOTSTRAP_SEED))
            meets_rule = scheme == "balanced_spatial" && baseline == AGREEMENT &&
                row.rel_RMSE <= -MIN_REL_GAIN && row.rel_RMSE_hi < 0 &&
                row.delta_CSI_8_hi >= 0 && row.delta_FBI_8_dev_lo <= 0 &&
                row.delta_abs_bias <= BIAS_TOLERANCE
            push!(rows, merge(row, (; meets_rule)))
        end
        println("  $scheme $product scored")
    end

    write_csv_atomic(joinpath(OUTDIR, "paired.csv"), DataFrame(rows))
    write_csv_atomic(joinpath(OUTDIR, "dry_weights.csv"), DataFrame(weight_rows))
    CSV.write(joinpath(OUTDIR, "provenance.csv"), finish_provenance!(provenance, ROOT))

    println("\nRMSE change, 95% day-block interval:")
    for row in rows
        @printf("%-17s %-11s %-23s vs %-23s RMSE %.4f vs %.4f  %+.2f%% (%+.2f..%+.2f)  bias %+.4f  CSI8 %+.4f (%+.4f..%+.4f)  FBI8 %+.3f  rule=%s\n",
            row.scheme, row.product, row.candidate, row.baseline, row.candidate_RMSE,
            row.baseline_RMSE, 100 * row.rel_RMSE, 100 * row.rel_RMSE_lo,
            100 * row.rel_RMSE_hi, row.delta_bias, row.delta_CSI_8, row.delta_CSI_8_lo,
            row.delta_CSI_8_hi, row.delta_FBI_8, row.meets_rule)
    end
    println("Wrote $(relpath(OUTDIR, ROOT))")
end

if abspath(PROGRAM_FILE) == @__FILE__
    compare_main()
end

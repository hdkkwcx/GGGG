"""
Diagnostic fields describing whether a fold/method's prediction used a DEM or joint-covariate
model, and (when it used joint covariates) which variable roles were selected vs. actually used.
Factored out because the per-method fold loop in `run_interpolation_benchmark` needs the exact
same computation whether the method succeeded or raised — see the `try`/`catch` call sites below.
"""
function _run_status_role_fields(
    dem_context, joint_context, mode::String, method::String, output_method::String,
)
    uses_dem = dem_context !== nothing && mode == "residual" &&
        method in ("gwr", "mixed_gwr", "mgwr")
    uses_joint = joint_context !== nothing && mode == "residual" &&
        method in ("gwr", "mixed_gwr", "mgwr")
    selected_roles = uses_joint ? join([
        "$group=$(joint_context.roles[group])" for group in joint_context.variables
    ], ";") : ""
    effective_role_map = uses_joint ? joint_effective_roles(
        joint_context, output_method,
    ) : Dict{String,String}()
    effective_roles = uses_joint ? join([
        "$group=$(effective_role_map[group])" for group in joint_context.variables
    ], ";") : ""
    # `mixed_gwr` collapses onto `residual_gwr` whenever this fold's role map has no "global":
    # identical design, identical grid, identical numbers (see `joint_models_coincide`). Recorded
    # as a column so a reader of `run_status.csv` can see the two rows are one model, instead of
    # having to notice that their RMSEs match to the last digit.
    duplicate_of = uses_joint && output_method == "mixed_gwr" &&
        joint_models_coincide(joint_context, "mixed_gwr", "residual_gwr") ? "residual_gwr" : ""
    return (;
        uses_dem, uses_joint, selected_roles, effective_role_map, effective_roles, duplicate_of,
    )
end

"""Non-NaN share of the held-out cells a method could have predicted."""
_prediction_coverage(eligible, prediction) =
    count(eligible) == 0 ? 0.0 :
        count(eligible .& .!isnan.(prediction)) / count(eligible)

"""
`auto`'s `(status, error)` for one fold.

"skipped" rather than "failed" where `auto` was never applicable, so a genuine breakage stays
visible instead of being lost among expected non-runs.

A fold that did choose a contender is held to the same coverage gate as the seven fitted runs,
with the same reason string. Without that the `status` column meant two different things in one
table: `auto` read "success" at the very coverage that makes the model it selected read
"partial", so filtering `run_status.csv` on `status == "success"` kept the `auto` row and dropped
the joint-model rows describing the same fold at the same coverage. Measured on the 2026-09-02
full nested run, that was three rows — balanced_spatial/fold 4 for each product, at coverage
0.9028 against a 0.95 bar.
"""
function _auto_status(
    chose_contender::Bool, coverage::Float64, min_coverage::Float64,
    dem_enabled::Bool, has_selection_split::Bool, failures::Vector{String},
)
    chose_contender && return coverage >= min_coverage ? ("success", "") :
        ("partial", "prediction coverage below minimum")
    dem_enabled && return ("skipped", "auto is not run on the legacy DEM path")
    has_selection_split || return ("skipped", "fold has no inner selection split to choose on")
    return ("failed", "no auto contender scored: $(join(failures, " | "))")
end

"""
One row of `run_status.csv`.

The fold loop assembles this row in three places - a method that predicted, a method that raised,
and `auto` - which differ only in `status`, `error` and `prediction_coverage`. Passing `nothing`
for both contexts gives the all-"not_applicable" provenance `auto` reports: it fits no covariate
model of its own and delegates provenance to whichever method it chose, whose own row carries it.
Keeping that a `nothing` argument rather than a separate literal is what stops the three copies
drifting apart, which is also why `covariate_selection_mode` still reads "not_applicable" there
and `auto` stays out of `covariate_model_status.csv`.
"""
function _benchmark_status_row(
    dem_context, joint_context, joint_inputs, nested_joint::Bool;
    scheme, product, fold, repeat::Int, seed::Int, mode::String, method::String,
    output_method::String, status::String, error::String, prediction_coverage::Float64,
)
    (; uses_dem, uses_joint, selected_roles, effective_roles, duplicate_of) =
        _run_status_role_fields(dem_context, joint_context, mode, method, output_method)
    return (;
        scheme, product, fold, repeat, seed, method=output_method, status, error,
        prediction_coverage,
        dem_variable_count=uses_dem ? dem_context.dem_variable_count : 0,
        dem_selection_status=uses_dem ? dem_context.selection_status : "not_applicable",
        dem_variables=uses_dem ? join(dem_context.variables, ",") : "",
        dem_roles=uses_dem ? dem_context.dem_roles : "",
        covariate_selection_mode=uses_joint ?
            (nested_joint ? "nested_per_fold" : "fixed_full_data") : "not_applicable",
        covariate_variable_count=uses_joint ? length(joint_context.variables) : 0,
        covariate_variables=uses_joint ? join(joint_context.variables, ",") : "",
        covariate_selected_roles=selected_roles,
        covariate_effective_roles=effective_roles,
        covariate_spec_sha256=uses_joint ? something(joint_inputs.spec_sha256, "") : "",
        # Appended rather than grouped with the other status fields: existing readers of this
        # table index it by name, but the ad-hoc awk/pandas kind does not, and a new column at
        # the end cannot shift anything that already exists.
        duplicate_of,
    )
end

"""
Provenance and input-QC tables for a joint-covariate run.

`joint_spec_provenance.csv` is written on both paths so a reader can always tell which selection
mode produced the run: nested per-fold selection is confirmatory, a fixed full-data spec is not.
"""
function _write_joint_provenance(cfg::InterpolationBenchmarkConfig, joint_inputs, nested_joint::Bool)
    if nested_joint
        CSV.write(joinpath(cfg.mger.outdir, "joint_spec_provenance.csv"), DataFrame([(
            source_path="", sha256="", selection_mode="nested_per_fold", confirmatory=true,
        )]))
    else
        CSV.write(
            joinpath(cfg.mger.outdir, "joint_variable_spec_used.csv"),
            joint_inputs.specification.table,
        )
        CSV.write(joinpath(cfg.mger.outdir, "joint_spec_provenance.csv"), DataFrame([(
            source_path=abspath(something(cfg.joint_covariates).spec_path),
            sha256=joint_inputs.spec_sha256,
            selection_mode="fixed_full_data",
            confirmatory=false,
        )]))
    end
    CSV.write(joinpath(cfg.mger.outdir, "joint_era5_input_qc.csv"), joint_inputs.era5.qc)
    if joint_inputs.ndvi !== nothing
        CSV.write(joinpath(cfg.mger.outdir, "joint_ndvi_alignment_qc.csv"),
            joint_inputs.ndvi.aligned.qc)
    end
    return nothing
end

"""
Screen the DEM covariates once over every station, before any fold is drawn.

Reported alongside the per-fold screens so a reader can see which variables the full station set
supports; the `"full_data"` scheme label keeps it out of the cross-fold stability tables.
"""
function _screen_dem_full_data!(dem_store, cfg::InterpolationBenchmarkConfig, terrain, lonlat, products, product_data)
    dem = something(cfg.dem)
    for product in products
        data = product_data[product]
        selection = screen_dem_subset(
            terrain, lonlat, data.times, data.Y_obs, data.Y_sat, dem;
            scheme="full_data", product, fold=0, phase="full_data",
            seed=cfg.seed + 100_000 + Int(sum(codeunits(product))),
        )
        _store_dem_selection!(dem_store, selection)
    end
    return dem_store
end

"""
Every table a run accumulates, one field per output.

The row vectors become the CSVs `_write_benchmark_outputs` writes. `dem_store` and `joint_store`
hold the per-fold selection tables, `joint_scaling` and `joint_qc` each fold's scaling and
quality-control table. `selection_fallback_cells` names the cells whose fold was too small for an
inner selection split and fell back to leave-one-out; `benchmark_scope.csv` reports them, so the
fallback is never silent.
"""
_benchmark_rows() = (;
    metrics=NamedTuple[], scans=NamedTuple[], bootstrap=NamedTuple[], status=NamedTuple[],
    auto_selection=NamedTuple[], blend_selection=NamedTuple[], stack_selection=NamedTuple[],
    fused_anchor=NamedTuple[], fused_grouped=NamedTuple[],
    fy4b_shift=NamedTuple[], fy4b_map=NamedTuple[], fy4b_strength=NamedTuple[],
    fy4b_gate=NamedTuple[], fy4b_twostage=NamedTuple[], hurdle=NamedTuple[],
    selection_fallback_cells=String[],
    dem_store=_empty_dem_store(), joint_store=_empty_joint_store(),
    joint_scaling=DataFrame[], joint_qc=DataFrame[],
)

"""
Assemble every result table from the accumulated `rows` (`_benchmark_rows`) and write the run's
CSVs. `setup` is the run's `_benchmark_setup`.

Split out of `run_interpolation_benchmark` so the orchestration reads as
setup -> repeat/scheme/product/fold loops -> outputs.
"""
function _write_benchmark_outputs(setup, seeds, rows; provenance::DataFrame=_git_provenance())
    (; cfg, products, nested_joint, joint_inputs) = setup
    # One row per fold x band x term; `group` is the agreement-envelope band (0 = no product
    # wet). `fell_back` marks a band too sparse or singular to fit alone, which took the pooled fit.
    isempty(rows.fused_grouped) ||
        CSV.write(joinpath(cfg.mger.outdir, "fused_anchor_grouped_coefficients.csv"),
            DataFrame(rows.fused_grouped))
    # One row per fold x stack band: the band's inner-split cell count, whether it took the pooled
    # fit, and its weight on every `STACK_FEATURES` member.
    isempty(rows.stack_selection) ||
        CSV.write(joinpath(cfg.mger.outdir, "stack_selection.csv"), DataFrame(rows.stack_selection))
    # One row per fold x FY4B regime: the pixel shift chosen on the training stations, with the
    # pooled training r it reached and the unshifted r it is measured against.
    isempty(rows.fy4b_shift) ||
        CSV.write(joinpath(cfg.mger.outdir, "fy4b_qm_shift_selection.csv"),
            DataFrame(rows.fy4b_shift))
    # One row per fold x season x regime: the wet threshold and the map's pool and anchor knots.
    isempty(rows.fy4b_map) ||
        CSV.write(joinpath(cfg.mger.outdir, "fy4b_qm_maps.csv"), DataFrame(rows.fy4b_map))
    # One row per fold x regime x month: how strongly the map was applied, and the training error
    # with no adjustment, with the full map and with what was applied.
    isempty(rows.fy4b_strength) ||
        CSV.write(joinpath(cfg.mger.outdir, "fy4b_qm_strengths.csv"),
            DataFrame(rows.fy4b_strength))
    # One row per fold x time fold: the humidity cut (-Inf for no gate) chosen on the fold's training
    # stations over the other time folds' hours, with the training CSI before and after it.
    isempty(rows.fy4b_gate) ||
        CSV.write(joinpath(cfg.mger.outdir, "fy4b_qm_rh_gate_cuts.csv"), DataFrame(rows.fy4b_gate))
    # One row per fold x season-by-regime cell x context class: the training counts, the level they
    # were pooled at, the chance of rain, the mean wet amount and their product, the value applied.
    isempty(rows.fy4b_twostage) ||
        CSV.write(joinpath(cfg.mger.outdir, "fy4b_twostage_tables.csv"),
            DataFrame(rows.fy4b_twostage))

    metrics = DataFrame(rows.metrics)
    scans = DataFrame(rows.scans)
    bootstrap = DataFrame(rows.bootstrap)
    status = DataFrame(rows.status)
    claim = if :balanced_spatial in cfg.cv_schemes && cfg.bootstrap_reps > 0
        assess_gwr_claim(metrics, bootstrap, products)
    else
        DataFrame()
    end
    repeat_summary = summarize_repeats(metrics)
    fold_summary = summarize_fold_spread(metrics)
    rank_stability = method_rank_stability(metrics)
    CSV.write(joinpath(cfg.mger.outdir, "metrics_stratified.csv"), metrics)
    CSV.write(joinpath(cfg.mger.outdir, "metrics_folds.csv"), filter(:fold => (x -> !ismissing(x)), metrics))
    CSV.write(joinpath(cfg.mger.outdir, "metrics_pooled.csv"), filter(:fold => ismissing, metrics))
    nrow(repeat_summary) > 0 &&
        CSV.write(joinpath(cfg.mger.outdir, "metrics_repeat_summary.csv"), repeat_summary)
    # Dispersion beside the pooled point estimate, so a headline gap can be read against how much
    # the methods move between folds. See `summarize_fold_spread` for why `_std` is not a standard
    # error.
    nrow(fold_summary) > 0 &&
        CSV.write(joinpath(cfg.mger.outdir, "metrics_fold_summary.csv"), fold_summary)
    nrow(rank_stability) > 0 &&
        CSV.write(joinpath(cfg.mger.outdir, "method_rank_stability.csv"), rank_stability)
    CSV.write(joinpath(cfg.mger.outdir, "parameter_scan.csv"), scans)
    # How much of each hurdle prediction actually came from a local fit. Without this a
    # globally-degenerate model is indistinguishable from a local one in the metrics.
    isempty(rows.hurdle) ||
        CSV.write(joinpath(cfg.mger.outdir, "hurdle_diagnostics.csv"), DataFrame(rows.hurdle))
    nrow(bootstrap) > 0 && CSV.write(joinpath(cfg.mger.outdir, "paired_comparisons.csv"), bootstrap)
    CSV.write(joinpath(cfg.mger.outdir, "run_status.csv"), status)
    # Which method `auto` picked in each fold, and by how much. A `chosen` column that changes
    # from fold to fold is the honest reading of "no single GWR variant is best here", and the
    # gap to `runner_up_rmse` says whether the choice was decisive or a coin flip.
    isempty(rows.auto_selection) ||
        CSV.write(joinpath(cfg.mger.outdir, "auto_selection.csv"), DataFrame(rows.auto_selection))
    isempty(rows.blend_selection) ||
        CSV.write(joinpath(cfg.mger.outdir, "blend_selection.csv"), DataFrame(rows.blend_selection))
    # The fusion coefficients each derived product's anchor was built from, per fold. Worth
    # reporting for the same reason `auto_selection.csv` is: they are fitted quantities the
    # published numbers rest on, and their stability across folds is the evidence that the anchor
    # is a property of the products rather than of one partition. `fell_back` marks a fold whose
    # normal equations were singular and which therefore used equal weights instead.
    isempty(rows.fused_anchor) ||
        CSV.write(joinpath(cfg.mger.outdir, "fused_anchor_selection.csv"),
            DataFrame(rows.fused_anchor))
    _dem_enabled(cfg) && _write_dem_outputs(cfg.mger.outdir, rows.dem_store, scans, status)
    if joint_inputs !== nothing
        scaling = isempty(rows.joint_scaling) ? DataFrame() :
            vcat(rows.joint_scaling...; cols=:union)
        quality = isempty(rows.joint_qc) ? DataFrame() :
            vcat(rows.joint_qc...; cols=:union)
        CSV.write(joinpath(cfg.mger.outdir, "joint_fold_scaling.csv"), scaling)
        CSV.write(joinpath(cfg.mger.outdir, "joint_fold_quality_control.csv"), quality)
        bandwidths = filter([:mode, :method, :selected] =>
            (mode, method, selected) -> mode == "residual" &&
                method in ("gwr", "mixed_gwr", "mgwr") && selected, scans)
        CSV.write(joinpath(cfg.mger.outdir, "joint_bandwidths.csv"), bandwidths)
        CSV.write(joinpath(cfg.mger.outdir, "covariate_model_status.csv"), filter(
            :covariate_selection_mode => in(("fixed_full_data", "nested_per_fold")), status,
        ))
        nested_joint && _write_joint_selection_outputs(cfg.mger.outdir, rows.joint_store)
    end
    if ncol(claim) > 0
        CSV.write(joinpath(cfg.mger.outdir, "claim_assessment.csv"), claim)
        # One row per (repeat, product) above; this collapses them so "does the claim hold" can be
        # separated from "did it hold in one partition". Trivial for a single-partition run.
        agreement = claim_agreement(claim)
        nrow(agreement) > 0 &&
            CSV.write(joinpath(cfg.mger.outdir, "claim_agreement.csv"), agreement)
    end
    scope = _benchmark_scope_table(
        cfg, seeds, nested_joint, joint_inputs, rows.selection_fallback_cells)
    append!(scope, provenance)
    # A tree edited while the run was going: `git_dirty` above is still the code that ran.
    end_git = _git_provenance()
    append!(scope, DataFrame(key=["git_dirty_at_end"], value=end_git.value[end_git.key .== "git_dirty"]))
    CSV.write(joinpath(cfg.mger.outdir, "benchmark_scope.csv"), scope)
    return (; metrics, scans, bootstrap, status, claim, repeat_summary, fold_summary,
        rank_stability, auto_selection=DataFrame(rows.auto_selection))
end

# Two levels up: this file lives in `src/benchmark/`. Derived rather than configured - `src/`
# must not read fixed absolute paths, and the process working directory is whatever the calling
# script was launched from.
const _REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))

"""
Which code produced this run, as `scope` rows: `RunProvenance.git_provenance` for this checkout.

`benchmark_scope.csv` recorded the CV scope and nothing about provenance, so a published run
could be dated but not attributed. `git_dirty` is the load-bearing one: a clean tree means
`git_commit` fully determines the code, and a dirty tree means it does not, which is the
difference between a citable baseline and a plausible one. `run_interpolation_benchmark` reads
it at the start, through `_start_provenance`: read at the end, a run launched on uncommitted
edits that were committed while it ran recorded itself as clean. `git_dirty_at_end` is the
second reading.
"""
_git_provenance() = DataFrame(git_provenance(_REPO_ROOT))

"""Every input file `cfg` reads: the gauges, the satellite products, and the covariate tables."""
function _benchmark_input_paths(cfg::InterpolationBenchmarkConfig)
    mger = cfg.mger
    paths = String[mger.station_meta_path, mger.obs_hourly_wide_path]
    append!(paths, [mger.sat_paths[product] for product in sort(collect(keys(mger.sat_paths)))])
    cfg.terrain_path === nothing || push!(paths, cfg.terrain_path)
    joint = cfg.joint_covariates
    if joint !== nothing
        joint.spec_path === nothing || push!(paths, joint.spec_path)
        push!(paths, joint.terrain_path)
        append!(paths, [joint.era5_paths[year_value] for year_value in sort(collect(keys(joint.era5_paths)))])
        joint.ndvi_path === nothing || push!(paths, joint.ndvi_path)
    end
    return paths
end

"""
The provenance rows a run records at its start (`RunProvenance.start_provenance`): the git state,
the Julia version, the SHA-256 of the active `Manifest.toml` (copied into the run directory, since
it is gitignored), and the SHA-256 of every file in `_benchmark_input_paths`, since `data/` is
outside git.
"""
_start_provenance(cfg::InterpolationBenchmarkConfig) = DataFrame(
    start_provenance(_REPO_ROOT, cfg.mger.outdir; inputs=_benchmark_input_paths(cfg)))

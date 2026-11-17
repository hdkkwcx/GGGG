# What a benchmark config must satisfy before a run starts. `_validate_benchmark_config` calls the
# groups below in a fixed order and each group keeps its checks in order, so the first rule a
# config breaks is the one reported: the order the checks had as one 176-line function.

"""`true` when the run loads FY4B, which several options change the use of."""
_loads_fy4b(cfg::InterpolationBenchmarkConfig) = haskey(cfg.mger.sat_paths, "FY4B")

"""
`true` when something reads the fusion's sources: a fused anchor, or a blend on a non-constant
axis. An option that changes what the fusion reads of FY4B needs one of them to act on.
"""
_has_fusion_consumer(cfg::InterpolationBenchmarkConfig) =
    !isempty(cfg.fused_anchor_variants) ||
    (cfg.satellite_wet_blend && any(!=(:constant), cfg.blend_axes))

"""The folds, the candidate grids and the tuning options."""
function _validate_cv_and_tuning(cfg::InterpolationBenchmarkConfig, n_station::Int)
    2 <= cfg.k <= n_station || throw(ArgumentError("k must be between 2 and station count"))
    length(unique(benchmark_seeds(cfg))) == length(benchmark_seeds(cfg)) ||
        throw(ArgumentError("seeds must be unique"))
    cfg.fold_center_init in (:hilbert, :kmeanspp, :farthest) ||
        throw(ArgumentError("fold_center_init must be :hilbert, :kmeanspp, or :farthest"))
    all(s -> s in (:balanced_spatial, :random, :strip), cfg.cv_schemes) ||
        throw(ArgumentError("cv_schemes must contain :balanced_spatial, :random, or :strip"))
    length(unique(cfg.cv_schemes)) == length(cfg.cv_schemes) ||
        throw(ArgumentError("cv_schemes must be unique"))
    all(>(0), cfg.idw_powers) || throw(ArgumentError("IDW/ADW powers must be positive"))
    all(x -> x === nothing || x > 0, cfg.neighbor_candidates) ||
        throw(ArgumentError("neighbor candidates must be positive or nothing"))
    all(>(0), cfg.tps_smooth_candidates) ||
        throw(ArgumentError("TPS smoothing candidates must be positive"))
    0 < cfg.min_tuning_coverage <= 1 ||
        throw(ArgumentError("min_tuning_coverage must be in (0, 1]"))
    cfg.bootstrap_reps >= 0 || throw(ArgumentError("bootstrap_reps must be non-negative"))
    cfg.tuning_max_times >= 0 || throw(ArgumentError("tuning_max_times must be non-negative"))
    cfg.tuning_time_weighting in (:stratified, :uniform) ||
        throw(ArgumentError("tuning_time_weighting must be :stratified or :uniform"))
    cfg.tuning_geometry in (:inner_spatial, :loocv) ||
        throw(ArgumentError("tuning_geometry must be :inner_spatial or :loocv"))
    cfg.tuning_inner_k == 0 || cfg.tuning_inner_k >= 2 ||
        throw(ArgumentError("tuning_inner_k must be 0 (use cfg.k) or at least 2"))
    cfg.mgwr_max_tuning_iterations > 0 ||
        throw(ArgumentError("mgwr_max_tuning_iterations must be positive"))
    isempty(cfg.residual_shrinkage_candidates) &&
        throw(ArgumentError("residual_shrinkage_candidates must not be empty"))
    all(s -> 0 < s <= 1, cfg.residual_shrinkage_candidates) ||
        throw(ArgumentError("residual shrinkage candidates must be in (0, 1]"))
    return nothing
end

"""The covariate path: the legacy DEM screen and the joint covariates exclude each other."""
function _validate_covariate_paths(cfg::InterpolationBenchmarkConfig)
    xor(cfg.terrain_path === nothing, cfg.dem === nothing) && throw(ArgumentError(
        "terrain_path and dem must either both be configured or both be omitted",
    ))
    cfg.joint_covariates !== nothing && cfg.terrain_path !== nothing && throw(ArgumentError(
        "joint covariates and legacy DEM screening cannot be enabled together",
    ))
    if cfg.dem !== nothing
        cfg.dem.min_wet_hours > 0 || throw(ArgumentError("DEM min_wet_hours must be positive"))
        cfg.dem.screen_permutations > 0 ||
            throw(ArgumentError("DEM screen_permutations must be positive"))
        cfg.dem.spatial_permutations > 0 ||
            throw(ArgumentError("DEM spatial_permutations must be positive"))
        all(>(1), cfg.dem.bandwidth_candidates) ||
            throw(ArgumentError("DEM bandwidth candidates must exceed one neighbor"))
    end
    return nothing
end

"""The blends and the stack, which need the inner selection split, and the blend axes."""
function _validate_blends(cfg::InterpolationBenchmarkConfig)
    # Refused on the legacy DEM path for the same reason `auto` is skipped there:
    # `inner_selection_prediction` deliberately carries no `dem_context`, so there is no way to
    # choose a blending weight inside the fold.
    cfg.satellite_wet_blend && cfg.dem !== nothing && throw(ArgumentError(
        "satellite_wet_blend needs an inner selection split to choose its weight on, which the " *
        "legacy DEM path does not provide",
    ))
    cfg.stacked_blend && cfg.dem !== nothing && throw(ArgumentError(
        "stacked_blend fits its weights on the inner selection split, which the legacy DEM " *
        "path does not provide",
    ))
    isempty(cfg.blend_axes) && throw(ArgumentError("blend_axes must not be empty"))
    all(axis -> axis in BLEND_AXES, cfg.blend_axes) ||
        throw(ArgumentError("blend_axes must be drawn from $BLEND_AXES"))
    length(unique(cfg.blend_axes)) == length(cfg.blend_axes) ||
        throw(ArgumentError("blend_axes must be unique"))
    return nothing
end

"""The fused anchors: known variants, and nested selection for a product with no file."""
function _validate_fused_anchors(cfg::InterpolationBenchmarkConfig)
    all(variant -> variant in FUSION_VARIANTS || variant in GROUPED_FUSION_VARIANTS,
        cfg.fused_anchor_variants) || throw(ArgumentError(
        "fused_anchor_variants must be drawn from $FUSION_VARIANTS or $GROUPED_FUSION_VARIANTS",
    ))
    length(unique(cfg.fused_anchor_variants)) == length(cfg.fused_anchor_variants) ||
        throw(ArgumentError("fused_anchor_variants must be unique"))
    # A derived product has no file, so a fixed full-data specification has no role map for it and
    # `joint_inputs.specification.role_maps[product]` would throw mid-run. Nested per-fold
    # selection screens whatever product it is handed, which is what these need.
    isempty(cfg.fused_anchor_variants) || cfg.joint_covariates === nothing ||
        cfg.joint_selection !== nothing || throw(ArgumentError(
            "fused anchor products need nested per-fold covariate selection: a product built " *
            "from other products has no entry in a fixed full-data specification, so " *
            "joint_inputs.specification.role_maps would have nothing to hand it",
        ))
    return nothing
end

"""
FY4B's options: the calibrated product and what it needs, then the options that change what the
fusion reads of FY4B, which exclude each other.
"""
function _validate_fy4b_options(cfg::InterpolationBenchmarkConfig)
    if cfg.fy4b_qm_grid !== nothing
        _loads_fy4b(cfg) || throw(ArgumentError(
            "fy4b_qm_grid calibrates FY4B, which this run does not load"))
        for path in (cfg.fy4b_qm_grid, _fy4b_hour_meta_path(cfg))
            isfile(path) || throw(ArgumentError("fy4b_qm_grid needs $path"))
        end
        # Same reason as the fused products: it has no file, so no fixed specification entry.
        cfg.joint_covariates === nothing || cfg.joint_selection !== nothing ||
            throw(ArgumentError("$FY4B_QM_PRODUCT needs nested per-fold covariate selection"))
    end
    (cfg.fy4b_qm_adaptive || cfg.fy4b_qm_fusion) && cfg.fy4b_qm_grid === nothing &&
        throw(ArgumentError("fy4b_qm_adaptive and fy4b_qm_fusion need fy4b_qm_grid"))
    cfg.fy4b_qm_crossed && cfg.fy4b_qm_grid === nothing && throw(ArgumentError(
        "fy4b_qm_crossed fits $FY4B_QM_PRODUCT out of time, which needs fy4b_qm_grid"))
    if cfg.fy4b_qm_fusion
        # The grouped variants and the stack build their designs and bands once from the sources,
        # so they cannot take a source that changes with the fold.
        any(in(GROUPED_FUSION_VARIANTS), cfg.fused_anchor_variants) && throw(ArgumentError(
            "fy4b_qm_fusion cannot feed the grouped fusion variants $GROUPED_FUSION_VARIANTS"))
        cfg.stacked_blend && throw(ArgumentError("fy4b_qm_fusion cannot feed the stacked blend"))
        _has_fusion_consumer(cfg) ||
            throw(ArgumentError("fy4b_qm_fusion needs a fused anchor or a non-constant blend " *
                "axis to feed"))
    end
    if cfg.fy4b_qm_rh_gate
        cfg.fy4b_qm_grid === nothing &&
            throw(ArgumentError("fy4b_qm_rh_gate gates $FY4B_QM_PRODUCT, which needs fy4b_qm_grid"))
        cfg.joint_covariates === nothing && throw(ArgumentError(
            "fy4b_qm_rh_gate reads ERA5-Land t2m and d2m from the joint covariates, which this " *
            "run does not load"))
    end
    if cfg.fy4b_twostage
        cfg.fy4b_qm_grid === nothing && throw(ArgumentError(
            "fy4b_twostage reads the offset grid at each fold's pixel shift, which needs fy4b_qm_grid"))
        cfg.fy4b_qm_crossed && throw(ArgumentError(
            "fy4b_twostage counts its tables on every hour, so beside a fy4b_qm_crossed " *
            "$FY4B_QM_PRODUCT the two would not be fitted alike"))
    end
    if cfg.fy4b_out_of_fusion
        _loads_fy4b(cfg) || throw(ArgumentError(
            "fy4b_out_of_fusion leaves FY4B out of the fusion, but this run does not load FY4B"))
        cfg.fy4b_qm_fusion && throw(ArgumentError(
            "fy4b_out_of_fusion and fy4b_qm_fusion contradict each other"))
        _has_fusion_consumer(cfg) ||
            throw(ArgumentError("fy4b_out_of_fusion needs a fused anchor or a non-constant " *
                "blend axis to leave FY4B out of"))
    end
    if cfg.fy4b_confirm_fusion
        _loads_fy4b(cfg) || throw(ArgumentError(
            "fy4b_confirm_fusion confirms FY4B's rain, but this run does not load FY4B"))
        length(cfg.mger.sat_paths) >= 2 || throw(ArgumentError(
            "fy4b_confirm_fusion needs another product to confirm FY4B's rain with"))
        cfg.fy4b_qm_fusion && throw(ArgumentError(
            "fy4b_confirm_fusion and fy4b_qm_fusion both replace FY4B in the fusion"))
        cfg.fy4b_out_of_fusion && throw(ArgumentError(
            "fy4b_confirm_fusion and fy4b_out_of_fusion contradict each other"))
        _has_fusion_consumer(cfg) ||
            throw(ArgumentError("fy4b_confirm_fusion needs a fused anchor or a non-constant " *
                "blend axis to feed"))
    end
    return nothing
end

"""The blend weight grid is a constant; this guards an edit to it."""
function _validate_blend_lambdas()
    # 0.0 is load-bearing twice over: it is the "do not blend" option, so a fold that cannot be
    # helped is not forced to blend anyway, and it is the reference `unblended_RMSE` that
    # `blend_selection.csv` reports the chosen weight against.
    0.0 in BLEND_LAMBDAS ||
        throw(ArgumentError("BLEND_LAMBDAS must contain 0.0, the unblended identity"))
    all(lambda -> 0.0 <= lambda <= 1.0, BLEND_LAMBDAS) ||
        throw(ArgumentError("BLEND_LAMBDAS must lie in [0, 1]"))
    return nothing
end

"""The joint covariates: one selection mode, its files, and the back-fit settings."""
function _validate_joint(cfg::InterpolationBenchmarkConfig)
    if cfg.joint_covariates !== nothing
        joint = cfg.joint_covariates
        xor(joint.spec_path === nothing, cfg.joint_selection === nothing) || throw(ArgumentError(
            "joint covariates require exactly one of spec_path (fixed full-data selection) or " *
            "joint_selection (nested per-fold selection)",
        ))
        joint.spec_path === nothing || cfg.exploratory_only || throw(ArgumentError(
            "joint covariates selected once on the full station set (spec_path) let the outer " *
            "held-out fold influence which covariates and local/global roles its own predictions " *
            "use. Use joint_selection for nested per-fold selection, or set exploratory_only=true " *
            "to acknowledge that this run's metrics are not a reportable result",
        ))
        joint.spec_path === nothing || isfile(joint.spec_path) ||
            throw(ArgumentError("joint specification does not exist"))
        isfile(joint.terrain_path) || throw(ArgumentError("joint terrain table does not exist"))
        all(>(1), joint.bandwidth_candidates) ||
            throw(ArgumentError("joint bandwidth candidates must exceed one neighbor"))
        joint.max_iterations > 0 || throw(ArgumentError("joint max_iterations must be positive"))
        joint.tolerance > 0 || throw(ArgumentError("joint tolerance must be positive"))
        # Successive over-relaxation diverges outside (0, 2); see the field comment in
        # `JointCovariateBenchmarkConfig` for why the default sits above 1.
        0 < joint.relaxation < 2 ||
            throw(ArgumentError("joint relaxation must be in (0, 2)"))
        joint.mgwr_spatial_grouping in (:split, :shared, :intercept_only) || throw(ArgumentError(
            "mgwr_spatial_grouping must be :split, :shared, or :intercept_only",
        ))
    end
    return nothing
end

function _validate_benchmark_config(cfg::InterpolationBenchmarkConfig, n_station::Int)
    _validate_cv_and_tuning(cfg, n_station)
    _validate_covariate_paths(cfg)
    _validate_blends(cfg)
    _validate_fused_anchors(cfg)
    _validate_fy4b_options(cfg)
    _validate_blend_lambdas()
    _validate_joint(cfg)
    return nothing
end

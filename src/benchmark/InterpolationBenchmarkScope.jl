"""
The rows of `benchmark_scope.csv` that say what a run did, as a `key, value` table: what was held
out and what the result supports, how candidates were tuned and a method chosen, the covariate
path, and one or more rows for each option that is on. `_write_benchmark_outputs` appends the
provenance rows and writes the file.

`selection_fallback_cells` names the cells whose inner selection split fell back to leave-one-out
(`_benchmark_rows`), which the `tuning_geometry` row reports.
"""
function _benchmark_scope_table(
    cfg::InterpolationBenchmarkConfig, seeds, nested_joint::Bool, joint_inputs,
    selection_fallback_cells,
)
    scope = DataFrame(
        key=[
            "repeated_cv_partitions", "repeated_cv_seeds", "fold_center_init",
            "fold_rotations",
            "validation_target", "training_signal", "temporal_holdout", "primary_cv",
            "secondary_cv", "common_evaluation_mask", "tuning_time_limit",
            "tuning_time_weighting", "tuning_geometry",
            "model_selection", "per_repeat_oof_tables", "results_admissible",
            "dem_variable_selection", "dem_role_assignment", "dem_leakage_control",
            "dem_empty_selection", "supported_claim", "unsupported_claim",
        ],
        value=[
            string(length(seeds)), join(seeds, ","), string(cfg.fold_center_init),
            cfg.fold_center_init === :hilbert ?
                "seed-free: partition i is rotation i-1 of the Hilbert frame (" *
                    join(0:(length(seeds) - 1), ",") * "); rotation 0 is the canonical one" :
                "not applicable: $(cfg.fold_center_init) draws its initial centers from the seed",
            "held-out stations at matched observation/satellite timestamps",
            "concurrent training-station observations for direct methods; observation-minus-satellite residuals plus fold-selected DEM variables for residual_gwr, mixed_gwr, and mgwr",
            "false", "balanced two-dimensional spatial 5-fold", "random station 5-fold",
            "true across all eight methods", string(cfg.tuning_max_times),
            cfg.tuning_time_weighting === :stratified ?
                "stratified: wettest eighth taken with certainty, remaining hours systematically subsampled and inverse-probability weighted so the tuning RMSE estimates the reported pooled RMSE" :
                "uniform (default): unweighted RMSE over a wet-oversampled subsample, so the tuning RMSE runs ~3x the reported pooled RMSE it stands in for; the level error is close to a constant multiplier and cancels in the candidate ranking",
            cfg.tuning_geometry === :inner_spatial ?
                "inner_spatial (default): candidates scored by out-of-fold prediction onto an inner $(cfg.tuning_inner_k == 0 ? cfg.k : cfg.tuning_inner_k)-group split of the training stations, built with the same splitter and scheme as the outer partition, so selection and reporting are the same estimand" *
                (isempty(selection_fallback_cells) ? "" :
                    "; FELL BACK to leave-one-out in $(length(selection_fallback_cells)) cell(s) too small to split: $(join(selection_fallback_cells, ", "))") :
                "loocv (legacy): candidates scored by leave-one-out at training stations, which measures interpolation next to a retained gauge while the reported metric measures extrapolation to a station 20+ km from any gauge",
            "auto: one GWR-family method chosen per fold on the inner selection split ($(join([_output_method(mode, method) for (mode, method) in AUTO_CANDIDATE_RUNS], ", "))), contenders re-predicted and compared on the intersection of their masks; every method is also reported on its own, but naming a per-method winner from those columns is selection on the held-out fold, which is what auto exists to avoid",
            # `oof_*.csv` and `common_evaluation_mask.csv` are large, so only the first partition
            # writes them. Stated here because the claim path no longer assumes repeat 1 exists.
            length(seeds) == 1 ? "written for the single partition" :
                "written for repeat_01 only; later partitions report metrics but not per-station OOF tables",
            cfg.exploratory_only ?
                "NO - exploratory_only=true: joint covariates and roles were selected once over every station, so the outer held-out fold helped choose the covariates its own predictions use" :
                "yes - every selection step ran inside its training fold",
            _dem_enabled(cfg) ? "training-fold Pearson/Spearman direction check, joint aspect F test, BH q<0.05, and grouped VIF<5" : "disabled",
            _dem_enabled(cfg) ? "training-fold GWR Monte Carlo spatial nonstationarity test with BH q<0.05" : "disabled",
            _dem_enabled(cfg) ? "all DEM screening, scaling, role tests, and bandwidth selection use training stations only" : "not applicable",
            _dem_enabled(cfg) ? "fall back to the original spatial-only residual design and record dem_variable_count=0" : "not applicable",
            "gauge-network-assisted interpolation/correction at stations not used for fitting",
            "temporal forecast or correction without concurrent gauge observations",
        ],
    )
    if joint_inputs !== nothing
        joint = something(cfg.joint_covariates)
        append!(scope, DataFrame(
            key=["mgwr_spatial_grouping", "backfit_relaxation", "backfit_max_iterations",
                "residual_shrinkage"],
            value=[
                joint.mgwr_spatial_grouping === :split ?
                    "split (default): the intercept, longitude and latitude each form their own single-column back-fitting group" :
                    joint.mgwr_spatial_grouping === :shared ?
                        "shared: intercept, longitude and latitude solved as one weighted least squares, so only the covariates carry separate bandwidths" :
                        "intercept_only: the coordinate columns are dropped; a locally varying intercept plus one group per local covariate",
                "$(joint.relaxation) (successive over-relaxation on the back-fitting sweeps; converges to the same fixed point for any admissible value, so this trades rate only — plain Gauss-Seidel at 1.0 fails to converge on most hours)",
                string(joint.max_iterations),
                length(cfg.residual_shrinkage_candidates) == 1 &&
                    only(cfg.residual_shrinkage_candidates) == 1.0 ?
                    "disabled: the residual correction is added back unshrunk, as GWR produces it" :
                    "enabled: a scale in $(minimum(cfg.residual_shrinkage_candidates))-$(maximum(cfg.residual_shrinkage_candidates)) is selected with the bandwidth on the inner spatial split and applied to the joint dynamic models' residual correction",
            ],
        ))
        training_row = findfirst(==("training_signal"), scope.key)
        scope.value[training_row] = nested_joint ?
            "concurrent training-station observations for direct methods; observation-minus-satellite residuals plus a per-fold, training-station-only product-specific DEM/ERA5 covariate specification for residual_gwr, mixed_gwr, and mgwr" :
            "concurrent training-station observations for direct methods; observation-minus-satellite residuals plus a fixed full-data product-specific DEM/ERA5 covariate specification for residual_gwr, mixed_gwr, and mgwr"
        append!(scope, DataFrame(
            key=[
                "covariate_selection_mode", "covariate_selection_leakage",
                "covariate_parameter_leakage", "covariate_temporal_fit",
                "inference_policy",
            ],
            value=nested_joint ? [
                "nested_per_fold",
                "variable and role selection re-run inside every training fold; validation stations never seen by selection",
                "scaling, bandwidth selection, coefficients, and predictions use training stations only",
                "hour-specific cross-sectional spatial fitting with matched ERA5 state",
                "eligible for the same paired significance test and claim assessment as the legacy DEM path",
            ] : [
                "fixed_full_data",
                "non-nested variable selection used all 237 stations before spatial cross-validation",
                "scaling, bandwidth selection, coefficients, and predictions use training stations only",
                "hour-specific cross-sectional spatial fitting with matched ERA5 state",
                "exploratory performance metrics only; no paired significance or claim assessment",
            ],
        ))
    end
    if cfg.satellite_wet_blend
        append!(scope, DataFrame(
            key=["satellite_wet_blend", "satellite_wet_blend_selection", "blend_axes"],
            value=[
                "blended counterparts of $(join(BLEND_SOURCE_METHODS, ", ")) reported as " *
                "$(join([_blend_method(m, a) for a in cfg.blend_axes
                         for m in BLEND_SOURCE_METHODS], ", ")): where the " *
                "satellite reports at least $(cfg.mger.rain_threshold) mm, the prediction is " *
                "blended toward $(BLEND_FALLBACK_METHOD)",
                "blending weight chosen per fold on the inner selection split over " *
                "$(BLEND_LAMBDAS), never on held-out cells; see blend_selection.csv. The " *
                "blended methods are scored on the shared mask but do not define it, so adding " *
                "them leaves every other method's denominator unchanged",
                "$(join(cfg.blend_axes, ", ")). `constant` is one weight for every " *
                "satellite-wet cell. `agreement_envelope` gives one weight per band of (how " *
                "many source products call the cell wet) x (their maximum, on " *
                "$(BLEND_ENVELOPE_EDGES)), which widens the intervention set: a cell this " *
                "product calls dry but another calls wet is blended, where the constant axis " *
                "would not have touched it" *
                (any(in(BLEND_DRY_AXES), cfg.blend_axes) ? ". `agreement_envelope_dry` adds a " *
                    "weight for the dry band (no source product wet), blended toward " *
                    "$(BLEND_FALLBACK_METHOD) and chosen on the evenly spaced tuning hours " *
                    "only; `agreement_envelope_dry_raw` blends that band toward the product's " *
                    "own anchor instead" : ""),
            ],
        ))
    end
    if cfg.stacked_blend
        append!(scope, DataFrame(
            key=["stacked_blend", "stacked_blend_selection"],
            value=[
                "$(STACK_METHOD): non-negative combination of $(join(STACK_FEATURES, ", ")), " *
                "one weight vector per band of (agreement-envelope band including no-product-wet) " *
                "x (nearest training gauge at $(STACK_DISTANCE_EDGES) km)",
                "weights fitted per fold on the inner selection split over the full training " *
                "record, never on held-out cells; see stack_selection.csv. Scored on the shared " *
                "mask without defining it",
            ],
        ))
    end
    if !isempty(cfg.fused_anchor_variants)
        append!(scope, DataFrame(
            key=["fused_anchor", "fused_anchor_leakage"],
            value=[
                "derived products $(join([_fused_product_name(v) for v in
                    sort(cfg.fused_anchor_variants)], ", ")): the anchor is a combination of the " *
                "products that have files, not a satellite field as shipped. Their `raw` row is " *
                "therefore a gauge-trained predictor that uses only satellites at prediction " *
                "time, and is reported so the fusion's contribution can be read separately from " *
                "the correction built on top of it",
                "fusion coefficients fitted inside every fold on that fold's training stations " *
                "alone and applied to every station, so a held-out gauge never reaches the " *
                "anchor its own prediction is built on; see fused_anchor_selection.csv. A " *
                "derived product's evaluation mask needs every source product finite, so it is " *
                "a different cell population from any of them",
            ],
        ))
    end
    if cfg.fy4b_qm_grid !== nothing
        append!(scope, DataFrame(
            key=["fy4b_qm", "fy4b_qm_leakage"],
            value=[
                "derived product $(FY4B_QM_PRODUCT): FY4B from $(cfg.fy4b_qm_grid), its sample " *
                "pixel shifted per regime (satellite longitude x software version) and then " *
                "quantile-mapped per season and regime, with a wet threshold matching the " *
                "gauges' wet frequency. Its `raw` row is a gauge-trained predictor that uses only " *
                "FY4B at prediction time",
                "shift and map chosen inside every fold on that fold's training stations alone " *
                "and applied to every station, so a held-out gauge never reaches the anchor its " *
                "own prediction is built on; see fy4b_qm_shift_selection.csv and fy4b_qm_maps.csv",
            ],
        ))
    end
    if cfg.fy4b_qm_adaptive
        append!(scope, DataFrame(key=["fy4b_qm_adaptive"], value=[
            "each regime's season map applied, month by month (regime x year-month), only as " *
            "strongly as the fold's training stations support (strength " *
            "$(join(QM_STRENGTHS, "/")) after the wet threshold, or no adjustment), chosen by " *
            "training squared error; see fy4b_qm_strengths.csv",
        ]))
    end
    if cfg.fy4b_qm_crossed
        append!(scope, DataFrame(key=["fy4b_qm_crossed"], value=[
            "$(FY4B_QM_PRODUCT) fitted out of time as well: the hours of each 7-day time fold (5 " *
            "folds, balanced within season x regime, seed $(FY4B_QM_CROSSED_TIME_FOLD_SEED)) take " *
            "the shift, map and strength fitted on the fold's training stations at the other time " *
            "folds' hours, a day either side left out, so no week's anchor saw that week at any " *
            "gauge; see fy4b_qm_time_folds.csv and the time_fold column of " *
            "fy4b_qm_shift_selection.csv, fy4b_qm_maps.csv and fy4b_qm_strengths.csv",
        ]))
    end
    if cfg.fy4b_qm_fusion
        append!(scope, DataFrame(key=["fy4b_qm_fusion"], value=[
            "the fused anchors and the blend's non-constant bands take each fold's " *
            "$(FY4B_QM_PRODUCT) anchor in place of the plain FY4B; the products with files keep " *
            "the raw and unblended rows of a run without this option",
        ]))
    end
    if cfg.fy4b_qm_rh_gate
        append!(scope, DataFrame(key=["fy4b_qm_rh_gate"], value=[
            "derived product $(FY4B_QM_RH_PRODUCT): each fold's $(FY4B_QM_PRODUCT) set to 0 where " *
            "ERA5-Land relative humidity at the gauge is below a cut in $(first(GATE_CUTS))-" *
            "$(last(GATE_CUTS))% and the value is below $(GATE_SPARED_MM) mm/h. The cut raises " *
            "training CSI the most, chosen per fold on training stations and per 7-day time fold on " *
            "the other blocks' hours (seed $(FY4B_QM_RH_TIME_FOLD_SEED)); see " *
            "fy4b_qm_rh_gate_cuts.csv. It stays out of the fusion. ERA5-Land is reanalysis, so this " *
            "is FY4B + reanalysis",
        ]))
    end
    if cfg.fy4b_twostage
        append!(scope, DataFrame(key=["fy4b_twostage"], value=[
            "derived product $(FY4B_TS_PRODUCT): FY4B replaced by P(gauge wet) x the mean wet gauge " *
            "amount of its context class (FY4B's intensity, rain in the surrounding shifts, rain " *
            "at t-1 and t+1, within season x regime), counted per fold on that fold's training " *
            "stations at its pixel shift and applied to every station; a class with fewer than " *
            "$(TS_MIN_PAIRS) training pairs is pooled. A conditional mean, so it flattens heavy " *
            "rain. The tables are not held out in time. It stays out of the fusion; see " *
            "fy4b_twostage_tables.csv",
        ]))
    end
    if cfg.fy4b_out_of_fusion
        append!(scope, DataFrame(key=["fy4b_out_of_fusion"], value=[
            "FY4B is scored but left out of the fused anchors and the blend's non-constant bands; " *
            "the product list, hours and selection seeds stay those of a run with FY4B in the " *
            "fusion, so the products with files keep that run's raw and unblended rows",
        ]))
    end
    if cfg.fy4b_confirm_fusion
        append!(scope, DataFrame(key=["fy4b_confirm_fusion"], value=[
            "the fused anchors and the blend's non-constant bands take FY4B only where another " *
            "product calls the cell wet (>= $(cfg.mger.rain_threshold) mm/h), and 0 where FY4B is " *
            "wet alone; gauge-free and fitted on nothing, so the same in every fold. The products " *
            "with files keep the raw and unblended rows of a run without this option",
        ]))
    end
    grouped = [v for v in sort(cfg.fused_anchor_variants) if v in GROUPED_FUSION_VARIANTS]
    if !isempty(grouped)
        append!(scope, DataFrame(
            key=["fused_anchor_grouped"],
            value=[
                "$(join([_fused_product_name(v) for v in grouped], ", ")): regressed on every " *
                "source product at t, t-1 and t+1 and on each product's mean over the 8 nearest " *
                "other stations, one coefficient set per agreement-envelope band; features read " *
                "satellite values only. Coefficients per fold and band in " *
                "fused_anchor_grouped_coefficients.csv, not fused_anchor_selection.csv",
            ],
        ))
    end
    return scope
end

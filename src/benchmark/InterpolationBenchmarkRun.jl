function run_interpolation_benchmark(cfg::InterpolationBenchmarkConfig)
    mkpath(cfg.mger.outdir)
    # Read before any work: the code that runs is the code checked out now, not at the end.
    provenance = _start_provenance(cfg)
    station_meta = load_station_meta(cfg.mger.station_meta_path;
        station_id_col=cfg.mger.station_id_col, lon_col=cfg.mger.lon_col, lat_col=cfg.mger.lat_col)
    source_products, ids, product_data = load_global_common_product_data(cfg.mger)
    lonlat = build_X_lonlat(station_meta, ids)
    _validate_benchmark_config(cfg, length(ids))
    # Derived products are appended after validation, and the joint inputs below are loaded for the
    # source products only: `load_joint_benchmark_inputs` reads a fixed specification keyed by
    # product, which a product with no file cannot appear in, and which the config validation has
    # already ruled out alongside a fused anchor.
    # `fy4b_out_of_fusion` narrows only what the fusion reads. `benchmark_products` and the joint
    # inputs below keep the full list: FY4B is still scored, and a product's position in that list
    # sets its selection seed (`_run_benchmark_fold!`), so the list has to stay as it is for this
    # run to differ from one with FY4B in the fusion in the fusion alone.
    fusion_products = cfg.fy4b_out_of_fusion ? filter(!=("FY4B"), source_products) : source_products
    _, fused_specifications = _build_fused_products!(cfg, fusion_products, product_data, lonlat)
    fy4b_calibration = _prepare_fy4b_calibration!(cfg, ids, product_data)
    # Each hour's 7-day block and time fold, for a crossed FY4B_QM.
    fy4b_calibration === nothing || fy4b_calibration.time_folds === nothing ||
        CSV.write(joinpath(cfg.mger.outdir, "fy4b_qm_time_folds.csv"), DataFrame(
            time=Dates.format.(product_data[FY4B_QM_PRODUCT].times, dateformat"yyyy-mm-ddTHH:MM:SS"),
            block=fy4b_calibration.time_folds.block_of_time,
            time_fold=fy4b_calibration.time_folds.fold_of_time))
    # The gated product shares FY4B_QM's grid and gauges; its anchor is built per fold below.
    cfg.fy4b_qm_rh_gate && (product_data[FY4B_QM_RH_PRODUCT] = product_data[FY4B_QM_PRODUCT])
    # So does the two-stage product. `twostage_nearby` keeps each centre shift's neighbour counts,
    # which depend on neither the fold nor the scheme.
    cfg.fy4b_twostage && (product_data[FY4B_TS_PRODUCT] = product_data[FY4B_QM_PRODUCT])
    twostage_nearby = Dict{Int,Tuple{Matrix{UInt8},Int}}()
    products = benchmark_products(cfg, source_products)
    terrain = _dem_enabled(cfg) ? load_aligned_terrain(something(cfg.terrain_path), ids) : nothing
    common_times = product_data[first(products)].times
    joint_inputs = _joint_enabled(cfg) ? load_joint_benchmark_inputs(
        something(cfg.joint_covariates), source_products, ids, common_times,
    ) : nothing
    nested_joint = _joint_enabled(cfg) && cfg.joint_selection !== nothing
    # The humidity gate's inputs: RH at every gauge-hour from the joint covariates' ERA5-Land (read
    # at the benchmark's own alignment), and 7-day time folds balanced within FY4B's season x regime
    # cells, over which each fold's cut is cross-fitted.
    gate_inputs = cfg.fy4b_qm_rh_gate ? (;
        humidity=station_humidity(joint_inputs.era5.values[:t2m_c], joint_inputs.era5.values[:d2m_c]),
        time_folds=time_block_folds(product_data[FY4B_QM_PRODUCT].times, fy4b_calibration.cells;
            seed=FY4B_QM_RH_TIME_FOLD_SEED),
    ) : nothing
    # Band matrices for the non-constant blend axes. They are read off the source products alone,
    # so they depend on neither the fold nor the scheme nor which product is being corrected, and
    # are built once here rather than per fold. The constant axis keeps its own code path in
    # `_run_fold_blend!` so that turning another axis on cannot perturb what it reports.
    blend_bands = Dict{Symbol,Matrix{Int}}()
    band_sources = Matrix{Float64}[]
    if cfg.satellite_wet_blend
        band_sources = _fusion_source_matrices(cfg, fusion_products, product_data)
        for axis in cfg.blend_axes
            axis === :constant && continue
            blend_bands[axis] = blend_band_matrix(
                axis, first(band_sources), band_sources, cfg.mger.rain_threshold,
            )
        end
    end
    # The stack's agreement bands come from the same source products, so they too are built once.
    stack_agreement = if cfg.stacked_blend
        stack_sources = _fusion_source_matrices(cfg, fusion_products, product_data)
        blend_band_matrix(
            :agreement_envelope, first(stack_sources), stack_sources, cfg.mger.rain_threshold,
        )
    else
        nothing
    end
    joint_inputs === nothing || _write_joint_provenance(cfg, joint_inputs, nested_joint)

    rows = _benchmark_rows()

    _dem_enabled(cfg) && _screen_dem_full_data!(rows.dem_store, cfg, terrain, lonlat, products, product_data)

    seeds = benchmark_seeds(cfg)
    # Repeated cross-validation: one independent fold partition per repeat.
    #
    # `repeat_seed` and the partition index are deliberately separate. Under the default
    # `:hilbert` initialisation the partition comes from `rotation = repeat_index - 1` and no RNG
    # is involved at all; `repeat_seed` still seeds the genuinely stochastic sub-processes below
    # (paired bootstrap, DEM permutation tests, joint variable selection) and the `:random`
    # scheme, which is random by definition.
    for (repeat_index, repeat_seed) in enumerate(seeds)
        repeat_root = _repeat_dir(cfg.mger.outdir, repeat_index, length(seeds))
        for scheme_symbol in cfg.cv_schemes
            scheme = string(scheme_symbol)
            scheme_dir = joinpath(repeat_root, scheme)
            mkpath(scheme_dir)
            folds = benchmark_folds(
                scheme_symbol, ids, lonlat; k=cfg.k, seed=repeat_seed,
                center_init=cfg.fold_center_init, rotation=repeat_index - 1,
            )
            _write_split(joinpath(scheme_dir, "split_common.csv"), ids, folds, scheme_symbol)
            id_map = Dict(id => index for (index, id) in enumerate(ids))
            # A fold's calibrated FY4B depends only on its training stations, so it is built once
            # per fold here: for FY4B_QM itself and, under `fy4b_qm_fusion`, for the fused anchors
            # and blend bands it feeds. Under `fy4b_qm_crossed` each time fold's hours come from
            # a fit on the other time folds' hours alone.
            fold_calibration = fy4b_calibration === nothing ? nothing : map(1:cfg.k) do fold
                _, _, train_idx, _ = _fold_station_indices(folds, id_map, cfg.k, fold)
                if cfg.fy4b_qm_crossed
                    calibrate_crossed(
                        fy4b_calibration.grid, fy4b_calibration.crossed_sums,
                        product_data[FY4B_QM_PRODUCT].Y_obs, train_idx, fy4b_calibration.regimes,
                        fy4b_calibration.regime_of_time, fy4b_calibration.cells,
                        fy4b_calibration.time_folds.fold_of_time, product_data[FY4B_QM_PRODUCT].times;
                        cache=fy4b_calibration.cache, wet=cfg.mger.rain_threshold,
                        adaptive=cfg.fy4b_qm_adaptive, choice_of_time=fy4b_calibration.month_of_time,
                        choice_labels=fy4b_calibration.month_labels,
                    )
                else
                    calibrate_fold(
                        fy4b_calibration.grid, fy4b_calibration.sums,
                        product_data[FY4B_QM_PRODUCT].Y_obs, train_idx, fy4b_calibration.regimes,
                        fy4b_calibration.regime_of_time, fy4b_calibration.cells;
                        cache=fy4b_calibration.cache, wet=cfg.mger.rain_threshold,
                        adaptive=cfg.fy4b_qm_adaptive, choice_of_time=fy4b_calibration.month_of_time,
                        choice_labels=fy4b_calibration.month_labels,
                    )
                end
            end
            # Each fold's gated FY4B_QM: the cut chosen on that fold's training stations, per time
            # fold on the other blocks' hours, and applied at every station.
            fold_gate = gate_inputs === nothing ? nothing : map(1:cfg.k) do fold
                _, _, train_idx, _ = _fold_station_indices(folds, id_map, cfg.k, fold)
                time_crossed_gate(
                    first(fold_calibration[fold]), product_data[FY4B_QM_PRODUCT].Y_obs,
                    gate_inputs.humidity, train_idx, gate_inputs.time_folds.fold_of_time,
                    product_data[FY4B_QM_PRODUCT].times; wet=cfg.mger.rain_threshold,
                )
            end
            # Each fold's two-stage FY4B: the tables counted on that fold's training stations at its
            # pixel shift (the cache it shares with the calibration above), applied at every station.
            fold_twostage = !cfg.fy4b_twostage ? nothing : map(1:cfg.k) do fold
                _, _, train_idx, _ = _fold_station_indices(folds, id_map, cfg.k, fold)
                expected, _, table_rows = two_stage_fold(
                    fy4b_calibration.grid, fy4b_calibration.sums,
                    product_data[FY4B_TS_PRODUCT].Y_obs, train_idx, fy4b_calibration.regimes,
                    fy4b_calibration.regime_of_time, fy4b_calibration.cells,
                    product_data[FY4B_TS_PRODUCT].times, first(fold_calibration[fold]);
                    cache=fy4b_calibration.cache, nearby_cache=twostage_nearby,
                    wet=cfg.mger.rain_threshold,
                )
                (expected, table_rows)
            end
            fold_bands = if cfg.fy4b_qm_fusion && !isempty(blend_bands)
                map(1:cfg.k) do fold
                    sources = _fold_sources(cfg, source_products, band_sources, fold_calibration, fold)
                    Dict(axis => blend_band_matrix(
                        axis, first(sources), sources, cfg.mger.rain_threshold,
                    ) for axis in keys(blend_bands))
                end
            else
                nothing
            end

            for product in products
                data = product_data[product]
                y_obs = data.Y_obs
                y_sat = data.Y_sat
                fusion = get(fused_specifications, product, nothing)
                calibrated = product == FY4B_QM_PRODUCT
                gated = product == FY4B_QM_RH_PRODUCT
                two_stage = product == FY4B_TS_PRODUCT
                predictions = Dict(
                    method => fill(NaN, size(y_obs)) for method in benchmark_methods(cfg)
                )
                # A derived product has no anchor until a fold says which stations may be used to
                # fit one, so `raw` is filled per fold below instead of in one assignment here.
                fusion === nothing && !calibrated && !gated && !two_stage &&
                    (predictions["raw"] .= y_sat)
                nearest_train_distance = fill(NaN, length(ids))

                for fold in 1:cfg.k
                    fold_sat = y_sat
                    if calibrated
                        _, _, calibration_train_idx, calibration_val_idx =
                            _fold_station_indices(folds, id_map, cfg.k, fold)
                        # Shift and map fitted on the training stations, applied everywhere, for
                        # the same reason as the fused anchors below.
                        fold_sat, shift_rows, map_rows, strength_rows = fold_calibration[fold]
                        predictions["raw"][calibration_val_idx, :] =
                            fold_sat[calibration_val_idx, :]
                        context = (; scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                            n_train_station=length(calibration_train_idx))
                        append!(rows.fy4b_shift, [merge(context, row) for row in shift_rows])
                        append!(rows.fy4b_map, [merge(context, row) for row in map_rows])
                        append!(rows.fy4b_strength,
                            [merge(context, row) for row in strength_rows])
                    elseif gated
                        _, _, gate_train_idx, gate_val_idx =
                            _fold_station_indices(folds, id_map, cfg.k, fold)
                        fold_sat, cut_rows = fold_gate[fold]
                        predictions["raw"][gate_val_idx, :] = fold_sat[gate_val_idx, :]
                        context = (; scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                            n_train_station=length(gate_train_idx))
                        append!(rows.fy4b_gate, [merge(context, row) for row in cut_rows])
                    elseif two_stage
                        _, _, twostage_train_idx, twostage_val_idx =
                            _fold_station_indices(folds, id_map, cfg.k, fold)
                        fold_sat, table_rows = fold_twostage[fold]
                        predictions["raw"][twostage_val_idx, :] = fold_sat[twostage_val_idx, :]
                        context = (; scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                            n_train_station=length(twostage_train_idx))
                        append!(rows.fy4b_twostage, [merge(context, row) for row in table_rows])
                    elseif fusion !== nothing && fusion.design !== nothing
                        _, _, fusion_train_idx, fusion_val_idx =
                            _fold_station_indices(folds, id_map, cfg.k, fold)
                        design = fusion.design
                        coefficients, used, fell_back = fit_grouped_fusion(
                            y_obs, fusion.sources, design.features, design.groups,
                            design.n_group, fusion_train_idx,
                        )
                        coefficients === nothing && error(
                            "$product fold $fold: the pooled fusion fit is singular",
                        )
                        fold_sat = apply_grouped_fusion(
                            fusion.sources, design.features, coefficients, design.groups,
                        )
                        predictions["raw"][fusion_val_idx, :] = fold_sat[fusion_val_idx, :]
                        lower = lowercase.(fusion.source_names)
                        terms = vcat("intercept", lower, lower .* "_lag1",
                            lower .* "_lead1", lower .* "_nbr8")
                        for group in 1:design.n_group, (row, term) in enumerate(terms)
                            push!(rows.fused_grouped, (;
                                scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                                variant=String(fusion.variant), group=group - 1,
                                n_cell=used[group], fell_back=fell_back[group], term,
                                coefficient=coefficients[row, group],
                            ))
                        end
                    elseif fusion !== nothing
                        _, _, fusion_train_idx, fusion_val_idx =
                            _fold_station_indices(folds, id_map, cfg.k, fold)
                        sources = _fold_sources(
                            cfg, fusion.source_names, fusion.sources, fold_calibration, fold,
                        )
                        coefficients, used, fell_back = fusion_coefficients(
                            fusion.variant, y_obs, sources, fusion_train_idx,
                        )
                        # Fitted on the training stations, applied everywhere: the fold's own
                        # training rows need the same anchor its held-out rows get, or the
                        # residual target would mean two different things inside one fit.
                        fold_sat = apply_satellite_fusion(sources, coefficients)
                        predictions["raw"][fusion_val_idx, :] = fold_sat[fusion_val_idx, :]
                        push!(rows.fused_anchor, merge(
                            (; scheme, product, fold, repeat=repeat_index, seed=repeat_seed,
                                variant=String(fusion.variant),
                                n_train_station=length(fusion_train_idx), n_train_cell=used,
                                fell_back, intercept=coefficients[1]),
                            NamedTuple{Tuple(Symbol("beta_", lowercase(name))
                                             for name in fusion.source_names)}(
                                Tuple(coefficients[2:end])),
                        ))
                    end
                    _run_benchmark_fold!(
                        cfg, fold, folds, id_map, ids, products, product,
                        data, lonlat, y_obs, fold_sat, terrain, joint_inputs, nested_joint,
                        scheme, scheme_symbol, repeat_index, repeat_seed,
                        predictions, nearest_train_distance,
                        fold_bands === nothing ? blend_bands : fold_bands[fold],
                        stack_agreement, rows,
                    )
                end

                common_mask = _common_method_mask(y_obs, predictions)
                if !any(common_mask)
                    failures = ["fold=$(row.fold) method=$(row.method): $(row.error)" for
                        row in rows.status if row.scheme == scheme && row.product == product &&
                        row.status != "success"]
                    error("[$scheme/$product] no common valid OOF samples across all methods; " *
                        join(failures, " | "))
                end
                product_dir = joinpath(scheme_dir, lowercase(product))
                mkpath(product_dir)
                for method in benchmark_methods(cfg)
                    # The per-station OOF tables are large; only the first repeat writes them.
                    repeat_index == 1 && write_wide(
                        joinpath(product_dir, "oof_$(method).csv"), data.times, ids, predictions[method],
                    )
                    append_stratified_metrics!(
                        rows.metrics, scheme, product, method, data.times, y_obs,
                        predictions[method], common_mask, nearest_train_distance, cfg.event_thresholds;
                        repeat=repeat_index, seed=repeat_seed,
                    )
                end
                if repeat_index == 1
                    mask_df = DataFrame(time=Dates.format.(data.times, dateformat"yyyy-mm-ddTHH:MM:SS"))
                    for (station_index, station_id) in enumerate(ids)
                        mask_df[!, Symbol(station_id)] = common_mask[station_index, :]
                    end
                    CSV.write(joinpath(product_dir, "common_evaluation_mask.csv"), mask_df)
                end
                if scheme_symbol == :balanced_spatial && cfg.bootstrap_reps > 0
                    append!(rows.bootstrap, paired_bootstrap_rows(
                        cfg, scheme, product, data.times, y_obs, predictions, common_mask;
                        repeat=repeat_index, seed=repeat_seed,
                    ))
                end
            end
        end
    end # repeat loop

    return _write_benchmark_outputs(
        cfg, products, seeds, nested_joint, joint_inputs, rows; provenance,
    )
end

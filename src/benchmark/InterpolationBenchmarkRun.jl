"""
What is fixed for the whole run, as one NamedTuple: the stations and products the benchmark reads,
and everything built from them once.

- `ids`, `lonlat`: the common stations. `source_products` are the products with files; `products`
  adds the derived ones (`benchmark_products`); `product_data` holds each one's `times`, `Y_obs`
  and `Y_sat`.
- `terrain`, `joint_inputs`, `nested_joint`: the covariates of the legacy DEM path or of the joint
  path.
- `fused_specifications`, `fy4b_calibration`, `gate_inputs`, `twostage_nearby`: what the derived
  products need that does not depend on the fold (`InterpolationBenchmarkDerived.jl`).
- `blend_bands`, `band_sources`, `stack_agreement`: the bands of the non-constant blend axes and of
  the stack, read off the source products.

Validates `cfg` once the stations are known. Writes `global_common_time_qc.csv`, for a crossed
FY4B_QM `fy4b_qm_time_folds.csv`, and the joint provenance tables.
"""
function _benchmark_setup(cfg::InterpolationBenchmarkConfig)
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
    # The gated product shares FY4B_QM's grid and gauges; its anchor is built per fold
    # (`_benchmark_partition`).
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
    return (; cfg, ids, lonlat, products, source_products, product_data, terrain, joint_inputs,
        nested_joint, fused_specifications, fy4b_calibration, gate_inputs, twostage_nearby,
        blend_bands, band_sources, stack_agreement)
end

"""
One fold partition of the stations, a repeat of a scheme, as a NamedTuple: its `folds`, the
folder its tables go to, and each fold's derived FY4B series (`fold_calibration`, `fold_gate`,
`fold_twostage`, `fold_bands`; `nothing` when the run does not use one). Those series depend on
the fold's training stations alone, so they are built here once for every product that reads them.

Writes `split_common.csv`.
"""
function _benchmark_partition(
    setup, repeat_index::Int, repeat_seed::Int, scheme_symbol::Symbol, repeat_root,
)
    (; cfg, ids, lonlat) = setup
    scheme = string(scheme_symbol)
    scheme_dir = joinpath(repeat_root, scheme)
    mkpath(scheme_dir)
    folds = benchmark_folds(
        scheme_symbol, ids, lonlat; k=cfg.k, seed=repeat_seed,
        center_init=cfg.fold_center_init, rotation=repeat_index - 1,
    )
    _write_split(joinpath(scheme_dir, "split_common.csv"), ids, folds, scheme_symbol)
    id_map = Dict(id => index for (index, id) in enumerate(ids))
    fold_calibration = _fold_calibrations(setup, folds, id_map)
    fold_gate = _fold_gates(setup, folds, id_map, fold_calibration)
    fold_twostage = _fold_two_stage(setup, folds, id_map, fold_calibration)
    fold_bands = _fold_blend_bands(setup, fold_calibration)
    return (; scheme, scheme_symbol, scheme_dir, repeat_index, repeat_seed, folds, id_map,
        fold_calibration, fold_gate, fold_twostage, fold_bands)
end

"""
One product on one partition: every fold, then the product's out-of-fold tables, its scores on the
shared mask and its paired bootstrap.

A derived product (a fusion, or one of FY4B's corrected series) has no anchor until a fold says
which stations may be used to fit one, so its `raw` is filled fold by fold from `_fold_anchor!`
rather than in one assignment.
"""
function _run_benchmark_product!(setup, partition, product::String, rows)
    (; cfg, ids, product_data, fused_specifications) = setup
    (; scheme, scheme_symbol, scheme_dir, folds, id_map, repeat_index, repeat_seed) = partition
    data = product_data[product]
    y_obs = data.Y_obs
    y_sat = data.Y_sat
    derived = haskey(fused_specifications, product) ||
        product in (FY4B_QM_PRODUCT, FY4B_QM_RH_PRODUCT, FY4B_TS_PRODUCT)
    predictions = Dict(
        method => fill(NaN, size(y_obs)) for method in benchmark_methods(cfg)
    )
    derived || (predictions["raw"] .= y_sat)
    nearest_train_distance = fill(NaN, length(ids))

    for fold in 1:cfg.k
        train_ids, val_ids, train_idx, val_idx =
            _fold_station_indices(folds, id_map, cfg.k, fold)
        stations = (; train_ids, val_ids, train_idx, val_idx)
        # The fold's identity: the first five columns of every row it writes.
        cell = (; scheme, product, fold, repeat=repeat_index, seed=repeat_seed)
        fold_sat = y_sat
        if derived
            fold_sat = _fold_anchor!(rows, setup, partition, cell, y_obs, train_idx)
            predictions["raw"][val_idx, :] = fold_sat[val_idx, :]
        end
        _run_benchmark_fold!(
            setup, partition, cell, stations, fold_sat, predictions,
            nearest_train_distance, rows,
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
    return nothing
end

"""
Run the benchmark `cfg` describes and write its tables to `cfg.mger.outdir`: every product, on
every fold of every scheme and repeat, then the pooled scores (`_write_benchmark_outputs`, whose
result this returns).
"""
function run_interpolation_benchmark(cfg::InterpolationBenchmarkConfig)
    mkpath(cfg.mger.outdir)
    # Read before any work: the code that runs is the code checked out now, not at the end.
    provenance = _start_provenance(cfg)
    setup = _benchmark_setup(cfg)
    rows = _benchmark_rows()

    _dem_enabled(cfg) && _screen_dem_full_data!(
        rows.dem_store, cfg, setup.terrain, setup.lonlat, setup.products, setup.product_data)

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
            partition = _benchmark_partition(
                setup, repeat_index, repeat_seed, scheme_symbol, repeat_root)
            for product in setup.products
                _run_benchmark_product!(setup, partition, product, rows)
            end
        end
    end # repeat loop

    return _write_benchmark_outputs(setup, seeds, rows; provenance)
end

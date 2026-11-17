"""
Products derived by fusing the ones that have files, keyed by product name.

Each entry carries the variant and the source matrices, which are the real products' own `Y_sat`
in `products` order - so the fusion's design is the same one everywhere, and
`fused_anchor_selection.csv` can name its coefficients.

The derived product shares the source products' `times`, `ids` and `Y_obs`:
`load_global_common_product_data` has already put every product on one common station set and one
common time grid, so there is nothing to align. Its `Y_sat` starts as NaN and is filled per fold,
because the anchor depends on which stations were held out - which is exactly what keeps a held-out
gauge out of the anchor its own prediction is built on.
"""
function _build_fused_products!(
    cfg::InterpolationBenchmarkConfig, products::Vector{String}, product_data::Dict{String,Any},
    lonlat::AbstractMatrix,
)
    isempty(cfg.fused_anchor_variants) && return String[], Dict{String,Any}()
    source_names = copy(products)
    sources = _fusion_source_matrices(cfg, source_names, product_data)
    template = product_data[first(source_names)]
    derived = String[]
    specifications = Dict{String,Any}()
    for variant in sort(cfg.fused_anchor_variants)
        name = _fused_product_name(variant)
        haskey(product_data, name) &&
            throw(ArgumentError("derived product $name collides with a configured product"))
        product_data[name] = (;
            template.times, template.ids, template.Y_obs,
            Y_sat=fill(NaN, size(template.Y_obs)),
        )
        # A grouped variant's design and bands depend only on the source products, so they are
        # built once here; only the coefficients are refitted per fold.
        design = if variant in GROUPED_FUSION_VARIANTS
            distance = haversine_distance_matrix(lonlat, lonlat)
            neighbours = [filter(!=(station), sortperm(distance[station, :]))
                          for station in axes(distance, 1)]
            groups = blend_band_matrix(
                :agreement_envelope, first(sources), sources, cfg.mger.rain_threshold,
            ) .+ 1
            (; features=fusion_features(sources, template.times, neighbours), groups,
                n_group=blend_band_count(:agreement_envelope) + 1)
        else
            nothing
        end
        specifications[name] = (; variant, source_names, sources, design)
        push!(derived, name)
    end
    return derived, specifications
end

"""
What the fusion reads of the products `names`, as fresh copies in that order: each product's own
`Y_sat`, except that under `fy4b_confirm_fusion` FY4B is `confirmed_by_others` the rest of `names`
at the bands' wet threshold. That reads no gauge, so one matrix serves every fold and scheme, and
the fused anchors, the blend bands and the stack's bands all see the same FY4B.
"""
function _fusion_source_matrices(cfg::InterpolationBenchmarkConfig, names, product_data)
    sources = [Matrix{Float64}(product_data[name].Y_sat) for name in names]
    cfg.fy4b_confirm_fusion || return sources
    position = something(findfirst(==("FY4B"), names))
    others = sources[setdiff(eachindex(sources), position)]
    sources[position] = confirmed_by_others(sources[position], others, cfg.mger.rain_threshold)
    return sources
end

"""
Everything `FY4B_QM_PRODUCT` needs that does not depend on the fold, or `nothing` when the
run does not score it: the offset grid aligned to the run's stations and hours, each hour's
regime and season-by-regime cell, and the correlation sums from which a fold picks its shift.
Registers the derived product's data, with an anchor that only a fold can fill.

Under `fy4b_qm_crossed` it also holds the 7-day time folds and, per time fold, the sums over that
fold's training hours alone (`crossed_shift_sums`); both are `nothing` otherwise.
"""
function _prepare_fy4b_calibration!(cfg::InterpolationBenchmarkConfig, ids, product_data)
    cfg.fy4b_qm_grid === nothing && return nothing
    source = product_data["FY4B"]
    regimes, regime_of_time = product_regimes(_fy4b_hour_meta_path(cfg), source.times)
    grid = open_offset_grid(something(cfg.fy4b_qm_grid), ids, source.times)
    sums = shift_sums(grid, source.Y_obs, regime_of_time, length(regimes))
    product_data[FY4B_QM_PRODUCT] = (;
        source.times, source.ids, source.Y_obs, Y_sat=fill(NaN, size(source.Y_obs)),
    )
    month_labels, month_of_time = strength_cells(source.times, regime_of_time)
    cells = calibration_cells(source.times, regime_of_time)
    time_folds = cfg.fy4b_qm_crossed ?
        time_block_folds(source.times, cells; seed=FY4B_QM_CROSSED_TIME_FOLD_SEED) : nothing
    crossed_sums = time_folds === nothing ? nothing : crossed_shift_sums(
        grid, source.Y_obs, regime_of_time, length(regimes), time_folds.fold_of_time, source.times)
    return (; grid, sums, regimes, regime_of_time, cells, month_labels, month_of_time,
        cache=Dict{Int,Matrix{Float64}}(), time_folds, crossed_sums)
end

"""
`sources` (named by `names`) as one fold feeds them to a fused anchor or the blend bands: as they
are, or under `fy4b_qm_fusion` with FY4B swapped for that fold's calibrated anchor.
"""
function _fold_sources(cfg::InterpolationBenchmarkConfig, names, sources, fold_calibration, fold)
    cfg.fy4b_qm_fusion || return sources
    swapped = copy(sources)
    swapped[something(findfirst(==("FY4B"), names))] = first(fold_calibration[fold])
    return swapped
end

"""
Each fold's calibrated FY4B, as `(anchor, shift_rows, map_rows, strength_rows)` per fold, or
`nothing` when the run does not score `FY4B_QM_PRODUCT`.

A fold's calibrated FY4B depends only on its training stations, so it is built once per fold of a
partition: for FY4B_QM itself and, under `fy4b_qm_fusion`, for the fused anchors and blend bands
it feeds. Under `fy4b_qm_crossed` each time fold's hours come from a fit on the other time folds'
hours alone.
"""
function _fold_calibrations(setup, folds, id_map)
    (; cfg, product_data, fy4b_calibration) = setup
    fy4b_calibration === nothing && return nothing
    return map(1:cfg.k) do fold
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
end

"""
Each fold's gated FY4B_QM, as `(anchor, cut_rows)` per fold, or `nothing` without
`fy4b_qm_rh_gate`: the cut chosen on that fold's training stations, per time fold on the other
blocks' hours, and applied at every station.
"""
function _fold_gates(setup, folds, id_map, fold_calibration)
    (; cfg, product_data, gate_inputs) = setup
    gate_inputs === nothing && return nothing
    return map(1:cfg.k) do fold
        _, _, train_idx, _ = _fold_station_indices(folds, id_map, cfg.k, fold)
        time_crossed_gate(
            first(fold_calibration[fold]), product_data[FY4B_QM_PRODUCT].Y_obs,
            gate_inputs.humidity, train_idx, gate_inputs.time_folds.fold_of_time,
            product_data[FY4B_QM_PRODUCT].times; wet=cfg.mger.rain_threshold,
        )
    end
end

"""
Each fold's two-stage FY4B, as `(anchor, table_rows)` per fold, or `nothing` without
`fy4b_twostage`: the tables counted on that fold's training stations at its pixel shift (the cache
it shares with the calibration), applied at every station.
"""
function _fold_two_stage(setup, folds, id_map, fold_calibration)
    (; cfg, product_data, fy4b_calibration, twostage_nearby) = setup
    cfg.fy4b_twostage || return nothing
    return map(1:cfg.k) do fold
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
end

"""
Each fold's blend bands under `fy4b_qm_fusion`, read off the sources with FY4B swapped for the
fold's calibrated anchor (`_fold_sources`). `nothing` otherwise: one set of bands, `setup`'s, then
serves every fold.
"""
function _fold_blend_bands(setup, fold_calibration)
    (; cfg, source_products, band_sources, blend_bands) = setup
    (cfg.fy4b_qm_fusion && !isempty(blend_bands)) || return nothing
    return map(1:cfg.k) do fold
        sources = _fold_sources(cfg, source_products, band_sources, fold_calibration, fold)
        Dict(axis => blend_band_matrix(
            axis, first(sources), sources, cfg.mger.rain_threshold,
        ) for axis in keys(blend_bands))
    end
end

"""
The anchor a derived product is corrected from in one fold, at every station, with the rows that
say how it was built appended to `rows`.

There are five kinds: FY4B's three per-fold series, which `partition` already holds
(`FY4B_QM_PRODUCT`, `FY4B_QM_RH_PRODUCT`, `FY4B_TS_PRODUCT`), and the two fusions, grouped or not,
whose coefficients are fitted here. Every kind is fitted on the fold's training stations
(`train_idx`) and applied everywhere: the fold's own training rows need the same anchor its
held-out rows get, or the residual target would mean two different things inside one fit.
"""
function _fold_anchor!(rows, setup, partition, cell, y_obs, train_idx)
    (; cfg, fused_specifications) = setup
    (; product, fold) = cell
    context = merge(cell, (; n_train_station=length(train_idx)))
    if product == FY4B_QM_PRODUCT
        fold_sat, shift_rows, map_rows, strength_rows = partition.fold_calibration[fold]
        append!(rows.fy4b_shift, [merge(context, row) for row in shift_rows])
        append!(rows.fy4b_map, [merge(context, row) for row in map_rows])
        append!(rows.fy4b_strength, [merge(context, row) for row in strength_rows])
        return fold_sat
    elseif product == FY4B_QM_RH_PRODUCT
        fold_sat, cut_rows = partition.fold_gate[fold]
        append!(rows.fy4b_gate, [merge(context, row) for row in cut_rows])
        return fold_sat
    elseif product == FY4B_TS_PRODUCT
        fold_sat, table_rows = partition.fold_twostage[fold]
        append!(rows.fy4b_twostage, [merge(context, row) for row in table_rows])
        return fold_sat
    end
    fusion = fused_specifications[product]
    if fusion.design !== nothing
        design = fusion.design
        coefficients, used, fell_back = fit_grouped_fusion(
            y_obs, fusion.sources, design.features, design.groups,
            design.n_group, train_idx,
        )
        coefficients === nothing && error(
            "$product fold $fold: the pooled fusion fit is singular",
        )
        fold_sat = apply_grouped_fusion(
            fusion.sources, design.features, coefficients, design.groups,
        )
        lower = lowercase.(fusion.source_names)
        terms = vcat("intercept", lower, lower .* "_lag1",
            lower .* "_lead1", lower .* "_nbr8")
        for group in 1:design.n_group, (row, term) in enumerate(terms)
            push!(rows.fused_grouped, merge(cell, (;
                variant=String(fusion.variant), group=group - 1,
                n_cell=used[group], fell_back=fell_back[group], term,
                coefficient=coefficients[row, group],
            )))
        end
        return fold_sat
    end
    sources = _fold_sources(
        cfg, fusion.source_names, fusion.sources, partition.fold_calibration, fold,
    )
    coefficients, used, fell_back = fusion_coefficients(
        fusion.variant, y_obs, sources, train_idx,
    )
    fold_sat = apply_satellite_fusion(sources, coefficients)
    push!(rows.fused_anchor, merge(
        cell,
        (; variant=String(fusion.variant),
            n_train_station=length(train_idx), n_train_cell=used,
            fell_back, intercept=coefficients[1]),
        NamedTuple{Tuple(Symbol("beta_", lowercase(name))
                         for name in fusion.source_names)}(
            Tuple(coefficients[2:end])),
    ))
    return fold_sat
end

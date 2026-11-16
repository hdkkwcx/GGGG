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

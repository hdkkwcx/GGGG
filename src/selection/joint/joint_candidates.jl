function _candidate_rows(
    product::String, Yobs::Matrix{Float64}, Ysat::Matrix{Float64},
    terrain::DataFrame, era5::AbstractDict, ndvi_aligned,
    train_indices::Vector{Int}, station_ids::Vector{String},
    times::Vector{DateTime}, cfg::JointSelectionConfig,
    run_seed::Int, dem_seed::Int,
)
    rows = NamedTuple[]
    active = String[]
    # DEM screening uses station mean wet residuals and grouped terrain VIF.
    response, _ = mean_wet_residual(
        Yobs[train_indices, :], Ysat[train_indices, :];
        threshold=cfg.wet_threshold, min_hours=cfg.min_wet_hours,
    )
    dem_screen, _, dem_selected = terrain_screen(
        terrain[train_indices, :], response; product,
        permutations=cfg.independent_permutations,
        q_threshold=cfg.q_threshold, vif_threshold=cfg.vif_threshold,
        seed=dem_seed,
    )
    for group in ("elevation", "slope", "aspect")
        row = only(eachrow(filter(:variable_group => ==(group), dem_screen)))
        selected = group in dem_selected
        push!(rows, (;
            family="dem", variable_group=group,
            predictor_columns=_predictor_label(group),
            independent_qvalue=Float64(row.qvalue),
            direction_stable=Bool(row.direction_stable),
            independent_selected=selected,
            independent_status=selected ? "selected" : String(row.exclusion_reason),
        ))
        selected && push!(active, group)
    end

    era_cfg, ndvi_cfg = _family_configs(cfg)
    era_panel = prepare_dynamic_panel(Yobs, Ysat, era5, train_indices, times, era_cfg)
    if era_panel.qc.status[1] == "ok"
        era_screen = dynamic_panel_screen(
            era_panel, times, era_cfg; rng=MersenneTwister(20260816 + run_seed),
        )
        for variable in ERA5_VARIABLES
            row = only(eachrow(filter(:variable => ==(String(variable)), era_screen.association)))
            selected = variable in era_screen.selected
            push!(rows, (;
                family="era5", variable_group=String(variable),
                predictor_columns=String(variable), independent_qvalue=Float64(row.qvalue),
                direction_stable=Bool(row.direction_stable), independent_selected=selected,
                independent_status=selected ? "selected" :
                    (!row.direction_stable ? "direction_conflict" : "q_not_significant_or_family_vif"),
            ))
            selected && push!(active, String(variable))
        end
    else
        for variable in ERA5_VARIABLES
            push!(rows, (;
                family="era5", variable_group=String(variable), predictor_columns=String(variable),
                independent_qvalue=NaN, direction_stable=false, independent_selected=false,
                independent_status="insufficient_stations",
            ))
        end
    end

    ndvi_panel = prepare_ndvi_panel(
        Yobs, Ysat, ndvi_aligned, train_indices,
        station_ids, ndvi_cfg,
    )
    if ndvi_panel.qc.status[1] == "ok"
        ndvi_screen = screen_ndvi(
            ndvi_panel, ndvi_aligned.periods, ndvi_cfg;
            rng=MersenneTwister(20260817 + run_seed),
        )
        row = ndvi_screen.association[1, :]
        selected = !isempty(ndvi_screen.selected)
        push!(rows, (;
            family="ndvi", variable_group="ndvi", predictor_columns="ndvi_qc",
            independent_qvalue=Float64(row.qvalue),
            direction_stable=Bool(row.direction_stable), independent_selected=selected,
            independent_status=selected ? "selected" :
                (!row.direction_stable ? "direction_conflict" : "q_not_significant"),
        ))
        selected && push!(active, "ndvi")
    else
        push!(rows, (;
            family="ndvi", variable_group="ndvi", predictor_columns="ndvi_qc",
            independent_qvalue=NaN, direction_stable=false, independent_selected=false,
            independent_status="insufficient_stations_or_variation",
        ))
    end
    active = [group for group in JOINT_GROUP_ORDER if group in Set(active)]
    return DataFrame(rows), active
end

"""
Select joint covariates for one (product, `train_indices`) cell: independent per-family
screening → cross-family joint VIF pruning → local/global spatial-role test. Used both by
`run_joint_variable_selection`'s own full-data/stability reporting and, per training fold, by
the interpolation benchmark's nested selection mode — the only difference between the two
callers is which `train_indices` they pass in.
"""
function select_joint_covariates(
    product::String, Yobs::Matrix{Float64}, Ysat::Matrix{Float64}, terrain::DataFrame,
    era5::AbstractDict, ndvi_aligned, train_indices::Vector{Int}, station_ids::Vector{String},
    times::Vector{DateTime}, lonlat::Matrix{Float64}, cfg::JointSelectionConfig,
    run_seed::Int, dem_seed::Int,
)
    candidates, active = _candidate_rows(
        product, Yobs, Ysat, terrain, era5, ndvi_aligned,
        train_indices, station_ids, times, cfg, run_seed, dem_seed,
    )
    if isempty(active)
        qc = DataFrame(
            training_station_count=[length(train_indices)], eligible_station_count=[0],
            eligible_time_count=[0], used_station_hours=[0], group_count=[0],
            predictor_column_count=[0], uses_ndvi_hierarchical_weight=[false],
            minimum_station_weight=[NaN], maximum_station_weight=[NaN],
            status=["no_independent_candidates"],
        )
        return (;
            candidates, active, qc, scales=DataFrame(), weight_audit=DataFrame(),
            vif=DataFrame(), retained=String[], spatial_scan=DataFrame(),
            spatial_result=DataFrame(), role_map=Dict{String,String}(),
            status="ok", status_reason="no_independent_candidates",
        )
    end
    independent_q = Dict(String(row.variable_group) => Float64(row.independent_qvalue)
        for row in eachrow(candidates))
    panel = prepare_joint_panel(Yobs, Ysat, terrain, era5, ndvi_aligned, train_indices, active, cfg)
    if panel.qc.status[1] != "ok"
        retained = String[]; vif = DataFrame()
    else
        vif, retained = joint_vif(panel, active, independent_q, cfg)
    end
    spatial = panel.qc.status[1] == "ok" ? joint_spatial_variability_test(
        panel, retained, lonlat[train_indices, :], cfg; rng=MersenneTwister(cfg.seed + run_seed),
    ) : joint_spatial_variability_test(panel, String[], lonlat[train_indices, :], cfg)
    role_map = Dict(String(row.variable_group) => String(row.role)
        for row in eachrow(spatial.result) if row.role in ("local", "global"))
    status = panel.qc.status[1] == "ok" ? "ok" : "failed"
    status_reason = status == "ok" ? (isempty(retained) ? "all_removed_by_joint_vif" : "") :
        String(panel.qc.status[1])
    return (;
        candidates, active, qc=panel.qc, scales=panel.scales, weight_audit=panel.weight_audit,
        vif, retained, spatial_scan=spatial.scan, spatial_result=spatial.result,
        role_map, status, status_reason,
    )
end

function run_joint_variable_selection(
    cfg::JointSelectionConfig, products::Vector{String}, station_ids::Vector{String},
    times::Vector{DateTime}, Yobs::Matrix{Float64}, satellite::AbstractDict,
    terrain::DataFrame, era5::AbstractDict, ndvi_table::DataFrame,
    lonlat::Matrix{Float64}; prerequisite_audit::Union{Nothing,DataFrame}=nothing,
)
    mkpath(cfg.outdir)
    n = length(station_ids)
    length(times) == size(Yobs, 2) || throw(DimensionMismatch("observation time dimension differs"))
    size(Yobs, 1) == n || throw(DimensionMismatch("observation station dimension differs"))
    nrow(terrain) == n || throw(DimensionMismatch("terrain station dimension differs"))
    String.(terrain.station_id) == station_ids ||
        throw(ArgumentError("terrain rows must follow common station order"))
    size(lonlat) == (n, 2) || throw(DimensionMismatch("coordinate dimensions differ"))
    all(size(satellite[product]) == size(Yobs) for product in products) ||
        throw(DimensionMismatch("satellite dimensions differ"))
    all(size(era5[variable]) == size(Yobs) for variable in ERA5_VARIABLES) ||
        throw(DimensionMismatch("ERA5 dimensions differ"))
    folds = balanced_spatial_folds(station_ids, lonlat; k=cfg.k)
    CSV.write(joinpath(cfg.outdir, "spatial_folds.csv"), DataFrame(
        station_id=station_ids, fold=folds, lon=lonlat[:, 1], lat=lonlat[:, 2],
    ))
    prerequisite_audit !== nothing && CSV.write(
        joinpath(cfg.outdir, "prerequisite_full_audit.csv"), prerequisite_audit,
    )
    ndvi_aligned = align_ndvi_asof(
        ndvi_table, station_ids, times; max_age_days=cfg.max_ndvi_age_days,
    )
    candidate_all, qc_all, scale_all, weight_all, vif_all =
        DataFrame(), DataFrame(), DataFrame(), DataFrame(), DataFrame()
    bandwidth_all, spatial_all, role_all, status_all = DataFrame(), DataFrame(), DataFrame(), DataFrame()
    final_spec = DataFrame(
        product=String[], family=String[], variable_group=String[], predictor_columns=String[],
        independent_selected=Bool[], independent_qvalue=Float64[], joint_vif_retained=Bool[],
        role=String[], final_included=Bool[], exclusion_reason=String[],
    )
    schemes = selection_schemes(n, folds, cfg.k)
    for (product_index, product) in enumerate(products)
        Ysat = Float64.(satellite[product])
        for (scheme, fold, train_indices) in schemes
            run_seed = 1000 * product_index + 10 * fold
            dem_seed = 20260815 + Int(sum(codeunits(product))) +
                (fold == 0 ? 0 : 10_000 + 10 * fold)
            selection = select_joint_covariates(
                product, Yobs, Ysat, terrain, era5, ndvi_aligned,
                train_indices, station_ids, times, lonlat, cfg, run_seed, dem_seed,
            )
            candidates = selection.candidates
            annotate_selection!(candidates, product, scheme, fold); append_selection!(candidate_all, candidates)
            if isempty(selection.active)
                qc = selection.qc
                annotate_selection!(qc, product, scheme, fold); append_selection!(qc_all, qc)
                push!(status_all, (product=product, scheme=scheme, fold=fold,
                    status="ok", reason="no_independent_candidates"))
                for row in eachrow(candidates)
                    push!(role_all, (product=product, scheme=scheme, fold=fold,
                        family=row.family, variable_group=row.variable_group,
                        independent_selected=false, joint_vif_retained=false,
                        role="not_selected", final_included=false,
                        exclusion_reason=row.independent_status))
                    if scheme == "full_data"
                        push!(final_spec, (product, row.family, row.variable_group,
                            row.predictor_columns, false, row.independent_qvalue, false,
                            "not_selected", false, row.independent_status))
                    end
                end
                continue
            end
            qc = copy(selection.qc); annotate_selection!(qc, product, scheme, fold); append_selection!(qc_all, qc)
            scales = copy(selection.scales)
            annotate_selection!(scales, product, scheme, fold); append_selection!(scale_all, scales)
            weight_audit = copy(selection.weight_audit)
            weight_audit.station_id = station_ids[weight_audit.global_station_index]
            annotate_selection!(weight_audit, product, scheme, fold); append_selection!(weight_all, weight_audit)
            vif = selection.vif; retained = selection.retained
            annotate_selection!(vif, product, scheme, fold); append_selection!(vif_all, vif)
            annotate_selection!(selection.spatial_scan, product, scheme, fold)
            append_selection!(bandwidth_all, selection.spatial_scan)
            annotate_selection!(selection.spatial_result, product, scheme, fold)
            append_selection!(spatial_all, selection.spatial_result)
            role_status_map = Dict(String(row.variable_group) => (role=String(row.role), status=String(row.status))
                for row in eachrow(selection.spatial_result))
            for row in eachrow(candidates)
                group = String(row.variable_group)
                selected = Bool(row.independent_selected)
                vif_retained = group in retained
                if !selected
                    role, included, reason = "not_selected", false, String(row.independent_status)
                elseif !vif_retained
                    role, included, reason = "not_selected", false,
                        selection.status == "ok" ? "removed_by_joint_vif" : selection.status_reason
                else
                    info = get(role_status_map, group, (role="uncertain", status="test_failed"))
                    role = info.role; included = role in ("local", "global")
                    reason = included ? "" : info.status
                end
                push!(role_all, (product=product, scheme=scheme, fold=fold,
                    family=row.family, variable_group=group, independent_selected=selected,
                    joint_vif_retained=vif_retained, role, final_included=included,
                    exclusion_reason=reason))
                if scheme == "full_data"
                    push!(final_spec, (product, String(row.family), group,
                        String(row.predictor_columns), selected, Float64(row.independent_qvalue),
                        vif_retained, role, included, reason))
                end
            end
            push!(status_all, (product=product, scheme=scheme, fold=fold,
                status=selection.status, reason=selection.status_reason))
        end
    end
    stability = DataFrame(
        product=String[], family=String[], variable_group=String[], fold_count=Int[],
        independent_selected_count=Int[], vif_retained_count=Int[],
        local_count=Int[], global_count=Int[], uncertain_count=Int[], final_included_count=Int[],
    )
    for group in groupby(filter(:scheme => ==("spatial_cv"), role_all),
        [:product, :family, :variable_group])
        push!(stability, (
            first(group.product), first(group.family), first(group.variable_group), nrow(group),
            count(group.independent_selected), count(group.joint_vif_retained),
            count(==("local"), group.role), count(==("global"), group.role),
            count(==("uncertain"), group.role), count(group.final_included),
        ))
    end
    consensus = DataFrame(
        family=String[], variable_group=String[], included_product_count=Int[],
        local_product_count=Int[], global_product_count=Int[],
        consensus_included=Bool[], consensus_role=String[],
    )
    for variable_group in JOINT_GROUP_ORDER
        rows = filter(:variable_group => ==(variable_group), final_spec)
        included_count = count(rows.final_included)
        local_count = count(==("local"), rows.role); global_count = count(==("global"), rows.role)
        role = included_count == 0 ? "not_selected" : local_count >= 2 ? "local" :
            global_count >= 2 ? "global" : included_count >= 2 ? "role_unstable" : "product_specific"
        push!(consensus, (GROUP_FAMILY[variable_group], variable_group, included_count,
            local_count, global_count, included_count >= 2, role))
    end
    outputs = Dict(
        "joint_independent_candidates.csv" => candidate_all,
        "joint_panel_quality_control.csv" => qc_all,
        "joint_scaling.csv" => scale_all,
        "joint_weight_audit.csv" => weight_all,
        "joint_vif.csv" => vif_all,
        "joint_bandwidth_scan.csv" => bandwidth_all,
        "joint_spatial_variability.csv" => spatial_all,
        "joint_fold_roles.csv" => role_all,
        "joint_role_stability.csv" => stability,
        "joint_final_full_data_spec.csv" => final_spec,
        "joint_cross_product_consensus.csv" => consensus,
        "run_status.csv" => status_all,
    )
    for (filename, table) in outputs
        CSV.write(joinpath(cfg.outdir, filename), table)
    end
    return (; final_spec, stability, consensus, status=status_all)
end

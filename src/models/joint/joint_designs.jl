function joint_effective_roles(context::JointFoldContext, method::String)
    method == "residual_gwr" && return Dict(group => "local" for group in context.variables)
    method in ("mixed_gwr", "mgwr") || throw(ArgumentError("unsupported joint model $method"))
    return copy(context.roles)
end

"""
Whether two joint models reduce to the same fit on this fold's covariate roles.

`residual_gwr` forces every covariate local; `mixed_gwr` honours the fold's role map. When that
map contains no `"global"` the two build the identical design — `_design_at` then differs only in
the group *name* (`"shared_all_local"` vs `"shared_local"`) — and they search the identical grid
through the identical scorer, so they are one model reported under two names. On the 2026-09-02
full nested run that held at 19 of 30 fold-cells, with bit-identical selected kernel, bandwidth,
shrinkage and inner RMSE; `joint_role_stability.csv` shows why, a `"global"` role appearing in
only 6 of 54 group-cells.

`mgwr` is never a duplicate of either: it gives every group its own bandwidth whatever the roles
say, and under `:intercept_only` it also drops the coordinate columns.

This is what `run_status.csv`'s `duplicate_of` column and `select_auto_method`'s contender
collapsing are keyed on, so a reader does not have to infer "same model" from matching RMSEs.
"""
function joint_models_coincide(
    context::JointFoldContext, first_method::AbstractString, second_method::AbstractString,
)
    first_method == second_method && return true
    ("mgwr" in (first_method, second_method)) && return false
    return joint_effective_roles(context, String(first_method)) ==
        joint_effective_roles(context, String(second_method))
end

function joint_group_names(context::JointFoldContext, method::String)
    design = _design_at(context, method, 1; target=false)
    return design.group_names
end

function _design_at(context::JointFoldContext, method::String, time::Int; target::Bool=false)
    spatial = target ? context.spatial_target : context.spatial_train
    predictors = target ? context.predictor_target : context.predictor_train
    roles = joint_effective_roles(context, method)
    local_groups = Matrix{Float64}[]
    group_names = String[]
    if method == "mgwr"
        # `mgwr_spatial_grouping` controls only the spatial block; covariates always get their
        # own group, and therefore their own bandwidth, in every variant.
        grouping = context.config.mgwr_spatial_grouping
        if grouping === :split
            for index in axes(spatial, 2)
                push!(local_groups, spatial[:, index:index])
                push!(group_names, ("intercept", "longitude", "latitude")[index])
            end
        elseif grouping === :shared
            push!(local_groups, copy(spatial))
            # Deliberately not "shared_local" (mixed_gwr's name) or "all" (`_scan_row`'s default
            # for ungrouped methods): the scan-row lookups match on the group name.
            push!(group_names, "shared_spatial")
        elseif grouping === :intercept_only
            push!(local_groups, spatial[:, 1:1])
            push!(group_names, "intercept")
        else
            throw(ArgumentError("unsupported mgwr_spatial_grouping: $grouping"))
        end
        for group in context.variables
            roles[group] == "local" || continue
            push!(local_groups, Matrix{Float64}(predictors[group][:, :, time]))
            push!(group_names, group)
        end
    else
        local_design = copy(spatial)
        for group in context.variables
            roles[group] == "local" || continue
            local_design = hcat(local_design, predictors[group][:, :, time])
        end
        push!(local_groups, local_design)
        push!(group_names, method == "residual_gwr" ? "shared_all_local" : "shared_local")
    end
    global_design = zeros(Float64, size(spatial, 1), 0)
    for group in context.variables
        roles[group] == "global" || continue
        global_design = hcat(global_design, predictors[group][:, :, time])
    end
    return (; local_groups, global_design, group_names)
end

function _weighted_rows(panel, active::Vector{String})
    columns = Symbol[]; column_groups = String[]
    for group in active, column in GROUP_COLUMNS[group]
        push!(columns, column); push!(column_groups, group)
    end
    rows = Vector{Vector{Float64}}(); weights = Float64[]
    for i in axes(panel.y, 1), time in axes(panel.y, 2)
        panel.weights[i, time] > 0 && isfinite(panel.y[i, time]) &&
            all(isfinite(panel.x[column][i, time]) for column in columns) || continue
        push!(rows, [panel.x[column][i, time] for column in columns])
        push!(weights, panel.weights[i, time])
    end
    X = isempty(rows) ? zeros(0, length(columns)) : reduce(vcat, permutedims.(rows))
    return X, weights, columns, column_groups
end

function _column_vifs(X::Matrix{Float64}, weights::Vector{Float64})
    p = size(X, 2)
    p == 0 && return Float64[]
    p == 1 && return [1.0]
    sw = sqrt.(weights ./ sum(weights))
    result = fill(Inf, p)
    for j in 1:p
        others = setdiff(1:p, [j])
        design = hcat(ones(size(X, 1)), X[:, others])
        beta = (design .* sw) \ (X[:, j] .* sw)
        prediction = design * beta
        mean_value = dot(X[:, j], weights) / sum(weights)
        rss = sum(weights .* (X[:, j] .- prediction) .^ 2)
        tss = sum(weights .* (X[:, j] .- mean_value) .^ 2)
        r2 = tss > 0 ? clamp(1 - rss / tss, 0.0, 1.0) : 1.0
        result[j] = r2 < 1 ? 1 / (1 - r2) : Inf
    end
    return result
end

function joint_vif(
    panel, candidates::Vector{String}, independent_q::Dict{String,Float64},
    cfg::JointSelectionConfig,
)
    active = copy(candidates)
    audit = DataFrame(
        iteration=Int[], variable_group=String[], predictor_columns=String[],
        group_vif=Float64[], removed=Bool[], reason=String[],
    )
    iteration = 1
    while !isempty(active)
        X, weights, columns, column_groups = _weighted_rows(panel, active)
        column_vifs = _column_vifs(X, weights)
        group_values = Dict(group => maximum(column_vifs[column_groups .== group]) for group in active)
        maximum_value = maximum(Base.values(group_values))
        if maximum_value < cfg.vif_threshold
            for group in active
                push!(audit, (iteration, group, join(String.(GROUP_COLUMNS[group]), "+"),
                    group_values[group], false, "retained"))
            end
            break
        end
        tied = [group for group in active if isapprox(group_values[group], maximum_value;
            rtol=1e-10, atol=1e-10) || (isinf(group_values[group]) && isinf(maximum_value))]
        remove_group = sort(tied; by=group ->
            (-get(independent_q, group, Inf), findfirst(==(group), JOINT_GROUP_ORDER)))[1]
        for group in active
            push!(audit, (iteration, group, join(String.(GROUP_COLUMNS[group]), "+"),
                group_values[group], group == remove_group,
                group == remove_group ? "joint_vif_ge_threshold" : "pending"))
        end
        filter!(!=(remove_group), active)
        iteration += 1
    end
    return audit, active
end

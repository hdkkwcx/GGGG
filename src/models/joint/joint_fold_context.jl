function _raw_group(
    group::String, indices::Vector{Int}, nt::Int, terrain::DataFrame,
    era5::AbstractDict, ndvi_aligned,
)
    columns = GROUP_COLUMNS[group]
    values = Array{Float64}(undef, length(indices), length(columns), nt)
    for (column_index, column) in enumerate(columns)
        if GROUP_FAMILY[group] == "dem"
            values[:, column_index, :] .= reshape(
                Float64.(terrain[indices, column]), length(indices), 1,
            )
        elseif GROUP_FAMILY[group] == "era5"
            values[:, column_index, :] .= Float64.(era5[column][indices, :])
        else
            ndvi_aligned === nothing && throw(ArgumentError("NDVI is selected but not loaded"))
            values[:, column_index, :] .= Float64.(ndvi_aligned.values[indices, :])
        end
    end
    return values
end

function _spatial_design(train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64})
    center = vec(mean(train_lonlat, dims=1))
    scale = vec(std(train_lonlat, dims=1))
    all(isfinite, scale) && all(>(0), scale) ||
        throw(ArgumentError("training coordinates have insufficient variation"))
    train = hcat(ones(size(train_lonlat, 1)),
        (train_lonlat[:, 1] .- center[1]) ./ scale[1],
        (train_lonlat[:, 2] .- center[2]) ./ scale[2])
    target = hcat(ones(size(target_lonlat, 1)),
        (target_lonlat[:, 1] .- center[1]) ./ scale[1],
        (target_lonlat[:, 2] .- center[2]) ./ scale[2])
    return train, target
end

function _hierarchical_weights(mask::BitMatrix, train_indices::Vector{Int}, ndvi_aligned)
    n, nt = size(mask)
    weights = zeros(n, nt)
    uses_ndvi = ndvi_aligned !== nothing
    for i in 1:n
        valid_times = findall(mask[i, :])
        isempty(valid_times) && continue
        if uses_ndvi
            periods = ndvi_aligned.source_period[train_indices[i], valid_times]
            valid_periods = sort(unique(filter(>(0), periods)))
            isempty(valid_periods) && continue
            for period in valid_periods
                period_times = valid_times[periods .== period]
                weights[i, period_times] .= 1 / (length(valid_periods) * length(period_times))
            end
        else
            weights[i, valid_times] .= 1 / length(valid_times)
        end
    end
    return weights
end

function build_joint_fold_context(
    product::String, roles::Dict{String,String}, train_indices::Vector{Int},
    target_indices::Vector{Int}, lonlat::Matrix{Float64}, Yobs::Matrix{Float64},
    Ysat::Matrix{Float64}, terrain::DataFrame, era5::AbstractDict, ndvi_aligned,
    cfg::JointCovariateBenchmarkConfig,
)
    variables = [group for group in JOINT_GROUP_ORDER if haskey(roles, group)]
    nt = size(Yobs, 2)
    train_lonlat = Matrix{Float64}(lonlat[train_indices, :])
    target_lonlat = Matrix{Float64}(lonlat[target_indices, :])
    spatial_train, spatial_target = _spatial_design(train_lonlat, target_lonlat)
    raw_train = Dict(group => _raw_group(
        group, train_indices, nt, terrain, era5,
        group == "ndvi" ? ndvi_aligned : nothing,
    ) for group in variables)
    raw_target = Dict(group => _raw_group(
        group, target_indices, nt, terrain, era5,
        group == "ndvi" ? ndvi_aligned : nothing,
    ) for group in variables)

    mask = isfinite.(Yobs[train_indices, :]) .& isfinite.(Ysat[train_indices, :]) .&
        (Yobs[train_indices, :] .>= cfg.wet_threshold)
    for group in variables, column in axes(raw_train[group], 2)
        mask .&= isfinite.(raw_train[group][:, column, :])
    end
    weights = _hierarchical_weights(
        mask, train_indices, "ndvi" in variables ? ndvi_aligned : nothing,
    )
    predictor_train = Dict{String,Array{Float64,3}}()
    predictor_target = Dict{String,Array{Float64,3}}()
    scaling_rows = NamedTuple[]
    for group in variables
        train = copy(raw_train[group]); target = copy(raw_target[group])
        for column in axes(train, 2)
            for time in 1:nt
                finite_train = isfinite.(@view train[:, column, time])
                if any(finite_train)
                    hour_mean = mean(@view train[finite_train, column, time])
                    train[finite_train, column, time] .-= hour_mean
                    finite_target = isfinite.(@view target[:, column, time])
                    target[finite_target, column, time] .-= hour_mean
                end
            end
            weighted_values = Float64[]; value_weights = Float64[]
            for i in axes(train, 1), time in 1:nt
                weight = weights[i, time]
                weight > 0 && isfinite(train[i, column, time]) || continue
                push!(weighted_values, train[i, column, time])
                push!(value_weights, weight)
            end
            total = sum(value_weights)
            center = total > 0 ? dot(weighted_values, value_weights) / total : NaN
            scale = total > 0 ? sqrt(dot(
                (weighted_values .- center) .^ 2, value_weights,
            ) / total) : NaN
            isfinite(scale) && scale > 0 || throw(ArgumentError(
                "training-fold scale is unavailable for $product/$group column $column",
            ))
            for index in eachindex(train)
                isfinite(train[index]) && (train[index] = (train[index] - center) / scale)
            end
            for index in eachindex(target)
                isfinite(target[index]) && (target[index] = (target[index] - center) / scale)
            end
            push!(scaling_rows, (;
                product, variable_group=group,
                predictor_column=String(GROUP_COLUMNS[group][column]),
                center, scale, wet_station_hours=count(weights .> 0),
            ))
        end
        predictor_train[group] = train
        predictor_target[group] = target
    end
    valid_station_weights = vec(sum(weights, dims=2))
    positive = valid_station_weights[valid_station_weights .> 0]
    qc = DataFrame([(
        product=product, training_station_count=length(train_indices),
        target_station_count=length(target_indices), time_count=nt,
        variable_count=length(variables), variables=join(variables, ","),
        roles=join(["$group=$(roles[group])" for group in variables], ";"),
        wet_station_hours=count(mask),
        minimum_station_weight=isempty(positive) ? NaN : minimum(positive),
        maximum_station_weight=isempty(positive) ? NaN : maximum(positive),
        status="ok",
    )])
    return JointFoldContext(
        product, variables, copy(roles), spatial_train, spatial_target,
        predictor_train, predictor_target, train_lonlat, target_lonlat,
        filter(<(length(train_indices)), unique(cfg.bandwidth_candidates)),
        DataFrame(scaling_rows), qc, cfg,
    )
end

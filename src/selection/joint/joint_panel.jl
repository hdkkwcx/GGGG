function _raw_predictors(
    groups::Vector{String}, terrain::DataFrame, era5::AbstractDict,
    ndvi_aligned, train_indices::Vector{Int}, time_count::Int,
)
    predictors = Dict{Symbol,Matrix{Float64}}()
    column_groups = String[]
    columns = Symbol[]
    for group in groups
        for column in GROUP_COLUMNS[group]
            matrix = if GROUP_FAMILY[group] == "dem"
                repeat(reshape(Float64.(terrain[train_indices, column]), :, 1), 1, time_count)
            elseif GROUP_FAMILY[group] == "era5"
                Float64.(era5[column][train_indices, :])
            else
                Float64.(ndvi_aligned.values[train_indices, :])
            end
            predictors[column] = matrix
            push!(columns, column); push!(column_groups, group)
        end
    end
    return predictors, columns, column_groups
end

function prepare_joint_panel(
    Yobs::Matrix{Float64}, Ysat::Matrix{Float64}, terrain::DataFrame,
    era5::AbstractDict, ndvi_aligned, train_indices::Vector{Int},
    groups::Vector{String}, cfg::JointSelectionConfig,
)
    n, nt = length(train_indices), size(Yobs, 2)
    predictors, columns, column_groups = _raw_predictors(
        groups, terrain, era5, ndvi_aligned, train_indices, nt,
    )
    residual = Yobs[train_indices, :] .- Ysat[train_indices, :]
    mask = falses(n, nt)
    for i in 1:n, time in 1:nt
        global_station = train_indices[i]
        mask[i, time] = isfinite(Yobs[global_station, time]) &&
            isfinite(Ysat[global_station, time]) &&
            Yobs[global_station, time] >= cfg.wet_threshold &&
            all(isfinite(predictors[column][i, time]) for column in columns)
    end
    previous = -1
    while count(mask) != previous
        previous = count(mask)
        if "ndvi" in groups
            for i in 1:n
                global_station = train_indices[i]
                period_values = ndvi_aligned.source_period[global_station, :]
                for period in unique(filter(>(0), period_values[mask[i, :]]))
                    period_times = findall(mask[i, :] .& (period_values .== period))
                    length(period_times) >= cfg.min_ndvi_wet_hours_per_period ||
                        (mask[i, period_times] .= false)
                end
                retained_periods = unique(filter(>(0), period_values[mask[i, :]]))
                length(retained_periods) >= cfg.min_ndvi_periods_per_station ||
                    (mask[i, :] .= false)
            end
        end
        station_ok = vec(sum(mask, dims=2)) .>= cfg.min_wet_hours
        mask[.!station_ok, :] .= false
        time_ok = vec(sum(mask, dims=1)) .>= cfg.min_stations_per_time
        mask[:, .!time_ok] .= false
    end
    eligible_station = vec(sum(mask, dims=2)) .>= cfg.min_wet_hours
    y = fill(NaN, n, nt)
    x = Dict(column => fill(NaN, n, nt) for column in columns)
    for time in 1:nt
        stations = findall(mask[:, time])
        length(stations) >= cfg.min_stations_per_time || continue
        y[stations, time] .= residual[stations, time] .- mean(residual[stations, time])
        for column in columns
            values = predictors[column][stations, time]
            x[column][stations, time] .= values .- mean(values)
        end
    end
    weights = zeros(n, nt)
    uses_ndvi = "ndvi" in groups
    for i in 1:n
        valid_times = findall(mask[i, :])
        isempty(valid_times) && continue
        if uses_ndvi
            global_station = train_indices[i]
            period_values = ndvi_aligned.source_period[global_station, valid_times]
            valid_periods = sort(unique(filter(>(0), period_values)))
            isempty(valid_periods) && continue
            for period in valid_periods
                period_times = valid_times[period_values .== period]
                weights[i, period_times] .= 1 / (length(valid_periods) * length(period_times))
            end
        else
            weights[i, valid_times] .= 1 / length(valid_times)
        end
    end
    scale_rows = NamedTuple[]
    for column in columns
        values = Float64[]; value_weights = Float64[]
        for index in eachindex(y)
            isfinite(y[index]) && isfinite(x[column][index]) && weights[index] > 0 || continue
            push!(values, x[column][index]); push!(value_weights, weights[index])
        end
        total = sum(value_weights)
        mean_value = total > 0 ? dot(values, value_weights) / total : NaN
        scale = total > 0 ? sqrt(dot((values .- mean_value) .^ 2, value_weights) / total) : NaN
        if isfinite(scale) && scale > 0
            for index in eachindex(x[column])
                isfinite(x[column][index]) &&
                    (x[column][index] = (x[column][index] - mean_value) / scale)
            end
        end
        push!(scale_rows, (; column=String(column), mean=mean_value, scale))
    end
    station_weight_sums = vec(sum(weights, dims=2))
    positive = station_weight_sums[station_weight_sums .> 0]
    qc = DataFrame(
        training_station_count=[n], eligible_station_count=[count(eligible_station)],
        eligible_time_count=[count(vec(sum(mask, dims=1)) .>= cfg.min_stations_per_time)],
        used_station_hours=[count(mask)], group_count=[length(groups)],
        predictor_column_count=[length(columns)], uses_ndvi_hierarchical_weight=[uses_ndvi],
        minimum_station_weight=[isempty(positive) ? NaN : minimum(positive)],
        maximum_station_weight=[isempty(positive) ? NaN : maximum(positive)],
        status=[count(eligible_station) >= cfg.min_stations_per_time &&
            all(isfinite(row.scale) && row.scale > 0 for row in scale_rows) ? "ok" :
            "insufficient_stations_or_variation"],
    )
    weight_rows = NamedTuple[]
    for i in 1:n
        global_station = train_indices[i]
        valid_times = findall(weights[i, :] .> 0)
        period_weights = Float64[]
        period_count = 0
        if uses_ndvi && !isempty(valid_times)
            periods = ndvi_aligned.source_period[global_station, valid_times]
            valid_periods = sort(unique(filter(>(0), periods)))
            period_count = length(valid_periods)
            period_weights = [sum(weights[i, valid_times[periods .== period]])
                for period in valid_periods]
        end
        push!(weight_rows, (;
            local_station_index=i, global_station_index=global_station,
            eligible=eligible_station[i], used_hour_count=length(valid_times),
            ndvi_period_count=period_count, total_weight=station_weight_sums[i],
            minimum_period_weight=isempty(period_weights) ? NaN : minimum(period_weights),
            maximum_period_weight=isempty(period_weights) ? NaN : maximum(period_weights),
        ))
    end
    return (; y, x, weights, mask, eligible_station, columns, column_groups,
        groups, scales=DataFrame(scale_rows), weight_audit=DataFrame(weight_rows), qc)
end

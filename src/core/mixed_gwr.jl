"""
Mixed GWR (one local group plus global columns) and multiscale GWR (one bandwidth per local
group), fitted by back-fitting the local hat matrices of `local_hat.jl`.
"""

export select_mixed_bandwidth, select_multiscale_bandwidths
export mixed_gwr_predict, multiscale_gwr_predict

function _backfit_components(
    y::Vector{Float64}, local_hats::Vector{Matrix{Float64}},
    global_hat::Union{Nothing,Matrix{Float64}};
    tolerance::Float64=1e-5, max_iterations::Int=200,
)
    local_components = [zeros(Float64, length(y)) for _ in local_hats]
    global_component = _apply_global(global_hat, y)
    previous = Inf
    for iteration in 1:max_iterations
        for j in eachindex(local_hats)
            partial = y - global_component
            for k in eachindex(local_components)
                k == j || (partial .-= local_components[k])
            end
            local_components[j] = local_hats[j] * partial
        end
        local_sum = isempty(local_components) ? zeros(Float64, length(y)) : reduce(+, local_components)
        global_component = _apply_global(global_hat, y - local_sum)
        fitted = local_sum + global_component
        rss = sum(abs2, y - fitted)
        change = isfinite(previous) ? abs(previous - rss) / max(previous, eps()) : Inf
        change < tolerance && return (; local_components, global_component, fitted, converged=true, iterations=iteration)
        previous = rss
    end
    local_sum = isempty(local_components) ? zeros(Float64, length(y)) : reduce(+, local_components)
    fitted = local_sum + global_component
    return (; local_components, global_component, fitted, converged=false, iterations=max_iterations)
end

function _linear_loocv_rmse(hat::Matrix{Float64}, y::Vector{Float64})
    fitted = hat * y
    denominator = 1 .- diag(hat)
    valid = abs.(denominator) .> sqrt(eps())
    any(valid) || return Inf
    errors = (y[valid] - fitted[valid]) ./ denominator[valid]
    return sqrt(mean(abs2, errors))
end

function select_mixed_bandwidth(
    Xlocal::Matrix{Float64}, Xglobal::Matrix{Float64}, y::Vector{Float64},
    lonlat::Matrix{Float64}, candidates::Vector{Float64}, kernel::Function;
    adaptive::Bool=true, ridge::Float64=1e-8,
    tolerance::Float64=1e-5, max_iterations::Int=200,
)
    distances = haversine_distance_matrix(lonlat, lonlat)
    global_hat = global_projection(Xglobal; ridge)
    rows = NamedTuple[]
    best_bandwidth = 0.0
    best_rmse = Inf
    for bandwidth in candidates
        # A fixed-km bandwidth has no relationship to the station count, so the sanity check
        # against `length(y)` only applies to the adaptive neighbor-count case.
        (!adaptive || bandwidth < length(y)) || continue
        try
            local_hat = local_hat_matrix(Xlocal, Xlocal, distances, bandwidth, kernel; adaptive, ridge)
            identity_fit = _backfit_components(
                Matrix{Float64}(I, length(y), length(y)), [local_hat], global_hat;
                tolerance, max_iterations,
            )
            hat = identity_fit.fitted
            rmse = _linear_loocv_rmse(hat, y)
            push!(rows, (; bandwidth, RMSE=rmse, converged=identity_fit.converged,
                iterations=identity_fit.iterations, status="success", error=""))
            if identity_fit.converged && rmse < best_rmse
                best_bandwidth, best_rmse = bandwidth, rmse
            end
        catch error
            push!(rows, (; bandwidth, RMSE=Inf, converged=false, iterations=0,
                status="failed", error=sprint(showerror, error)))
        end
    end
    best_bandwidth > 0 || throw(ArgumentError("all Mixed GWR bandwidth candidates failed"))
    return best_bandwidth, DataFrame(rows)
end

function _backfit_components(
    Y::Matrix{Float64}, local_hats::Vector{Matrix{Float64}}, global_hat::Matrix{Float64};
    tolerance::Float64=1e-5, max_iterations::Int=200,
)
    local_components = [zeros(Float64, size(Y)) for _ in local_hats]
    global_component = global_hat * Y
    previous = Inf
    for iteration in 1:max_iterations
        for j in eachindex(local_hats)
            partial = Y - global_component
            for k in eachindex(local_components)
                k == j || (partial .-= local_components[k])
            end
            local_components[j] = local_hats[j] * partial
        end
        local_sum = isempty(local_components) ? zeros(Float64, size(Y)) : reduce(+, local_components)
        global_component = global_hat * (Y - local_sum)
        fitted = local_sum + global_component
        rss = sum(abs2, Y - fitted)
        change = isfinite(previous) ? abs(previous - rss) / max(previous, eps()) : Inf
        change < tolerance && return (; local_components, global_component, fitted, converged=true, iterations=iteration)
        previous = rss
    end
    local_sum = isempty(local_components) ? zeros(Float64, size(Y)) : reduce(+, local_components)
    fitted = local_sum + global_component
    return (; local_components, global_component, fitted, converged=false, iterations=max_iterations)
end

function select_multiscale_bandwidths(
    local_groups::Vector{Matrix{Float64}}, Xglobal::Matrix{Float64}, y::Vector{Float64},
    lonlat::Matrix{Float64}, candidates::Vector{Float64}, kernel::Function;
    adaptive::Bool=true, ridge::Float64=1e-8,
    tolerance::Float64=1e-5, max_iterations::Int=200,
)
    isempty(local_groups) && return Float64[], DataFrame(), true
    # A fixed-km bandwidth has no relationship to the station count, so the sanity check
    # against `length(y)` only applies to the adaptive neighbor-count case.
    usable = adaptive ? filter(<(length(y)), candidates) : candidates
    isempty(usable) && throw(ArgumentError("no usable multiscale bandwidth candidates"))
    distances = haversine_distance_matrix(lonlat, lonlat)
    global_hat = global_projection(Xglobal; ridge)
    bandwidths = fill(last(usable), length(local_groups))
    components = [zeros(Float64, length(y)) for _ in local_groups]
    global_component = global_hat * y
    rows = NamedTuple[]
    previous_rss = Inf
    for iteration in 1:max_iterations
        changed = false
        for group_index in eachindex(local_groups)
            partial = y - global_component
            for other in eachindex(components)
                other == group_index || (partial .-= components[other])
            end
            best_bw, best_rmse, best_hat = 0.0, Inf, nothing
            for bandwidth in usable
                try
                    hat = local_hat_matrix(
                        local_groups[group_index], local_groups[group_index], distances,
                        bandwidth, kernel; adaptive, ridge,
                    )
                    rmse = _linear_loocv_rmse(hat, partial)
                    push!(rows, (; iteration, group_index, bandwidth, RMSE=rmse,
                        status="success", error="", selected=false))
                    if rmse < best_rmse
                        best_bw, best_rmse, best_hat = bandwidth, rmse, hat
                    end
                catch error
                    push!(rows, (; iteration, group_index, bandwidth, RMSE=Inf,
                        status="failed", error=sprint(showerror, error), selected=false))
                end
            end
            best_hat === nothing && return bandwidths, DataFrame(rows), false
            selected_index = findlast(row -> row.iteration == iteration &&
                row.group_index == group_index && row.bandwidth == best_bw && row.status == "success", rows)
            rows[selected_index] = merge(rows[selected_index], (; selected=true))
            changed |= bandwidths[group_index] != best_bw
            bandwidths[group_index] = best_bw
            components[group_index] = best_hat * partial
        end
        local_sum = reduce(+, components)
        global_component = global_hat * (y - local_sum)
        fitted = local_sum + global_component
        rss = sum(abs2, y - fitted)
        relative_change = isfinite(previous_rss) ?
            abs(previous_rss - rss) / max(previous_rss, eps()) : Inf
        (!changed && relative_change < tolerance) && return bandwidths, DataFrame(rows), true
        previous_rss = rss
    end
    return bandwidths, DataFrame(rows), false
end

"""
`train_distances`/`target_distances` let a caller that fits many response columns over the same
station geometry hand in the haversine matrices instead of having them rebuilt here.

The contract is "the distances of the coordinate arguments exactly as passed":
`train_distances[i, j]` is the distance between `train_lonlat[i, :]` and `train_lonlat[j, :]`, and
`target_distances[i, j]` between `train_lonlat[i, :]` and `target_lonlat[j, :]`. The per-column
`indices`/`target_indices` subsetting below is then applied to the matrix the same way it is
applied to the coordinates, so the submatrix is elementwise what `haversine_distance_matrix`
would have returned for the subset. `nothing` (the default) computes them here.

`unsupported` is forwarded to the *target* hat only, never to the back-fit hat - see
[`local_hat_matrix`](@ref) for why that distinction matters. `:missing` is what lets a target the
local fit could not support surface as a `NaN` prediction instead of a silent zero correction.

`train_weights`/`target_weights` are the [`local_weight_matrix`](@ref) of `train_distances` and
`target_distances` at `bandwidth` (the train one without `exclude_self`, the target one with the
caller's `exclude_self`). They describe the full station set, so they are used only for a column
in which every training row is valid; any other column computes its own.
"""
function mixed_gwr_predict(
    Xlocal_train::Matrix{Float64}, Xglobal_train::Matrix{Float64}, Ytrain::Matrix{Float64},
    train_lonlat::Matrix{Float64}, Xlocal_target::Matrix{Float64},
    Xglobal_target::Matrix{Float64}, target_lonlat::Matrix{Float64}, bandwidth::Float64,
    kernel::Function; adaptive::Bool=true,
    ridge::Float64=1e-8, tolerance::Float64=1e-5, max_iterations::Int=200,
    exclude_self::Bool=false, unsupported::Symbol=:zero,
    train_distances::Union{Nothing,Matrix{Float64}}=nothing,
    target_distances::Union{Nothing,Matrix{Float64}}=nothing,
    train_weights::Union{Nothing,Matrix{Float64}}=nothing,
    target_weights::Union{Nothing,Matrix{Float64}}=nothing,
)
    if exclude_self
        size(train_lonlat) == size(target_lonlat) ||
            throw(DimensionMismatch("exclude_self requires matching train and target coordinates"))
        train_lonlat == target_lonlat ||
            throw(ArgumentError("exclude_self requires train and target coordinates in the same order"))
    end
    prediction = fill(NaN, size(Xlocal_target, 1), size(Ytrain, 2))
    converged = trues(size(Ytrain, 2))
    # Operators per set of valid stations; `Vector{Int}` keys hash by content.
    cache = Dict{Vector{Int},NamedTuple}()
    for time in axes(Ytrain, 2)
        valid = isfinite.(@view Ytrain[:, time])
        if count(valid) <= size(Xlocal_train, 2) + size(Xglobal_train, 2) + 2
            converged[time] = false
            continue
        end
        indices = findall(valid)
        operators = get!(cache, indices) do
            local_train_valid = Xlocal_train[indices, :]
            global_train_valid = Xglobal_train[indices, :]
            lonlat_valid = train_lonlat[indices, :]
            target_indices = exclude_self ? indices : collect(axes(Xlocal_target, 1))
            local_target_valid = Xlocal_target[target_indices, :]
            global_target_valid = Xglobal_target[target_indices, :]
            target_lonlat_valid = target_lonlat[target_indices, :]
            # Clamping to the valid station count only makes sense for an adaptive neighbor
            # count; a fixed-km bandwidth must be applied as given.
            adjusted_bandwidth = adaptive ? min(bandwidth, length(indices) - 1) : bandwidth
            train_distances_valid = train_distances === nothing ?
                haversine_distance_matrix(lonlat_valid, lonlat_valid) :
                distance_subset(train_distances, indices, indices)
            target_distances_valid = target_distances === nothing ?
                haversine_distance_matrix(lonlat_valid, target_lonlat_valid) :
                distance_subset(target_distances, indices, target_indices)
            all_rows = length(indices) == size(Ytrain, 1)
            local_hat = local_hat_matrix(
                local_train_valid, local_train_valid,
                train_distances_valid, adjusted_bandwidth, kernel;
                adaptive, ridge, weights=all_rows ? train_weights : nothing,
            )
            target_hat = local_hat_matrix(
                local_train_valid, local_target_valid,
                target_distances_valid, adjusted_bandwidth, kernel;
                adaptive, ridge, exclude_self, unsupported,
                weights=all_rows ? target_weights : nothing,
            )
            # `nothing` without global columns; see `_apply_global`.
            global_hat = isempty(global_train_valid) ? nothing :
                global_projection(global_train_valid; ridge)
            (; local_train_valid, global_train_valid, global_target_valid,
                target_indices, local_hat, target_hat, global_hat)
        end
        y = Vector{Float64}(Ytrain[indices, time])
        fitted = _backfit_components(
            y, [operators.local_hat], operators.global_hat; tolerance, max_iterations,
        )
        if !fitted.converged
            converged[time] = false
            continue
        end
        global_beta = isempty(operators.global_train_valid) ? Float64[] :
            (operators.global_train_valid' * operators.global_train_valid + ridge * I) \
            (operators.global_train_valid' * (y - fitted.local_components[1]))
        partial = y - (isempty(operators.global_train_valid) ? zeros(length(y)) :
            operators.global_train_valid * global_beta)
        target_prediction = operators.target_hat * partial
        isempty(operators.global_target_valid) ||
            (target_prediction .+= operators.global_target_valid * global_beta)
        prediction[operators.target_indices, time] = target_prediction
    end
    return prediction, converged
end

function multiscale_gwr_predict(
    local_train::Vector{Matrix{Float64}}, Xglobal_train::Matrix{Float64},
    Ytrain::Matrix{Float64}, train_lonlat::Matrix{Float64},
    local_target::Vector{Matrix{Float64}}, Xglobal_target::Matrix{Float64},
    target_lonlat::Matrix{Float64}, bandwidths::Vector{Float64}, kernel::Function;
    adaptive::Bool=true, ridge::Float64=1e-8,
    tolerance::Float64=1e-5, max_iterations::Int=200, exclude_self::Bool=false,
    unsupported::Symbol=:zero,
    train_distances::Union{Nothing,Matrix{Float64}}=nothing,
    target_distances::Union{Nothing,Matrix{Float64}}=nothing,
)
    if exclude_self
        isempty(Xglobal_train) ||
            throw(ArgumentError("multiscale exclude_self does not support global variables"))
        size(train_lonlat) == size(target_lonlat) ||
            throw(DimensionMismatch("exclude_self requires matching train and target coordinates"))
        train_lonlat == target_lonlat ||
            throw(ArgumentError("exclude_self requires train and target coordinates in the same order"))
    end
    prediction = fill(NaN, size(target_lonlat, 1), size(Ytrain, 2))
    converged = trues(size(Ytrain, 2))
    # Operators per set of valid stations; `Vector{Int}` keys hash by content.
    cache = Dict{Vector{Int},NamedTuple}()
    for time in axes(Ytrain, 2)
        valid = isfinite.(@view Ytrain[:, time])
        min_required = sum(size(X, 2) for X in local_train) + size(Xglobal_train, 2) + 2
        if count(valid) <= min_required
            converged[time] = false
            continue
        end
        indices = findall(valid)
        operators = get!(cache, indices) do
            local_valid = [X[indices, :] for X in local_train]
            global_valid = Xglobal_train[indices, :]
            lonlat_valid = train_lonlat[indices, :]
            target_indices = exclude_self ? indices : collect(axes(target_lonlat, 1))
            local_target_valid = [X[target_indices, :] for X in local_target]
            global_target_valid = Xglobal_target[target_indices, :]
            target_lonlat_valid = target_lonlat[target_indices, :]
            # Clamping to the valid station count only makes sense for an adaptive neighbor
            # count; a fixed-km bandwidth must be applied as given.
            adjusted_bandwidths = adaptive ? min.(bandwidths, length(indices) - 1) : bandwidths
            # Same `train_distances`/`target_distances` contract as `mixed_gwr_predict`; see its
            # docstring. Local names differ from the keywords so the `do` block does not capture.
            train_distances_valid = train_distances === nothing ?
                haversine_distance_matrix(lonlat_valid, lonlat_valid) :
                distance_subset(train_distances, indices, indices)
            train_hats = [local_hat_matrix(X, X, train_distances_valid, bw, kernel; adaptive, ridge)
                for (X, bw) in zip(local_valid, adjusted_bandwidths)]
            global_hat = global_projection(global_valid; ridge)
            target_hats = if exclude_self
                [local_hat_matrix(X, X, train_distances_valid, bw, kernel;
                    adaptive, ridge, exclude_self=true, unsupported)
                    for (X, bw) in zip(local_valid, adjusted_bandwidths)]
            else
                target_distances_valid = target_distances === nothing ?
                    haversine_distance_matrix(lonlat_valid, target_lonlat_valid) :
                    distance_subset(target_distances, indices, target_indices)
                [local_hat_matrix(X, Xt, target_distances_valid, bw, kernel;
                    adaptive, ridge, unsupported)
                    for (X, Xt, bw) in zip(local_valid, local_target_valid, adjusted_bandwidths)]
            end
            (; local_valid, global_valid, global_target_valid, target_indices,
                train_hats, target_hats, global_hat)
        end
        y = Vector{Float64}(Ytrain[indices, time])
        fitted = _backfit_components(
            y, operators.train_hats, operators.global_hat; tolerance, max_iterations,
        )
        if !fitted.converged
            converged[time] = false
            continue
        end
        local_sum = reduce(+, fitted.local_components)
        global_beta = isempty(operators.global_valid) ? Float64[] :
            (operators.global_valid' * operators.global_valid + ridge * I) \
            (operators.global_valid' * (y - local_sum))
        global_train_component = isempty(operators.global_valid) ? zeros(length(y)) :
            operators.global_valid * global_beta
        target_prediction = isempty(operators.global_target_valid) ?
            zeros(length(operators.target_indices)) : operators.global_target_valid * global_beta
        for group_index in eachindex(operators.local_valid)
            partial = y - global_train_component
            for other in eachindex(fitted.local_components)
                other == group_index || (partial .-= fitted.local_components[other])
            end
            target_prediction .+= operators.target_hats[group_index] * partial
        end
        prediction[operators.target_indices, time] = target_prediction
    end
    return prediction, converged
end

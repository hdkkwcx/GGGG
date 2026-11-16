function _local_predict(
    Xtrain::Matrix{Float64}, y::Vector{Float64}, Xtarget::Matrix{Float64},
    weights::Matrix{Float64}, ridge::Float64,
)
    prediction = fill(NaN, size(Xtarget, 1))
    p = size(Xtrain, 2)
    for target in axes(Xtarget, 1)
        w = @view weights[target, :]
        valid = (w .> 0) .& isfinite.(y) .& vec(all(isfinite, Xtrain; dims=2))
        count(valid) >= p + 1 || continue
        X = Xtrain[valid, :]; values = y[valid]; selected_weights = w[valid]
        beta = try
            (X' * (selected_weights .* X) + ridge * I) \
                (X' * (selected_weights .* values))
        catch
            continue
        end
        all(isfinite, @view Xtarget[target, :]) &&
            (prediction[target] = dot(@view(Xtarget[target, :]), beta))
    end
    return prediction
end

function _global_predict(
    Xtrain::Matrix{Float64}, y::Vector{Float64}, Xtarget::Matrix{Float64},
    ridge::Float64; leave_one_out::Bool=false,
)
    isempty(Xtrain) && return zeros(size(Xtarget, 1))
    if !leave_one_out
        beta = (Xtrain' * Xtrain + ridge * I) \ (Xtrain' * y)
        return Xtarget * beta
    end
    size(Xtrain, 1) == size(Xtarget, 1) ||
        throw(DimensionMismatch("global leave-one-out requires matching rows"))
    prediction = fill(NaN, size(Xtrain, 1))
    for target in axes(Xtrain, 1)
        keep = collect(axes(Xtrain, 1)) .!= target
        X = Xtrain[keep, :]; values = y[keep]
        beta = try (X' * X + ridge * I) \ (X' * values) catch; continue end
        prediction[target] = dot(@view(Xtrain[target, :]), beta)
    end
    return prediction
end

function _multiscale_predict_damped(
    local_train::Vector{Matrix{Float64}}, global_train::Matrix{Float64},
    response::Matrix{Float64}, train_lonlat::Matrix{Float64},
    local_target::Vector{Matrix{Float64}}, global_target::Matrix{Float64},
    target_lonlat::Matrix{Float64}, bandwidths::Vector{Float64}, kernel::GWKernel, cfg;
    adaptive::Bool=true, exclude_self::Bool=false, unsupported::Symbol=:zero,
    train_distances::Union{Nothing,Matrix{Float64}}=nothing,
    target_distances::Union{Nothing,Matrix{Float64}}=nothing,
    train_weights::Union{Nothing,Vector{Matrix{Float64}}}=nothing,
    target_weights::Union{Nothing,Vector{Matrix{Float64}}}=nothing,
)
    length(local_train) == length(local_target) == length(bandwidths) ||
        throw(DimensionMismatch(
            "local group count ($(length(local_train)) train, $(length(local_target)) target) " *
            "must match the bandwidth count ($(length(bandwidths)))",
        ))
    y = vec(response)
    # Supplied by `dynamic_covariate_predict`, which holds the station geometry fixed across every
    # hour it fits and so builds these once instead of once per hour. Same contract as
    # `mixed_gwr_predict`: the distances of `train_lonlat`/`target_lonlat` exactly as passed.
    train_d::Matrix{Float64} = train_distances === nothing ?
        haversine_distance_matrix(train_lonlat, train_lonlat) : train_distances
    # `train_weights`/`target_weights`, when given, hold one `local_weight_matrix` per group of
    # `train_d` (with this call's `exclude_self`) and `target_d` (without it) at that group's
    # bandwidth, so the hats below read the weights instead of recomputing them.
    train_hats = [local_hat_matrix(
        X, X, train_d, bandwidth, kernel; adaptive, ridge=cfg.ridge, exclude_self,
        weights=train_weights === nothing ? nothing : train_weights[group],
    ) for (group, (X, bandwidth)) in enumerate(zip(local_train, bandwidths))]
    # With no global columns the projection is an n x n zero matrix, and multiplying by it every
    # sweep was a whole wasted gemv per sweep. BLAS writes that product as +0.0 (it zeroes the
    # output for `beta = 0` and then adds signed-zero terms), so a component held at `zeros`
    # is the same value, and a non-finite iterate still leaves the reals through the local sum.
    has_global = !isempty(global_train)
    global_hat = has_global ?
        global_projection(global_train; ridge=cfg.ridge) : nothing
    local_components = [zeros(length(y)) for _ in local_train]
    global_component = has_global ? global_hat * y : zeros(length(y))
    previous = fill(NaN, length(y))
    converged = false
    relaxation = cfg.relaxation
    # Scratch for the sweep loop, allocated once instead of once per group per sweep. The sweep
    # count runs to the hundreds and this loop is the benchmark's innermost hot spot, so the
    # allocations were thousands of short-lived vectors per hour-fit, inside a threaded region
    # where the resulting GC pauses stop every thread at once.
    #
    # The arithmetic is deliberately left in exactly its original order and grouping, including
    # `component_sum` starting from zero (which is what `reduce(+, ...; init=zeros(...))` did) and
    # the `norm` of a materialised difference, so every value is bit-for-bit unchanged.
    n_obs = length(y)
    partial = Vector{Float64}(undef, n_obs)
    update = Vector{Float64}(undef, n_obs)
    component_sum = Vector{Float64}(undef, n_obs)
    global_update = Vector{Float64}(undef, n_obs)
    fitted = Vector{Float64}(undef, n_obs)
    difference = Vector{Float64}(undef, n_obs)
    for _ in 1:cfg.max_iterations
        for group in eachindex(local_train)
            partial .= y .- global_component
            for other in eachindex(local_components)
                other == group || (partial .-= local_components[other])
            end
            mul!(update, train_hats[group], partial)
            local_components[group] .= relaxation .* update .+
                (1 - relaxation) .* local_components[group]
        end
        fill!(component_sum, 0.0)
        for component in local_components
            component_sum .+= component
        end
        if has_global
            partial .= y .- component_sum
            mul!(global_update, global_hat, partial)
            global_component .= relaxation .* global_update .+
                (1 - relaxation) .* global_component
        end
        fitted .= component_sum .+ global_component
        # A diverging sweep makes `change` NaN, which fails the tolerance test silently and then
        # burns the whole iteration budget. Bail as soon as the iterate leaves the reals, the way
        # the undamped sibling does.
        all(isfinite, fitted) || return fill(NaN, size(target_lonlat, 1), 1), falses(1)
        change = if all(isfinite, previous)
            difference .= fitted .- previous
            norm(difference) / max(norm(previous), eps())
        else
            Inf
        end
        previous .= fitted
        if change < cfg.tolerance
            converged = true
            break
        end
    end
    converged || return fill(NaN, size(target_lonlat, 1), 1), falses(1)
    if exclude_self
        # `unsupported` cannot be honoured here. Unlike `mixed_gwr_predict` and
        # `multiscale_gwr_predict`, this branch has no separate target hat to mark: the
        # leave-one-out prediction *is* the fitted vector of the self-excluding `train_hats`
        # above, and a NaN row in those would propagate into `fitted`, trip the
        # `all(isfinite, fitted)` guard and discard the whole hour rather than one station.
        # Only `tuning_geometry = :loocv` and the no-inner-split joint scoring path reach it;
        # the default `:inner_spatial` never does.
        return reshape(previous, :, 1), trues(1)
    end
    local_sum = reduce(+, local_components; init=zeros(length(y)))
    global_beta = isempty(global_train) ? Float64[] :
        (global_train' * global_train + cfg.ridge * I) \
            (global_train' * (y - local_sum))
    prediction = isempty(global_target) ? zeros(size(target_lonlat, 1)) :
        global_target * global_beta
    target_d::Matrix{Float64} = target_distances === nothing ?
        haversine_distance_matrix(train_lonlat, target_lonlat) : target_distances
    for group in eachindex(local_train)
        partial = y - (isempty(global_train) ? zeros(length(y)) : global_train * global_beta)
        for other in eachindex(local_components)
            other == group || (partial .-= local_components[other])
        end
        target_hat = local_hat_matrix(
            local_train[group], local_target[group], target_d,
            bandwidths[group], kernel; adaptive, ridge=cfg.ridge, unsupported,
            weights=target_weights === nothing ? nothing : target_weights[group],
        )
        prediction .+= target_hat * partial
    end
    return reshape(prediction, :, 1), trues(1)
end

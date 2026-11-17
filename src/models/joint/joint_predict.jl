function dynamic_covariate_predict(
    context::JointFoldContext, residuals::Matrix{Float64}, method::String,
    bandwidths::Vector{Float64}, kernel::GWKernel;
    adaptive::Bool=true, time_indices::Vector{Int}=collect(axes(residuals, 2)),
    leave_one_out::Bool=false,
)
    method in ("residual_gwr", "mixed_gwr", "mgwr") ||
        throw(ArgumentError("unsupported joint model $method"))
    # Forwarded to every *target* hat below. It is a documented no-op on
    # `_multiscale_predict_damped`'s `exclude_self` path, which has no separate target hat;
    # see the comment at that early return.
    unsupported = context.config.unsupported_local_target
    target_count = leave_one_out ? size(residuals, 1) : size(context.target_lonlat, 1)
    prediction = fill(NaN, target_count, length(time_indices))
    # `Vector{Bool}`, not a `BitVector`. The loop below writes one entry per iteration from
    # many threads, and a BitVector packs 64 entries per word, so those writes would be
    # read-modify-writes of a shared word and would drop each other. That race predates the
    # scheduling change below but was mostly hidden by contiguous per-thread chunks; interleaved
    # scheduling would make it routine. One byte per entry makes each write independent.
    converged = fill(false, length(time_indices))
    # The station geometry is a property of the fold context, not of the hour: only which rows are
    # `valid` changes below. Building the two haversine matrices here rather than inside the loop
    # replaces ~2*N^2 transcendental evaluations per hour with one row/column gather. The
    # submatrix of a distance matrix is elementwise the distance matrix of the submatrix, so the
    # values handed to `local_hat_matrix` are unchanged.
    full_train_distances =
        haversine_distance_matrix(context.train_lonlat, context.train_lonlat)
    full_target_distances = leave_one_out ? full_train_distances :
        haversine_distance_matrix(context.train_lonlat, context.target_lonlat)
    # The local kernel weights depend on the distances, bandwidth, kernel and `exclude_self`
    # alone. On an hour where every training station is valid those are the same for every hour,
    # every MGWR group sharing a bandwidth, and every hat - yet the per-target k-th-neighbour
    # search was a third of each hat's cost. Built once here, keyed (distances, exclude_self,
    # bandwidth), and read-only inside the loop; hours that drop a station compute their own.
    ntrain = size(context.train_lonlat, 1)
    all_valid_bandwidths = unique(adaptive ? min.(bandwidths, ntrain - 1) : bandwidths)
    weight_kinds = leave_one_out ? ((:train, false), (:train, true)) :
        ((:train, false), (:target, false))
    full_weights = Dict(
        (kind, exclude, bw) => local_weight_matrix(
            kind === :train ? full_train_distances : full_target_distances, bw, kernel;
            adaptive, exclude_self=exclude,
        ) for (kind, exclude) in weight_kinds for bw in all_valid_bandwidths
    )
    # `:greedy`, not the default chunked schedule. Hour cost here spans orders of magnitude - a
    # dry or underdetermined hour falls straight through the `count(valid) > min_required` guard
    # below, while a wet one runs hundreds of back-fitting sweeps - and precipitation is strongly
    # autocorrelated in time, so splitting the range into contiguous per-thread chunks drops whole
    # storms onto single threads. Iterations write only their own column of `prediction` and
    # `converged` and read nothing another iteration writes, so the schedule cannot change a value.
    Threads.@threads :greedy for output_time in eachindex(time_indices)
        time = time_indices[output_time]
        train_design = _design_at(context, method, time; target=false)
        target_design = leave_one_out ? train_design :
            _design_at(context, method, time; target=true)
        y = Vector{Float64}(residuals[:, time])
        valid = isfinite.(y)
        for X in train_design.local_groups
            valid .&= vec(all(isfinite, X; dims=2))
        end
        isempty(train_design.global_design) ||
            (valid .&= vec(all(isfinite, train_design.global_design; dims=2)))
        min_required = sum(size(X, 2) for X in train_design.local_groups) +
            size(train_design.global_design, 2) + 2
        count(valid) > min_required || continue
        local_train = [Matrix{Float64}(X[valid, :]) for X in train_design.local_groups]
        global_train = Matrix{Float64}(train_design.global_design[valid, :])
        train_lonlat = context.train_lonlat[valid, :]
        # Subset the fold's distance matrices the same way the coordinates are subset, so the
        # models below receive exactly the matrices they would have rebuilt for themselves - and
        # skip the copy entirely on the common hour where nothing is dropped. See
        # `MixedGWR.distance_subset` for why these copies are worth avoiding.
        all_valid = all(valid)
        train_distances = all_valid ? full_train_distances : full_train_distances[valid, valid]
        if leave_one_out
            local_target = local_train
            global_target = global_train
            target_lonlat = train_lonlat
            target_distances = train_distances
        else
            local_target = target_design.local_groups
            global_target = target_design.global_design
            target_lonlat = context.target_lonlat
            target_distances = all_valid ? full_target_distances : full_target_distances[valid, :]
        end
        # Clamping to the valid station count only makes sense for an adaptive neighbor count; a
        # fixed-km bandwidth must be applied as given.
        adjusted = adaptive ? min.(bandwidths, count(valid) - 1) : bandwidths
        # `nothing` on an hour that dropped a station: the precomputed weights describe the
        # full station set, and the hats then compute their own as they always have.
        hat_weights(kind, exclude, bw) = all_valid ? full_weights[(kind, exclude, bw)] : nothing
        group_weights(kind, exclude) = all_valid ?
            [full_weights[(kind, exclude, bw)] for bw in adjusted] : nothing
        response = reshape(y[valid], :, 1)
        values_matrix, ok_vector = if leave_one_out && !isempty(global_train)
            # The held station is excluded analytically from the global coefficient;
            # local components are then tuned with a self-excluding spatial fit.
            global_loo = _global_predict(
                global_train, vec(response), global_train,
                context.config.ridge; leave_one_out=true,
            )
            partial = reshape(vec(response) - global_loo, :, 1)
            if method == "mgwr"
                local_prediction, ok = _multiscale_predict_damped(
                    local_train, zeros(Float64, count(valid), 0), partial,
                    train_lonlat, local_train, zeros(Float64, count(valid), 0),
                    train_lonlat, adjusted, kernel, context.config;
                    adaptive, exclude_self=true, unsupported,
                    train_distances, target_distances=train_distances,
                    train_weights=group_weights(:train, true),
                )
                local_prediction .+= global_loo
                local_prediction, ok
            else
                local_prediction, ok = mixed_gwr_predict(
                    only(local_train), zeros(Float64, count(valid), 0), partial,
                    train_lonlat, only(local_train), zeros(Float64, count(valid), 0),
                    train_lonlat, only(adjusted), kernel; adaptive, ridge=context.config.ridge,
                    tolerance=context.config.tolerance,
                    max_iterations=context.config.max_iterations,
                    exclude_self=true, unsupported,
                    train_distances, target_distances=train_distances,
                    train_weights=hat_weights(:train, false, only(adjusted)),
                    target_weights=hat_weights(:train, true, only(adjusted)),
                )
                local_prediction .+= global_loo
                local_prediction, ok
            end
        elseif method == "mgwr"
            _multiscale_predict_damped(
                local_train, global_train, response, train_lonlat,
                local_target, global_target, target_lonlat, adjusted, kernel, context.config;
                adaptive, exclude_self=leave_one_out, unsupported,
                train_distances, target_distances,
                train_weights=group_weights(:train, leave_one_out),
                # The self-excluding branch returns before it builds a target hat.
                target_weights=leave_one_out ? nothing : group_weights(:target, false),
            )
        else
            mixed_gwr_predict(
                only(local_train), global_train, response, train_lonlat,
                only(local_target), global_target, target_lonlat, only(adjusted), kernel;
                adaptive, ridge=context.config.ridge, tolerance=context.config.tolerance,
                max_iterations=context.config.max_iterations,
                exclude_self=leave_one_out, unsupported,
                train_distances, target_distances,
                train_weights=hat_weights(:train, false, only(adjusted)),
                target_weights=leave_one_out ? hat_weights(:train, true, only(adjusted)) :
                    hat_weights(:target, false, only(adjusted)),
            )
        end
        values = vec(values_matrix)
        ok = Bool(ok_vector[1])
        if leave_one_out
            prediction[findall(valid), output_time] = values
        else
            prediction[:, output_time] = values
        end
        converged[output_time] = ok
    end
    return prediction, converged
end

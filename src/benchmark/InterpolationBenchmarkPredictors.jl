function _gwr_predict(
    train_lonlat::Matrix{Float64}, values::Matrix{Float64}, target_lonlat::Matrix{Float64};
    kernel::Int, adaptive::Bool, bw::Float64, exclude_self::Bool=false,
)
    center = (mean(train_lonlat[:, 1]), mean(train_lonlat[:, 2]))
    X_train = build_X_intercept_centered(train_lonlat; center=center)
    X_target = build_X_intercept_centered(target_lonlat; center=center)
    distances = haversine_distance_matrix(train_lonlat, target_lonlat)
    if exclude_self
        size(train_lonlat, 1) == size(target_lonlat, 1) ||
            throw(DimensionMismatch("GWR exclude_self requires matching rows"))
        # An infinite distance leaves each station out of its own fit, at any bandwidth.
        distances = loocv_distances(distances)
    end
    weights = gw_weight(distances, bw; kernel=kernel, adaptive=adaptive)
    return st_gwr_predict_nanaware(X_train, values, weights; Xpred=X_target, min_obs=3)
end

"""
Hyperparameter-free reference predictors for one fold, estimated from training stations only.

- `zero`: constant 0. Hourly precipitation is ~93% dry, so this is a surprisingly strong
  RMSE competitor and the floor any method must clear.
- `train_clim`: the training stations' overall mean. A per-station climatology is not
  estimable for a held-out station, so the pooled mean is the station-free analogue.
- `hour_field_mean`: each hour's spatial mean over training stations - "is it raining
  anywhere in the domain right now", carrying no spatial structure at all. An interpolator
  that does not beat this is not doing spatial work.

Reported alongside the real methods so `metrics_*.csv` can be read as skill rather than as a
bare RMSE. See `NULL_METHODS`.
"""
function _null_fold_predictions(y_obs_train::Matrix{Float64}, n_val::Int)
    n_time = size(y_obs_train, 2)
    finite = filter(isfinite, vec(y_obs_train))
    climatology = isempty(finite) ? NaN : mean(finite)
    field_mean = Matrix{Float64}(undef, n_val, n_time)
    @inbounds for time in 1:n_time
        total = 0.0
        count = 0
        for station in axes(y_obs_train, 1)
            value = y_obs_train[station, time]
            if !isnan(value)
                total += value
                count += 1
            end
        end
        field_mean[:, time] .= count > 0 ? total / count : NaN
    end
    return Dict{String,Matrix{Float64}}(
        "zero" => zeros(Float64, n_val, n_time),
        "train_clim" => fill(climatology, n_val, n_time),
        "hour_field_mean" => field_mean,
    )
end

"""Build the provisional mixed-GWR design; replace this when variable roles are finalized."""
function build_mixed_gwr_designs(
    train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64},
)
    center = (mean(train_lonlat[:, 1]), mean(train_lonlat[:, 2]))
    local_train = build_X_intercept_centered(train_lonlat; center)
    local_target = build_X_intercept_centered(target_lonlat; center)
    global_train = zeros(Float64, size(train_lonlat, 1), 0)
    global_target = zeros(Float64, size(target_lonlat, 1), 0)
    return (; local_train, local_target, global_train, global_target)
end

function _mixed_gwr_predict(
    train_lonlat::Matrix{Float64}, values::Matrix{Float64},
    target_lonlat::Matrix{Float64}; bw::Float64, kernel::Int, adaptive::Bool=true,
    exclude_self::Bool=false,
)
    designs = build_mixed_gwr_designs(train_lonlat, target_lonlat)
    prediction, _ = mixed_gwr_predict(
        designs.local_train, designs.global_train, values, train_lonlat,
        designs.local_target, designs.global_target, target_lonlat, bw,
        kernel; adaptive, exclude_self,
    )
    return prediction
end

"""Build the provisional MGWR groups; replace these when covariates are finalized."""
function build_mgwr_designs(
    train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64},
)
    center = (mean(train_lonlat[:, 1]), mean(train_lonlat[:, 2]))
    X_train = build_X_intercept_centered(train_lonlat; center)
    X_target = build_X_intercept_centered(target_lonlat; center)
    local_train = [X_train[:, index:index] for index in axes(X_train, 2)]
    local_target = [X_target[:, index:index] for index in axes(X_target, 2)]
    global_train = zeros(Float64, size(train_lonlat, 1), 0)
    global_target = zeros(Float64, size(target_lonlat, 1), 0)
    group_names = ["intercept", "longitude", "latitude"]
    return (; local_train, local_target, global_train, global_target, group_names)
end

function _mgwr_predict(
    train_lonlat::Matrix{Float64}, values::Matrix{Float64},
    target_lonlat::Matrix{Float64}; bandwidths::Vector{Float64}, kernel::Int,
    adaptive::Bool=true, exclude_self::Bool=false,
)
    designs = build_mgwr_designs(train_lonlat, target_lonlat)
    prediction, _ = multiscale_gwr_predict(
        designs.local_train, designs.global_train, values, train_lonlat,
        designs.local_target, designs.global_target, target_lonlat, bandwidths,
        kernel; adaptive, exclude_self,
    )
    return prediction
end

function predict_selected(
    selected, method::String, mode::String,
    train_lonlat::Matrix{Float64}, target_lonlat::Matrix{Float64},
    y_obs_train::Matrix{Float64}, y_sat_train::Matrix{Float64}, y_sat_target::Matrix{Float64};
    dem_context=nothing, joint_context=nothing, hurdle_context=nothing,
    time_indices::Union{Nothing,Vector{Int}}=nothing,
)
    (mode, method) in BENCHMARK_RUNS ||
        throw(ArgumentError("unsupported benchmark method/mode pair: $method/$mode"))
    values = mode == "direct" ? y_obs_train : y_obs_train .- y_sat_train
    is_dem_model = dem_context !== nothing && hasproperty(selected, :dem_model) && selected.dem_model
    is_joint_model = joint_context !== nothing &&
        hasproperty(selected, :joint_model) && selected.joint_model
    # `time_indices` picks the hours (columns, in that order) a joint model predicts. A
    # `JointFoldContext` indexes hours itself, so the joint model takes the full matrices and the
    # hours separately; every other method is handed pre-sliced inputs instead.
    time_indices === nothing || (is_joint_model && mode == "residual") ||
        throw(ArgumentError("time_indices applies only to a joint residual model"))
    interpolated = if is_joint_model && mode == "residual"
        prediction, converged = dynamic_covariate_predict(
            joint_context, values, selected.joint_method,
            Float64.(selected.bandwidths), selected.kernel;
            adaptive=selected.adaptive,
            time_indices=something(time_indices, collect(axes(values, 2))),
        )
        any(.!converged) && @warn(
            "joint dynamic model did not converge for some hours",
            product=joint_context.product, method=selected.joint_method,
            failed_hours=count(.!converged),
        )
        # Folded in here rather than at the shared `y_sat_target .+ interpolated` tail below,
        # which serves every other method and must keep its historical behaviour.
        shrink = hasproperty(selected, :shrink) ? selected.shrink : 1.0
        isfinite(shrink) ? prediction .* shrink : prediction
    elseif is_dem_model && method == "gwr" && mode == "residual"
        designs = dem_context.all_local
        prediction, _ = mixed_gwr_predict(
            designs.mixed_local_train, zeros(Float64, size(train_lonlat, 1), 0), values,
            train_lonlat, designs.mixed_local_target,
            zeros(Float64, size(target_lonlat, 1), 0), target_lonlat,
            selected.bw, selected.kernel; adaptive=selected.adaptive,
            ridge=dem_context.dem.ridge,
            tolerance=dem_context.dem.tolerance,
            max_iterations=dem_context.dem.max_iterations,
        )
        prediction
    elseif is_dem_model && method == "mixed_gwr"
        designs = dem_context.mixed
        prediction, _ = mixed_gwr_predict(
            designs.mixed_local_train, designs.global_train, values, train_lonlat,
            designs.mixed_local_target, designs.global_target, target_lonlat,
            selected.bw, selected.kernel; adaptive=selected.adaptive,
            ridge=dem_context.dem.ridge,
            tolerance=dem_context.dem.tolerance,
            max_iterations=dem_context.dem.max_iterations,
        )
        prediction
    elseif is_dem_model && method == "mgwr"
        designs = dem_context.mixed
        prediction, _ = multiscale_gwr_predict(
            designs.multiscale_train, designs.global_train, values, train_lonlat,
            designs.multiscale_target, designs.global_target, target_lonlat,
            selected.bandwidths, selected.kernel;
            adaptive=selected.adaptive,
            ridge=dem_context.dem.ridge,
            tolerance=dem_context.dem.tolerance,
            max_iterations=dem_context.dem.max_iterations,
        )
        prediction
    elseif method in ("idw", "adw")
        selected_neighbors = ismissing(selected.neighbors) ? nothing :
            (selected.neighbors == 0 ? nothing : Int(selected.neighbors))
        predictor = method == "idw" ? idw_predict : adw_predict
        predictor(train_lonlat, values, target_lonlat;
            power=selected.power, neighbors=selected_neighbors)
    elseif method == "tps"
        tps_predict(train_lonlat, values, target_lonlat; smooth=selected.smooth)
    elseif method == "gwr"
        _gwr_predict(train_lonlat, values, target_lonlat;
            kernel=selected.kernel, adaptive=selected.adaptive, bw=selected.bw)
    elseif method == "hurdle_gwr"
        _predict_hurdle_selected(selected, train_lonlat, target_lonlat,
            y_obs_train, y_sat_train, y_sat_target, hurdle_context)
    elseif method == "mixed_gwr"
        _mixed_gwr_predict(
            train_lonlat, values, target_lonlat;
            bw=selected.bw, kernel=selected.kernel, adaptive=selected.adaptive,
        )
    elseif method == "mgwr"
        _mgwr_predict(
            train_lonlat, values, target_lonlat; bandwidths=Float64.(selected.bandwidths),
            kernel=selected.kernel, adaptive=selected.adaptive,
        )
    else
        throw(ArgumentError("unknown method: $method"))
    end
    time_indices === nothing || (y_sat_target = y_sat_target[:, time_indices])
    return mode == "direct" ? max.(interpolated, 0.0) : max.(y_sat_target .+ interpolated, 0.0)
end

"""
Refit a method's already-selected hyperparameters across the inner selection split and stitch the
held-out inner predictions into one `n_train × n_time` matrix.

This is what makes `auto` a fair comparison rather than a scan-row lookup. Each method's
`selected` scan row already carries an inner-split RMSE, but those numbers are scored on slightly
different cell sets: a method is only required to reach `min_tuning_coverage`, and in the full run
the selected candidates' coverage ranges 0.98-1.00, with `mgwr` — the method most likely to win —
sitting lowest. Comparing RMSEs computed over different denominators would quietly reward whichever
method dropped the hardest cells, on margins of about a percent. Re-predicting here lets the caller
score every contender on the intersection of their masks.

Deliberately no `dem_context`: the legacy DEM path is mutually exclusive with the joint path and is
not part of the reported comparison, so `auto` does not run there. `joint_contexts` is
`joint_selection_contexts` — one `JointFoldContext` per inner group, already built for the scan, and
matched here by `target_positions` rather than by position so a reordering cannot silently misalign
a context with the wrong stations.

`time_indices` restricts the result to the tuning hours, so `auto` compares methods over exactly
the hours the scan scored them on. The two model families need it applied at different points: a
`JointFoldContext` indexes hours internally, so the joint methods get the full matrices plus
`predict_selected`'s `time_indices` and predict only those hours, while everything else is given
pre-sliced inputs. Neither ever touches the other hours. `dynamic_covariate_predict` fits each
hour on its own, so an hour's prediction does not depend on which other hours are predicted with
it. On the 26,304-hour grid, predicting every hour and slicing cost the joint methods about twice
their whole scan, a third of each fold's time.
"""
function inner_selection_prediction(
    selected, method::String, mode::String, train_lonlat::Matrix{Float64},
    y_obs_train::Matrix{Float64}, y_sat_train::Matrix{Float64},
    selection_groups::Vector{Vector{Int}}; joint_contexts=nothing,
    time_indices::Union{Nothing,Vector{Int}}=nothing,
)
    is_joint = joint_contexts !== nothing &&
        hasproperty(selected, :joint_model) && selected.joint_model
    columns = time_indices === nothing ? collect(axes(y_obs_train, 2)) : time_indices
    obs = is_joint ? y_obs_train : y_obs_train[:, columns]
    sat = is_joint ? y_sat_train : y_sat_train[:, columns]
    n_train = size(obs, 1)
    out_of_fold = _selection_oof(selection_groups, n_train, length(columns), function (inner_train, group)
        joint_context = if joint_contexts === nothing
            nothing
        else
            entry = findfirst(candidate -> candidate.target_positions == group, joint_contexts)
            entry === nothing && throw(ArgumentError(
                "no joint selection context matches inner group $(group)",
            ))
            joint_contexts[entry].context
        end
        return predict_selected(
            selected, method, mode,
            train_lonlat[inner_train, :], train_lonlat[group, :],
            obs[inner_train, :], sat[inner_train, :], sat[group, :];
            joint_context, time_indices=is_joint ? time_indices : nothing,
        )
    end)
    return out_of_fold
end

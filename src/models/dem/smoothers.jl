"""Adaptive bisquare weights; a thin, behaviour-preserving wrapper over [`_gw_local_weights!`](@ref)."""
function _adaptive_bisquare!(
    weights::Vector{Float64}, distances::Vector{Float64}, neighbors::Int,
    buffer::Vector{Float64},
)
    _gw_local_weights!(weights, distances, Float64(neighbors), gw_weight_bisquare, buffer; adaptive=true)
end

function _adaptive_bisquare(distances::Vector{Float64}, neighbors::Int)
    n = length(distances)
    return _adaptive_bisquare!(
        Vector{Float64}(undef, n), distances, neighbors, Vector{Float64}(undef, n),
    )
end

"""
Per-target local weights for [`_gwr_smoothers`](@ref): the surviving neighbour indices and their
weights, for every target.

Split out because these depend on `distances`, `neighbors` and `exclude_self` alone — never on
the design matrix. `spatial_variability_test` rebuilds a smoother a thousand times over the same
geometry with only `X` permuted, and used to recompute this, including a sort per target, every
single time.

The `p + 1` sufficiency check lives here too, so it is raised from the caller's own task rather
than from inside a threaded permutation loop, where it would surface wrapped.
"""
function _gwr_local_weights(
    distances::Matrix{Float64}, neighbors::Int, p::Int; exclude_self::Bool=true,
)
    n = size(distances, 1)
    size(distances) == (n, n) || throw(DimensionMismatch("GWR smoother requires square distances"))
    entries = Vector{Tuple{Vector{Int},Vector{Float64}}}(undef, n)
    d = Vector{Float64}(undef, n)
    buffer = Vector{Float64}(undef, n)
    w = Vector{Float64}(undef, n)
    for target in 1:n
        copyto!(d, view(distances, :, target))
        exclude_self && (d[target] = Inf)
        _adaptive_bisquare!(w, d, neighbors, buffer)
        valid = w .> 0
        count(valid) >= p + 1 || throw(ArgumentError("insufficient local observations for bandwidth $neighbors"))
        indices = findall(valid)
        entries[target] = (indices, w[indices])
    end
    return entries
end

function _gwr_smoothers(
    X::Matrix{Float64}, distances::Matrix{Float64}, neighbors::Int;
    ridge::Float64=1e-8, exclude_self::Bool=true, local_weights=nothing,
)
    n, p = size(X)
    size(distances) == (n, n) || throw(DimensionMismatch("GWR smoother requires square distances"))
    entries = local_weights === nothing ?
        _gwr_local_weights(distances, neighbors, p; exclude_self) : local_weights
    smoothers = [zeros(Float64, n, n) for _ in 1:p]
    for target in 1:n
        indices, wv = entries[target]
        Xv = X[indices, :]
        A = Xv' * (wv .* Xv) + ridge * I
        B = A \ (Xv' .* wv')
        for coefficient in 1:p
            smoothers[coefficient][target, indices] = B[coefficient, :]
        end
    end
    return smoothers
end

function _select_gwr_bandwidth(
    X::Matrix{Float64}, y::Vector{Float64}, distances::Matrix{Float64}, candidates::Vector{Int};
    ridge::Float64=1e-8,
)
    rows = NamedTuple[]
    best = nothing
    for bandwidth in candidates
        bandwidth < size(X, 1) || continue
        try
            smoothers = _gwr_smoothers(X, distances, bandwidth; ridge, exclude_self=true)
            prediction = zeros(Float64, length(y))
            for j in axes(X, 2)
                prediction .+= X[:, j] .* (smoothers[j] * y)
            end
            rmse = sqrt(mean(abs2, prediction - y))
            push!(rows, (; bandwidth, RMSE=rmse, status="success", error=""))
            if best === nothing || rmse < best.RMSE
                best = (; bandwidth, RMSE=rmse, smoothers)
            end
        catch error
            push!(rows, (; bandwidth, RMSE=Inf, status="failed", error=sprint(showerror, error)))
        end
    end
    best === nothing && throw(ArgumentError("all GWR bandwidth candidates failed"))
    return best, DataFrame(rows)
end

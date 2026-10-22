"""
Local GWR operators shared by the mixed and multiscale estimators and the joint covariate
models: per-target kernel weights, the local hat matrix that maps a response onto its local
fit, and the global projection.
"""

export local_hat_matrix, local_weight_matrix, global_projection, distance_subset

"""
Rows `rows` and columns `columns` of a supplied distance matrix, without copying when the
selection is the whole matrix.

Nothing is dropped on most hours, and in the threaded hour loop of `dynamic_covariate_predict`
the large gather allocations, not the arithmetic, are what limit scaling. A view would still
allocate a wrapper and push a non-`Matrix` type through `local_hat_matrix`.
"""
function distance_subset(distances::Matrix{Float64}, rows, columns)
    length(rows) == size(distances, 1) && length(columns) == size(distances, 2) &&
        return distances
    return distances[rows, columns]
end

"""
`unsupported` decides how a target that cannot be fitted locally - fewer than `p + 1` stations
with positive weight - is recorded.

`:zero` (the default) leaves the all-zero row this matrix was allocated with, so the target
contributes exactly 0 to `H * y`, indistinguishable from a genuine "no local correction".

`:missing` writes `NaN` across the row, so `H * y` is `NaN` there and the benchmark's coverage
machinery (`min_tuning_coverage`, `_prediction_coverage`, `_common_method_mask`) can see it. Pass
it only for a hat used to *predict at targets*: a `NaN` row in a hat fed to a back-fit propagates
into `fitted`, trips its `all(isfinite, ...)` guard and discards the whole response column.

`weights`, when given, is [`local_weight_matrix`](@ref) of the same `distances`, `bw`, `kernel`,
`adaptive` and `exclude_self`: column `target` is read instead of recomputed. The weights never
depend on `X`, so a caller fitting many designs over one geometry builds them once.
"""
function local_hat_matrix(
    Xtrain::Matrix{Float64}, Xtarget::Matrix{Float64}, distances::Matrix{Float64},
    bw::Float64, kernel::GWKernel; adaptive::Bool=true, ridge::Float64=1e-8,
    exclude_self::Bool=false, unsupported::Symbol=:zero,
    weights::Union{Nothing,Matrix{Float64}}=nothing,
)
    unsupported in (:zero, :missing) ||
        throw(ArgumentError("unsupported must be :zero or :missing, got $unsupported"))
    ntrain, p = size(Xtrain)
    ntarget = size(Xtarget, 1)
    size(Xtarget, 2) == p || throw(DimensionMismatch("local design columns differ"))
    size(distances) == (ntrain, ntarget) || throw(DimensionMismatch("distance dimensions differ"))
    weights === nothing || size(weights) == size(distances) ||
        throw(DimensionMismatch("weight dimensions differ"))
    exclude_self && ntrain != ntarget &&
        throw(DimensionMismatch("exclude_self requires matching train and target rows"))
    H = zeros(Float64, ntarget, ntrain)
    # Allocated once per call, not per target. All buffers are call-local, so threaded callers
    # stay safe.
    d = Vector{Float64}(undef, ntrain)
    w = Vector{Float64}(undef, ntrain)
    buffer = Vector{Float64}(undef, ntrain)
    indices = Vector{Int}(undef, ntrain)
    wv = Vector{Float64}(undef, ntrain)
    hv = Vector{Float64}(undef, ntrain)
    Xv = Matrix{Float64}(undef, ntrain, p)
    WX = Matrix{Float64}(undef, ntrain, p)
    A = Matrix{Float64}(undef, p, p)
    xt = Vector{Float64}(undef, p)
    @inbounds for target in 1:ntarget
        if weights === nothing
            _target_local_weights!(w, d, buffer, distances, target, bw, kernel, adaptive, exclude_self)
        else
            copyto!(w, view(weights, :, target))
        end
        n_valid = 0
        for i in 1:ntrain
            w[i] > 0 || continue
            n_valid += 1
            indices[n_valid] = i
            wv[n_valid] = w[i]
            for j in 1:p
                Xv[n_valid, j] = Xtrain[i, j]
                WX[n_valid, j] = w[i] * Xtrain[i, j]
            end
        end
        if n_valid < p + 1
            unsupported === :missing && (H[target, :] .= NaN)
            continue
        end
        Xvalid = view(Xv, 1:n_valid, :)
        mul!(A, transpose(Xvalid), view(WX, 1:n_valid, :))
        for j in 1:p
            A[j, j] += ridge
        end
        for j in 1:p
            xt[j] = Xtarget[target, j]
        end
        # `A` is symmetric, so `xt' * (A \ (Xv' .* wv')) == (A \ xt)' * (Xv' .* wv')`: one
        # right-hand side instead of `n_valid`. A plain backslash rather than a Cholesky, so an
        # ill-conditioned or non-finite `A` behaves like a general solve.
        coefficients = A \ xt
        hvalid = view(hv, 1:n_valid)
        mul!(hvalid, Xvalid, coefficients)
        for t in 1:n_valid
            H[target, indices[t]] = hvalid[t] * wv[t]
        end
    end
    return H
end

"""One target's local weights for `local_hat_matrix`; `d` and `buffer` are scratch."""
function _target_local_weights!(
    w::Vector{Float64}, d::Vector{Float64}, buffer::Vector{Float64},
    distances::Matrix{Float64}, target::Int, bw::Float64, kernel::GWKernel,
    adaptive::Bool, exclude_self::Bool,
)
    copyto!(d, view(distances, :, target))
    exclude_self && (d[target] = Inf)
    return gw_weight!(w, d, bw; kernel, adaptive, buffer)
end

"""
Every target's local weights for [`local_hat_matrix`](@ref)'s `weights` keyword: column `target` is
exactly the `w` that `local_hat_matrix` would compute for it with the same arguments.
"""
function local_weight_matrix(
    distances::Matrix{Float64}, bw::Float64, kernel::GWKernel;
    adaptive::Bool=true, exclude_self::Bool=false,
)
    ntrain, ntarget = size(distances)
    exclude_self && ntrain != ntarget &&
        throw(DimensionMismatch("exclude_self requires matching train and target rows"))
    weights = Matrix{Float64}(undef, ntrain, ntarget)
    d = Vector{Float64}(undef, ntrain)
    w = Vector{Float64}(undef, ntrain)
    buffer = Vector{Float64}(undef, ntrain)
    for target in 1:ntarget
        _target_local_weights!(w, d, buffer, distances, target, bw, kernel, adaptive, exclude_self)
        copyto!(view(weights, :, target), w)
    end
    return weights
end

function global_projection(X::Matrix{Float64}; ridge::Float64=1e-8)
    isempty(X) && return zeros(Float64, size(X, 1), size(X, 1))
    return X * ((X' * X + ridge * I) \ X')
end

# `global_hat = nothing` stands for a model with no global columns, whose projection is an n x n
# zero matrix: the component is held at `zeros` instead of multiplied out (BLAS writes that
# product as +0.0, so the values are the same).
_apply_global(global_hat::Matrix{Float64}, v::Vector{Float64}) = global_hat * v
_apply_global(::Nothing, v::Vector{Float64}) = zeros(Float64, length(v))

"""
    GWR!(β, X, Y, wMat)

Local regression coefficients from precomputed weights: `β[i, :]` is the weighted least-squares
fit of `Y` on `X` with weights `wMat[:, i]`.

- `X`: [n_control, k_local]
- `Y`: [n_control]
- `wMat`: [n_control, n_target]
- `β`: [n_target, k_local]
"""
function GWR!(β::AbstractMatrix{T},
  X::AbstractMatrix{T}, Y::AbstractVector{T},
  wMat::AbstractMatrix{T})::Matrix{T} where {T<:Real}
  n_target = size(wMat, 2)
  @inbounds @threads for i in 1:n_target
    w = @view wMat[:, i]
    β[i, :] = solve_reg(X, Y, w) # [k_local, 1]
  end
  return β
end

function GWR(X::AbstractMatrix{T}, Y::AbstractVector{T},
  wMat::AbstractMatrix{T})::Matrix{T} where {T<:Real}
  β = zeros(T, size(wMat, 2), size(X, 2))
  GWR!(β, X, Y, wMat)
end

"""GWR from a distance matrix `dMat` [n_control, n_target], weighted by [`gw_weight`](@ref)."""
GWR(x::Matrix{T}, y::Vector{T}, dMat::AbstractMatrix{T}, bw::T;
  kernel::Int, adaptive::Bool=false) where {T<:Real} =
  GWR(x, y, gw_weight(dMat, bw; kernel, adaptive))

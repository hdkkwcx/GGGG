export solve_reg

"""
    solve_reg(X, Y, w)

Weighted least squares at one location for every column of `Y` at once,
`β = (X'WX) \\ X'Wy`.

- `X`: [n, p]
- `Y`: [n, ntime], or a vector
- `w`: [n]

Returns `β`, [p, ntime].
"""
function solve_reg(X::AbstractMatrix{T}, Y::AbstractMatrix{T}, w::AbstractVector{T}) where {T<:Real}
  w_sqrt = sqrt.(w)
  xw = X .* w_sqrt
  XtWX = xw' * xw
  p = size(X, 2)
  ntime = size(Y, 2)
  β = zeros(T, p, ntime)

  @inbounds @threads for i in 1:ntime
    yw = @view(Y[:, i]) .* w_sqrt
    β[:, i] .= XtWX \ (xw' * yw)
  end
  return β
end

solve_reg(X::AbstractMatrix{T}, y::AbstractVector{T}, w::AbstractVector{T}) where {T<:Real} =
  solve_reg(X, reshape(y, :, 1), w)

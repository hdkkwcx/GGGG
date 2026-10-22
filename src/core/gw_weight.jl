"""
Apply one kernel across a distance vector; a non-finite distance gets weight 0.

A function barrier: `kernel_function(::Integer)` is only known to return some `Function`, and
passing it in as an argument lets Julia specialise the loop on its concrete type.
"""
function _apply_kernel!(wv::AbstractVector{T}, vdist::AbstractVector{T}, bw::T,
  kerf::F) where {T<:Real,F}
  @inbounds for i in eachindex(wv, vdist)
    d = vdist[i]
    wv[i] = isfinite(d) ? kerf(d, bw) : zero(T)
  end
  return wv
end

"""
    gw_weight!(wv, vdist, bw; kernel=GAUSSIAN, adaptive=true, buffer=nothing)

GW weights of one target from its distances `vdist` to the control points: the only GW
weighting in the project.

`kernel` is a kernel index (`GAUSSIAN` ... `BOXCAR`) or a weight function `kernel(dist, bw)`.
A non-finite distance gets weight 0; that is how a point is left out, e.g. the target itself in
leave-one-out ([`loocv_distances`](@ref)).

With `adaptive`, `bw` is a neighbour count k and the bandwidth is the k-th smallest finite
distance. With fewer than k finite distances it is GWmodel's `k/m * maximum` over the m finite
ones. A zero bandwidth (k coincident points) gives weight 1 at distance 0 and 0 elsewhere.

`buffer`, when given, is scratch of at least `length(vdist)` for the order statistic. `vdist`
itself is never modified.
"""
function gw_weight!(wv::AbstractVector{T}, vdist::AbstractVector{T}, bw::T;
  kernel::GWKernel=GAUSSIAN, adaptive::Bool=true,
  buffer::Union{Nothing,AbstractVector{T}}=nothing) where {T<:Real}

  if adaptive
    scratch = buffer === nothing ? Vector{T}(undef, length(vdist)) : buffer
    m = 0
    @inbounds for d in vdist
      isfinite(d) || continue
      m += 1
      scratch[m] = d
    end
    m == 0 && return fill!(wv, zero(T))
    k = round(Int, bw)
    k >= 1 || throw(ArgumentError("an adaptive bandwidth needs at least one neighbour, got $bw"))
    finite = view(scratch, 1:m)
    bw = k <= m ? partialsort!(finite, k) : T(k) / m * maximum(finite)
    if iszero(bw)
      @inbounds for i in eachindex(wv, vdist)
        wv[i] = iszero(vdist[i]) ? one(T) : zero(T)
      end
      return wv
    end
  end

  return _apply_kernel!(wv, vdist, bw, kernel_function(kernel))
end

gw_weight(vdist::AbstractVector{T}, bw::T; kwargs...) where {T<:Real} =
  gw_weight!(Vector{T}(undef, length(vdist)), vdist, bw; kwargs...)

# Matrix form: one column of weights per target.
function gw_weight!(ws::AbstractMatrix{T}, dist::AbstractMatrix{T}, bw::T;
  kernel::GWKernel=GAUSSIAN, adaptive::Bool=true) where {T<:Real}
  n_control, n_target = size(dist)

  # One scratch vector for the whole matrix instead of one per column.
  buffer = adaptive ? Vector{T}(undef, n_control) : nothing
  @inbounds for j in 1:n_target
    gw_weight!(view(ws, :, j), view(dist, :, j), bw; kernel, adaptive, buffer)
  end
  return ws
end

function gw_weight(dist::AbstractMatrix{T}, bw::T; kernel::GWKernel=GAUSSIAN,
  adaptive::Bool=true) where {T<:Real}
  wMat = zeros(T, size(dist))
  gw_weight!(wMat, dist, bw; kernel, adaptive)
end


export gw_weight, gw_weight!

# ---------------------------------------------------------------------------------------------
# Where a dry hour sits

"""
Hours from each hour of one gauge series to its nearest wet hour (`>= wet`), looking back (`before`)
and ahead (`after`).

A distance is positive when a wet hour is reached, `Inf` when the series ends first, and negative
(minus the distance) when a missing hour is reached first - the rain could be hiding there, so the
true distance is only bounded. Wet hours are at distance 0; missing hours carry `NaN`. The matrix
method applies this to every row.
"""
function hours_to_rain(y::AbstractVector{<:Real}; wet::Real=WET_MM)
    n = length(y)
    before, after = fill(NaN, n), fill(NaN, n)
    for (range, out) in ((1:n, before), (n:-1:1, after))
        last, last_missing = 0, false
        for t in range
            v = y[t]
            if !isfinite(v)
                last, last_missing = t, true
            elseif v >= wet
                out[t] = 0.0
                last, last_missing = t, false
            else
                out[t] = last == 0 ? Inf : last_missing ? -abs(t - last) : abs(t - last)
            end
        end
    end
    return before, after
end

function hours_to_rain(Y::AbstractMatrix{<:Real}; wet::Real=WET_MM)
    before, after = fill(NaN, size(Y)), fill(NaN, size(Y))
    for i in axes(Y, 1)
        b, a = hours_to_rain(view(Y, i, :); wet)
        before[i, :] = b
        after[i, :] = a
    end
    return before, after
end

"""
Class of the distance from a dry hour to the station's nearest rain (`PROXIMITY_LEVELS`), 0 for wet or
missing hours and for hours whose class a nearby missing stretch leaves undecided.

With the nearest known wet hour at `dw` and the nearest missing hour at `dm`, the true distance lies in
`[min(dw, dm), dw]`; the class is kept only when both ends fall in the same class.
"""
function proximity_labels(before::AbstractMatrix, after::AbstractMatrix; edges=PROXIMITY_EDGES)
    size(before) == size(after) || throw(DimensionMismatch("before and after must have the same size"))
    labels = zeros(Int8, size(before))
    for k in eachindex(before)
        b, a = before[k], after[k]
        (isnan(b) || isnan(a) || b == 0 || a == 0) && continue
        dw = min(b > 0 ? b : Inf, a > 0 ? a : Inf)
        dm = min(b < 0 ? -b : Inf, a < 0 ? -a : Inf)
        known = searchsortedfirst(edges, dw)
        (dm >= dw || searchsortedfirst(edges, dm) == known) && (labels[k] = known)
    end
    return labels
end

"""
Which side of the station's rain a dry hour is on (`SIDE_LEVELS`), looking `window` hours each way: rain
on both sides (a gap inside an event), rain only ahead (before the onset), rain only behind (after the
end), or neither. 0 for wet or missing hours, and when a missing hour within the window leaves a side
undecided.
"""
function side_labels(before::AbstractMatrix, after::AbstractMatrix; window::Int=3)
    size(before) == size(after) || throw(DimensionMismatch("before and after must have the same size"))
    labels = zeros(Int8, size(before))
    side(v) = 0 < v <= window ? 1 : -window <= v < 0 ? -1 : 0   # rain, undecided, no rain in the window
    for k in eachindex(before)
        b, a = before[k], after[k]
        (isnan(b) || isnan(a) || b == 0 || a == 0) && continue
        sb, sa = side(b), side(a)
        (sb == -1 || sa == -1) && continue
        labels[k] = sb == 1 && sa == 1 ? 1 : sa == 1 ? 2 : sb == 1 ? 3 : 4
    end
    return labels
end

"""
Share of each station's reporting neighbours within `radius_km` that are wet at each hour; `NaN` where
a station has no neighbour in range or none of them reports.
"""
function neighbour_wet_share(Y::AbstractMatrix{<:Real}, lonlat::AbstractMatrix{<:Real};
        radius_km::Real=15.0, wet::Real=WET_MM)
    size(lonlat, 1) == size(Y, 1) || throw(DimensionMismatch("one lon/lat row per station"))
    distance = haversine_distance_matrix(lonlat, lonlat)
    neighbours = [[k for k in axes(Y, 1) if k != i && distance[k, i] <= radius_km] for i in axes(Y, 1)]
    share = fill(NaN, size(Y))
    for j in axes(Y, 2), i in axes(Y, 1)
        reporting, raining = 0, 0
        for k in neighbours[i]
            v = Y[k, j]
            isfinite(v) || continue
            reporting += 1
            raining += v >= wet
        end
        reporting > 0 && (share[i, j] = raining / reporting)
    end
    return share, neighbours
end

"""`NEIGHBOUR_LEVELS` of the gauge-dry cells from `neighbour_wet_share`; 0 elsewhere."""
function neighbour_labels(share::AbstractMatrix, Y_obs::AbstractMatrix; wet::Real=WET_MM)
    labels = zeros(Int8, size(share))
    for k in eachindex(share)
        s, o = share[k], Y_obs[k]
        (isfinite(o) && o < wet && isfinite(s)) || continue
        labels[k] = s == 0 ? 1 : s < 0.5 ? 2 : 3
    end
    return labels
end

"""
Share of reporting gauges that are wet, per hour; `NaN` for hours where fewer than `min_reporting` of
the stations report, which cannot be called network-dry.
"""
function network_wet_share(Y::AbstractMatrix{<:Real}; min_reporting::Real=0.9, wet::Real=WET_MM)
    share = fill(NaN, size(Y, 2))
    for j in axes(Y, 2)
        reporting = count(isfinite, view(Y, :, j))
        reporting >= min_reporting * size(Y, 1) || continue
        share[j] = count(v -> isfinite(v) && v >= wet, view(Y, :, j)) / reporting
    end
    return share
end

"""`NETWORK_LEVELS` of the gauge-dry cells from a per-hour `network_wet_share`; 0 elsewhere."""
function network_labels(share::AbstractVector, Y_obs::AbstractMatrix; wet::Real=WET_MM, edges=NETWORK_EDGES)
    length(share) == size(Y_obs, 2) || throw(DimensionMismatch("one share per hour"))
    labels = zeros(Int8, size(Y_obs))
    for j in axes(Y_obs, 2)
        s = share[j]
        isfinite(s) || continue
        level = s == 0 ? 1 : searchsortedfirst(edges, s) + 1
        for i in axes(Y_obs, 1)
            o = Y_obs[i, j]
            isfinite(o) && o < wet && (labels[i, j] = level)
        end
    end
    return labels
end

"""
How widely an estimate rains on the hours every reporting gauge is dry (`network_wet_share == 0`).

For each such hour, the share of the estimate's scored cells (inside `mask`) that are `>= wet`; returns
`n_hours` and, for each cutoff `c`, the share of those hours on which that share exceeds `c` (`c = 0`:
at least one gauge pixel rains). With no gauge wet anywhere, rain at a gauge pixel is either rain that
fell between the gauges or rain that never reached the ground.
"""
function phantom_extent(network_share::AbstractVector, Y_est::AbstractMatrix, mask::AbstractMatrix;
        wet::Real=WET_MM, cutoffs=[0.0, 0.05, 0.25, 0.5], key::NamedTuple=(;))
    length(network_share) == size(Y_est, 2) || throw(DimensionMismatch("one share per hour"))
    size(Y_est) == size(mask) || throw(DimensionMismatch("estimate and mask must match"))
    shares = Float64[]
    for j in axes(Y_est, 2)
        network_share[j] == 0 || continue
        scored = [Y_est[i, j] for i in axes(Y_est, 1) if mask[i, j] && isfinite(Y_est[i, j])]
        isempty(scored) || push!(shares, count(>=(wet), scored) / length(scored))
    end
    exceed = (; (Symbol("share_hours_gt_", round(Int, 100c)) => _ratio(count(>(c), shares), length(shares))
        for c in cutoffs)...)
    return merge(key, (; n_hours=length(shares), mean_wet_share=isempty(shares) ? NaN : mean(shares)), exceed)
end

"""Label `k` where `edges[k] <= x < edges[k + 1]`; 0 for `NaN` or out of range."""
function bin_labels(X::AbstractArray{<:Real}, edges::AbstractVector{<:Real})
    issorted(edges) || throw(ArgumentError("edges must be sorted"))
    labels = zeros(Int8, size(X))
    for k in eachindex(X)
        x = X[k]
        isnan(x) && continue
        b = searchsortedlast(edges, x)
        1 <= b < length(edges) && (labels[k] = b)
    end
    return labels
end

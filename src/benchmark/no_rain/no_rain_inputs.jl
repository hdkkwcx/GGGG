# ---------------------------------------------------------------------------------------------
# Inputs and sensitivity

"""
A long table (`time`, `station`, `value`) as a `[station, time]` matrix on `ref_ids` x `ref_times`.

`time + Hour(offset_hours)` is matched to `ref_times`: 9 turns a UTC hour-start label into the
hour-ending Beijing-time label the gauges use. Rows for other stations or hours are skipped; a
duplicate cell is an error. Returns `(; Y, matched)`.
"""
function long_to_wide_hourly(
    time::AbstractVector{DateTime}, station::AbstractVector{<:AbstractString}, value::AbstractVector{<:Real},
    ref_times::AbstractVector{DateTime}, ref_ids::AbstractVector{<:AbstractString}; offset_hours::Int=9,
)
    length(time) == length(station) == length(value) || throw(DimensionMismatch("time, station and value must match"))
    row_of = Dict(String(id) => i for (i, id) in enumerate(ref_ids))
    column_of = Dict(t => j for (j, t) in enumerate(ref_times))
    Y = fill(NaN, length(ref_ids), length(ref_times))
    seen = falses(size(Y))
    matched = 0
    for k in eachindex(time)
        i = get(row_of, String(station[k]), 0)
        j = get(column_of, time[k] + Hour(offset_hours), 0)
        (i == 0 || j == 0) && continue
        seen[i, j] && throw(ArgumentError("duplicate value for $(station[k]) at $(time[k])"))
        seen[i, j] = true
        Y[i, j] = value[k]
        matched += 1
    end
    return (; Y, matched)
end

"""
Pooled Pearson r between `A[:, j]` and `B[:, j + lag]` for each lag, over cells both finite and inside
`mask` at their own column. A clock offset between two series of the same field shows as a peak away
from lag 0.
"""
function lagged_correlation(A::AbstractMatrix, B::AbstractMatrix, mask::AbstractMatrix; lags=-3:3)
    size(A) == size(B) == size(mask) || throw(DimensionMismatch("A, B and mask must match"))
    rows = NamedTuple[]
    for lag in lags
        sx = sy = sxx = syy = sxy = 0.0
        n = 0
        for j in axes(A, 2)
            jj = j + lag
            1 <= jj <= size(B, 2) || continue
            for i in axes(A, 1)
                (mask[i, j] && mask[i, jj]) || continue
                x, y = A[i, j], B[i, jj]
                (isfinite(x) && isfinite(y)) || continue
                n += 1
                sx += x; sy += y; sxx += x^2; syy += y^2; sxy += x * y
            end
        end
        cov = sxy / n - (sx / n) * (sy / n)
        push!(rows, (; lag, n, r=cov / sqrt((sxx / n - (sx / n)^2) * (syy / n - (sy / n)^2))))
    end
    return DataFrame(rows)
end

"""
What setting every value below `tau` to zero would change, for each `tau` in `taus` (0 is the estimate
as it is): overall and dry-hour RMSE, the mean on gauge-dry cells, the volume ratio to the gauges, and
POD / POFD / FAR / CSI at `wet`. Descriptive only - a `tau` read off this table is chosen on held-out
cells, so it is not an out-of-sample result.
"""
function zero_threshold_sensitivity(Y_obs::AbstractMatrix, Y_est::AbstractMatrix, mask::AbstractMatrix;
        taus=[0.0, 0.05, 0.1, 0.2, 0.3, 0.5], wet::Real=WET_MM, key::NamedTuple=(;))
    size(Y_obs) == size(Y_est) == size(mask) || throw(DimensionMismatch("obs, estimate and mask must match"))
    rows = NamedTuple[]
    for tau in taus
        n = n_dry = hits = misses = false_alarms = 0
        sse = sse_dry = vol = vol_dry = vol_obs = 0.0
        for k in eachindex(Y_obs)
            mask[k] || continue
            o, e = Y_obs[k], Y_est[k]
            (isfinite(o) && isfinite(e)) || continue
            e = e < tau ? 0.0 : e
            n += 1
            sse += (e - o)^2
            vol += e
            vol_obs += o
            if o >= wet
                e >= wet ? (hits += 1) : (misses += 1)
            else
                n_dry += 1
                sse_dry += (e - o)^2
                vol_dry += e
                e >= wet && (false_alarms += 1)
            end
        end
        push!(rows, merge(key, (; tau=Float64(tau), n, RMSE=sqrt(_ratio(sse, n)), RMSE_dry=sqrt(_ratio(sse_dry, n_dry)),
            mean_dry=_ratio(vol_dry, n_dry), volume_ratio=_ratio(vol, vol_obs), POD=_ratio(hits, hits + misses),
            POFD=_ratio(false_alarms, n_dry), FAR=_ratio(false_alarms, hits + false_alarms),
            CSI=_ratio(hits, hits + misses + false_alarms))))
    end
    return DataFrame(rows)
end

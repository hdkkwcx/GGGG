# ---------------------------------------------------------------------------------------------
# Occurrence and dry-hour amounts

# Per-day accumulator columns; threshold t adds hits at BASE + 2t - 1 and false alarms at BASE + 2t.
const N_DRY, N_WET, VOL_DRY, VOL_ALL, ZERO_DRY, SSE_DRY, SSE_ALL, VOL_OBS, VOL_OBS_DRY = 1:9
const BASE = 9

"""Scores from one accumulator total; wet-side scores are `NaN` where the level has no gauge-wet cell."""
function _occurrence_scores(t::AbstractVector, n_thresholds::Int)
    n_dry, n_wet = t[N_DRY], t[N_WET]
    out = [_ratio(t[VOL_DRY], n_dry), sqrt(_ratio(t[SSE_DRY], n_dry)), _ratio(t[SSE_DRY], t[SSE_ALL])]
    for k in 1:n_thresholds
        a, b = t[BASE + 2k - 1], t[BASE + 2k]
        c, d = n_wet - a, n_dry - b
        pofd = _ratio(b, n_dry)
        if n_wet > 0
            hss_den = (a + c) * (c + d) + (a + b) * (b + d)
            append!(out, [pofd, _ratio(b, a + b), _ratio(a + b, n_wet), _ratio(d, c + d), _ratio(a, n_wet) - pofd,
                hss_den > 0 ? 2 * (a * d - b * c) / hss_den : NaN])
        else
            append!(out, [pofd, NaN, NaN, NaN, NaN, NaN])
        end
    end
    return out
end

const _THRESHOLD_SCORES = (:POFD, :FAR, :freq_bias, :NPV, :PSS, :HSS)

"""Rows of an occurrence table from per-level, per-day accumulators `acc[level, day, column]`."""
function _occurrence_rows(acc::Array{Float64,3}, levels, thresholds; reps::Int, seed::Int, key::NamedTuple)
    T = length(thresholds)
    rows = NamedTuple[]
    for (l, level) in enumerate(levels)
        per_day = acc[l, :, :]
        used = [d for d in axes(per_day, 1) if per_day[d, N_DRY] + per_day[d, N_WET] > 0]
        total = vec(sum(per_day; dims=1))
        n_dry, n_wet = total[N_DRY], total[N_WET]
        n_cells = n_dry + n_wet
        point = _occurrence_scores(total, T)
        ci = reps > 0 && !isempty(used) ?
            day_block_bootstrap(collect(eachindex(used)), per_day[used, :], t -> _occurrence_scores(t, T);
                reps, rng=MersenneTwister(seed + l)) : nothing
        lo(k) = ci === nothing ? NaN : ci.lo[k]
        hi(k) = ci === nothing ? NaN : ci.hi[k]
        base = (;
            level, n_cells=Int(n_cells), n_dry=Int(n_dry), n_wet=Int(n_wet), dry_share=_ratio(n_dry, n_cells),
            n_days=length(used),
            mean_dry=point[1], mean_dry_lo=lo(1), mean_dry_hi=hi(1),
            RMSE_dry=point[2], RMSE_dry_lo=lo(2), RMSE_dry_hi=hi(2),
            mean_obs_dry=_ratio(total[VOL_OBS_DRY], n_dry),
            zero_share_dry=_ratio(total[ZERO_DRY], n_dry), dry_volume_share=_ratio(total[VOL_DRY], total[VOL_ALL]),
            dry_sse_share=point[3],
            spurious_mm_per_year=_ratio(total[VOL_DRY], n_cells) * HOURS_PER_YEAR,
            estimate_mm_per_year=_ratio(total[VOL_ALL], n_cells) * HOURS_PER_YEAR,
            gauge_mm_per_year=_ratio(total[VOL_OBS], n_cells) * HOURS_PER_YEAR,
        )
        for (k, threshold) in enumerate(thresholds)
            offset = 3 + 6 * (k - 1)
            scores = Pair{Symbol,Any}[:threshold => Float64(threshold), :POD => _ratio(total[BASE + 2k - 1], n_wet),
                :specificity => 1 - point[offset + 1]]
            for (s, name) in enumerate(_THRESHOLD_SCORES)
                append!(scores, [name => point[offset + s], Symbol(name, "_lo") => lo(offset + s),
                    Symbol(name, "_hi") => hi(offset + s)])
            end
            push!(rows, merge(key, base, NamedTuple{Tuple(first.(scores))}(Tuple(last.(scores)))))
        end
    end
    return DataFrame(rows)
end

"""
Occurrence and dry-hour skill of an estimate against the gauges, per level of a stratifier.

`labels` (a matrix, or a function `(i, j) -> label`) puts cell `(i, j)` in level 1..length(levels); 0,
or a cell outside `mask` or with a non-finite value, leaves it out. A cell is gauge-wet at
`obs >= wet`, otherwise gauge-dry; the estimate is wet at `est >= threshold` for each of `thresholds`.
One row per level and threshold, prefixed by the fields of `key`:

- counts: `n_cells`, `n_dry`, `n_wet`, `dry_share`, `n_days`
- occurrence at the threshold: `POFD` = P(est wet | gauge dry), the false-alarm *rate*;
  `specificity` = 1 - POFD; `FAR` = P(gauge dry | est wet), the false-alarm *ratio*;
  `NPV` = P(gauge dry | est dry); `POD`; `freq_bias` = est-wet over gauge-wet cells; `PSS` = POD - POFD;
  `HSS`. Wet-side scores are `NaN` for a level without gauge-wet cells.
- dry-hour amounts, the same on every threshold row: `mean_dry` (what the estimate reports on gauge-dry
  cells, mm/h), `mean_obs_dry` (the gauge's own mean there: a few readings are non-zero but below 0.1),
  `RMSE_dry`, `zero_share_dry` (estimate exactly 0), `dry_volume_share` (estimate volume on gauge-dry
  cells over its total), `dry_sse_share` (squared error on gauge-dry cells over the level's total), and
  per year of record `spurious_mm_per_year` (estimate volume on gauge-dry cells), `estimate_mm_per_year`
  and `gauge_mm_per_year`.

`_lo` / `_hi` are 95% day-block bootstrap bounds; `reps = 0` leaves them `NaN`.
"""
function occurrence_table(
    labels, levels::AbstractVector, Y_obs::AbstractMatrix, Y_est::AbstractMatrix, mask::AbstractMatrix,
    day::AbstractVector{<:Integer}; thresholds=DRY_THRESHOLDS, wet::Real=WET_MM, reps::Int=1000,
    seed::Int=20260919, key::NamedTuple=(;),
)
    size(Y_obs) == size(Y_est) == size(mask) || throw(DimensionMismatch("obs, estimate and mask must match"))
    length(day) == size(Y_obs, 2) || throw(DimensionMismatch("one day per column"))
    labels isa AbstractMatrix && size(labels) != size(Y_obs) && throw(DimensionMismatch("labels must match obs"))
    T = length(thresholds)
    acc = zeros(length(levels), maximum(day), BASE + 2T)
    for j in axes(Y_obs, 2)
        d = day[j]
        for i in axes(Y_obs, 1)
            mask[i, j] || continue
            l = _label(labels, i, j)
            l == 0 && continue
            o, e = Y_obs[i, j], Y_est[i, j]
            (isfinite(o) && isfinite(e)) || continue
            err2 = (e - o)^2
            acc[l, d, VOL_ALL] += e
            acc[l, d, SSE_ALL] += err2
            acc[l, d, VOL_OBS] += o
            if o >= wet
                acc[l, d, N_WET] += 1
                for k in 1:T
                    e >= thresholds[k] && (acc[l, d, BASE + 2k - 1] += 1)
                end
            else
                acc[l, d, N_DRY] += 1
                acc[l, d, VOL_DRY] += e
                acc[l, d, VOL_OBS_DRY] += o
                acc[l, d, ZERO_DRY] += e == 0
                acc[l, d, SSE_DRY] += err2
                for k in 1:T
                    e >= thresholds[k] && (acc[l, d, BASE + 2k] += 1)
                end
            end
        end
    end
    return _occurrence_rows(acc, levels, thresholds; reps, seed, key)
end

"""
Squared error of `Y_method` against `Y_ref` on the same cells, per level of a stratifier.

Per level: `n`, `RMSE_method`, `RMSE_ref`, `improvement` = 1 - RMSE_method / RMSE_ref (positive: the
method is better), the mean of each (`mean_method`, `mean_ref`, `mean_obs`), and `sse_gap` = SSE_method
- SSE_ref with `sse_gap_share`, its share of the summed gap over all levels - meaningful when the levels
partition the cells, as the gauge x anchor quadrants do. `_lo` / `_hi`: paired day-block bootstrap
bounds of `improvement` and of `mean_diff` = mean_method - mean_ref.
"""
function paired_error_table(
    labels, levels::AbstractVector, Y_obs::AbstractMatrix, Y_method::AbstractMatrix, Y_ref::AbstractMatrix,
    mask::AbstractMatrix, day::AbstractVector{<:Integer}; reps::Int=1000, seed::Int=20260919,
    key::NamedTuple=(;),
)
    size(Y_obs) == size(Y_method) == size(Y_ref) == size(mask) || throw(DimensionMismatch("matrices must match"))
    acc = zeros(length(levels), maximum(day), 6)     # n, SSE method, SSE ref, sum method, sum ref, sum obs
    for j in axes(Y_obs, 2)
        d = day[j]
        for i in axes(Y_obs, 1)
            mask[i, j] || continue
            l = _label(labels, i, j)
            l == 0 && continue
            o, m, r = Y_obs[i, j], Y_method[i, j], Y_ref[i, j]
            (isfinite(o) && isfinite(m) && isfinite(r)) || continue
            acc[l, d, 1] += 1
            acc[l, d, 2] += (m - o)^2
            acc[l, d, 3] += (r - o)^2
            acc[l, d, 4] += m
            acc[l, d, 5] += r
            acc[l, d, 6] += o
        end
    end
    scores(t) = (1 - sqrt(_ratio(t[2], t[1])) / sqrt(_ratio(t[3], t[1])), _ratio(t[4] - t[5], t[1]))
    gaps = [sum(acc[l, :, 2]) - sum(acc[l, :, 3]) for l in eachindex(levels)]
    rows = NamedTuple[]
    for (l, level) in enumerate(levels)
        per_day = acc[l, :, :]
        used = [d for d in axes(per_day, 1) if per_day[d, 1] > 0]
        t = vec(sum(per_day; dims=1))
        point = scores(t)
        ci = reps > 0 && !isempty(used) ?
            day_block_bootstrap(collect(eachindex(used)), per_day[used, :], scores; reps, rng=MersenneTwister(seed + l)) :
            nothing
        push!(rows, merge(key, (;
            level, n=Int(t[1]), RMSE_method=sqrt(_ratio(t[2], t[1])), RMSE_ref=sqrt(_ratio(t[3], t[1])),
            improvement=point[1], improvement_lo=ci === nothing ? NaN : ci.lo[1],
            improvement_hi=ci === nothing ? NaN : ci.hi[1],
            mean_method=_ratio(t[4], t[1]), mean_ref=_ratio(t[5], t[1]), mean_obs=_ratio(t[6], t[1]),
            mean_diff=point[2], mean_diff_lo=ci === nothing ? NaN : ci.lo[2], mean_diff_hi=ci === nothing ? NaN : ci.hi[2],
            sse_gap=gaps[l], sse_gap_share=sum(gaps) != 0 ? gaps[l] / sum(gaps) : NaN,
        )))
    end
    return DataFrame(rows)
end

"""
The false-alarm rate of a *gauge* as an estimate for its neighbour, by separation.

Every ordered pair of stations `(i, k)` whose separation falls in `[edges_km[b], edges_km[b + 1])` scores
gauge `k` as the estimate at station `i`, on the hours both report inside `mask`; a last level pairs every
station with its nearest neighbour whatever the distance. Scored with `occurrence_table`'s columns at
`threshold = wet`, plus `n_pairs` and `median_km`. Two gauges a few km apart sample the same rain field
at two points, so this is the floor a point-to-area comparison cannot get under: part of what a
satellite pixel is charged as a false alarm is rain the gauge simply was not under.
"""
function neighbour_gauge_baseline(
    Y_obs::AbstractMatrix, lonlat::AbstractMatrix, mask::AbstractMatrix, day::AbstractVector{<:Integer};
    edges_km=[0.0, 5.0, 10.0, 15.0, 20.0, 30.0], wet::Real=WET_MM, reps::Int=1000, seed::Int=20260919,
    key::NamedTuple=(;),
)
    n = size(Y_obs, 1)
    distance = haversine_distance_matrix(lonlat, lonlat)
    for i in 1:n
        distance[i, i] = Inf
    end
    B = length(edges_km) - 1
    levels = vcat(["$(Int(edges_km[b]))-$(Int(edges_km[b + 1])) km" for b in 1:B], ["nearest gauge"])
    pairs = [Tuple{Int,Int}[] for _ in 1:(B + 1)]
    for i in 1:n, k in 1:n
        b = searchsortedlast(edges_km, distance[k, i])
        1 <= b <= B && push!(pairs[b], (i, k))
    end
    for i in 1:n
        push!(pairs[B + 1], (i, argmin(view(distance, :, i))))
    end
    acc = zeros(length(levels), maximum(day), BASE + 2)
    for (l, level_pairs) in enumerate(pairs), (i, k) in level_pairs, j in axes(Y_obs, 2)
        (mask[i, j] && mask[k, j]) || continue
        o, e = Y_obs[i, j], Y_obs[k, j]
        (isfinite(o) && isfinite(e)) || continue
        d = day[j]
        acc[l, d, VOL_ALL] += e
        acc[l, d, SSE_ALL] += (e - o)^2
        acc[l, d, VOL_OBS] += o
        if o >= wet
            acc[l, d, N_WET] += 1
            e >= wet && (acc[l, d, BASE + 1] += 1)
        else
            acc[l, d, N_DRY] += 1
            acc[l, d, VOL_DRY] += e
            acc[l, d, VOL_OBS_DRY] += o
            acc[l, d, ZERO_DRY] += e == 0
            acc[l, d, SSE_DRY] += (e - o)^2
            e >= wet && (acc[l, d, BASE + 2] += 1)
        end
    end
    table = _occurrence_rows(acc, levels, [wet]; reps, seed, key)
    table.n_pairs = length.(pairs)
    table.median_km = [isempty(p) ? NaN : median(distance[k, i] for (i, k) in p) for p in pairs]
    return table
end

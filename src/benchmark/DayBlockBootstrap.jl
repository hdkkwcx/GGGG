"""
Pooled scores and paired day-block bootstrap intervals from per-group sums, for the scripts that
score a benchmark run's out-of-fold (OOF) files without refitting anything.

- `read_wide_values` / `read_wide_mask`: an OOF table (`time`, then one column per station) as a
  stations × hours matrix, checked against the grid's hours.
- `group_sums`: per group, the sums every score is built from, in one pass over the cells.
- `sum_stats`: n, RMSE, bias, |bias|, r, and POD, FAR, CSI, FBI and |FBI − 1| per threshold, from
  one row of those sums. Events use `>=` on both sides, as `metric_event` does.
- `day_groups` / `compare_days`: a candidate against a baseline on one shared mask, with whole days
  resampled for both series at once, so the pairing is kept and one weather system stays in one
  block.

Matrices are `[station, time]`, with `NaN` for missing.
"""
module DayBlockBootstrap

using CSV, DataFrames, Dates, Random, Statistics

export read_wide_values, read_wide_mask, group_sums, sum_stats, day_groups, compare_days,
    threshold_label

"""`2.0` -> `"2"`, `0.1` -> `"0.1"`: the suffix a threshold's columns carry."""
threshold_label(threshold) = isinteger(threshold) ? string(Int(threshold)) : string(threshold)

function _read_wide(path, ids, times, convert)
    table = CSV.read(path, DataFrame)
    DateTime.(string.(table.time)) == times || error("$path: its hours differ from the grid's")
    return permutedims(reduce(hcat, [convert(table[!, id]) for id in ids]))
end

"""An OOF value table as a stations × hours matrix; a missing value becomes `NaN`."""
read_wide_values(path, ids, times) =
    _read_wide(path, ids, times, column -> Float64.(coalesce.(column, NaN)))

"""An OOF mask table (`true`/`false`) as a stations × hours `BitMatrix`; missing is `false`."""
read_wide_mask(path, ids, times) =
    BitMatrix(_read_wide(path, ids, times, column -> coalesce.(column, false) .== true))

_n_sum(thresholds) = 8 + 3 * length(thresholds)

"""
Sums per group over the masked cells:
- n, Σe, Σe², Σx, Σy, Σx², Σy² and Σxy, where e = x − y;
- then hits, misses and false alarms per threshold.

`group[i, t]` is the cell's group, and 0 skips the cell.
"""
function group_sums(x, y, mask, group::AbstractMatrix{<:Integer}, n_group::Int,
                    thresholds::AbstractVector{<:Real})
    sums = zeros(n_group, _n_sum(thresholds))
    @inbounds for index in eachindex(y)
        mask[index] || continue
        g = group[index]
        g == 0 && continue
        xv, yv = x[index], y[index]
        e = xv - yv
        sums[g, 1] += 1
        sums[g, 2] += e
        sums[g, 3] += e * e
        sums[g, 4] += xv
        sums[g, 5] += yv
        sums[g, 6] += xv * xv
        sums[g, 7] += yv * yv
        sums[g, 8] += xv * yv
        for (j, threshold) in enumerate(thresholds)
            observed = yv >= threshold
            estimated = xv >= threshold
            column = 8 + 3 * (j - 1)
            sums[g, column + 1] += observed & estimated
            sums[g, column + 2] += observed & !estimated
            sums[g, column + 3] += !observed & estimated
        end
    end
    return sums
end

"""The scores one row of `group_sums` gives, the event ones at every threshold."""
function sum_stats(total, thresholds::AbstractVector{<:Real})
    n = total[1]
    mx, my = total[4] / n, total[5] / n
    r = (total[8] / n - mx * my) / sqrt((total[6] / n - mx^2) * (total[7] / n - my^2))
    bias = total[2] / n
    stats = Pair{Symbol,Float64}[:n => n, :RMSE => sqrt(total[3] / n), :bias => bias,
        :abs_bias => abs(bias), :r => r]
    for (j, threshold) in enumerate(thresholds)
        column = 8 + 3 * (j - 1)
        hits, misses, false_alarms = total[column + 1], total[column + 2], total[column + 3]
        tag = threshold_label(threshold)
        fbi = (hits + false_alarms) / (hits + misses)
        append!(stats, [Symbol("POD_", tag) => hits / (hits + misses),
            Symbol("FAR_", tag) => false_alarms / (hits + false_alarms),
            Symbol("CSI_", tag) => hits / (hits + misses + false_alarms),
            Symbol("FBI_", tag) => fbi, Symbol("FBI_", tag, "_dev") => abs(fbi - 1)])
    end
    return (; stats...)
end

"""Each cell's day as a group index (a stations × hours matrix), and the number of days."""
function day_groups(times::AbstractVector{DateTime}, n_station::Int)
    days = Date.(times)
    unique_days = sort(unique(days))
    day_of = [searchsortedfirst(unique_days, day) for day in days]
    return repeat(permutedims(day_of), n_station), length(unique_days)
end

_prefixed(stats, prefix) =
    NamedTuple{Tuple(Symbol(prefix, name) for name in keys(stats))}(values(stats))

"""
Candidate against baseline, from per-day `group_sums` over one shared mask.

Returns both series' scores (`candidate_*`, `baseline_*`), the relative change in RMSE, and the
point change in every `ci_stats` entry, each with a 95% paired day-block bootstrap interval.
"""
function compare_days(candidate::Matrix{Float64}, baseline::Matrix{Float64},
                      thresholds::AbstractVector{<:Real}; ci_stats, reps::Int, seed::Int)
    keep = findall(>(0), candidate[:, 1])
    keep == findall(>(0), baseline[:, 1]) || error("the two series are not on one mask")
    c = sum_stats(vec(sum(candidate; dims=1)), thresholds)
    b = sum_stats(vec(sum(baseline; dims=1)), thresholds)
    rng = MersenneTwister(seed)
    draws = Dict(stat => Vector{Float64}(undef, reps) for stat in ci_stats)
    for rep in 1:reps
        sampled = rand(rng, keep, length(keep))
        cs = sum_stats(vec(sum(view(candidate, sampled, :); dims=1)), thresholds)
        bs = sum_stats(vec(sum(view(baseline, sampled, :); dims=1)), thresholds)
        for stat in ci_stats
            draws[stat][rep] = getproperty(cs, stat) - getproperty(bs, stat)
        end
    end
    fields = Pair{Symbol,Float64}[:rel_RMSE => (c.RMSE - b.RMSE) / b.RMSE]
    if :RMSE in ci_stats
        append!(fields, [:rel_RMSE_lo => quantile(draws[:RMSE], 0.025) / b.RMSE,
            :rel_RMSE_hi => quantile(draws[:RMSE], 0.975) / b.RMSE])
    end
    for stat in ci_stats
        append!(fields, [Symbol("delta_", stat) => getproperty(c, stat) - getproperty(b, stat),
            Symbol("delta_", stat, "_lo") => quantile(draws[stat], 0.025),
            Symbol("delta_", stat, "_hi") => quantile(draws[stat], 0.975)])
    end
    return merge(_prefixed(c, "candidate_"), _prefixed(b, "baseline_"), (; fields...))
end

end # module

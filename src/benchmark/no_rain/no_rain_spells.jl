# ---------------------------------------------------------------------------------------------
# Dry spells and dry days

"""
Runs of consecutive `true` in `flags` over the positions where `valid` holds.

`breaks[t]` marks a gap in time before position `t` (a non-contiguous axis), which ends a run as a
missing position does. Returns `(; lengths, censored)`: a run is censored when it touches the ends of
the axis, a missing position or a break, since its true length is then unknown.
"""
function run_lengths(flags::AbstractVector{Bool}, valid::AbstractVector{Bool};
        breaks::AbstractVector{Bool}=falses(length(flags)))
    length(flags) == length(valid) == length(breaks) || throw(DimensionMismatch("flags, valid and breaks must match"))
    lengths, censored = Int[], Bool[]
    n = length(flags)
    t = 1
    while t <= n
        if valid[t] && flags[t]
            start = t
            while t < n && valid[t + 1] && flags[t + 1] && !breaks[t + 1]
                t += 1
            end
            left_open = start == 1 || breaks[start] || !valid[start - 1]
            right_open = t == n || breaks[t + 1] || !valid[t + 1]
            push!(lengths, t - start + 1)
            push!(censored, left_open || right_open)
        end
        t += 1
    end
    return (; lengths, censored)
end

"""
Dry spells of one source: runs of cells below `threshold` inside `valid`.

Returns `(; summary, histogram)`. `summary` has a pooled `all` row and one row per station: dry share
of the valid cells, the number of complete and censored spells, and the mean, median, 90th percentile
and maximum length of the complete ones, plus the share of all valid cells lying in complete spells of
at least 24 and 72 steps. `histogram` bins the pooled complete spells by `bin_edges` (lower bounds):
spell count and the cells they hold. On a daily axis the same numbers are in days.
"""
function dry_spell_table(
    Y::AbstractMatrix, valid::AbstractMatrix; threshold::Real, station_ids::AbstractVector,
    breaks::AbstractVector{Bool}=falses(size(Y, 2)), bin_edges=SPELL_BIN_EDGES, key::NamedTuple=(;),
)
    size(Y) == size(valid) || throw(DimensionMismatch("Y and valid must match"))
    summarize(level, lengths, censored, n_valid, n_dry) = begin
        complete = lengths[.!censored]
        (; level, n_valid, dry_share=_ratio(n_dry, n_valid), n_spells=length(complete), n_censored=count(censored),
            mean_length=isempty(complete) ? NaN : mean(complete), median_length=isempty(complete) ? NaN : median(complete),
            p90_length=isempty(complete) ? NaN : quantile(complete, 0.9), max_length=isempty(complete) ? 0 : maximum(complete),
            share_in_ge24=_ratio(sum(filter(>=(24), complete); init=0), n_valid),
            share_in_ge72=_ratio(sum(filter(>=(72), complete); init=0), n_valid))
    end
    rows = NamedTuple[]
    pooled_lengths, pooled_censored = Int[], Bool[]
    pooled_valid, pooled_dry = 0, 0
    for i in axes(Y, 1)
        v = BitVector(view(valid, i, :) .& isfinite.(view(Y, i, :)))
        flags = BitVector([v[j] && Y[i, j] < threshold for j in axes(Y, 2)])
        runs = run_lengths(flags, v; breaks)
        n_valid, n_dry = count(v), count(flags)
        push!(rows, merge(key, summarize(String(station_ids[i]), runs.lengths, runs.censored, n_valid, n_dry)))
        append!(pooled_lengths, runs.lengths)
        append!(pooled_censored, runs.censored)
        pooled_valid += n_valid
        pooled_dry += n_dry
    end
    pushfirst!(rows, merge(key, summarize("all", pooled_lengths, pooled_censored, pooled_valid, pooled_dry)))
    complete = pooled_lengths[.!pooled_censored]
    hist = NamedTuple[]
    for (b, lower) in enumerate(bin_edges)
        upper = b < length(bin_edges) ? bin_edges[b + 1] : typemax(Int)
        in_bin = filter(x -> lower <= x < upper, complete)
        push!(hist, merge(key, (; bin_lower=lower, bin_upper=upper == typemax(Int) ? -1 : upper,
            n_spells=length(in_bin), cells_in_spells=sum(in_bin; init=0),
            spell_share=_ratio(length(in_bin), length(complete)),
            cell_share=_ratio(sum(in_bin; init=0), sum(complete; init=0)))))
    end
    return (; summary=DataFrame(rows), histogram=DataFrame(hist))
end

"""
Daily sums over the `valid` cells of each met-day: `totals[station, day]` for `day = 1:maximum(day)`,
`NaN` where fewer than `min_hours` of that station-day's hours are valid. Summing gauge and estimate over
one shared `valid` mask keeps a partial day fair: both sides then cover the same hours.
"""
function daily_sums(day::AbstractVector{<:Integer}, Y::AbstractMatrix, valid::AbstractMatrix; min_hours::Int=24)
    size(Y) == size(valid) || throw(DimensionMismatch("Y and valid must match"))
    length(day) == size(Y, 2) || throw(DimensionMismatch("one day per column"))
    D = maximum(day)
    totals, hours = zeros(size(Y, 1), D), zeros(Int, size(Y, 1), D)
    for j in axes(Y, 2), i in axes(Y, 1)
        (valid[i, j] && isfinite(Y[i, j])) || continue
        totals[i, day[j]] += Y[i, j]
        hours[i, day[j]] += 1
    end
    totals[hours .< min_hours] .= NaN
    return totals
end

"""
Runs of dry days at a gauge during which its neighbours were clearly wet: the signature of a gauge that
reads 0.0 while it is not measuring (a blocked funnel, a stopped logger), not of a dry spell.

`D[station, day]` holds daily totals (`NaN` when incomplete). Every run of at least `min_days` consecutive
complete days below `dry_mm` is compared, day by day, with the median of the `k` nearest stations reporting
that day; the run is flagged when that median reaches `wet_mm` on at least `min_wet_days` of its days.
`p_dry_given_wet` - how often any gauge stays below `dry_mm` on a day its neighbours' median reaches
`wet_mm` - is returned alongside, as the single-day rate the rule rests on: several such days inside one
dry run are very unlikely to be weather.

Returns `(; spells, flagged, p_dry_given_wet)`: one row per run (`station_id`, `first_day`, `last_day`
as indices into the day axis, `n_days`, `gauge_total`, `neighbour_total`, `neighbour_wet_days`,
`flagged`), and `flagged[station, day]`, true on every day of a flagged run.
"""
function silent_gauge_spells(D::AbstractMatrix, lonlat::AbstractMatrix, station_ids::AbstractVector;
        dry_mm::Real=1.0, wet_mm::Real=5.0, min_days::Int=7, min_wet_days::Int=3, k::Int=3)
    n, n_days = size(D)
    size(lonlat, 1) == n == length(station_ids) || throw(DimensionMismatch("one lon/lat row and id per station"))
    distance = haversine_distance_matrix(lonlat, lonlat)
    order = [filter(!=(i), sortperm(view(distance, :, i))) for i in 1:n]
    neighbour = fill(NaN, n, n_days)
    for d in 1:n_days, i in 1:n
        values = Float64[]
        for m in order[i]
            isfinite(D[m, d]) && push!(values, D[m, d])
            length(values) == k && break
        end
        isempty(values) || (neighbour[i, d] = median(values))
    end
    paired = [(D[c] < dry_mm) for c in eachindex(D) if isfinite(D[c]) && isfinite(neighbour[c]) && neighbour[c] >= wet_mm]
    flagged = falses(n, n_days)
    rows = NamedTuple[]
    for i in 1:n
        d = 1
        while d <= n_days
            if isfinite(D[i, d]) && D[i, d] < dry_mm
                first_day = d
                while d < n_days && isfinite(D[i, d + 1]) && D[i, d + 1] < dry_mm
                    d += 1
                end
                if d - first_day + 1 >= min_days
                    around = filter(isfinite, neighbour[i, first_day:d])
                    wet_days = count(>=(wet_mm), around)
                    flag = wet_days >= min_wet_days
                    flag && (flagged[i, first_day:d] .= true)
                    push!(rows, (; station_id=String(station_ids[i]), first_day, last_day=d, n_days=d - first_day + 1,
                        gauge_total=sum(D[i, first_day:d]), neighbour_total=sum(around; init=0.0),
                        neighbour_wet_days=wet_days, flagged=flag))
                end
            end
            d += 1
        end
    end
    return (; spells=DataFrame(rows), flagged, p_dry_given_wet=isempty(paired) ? NaN : mean(paired))
end

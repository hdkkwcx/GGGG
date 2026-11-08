"""
Cross-validation that holds out whole blocks of time as well as gauges.

Holding out gauges alone leaves one weather system on both sides of a split: a gauge in training and
a held-out gauge share the same hour, and rain or false rain often covers the network at once. A
model that can learn one episode's signature is then "confirmed" by that episode's other gauges.

- `time_block_folds`: every hour's block and time fold. Blocks are whole days, and a block's fold is
  drawn so each stratum (say regime × season) spreads evenly over the folds.
- `training_hours`: the hours a time fold may train on, leaving out a guard around its held-out days,
  so an episode that crosses a block's edge does not stay on both sides.
- `crossed_splits`: every pair of gauge fold and time fold. Each scores its held-out gauges in its
  held-out hours, so each gauge-hour is scored once, by a model that saw neither.

Hour labels are hour ends, so an hour's day is the day of `label - 1 h`.

- `block_times` / `block_days_option`: for the day-block bootstraps, which group hours by `Date(t)`.
  Handed `block_times` in place of the hours, they resample whole blocks of days instead.
- `from_option` / `from_suffix`: a script's `--from YYYY-MM`, the month its record starts.
"""
module TimeBlockFolds

using Dates, Random

export label_day, time_block_folds, training_hours, crossed_splits, block_times, block_days_option,
    from_option, from_suffix

"""The calendar day of an hour labelled by its end."""
label_day(time::DateTime) = Date(time - Hour(1))

"""
Each hour's block and time fold, as `(; block_of_time, fold_of_block, fold_of_time)`.

Blocks are `block_days` long, counted from the first day, and those with no hours are skipped. A
block's stratum is that of its first hour. Within each stratum, in sorted order, the shuffled blocks
are dealt to the folds in turn, with the deal carrying over between strata so the totals stay even.
"""
function time_block_folds(times::AbstractVector{DateTime}, strata::AbstractVector;
                          k::Integer=5, block_days::Integer=7, seed::Integer)
    length(times) == length(strata) || throw(DimensionMismatch("one stratum per hour is needed"))
    days = label_day.(times)
    first_day = minimum(days)
    raw_block = [Dates.value(day - first_day) ÷ block_days for day in days]
    blocks = sort(unique(raw_block))
    block_of_time = [searchsortedfirst(blocks, block) for block in raw_block]
    stratum_of_block = [strata[findfirst(==(b), block_of_time)] for b in eachindex(blocks)]
    fold_of_block = zeros(Int, length(blocks))
    rng = MersenneTwister(seed)
    dealt = 0
    for stratum in sort(unique(stratum_of_block))
        for block in shuffle(rng, findall(==(stratum), stratum_of_block))
            dealt += 1
            fold_of_block[block] = mod1(dealt, k)
        end
    end
    return (; block_of_time, fold_of_block, fold_of_time=fold_of_block[block_of_time])
end

"""
Whether each hour may train a model scored on time fold `fold`: not on a held-out day, nor within
`guard_days` of one.
"""
function training_hours(fold_of_time::AbstractVector{<:Integer}, times::AbstractVector{DateTime},
                        fold::Integer; guard_days::Integer=1)
    days = label_day.(times)
    excluded = Set(day + Day(gap) for day in days[fold_of_time .== fold]
        for gap in -guard_days:guard_days)
    return BitVector([!(day in excluded) for day in days])
end

"""
Every pair of gauge fold and time fold, as NamedTuples of `gauge_fold`, `time_fold`,
`train_stations`, `test_stations`, `train_hours` and `test_hours` (all indices).
"""
function crossed_splits(fold_of_station::AbstractVector{<:Integer},
                        fold_of_time::AbstractVector{<:Integer},
                        times::AbstractVector{DateTime}; guard_days::Integer=1)
    return [(; gauge_fold, time_fold,
        train_stations=findall(!=(gauge_fold), fold_of_station),
        test_stations=findall(==(gauge_fold), fold_of_station),
        train_hours=findall(training_hours(fold_of_time, times, time_fold; guard_days)),
        test_hours=findall(==(time_fold), fold_of_time))
        for gauge_fold in sort(unique(fold_of_station)) for time_fold in sort(unique(fold_of_time))]
end

"""
Each hour as the start of its `block_days`-day block, for the code that groups hours by `Date(t)`:
`_daily_bootstrap_delta`, `daily_bootstrap_event_delta` and `DayBlockBootstrap.day_groups`. Handed
these in place of the hours, that code resamples whole blocks instead of days.

Blocks count `Date(t)` from the first hour's date, as that code does, not `label_day`. With
`block_days = 1` every hour maps to the start of its own date, so the grouping is unchanged.
"""
function block_times(times::AbstractVector{DateTime}; block_days::Integer=7)
    block_days >= 1 || throw(ArgumentError("block_days must be at least 1"))
    first_day = minimum(Date.(times))
    return [DateTime(first_day + Day(fld(Dates.value(Date(t) - first_day), block_days) * block_days))
        for t in times]
end

"""
`--block-days N` read out of command-line `args`, as `(N, the other args)`. `N` is `nothing` without
the flag, so a script can keep its default path then, and still run `--block-days 1` through
`block_times` as a check that it reproduces that default.
"""
function block_days_option(args)
    rest = collect(String, args)
    index = findfirst(==("--block-days"), rest)
    index === nothing && return nothing, rest
    index < length(rest) || throw(ArgumentError("--block-days needs a number of days"))
    block_days = parse(Int, rest[index + 1])
    block_days >= 1 || throw(ArgumentError("--block-days must be at least 1"))
    return block_days, deleteat!(rest, index:(index + 1))
end

"""
`--from YYYY-MM` read out of command-line `args`, as `(first day of that month, the other args)`,
or `(nothing, args)` without the flag. A script starts its 08-08 days there instead of at its
default start, so its day blocks and per-record counts start there too; `from_suffix` names the
output folder after it.
"""
function from_option(args)
    rest = collect(String, args)
    index = findfirst(==("--from"), rest)
    index === nothing && return nothing, rest
    text = get(rest, index + 1, "")
    parsed = match(r"^(\d{4})-(\d{2})$", text)
    parsed === nothing && throw(ArgumentError("--from expects YYYY-MM, got \"$text\""))
    return Date(parse(Int, parsed[1]), parse(Int, parsed[2]), 1), deleteat!(rest, index:(index + 1))
end

"""`_from<yyyymm>` for a `from_option` date, `""` for `nothing`."""
from_suffix(from::Union{Nothing,Date}) =
    from === nothing ? "" : "_from" * Dates.format(from, "yyyymm")

end # module

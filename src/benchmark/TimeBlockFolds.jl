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
"""
module TimeBlockFolds

using Dates, Random

export label_day, time_block_folds, training_hours, crossed_splits

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

end # module

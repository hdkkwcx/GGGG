"""
A two-stage FY4B correction by counting on training gauges: first whether it rains, then how much.

`SatelliteCalibration`'s quantile map fixes how much FY4B rains. It is monotone within each
season-by-regime cell, though, so it barely changes which hours rank as wetter, and hourly r moved
little (0.223 plain, 0.246 adaptive). This module re-ranks hours by their context. Each
station-hour gets a class from FY4B alone, at the fold's pixel shift (`shift_fold`):

- **intensity:** the shifted 3x3 mean, on `TS_INTENSITY_EDGES` (8 classes, dry below 0.1);
- **rain nearby:** how many of the 8 shifts around the chosen one are wet (none, some or all),
  read from the same offset grid. Together the nine 3x3 boxes cover a 5x5-pixel area.
- **rain in the adjacent hours:** FY4B wet at t-1 and at t+1 (neither, one or both). Only an hour
  exactly 1 h away counts, and a missing hour counts as dry. Using the next hour is fine for a
  retrospective correction; the fused anchor's lead term does the same.

Then, per season-by-regime cell and class, counted on the fold's training gauges only:

1. **occurrence:** P(gauge >= wet | class), the share of wet training gauge-hours;
2. **amount:** the class's mean gauge rain over its wet hours.

A class with fewer than `TS_MIN_PAIRS` training pairs in its cell takes, in turn, the same class
pooled over its regime, over every cell, then its intensity class alone (`POOL_LEVELS`).

There are two predictions:

- `expected_anchor` gives P x amount, the conditional mean. It serves r and RMSE, and it flattens
  heavy rain.
- `occurrence_anchor` calls the classes with the highest P wet, in each cell, until they cover as
  many training gauge-hours as the gauges' own wet count, stopping before a class that would
  leave the count further from it. This is the quantile map's frequency rule. A wet hour keeps the calibrated amount it already has, or takes the class's mean wet
  amount when that is 0; every other hour is 0. That removes isolated false rain and keeps the
  heavy tail.

Every choice above (the edges, the classes, `TS_MIN_PAIRS`) is fixed in advance and not tuned on
gauges, and nothing is fitted: the tables are counts.

`hurdle_gwr` (`HurdleGWR`) is the repository's other two-part model. Its occurrence part is a
spatially weighted wet frequency per satellite-intensity class and month, with no context to
re-rank by, so none of it is reused here.

Leakage is guarded as in `SatelliteCalibration`: the shift and the tables come from the `train`
stations and are applied to every station, so no held-out gauge shapes the value it is scored
against.
"""
module SatelliteTwoStage

using Dates
using Main.SatelliteCalibration
using Main: SatelliteCalibration

export TS_INTENSITY_EDGES, TS_MIN_PAIRS, N_CONTEXT, POOL_LEVELS
export fy4b_intensity_class, context_class, context_label
export neighbour_wet_count, fold_nearby, context_classes, fit_context_tables
export expected_anchor, occurrence_anchor, context_table_rows, two_stage_fold

"""Upper edges (mm/h) of the intensity classes, below 0.1 being dry; fixed in advance."""
const TS_INTENSITY_EDGES = [0.1, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]
const N_INTENSITY = length(TS_INTENSITY_EDGES) + 1
const NEARBY_LABELS = ("none", "some", "all")
const ADJACENT_LABELS = ("neither", "one", "both")
"""Context classes per season-by-regime cell: intensity x rain nearby x rain at t-1/t+1."""
const N_CONTEXT = N_INTENSITY * length(NEARBY_LABELS) * length(ADJACENT_LABELS)

"""Training pairs a class needs at one pooling level to be counted there; fixed in advance."""
const TS_MIN_PAIRS = 200

"""Where a class's counts came from, finest first; `pooled` is every training pair at once."""
const POOL_LEVELS = ("cell", "regime", "all", "intensity", "pooled")

"""Intensity class of a FY4B value: 1 below 0.1 mm/h (dry), up to 8 at 16 mm/h and above."""
fy4b_intensity_class(x::Real) = searchsortedlast(TS_INTENSITY_EDGES, x) + 1

"""Context class (1 to `N_CONTEXT`) of intensity class `i`, nearby class `j` and `m` wet adjacent hours."""
context_class(i::Integer, j::Integer, m::Integer) = (i - 1) * 9 + (j - 1) * 3 + m + 1

_number(x) = isinteger(x) ? string(Int(x)) : string(x)

"""What context class `q` means, as `(; intensity, nearby, adjacent)` labels."""
function context_label(q::Integer)
    i = (q - 1) ÷ 9 + 1
    j = ((q - 1) % 9) ÷ 3 + 1
    m = (q - 1) % 3
    intensity = i == 1 ? "dry" : i == N_INTENSITY ? ">=$(_number(TS_INTENSITY_EDGES[end]))" :
        "$(_number(TS_INTENSITY_EDGES[i - 1]))-$(_number(TS_INTENSITY_EDGES[i]))"
    return (; intensity, nearby=NEARBY_LABELS[j], adjacent=ADJACENT_LABELS[m + 1])
end

"""
How many of the 8 shifts around shift `k` are wet at every station and hour, as
`(count, n_neighbour)`: a `UInt8` `[station, time]` matrix, and how many of the 8 lie inside the
grid. A shift's 3x3 mean is wet at `>= wet`, and a gap counts as dry.
"""
function neighbour_wet_count(grid::OffsetGrid, k::Integer; wet::Real=0.1)
    dx, dy = grid.shifts[k]
    index = Dict(shift => i for (i, shift) in enumerate(grid.shifts))
    count = zeros(UInt8, length(grid.station_rows), length(grid.time_cols))
    n_neighbour = 0
    for b in -1:1, a in -1:1
        (a, b) == (0, 0) && continue
        neighbour = get(index, (dx + a, dy + b), nothing)
        neighbour === nothing && continue
        n_neighbour += 1
        series = read_shift(grid, neighbour)
        @inbounds for i in eachindex(series)
            series[i] >= wet && (count[i] += 0x01)
        end
    end
    return count, n_neighbour
end

"""
The nearby class of every station and hour under a fold's chosen shifts (`shift_fold`): 1 when
none of the surrounding shifts is wet, 3 when all are, 2 otherwise, and 0 in a regime with no
chosen shift. `cache` keeps each centre shift's counts, so folds that agree read nothing twice.
"""
function fold_nearby(grid::OffsetGrid, chosen::AbstractVector{<:Integer},
                     regime_of_time::AbstractVector{<:Integer}; wet::Real=0.1,
                     cache::AbstractDict=Dict{Int,Tuple{Matrix{UInt8},Int}}())
    nearby = zeros(UInt8, length(grid.station_rows), length(regime_of_time))
    for (g, k) in enumerate(chosen)
        k == 0 && continue
        count, n_neighbour = get!(() -> neighbour_wet_count(grid, k; wet), cache, k)
        for t in findall(==(g), regime_of_time), s in axes(nearby, 1)
            c = count[s, t]
            nearby[s, t] = c == 0 ? 1 : c >= n_neighbour ? 3 : 2
        end
    end
    return nearby
end

"""
The context class of every station and hour, as an `Int` `[station, time]` matrix, 0 where the
shifted value is missing. `times` must be sorted. An adjacent hour counts only when it is exactly
1 h away, and then only when its shifted value is `>= wet`.
"""
function context_classes(shifted::AbstractMatrix{<:Real}, nearby::AbstractMatrix{<:Integer},
                         times::AbstractVector{DateTime}; wet::Real=0.1)
    size(shifted) == size(nearby) && size(shifted, 2) == length(times) ||
        throw(DimensionMismatch("shifted, nearby and times must describe the same stations x hours"))
    issorted(times) || throw(ArgumentError("times must be sorted"))
    classes = zeros(Int, size(shifted))
    n_time = length(times)
    for t in 1:n_time
        has_prev = t > 1 && times[t - 1] == times[t] - Hour(1)
        has_next = t < n_time && times[t + 1] == times[t] + Hour(1)
        for s in axes(shifted, 1)
            x = shifted[s, t]
            isnan(x) && continue
            adjacent = (has_prev && shifted[s, t - 1] >= wet) + (has_next && shifted[s, t + 1] >= wet)
            classes[s, t] = context_class(fy4b_intensity_class(x), nearby[s, t], adjacent)
        end
    end
    return classes
end

_regime(cell) = (cell - 1) ÷ 4 + 1

"""
Occurrence and amount per season-by-regime cell (`calibration_cells`) and context class, counted
on the `train` stations only.

Returns `[cell, class]` matrices: `p` (P(gauge `>= wet`)), `amount` (mean gauge rain over the wet
pairs, NaN where there are none), `pool` (index into `POOL_LEVELS`), the counts `n` and `n_wet` at
that level and `n_own`, `n_wet_own` in the cell itself. Also `p_star`, per cell, the lowest P that
`occurrence_anchor` still calls wet, and `cells`, the cells present.
"""
function fit_context_tables(y_obs::AbstractMatrix{<:Real}, classes::AbstractMatrix{<:Integer},
                            cell_of_time::AbstractVector{<:Integer},
                            train::AbstractVector{<:Integer}; wet::Real=0.1,
                            min_pairs::Integer=TS_MIN_PAIRS)
    size(y_obs) == size(classes) && size(y_obs, 2) == length(cell_of_time) ||
        throw(DimensionMismatch("y_obs, classes and cell_of_time must describe the same grid"))
    n_cell = maximum(cell_of_time)
    n_own = zeros(Int, n_cell, N_CONTEXT)
    n_wet_own = zeros(Int, n_cell, N_CONTEXT)
    sum_wet = zeros(n_cell, N_CONTEXT)
    for t in axes(y_obs, 2)
        cell = cell_of_time[t]
        for s in train
            q = classes[s, t]
            y = y_obs[s, t]
            (q == 0 || isnan(y)) && continue
            n_own[cell, q] += 1
            if y >= wet
                n_wet_own[cell, q] += 1
                sum_wet[cell, q] += y
            end
        end
    end
    sum(n_own) > 0 || error("no training pairs with a context class to count")
    cells = sort(unique(cell_of_time))
    same_intensity(q) = [r for r in 1:N_CONTEXT if (r - 1) ÷ 9 == (q - 1) ÷ 9]
    p = fill(NaN, n_cell, N_CONTEXT)
    amount = fill(NaN, n_cell, N_CONTEXT)
    pool = zeros(Int, n_cell, N_CONTEXT)
    n = zeros(Int, n_cell, N_CONTEXT)
    n_wet = zeros(Int, n_cell, N_CONTEXT)
    for cell in cells, q in 1:N_CONTEXT
        regime_cells = [c for c in 1:n_cell if _regime(c) == _regime(cell)]
        levels = ((cell:cell, q:q), (regime_cells, q:q), (1:n_cell, q:q),
            (1:n_cell, same_intensity(q)), (1:n_cell, 1:N_CONTEXT))
        level = something(findfirst(l -> sum(view(n_own, l...)) >= min_pairs, levels),
            length(levels))
        rows, cols = levels[level]
        pool[cell, q] = level
        n[cell, q] = sum(view(n_own, rows, cols))
        n_wet[cell, q] = sum(view(n_wet_own, rows, cols))
        p[cell, q] = n_wet[cell, q] / n[cell, q]
        n_wet[cell, q] > 0 && (amount[cell, q] = sum(view(sum_wet, rows, cols)) / n_wet[cell, q])
    end
    # Per cell, call classes wet in order of falling P until the training pairs they cover reach
    # the gauges' own wet count there; the P of the last class called wet is the cell's cut. A
    # class that would leave the count further from the target than it already is stays dry: in
    # winter FY4B's wet classes cover fewer hours than the gauges are wet, and calling the next
    # class (dry, nothing nearby, a third of a million pairs) wet to close that gap would call
    # every dry hour wet.
    p_star = fill(Inf, n_cell)
    for cell in cells
        target = sum(view(n_wet_own, cell, :))
        target == 0 && continue
        candidates = [q for q in 1:N_CONTEXT if n_own[cell, q] > 0 && p[cell, q] > 0]
        sort!(candidates; by=q -> (-p[cell, q], -q))
        covered = 0
        for q in candidates
            covered + n_own[cell, q] - target > target - covered && break
            covered += n_own[cell, q]
            p_star[cell] = p[cell, q]
            covered >= target && break
        end
    end
    return (; p, amount, pool, n, n_wet, n_own, n_wet_own, p_star, cells)
end

_called_wet(tables, cell, q) = (p = tables.p[cell, q]; p > 0 && p >= tables.p_star[cell])

_expected(tables, cell, q) = (p = tables.p[cell, q]; p == 0 ? 0.0 : p * tables.amount[cell, q])

"""P x amount at every station and hour; NaN where there is no context class."""
expected_anchor(tables, classes::AbstractMatrix{<:Integer}, cell_of_time) =
    [classes[s, t] == 0 ? NaN : _expected(tables, cell_of_time[t], classes[s, t])
     for s in axes(classes, 1), t in axes(classes, 2)]

"""
The occurrence rule at every station and hour: in a class called wet, the `calibrated` value
(the fold's `FY4B_QM`), or the class's mean wet amount where that is 0; elsewhere 0. NaN where
there is no context class or no calibrated value.
"""
function occurrence_anchor(tables, classes::AbstractMatrix{<:Integer}, cell_of_time,
                           calibrated::AbstractMatrix{<:Real})
    size(calibrated) == size(classes) ||
        throw(DimensionMismatch("calibrated must have the classes' stations x hours"))
    return [begin
            q = classes[s, t]
            cell = cell_of_time[t]
            x = calibrated[s, t]
            q == 0 || isnan(x) ? NaN :
                !_called_wet(tables, cell, q) ? 0.0 :
                x > 0 ? Float64(x) : tables.amount[cell, q]
        end for s in axes(classes, 1), t in axes(classes, 2)]
end

"""One row per cell present and context class: its labels, counts, pooling level, P, amount and
whether the occurrence rule calls it wet."""
function context_table_rows(tables, regimes::AbstractVector{<:AbstractString})
    rows = NamedTuple[]
    for cell in tables.cells, q in 1:N_CONTEXT
        push!(rows, (; regime=regimes[_regime(cell)],
            season=SatelliteCalibration.SEASONS[(cell - 1) % 4 + 1], class=q, context_label(q)...,
            n_own=tables.n_own[cell, q], n_wet_own=tables.n_wet_own[cell, q],
            pool=POOL_LEVELS[tables.pool[cell, q]], n=tables.n[cell, q],
            n_wet=tables.n_wet[cell, q], p_wet=tables.p[cell, q],
            wet_amount=tables.amount[cell, q], expected=_expected(tables, cell, q),
            p_star=tables.p_star[cell], called_wet=_called_wet(tables, cell, q)))
    end
    return rows
end

"""
One fold's two-stage anchors for every station, as `(expected, occurrence, table_rows)`.

The shift is chosen per regime on `train` (`shift_fold`, sharing `cache` with `calibrate_fold`),
the classes are read at it, and the tables are counted on `train`. `calibrated` is the fold's
`FY4B_QM` anchor, whose amounts the occurrence rule keeps. `nearby_cache` keeps each centre
shift's neighbour counts.
"""
function two_stage_fold(grid::OffsetGrid, sums, y_obs::AbstractMatrix{<:Real},
                        train::AbstractVector{<:Integer}, regimes::AbstractVector{<:AbstractString},
                        regime_of_time::AbstractVector{<:Integer}, cell_of_time,
                        times::AbstractVector{DateTime}, calibrated::AbstractMatrix{<:Real};
                        cache::AbstractDict=Dict{Int,Matrix{Float64}}(),
                        nearby_cache::AbstractDict=Dict{Int,Tuple{Matrix{UInt8},Int}}(),
                        wet::Real=0.1, min_pairs::Integer=TS_MIN_PAIRS)
    shifted, chosen, _ = shift_fold(grid, sums, train, regimes, regime_of_time; cache)
    nearby = fold_nearby(grid, chosen, regime_of_time; wet, cache=nearby_cache)
    classes = context_classes(shifted, nearby, times; wet)
    tables = fit_context_tables(y_obs, classes, cell_of_time, train; wet, min_pairs)
    return expected_anchor(tables, classes, cell_of_time),
        occurrence_anchor(tables, classes, cell_of_time, calibrated),
        context_table_rows(tables, regimes)
end

end # module

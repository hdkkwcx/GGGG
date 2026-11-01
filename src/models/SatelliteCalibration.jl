"""
Calibrate one satellite product against gauges, inside a cross-validation fold: a pixel shift
and a wet-threshold quantile map, both fitted on the fold's training stations only.

Built for FY4B, whose rebuilt series (`FY4BPreprocessing.aggregate_fy4b_variants`) still carries
two faults that no gauge-free step can remove:

- **A displacement.** Moving the sample pixel and scoring against the gauges shows a systematic
  shift that is not parallax-shaped: 2-5 px west in 2022 (software V1.0, 133E), none consistent
  in 2023 (V1.0.0, 133E), 1-4 px east and mostly north in 2024 (105E). A fixed-height parallax
  shift made June 2022 worse. The shift can only be measured against gauges, so it is chosen per
  regime from the rebuild's `fy4b_offset_grid.nc`, by the highest pooled correlation over the
  training stations.
- **Amount and occurrence errors.** Biases of +170% (August 2022) to -50% (June 2024), and
  false rain in dry hours. Per season and regime, a threshold sets FY4B's wet frequency to the
  training gauges' and a quantile map gives its wet amounts their distribution. Above the top
  knot values are scaled by the top knots' ratio, which keeps the product's 20 mm/h cap a lower
  bound instead of a ceiling. The full map helped June but raised July and August's error, so an
  `adaptive` calibration applies the map, month by month, only as strongly as the training
  gauges support.

Leakage is the point of the `train` arguments, as in `SatelliteFusion`: the benchmark passes a
fold's training stations, applies the result to every station, and so no held-out gauge shapes
the anchor it is scored against. Fitting on all 237 gauges instead moved r by about 0.001
(measured 2026-10-30), so the guard costs nothing.
"""
module SatelliteCalibration

using CSV, DataFrames, Dates, NCDatasets, Statistics

export OffsetGrid, open_offset_grid, read_shift, shift_sums, choose_shift, shift_fold
export product_regimes, calibration_cells, strength_cells, fit_quantile_map, apply_quantile_map
export QM_STRENGTHS, choose_strength, calibrate_fold

"""Quantile knots of the map; the outermost 0.5% on each side is left to the tail rules."""
const QM_PROBABILITIES = collect(0.005:0.005:0.995)

"""
Wet pairs a season-by-regime cell needs on each side to get its own map, fixed in advance.
A sparser cell takes its regime's map, and a sparse regime the pooled one.
"""
const QM_MIN_WET = 500

season_index(t) = month(t) in (12, 1, 2) ? 1 : month(t) in (3, 4, 5) ? 2 :
    month(t) in (6, 7, 8) ? 3 : 4
const SEASONS = ("DJF", "MAM", "JJA", "SON")

_parse_label(label) = DateTime(first(string(label), 19), dateformat"yyyy-mm-ddTHH:MM:SS")

"""
Regime of every time, as `(labels, regime_of_time)`, from the rebuild's `fy4b_hour_meta.csv`.

One regime per distinct (satellite longitude, software version), so either changing starts a new
one; read from the files rather than hard-coded dates.
"""
function product_regimes(meta_path::AbstractString, times::AbstractVector{DateTime})
    meta = CSV.read(meta_path, DataFrame; types=Dict(:time => String, :software => String))
    label_of = Dict(_parse_label(row.time) => "$(row.sat_lon)E_$(row.software)"
        for row in eachrow(meta))
    missing_times = [t for t in times if !haskey(label_of, t)]
    isempty(missing_times) || error(
        "$(length(missing_times)) benchmark hours have no regime in $meta_path, " *
        "e.g. $(first(missing_times))")
    labels = sort(unique(label_of[t] for t in times))
    index = Dict(label => i for (i, label) in enumerate(labels))
    return labels, [index[label_of[t]] for t in times]
end

"""Season-by-regime cell of every time: `(regime - 1) * 4 + season`."""
calibration_cells(times, regime_of_time) =
    [(regime - 1) * 4 + season_index(t) for (t, regime) in zip(times, regime_of_time)]

"""
Month cell of every time, as `(labels, cell_of_time)`: one per regime and calendar month of one
year, labelled by the month (`"2022-07"`); a month split between two regimes gets two cells.

An adaptive map chooses its strength per month cell, finer than the season cell its map is fitted
on: in one summer the full map helped June and hurt July and August.
"""
function strength_cells(times::AbstractVector{DateTime}, regime_of_time::AbstractVector{<:Integer})
    keys = [(regime, year(t), month(t)) for (t, regime) in zip(times, regime_of_time)]
    unique_keys = sort(unique(keys))
    index = Dict(key => i for (i, key) in enumerate(unique_keys))
    labels = ["$(y)-$(lpad(m, 2, '0'))" for (_, y, m) in unique_keys]
    return labels, [index[key] for key in keys]
end

"""
The rebuild's shifted grid, aligned to a benchmark's stations and hours.

Shift `k` runs over `offset_x` fastest, as `read_shift` and `shifts` order them.
"""
struct OffsetGrid
    path::String
    shifts::Vector{Tuple{Int,Int}}
    n_x::Int
    station_rows::Vector{Int}
    time_cols::Vector{Int}
end

function open_offset_grid(path::AbstractString, ids::AbstractVector{<:AbstractString},
                          times::AbstractVector{DateTime})
    NCDataset(path) do ds
        offset_x = Int.(ds["offset_x"][:])
        offset_y = Int.(ds["offset_y"][:])
        station_rows = indexin(ids, string.(ds["station_id"][:]))
        any(isnothing, station_rows) && error("stations missing from $path: " *
            join(ids[isnothing.(station_rows)], ", "))
        time_cols = indexin(times, _parse_label.(ds["time_label"][:]))
        any(isnothing, time_cols) && error("$(count(isnothing, time_cols)) benchmark hours " *
            "missing from $path")
        shifts = [(dx, dy) for dy in offset_y for dx in offset_x]
        return OffsetGrid(String(path), shifts, length(offset_x),
            Int.(station_rows), Int.(time_cols))
    end
end

"""Shift `k`'s series as a benchmark-aligned `[station, time]` matrix."""
function read_shift(grid::OffsetGrid, k::Integer)
    a = (k - 1) % grid.n_x + 1
    b = (k - 1) ÷ grid.n_x + 1
    NCDataset(grid.path) do ds
        values = ds["fy4b"][a, b, :, :]
        return Float64.(values[grid.station_rows, grid.time_cols])
    end
end

"""
Per shift, station and regime, the sums a pooled correlation needs, as `[k, station, regime, s]`
with `s` = n, Σx, Σy, Σx², Σy², Σxy over hours where both values are finite.

Computed once for all stations; a fold then adds up its training stations' rows, so choosing a
shift costs no further reads.
"""
function shift_sums(grid::OffsetGrid, y_obs::AbstractMatrix{<:Real},
                    regime_of_time::AbstractVector{<:Integer}, n_regime::Integer)
    n_station, n_time = size(y_obs)
    sums = zeros(length(grid.shifts), n_station, n_regime, 6)
    for k in eachindex(grid.shifts)
        x = read_shift(grid, k)
        @inbounds for time in 1:n_time, station in 1:n_station
            xv = x[station, time]
            yv = y_obs[station, time]
            (isnan(xv) || isnan(yv)) && continue
            g = regime_of_time[time]
            sums[k, station, g, 1] += 1
            sums[k, station, g, 2] += xv
            sums[k, station, g, 3] += yv
            sums[k, station, g, 4] += xv * xv
            sums[k, station, g, 5] += yv * yv
            sums[k, station, g, 6] += xv * yv
        end
    end
    return sums
end

function _pooled_r(s)
    n = s[1]
    n < 2 && return NaN
    mx, my = s[2] / n, s[3] / n
    vx, vy = s[4] / n - mx^2, s[5] / n - my^2
    (vx <= 0 || vy <= 0) && return NaN
    return (s[6] / n - mx * my) / sqrt(vx * vy)
end

"""
The shift with the highest pooled correlation over `train` stations in regime `g`, as
`(k, r, r_unshifted)`. Shifts are tried nearest first, so a tie keeps the smaller one, and a
regime with no usable pairs keeps no shift.
"""
function choose_shift(grid::OffsetGrid, sums::AbstractArray{<:Real,4},
                      train::AbstractVector{<:Integer}, g::Integer)
    order = sortperm(grid.shifts; by=s -> (s[1]^2 + s[2]^2, s[2], s[1]))
    zero_k = findfirst(==((0, 0)), grid.shifts)
    pooled(k) = _pooled_r(vec(sum(view(sums, k, train, g, :); dims=1)))
    best_k, best_r = zero_k, -Inf
    for k in order
        r = pooled(k)
        isnan(r) && continue
        r > best_r && ((best_k, best_r) = (k, r))
    end
    return best_k, best_r, pooled(zero_k)
end

"""
The fold's shifted series for every station, as `(shifted, chosen, shift_rows)`.

Per regime, the shift `choose_shift` picks on `train` is read from the grid and used for that
regime's hours; `chosen[g]` is regime `g`'s shift index, `0` for a regime with no hours here.
`cache` keeps the shifts already read, since folds tend to agree.
"""
function shift_fold(grid::OffsetGrid, sums, train::AbstractVector{<:Integer},
                    regimes::AbstractVector{<:AbstractString},
                    regime_of_time::AbstractVector{<:Integer};
                    cache::AbstractDict=Dict{Int,Matrix{Float64}}())
    shifted = fill(NaN, length(grid.station_rows), length(regime_of_time))
    chosen = zeros(Int, length(regimes))
    shift_rows = NamedTuple[]
    for (g, label) in enumerate(regimes)
        cols = findall(==(g), regime_of_time)
        isempty(cols) && continue
        k, r, r_unshifted = choose_shift(grid, sums, train, g)
        series = get!(() -> read_shift(grid, k), cache, k)
        shifted[:, cols] = series[:, cols]
        chosen[g] = k
        dx, dy = grid.shifts[k]
        push!(shift_rows, (; regime=label, dx, dy, r_train=r, r_train_unshifted=r_unshifted,
            n_hours=length(cols)))
    end
    return shifted, chosen, shift_rows
end

function _fit_cell(y_obs, y_sat, train, cols, wet)
    x = Float64[]
    y = Float64[]
    for time in cols, station in train
        xv = y_sat[station, time]
        yv = y_obs[station, time]
        (isnan(xv) || isnan(yv)) && continue
        push!(x, xv)
        push!(y, yv)
    end
    isempty(x) && return nothing
    wet_fraction = mean(y .>= wet)
    tau = wet_fraction > 0 ? quantile(x, 1 - wet_fraction) : Inf
    x_wet = sort(x[x .> tau])
    y_wet = sort(y[y .>= wet])
    (isempty(x_wet) || isempty(y_wet)) && return nothing
    return (; tau, qx=quantile(x_wet, QM_PROBABILITIES; sorted=true),
        qy=quantile(y_wet, QM_PROBABILITIES; sorted=true),
        n_pairs=length(x), n_wet_sat=length(x_wet), n_wet_obs=length(y_wet))
end

"""
A wet-threshold quantile map per cell, fitted on the `train` stations only, as
`Dict(cell => map)`; each map also says which pool it was fitted on.

`tau` gives FY4B the training gauges' frequency of `>= wet`; above it, FY4B's wet values are
mapped onto the gauges' wet values knot by knot. A cell with fewer than `min_wet` wet pairs on
either side takes its regime's map, and a sparse regime the pooled map.
"""
function fit_quantile_map(y_obs::AbstractMatrix{<:Real}, y_sat::AbstractMatrix{<:Real},
                          train::AbstractVector{<:Integer}, cell_of_time::AbstractVector{<:Integer};
                          wet::Real=0.1, min_wet::Integer=QM_MIN_WET)
    enough(m) = m !== nothing && min(m.n_wet_sat, m.n_wet_obs) >= min_wet
    pooled = _fit_cell(y_obs, y_sat, train, eachindex(cell_of_time), wet)
    pooled === nothing && error("no finite training pairs to fit a quantile map on")
    maps = Dict{Int,Any}()
    for cell in sort(unique(cell_of_time))
        regime = (cell - 1) ÷ 4 + 1
        own = _fit_cell(y_obs, y_sat, train, findall(==(cell), cell_of_time), wet)
        if enough(own)
            maps[cell] = merge(own, (; pool="cell"))
            continue
        end
        regime_cols = findall(c -> (c - 1) ÷ 4 + 1 == regime, cell_of_time)
        regime_map = _fit_cell(y_obs, y_sat, train, regime_cols, wet)
        maps[cell] = enough(regime_map) ? merge(regime_map, (; pool="regime")) :
            merge(pooled, (; pool="pooled"))
    end
    return maps
end

function _map_value(m, x)
    isnan(x) && return NaN
    x <= m.tau && return 0.0
    qx, qy = m.qx, m.qy
    x >= qx[end] && return qy[end] * x / max(qx[end], eps())
    x <= qx[1] && return qy[1]
    k = clamp(searchsortedlast(qx, x), 1, length(qx) - 1)
    width = qx[k + 1] - qx[k]
    return width > 0 ? qy[k] + (x - qx[k]) / width * (qy[k + 1] - qy[k]) : qy[k]
end

"""Apply `fit_quantile_map`'s maps to every station and hour of `y_sat`."""
apply_quantile_map(maps, y_sat::AbstractMatrix{<:Real}, cell_of_time) =
    [_map_value(maps[cell_of_time[time]], y_sat[station, time])
     for station in axes(y_sat, 1), time in axes(y_sat, 2)]

"""
How far an adaptive calibration may move the wet amounts towards the map, smallest first.

Candidate 0 leaves the value alone; candidate `c >= 1` zeroes values at or below `tau` and moves
the rest `QM_STRENGTHS[c]` of the way to the mapped value, so the first is the threshold alone and
the last is the full map. Every candidate is monotone in the satellite value.
"""
const QM_STRENGTHS = (0.0, 0.25, 0.5, 0.75, 1.0)

function _adjusted_value(x, mapped, tau, candidate)
    candidate == 0 && return x
    isnan(x) && return NaN
    x <= tau && return 0.0
    strength = QM_STRENGTHS[candidate]
    strength == 1 && return mapped
    return (1 - strength) * x + strength * mapped
end

"""
Per choice cell, the candidate (see `QM_STRENGTHS`) with the lowest squared error on the `train`
stations, as `Dict(choice => (; candidate, n_pairs, rmse_none, rmse_full, rmse_chosen))`.

Each pair is adjusted with the map of its own cell (`cell_of_time`) and scored in its choice cell
(`choice_of_time`), so the choice can be finer than the map: per month (`strength_cells`) while
the map is per season. The full map restores the gauges' distribution, which raises the error
wherever FY4B places rain badly (July and August, where it put the gauge tail onto misplaced
convective cells). A tie goes to the smaller adjustment, and a cell with no training pairs keeps
the full map.
"""
function choose_strength(y_obs::AbstractMatrix{<:Real}, y_sat::AbstractMatrix{<:Real}, maps,
                         train::AbstractVector{<:Integer}, cell_of_time::AbstractVector{<:Integer},
                         choice_of_time::AbstractVector{<:Integer}=cell_of_time)
    n_candidate = length(QM_STRENGTHS) + 1
    choices = unique(choice_of_time)
    sse = Dict(choice => zeros(n_candidate) for choice in choices)
    n = Dict(choice => 0 for choice in choices)
    for time in axes(y_sat, 2)
        choice = choice_of_time[time]
        m = maps[cell_of_time[time]]
        s = sse[choice]
        for station in train
            x = y_sat[station, time]
            y = y_obs[station, time]
            (isnan(x) || isnan(y)) && continue
            n[choice] += 1
            mapped = _map_value(m, x)
            for c in 0:(n_candidate - 1)
                s[c + 1] += (_adjusted_value(x, mapped, m.tau, c) - y)^2
            end
        end
    end
    full = length(QM_STRENGTHS)
    return Dict(map(choices) do choice
        rmse = sqrt.(sse[choice] ./ n[choice])
        candidate = n[choice] == 0 ? full : argmin(rmse) - 1
        choice => (; candidate, n_pairs=n[choice], rmse_none=rmse[1], rmse_full=rmse[end],
            rmse_chosen=rmse[candidate + 1])
    end)
end

"""
One fold's calibrated anchor for every station, as `(anchor, shift_rows, map_rows, strength_rows)`.

Per regime, the shift is chosen on `train` (`shift_fold`) and that shift's series is used for
the regime's hours at every station; the quantile map is then fitted on `train` and applied
everywhere. `cache` keeps the shifts already read, since folds tend to agree.

`adaptive` applies the map only as strongly as `choose_strength` finds best on `train` in each
choice cell (`choice_of_time`, labelled by `choice_labels`; the benchmark passes month cells from
`strength_cells`), down to leaving the shifted series alone; without it every hour gets the full
map. Either way the strength rows report each choice cell's training error with no adjustment,
with the full map and with what was applied.
"""
function calibrate_fold(grid::OffsetGrid, sums, y_obs::AbstractMatrix{<:Real},
                        train::AbstractVector{<:Integer}, regimes::AbstractVector{<:AbstractString},
                        regime_of_time::AbstractVector{<:Integer}, cell_of_time;
                        cache::AbstractDict=Dict{Int,Matrix{Float64}}(), wet::Real=0.1,
                        adaptive::Bool=false, choice_of_time::AbstractVector{<:Integer}=cell_of_time,
                        choice_labels=nothing)
    size(y_obs) == (length(grid.station_rows), length(regime_of_time)) ||
        throw(DimensionMismatch("y_obs must be the grid's stations x the run's hours"))
    shifted, _, shift_rows = shift_fold(grid, sums, train, regimes, regime_of_time; cache)
    maps = fit_quantile_map(y_obs, shifted, train, cell_of_time; wet)
    map_rows = [(; regime=regimes[(cell - 1) ÷ 4 + 1], season=SEASONS[(cell - 1) % 4 + 1],
        pool=m.pool, tau=m.tau, n_pairs=m.n_pairs, n_wet_sat=m.n_wet_sat, n_wet_obs=m.n_wet_obs,
        q50_sat=m.qx[100], q50_obs=m.qy[100], q99_sat=m.qx[198], q99_obs=m.qy[198])
        for (cell, m) in sort(collect(maps); by=first)]
    strengths = choose_strength(y_obs, shifted, maps, train, cell_of_time, choice_of_time)
    full = length(QM_STRENGTHS)
    applied(choice) = adaptive ? strengths[choice].candidate : full
    strength_rows = map(sort(collect(keys(strengths)))) do choice
        first_time = findfirst(==(choice), choice_of_time)
        s = strengths[choice]
        (; regime=regimes[regime_of_time[first_time]],
            cell=choice_labels === nothing ? SEASONS[(choice - 1) % 4 + 1] : choice_labels[choice],
            n_pairs=s.n_pairs, threshold_applied=applied(choice) >= 1,
            strength=applied(choice) == 0 ? 0.0 : QM_STRENGTHS[applied(choice)],
            train_rmse_none=s.rmse_none, train_rmse_full=s.rmse_full,
            train_rmse_chosen=adaptive ? s.rmse_chosen : s.rmse_full)
    end
    adaptive || return apply_quantile_map(maps, shifted, cell_of_time), shift_rows, map_rows,
        strength_rows
    anchor = [begin
            m = maps[cell_of_time[time]]
            x = shifted[station, time]
            _adjusted_value(x, _map_value(m, x), m.tau, strengths[choice_of_time[time]].candidate)
        end for station in axes(shifted, 1), time in axes(shifted, 2)]
    return anchor, shift_rows, map_rows, strength_rows
end

end # module

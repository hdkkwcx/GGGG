"""
A near-surface humidity gate for a satellite anchor: set its light rain to 0 where the air at the gauge
is too dry for rain to reach the ground.

FY4B's rain is a lookup on one infrared channel, so cold cloud over dry air reads as rain. The screen
(`scripts/screen_fy4b_humidity.jl`) found that ERA5-Land relative humidity at the gauge separates much
of that false rain from real rain. This is its gate "B8":
- the cut is the RH in `GATE_CUTS` below which zeroing the anchor most raises CSI at the wet threshold
  on the training cells, or no gate (`-Inf`) when no cut raises it;
- values of at least `GATE_SPARED_MM` are never touched, so detection at 8 mm/h cannot fall.

`time_crossed_gate` fits the cut per time fold on the other blocks' hours (`TimeBlockFolds`), so a week's
anchor never comes from a cut that saw that week. ERA5-Land is reanalysis, not FY4B: a gated product is
"FY4B + reanalysis".

Matrices are `[station, hour]`, with `NaN` for missing.
"""
module HumidityGate

using Dates
using Main.TimeBlockFolds: training_hours
using Main.ERA5LandProcessing: relative_humidity

export GATE_CUTS, GATE_SPARED_MM, station_humidity, choose_gate_cut, apply_gate, time_crossed_gate

const GATE_CUTS = 40.0:5.0:90.0
const GATE_SPARED_MM = 8.0

"""Relative humidity (%) from 2 m temperature and dew point in °C, element by element."""
station_humidity(t2m_c::AbstractMatrix, d2m_c::AbstractMatrix) =
    relative_humidity.(t2m_c .+ 273.15, d2m_c .+ 273.15)

"""Whether the gate zeroes `x` at humidity `rh` under `cut`. `NaN` in either never gates."""
_gated(x, rh, cut, spared) = rh < cut && x < spared

"""
The cut in `cuts` that most raises CSI at `wet` over the cells `stations` x `hours` where `x`, `y`
and `rh` are finite, when `x` is zeroed wherever `_gated`. Returns `(; cut, csi_before, csi_after,
n_cells, n_wet, removed)`, with `cut = -Inf` (no gate) when no cut raises CSI; a tie keeps the lower cut.
"""
function choose_gate_cut(x::AbstractMatrix, y::AbstractMatrix, rh::AbstractMatrix,
                         stations::AbstractVector{<:Integer}, hours::AbstractVector{<:Integer};
                         wet::Real=0.1, spared::Real=GATE_SPARED_MM, cuts=GATE_CUTS)
    hits = false_alarms = misses = n_cells = 0
    removable = Tuple{Float64,Bool}[]           # (rh, is a hit) for wet cells the gate could zero
    for t in hours, i in stations
        xv, yv, rv = x[i, t], y[i, t], rh[i, t]
        (isfinite(xv) && isfinite(yv) && isfinite(rv)) || continue
        n_cells += 1
        observed, estimated = yv >= wet, xv >= wet
        if estimated
            observed ? (hits += 1) : (false_alarms += 1)
            xv < spared && push!(removable, (rv, observed))
        elseif observed
            misses += 1
        end
    end
    denominator = hits + misses + false_alarms
    before = denominator > 0 ? hits / denominator : NaN
    best_cut, best_csi, best_removed = -Inf, before, 0
    for cut in cuts
        removed_hits = count(cell -> cell[1] < cut && cell[2], removable)
        removed_false = count(cell -> cell[1] < cut && !cell[2], removable)
        remaining = denominator - removed_false
        remaining > 0 || continue
        csi = (hits - removed_hits) / remaining
        csi > best_csi && ((best_cut, best_csi, best_removed) = (cut, csi, removed_hits + removed_false))
    end
    return (; cut=Float64(best_cut), csi_before=before, csi_after=best_csi, n_cells,
        n_wet=hits + false_alarms, removed=best_removed)
end

"""`x` with every cell `_gated` under `cut` set to 0; with `cut = -Inf`, an unchanged copy."""
function apply_gate(x::AbstractMatrix, rh::AbstractMatrix, cut::Real; spared::Real=GATE_SPARED_MM)
    gated = Matrix{Float64}(x)
    for index in eachindex(gated)
        _gated(gated[index], rh[index], cut, spared) && (gated[index] = 0.0)
    end
    return gated
end

"""
`x` gated per time fold, at every station: each fold's hours take the cut `choose_gate_cut` picks on
`train_stations` over the hours `training_hours` leaves that fold (its own days and `guard_days` either
side held out). Returns `(gated, rows)`, one row per time fold with its cut and training counts.
"""
function time_crossed_gate(x::AbstractMatrix, y::AbstractMatrix, rh::AbstractMatrix,
                           train_stations::AbstractVector{<:Integer},
                           fold_of_time::AbstractVector{<:Integer}, times::AbstractVector{DateTime};
                           guard_days::Integer=1, wet::Real=0.1, spared::Real=GATE_SPARED_MM,
                           cuts=GATE_CUTS)
    size(x) == size(y) == size(rh) || throw(DimensionMismatch("x, y and rh must have one shape"))
    length(fold_of_time) == length(times) == size(x, 2) ||
        throw(DimensionMismatch("one time fold and one time per hour are needed"))
    gated = Matrix{Float64}(x)
    rows = NamedTuple[]
    for fold in sort(unique(fold_of_time))
        train_hours = findall(training_hours(fold_of_time, times, fold; guard_days))
        choice = choose_gate_cut(x, y, rh, train_stations, train_hours; wet, spared, cuts)
        test_hours = findall(==(fold), fold_of_time)
        for t in test_hours, i in axes(gated, 1)
            _gated(gated[i, t], rh[i, t], choice.cut, spared) && (gated[i, t] = 0.0)
        end
        push!(rows, merge((; time_fold=fold, n_train_hours=length(train_hours),
            n_test_hours=length(test_hours)), choice))
    end
    return gated, rows
end

end # module

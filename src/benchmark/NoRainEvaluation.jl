"""
No-rain evaluation: what the satellite products, and the interpolation and fusion methods built on
them, report on the station-hours where the gauge records no rain.

About 93% of gauge station-hours are dry, so dry-hour behaviour decides most pooled scores, yet
`SatelliteTemporalEvaluation.false_alarm_table` is the only place it has been looked at on its own.
This module widens that view:

- `occurrence_table`: the full occurrence contingency and the dry-hour amount scores of one estimate,
  per level of any stratifier, with day-block bootstrap intervals. Every stratification below (season,
  hour, station, distance to gauge rain, neighbouring gauges, the network's state, meteorology, product
  agreement) goes through this one function, so they are all scored the same way.
- Stratifiers of the gauge-dry hours: `hours_to_rain` / `proximity_labels` / `side_labels` (where a dry
  hour sits relative to the station's own rain), `neighbour_wet_share` / `neighbour_labels` (whether
  nearby gauges are wet), `network_wet_share` / `network_labels` (whether the whole network is dry),
  and `bin_labels` for continuous covariates.
- `neighbour_gauge_baseline`: how often a *neighbouring gauge* reports rain when a station is dry, by
  separation - the false-alarm rate a perfect areal estimate would still show against a point gauge.
- `run_lengths` / `dry_spell_table` / `daily_sums`: dry spells, hourly and daily.
- `silent_gauge_spells`: the other side of a false alarm - long runs of gauge zeros during which the
  neighbouring gauges were clearly wet, i.e. a gauge reporting 0.0 while it was not measuring.
- `paired_error_table` / `zero_threshold_sensitivity`: dry-hour error of one prediction against a
  reference, and what zeroing small values would change.
- `long_to_wide_hourly` / `lagged_correlation`: bring a long-format product onto the gauge grid and
  check its clock.

Conventions, shared with `SatelliteTemporalEvaluation`:

- Matrices are `[station, time]`, `NaN` for missing, every estimate aligned to the gauges.
- `mask` names the cells a sample scores; a cell outside it is ignored everywhere.
- A gauge-dry cell is `gauge < WET_MM` (0.1 mm/h). The gauge record is quantised to 0.5 mm, so that is a
  0.0 reading, and a satellite value between 0.1 and 0.5 mm/h may be rain below that resolution.
- A stratifier is an integer label per cell - a matrix, or a function `(i, j) -> label` for labels that
  only depend on the station or the hour. Labels run 1..L, and 0 leaves a cell out.
- `day[j]` is column `j`'s met-day number (positive integers); intervals resample whole days.
"""
module NoRainEvaluation

using DataFrames, Dates, Random, Statistics
using Main.SatelliteTemporalEvaluation: WET_MM, GAUGE_RESOLUTION_MM, day_block_bootstrap
using Main.TraditionalInterpolation: haversine_distance_matrix

export DRY_THRESHOLDS, EXCEEDANCE_THRESHOLDS, HOURS_PER_YEAR
export PROXIMITY_EDGES, PROXIMITY_LEVELS, SIDE_LEVELS, NEIGHBOUR_LEVELS, NETWORK_EDGES, NETWORK_LEVELS
export SPELL_BIN_EDGES
export hours_to_rain, proximity_labels, side_labels, neighbour_wet_share, neighbour_labels
export network_wet_share, network_labels, bin_labels, phantom_extent
export occurrence_table, paired_error_table, neighbour_gauge_baseline
export run_lengths, dry_spell_table, daily_sums, silent_gauge_spells
export long_to_wide_hourly, lagged_correlation, zero_threshold_sensitivity

"""Estimate thresholds for occurrence: the wet threshold and the gauge's 0.5 mm resolution."""
const DRY_THRESHOLDS = [WET_MM, GAUGE_RESOLUTION_MM]
"""Thresholds of the exceedance curve P(estimate >= t | gauge dry)."""
const EXCEEDANCE_THRESHOLDS = [0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0, 4.0, 8.0]
const HOURS_PER_YEAR = 8766.0
"""Upper bounds (h) of the distance-to-rain classes; the last class is open."""
const PROXIMITY_EDGES = [1, 3, 6, 24]
const PROXIMITY_LEVELS = ["1 h", "2-3 h", "4-6 h", "7-24 h", "> 24 h"]
const SIDE_LEVELS = ["inside an event gap", "before rain", "after rain", "away from rain"]
const NEIGHBOUR_LEVELS = ["all neighbours dry", "under half wet", "half or more wet"]
"""Upper bounds of the network wet share for the classes after `network dry`."""
const NETWORK_EDGES = [0.05, 0.25]
const NETWORK_LEVELS = ["network dry", "isolated (<= 5% wet)", "scattered (5-25% wet)", "widespread (> 25% wet)"]
"""Lower bounds (h, or days on a daily axis) of the spell-length histogram bins."""
const SPELL_BIN_EDGES = [1, 2, 3, 4, 6, 8, 12, 16, 24, 36, 48, 72, 96, 144, 192, 288, 384, 576, 768]

_ratio(a, b) = b > 0 ? a / b : NaN
_label(labels::AbstractMatrix, i, j) = labels[i, j]
_label(labels::Function, i, j) = labels(i, j)

# The evaluation, one file per section in `no_rain/`; each fragment opens with its own banner.
include("no_rain/no_rain_context.jl")
include("no_rain/no_rain_occurrence.jl")
include("no_rain/no_rain_spells.jl")
include("no_rain/no_rain_inputs.jl")

end

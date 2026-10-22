module DEMTerrainExperiment

using CSV
using DataFrames
using Dates
using LinearAlgebra
using Random
using Statistics
using MixedGWR: INTENSITY_RANGES, haversine_distance_matrix, tied_ranks, gw_weight!
using MixedGWR: BISQUARE
# The mixed/multiscale GWR estimators live in MixedGWR; re-exported for the existing callers.
using MixedGWR: mixed_gwr_predict, multiscale_gwr_predict
using MixedGWR: select_mixed_bandwidth, select_multiscale_bandwidths

export DEMExperimentConfig, bh_adjust, terrain_screen, spatial_variability_test
export mixed_gwr_predict, multiscale_gwr_predict, select_mixed_bandwidth
export select_multiscale_bandwidths, run_dem_experiment
export mean_wet_residual, monthly_correlation_rows, terrain_model_designs
export terrain_groups, terrain_columns

const TERRAIN_GROUPS = ("elevation", "slope", "aspect")
const TERRAIN_COLUMNS = Dict(
    "elevation" => [:elevation_m],
    "slope" => [:slope_deg],
    "aspect" => [:aspect_sin, :aspect_cos],
)

terrain_groups() = collect(TERRAIN_GROUPS)
terrain_columns(group::String) = copy(TERRAIN_COLUMNS[group])

Base.@kwdef struct DEMExperimentConfig
    outdir::String
    wet_threshold::Float64 = 0.1
    min_wet_hours::Int = 100
    k::Int = 5
    seed::Int = 20260815
    bandwidth_candidates::Vector{Int} = [30, 50, 80, 120, 160]
    screen_permutations::Int = 999
    spatial_permutations::Int = 999
    q_threshold::Float64 = 0.05
    vif_threshold::Float64 = 5.0
    ridge::Float64 = 1e-8
    tolerance::Float64 = 1e-5
    max_iterations::Int = 200
end

# The terrain study, split by concern in `dem/`. `experiment.jl` is the study orchestrator that
# writes the fixed-name CSVs; the files before it are the screening and smoother-based spatial
# variability tests it is built on.

include("dem/screening.jl")
include("dem/smoothers.jl")
include("dem/spatial_variability.jl")
include("dem/summaries.jl")
include("dem/experiment.jl")

end

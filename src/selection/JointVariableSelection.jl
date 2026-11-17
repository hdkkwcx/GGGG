module JointVariableSelection

using CSV
using DataFrames
using Dates
using LinearAlgebra
using Random
using Statistics
using MixedGWR: gw_weight, haversine_distance_matrix, loocv_distances, BISQUARE

# Shared bookkeeping for every variable-selection path; see src/selection/SelectionScaffolding.jl.
using Main.SelectionScaffolding

# All three are loaded through src/load_modules.jl before this file, so each is shared
# rather than recompiled into this module. See src/load_modules.jl.
using Main.DEMTerrainExperiment
using Main.ERA5VariableSelection
using Main.NDVIVariableSelection

# The covariate design vocabulary, shared with `JointCovariateModels` so the two cannot drift.
using Main.CovariateGroups: JOINT_GROUP_ORDER, GROUP_COLUMNS, GROUP_FAMILY

export JointSelectionConfig, JOINT_GROUP_ORDER, prepare_joint_panel, joint_vif
export joint_spatial_variability_test, run_joint_variable_selection, select_joint_covariates


_predictor_label(group::String) = group == "ndvi" ? "ndvi_qc" :
    join(String.(GROUP_COLUMNS[group]), "+")

Base.@kwdef struct JointSelectionConfig
    outdir::String
    wet_threshold::Float64 = 0.1
    min_wet_hours::Int = 100
    min_stations_per_time::Int = 12
    min_ndvi_wet_hours_per_period::Int = 5
    min_ndvi_periods_per_station::Int = 8
    max_ndvi_age_days::Int = 32
    k::Int = 5
    seed::Int = 20260818
    bandwidth_candidates::Vector{Int} = [30, 50, 80, 120, 160]
    independent_permutations::Int = 999
    spatial_permutations::Int = 999
    q_threshold::Float64 = 0.05
    vif_threshold::Float64 = 5.0
    ridge::Float64 = 1e-8
end

function _family_configs(cfg::JointSelectionConfig)
    era = ERA5SelectionConfig(
        outdir=cfg.outdir, wet_threshold=cfg.wet_threshold,
        min_wet_hours=cfg.min_wet_hours,
        min_stations_per_time=cfg.min_stations_per_time,
        k=cfg.k, seed=20260816,
        bandwidth_candidates=cfg.bandwidth_candidates,
        association_permutations=cfg.independent_permutations,
        spatial_permutations=cfg.spatial_permutations,
        q_threshold=cfg.q_threshold, vif_threshold=cfg.vif_threshold,
        ridge=cfg.ridge,
    )
    ndvi = NDVISelectionConfig(
        outdir=cfg.outdir, wet_threshold=cfg.wet_threshold,
        min_wet_hours_per_period=cfg.min_ndvi_wet_hours_per_period,
        min_periods_per_station=cfg.min_ndvi_periods_per_station,
        min_stations_per_period=cfg.min_stations_per_time,
        max_age_days=cfg.max_ndvi_age_days,
        k=cfg.k, seed=20260817,
        bandwidth_candidates=cfg.bandwidth_candidates,
        association_permutations=cfg.independent_permutations,
        spatial_permutations=cfg.spatial_permutations,
        q_threshold=cfg.q_threshold, ridge=cfg.ridge,
    )
    return era, ndvi
end

# The selection, split by step in `joint/`: the independent screen of each family
# (`joint_candidates.jl`), the panel the screened groups share (`joint_panel.jl`), the
# cross-family VIF (`joint_vif.jl`), the local/global role test (`joint_panel_variability.jl`), and
# `select_joint_covariates` with the study run that writes the CSVs (`joint_selection_run.jl`).
include("joint/joint_candidates.jl")
include("joint/joint_panel.jl")
include("joint/joint_vif.jl")
include("joint/joint_panel_variability.jl")
include("joint/joint_selection_run.jl")

end

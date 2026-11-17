module JointCovariateModels

using CSV
using DataFrames
using Dates
using LinearAlgebra
using SHA
using Statistics

using MixedGWR: GWKernel, haversine_distance_matrix, local_hat_matrix, local_weight_matrix
using MixedGWR: global_projection
using MixedGWR: mixed_gwr_predict

# The covariate design vocabulary, shared with `JointVariableSelection` so the two cannot drift.
using Main.CovariateGroups: JOINT_GROUP_ORDER, GROUP_COLUMNS, GROUP_FAMILY

export JointCovariateBenchmarkConfig, JointFoldContext
export load_joint_covariate_spec, joint_spec_sha256, build_joint_fold_context
export dynamic_covariate_predict, joint_effective_roles, joint_group_names
export joint_models_coincide


Base.@kwdef struct JointCovariateBenchmarkConfig
    # `nothing` when the benchmark is configured for nested per-fold selection instead of a
    # fixed full-data specification (see `InterpolationBenchmarkConfig.joint_selection`).
    spec_path::Union{Nothing,String} = nothing
    terrain_path::String
    era5_paths::Dict{Int,String}
    ndvi_path::Union{Nothing,String} = nothing
    bandwidth_candidates::Vector{Int} = [30, 50, 80, 120, 160]
    feature_time_offset_hours::Int = 9
    max_ndvi_age_days::Int = 32
    wet_threshold::Float64 = 0.1
    ridge::Float64 = 1e-8
    tolerance::Float64 = 1e-5
    # The back-fit routinely needs several hundred sweeps once smooth covariate groups are in
    # play; at the old cap of 200 the median sweep count sat at 150-194 and MGWR discarded whole
    # hours for want of a few more iterations. Measured on the real `local_hat_matrix` and
    # `global_projection` operators (150 stations, 93%-dry residuals, two smooth ERA5-like
    # covariate groups): the `:split` layout converges on 60% of hours at 200 sweeps and 100% at
    # 500.
    max_iterations::Int = 1000
    # Successive over-relaxation factor for the back-fitting sweeps. It converges to the same
    # fixed point for any admissible value - `relaxation` cancels at the fixed point - so this
    # trades convergence *rate* only.
    #
    # DO NOT "fix" this to 1.0. The name `_multiscale_predict_damped` suggests damping and a
    # value above 1 looks like a bug, but plain Gauss-Seidel is far worse here: on the fixture
    # above, the share of hours converging within 200 sweeps runs
    #
    #   layout            w=1.0  w=1.1  w=1.2  w=1.3  w=1.5  w=1.7
    #   :split               4%    12%    40%    98%    60%     1%
    #   :shared              3%     8%    22%   100%   100%    52%
    #   :intercept_only    100%   100%   100%   100%   100%   100%
    #
    # i.e. over-relaxation is what makes the coupled layouts usable at all. 1.5 is kept as the
    # historical value; 1.3 looked better on that one synthetic fixture, which is not enough
    # evidence to move a default, and the layout that wins on accuracy is insensitive to it.
    relaxation::Float64 = 1.5
    # How MGWR splits the spatial trend into back-fitting groups. Covariate groups always get
    # their own bandwidth, which is what makes the model multiscale; this controls only the
    # `[1, z_lon, z_lat]` block.
    #
    #   :split           - one single-column group each for the intercept, longitude and latitude
    #                      (the historical layout).
    #   :shared          - one three-column group, solved by a single weighted least squares the
    #                      way `mixed_gwr` does.
    #   :intercept_only  - drop the coordinate columns; a locally varying intercept plus the
    #                      covariate groups.
    #
    # The coordinate groups are near-constant inside any local window, so they are near-collinear
    # with the intercept group *and* with any spatially smooth covariate. That coupling, not the
    # relaxation factor, is what governs how many sweeps the back-fit needs.
    #
    # `:intercept_only` is the default on measured evidence. Against a common-cell paired
    # comparison of the smoke run (`scripts/compare_benchmark_runs.jl --smoke`, where the
    # traditional methods move +0.0% and so confirm the pairing), MGWR's RMSE change is
    #
    #   layout            FY4B     GPM   GSMaP   coverage at 1.000   smoke wall clock
    #   :split           +0.1%   -6.8%   -5.1%             8 / 15             30-58 min
    #   :shared          +0.9%   -5.7%   -3.5%            13 / 15                16 min
    #   :intercept_only  -4.3%   -6.5%   -6.4%            15 / 15                14 min
    #
    # i.e. it wins on two products, is 0.3pp behind on the third against a layout scored on
    # fewer cells, and is the only one that never discards an hour. Dropping the coordinate plane
    # also frees the intercept to go genuinely local (8-80 neighbours) instead of pinning near the
    # grid maximum, which is the multiscale behaviour the method is supposed to have.
    mgwr_spatial_grouping::Symbol = :intercept_only
    # How a local group's fit at a target with too few supporting stations is recorded.
    #
    # `:missing` marks that target NaN, so the coverage machinery (`min_tuning_coverage`,
    # `_prediction_coverage`, `_common_method_mask`) can see it and reject the candidate.
    # `:zero` is the historical all-zero hat row, which reports a fabricated "no local
    # correction" as a successful prediction and is invisible to every one of those gates:
    # measured on the 2026-09-02 full nested run, the gate rejected 195/675 of direct `gwr`'s
    # candidates and 0/675, 0/675, 0/3939 of `residual_gwr`/`mixed_gwr`/`mgwr`'s on the
    # identical grid, and `mgwr`'s reported GPM own-coverage of 0.9985 was 0.793 once the
    # unfitted targets were counted honestly. Only a fixed-km bandwidth with a compact kernel
    # (bisquare/tricube/boxcar) can reach the branch at all.
    #
    # `:zero` is kept so a run can reproduce the pre-fix numbers exactly.
    unsupported_local_target::Symbol = :missing
end

struct JointFoldContext
    product::String
    variables::Vector{String}
    roles::Dict{String,String}
    spatial_train::Matrix{Float64}
    spatial_target::Matrix{Float64}
    predictor_train::Dict{String,Array{Float64,3}}
    predictor_target::Dict{String,Array{Float64,3}}
    train_lonlat::Matrix{Float64}
    target_lonlat::Matrix{Float64}
    bandwidth_candidates::Vector{Int}
    scaling::DataFrame
    quality_control::DataFrame
    config::JointCovariateBenchmarkConfig
end

# The models, split by step in `joint/`: reading a fixed specification (`joint_spec.jl`), a
# fold's scaled designs (`joint_fold_context.jl`), which columns each model fits locally or
# globally (`joint_designs.jl`), the per-hour fits (`joint_backfit.jl`), and the hour loop
# `dynamic_covariate_predict` (`joint_predict.jl`).
include("joint/joint_spec.jl")
include("joint/joint_fold_context.jl")
include("joint/joint_designs.jl")
include("joint/joint_backfit.jl")
include("joint/joint_predict.jl")

end

"""
Geographically weighted regression primitives: kernels, weights, local solvers, ST-GWR and the
accuracy metrics every evaluation path shares.
"""
module MixedGWR

using LoopVectorization
using ProgressMeter
using LinearAlgebra, Statistics
using Parameters
using Base.Threads

export GWR
export ST_GWR, ST_GWR_fast, ST_GWR_fast!, gwr_neighbors
export fitted, fitted!

export cor
export metric_continuous, metric_event, common_valid_mask, complete_time_mask
export INTENSITY_BOUNDS, INTENSITY_CLASSES, INTENSITY_RANGES, intensity_class

include("core/fitted.jl")
include("core/metrics.jl")
include("core/kernel.jl")
include("core/gw_weight.jl")
include("core/PrecipitationCorrection.jl")
include("core/solve_chol.jl")
include("core/solve_reg.jl")
include("core/GWR.jl")
include("core/ST_GWR.jl")

end # module MixedGWR

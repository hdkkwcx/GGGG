using MixedGWR, RTableTools, Distances, Test
include("main_pkgs.jl");

# These use the Shiyan station fixture, and test-gwmodel.jl also needs R. `data/` is gitignored,
# so the fixture is absent from a fresh checkout and from CI; they are skipped rather than
# aborting the suite. Everything below them is self-contained and always runs.
if HAS_SHIYAN_DATA
  include("test-gwmodel.jl")
  include("test-solver.jl")
  include("test-ST_GWR.jl")
else
  @warn "Skipping the Shiyan fixture testsets: $SHIYAN_DATA_PATH is missing"
end
include("test-gw-weight.jl")
include("test-study-area.jl")
include("test-mger-five-kernels.jl")
include("test-interpolation-benchmark.jl")
include("test-satellite-fusion.jl")
include("test-appeears-ndvi.jl")
include("test-era5-land-stations.jl")
include("test-era5-land-processing.jl")
include("test-terrain-features.jl")
include("test-era5-land-covariates.jl")
include("test-mod13a2-ndvi-processing.jl")
include("test-dem-terrain-experiment.jl")
include("test-era5-variable-selection.jl")
include("test-ndvi-variable-selection.jl")
include("test-joint-variable-selection.jl")
include("test-joint-covariate-models.jl")
include("test-benchmark-diagnostics.jl")
include("test-heavy-rain-events.jl")
include("test-heavy-core-diagnostics.jl")
include("test-landform-classification.jl")
include("test-satellite-temporal-evaluation.jl")
include("test-no-rain-evaluation.jl")
include("test-grid-support-diagnostics.jl")
include("test-run-provenance.jl")
include("test-imerg-uncal-stations.jl")

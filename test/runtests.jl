# test-interpolation-benchmark.jl is about half of the suite, so it runs in its own Julia process
# while the other files run here one after another. Compiling the code under test dominates the
# suite, and these files share what they compile; one process per file compiled it again for each
# and made CI slower (test step 381 s against 295 s serially, on 4 vCPUs).

# The child inherits this process's flags (Pkg.test's --check-bounds=yes) and environment. Its
# output is printed when it finishes.
const BENCHMARK_FILE = "test-interpolation-benchmark.jl"
benchmark_output = IOBuffer()
benchmark_cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) $(joinpath(@__DIR__, BENCHMARK_FILE))`
benchmark = run(pipeline(benchmark_cmd; stdout=benchmark_output, stderr=benchmark_output); wait=false)

try
  include("test-shiyan-fixture.jl")
  include("test-gw-weight.jl")
  include("test-study-area.jl")
  include("test-mger-five-kernels.jl")
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
finally
  # Wait for the benchmark even when a file above fails, so its output is never lost.
  status = success(benchmark) ? "passed" : "FAILED"
  print(string("\n── $BENCHMARK_FILE $status ──\n", String(take!(benchmark_output))))
end
success(benchmark) || error("$BENCHMARK_FILE failed")

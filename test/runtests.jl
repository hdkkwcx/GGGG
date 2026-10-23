# Runs each test file in its own Julia process, `Sys.CPU_THREADS` at a time (set JULIA_CPU_THREADS
# to change that). Compiling the code under test dominates the suite, and a process compiles on one
# core, so the files run side by side; each loads everything it needs. A file's output is printed
# when it finishes.

# Slowest first: test-interpolation-benchmark.jl alone is about half of the serial suite.
const TEST_FILES = [
  "test-interpolation-benchmark.jl",
  "test-mger-five-kernels.jl",
  "test-shiyan-fixture.jl",
  "test-gw-weight.jl",
  "test-study-area.jl",
  "test-satellite-fusion.jl",
  "test-appeears-ndvi.jl",
  "test-era5-land-stations.jl",
  "test-era5-land-processing.jl",
  "test-terrain-features.jl",
  "test-era5-land-covariates.jl",
  "test-mod13a2-ndvi-processing.jl",
  "test-dem-terrain-experiment.jl",
  "test-era5-variable-selection.jl",
  "test-ndvi-variable-selection.jl",
  "test-joint-variable-selection.jl",
  "test-joint-covariate-models.jl",
  "test-benchmark-diagnostics.jl",
  "test-heavy-rain-events.jl",
  "test-heavy-core-diagnostics.jl",
  "test-landform-classification.jl",
  "test-satellite-temporal-evaluation.jl",
  "test-no-rain-evaluation.jl",
  "test-grid-support-diagnostics.jl",
  "test-run-provenance.jl",
  "test-imerg-uncal-stations.jl",
]

# The child inherits this process's flags (Pkg.test's --check-bounds=yes) and environment.
function run_test_file(file)
  output = IOBuffer()
  cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) $(joinpath(@__DIR__, file))`
  seconds = @elapsed passed = success(pipeline(cmd; stdout=output, stderr=output))
  status = passed ? "passed" : "FAILED"
  print(string("\n── $file $status in $(round(seconds; digits=1)) s ──\n", String(take!(output))))
  return passed
end

passed = asyncmap(run_test_file, TEST_FILES; ntasks=Sys.CPU_THREADS)
failed = TEST_FILES[.!passed]
isempty(failed) || error("$(length(failed)) test files failed: $(join(failed, ", "))")

using MixedGWR, RTableTools, Distances, Test
include("main_pkgs.jl");

# These use the Shiyan station fixture, and test-gwmodel.jl also needs R. `data/` is gitignored,
# so the fixture is absent from a fresh checkout and from CI; they are skipped rather than
# aborting the suite.
if HAS_SHIYAN_DATA
  include("test-gwmodel.jl")
  include("test-solver.jl")
  include("test-ST_GWR.jl")
else
  @warn "Skipping the Shiyan fixture testsets: $SHIYAN_DATA_PATH is missing"
end

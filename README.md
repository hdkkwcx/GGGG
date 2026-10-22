# Mixed Geographically Weighted Regression (MGWR) in Julia

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jl-pkgs.github.io/MGWR.jl/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jl-pkgs.github.io/MGWR.jl/dev)
[![CI](https://github.com/jl-pkgs/MGWR.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/jl-pkgs/MGWR.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/jl-pkgs/MGWR.jl/branch/master/graph/badge.svg)](https://app.codecov.io/gh/jl-pkgs/MGWR.jl/tree/master)

> Translated from R package GWmodel

## Examples

- Spatial: the `"GWR"` testset in [test/runtests.jl](test/runtests.jl)

- SpatioTemporal: [test/test-ST_GWR.jl](test/test-ST_GWR.jl)

Both read the Shiyan fixture `data/prcp_st174_shiyan.csv` and compare against R's `GWmodel`.
The original `example/` scripts were removed in commit `3b06d88`; recover one with
`git show 3b06d88^:example/ex01_spatial.jl`.

## References

1. <https://github.com/GWmodel-Lab/GWmodel>

2. <https://github.com/rpkgs/GWmodel>

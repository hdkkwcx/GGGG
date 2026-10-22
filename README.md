# Mixed Geographically Weighted Regression (MGWR) in Julia

[![CI](https://github.com/hdkkwcx/GGGG/actions/workflows/CI.yml/badge.svg)](https://github.com/hdkkwcx/GGGG/actions/workflows/CI.yml)

> Translated from R package GWmodel

## Examples

- Spatial: the `"GWR"` testset in [test/runtests.jl](test/runtests.jl)

- SpatioTemporal: [test/test-ST_GWR.jl](test/test-ST_GWR.jl)

Both read the Shiyan fixture `data/prcp_st174_shiyan.csv` and compare against R's `GWmodel`.
The original `example/` scripts were removed in commit `3b06d88`; recover one with
`git show 3b06d88^:example/ex01_spatial.jl`.

## Project layout

The package is `MixedGWR`; the research project around it corrects FY4B, GPM and GSMaP satellite
precipitation against rain gauges.

| Folder | Holds |
|---|---|
| `src/` | Reusable code: the GWR library (`core/`: kernels, distances, local solvers, mixed/multiscale GWR, metrics), data sources, variable selection, models, evaluation (`benchmark/`), the MGER pipeline and shared utilities (`util/`). See `CLAUDE.md` for how each tier is loaded. |
| `scripts/` | Entry points, grouped by prefix: `download_*` / `gee_export_*` fetch raw data, `prepare_*` build `data/processed/`, `run_*` are the analyses, `verify_*` regression gates, `plot_*` the figures, and `make_paper.jl` rebuilds the paper. |
| `test/` | `julia --project=. -e 'using Pkg; Pkg.test()'` |
| `data/` | Research inputs (gitignored): `raw/`, `processed/`, and a few older files at the top level. |
| `output/` | Every program result, one folder per run, with its provenance (gitignored). |
| `paper/` | The manuscript `paper_plan.typ`, its figures and cited numbers, and `sources.toml`, which says where each comes from. |
| `Interim_results/` | Interim Typst reports for the advisor; the plot scripts refresh their figures and tables. |
| `reference_code/` | Git submodule with the earlier code this project grew from; kept for reference, never loaded. |
| `Projects_Pending_Handover/` | Handover notes from finished pieces of work. |
| `assets/` | Static images. |

## Reproduce the paper

```sh
julia --project=. scripts/make_paper.jl                # what is reused, what would rerun (changes nothing)
julia --project=. scripts/make_paper.jl --only fig06b  # one figure or table and the steps behind it
julia --project=. scripts/make_paper.jl --run          # on a clean tree: rerun stale steps, refresh paper/
```

`paper/sources.toml` is the record: every figure and cited number names the step that makes it,
and every step its exact command, its inputs and what it needs. The dry run lists all of them,
including the data steps (which it never runs) and their inputs. In short:

| Paper | Step chain |
|---|---|
| Fig. 1 | `plot_study_area_map.py` (from `data/dem_ShiYan_1km.tif`, see `scripts/download_dem_shiyan.jl`) |
| Figs. 2-3 | `run_heavy_rain_event_evaluation.jl --gpm-gsmap` → `plot_heavy_rain_events.py` |
| Fig. 4 | `run_satellite_temporal_evaluation.jl` → `plot_satellite_temporal_evaluation.py` |
| Figs. 6-7, Table 3 | `run_interpolation_benchmark.jl full …` → `run_claim_reassessment.jl` → `run_gauge_satellite_fusion.jl` → `plot_gauge_satellite_fusion.py --gpm-gsmap` |

Figure numbers are still a draft. A figure's number lives in its `paper/figures/` file name and its
`id` in `sources.toml`, not in the script names. `python scripts/check_paper_provenance.py` checks
that everything in the paper traces back to committed code.

A result made before its code changed is reported as CHANGED and never rerun on its own. After
changing analysis code, rerun the affected steps with `--run --force <step>`, or, if the change
cannot move the result, record that review as the step's `reviewed` commit in `sources.toml`.

## References

1. <https://github.com/GWmodel-Lab/GWmodel>

2. <https://github.com/rpkgs/GWmodel>

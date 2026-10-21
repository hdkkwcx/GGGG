# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Description

Satellite precipitation correction in Julia: GWR / MGWR (geographically weighted regression,
mixed/multiscale GWR) applied to FY4B, GPM and GSMaP satellite precipitation against rain gauges.

The Julia package is named `MixedGWR` (see `Project.toml`) — `using MixedGWR`, not the repo
directory name.

## Data safety

`data/` and `output/` are gitignored, so git commands that reach ignored files reach the dataset.
**Never run `git clean -x` or `git stash --all`**: `git clean -xfd` deletes all of `data/`
outright, and `git stash --all` strips it off disk. Use `git clean -fd` (no `-x`) or
`git stash -u` instead; `git clean -xdn` is a dry run and is fine.

`.claude/settings.json` denies `Edit`/`Write` into `data/` — results go to `output/`.

## Commands

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # resolve dependencies
julia --project=. -e 'using Pkg; Pkg.test()'          # full suite (test/runtests.jl)
julia --project=. test/test-study-area.jl             # one file; most are self-contained
julia --project=. scripts/run_mger_smoke_202206.jl    # scripts/run_*.jl, verify_*.jl are the entry points
```

`test/test-speed.jl` is a standalone benchmark, not in `runtests.jl`.

The first `@testset "GWR"` in `runtests.jl`, `test-solver.jl` and `test-ST_GWR.jl` check results
against R's `GWmodel` via `RCall`, using `data/prcp_st174_shiyan.csv`. They need R with `GWmodel`
installed, and they are skipped with a warning when that fixture is absent (`HAS_SHIYAN_DATA` in
`test/main_pkgs.jl`), e.g. on a fresh checkout or CI. The other test files don't need R.

### Threads and reproducibility

- Run the interpolation benchmark with `-t 4`, **not** `-t auto` and not single-threaded:
  `julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates`.
  Each hour-fit task allocates ~1.6 MB, and past a few threads that allocation contention makes
  the loop slower (on the 24-core box, 4 threads beat 8 and 24; single-threaded is ~3x slower).
  The output does not depend on the thread count. On new hardware, re-derive the best count with
  `scripts/profile_hour_fit.jl`.
- **Never call `BLAS.set_num_threads`.** `tps` factorizes a dense (n+3)×(n+3) system, and threaded
  LU changes its RMSE in the last digits, which propagates into `metrics_pooled.csv`,
  `paired_comparisons.csv` and `claim_assessment.csv`. The BLAS thread count is part of what makes
  published numbers reproducible.
- If `ci_high` in `paired_comparisons.csv` moves only in its last digit or two, and nothing else
  changes, that is benign floating-point non-associativity, not a regression. See the comment on
  `_daily_bootstrap_delta` in `src/benchmark/InterpolationBenchmarkBootstrap.jl`.

## Architecture

### `src/` layout: the directory decides how a file is loaded

```
src/
  MixedGWR.jl      package entry point (path fixed by Project.toml)
  load_modules.jl  the loader every script and test goes through
  core/            files `include`d into `module MixedGWR` — regression/kernel primitives
  sources/         one standalone module per data source or ingest stage
  selection/       standalone modules: the variable-selection searches
  models/          standalone modules: the estimators, plus CovariateGroups (shared vocabulary)
  benchmark/       the InterpolationBenchmark fragments + standalone evaluation/diagnostic modules
  mger/            the MGERPipeline fragments + MGERDataPrep
  util/            TableIO, RunProvenance — domain-free plumbing shared across tiers
```

Three loading disciplines:

1. **`core/`** files define no module of their own. `src/MixedGWR.jl` `include`s them, and they
   are reached via `using MixedGWR`. Never `include` them from anywhere else.
2. **Standalone modules** (every other `PascalCase.jl` in `sources/`, `selection/`, `models/`,
   `benchmark/`, `mger/`, `util/`) each define their own `module X`, and are loaded only through
   `load_standalone_modules`.
3. **`MGERPipeline.jl` and `InterpolationBenchmark.jl` are not modules.** They are top-level
   fragments evaluated into `Main` by `load_pipeline`, after `using MixedGWR`. Each one includes
   its own fragment files, which share one namespace, so include order only needs to put shared
   consts and structs first.

Lowercase subdirectories inside a tier (`benchmark/diagnostics/`, `models/dem/`,
`mger/pipeline/`, `selection/era5/`, …) hold fragments `include`d into their parent module or
pipeline. Never load them on their own. `MGERConfig` stays in `MGERPipeline.jl` because it is
`load_pipeline`'s sentinel.

Adding a file: a reusable regression/kernel primitive goes in `core/`, plus an `include` line in
`src/MixedGWR.jl`. Anything else becomes a standalone module in the matching directory, and the
loader finds it by name with no registration.

Do not create `src/data/` or `src/docs/`. `.gitignore` has bare `data/` and `docs` patterns that
match at any depth, so files there would be silently untracked.

### Loading `src/` from a script or test

There is exactly one idiom:

```julia
include(joinpath(ROOT, "src", "load_modules.jl"))
load_pipeline("InterpolationBenchmark")            # MGERPipeline / InterpolationBenchmark
load_standalone_modules("BenchmarkDiagnostics")    # any module; dependencies load first
using Main.BenchmarkDiagnostics
```

- Never `include` a standalone module file directly. A second `include` compiles a second,
  type-incompatible copy of the module.
- Do **not** add `pushfirst!(LOAD_PATH, joinpath(ROOT, "src"))`. It makes `using MixedGWR` load a
  second copy of the package, which is a hard error on Julia 1.11.
- `locate_src_file` resolves a file by basename across all of `src/` and refuses when a name
  matches zero files or two. So every basename, fragments included, must be unique. Watch for
  leftover copies after moving a file.
- Inside a standalone module, reach a sibling with `using Main.X`. If you also need qualified
  calls, add `using Main: X`, because `using Main.X: a, b` does not bind `X`.

### Interpolation benchmark invariants

- Opt-in flags (`--satellite-wet-blend`, `--blend-axis`, `--fused-anchor`,
  `--fused-anchor-lagnbr`, `--stacked-blend`) add methods or products without moving any existing
  metric row; the end-to-end test asserts this. Blended methods are scored on the shared mask but
  never join `MASK_METHODS`. The flags themselves are documented where they are parsed in
  `scripts/run_interpolation_benchmark.jl`.
- By default the CLI scores uncalibrated IMERG (`GPM_uncal`, output suffix `_gpmuncal`).
  `--gpm-cal` scores calibrated `GPM` instead. `benchmark_config`'s `gpm_uncal` keyword still
  defaults to `false`, so the verify scripts reproduce the calibrated baselines.
- By default the CLI also leaves FY4B out (output suffix `_nofy4b`, before `_gpmuncal`). That
  removes it from the merges and the agreement-envelope bands, and widens the common grid from
  FY4B's 13,471 hours to all 26,304 hours of 2022–2024. `--with-fy4b` restores it.
  `benchmark_config`'s `fy4b` keyword still defaults to `true`, for the verify scripts.
- `hurdle_gwr` (`InterpolationBenchmarkHurdle.jl`) is deliberately disabled and absent from
  `BENCHMARK_RUNS`.

## Working rules

- `scripts/` handles the process flow and calls into `src/`. Reusable algorithms live in `src/`.
  Do not put one-off experiment workflows in `src/`.
- `src/` takes data and configuration as parameters and never reads fixed local absolute paths.
- Directories: `data/raw/` holds raw data (do not modify), `data/processed/` holds cleaned or
  matched data, `output/` holds all program results (never write results back to `data/`), and
  `assets/` holds images and diagrams. Never create scratch scripts, data or images in the repo
  root.
- Only modify code relevant to the current task. Avoid unrelated refactoring, keep existing
  function names and parameters unless the task needs a change, and never delete existing
  functionality. If you notice other problems, point them out instead of fixing them.
- Prefer simple, concise implementations that match the surrounding code.
- A script whose outputs reach the paper records its provenance: `start_provenance` first thing in
  `main` (git state, Julia version, `Manifest.toml`, a SHA-256 for every input file) and
  `finish_provenance!` into `provenance.csv` at the end (`src/util/RunProvenance.jl`). Its plot
  script calls `write_plot_provenance` (`scripts/plot_provenance.py`). A figure copied into
  `paper/figures/` gets an entry in `paper/sources.toml`, and
  `python scripts/check_paper_provenance.py` checks the whole paper. Runs for the paper start from
  a clean tree. The paper never types a result: a cited number gets a `[[table]]` block in
  `paper/sources.toml`, `python scripts/export_paper_numbers.py` writes it to `paper/numbers/`, and
  `paper_plan.typ` reads it with `val` / `pct`. A number worked out from several rows (a range, a
  median, a largest gap) is a `[table.stats.<name>]` entry, read with `stat`, so the check covers it;
  the paper computes nothing itself.
- The RQ1/RQ2 scripts (`run_heavy_rain_event_evaluation.jl`, `run_satellite_temporal_evaluation.jl`,
  `run_heavy_core_diagnostics.jl`) follow the benchmark's defaults: IMERG Uncal (`GPM_uncal`, suffix
  `_gpmuncal`; `--gpm-cal` for IMERG Final), and RQ2 leaves FY4B out (suffix `_nofy4b`;
  `--with-fy4b` restores its `all_products` sample). `--gpm-cal --with-fy4b` reproduces RQ2's
  original run.

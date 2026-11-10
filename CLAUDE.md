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
julia --project=. test/test-study-area.jl             # one file; each is self-contained
julia --project=. scripts/run_mger_smoke_202206.jl    # scripts/run_*.jl, verify_*.jl are the entry points
julia --project=. scripts/make_paper.jl               # paper status; --run rebuilds what is stale (see its header)
```

One-time setup for the paper workflow:

```sh
python3.13 -m venv ~/.venvs/mixedgwr-plot && ~/.venvs/mixedgwr-plot/bin/pip install -r scripts/requirements-plot.txt
export MIXEDGWR_PYTHON=~/.venvs/mixedgwr-plot/bin/python  # make_paper.jl draws the figures with it (in ~/.zshrc)
ln -s ../../scripts/hooks/pre-commit .git/hooks/pre-commit # refuse a paper/ commit that fails the provenance check
```

That venv redraws the paper's figures byte for byte (checked on Fig. 4). `/paper-check` and
`/merge-branch` (`.claude/skills/`) are the read-only paper check and the repo's merge-and-push
routine.

`runtests.jl` runs `test-interpolation-benchmark.jl` in its own Julia process and `include`s the
other files one after another alongside it. A new test file must still load everything it uses,
and gets an `include` line there. `test/test-speed.jl` is a standalone benchmark, not in
`runtests.jl`.

`test-gwmodel.jl`, `test-solver.jl` and `test-ST_GWR.jl` use `data/prcp_st174_shiyan.csv`. They
run together from `test-shiyan-fixture.jl`, and are skipped with a warning when that fixture is
absent (`HAS_SHIYAN_DATA` in `test/main_pkgs.jl`), e.g. on a fresh checkout or CI. Only `test-gwmodel.jl` needs R: it checks
`GWR` against R's `GWmodel` via `RCall`, so it needs R with `GWmodel` installed. CI installs no R.

CI (`.github/workflows/CI.yml`) runs only when started by hand: `gh workflow run CI`, then
`gh run watch`. A push or a PR does not trigger it, so run the tests locally before committing.

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
  core/            files `include`d into `module MixedGWR` — kernels, weights, station distances,
                   local solvers, mixed/multiscale GWR, metrics
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

Adding a file: a reusable GWR primitive (anything two modules would otherwise each copy) goes in
`core/`, plus an `include` line in `src/MixedGWR.jl`. Anything else becomes a standalone module in
the matching directory, and the loader finds it by name with no registration.

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

- Opt-in flags (`--satellite-wet-blend`, `--blend-axis`, `--blend-dry-band`, `--fused-anchor`,
  `--fused-anchor-lagnbr`, `--stacked-blend`, `--fy4b-qm`) add methods or products without moving
  any existing metric row; the end-to-end test asserts this.
- `--fy4b-file <path>` (with `--with-fy4b`) scores a FY4B series rebuilt by
  `scripts/prepare_fy4b_rebuild.jl` instead of the shipped one. `--fy4b-qm` then adds `FY4B_QM`:
  that rebuild's FY4B, pixel-shifted and quantile-mapped per regime inside each fold on training
  stations only (`SatelliteCalibration`). The rebuild's 1-, 3- and 24-hour raw files are
  gauge-corrected and never read. Blended methods are scored on the shared mask but
  never join `MASK_METHODS`. The flags themselves are documented where they are parsed in
  `scripts/run_interpolation_benchmark.jl`.
- Three options of `FY4B_QM` change rows by design and so get their own suffix, as `--fy4b-file`
  does. `--fy4b-qm-adaptive` (`_fy4bqmadapt`) applies the season-by-regime map, month by month,
  only as strongly as the fold's training gauges support, and changes only `FY4B_QM`'s rows.
  `--fy4b-qm-crossed` (`_qmcrossed`) fits `FY4B_QM` out of time as well: each 7-day time fold's
  hours take the shift, map and strength fitted on the training gauges of the other time folds, a
  day either side left out (`SatelliteCalibration.calibrate_crossed`). It changes only `FY4B_QM`'s
  rows. `--fy4b-qm-fusion`
  (`_qmfusion`) feeds each fold's `FY4B_QM` to the fused anchors and the agreement bands instead
  of the plain FY4B, so `MERGED_*` and every blend row move; the real products' raw and
  unblended rows do not. With `--fy4b-qm-crossed` as well, the fusion fits gauge-hours that never
  shaped the correction at that hour; `scripts/compare_fy4b_fusion_same_hours.jl --crossed`
  compares that run against the paper run.
- A product's position in the run's product list sets its covariate-selection seed (`run_seed` in
  `InterpolationBenchmarkFold.jl`), so adding or removing a product moves the other products'
  selection. Runs are only directly comparable on the same product list. `--fy4b-out-of-fusion`
  (`_fy4bnofusion`, with `--with-fy4b`) keeps FY4B in the list and on the common hours but out of
  the fused anchors and the agreement bands, so it differs from `--fy4b-qm-fusion` in the fusion
  alone; `scripts/compare_fy4b_fusion_same_hours.jl` compares the two.
- `--fy4b-qm-rh-gate` (`_qmrhgate`, with `--fy4b-qm`) appends `FY4B_QM_RH`: each fold's `FY4B_QM`
  set to 0 where ERA5-Land humidity is too low (`HumidityGate`), with the cut cross-fitted over
  gauges and 7-day time folds. Appended last and kept out of the fusion, so no other row moves. It
  needs the joint covariates' ERA5, which the end-to-end test cannot load, so
  `scripts/compare_fy4b_humidity_gate.jl` checks that against the paper run cell by cell instead.
- `--fy4b-confirm-fusion` (`_fy4bconfirm`, with `--with-fy4b`) feeds the fused anchors and the
  agreement bands a FY4B kept only where GPM or GSMaP also calls the cell wet, and 0 where FY4B is
  wet alone (`SatelliteFusion.confirmed_by_others`). It reads no gauge and fits nothing, so only
  `MERGED_*` and the agreement blends move. `scripts/compare_fy4b_fusion_same_hours.jl --confirmed`
  judges it against the paper run by the last FY4B try's stricter rule. Every FY4B try is listed in
  `Projects_Pending_Handover/fy4b_tries.md`.
- By default the CLI scores uncalibrated IMERG (`GPM_uncal`, output suffix `_gpmuncal`).
  `--gpm-cal` scores calibrated `GPM` instead. `benchmark_config`'s `gpm_uncal` keyword still
  defaults to `false`, so the verify scripts reproduce the calibrated baselines.
- By default the CLI also leaves FY4B out (output suffix `_nofy4b`, before `_gpmuncal`). That
  removes it from the merges and the agreement-envelope bands, and widens the common grid from
  FY4B's 13,471 hours to all 26,304 hours of 2022–2024. `--with-fy4b` restores it.
  `benchmark_config`'s `fy4b` keyword still defaults to `true`, for the verify scripts.
- The paper's main comparison has three products on one matched sample: its `benchmark` step runs
  `--with-fy4b --fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv --fy4b-qm --fy4b-qm-adaptive`
  (plain FY4B in the merges and bands, FY4B_QM appended last). The run still scores each product on
  its own mask; `--matched-products` on `run_claim_reassessment.jl`, `run_gauge_satellite_fusion.jl`,
  `compare_fy4b_quality.jl`, `check_fy4b_qm_time_holdout.jl` and `screen_fy4b_clt.jl` rescores on
  the cells every product of the run is scored on (`BenchmarkDiagnostics.matched_evaluation_mask`,
  folder suffix `_matched`) and stops if a method failed a whole fold. Claims and fusion read a
  rebuilt-FY4B run only with the same `--fy4b-file`; the three FY4B scripts take the run folder.
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
  the paper computes nothing itself. Every `[[item]]` and `[[table]]` names a `step` in the
  `[[step]]` graph of `sources.toml`, which `scripts/make_paper.jl` walks and the check enforces: a
  new figure or table, or a new run behind one, gets its step there with the exact command. After
  changing analysis code, `make_paper.jl` reports the steps it touches as CHANGED: rerun them with
  `--run --force <step>`, or, only once the diff is checked and cannot move the result, bump the
  step's `reviewed` commit.
- The RQ1/RQ2 scripts (`run_heavy_rain_event_evaluation.jl`, `run_satellite_temporal_evaluation.jl`,
  `run_heavy_core_diagnostics.jl`) follow the benchmark's defaults: IMERG Uncal (`GPM_uncal`, suffix
  `_gpmuncal`; `--gpm-cal` for IMERG Final), and RQ2 leaves FY4B out (suffix `_nofy4b`;
  `--with-fy4b` restores its `all_products` sample). `--gpm-cal --with-fy4b` reproduces RQ2's
  original run. As in the benchmark, `--fy4b-file <path>` (with `--with-fy4b` in RQ2) scores a rebuilt
  FY4B series instead of the shipped one (suffix `_fy4bw3h0` for `fy4b_w3_h0.csv`); in RQ1 it selects
  the events on the hours the gauges and all three products share, and the core diagnostics then work
  on those windows. `--from YYYY-MM` starts the record at that month (suffix `_from<yyyymm>`, before
  `_block<N>d`). The paper runs all three, and `run_heavy_rain_event_fy4b.jl`, with
  `--fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv --from 2022-06` (RQ2 also `--with-fy4b
  --block-days 7`) and cites RQ2's `all_products` sample only.

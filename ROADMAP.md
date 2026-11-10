# Roadmap

Ideas that are planned but not built. What is already done for FY4B is in
`Projects_Pending_Handover/fy4b_quality_handover.md`.

## FY4B: a substantial gain in hourly correlation

**Goal.** Raise FY4B's hourly correlation with the gauges. On the 20,087-hour grid, held-out
gauges, balanced spatial folds:

| product | hourly r |
|---|---|
| FY4B (rebuilt, 3×3 mean) | 0.223 |
| FY4B_QM (per-fold shift + quantile map, monthly strength) | 0.246 (0.220 with weeks held out too) |
| GSMaP | 0.379 |
| IMERG Uncal | 0.436 |

**Status (2026-11-01): built as counting tables, screened, and stopped there by the user's
decision; `FY4B_QM` stays the FY4B quality step. See the handover's §0.**
- **Built:**
  - `src/models/SatelliteTwoStage.jl`, without machine learning at the user's request: per fold,
    on training gauges, it counts P(wet) and the wet amount per context class. The classes are
    intensity, rain nearby, and rain at t±1.
  - `scripts/screen_fy4b_two_stage.jl` scores it on held-out raw FY4B.
- **Results:**
  - P × amount raises r from 0.247 to 0.290 and cuts RMSE by 13%, but it never reaches 8 mm/h.
  - Most of that is not context. An amount-only table (the intensity class alone) already gives
    +0.032 of the +0.043 in r. The context adds +0.011, mainly by sorting light rain (CSI at
    0.1 mm/h 0.122 → 0.181).
  - The occurrence rule with the QM amounts changes nothing.
  - **These numbers held out gauges, not weeks.** Rerun with weeks held out as well (2026-11-10,
    `screen_fy4b_two_stage.jl --crossed`, as FY4B_QM's 11-07 check), the r gain does not survive:
    - TS_E against FY4B_QM goes from +0.043 to +0.022 (−0.013..+0.056);
    - against plain FY4B it is +0.019 (−0.019..+0.056);
    - TS_E's own r falls 0.290 → 0.243.
    - The RMSE gain stays (−26.5% against plain FY4B), but it is shrinkage toward a conditional mean, with heavy
      rain still lost (CSI_8 0).
    - The context adds +0.009, still small.
- **Not built:**
  - the 15-minute quarters, the product-quality features, hour of day and elevation;
  - a benchmark run, to score MGWR and the blend downstream.

**Why the current fixes stop short.**
- The pixel shift fixes where FY4B puts rain, and the quantile map fixes how much.
- A quantile map is monotone within each season-by-regime cell, so it barely changes which hours
  rank as wetter. Correlation therefore moves little.
- A larger gain needs a model that can re-rank hours: for example, one that recognises FY4B's
  false rain and its missed rain from context.

**Idea: a two-stage model.** First decide whether it rains, then estimate how much.
1. **Occurrence:** P(gauge ≥ 0.1 mm/h | features), e.g. logistic regression or boosted trees.
2. **Amount:** E(gauge | rain, features), fitted on wet training hours only.
3. **Prediction:** P × amount, or the amount where P passes a threshold chosen on training gauges.

**Leakage rule.** As with `FY4B_QM`, both stages are fitted per fold on that fold's training gauges
only, inside the benchmark (`src/models/SatelliteCalibration.jl` is the pattern).

**Features: FY4B only**, so this stays a FY4B quality step rather than a fusion with GPM or GSMaP:
- **From the offset grid** (`fy4b_offset_grid.nc`): the shifted 3×3 mean, plus the 5×5 max and
  wet-pixel fraction around the station.
- **The four 15-minute values inside the hour:** number of wet quarters, maximum, and trend. This
  needs a rebuild variant in `scripts/prepare_fy4b_rebuild.jl` that keeps quarter statistics; today
  the quarters are summed.
- **Product quality:** the 20 mm/h cap fraction and the share of DQF-flagged pixels.
- **Time context:** the previous and next hour (persistence), hour of day, season and regime
  (satellite position and software version).
- **Station:** elevation. It is gauge-free, so it is safe to use everywhere.

**Check first.** `hurdle_gwr` (`src/benchmark/InterpolationBenchmarkHurdle.jl`) is an existing
two-part model. It is kept but deliberately absent from the benchmark, and no reason is recorded
(`BENCHMARK_RUNS` in `InterpolationBenchmarkConfig.jl`). It corrects with GWR rather than rebuilding
FY4B, but its occurrence and amount pieces may be reusable. (Checked 2026-11-01: its occurrence
part is a wet frequency per satellite-intensity class and month, with no context to re-rank by,
so nothing of it was reused.)

**How to judge it.**
- r, as the target metric;
- CSI at 0.1, 2 and 8 mm/h;
- the frequency bias index, since P × amount compresses the distribution and can look good on r
  and RMSE while losing heavy rain;
- downstream: MGWR and the agreement blend, as for `FY4B_QM`.

No gain is promised: FY4B's own information, an infrared-based estimate, sets the ceiling.

## FY4B: confirmed FY4B in the merges only (exploratory)

**Status (2026-11-09): not built.** The user parked it as an exploratory approach while try 15 ran.
Try 15 put confirmed FY4B in the merges *and* the agreement bands (`--fy4b-confirm-fusion`; register
in `Projects_Pending_Handover/fy4b_tries.md`).

**Update (2026-11-10): try 15 failed, so this stays parked unless asked.** This variant would share try
15's merged anchor. That anchor's raw gain (MERGED_MEAN −0.97%) is gone after MGWR (+0.17%,
−0.04..+0.39), before any band is involved, so this variant is unlikely to clear the bar. It would keep
FY4B's own agreement blend whole, though: try 15 cost that blend 2.3%.

**What it is.** Confirmed FY4B is FY4B kept only where GPM_uncal or GSMaP also calls the cell wet, and 0
where FY4B is wet alone (`SatelliteFusion.confirmed_by_others`). Here it would go into MERGED_MEAN and
MERGED_OLS only. The agreement bands would keep plain FY4B.

**Why it might be worth a run.**
- **A cleaner test.** With the bands unchanged, the GPM_uncal, GSMaP, FY4B and FY4B_QM agreement blends
  are identical by construction, so they become guards. Only the merged rows can move.
- **FY4B's own Table 3 row stays safe.** In try 15's bands, the cells where FY4B is wet alone go to band
  0, so FY4B's and FY4B_QM's own agreement blends no longer blend them.
- **It isolates the anchor.** It tests the merged anchor's gain without the bands' share of it.

**How it would be built.**
- A scope option on the same flag, so that the bands keep the plain sources while the fused anchors
  take the confirmed one. `_fusion_source_matrices` in `InterpolationBenchmarkRun.jl` is where the two
  part.
- A `compare_fy4b_fusion_same_hours.jl` mode with those guards.

**How to judge it.**
- It would be try 16, so the primary's interval is 1 − 0.05/16, under try 15's rule otherwise.
- It is not independent of try 15: it uses the same series, and its merged anchors are the same. So
  report it as a follow-up of try 15, never as a fresh test.
- Run it only if try 15's result leaves a reason to: for example, the merges gain but FY4B's own blend
  pays for it.

## Later

- **Paper write-up of the FY4B work:**
  - a comparison script with provenance, replacing the scratch `output/fy4b_partB_analysis/paired.jl`;
  - a `[[step]]` plus `[[table]]` blocks in `paper/sources.toml`;
  - RQ2 on the rebuilt FY4B, which needs a `--fy4b-file` option in
    `run_satellite_temporal_evaluation.jl`;
  - a rewrite of the scope note at `paper/paper_plan.typ:174`.
- **Data gaps:** the official site has nothing more for 2022 Jan–May, Feb 2024, the daily 14:15 UTC
  file and about 2,400 scattered files (`output/fy4b_download_latest/`). Revisit only if a new
  source appears.

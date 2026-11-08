# Handover: FY4B quality work (raw data coming via the Samsung T9)

**Date:** 2026-11-08 · **State:** the FY4B work is merged (`805bad6`); two-stage counting tables
screened on raw FY4B and stopped there by the user's decision (§0), so `FY4B_QM` stays the FY4B
quality step. A blend screen and a nested run of two dry-band blends found nothing that clears
the decision rule (merged as `f720404`). Adding FY4B_QM to the fusion clears it for no product's
agreement blend (merged as `24e392e`). A whole-hour time shift does not help, 105°E is not shown
to be better, and the 2023 software-update gain is not shown (merged as `1dc0262` and
`3796bdb`). The cloud filter from FY4B's Cloud Type product (`CLT`) fails its pilot on July 2023
(merged as `1a61f2a`). The paper now has FY4B as "raising its quality" (Previous, merged as
`0982a32`), and its Scope scores the same rebuilt FY4B series as §2.5 and §3.4. A same-hours,
same-seeds test of FY4B in the fusion gives "keep FY4B out" (merged as `b5499c7`). An error budget shows
where FY4B's error sits, and a near-surface humidity gate passes its screen with a small gain. With time
held out as well as gauges, FY4B_QM's raw gain shrinks materially (Previous). Six AGRI L1 scans give a weak
"no" to the channel pilot. All of it is merged (`be63c2e`), and §3.4 now carries the time-held-out figures. The
paper's day-block intervals, rechecked with 7- and 14-day blocks, keep every verdict the text states, and the
paper now uses 7-day blocks (Previous, merged as `3d41498`). The humidity gate's benchmark run says stop: the gate
does not survive MGWR and the blend (Newest, branch `fy4b/humidity-gate-benchmark`, unmerged) · **Next step:** the
user decides whether to merge that branch; no FY4B quality step is pending

---

## Newest (2026-11-08): the humidity gate in the benchmark: stop

The user asked for the benchmark run of the humidity gate that passed its screen. ERA5-Land is reanalysis, so a kept
gate would have been "FY4B + reanalysis". No paper text or number changes.
- **Code** (branch `fy4b/humidity-gate-benchmark`):
  - `src/models/HumidityGate.jl` (`24f63ca`): the screen's gate B8, with the cut cross-fitted over gauge folds and
    7-day time folds;
  - `--fy4b-qm-rh-gate` (`d7a8cc7`, suffix `_qmrhgate`): appends `FY4B_QM_RH`, kept out of the fusion. It writes
    `fy4b_qm_rh_gate_cuts.csv`;
  - `scripts/compare_fy4b_humidity_gate.jl` (`423b2f4`), committed with its rule before the run;
  - review bump `682655a`.
- **Run:** clean at `682655a`, `-t 4`, 3 h 52 min, 12 GB. Log: `output/full_benchmark_fy4bw3h0_fy4bqm_fy4bqmadapt_qmfusion_qmrhgate.log`.
  Output: `..._qmfusion_qmrhgate_gpmuncal/`; comparison in `output/fy4b_humidity_gate_comparison/`.
- **Guards passed:**
  - every base product's `raw`, `adw`, `mgwr`, both blends and mask are identical, cell for cell, to the paper's
    `fy4b_benchmark` run;
  - `FY4B_QM_RH` raw differs from FY4B_QM only where the gate zeroed a value under 8 mm/h.
- **What the gate did:** it zeroed 26% of FY4B_QM's wet cells. The cut was RH < 70 in 37 fold × time-fold fits and
  RH < 75 in 13, like the screen.

### Result: the gain at raw does not survive MGWR or the blend

FY4B_QM_RH against FY4B_QM, balanced spatial, 95% intervals over 7-day blocks:

| Field | ΔRMSE | ΔCSI_0.1 | ΔCSI_8 |
|---|---|---|---|
| raw | −1.73% (−2.81..−0.84) | +0.014 | 0 |
| MGWR | +0.25% (−0.20..+0.71) | +0.024 | −0.006 (−0.010..−0.003) |
| ADW blend | +0.31% (−0.03..+0.76) | +0.021 | −0.001 |
| agreement blend | +0.03% (−0.13..+0.21) | +0.022 | −0.000 |

**Rule (agreement blend): fails.**

| Check | Result |
|---|---|
| ΔRMSE −0.3% or better, interval below 0 | FAIL: +0.03% |
| ΔCSI_8 not significantly lower | pass: −0.0002 (−0.0016..+0.0014) |
| \|FBI_8 − 1\| not significantly worse | FAIL: +0.005 (+0.000..+0.010) |
| \|bias\| up by at most 0.005 mm/h | pass: +0.004 |

- **The gap to GPM_uncal's agreement blend does not move:** +2.01% with the gate, +1.98% without.
- **Random split** (reported only): every blend is within ±0.3%.
- **Reading:** the gate does sort light rain (CSI at 0.1 mm/h +0.02 at every stage), but its RMSE gain is gone after
  MGWR. That fits the screen's warning: `d2m_c` is already a local MGWR covariate in most folds, so the residual model
  already knows the air is dry. At MGWR alone it even loses some 8 mm/h detection.
- **Dew point and the gate do the same job** (checked after the merge, in the run's `joint_fold_roles.csv`): MGWR
  picks `d2m_c` in 8 of 10 fold fits for FY4B_QM (3 of 5 balanced, 5 of 5 random) but in 1 of 10 for FY4B_QM_RH
  (1 of 5, 0 of 5). The other covariates hold: `u10` 10 of 10 in both, `sp_hpa` 9 → 7, `v10` 2 → 1. Once the gate
  has removed the dry-air rain, the residual model no longer needs humidity. FY4B_QM_RH comes last in the product
  list, so its selection seed differs from FY4B_QM's (`run_seed`): strong evidence, not a clean test.
- **So:** FY4B_QM stays FY4B's quality step, and no humidity gate goes into the paper. Merging the branch keeps the
  option and the module as the record; that is the user's call.

## Previous (2026-11-07): the paper's intervals with week blocks

A side agent noted that most of the paper's 95% intervals resample single days, but FY4B's errors come in
multi-day episodes. On §3.4's raw −16.9%, week blocks had widened the interval about 40%. No paper number changes.
- **Code** (branch `stats/week-block-intervals`):
  - `TimeBlockFolds.block_times` and `block_days_option`, with tests (`21cfc46`);
  - an opt-in `--block-days N` in `run_claim_reassessment.jl`, `run_gauge_satellite_fusion.jl` and
    `compare_fy4b_quality.jl`. It changes only the times the bootstrap groups by, and writes to `_block{N}d`
    folders beside the published ones;
  - review bump `763ea8d`;
  - `scripts/compare_interval_blocks.jl` (`1ceb0fe`), its reading committed before the runs. `c5a172a` loosened
    one check: `fy4b_quality`'s point estimates differ across block lengths by summation order alone (≤ 1.3e-12
    relative).
- **Runs:** all clean, about 2 min each for N = 1, 7, 14. Logs: `output/interval_block_runs.log`, `output/*_block{N}d.log`,
  `output/compare_interval_blocks.log`. Output: `output/interval_block_check/` (`verdicts.csv`, `summary.csv`).
- **Guard passed:** every `--block-days 1` file is byte-identical to the published one, across
  `paired_comparisons_by_method.csv`, the four fusion tables and `fy4b_quality/paired.csv`.

### Result: every verdict the text states holds; one Fig. 6 cell changes

| Family | Rows | Median widening, 7 d / 14 d | Flips at 7 d / 14 d |
|---|---|---|---|
| RQ3 fusion vs ADW (Holm) | 120 | 1.07× / 1.14× | 4 / 11 |
| RQ3 interpolation ranking | 15 | 0.99× / 1.06× | 0 / 1 |
| RQ3 anchors vs GPM | 15 | 1.15× / 1.09× | 1 / 1 |
| RQ3 detection skill | 480 | 1.09× / 1.11× | 15 / 18 |
| §3.4 FY4B quality | 18 | 1.13× / 1.17× | 2 / 3 |

**Every verdict the text prints survives both block lengths:**
- RQ3: the blend beats ADW on GPM and on MERGED_OLS (Holm p 0.000 at 1, 7 and 14 days), but not on GSMaP (0.090 →
  0.248 → 0.392). On heavy hours, MGWR and the blend beat ADW on GPM and GSMaP (all 0.000).
- §3.4: the intervals keep their verdicts. FY4B_QM against plain FY4B at raw, MGWR and blend still excludes 0; against
  GPM after MGWR and against GSMaP after the blend it still spans 0. The FY4B intervals widen most, by up to 1.39× at
  raw and MGWR, which fits its episodes.

**What does change at 7 days:**
- **Draft Fig. 6, one cell:** MERGED_OLS with MGWR alone, all hours, +2.2%. Holm p goes 0.033 → 0.165 (0.143 at 14 days),
  so the cell turns hollow. The text never cites this cell.
- **Not in the text:**
  - three MERGED_MEAN blends at 2.5–8 mm/h;
  - GSMaP's agreement blend against GPM's in the anchor comparison;
  - 15 of 480 detection-skill rows, all borderline (p 0.01–0.05);
  - in `fy4b_quality`, FY4B_QM against GSMaP raw (−12%), and the random split's blend against plain FY4B.
- **At 14 days only (sensitivity):** four more Fig. 6 cells turn hollow. These are GPM's residual and agreement-residual
  methods on heavy hours, and MERGED_MEAN's mixed and residual methods on all hours.

**The user's call:**
1. Keep the day-block intervals and add one methods sentence: week blocks were checked and change no stated verdict,
   with one Fig. 6 cell turning hollow.
2. Or switch the paper's RQ3 and §3.4 intervals to 7-day blocks. This is cheap: about 2 min per run, a Fig. 6 redraw,
   and re-exported numbers. RQ2's whiskers would then need the same option for consistency.

I recommend 2: the evidence for multi-day dependence is direct, and no text verdict moves.

**Done (2026-11-07, at the user's request), on branch `paper/week-block-intervals`:**
- **Wording:** RQ3's "only fused GPM beats ADW significantly" now reads "of the two [products] only fused GPM…",
  so it no longer contradicts the merged anchor's significant gain. A side agent's note.
- **RQ2 option:** `--block-days` added to RQ2 (`SatelliteTemporalEvaluation.bootstrap_blocks`). Only the two
  bootstraps resample 7-met-day blocks; counts and tables keep days. `--block-days 1` reproduced every RQ2 CSV byte
  for byte.
- **Plot labels:** the fusion and RQ2 plots name the block length from a `_block<N>d` folder.
- **Reruns:** claims, fusion, RQ2, §3.4's FY4B comparison (new step `fy4b_quality_week`) and both plots, clean at
  `a725000`. The paper reads the `_block7d` folders. `fy4b_quality` (day blocks) stays as the time-holdout check's
  input.
- **What moved:** only intervals and p-values.
  - GSMaP's p is 0.25.
  - Draft Fig. 6's MERGED_OLS + MGWR all-hours cell is hollow.
  - §3.4's raw interval is −27.6 to −7.0%.
- **Methods text:** says intervals resample whole 7-day blocks, and that they resample time only, with the gauges
  and five folds fixed (another side-agent note). The one-month cloud-type pilot keeps day blocks, as §3.4 says.
- **Checks:** provenance OK, `make_paper.jl` with nothing flagged, typst clean, `Pkg.test()` passes. The branch is
  unmerged.

## Previous (2026-11-07): FY4B_QM's raw gain depends on having seen the same weeks

A side agent warned that the channel pilot held out gauges but not days. False rain covers the network at
once, so a held-out gauge is not new weather. That holds, and it reaches a number the paper cites. No paper
text or number changes here: what §3.4 says is the user's call.
- **Code:**
  - `src/benchmark/TimeBlockFolds.jl` with `test/test-time-block-folds.jl` (`0839720`): 7-day blocks in 5 time
    folds, balanced within regime × season, ±1 day guard, crossed with the gauge folds. A synthetic test shows
    a table that memorises episodes scores 100% under gauges-only and under 65% crossed, and that a real
    signal survives both.
  - `scripts/check_fy4b_qm_time_holdout.jl` (`27f2278`), its reading committed before the run. Review bump
    `456b9cf`, for the new module only.
  - `scripts/explain_fy4b_qm_time_holdout.jl` (`3e45ae6`), written after that result, with no verdict of its own.
- **Runs:** both clean, 3 min each. Logs: `output/check_fy4b_qm_time_holdout.log`,
  `output/explain_fy4b_qm_time_holdout.log`. Output: `output/fy4b_qm_time_holdout/`,
  `output/fy4b_qm_time_holdout_explained/`.
- **Guards passed:** the gauges-only version equals the run's `oof_raw.csv` in all 4,760,619 cells and
  reproduces `fy4b_quality` (−16.86%, r 0.246 vs 0.223). The pixel shift is not held out in time.

### Result: material, by the reading fixed before the run

Raw FY4B_QM against plain FY4B (RMSE 1.1426, r 0.223), balanced spatial, held-out gauges, whole-week intervals:

| Version | RMSE | Change | r | Change in r |
|---|---|---|---|---|
| Gauges held out (the paper's) | 0.9500 | −16.9% (−27.6..−7.0) | 0.246 | +0.023 (+0.000..+0.043) |
| Gauges and weeks held out | 1.0338 | **−9.5% (−15.4..−3.5)** | 0.220 | **−0.004 (−0.021..+0.015)** |
| Same training hours, scattered at random | 0.9519 | −16.7% (−27.0..−7.2) | 0.248 | +0.025 (+0.002..+0.044) |
| Only the map held out | 1.0092 | −11.7% (−18.0..−5.3) | 0.224 | +0.000 (−0.020..+0.021) |
| Only the monthly strength held out | 0.9566 | −16.3% (−27.0..−6.6) | 0.245 | +0.022 (−0.000..+0.041) |

The RMSE gain moves 7.3 points (limit 2) and the r gain 0.027 (limit 0.01). The random split moves about the
same: −17.3% → −9.2%, r +0.027 → −0.001.

- **Not less data:** with as many training hours, scattered so every weather system stays on both sides, the
  gain is intact (+0.20%, −0.25..+0.72, against gauges only).
- **The map, not the monthly strength.** My plan suspected the strength; holding it out costs 0.6 points. The
  map's top end is what learns the weeks.
- **A few episodes carry most of it.** 8 days hold 82% of the net rise in squared error: 2023-09-18/19,
  2023-03-16/18/19 and 2022-09-22/23/24. On the March 2023 and September 2022 days the crossed version scores about
  like plain FY4B; on the September 2023 days it lies between plain and gauges only. Most of the loss is MAM
  (−26.7% → −9.8%) and SON (−36.2% → −20.3%). JJA hardly moves (−3.4% → −3.0%).
- **How the map does it,** gauge fold 1. These episodes supply much of a season's largest FY4B values, so the
  map learns to shrink those values. FY4B's 16 mm/h maps to:

  | Season, regime | All hours | Without the fold holding the episode |
  |---|---|---|
  | R1 MAM (Mar 2023) | 4.0 | 7.1 |
  | R1 SON (Sep 2022) | 9.0 | 15.8 |
  | R2 SON (Sep 2023) | 5.0 | 10.7 |

- **Not only the gross hours:** leaving them out, the gain still moves from −15.0% to −10.1%, and r's from +0.030
  to +0.007 (−0.012..+0.025).
- **Heavy rain:** POD_8 is 0.135 plain, 0.085 gauges only, 0.099 crossed; CSI_8 0.047, 0.054, 0.048. So CSI_8's
  gain goes too, and the crossed map dries heavy rain a little less. This is consistent with the map's tail
  drying the RQ1 storm days, though that is not tested.
- **No sign of an artefact:** the crossed version sits at or within 0.006 of plain FY4B's RMSE in every month, so
  its sparse-cell extrapolations (a few reach hundreds of mm/h in thin cells) do not drive the result.

### How to read it

**This is not a leak, and the −16.9% is not wrong.** A held-out gauge still never shapes its anchor, so §2.5's
sentence stands. The two numbers answer different questions:
- **−16.9%:** an ungauged site, with gauges elsewhere in the same weeks. That is the benchmark's question.
- **−9.5%, with no gain in r:** a correction that has not seen those weeks. That is nearer what "FY4B's quality
  was raised" suggests.

**Not measured:** MGWR (−5.4%) and the blend (−0.8%). They would need a benchmark run on a crossed anchor,
meaning a new flag and about 3.6 h. FY4B_QM vs GSMaP after MGWR (−6.8%) rests on them too.

**Recommendation (the user's call):** keep §3.4's numbers, and add that the calibration sees the same weeks, with
the time-held-out figures beside them (needs a `[[step]]` and `[[table]]`). Drop any wording that r improved as a
lasting property of the product. Run the crossed benchmark only if the tutor or a reviewer needs the blend number.

**Done (2026-11-07, at the user's request):** branch `paper/fy4b-time-holdout`.
- §3.4 keeps every number. After raw r's rise, it now says the correction is fitted in the same weeks as the
  gauges it is scored at, and gives −9.5% (−15.4 to −3.5%) against −16.9% (−27.6 to −7.0%), both with whole-week
  intervals. It says raw r no longer rises (−0.004, −0.021 to +0.015), that a few false-rain episodes carry the
  difference, and that MGWR and the blend hold for the same weeks only.
- `paper/sources.toml` gains step `fy4b_time_holdout` and table `fy4b_time_holdout`. The header comment names both
  scripts.

**For improving FY4B:** the map is partly absorbing a few false-rain episodes, at a cost to real heavy rain. A way
to flag those episodes would help more than a better map, which is the channel pilot's target.

### What changes elsewhere

- **The humidity gate's verdict stands** (exploratory run, not committed). With gauges and weeks held out it picks
  the same cuts and barely moves (CSI_0.1 0.1890 → 0.1888, RMSE 0.9328 → 0.9322). One cut cannot memorise an
  episode. It was not retested on the crossed FY4B_QM.
- **The channel pilot's design changes, before its rule is committed.** No extra download:
  1. Crossed hold-out with `TimeBlockFolds`, judged on that; gauges-only reported, to measure the optimism.
  2. A table cell needs at least 10 distinct days as well as its pairs: `TS_MIN_PAIRS = 200` can be one hour.
  3. Rule (proposed): at a cut losing ≤10% of hits, the extra channels remove ≥50% of held-out gross false-rain
     squared error, and ≥50% in at least half of the 81 runs of consecutive days, and beat IR1 alone by ≥20 points
     with the difference's interval above 0. Intervals resample whole weeks.
  4. Train on two regimes and score the third, reported only.
  5. The whole-record screen and any later product are cross-fitted in time as well.
- **Why the sample is small:** FY4B's 398 gross hours are 90,125 gauge-hours, but 113 days in 81 runs, and 10 days
  hold 73% of their false-rain squared error.

### Later the same day: a first look at two AGRI L1 scans

The user downloaded the two scans named in Previous to `/Volumes/T9/FY4B_L1`, linked as `data/FY4B_L1`. These
are read-only checks run inline, with nothing committed.

**The files:**
- 107 and 108 MB each. September 2022 was available, so the archive reaches back at least that far.
- HDF5, which `NCDatasets` opens. `Data/NOMChannel01`–`15` hold counts on the 2748 × 2748 nominal grid at 133°E.
  `Calibration/CALChannel01`–`15` are 4096-entry lookups from count to brightness temperature (or reflectance).
- Infrared channels: 9 is 6.25 µm, 10 is 6.95, 11 is 7.42, 12 is 8.55, 13 is 10.8, 14 is 12.0, 15 is 13.3; 7 and 8
  are 3.75 µm.

**Aligned with QPE as `[x, y]`:** of QPE's pixels wet at ≥ 1 mm/h, 0.7% and 0.2% are warmer than 250 K as read,
against 88% and 74% with the axes swapped.

**QPE is a one-channel lookup per region.**
- Within any 5° × 5° box, BT10.8 alone (0.5 K bins) explains R² 0.93–0.995 of QPE; over the whole disk only
  0.51–0.64, because the threshold varies by region.
- Hubei's box (30–35°N, 108–113°E) gives 0.982 and 0.931. At ≤ 220 K, 94–100% of pixels are wet.
- So the earlier description holds, with "per region" added.

**The two worst false-rain hours are thick, cold, high ice cloud, not thin cirrus.** Medians at the gauges, with
the same scan's QPE-wet pixels over the whole disk in brackets:

| Scan, gauges | BT10.8 | BT6.25 − BT10.8 | Split window BT10.8 − BT12.0 |
|---|---|---|---|
| 2023-03-18 23:30 UTC; gauges dry, QPE 16.2 mm/h | 213.7 K (221.3) | −2.0 (−4.8) | 1.18 (1.59) |
| 2022-09-23 07:30 UTC; gauges dry, QPE 8.6 mm/h | 226.1 K (220.1) | −6.4 (−4.1) | 1.53 (1.59) |

- Water vapour minus window near 0 means a high cloud top. A small split window means the cloud is optically
  thick: thin cirrus would show a larger one. So on these channels, these false-rain pixels look at least as
  rain-like as FY4B's typical rain pixels.
- Two scans cannot decide this, but they lower the chance that a table over these channels separates the worst
  episodes.

**Practicality:** the pilot's ~1,100 scans are about 118 GB, which fits on the T9 (786 GB free). The whole record
at one scan per hour is about 2.1 TB, which does not; it would need a download, extract, delete loop.

**Proposed next, about 430 MB:** four scans to contrast real and false rain in the same seasons and regimes. All
four QPE scans exist.

| Scan (UTC start) | Label (BJT) | Kind |
|---|---|---|
| `20230918203000` | 2023-09-19 05 | False: FY4B 2.5 mm/h, wet at 55% of gauges, gauges 0.01 |
| `20230919123000` | 2023-09-19 21 | Real, same day: gauges 3.4 mm/h, wet at 80% |
| `20221005213000` | 2022-10-06 06 | Real, R1 SON: gauges 2.1 mm/h, wet at 98% |
| `20230422223000` | 2023-04-23 07 | Real, R1 MAM: gauges 2.6 mm/h, wet at 89% |

- File names follow `FY4B-_AGRI--_N_DISK_1330E_L1-_FDI-_MULT_NOM_<start>_<start + 14:59>_4000M_V0001.HDF`.
- The September 2023 episode's large errors are mostly FY4B far too wet during real rain. At label 2023-09-19 20,
  FY4B shows 13.0 mm/h against the gauges' 3.1, so its pure false-rain hour is milder.
- **Reading, fixed before looking:** if the real-rain scans' cloud tops at the gauges fall inside the false ones'
  10–90% range on BT10.8, water vapour minus window, and the split window, stop the channel pilot. If they clearly
  differ, go to the ~1,100-scan sample.

### Result of the four contrast scans: stop the channel pilot

The user downloaded the four scans (85–103 MB each), and all six now sit in `/Volumes/T9/FY4B_L1`, linked as
`data/FY4B_L1`. They were read inline; nothing is committed. Medians at all gauges (10–90%):

| Scan, label (BJT) | Kind | Gauges | QPE | BT10.8 (K) | WV − IR | Split window |
|---|---|---|---|---|---|---|
| 2023-03-19 08 | false | 0.00, wet 0% | 16.2, wet 100% | 213.6 (211.9..221.9) | −2.0 (−5.4..−1.4) | 1.18 (0.90..1.46) |
| 2022-09-23 16 | false | 0.01, wet 1% | 8.6, wet 100% | 226.1 (218.0..239.2) | −6.4 (−15.1..−1.6) | 1.53 (1.18..2.24) |
| 2023-09-19 05 | false | 0.01, wet 2% | 2.8, wet 51% | 252.6 (222.2..281.6) | −22.5 (−46.4..−5.0) | 1.51 (0.87..2.55) |
| 2023-09-19 21 | real | 3.37, wet 80% | 11.3, wet 93% | 221.8 (206.1..244.8) | −5.4 (−18.4..−2.2) | 2.01 (1.43..3.05) |
| 2022-10-06 06 | real | 2.06, wet 98% | 1.4, wet 73% | 239.9 (228.2..251.9) | −12.5 (−19.5..−6.6) | 1.45 (0.85..2.09) |
| 2023-04-23 07 | real | 2.60, wet 89% | 1.2, wet 79% | 228.7 (223.8..237.6) | −4.7 (−9.4..−3.1) | 1.26 (1.04..1.47) |

- **The reading says stop.** All 9 real-scan medians fall inside the false scans' pooled 10–90% ranges.
- **The pooled range is wide partly from one scan.** The third false scan is partly clear (FY4B wet at 51% of
  gauges), which widens it.
- **The sharper comparison points the same way** (information only). At gauges where QPE is wet, pooled over the six
  scans, there are 515 hits and 635 false alarms. Single features barely separate them, and false alarms look *more*
  rain-like:

  | Feature | Hits | False alarms | AUC |
  |---|---|---|---|
  | BT10.8 | 227.9 K | 223.6 K (colder) | 0.58 |
  | WV − IR | −5.4 | −4.5 (higher top) | 0.43 |
  | Split window | 1.50 | 1.35 | 0.60 |
  | QPE amount | 2.4 | 9.4 | 0.35 |

  FY4B's false rain sits under its coldest, highest, thickest cloud tops. This fits the CLT pilot, where ice cloud
  held 83% of the false alarms and 88% of the hits.
- **Limits:** six scans from six episodes, and only three infrared features. Visible and near-infrared channels
  (daytime only) were not looked at.
  - The 1,150 gauge comparisons are not 1,150 cases. Gauges in one hour sit under the same cloud, so this is
    about six cases.
  - The stop rule was easy to trigger: the partly clear false scan widened the range real rain only had to fall
    inside.
- **So, a weak "no" from a small pilot, not a proven one.** The evidence leans no: false rain looked even more
  rain-like than real rain. No L1 sample download for now, and `FY4B_QM` stays the quality step.
  - Call it a small pilot when telling the tutor.
  - If a firmer answer is wanted, a few dozen more scans (about 100 MB each), spread over episodes and seasons,
    would give one.
  - The budget, the time-holdout check and this pilot suggest that FY4B's own infrared does not tell this false
    rain from real rain. Separating them would need information FY4B does not have, such as microwave, which the
    fusion tests already cover.

## Previous (2026-11-06): a near-surface humidity gate passes its screen, with a small gain

The budget's false rain sits in episodes FY4B alone cannot flag, so this asks whether humidity on disk can.
ERA5-Land is reanalysis, not FY4B: anything kept is "FY4B + reanalysis" and could not run in real time.
No paper number moves.
- **Script:** `scripts/screen_fy4b_humidity.jl` (`3f1c3cc`), committed with its rule before any run, run
  clean, 2 min. Log: `output/screen_fy4b_humidity.log`. Output: `output/fy4b_humidity_screen/`.
- **Setup:** as the CLT screen, on all grid hours. Held-out FY4B_QM equals the run's `oof_raw.csv` in all
  4,760,619 cells of both schemes.
- **Alignment:** read as the benchmark reads its covariates, `load_era5_panel` at `time_utc + 9 h`, so an
  hour takes the value at its start. t2m peaks at label 15 BJT. RH's anomaly best matches the gauge's wet
  call at lag +3 h (r 0.314; 0.296 at lag 0), within the allowed −2..+4: RH responds after rain.
- **Arms,** fitted per fold on training gauges:
  - A FY4B_QM;
  - B `FY4B_QM_RH`, set to 0 where RH is below the cut that most raises training CSI_0.1;
  - B8 `FY4B_QM_RH8`, the same but leaving values ≥ 8 mm/h alone;
  - C and C8, each one's amount-matched control;
  - D `GPM_uncal_RH`, the same gate on GPM_uncal, in no rule.
- **Rule:** the CLT screen's, for B and B8: ΔCSI_0.1 ≥ +0.01 over the control with its interval above 0, POD_8
  not significantly lower, and a gate in at least 3 of 5 balanced folds. GO if either passes.

### Result: GO for the gate that spares heavy rain, and a small gain

The folds pick RH < 70 (RH < 75 in balanced fold 4), and RH < 75 for GPM_uncal. Held out, balanced spatial:

| Series | r | RMSE | bias | CSI_0.1 | POD_0.1 | FAR_0.1 | FBI_0.1 | POD_8 |
|---|---|---|---|---|---|---|---|---|
| FY4B_QM | 0.247 | 0.949 | −0.0152 | 0.175 | 0.296 | 0.702 | 0.99 | 0.086 |
| B gate | 0.253 | 0.898 | −0.0400 | 0.188 | 0.275 | 0.626 | 0.73 | 0.064 |
| C control | 0.244 | 0.950 | −0.0202 | 0.163 | 0.250 | 0.681 | 0.78 | 0.086 |
| B8 gate | 0.249 | 0.933 | −0.0345 | 0.189 | 0.277 | 0.626 | 0.74 | 0.086 |
| C8 control | 0.245 | 0.950 | −0.0186 | 0.165 | 0.259 | 0.688 | 0.83 | 0.086 |
| GPM_uncal | 0.437 | 0.800 | −0.0119 | 0.338 | 0.602 | 0.565 | 1.38 | 0.102 |
| GPM_uncal gate | 0.419 | 0.799 | −0.0327 | 0.366 | 0.541 | 0.469 | 1.02 | 0.083 |

- **B fails:** ΔCSI_0.1 over C +0.0256 (+0.0189..+0.0326), but POD_8 −0.0216 (−0.0441..−0.0056).
- **B8 passes:** ΔCSI_0.1 over C8 +0.0243 (+0.0174..+0.0312). Its POD_8 cannot fall, by construction.
- **Against FY4B_QM, B8 gives:** CSI_0.1 +0.0143 (+0.0095..+0.0189), CSI_2 +0.0054 (+0.0013..+0.0098), RMSE −1.70%
  (−2.46..−1.01), r +0.0014 (−0.0040..+0.0074).

**How to read the GO:**
- **The gain is small and r does not move.** FY4B's real deficit, half of GPM_uncal's POD, gets worse:
  POD_0.1 0.296 → 0.277.
- **Bias worsens.** |bias| +0.0193 (+0.0142..+0.0237), so the bias goes to about 30% of mean gauge rain
  (0.114), and FBI_0.1 from 0.99 to 0.74. This rule has no bias test. The tolerance used in the blend screen
  (+0.005) would fail it.
- **Not specific to FY4B.** GPM_uncal gains more, ΔCSI_0.1 +0.0281 (+0.0209..+0.0353), though its r falls by
  0.0185. FY4B_QM's CSI stays far behind (0.189 vs 0.338).
- **The episodes mostly survive B8.** Of the false-rain squared error in FY4B_QM's gross hours, B8 removes
  20.8% and B 56.5%, because much of it sits at ≥ 8 mm/h. Controls remove under 1%. The worst hour's
  network RH is about 90%, so no cut near the chosen ones touches it.
- **B's RMSE gain is mostly R1,** which holds the dry 2022 season: −13.5% in R1, −0.8% in R2, −2.5% in R3. For
  B8: −2.8%, −0.7%, −2.3%.
- **By season,** ΔCSI_0.1 over the control is significant in JJA, SON and, barely, DJF, but not in MAM.
- **Downstream is untested and looks unpromising.** `d2m_c` is already a local MGWR covariate in 8 of
  FY4B_QM's 10 folds (3 of 5 balanced, 5 of 5 random), and a raw gain of 17% became 0.8% in the blend.
- **Correction to the first look:** an in-sample check with RH at the hour's end gave r +0.013. Per fold with
  RH at the hour's start, B gives +0.0053 (−0.0176..+0.0261), so most of that was RH responding to the rain.

**Recommendation:** by its rule the gate has earned a benchmark run, but that run needs a new derived product
and flag, and these numbers do not suggest it will survive. Do the channel pilot first. The user's call.

### Next: the channel pilot needs two files first

No AGRI L1 file is on the T9, and their size and the archive's reach back to June 2022 are unknown. So the
pilot stops at its first gate until these two exist, the scans of two of the worst gross hours, both in R1
(the matching QPE files exist):
```
FY4B-_AGRI--_N_DISK_1330E_L1-_FDI-_MULT_NOM_20230318233000_20230318234459_4000M_V0001.HDF
FY4B-_AGRI--_N_DISK_1330E_L1-_FDI-_MULT_NOM_20220923073000_20220923074459_4000M_V0001.HDF
```
- **Where to put them:** a new folder on the T9 linked as `data/FY4B_L1`, never under `data/FY4B`:
  `find_nc_files` keys on start time only, as for CLT.
- **They decide:** whether the sample's download is practical, and they settle the file layout and channels
  before any reader is written. They also give a first look at whether these hours are cold cloud.
- **After that, per the plan:** one scan per hour from gross false rain (`gross_hours.csv`), real widespread
  rain matched by month and hour of day, and missed rain; counting tables on training gauges over IR1, water
  vapour minus IR1, and the split-window difference; a rule fixed before the run. The hold-out, table cells
  and rule changed on 2026-11-07: see Newest.

## Earlier (2026-11-06): where FY4B's error sits

The user asked what makes FY4B hard to improve. No paper number moves, and nothing here is cited yet:
the script is not a step in `paper/sources.toml`.
- **Script:** `scripts/diagnose_fy4b_error_budget.jl` (`5a8e9d5`), run clean, 25 s. It fits nothing.
  Log: `output/fy4b_error_budget.log`. Output: `output/fy4b_error_budget/`.
- **Sample:** the raw fields of the `fy4b_benchmark` run, on the cells all four products score
  (4,452,108 of 4,760,619, balanced spatial).
- **Check:** FY4B_QM and FY4B on their shared mask reproduce `output/fy4b_quality/paired.csv` in both
  schemes (balanced: n 4,476,387, RMSE 0.95005 and 1.14265).

### What it finds

**FY4B's number holds one channel's information.** NSMC's slides
(`img.nsmc.org.cn/PORTAL/NSMC/VIDEO/FY4B/20230818/09-Atmosphere_Properties.pdf`) describe the
satellite-only estimate as a lookup from one infrared brightness temperature (IR1), matched offline to
polar-orbiter microwave rain. The 15-minute files hold only `Precipitation` and `DQF`, with an algorithm
dated 2016-10-10. Every fix so far reworks that number, or a class from the same imager (CLT). AGRI's
other infrared channels are untried.

**The gap is whether it rains, and misses count nearly as much as false rain** (0.1 mm/h; the last three
columns are shares of squared error):

| Raw field | r | RMSE | POD | FAR | False rain | Misses | Hit amounts |
|---|---|---|---|---|---|---|---|
| FY4B | 0.223 | 1.143 | 0.311 | 0.710 | 39.8% | 20.3% | 39.9% |
| FY4B_QM | 0.246 | 0.951 | 0.296 | 0.701 | 28.0% | 30.2% | 41.9% |
| GSMaP | 0.380 | 1.074 | 0.624 | 0.598 | 13.6% | 9.5% | 76.9% |
| GPM_uncal | 0.436 | 0.802 | 0.601 | 0.563 | 11.8% | 20.9% | 67.3% |

- FY4B_QM's MSE gap to GPM_uncal is +0.262 (+0.167..+0.370): false rain +0.177 (+0.095..+0.273), misses
  +0.139 (+0.101..+0.183), hit amounts −0.054 (−0.096..−0.010). The last is partly because it has half as
  many hits, so never read it as FY4B's amounts being better.
- With all false rain removed by oracle, FY4B_QM's r is 0.386, still below GPM_uncal as it stands. An
  oracle is a bound, not a method.
- §3.4's "What remains is false rain" holds for the blend (51% of the gap below), not the raw field.
  Not changed; the user's call.

**It is not mainly misplacement.** r at the stations and for the network mean:

| | 1 h station | 24 h station | 1 h network | 24 h network |
|---|---|---|---|---|
| FY4B | 0.223 | 0.378 | 0.405 | 0.458 |
| FY4B_QM | 0.246 | 0.417 | 0.451 | 0.514 |
| GSMaP | 0.380 | 0.658 | 0.668 | 0.803 |
| GPM_uncal | 0.436 | 0.721 | 0.833 | 0.894 |

Averaging does not close the gap, and rank correlation shows the same (network-day Spearman 0.53 vs 0.76),
as does dropping each product's worst 10 days (0.53 vs 0.92). So §4's "location error" row is not
supported as the reason point r is below areal r.

**False rain comes in episodes that FY4B alone cannot flag.**
- 30 of its worst days hold 80.8% of FY4B's false-rain squared error (GPM_uncal 43.6%).
- Gross false-rain hours, with the product wet at ≥50% of gauges and the gauges wet at ≤5%: 398 for FY4B
  (1.98% of hours, 113 days), holding 38.7% of its false-rain squared error and 15.8% of all of it. FY4B_QM:
  362 hours, 38.4%, 11.1%. GPM_uncal: 94 hours, 6.6%, 0.8%. The thresholds were chosen after a first look;
  under the others in `gross_summary.csv`, FY4B's share of false rain runs 27.6–48.4%.
- They fall in every month but January and February, most in September and October, at all hours of the day.
- Worst hour, label 2023-03-19 08 BJT: FY4B's network mean is 16.3 mm/h, wet at every gauge; the gauges
  read 0.00 and GPM_uncal 0.08.
- On the two worst days looked at, the files' quality flags read 0, and the first scan of one hour on each
  day shows a coherent cloud-like shape, not stripes. Cold non-raining cloud read as rain would fit, but
  it is **not shown**: there are no brightness temperatures for these hours.
- FY4B's own amounts overlap: median 0.80 mm/h in gross hours, 1.19 in real widespread rain. Of hours FY4B
  calls wet at ≥90% of gauges, 24% are false, against 3% for GPM_uncal.

**The product changes between regimes.** FY4B's June–August total against the gauges is ×1.76 in R1 (2022),
×0.89 in R2 (2023) and ×0.59 in R3 (2024). Each regime has one year, so software, position and weather
cannot be separated.

**FY4B adds little that the others lack.** P(gauge wet) is 0.07 when only FY4B is wet, 0.40 when only
GPM_uncal is, 0.49 when both are and 0.03 when neither is.

**Not the limit:**
- The 20 mm/h cap: FY4B reaches it in 743 cells, against 120 for GPM_uncal with no such cap. Gauges reach
  20 mm/h in 1,420 cells, and FY4B is at the cap in 27 of them.
- Gauge noise: gauges within 5 km correlate at r 0.76 (71 pairs), at 5–10 km 0.63 (246), at 10–20 km 0.47 (918).

## Earlier (2026-11-06): FY4B in the fusion, on the same hours and seeds

The user asked for FY4B in the comparison on the hours common to all three satellites. The first
check (below, 2026-11-01) was confounded, so this replaces it as the answer on the fusion. No paper
number moves.
- **Why the first check was confounded:** not by tuning hours alone. `run_seed` uses a product's
  position in the list (`InterpolationBenchmarkFold.jl`), so adding FY4B also shifted every product's
  covariate-selection seed.
- **Design:** a new arm with the qmfusion run's product list, hours and seeds, with FY4B scored but
  left out of the merges and the agreement bands. So the two arms differ in the fusion alone.
- **Code, branch `fy4b/same-hours-fusion`:**
  - `d273b46`: `--fy4b-out-of-fusion` (suffix `_fy4bnofusion`), off by default;
  - `e59f367`: `scripts/compare_fy4b_fusion_same_hours.jl`, committed with its verdict before the run;
  - `bc4f801`: `reviewed` bumped on the 12 steps flagged by the `src/` edit, none rerun.
- **Run:** clean at `bc4f801`, 3 h 37 min, 10 GB:
  ```sh
  julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates \
    --satellite-wet-blend --blend-axis agreement_envelope --fused-anchor --with-fy4b \
    --fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv --fy4b-qm --fy4b-qm-adaptive \
    --fy4b-out-of-fusion
  ```
  Logs: `output/full_benchmark_fy4bw3h0_fy4bqm_fy4bqmadapt_fy4bnofusion.log`,
  `output/compare_fy4b_fusion_same_hours.log`. Output: `output/fy4b_same_hours_comparison/`
  (`paired.csv`, `fold_status.csv`, `selection.csv`, `provenance.csv`).
- **Checks passed:** the same 237 stations, 20,087-hour grid, gauges and folds. Every `raw`, `adw`,
  `mgwr` and `blend_mgwr` of GPM_uncal, GSMaP, FY4B and FY4B_QM is identical in both arms, NaNs and
  masks included, against the qmfusion run from `b55a3a9`. So the later commits moved none of them.

### Result: the fixed verdict is "keep FY4B out"

Agreement blends (`blend_agrenv_mgwr`), balanced spatial, shared cells, FY4B_QM in the fusion against
out of it. "First check" is the same field from the confounded comparison.

| Product | RMSE | ΔRMSE (95%) | First check | ΔCSI_8 | FBI_8 | Fails on |
|---|---|---|---|---|---|---|
| GPM_uncal | 0.7243 vs 0.7248 | −0.06% (−0.16..+0.04) | −0.24% | −0.0005 | 0.239 → 0.235 | RMSE bar, FBI |
| GSMaP | 0.7346 vs 0.7353 | −0.10% (−0.24..+0.04) | −0.23% | +0.0014 | 0.337 → 0.327 | RMSE bar, FBI |
| MERGED_MEAN | 0.7212 vs 0.7244 | −0.45% (−0.74..−0.15) | −0.46% | −0.0038 | 0.307 → 0.281 | FBI |
| MERGED_OLS | 0.7183 vs 0.7199 | −0.22% (−0.34..−0.11) | −0.63% | −0.0000 | 0.231 → 0.228 | RMSE bar |

- **MERGED_OLS:** a real gain, negative in all five folds (−0.12..−0.38%), but under the 0.3% bar;
  8 mm/h skill is flat. Most of the first check's −0.63% was the confound, not its fold set: the clean
  effect on folds 3–5 alone is −0.19%. Its FBI failure was also the confound (+0.033 → +0.003).
- **MERGED_MEAN:** clears the RMSE bar in every fold (−0.05..−0.75%), but heavy rain gets worse:
  POD_8 0.110 → 0.104, FAR_8 0.641 → 0.630, |FBI_8 − 1| +0.026 (+0.004..+0.051). Unlike MERGED_OLS in
  the first check, this failure is not fewer false alarms alone.
- **Other fields:** MERGED_MEAN's `raw` (−4.21%) is the only field that meets the rule. Its `mgwr`
  (−3.91%) fails on CSI and FBI. MERGED_OLS's `mgwr`, which passed in the first check at −0.67%, is
  −0.21% and no longer does.
- **Random split:** every agreement blend is within ±0.11%, and every interval spans 0. Merged `adw`
  (control) is identical in balanced spatial and moves −0.18% in random, where the rule does not apply.
- **Selection:** a merged product's covariates are selected on its own anchor's residuals, which the
  intervals do not cover. MERGED_OLS selects the same in 5 of 5 folds in both schemes; MERGED_MEAN in
  3 of 5 balanced (folds 3 and 4 differ) and 5 of 5 random. It decides nothing here, since neither
  passes, and MERGED_MEAN's gain also appears in the folds that select the same.
- **FY4B's own blends (information, in no verdict):** FY4B −2.03%, FY4B_QM −1.27% with FY4B in the
  bands. Out of them FY4B_QM's blend is 0.7487, behind `adw` at 0.7459. So its agreement blend depends
  on FY4B's own signal in the bands. Never read this as FY4B helping the fusion.

### Fold coverage and a check on Table 3

Neither arm has a failed joint-method fold for any product (only fold 4 `partial`), so GPM_uncal and
MERGED_OLS score 4.5–4.6 M cells here against the paper run's 2.7 M on these hours. That gives the
paper's method over all five folds, on the common hours (blend vs `adw`, out-of-fusion arm):

| Product | All five folds | Folds 3–5 | Folds 1–2 | Paper run, all hours | Paper run, common hours |
|---|---|---|---|---|---|
| GPM_uncal | −2.79% | −2.95% | −2.56% | −2.78% | −2.72% |
| GSMaP | −1.25% | −1.40% | −1.04% | −1.30% | −1.12% |
| MERGED_MEAN | −3.01% | −2.97% | −3.07% | −3.17% | −3.00% |
| MERGED_OLS | −3.63% | −3.52% | −3.79% | −3.07% | −2.91% |

- Table 3's IMERG gain does not rest on its missing folds: folds 1–2 give a similar gain.
- MERGED_OLS's gain moves by about 0.6 points between runs on much the same folds, the others by
  about 0.2 or less. That is hours and seeds together; this run does not split them.
- The paper still does not state that Table 3's IMERG Uncal and merged rows rest mostly on folds 3–5.

**Not done, the user's call:** citing any of this. `compare_fy4b_fusion_same_hours.jl` is not a step in
`paper/sources.toml`, so a cited number needs a `[[step]]` and a `[[table]]` first.

## Earlier (2026-11-02): one FY4B series throughout the paper

The user asked to rerun RQ2 on the rebuilt series and update the paper (the follow-up offered in
Previous, last bullet).
- **Script** (`8c5709b`): `run_satellite_temporal_evaluation.jl` gains `--fy4b-file <path>`
  (with `--with-fy4b`), parsed as the benchmark's option. The run writes to
  `..._fy4bw3h0_gpmuncal`. Without the flag, the files, the folder and the outputs are unchanged.
  CLAUDE.md's RQ2 line names the flag.
- **Run** (clean, at `8c5709b`, about 40 s):
  `julia --project=. scripts/run_satellite_temporal_evaluation.jl --with-fy4b --fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv`
  writes `output/satellite_temporal_evaluation_fy4bw3h0_gpmuncal/` (log in `output/logs/rq2_fy4b_w3.log`).
  - The `all_products` sample: 19,818 hours and 298,653 gauge-wet cells, against 13,243 and
    214,803 on the first series.
  - Scope numbers (FY4B vs GPM_uncal, GSMaP), with the first series in brackets:
    - median hourly r at the gauges: 0.226 vs 0.450, 0.386 (0.229 vs 0.438, 0.374);
    - areal-mean hourly r: 0.405 vs 0.834, 0.667 (0.433 vs 0.833, 0.657);
    - wet-hour POD: 0.311 vs 0.602, 0.624 (0.344 vs 0.639, 0.656).
  - FY4B still scores worst on all three, so "scores worst" stays.
- **Checked:**
  - The `gpm_gsmap_full` rows of every output table equal the old `--with-fy4b` run's.
  - They also equal the `_nofy4b` run's, except the bootstrap intervals in
    `event_timing_summary.csv`. Those differ in the old `--with-fy4b` run too: the FY4B sample
    draws first, so the random stream shifts.
- **`paper/sources.toml`:**
  - `scope_fy4b_*` tables point at the new folder;
  - `rq2_fy4b`'s command gains the `--fy4b-file` argument, it needs `fy4b_rebuild`, and its
    `reviewed` is dropped (a fresh run);
  - `rq2`'s `reviewed` moves to `8c5709b`, with a note (the flag is opt-in, and the rows are
    identical).
- **Paper text:**
  - The Scope now says "In the hourly series rebuilt for this study (§2.5), on the hours it covers".
  - §2.5 says the rebuilt series is "the FY4B of the Scope paragraph and the plain FY4B of §3.4".
    The warning that its scores differ from the Scope's is gone.
  - The header comment and the `scope_fy4b_r` loader comment name the new run.
  - Then, at the user's request, §3.4 adds the MGWR stage against GSMaP: FY4B_QM is better,
    −6.8% (−11.2..−2.9), though its r stays lower (0.47 vs 0.49). After the blend it is +0.5%
    (the interval spans 0). The next sentence now says "the gain over plain FY4B after MGWR".
  - The Scope's r is still the median over gauges, and §3.4's raw r is still pooled over all
    gauge-hours. Each passage names its own.
- **Checks:**
  - `check_paper_provenance.py`: every entry traces.
  - `make_paper.jl`: 15 REUSE, 11 OK, 2 ABSENT, none CHANGED.
  - `typst compile`: clean, with no red "?".

## Previous (2026-11-02): FY4B written into the paper

The user asked for the FY4B paragraph, on the framing decided below ("raising FY4B's quality").
- **Paper text** (`paper/paper_plan.typ`):
  - Scope (Introduction) rewritten: RQ1–RQ3 compare GPM and GSMaP. FY4B is treated on its own;
    the uncorrected scores stay, and §3.4 asks how far its quality can be raised.
  - New §2.5 "Raising FY4B's quality" (method, no numbers): the rebuild, the per-fold shift per
    regime, and the adaptive quantile map.
  - New §3.4 "FY4B — how far can its quality be raised?", one paragraph:
    - gains against plain FY4B at every stage;
    - the remaining gap to GPM, and against GSMaP an interval spanning 0;
    - the random-split caveat;
    - the ice-cloud pilot;
    - "stays out of the merge and the agreement blend".
  - The Discussion's quantile-mapping sentence now says it is applied only to FY4B (§2.5), and
    not assessed on the peaks.
  - "What remains" gains the line "FY4B (§3.4): complete".
- **Numbers:** `scripts/compare_fy4b_quality.jl` (`9edc185`, run clean, 20 s), written to
  `output/fy4b_quality/paired.csv`.
  - It scores FY4B_QM against FY4B, GPM_uncal and GSMaP (raw, mgwr, agreement blend) in the FY4B
    benchmark run, with a paired day-block bootstrap (1000 reps, seed 20261101).
  - Its point estimates equal the scratch `rank_qmfusion.csv` (−16.86 / −5.36 / −0.80%, and
    +1.59 / +1.98% against GPM_uncal).
  - Balanced spatial:
    - against plain FY4B: raw −16.9% (−24.8..−10.0), mgwr −5.4% (−8.2..−2.7), blend −0.8%
      (−1.2..−0.4); raw r 0.22 → 0.25;
    - against GPM_uncal: blend +2.0% (+1.3..+2.6), mgwr +1.6% (−1.0..+4.5);
    - against GSMaP: blend +0.5% (−0.3..+1.3).
  - Random split, against plain FY4B: mgwr −0.1% (−0.6..+0.4).
- **`paper/sources.toml`:**
  - New steps:
    - `raw_fy4b` and `raw_fy4b_clt` (data);
    - `fy4b_rebuild`;
    - `fy4b_benchmark` (the qmfusion run's exact command);
    - `fy4b_quality`;
    - `fy4b_clt`.
  - New tables: `fy4b_quality`, `fy4b_clt` and `fy4b_clt_filter`.
  - `reviewed` bumped to `9edc185` after reading each diff:
    - `fy4b_benchmark`: `--blend-dry-band` is opt-in, `_tuning_spaced_hours` is an exact
      extraction, and `shift_fold` is a pure extraction;
    - `fy4b_rebuild`: none of its modules changed;
    - the seven existing steps: only new modules since `59e17d6`.
- **Corrected fact:** the processing switch at 133°E was on 2023-05-16 (last V1.0 hour 10:00 UTC
  label; `fy4b_hour_meta.csv`), not in June 2023 as §4.1 says. 105°E starts 2024-03-05.
- **Checked:**
  - `check_paper_provenance.py`: every entry traces to committed code.
  - `make_paper.jl` dry run: all steps REUSE or OK, none CHANGED.
  - `typst compile`: clean, with no red "?".
- **Two FY4B series in the paper (side-agent note, fixed in the text):**
  - **Scope's numbers:** they come from the first processed file. It has 13,471 hours and one pixel
    per gauge, and its r is the median over gauges.
  - **§3.4's "plain FY4B":** the rebuild, with 20,087 hours and a 3×3 mean. Its r is pooled over
    all gauge-hours.
  - **The fix:** each passage now names its series and its r. §2.5 says the rebuilt series is
    §3.4's plain FY4B and that its scores are not the Scope's. No numbers changed.
  - **If the tutor wants one series throughout:** rerun the RQ2 `--with-fy4b` evaluation on the
    rebuilt w3 file. That needs a file option in `run_satellite_temporal_evaluation.jl`. Done:
    see Newest.

## Earlier (2026-11-02): the CLT cloud-type filter, a one-month pilot — stop

The user downloaded July 2023 of `CLT` (2,921 quarter files, 5.8 GB, on the T9 at
`/Volumes/T9/CLT`, linked as `data/FY4B_CLT`) and asked me to pick the open decisions.

**Decisions taken (2026-11-02, picked by me at the user's request):**
- **Cloud filter:** try it on the July 2023 month (this section). The answer is no.
- **Main comparison:** FY4B stays out (`_nofy4b`). Even at face value, adding it cleared the rule
  for no blend, and its grid costs 24% of the hours.
- **Decision rule:** unchanged; past verdicts stand. Changing the FBI_8 part after seeing results
  would fit the rule to the data. New screens fix their own rule before the run, as this one did.
- **Paper framing:** "raising FY4B's quality".
  - What was done: FY4B rebuilt from the raw 15-minute QPE, then a per-fold pixel shift and an
    adaptive quantile map (FY4B_QM).
  - What it gained: raw RMSE −17%, MGWR −5.4%, blend −0.8% against plain FY4B, all significant
    (§0a).
  - What it still lacks: its blend trails IMERG Uncal's by about 2%, and it stays out of the main
    fusion.
  - No credit to the 2023 software update or the 105°E position.
  - Its false rain sits under ice cloud, where the cloud type cannot separate it (below).

**Code** (`73a0906`, branch `fy4b/clt-screen`, merged as `1a61f2a`):
- `src/sources/FY4BCloudType.jl` reads CLT at the gauges, reusing `FY4BPreprocessing`'s file-name
  parsing, projection and hour labels.
  - Hours are strict: all four quarters or no value, labelled by the Beijing hour end.
  - Per station, shift and hour, it counts the classes over the 3×3 box × 4 quarters (36
    pixel-quarters). The dominant class needs at least 18 valid; a tie goes to the smaller code.
  - CLT files are kept out of `data/FY4B`, because `find_nc_files` keys on start time only and
    would take them for duplicates of the QPE quarters.
- Test: `test/test-fy4b-cloud-type.jl` (in `runtests.jl`). The full suite passes.
- **Checked on the real files:**
  - CLT and QPE share the grid and axis order (`y, x` dims, read as `[x, y]`).
  - Where QPE is ≥ 1 mm/h, CLT calls 0.0% of pixels clear, against 27% with the axes swapped.
    So FY4B's QPE already applies the cloud mask: it never rains on a clear pixel.
  - The module's counts equal a direct full-array read at two stations and two shifts.
- **Screen:** `scripts/screen_fy4b_clt.jl`, run at `73a0906` from a clean tree, 2 min. Log:
  `output/screen_fy4b_clt.log`. Output: `output/fy4b_clt_screen/`.
  - It reuses the two-stage screen's setup: the §0a run's grid, folds and per-fold FY4B_QM.
    Held-out FY4B_QM equals that run's `oof_raw.csv` in all 4,760,619 cells of both schemes.
  - CLT: 694 of 744 July hours complete (the daily 14:15 UTC gap and a few others). 688 of them
    are on the FY4B grid, all regime R2 (133°E V1.0.0). That gives 162,150 scored cells over
    32 Beijing days.
  - The class is read at each fold's own FY4B_QM shift, (0,−1) or (1,−1).
- **Arms**, fitted on the fold's training gauges and scored on held-out gauges:
  - A: FY4B_QM.
  - B: FY4B_QM set to 0 in greedily chosen classes, adding classes while training CSI_0.1 rises.
  - C: an amount-matched control that zeroes as many of the lightest training wet hours.
- **Rule, fixed before the run:** go only if B beats C by ΔCSI_0.1 ≥ +0.01 with its interval
  above 0, and B does not lose POD_8.

### Result: the rule fails

- **The filter is always the same:** every fold, in both schemes, picks one class, *overlap*
  (thin ice over lower water cloud). It removes about 5.5% of the training wet hours, and the
  matched cut is 0.139 mm/h.
- Balanced spatial, held out:

| Pair | ΔCSI_0.1 (95%) | ΔPOD_8 | ΔRMSE |
|---|---|---|---|
| B (CLT) vs A (FY4B_QM) | +0.0062 (+0.0020..+0.0104) | 0 | −0.00% |
| B (CLT) vs C (amount) | **+0.0056 (+0.0011..+0.0090)** | 0 | −0.02% (sig) |
| C (amount) vs A | +0.0007 (−0.0034..+0.0046) | 0 | +0.02% |

  - B beats the amount control significantly, but by half the bar, so **stop**. FY4B_QM stays
    the quality step.
  - Random split (not judged): B vs C gives +0.0049 (+0.0008..+0.0084).
  - POD_8 cannot move: overlap hours never reach 8 mm/h.
- **Why so little: FY4B's false rain is under ice cloud, the same class as its real rain.**
  FY4B_QM held out, balanced, per dominant class:

| Class | Share of cells | FY4B_QM wet | FAR_0.1 | GPM_uncal FAR_0.1 | Share of FAs | Share of hits |
|---|---|---|---|---|---|---|
| clear | 31% | 1 | — | 0.99 | 0.0% | 0.0% |
| water | 27% | 137 | 0.89 | 0.89 | 1.1% | 0.4% |
| supercooled | 15% | 163 | 0.61 | 0.71 | 0.9% | 1.7% |
| mixed | 11% | 1,269 | 0.74 | 0.74 | 8.2% | 8.9% |
| ice | 10% | 12,742 | 0.75 | **0.57** | **82.7%** | **88.3%** |
| cirrus | 0.2% | 23 | 1.00 | 1.00 | 0.2% | 0.0% |
| overlap | 5.6% | 829 | 0.97 | 0.89 | 7.0% | 0.7% |

  - FY4B_QM calls rain in 80% of ice-dominated hours (12,742 of 15,938); GPM_uncal does in 44%
    (6,935). That excess is FY4B's false rain. It shares a class with 88% of the hits, so no class
    filter can remove it without the hits.
  - Within ice hours the amount already sorts it. FAR is 0.86 at 0.1–0.5 mm/h, 0.80 at 0.5–2,
    0.71 at 2–8 and 0.41 at ≥ 8 (`class_by_amount.csv`), and the QM works on the amount.
  - Overlap is the one easy class (FAR 0.97), but it holds only 7% of the false alarms.
  - Cirrus barely occurs as the dominant class (0.2% of cells), so the "cirrus false rain"
    hypothesis has nothing to act on at this scale.
  - FY4B misses warm and supercooled rain: it is wet in 163 supercooled hours against GPM_uncal's
    3,310. That is misses, not false rain.
- **Time of day** (`diurnal.csv`, Beijing hour end): false alarms peak 18–22 h (680–830 a
  slot, FAR 0.67–0.71, mostly ice), with the hits. FAR is highest at midday (11–13 h,
  0.91–0.98), on few wet hours. The evening tail is ice-topped rain hours, not a separate class.
- **This month is not typical:** in July 2023, FY4B_QM's FBI_0.1 is 1.72 on these cells,
  against 0.99 over the full grid. The map matches the wet frequency per season and regime, not
  per month, and July 2023 was wet.
- **For the paper:** at most one sentence. A one-month pilot found FY4B's false rain under ice
  cloud, like its hits, so a cloud-type filter gains little over the amount.
- **Not done, by the rule:** the full CLT download (about 180 GB), a benchmark flag, and a nested
  run.

## Earlier (2026-11-01): two cheap screens, a time shift and the satellite position

The user asked for ways to improve FY4B. These are the two that need no new data. Both are
exploratory screens on branch `fy4b/timing-position`, run from a clean tree at `cec3aab`, a minute
or two each. Each rule was fixed before its run. After two side-agent notes (the viewing angle,
and the 2022 drought), the position script was rerun from a clean tree at `3537881` with a
dry-year pair and the gauges' rain per year (branch `fy4b/drought-check`). Every earlier row is
unchanged.

### 1. A whole-hour time shift: no

- **Script:** `scripts/screen_fy4b_time_shift.jl`. Log: `output/screen_fy4b_time_shift.log`.
  Output: `output/fy4b_time_shift_screen/`. It scores the gauge at hour t against the rebuilt FY4B
  (w3, before shift and QM) at t+k, k = −3..+3, per season × regime.
- **Summer and autumn:** the best k is 0 in every regime. The summer composite peak is still late
  on the rebuilt file (diurnal phase +1.4 h, as RQ2 found on the old file). But hour by hour,
  FY4B matches best at the same hour. So the late peak is not shifted rain. It is consistent with
  FY4B keeping rain going into the evening, i.e. false rain, though that is not shown directly.
- **Spring and winter:** FY4B's previous hour matches slightly better (FY4B early). MAM k −1 gives
  Δr +0.006, and DJF k −1..−3 gives up to +0.02 within a regime.
- **Held-out test:** k is picked per fold on training stations and scored on held-out stations.
  - Pooled Δr is −0.0000 (−0.0013..+0.0013): **the rule fails**.
  - CSI_0.1 rises by +0.003, significant but tiny.
  - DJF alone: Δr +0.011 (−0.013..+0.047), not significant.
- **Control:** GPM_uncal's best k is 0 everywhere except a trace in DJF. So the method and the hour
  labels check out.

### 2. 105°E against 133°E: not shown

- **Script:** `scripts/compare_fy4b_by_position.jl`. Log: `output/compare_fy4b_by_position.log`.
  Output: `output/fy4b_by_position/`.
  - It reads the qmfusion run's out-of-fold files.
  - It compares regimes on JJA+SON only, so each regime contributes one year.
  - GPM_uncal's same field is the yardstick for each year's weather.
  - Check: FY4B_QM raw on its own mask reproduces `metrics_pooled.csv` (0.950043 balanced).
- **Position (R2 → R3, same software V1.0.0), balanced, FY4B_QM's gap to GPM_uncal:**

| Field | RMSE gap R2 → R3 | change (95%) | r gap R2 → R3 | change (95%) |
|---|---|---|---|---|
| `raw` | +15.8% → +8.3% | −7.5 pts (−17.3..+1.5) | −0.162 → −0.180 | −0.017 (−0.095..+0.065) |
| `blend_mgwr` | +1.5% → +1.8% | +0.3 pts (−1.5..+2.0) | −0.018 → −0.028 | −0.010 (−0.034..+0.014) |

  **The reading fails.** Plain FY4B's RMSE gap does shrink a lot, +45.8% → +15.2% (−30.6 pts,
  −62.1..−0.8), but its r gap does not. That is consistent with FY4B's amounts changing in 2024
  (the move, or that year's weather), which the QM already corrects.
- **Don't credit the viewing angle.** Over our gauges (109.6–111.5°E, 31.3–33.2°N) the viewing
  zenith only drops from 44–46° (133°E) to 37–39° (105°E). The shift of a 10 km cloud top goes
  from 10.0 to 7.8 km. 105°E is not "straight on".
- **Unexpected, and not the question set before the run:** from R1 to R2 (the May 2023 software
  update, both 133°E), FY4B_QM's own-blend gap to GPM_uncal narrows in balanced spatial:
  - RMSE gap +3.6% → +1.5%, change −2.1 pts (−4.1..−0.2);
  - r gap −0.049 → −0.018, change +0.030 (+0.008..+0.055);
  - `mgwr`'s r gap also changes by +0.058 (+0.004..+0.107).

  In the random split the same changes are not significant.
- **But R1 → R2 also swaps a dry year for a wet one.** Gauge rain, Jun–Nov, all gauge hours
  (`rain_years.csv`):

| Year | Wet hours (≥0.1 mm/h) | Hours ≥2 mm/h per 1000 | Mean, mm/h |
|---|---|---|---|
| 2014–2021 | 5.96–11.08% | 11.3–35.8 | 0.091–0.211 |
| 2022 (R1) | 6.08% | 15.6 | 0.112 |
| 2023 (R2) | 10.10% | 33.8 | 0.206 |
| 2024 (R3) | 6.45% | 17.6 | 0.119 |

  - 2022 was dry, and 2023 was the second-wettest Jun–Nov since 2014 by mean rain (after 2021).
    The scored
    cells show the same (`weather.csv`).
  - FAR_0.1 follows the weather for both products: FY4B_QM raw 0.72 / 0.61 / 0.73 and GPM_uncal
    raw 0.61 / 0.50 / 0.61 in R1 / R2 / R3.
  - So GPM_uncal as the yardstick takes out only part of each year's weather.
- **Dry-year check (rule fixed before the rerun).** R1 → R3 (2022 → 2024, both dry) holds the
  weather roughly fixed. It changes the software and the position. FY4B_QM `blend_mgwr` against
  GPM_uncal `blend_mgwr`, balanced:
  - RMSE gap +3.6% → +1.8%, change −1.8 pts (−4.0..+0.1);
  - r gap −0.049 → −0.028, change +0.021 (−0.005..+0.049).

  **The rule fails, so the update is not shown.**
  - The point estimates are nearly as large as from R1 to R2. So the data don't show that the
    gain was the weather either. With one year per regime, the two can't be told apart.
  - In the random split the change is +0.1 pts (−0.5..+0.7).
  - FY4B_QM `raw`'s RMSE gap does shrink between the dry years, −21.1 pts (−41.1..−2.5), but its
    r gap does not, +0.031 (−0.046..+0.125). So again this is amounts.
- **Also:** in 2023 (R2, the wet year), FY4B_QM's own blend beats ADW, −1.24% (−2.30..−0.21). In
  2022 and 2024 it does not.
- **For the paper:**
  - Don't credit the software update or the viewing angle.
  - At most, say that FY4B_QM's blend was closer to GPM_uncal in 2023 and 2024 than in 2022.
  - Say also that with one year per regime, and 2023 unusually wet, this can't be attributed to
    the update.

### What this leaves

- Timing is not a lever. The remaining target is false rain: FAR_0.1 is 0.61–0.73 for FY4B_QM raw
  against 0.50–0.61 for GPM_uncal raw (JJA+SON, by regime).
- **The idea aimed at it is a cloud filter using FY4B's own Cloud Type product, `CLT`.** The
  details are from NSMC's FY-4B AGRI L2 CLT product sheet V1.0.1, 2023-07.
  - Full disk, 4 km, every 15 min, NetCDF.
  - Same nominal grid and scans as our QPE files, so it lines up pixel for pixel, quarter by
    quarter.
  - File name: `FY4B-_AGRI--_N_DISK_1330E_L2-_CLT-_MULT_NOM_<start>_<end>_4000M_V0001.NC`, with
    `1050E` from March 2024.
  - Classes: 0 clear, 2 water, 3 supercooled water, 4 mixed, 5 ice, 6 cirrus, 7 overlap,
    9 uncertain (126 space, 127 fill).
  - Its `DQF` bits 1–2 carry the cloud-mask result (cloud / possible cloud / possible clear /
    clear), so `CLM` is not needed separately.
  - **Hypothesis, untested:** IR-based rain estimates often give false rain under thin cold ice
    cloud (cirrus, decaying anvils). That would fit the evening tail. Whether FY4B's QPE algorithm
    (version 2016-10-10 in the file) already uses cloud type is undocumented.
  - **Not checked:** whether NSMC's archive has CLT back to June 2022, and the file size. QPE
    files are about 0.9 MB each.
  - **Suggested first step:** the user downloads one month, e.g. July 2023. A small screen then
    checks whether FY4B's false-rain hours at the gauge pixels fall in some classes (cirrus,
    overlap) more often than its hits do. Only then is the full 2022–2024 download worth it. Any
    filter rule would be fitted per fold on training gauges.
- For the paper, the software-update result no longer supports "FY4B's quality is rising" (see
  the dry-year check above).

## Earlier (2026-11-01): does FY4B help the fusion?

**Superseded by Newest (2026-11-06).** Two things below were incomplete: the confound is hours *and*
selection seeds, and the run proposed at the end, restricted to FY4B's hours, would still differ by
seed. The numbers below are kept as the record of this check.

The user asked whether FY4B should go back into the main benchmark; the paper run leaves it out
(`_nofy4b`). This was checked from the two existing runs, with nothing refitted.
- **Script:** `scripts/compare_fy4b_in_fusion.jl` (`80287f2`, branch `fy4b/fusion-check`). It was
  run from a clean tree and takes under a minute. Log: `output/compare_fy4b_in_fusion.log`.
  Output: `output/fy4b_in_fusion_comparison/` (`paired.csv`, `provenance.csv`).
- **What it compares:** the paper run (`_nofy4b_gpmuncal`, `3ef424c`) against the qmfusion run
  (`b55a3a9`, FY4B_QM in the merges and bands). Each field is scored against the same field in the
  other run, on the cells both runs score: FY4B's 20,087 hours and both evaluation masks.
  - In balanced spatial, GPM_uncal and MERGED_OLS share 2.7 M cells. That is mostly folds 3–5,
    because of the paper run's joint-fold failures. The other products share 4.4–4.5 M.
- **Alignment checks passed:**
  - the same 237 stations;
  - the gauges are equal on the shared hours;
  - `split_common.csv` is equal for both schemes;
  - GPM_uncal's and GSMaP's `raw` fields are identical, so both runs read the same GSMaP file.
- **Sanity check:** GSMaP's agreement blend, scored on each run's own mask, reproduces
  `metrics_pooled.csv` (0.734575 and 0.687429) with the same n.

### Result: the fixed verdict is "confounded", and no agreement blend clears the rule anyway

- **Confound:** each run tunes ADW, the MGWR bandwidth and κ, and the blend weights on hours drawn
  from its own grid. GSMaP's `mgwr`, which FY4B never reaches, moved −0.56% (−1.08..−0.17). That
  is over the 0.3% that the verdict, set before the run, allows a control.
  - ADW is identical in balanced spatial.
  - The blends' own controls (`blend_mgwr`) moved −0.16% for GPM_uncal and −0.13% for GSMaP.

Agreement blends (`blend_agrenv_mgwr`), balanced spatial, with FY4B_QM against without:

| Product | RMSE | ΔRMSE (95%) | ΔCSI_8 | FBI_8 | Fails on |
|---|---|---|---|---|---|
| GPM_uncal | 0.7215 vs 0.7233 | −0.24% (−0.45..−0.04) | −0.0038 (sig) | 0.251 → 0.225 | RMSE bar, CSI, FBI |
| GSMaP | 0.7346 vs 0.7363 | −0.23% (−0.45..−0.04) | +0.0015 | 0.335 → 0.327 | RMSE bar |
| MERGED_MEAN | 0.7212 vs 0.7245 | −0.46% (−0.75..−0.17) | −0.0051 (sig) | 0.314 → 0.281 | CSI, FBI |
| MERGED_OLS | 0.7165 vs 0.7211 | −0.63% (−0.88..−0.39) | −0.0012 | 0.279 → 0.246 | FBI |

- **Taken at face value**, the verdict would be "keep FY4B out": no agreement blend meets the rule.
  MERGED_OLS fails only the FBI test. Its FBI_8 falls mostly through fewer false alarms
  (FAR 0.654 → 0.622), while POD goes 0.097 → 0.093.
- **Much of the change is the tuning hours, not FY4B:**
  - GPM_uncal's and GSMaP's controls move nearly as much as their agreement blends (−0.16% vs
    −0.24%, and −0.13% vs −0.23%). So FY4B's own part in those two blends is about 0.1%.
  - The controls' FBI_8 shifts by about as much as the treatments' (GPM_uncal's `blend_mgwr`
    0.252 → 0.223).
- **Where FY4B does help:** the merged anchors before blending.
  - MERGED_MEAN's raw merge improves by 4.2%, and its `mgwr` by 3.9%, but its POD_8 falls from
    0.193 to 0.149.
  - MERGED_OLS's `mgwr` improves by 0.67%, which meets the rule.
  - Blending toward ADW then takes most of that away.
- **Random split:** every agreement blend is within ±0.1%, and every interval spans 0.
- **Cost not in these scores:** the FY4B grid drops 24% of the hours (26,304 → 20,087).
- **So:** keep FY4B out of the main comparison. Separating FY4B from the tuning hours would need a
  new run of the paper setup restricted to FY4B's hours (new code, about 4 h). Nothing passes even
  at face value, so that run is worth it only if the paper needs a clean number.

## Earlier (2026-11-01): nested run of the dry-band blends

The user asked for a full nested run of both band-0 ideas from the screen below.
- **Code** (`8232154`, branch `fy4b/blend-headroom`): `--blend-dry-band` adds two opt-in axes on
  top of `--blend-axis agreement_envelope`. Both give the cells no source product calls wet (the
  agreement blend's never-blended band 0) a weight of their own, as band 10:
  - `blend_agrenvdry_*` leans band 10 toward adw, like every other band (the screen's K3);
  - `blend_agrenvdryraw_*` leans it toward the product's own anchor (the K5 idea).
- **Weights:** bands 1–9 repeat the agreement blend's weights exactly, because the bands separate
  in `select_blend_lambdas`. Band 10's weight is chosen on the inner split, on the representative
  tuning hours only, i.e. the evenly spaced half (`_tuning_representative`).
  - Why: the tuning sample's wettest half would make gauge rain in the dry band look commoner than
    it is. A side agent flagged this before the run; the stacked blend avoids it for the same
    reason.
  - So a null result would not come from the wet-heavy sample. The weight still rests on about
    168 hours.
- **Tests:** unit and end-to-end tests; the full suite passes. The end-to-end test checks that no
  existing row moves and that bands 1–9 repeat the agreement weights. The review bump is
  `53367fe`.
- **Run:** started at `8232154` from a clean tree with `-t 4` and the §0a command plus
  `--blend-dry-band`. Log: `output/full_benchmark_blenddry_qmfusion.log`. Output:
  `output/interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only_satwetblend_blendagrenv_blenddry_fusedanchor_fy4bw3h0_fy4bqm_fy4bqmadapt_qmfusion_gpmuncal/`.
- **Scored** with `julia --project=. scripts/compare_blend_dry_band.jl` (`fd59d5f`), run at
  `3e3a6f4` from a clean tree. Log: `output/compare_blend_dry_band.log`. Output:
  `output/blend_dry_band_comparison/` (`paired.csv`, `dry_weights.csv`, `provenance.csv`).
  - **Nothing else moved:** 7,344 pooled, 36,720 fold and 44,064 stratified reference rows are
    unchanged (largest relative difference 0), and so are the 360 `blend_selection.csv` rows.
    Bands 1–9 of both dry axes repeat the agreement weights.

### Result: neither dry-band blend clears the decision rule

FY4B_QM, balanced spatial, each against the agreement blend (95% day-block intervals):

| Blend | RMSE | ΔRMSE | Bias (mm/h) | ΔCSI_8 | FBI_8 | Fails on |
|---|---|---|---|---|---|---|
| agreement (`blend_agrenv_mgwr`) | 0.7391 | | +0.0003 | | 0.233 | |
| band 0 toward adw (`blend_agrenvdry_mgwr`) | 0.7376 | −0.20% (−0.28..−0.14) | +0.0017 | +0.0003 | 0.229 | RMSE bar, FBI |
| band 0 toward the anchor (`blend_agrenvdryraw_mgwr`) | 0.7339 | −0.71% (−0.90..−0.53) | −0.0124 | +0.0002 (ns) | 0.225 | bias, FBI |

- **Toward adw:** −0.20%, below the −0.3% bar. The screen's K3 gave −0.17%, so the cross-fitted
  screen predicted the nested result well.
- **Toward the anchor:** clears the RMSE bar, but Δ|bias| is +0.0121 (+0.0050..+0.0141), over the
  0.005 tolerance. The bias goes from about 0 to −0.0124 mm/h, about 11% of mean gauge rain
  (0.114 mm/h). Band 0 is where no product is wet, so the anchor there is near 0. Pulling toward
  it takes out the rain the gauges do get in those cells.
- **Each check on its own** (a side agent asked for this before dropping the idea):
  - toward adw: passes CSI_8 and bias, fails the RMSE bar and FBI;
  - toward the anchor: passes RMSE and CSI_8, fails bias and FBI.

  Without the FBI test, both would still fail.
- **FBI_8:** both variants fail "|FBI_8 − 1| not significantly worse". FBI_8 is already 0.23. The
  dry band cuts 8 mm/h false alarms (FAR 0.628 → 0.621 / 0.615) at the same POD (0.087), and that
  moves FBI further from 1. So the test counts fewer false alarms as worse here. For FY4B_QM this
  decides nothing, because each variant also fails another part of the rule.
- **Other products:** not FY4B_QM-specific; GSMaP gains most. Toward the anchor, plain FY4B
  (−0.96%) and GSMaP (−1.25%) fail only the FBI test. FY4B's bias flips sign (+0.0067 → −0.0073).
- **Edge over adw** (FY4B_QM, balanced): 0.91% with the agreement blend, 1.11% (−1.81..−0.48)
  toward adw, and 1.61% (−2.32..−0.93) toward the anchor, at the bias cost above.
- **Gap to GPM_uncal's blend** barely moves: 2.0% (agreement), 2.0% (toward adw), 1.7% (toward the
  anchor). Band 0 is not where FY4B_QM loses to GPM_uncal.
- **Band-0 weights:** FY4B_QM's balanced folds pick 0.2–0.5 toward adw and 0.6–0.7 toward the
  anchor. Over all products: 0–0.8 (mean 0.36) and 0.4–0.8 (mean 0.64). The weights are not stuck
  at 0, and they were picked on the representative hours, so the small gain toward adw does not
  come from the wet-heavy sample.
- **Random split** (not judged): both beat the agreement blend in RMSE, toward adw −0.44% and toward
  the anchor −0.75%. The anchor variant has bias −0.0095 and FBI_8 0.632 → 0.615.
- **So** the agreement blend stays the headline. The branch is not merged; whether to keep the
  `--blend-dry-band` code is the user's call.
- **Also on the branch:** `DayBlockBootstrap` (`c66893f`), the day-block scoring now shared by the
  screen and the comparison. The screen re-run gave byte-identical CSVs.

## Latest (2026-11-01): why the QM barely moves the blend, and a blend screen

- **Branch `fy4b/blend-headroom`:** `4d44e59` (the script), then this note. Not merged.
- **The question** (user): the QM cuts raw FY4B's RMSE by 17%, but the blend gains only 0.8%. Why?
  Is FY4B corrected more than once?
- **Answer, traced in code:**
  - **No.** The QM is applied once per fold. No later step maps, thresholds or rescales `FY4B_QM`.
  - **But the steps overlap.** Three gauge-fitted steps, on the same training gauges, go after the
    same amount error: the QM, then the MGWR residual with κ, then the blend weight toward ADW.
  - **MGWR redoes most of the QM's work.** Within an hour the QM is one monotone map, applied to
    every station alike. MGWR is a residual model, ĝ = max(0, s + κ·r̂), and the hourly residual
    surface from neighbouring gauges already fixes the spatially coherent amount errors.
  - **The blend leans on ADW**, because MGWR alone loses to ADW for every single-satellite product
    (FY4B, FY4B_QM, GPM_uncal, GSMaP). The merged anchors' MGWR does beat it: MERGED_MEAN 0.7412
    and MERGED_OLS 0.7244, against ADW ≈ 0.746.
  - **One subtlety:** the shift, QM and strength are fitted on all of a fold's training stations,
    including the inner hold-out group that κ and λ are tuned on.
- **Screen:** `scripts/screen_fy4b_blend_headroom.jl` → `output/fy4b_blend_headroom/` (log
  `output/fy4b_blend_headroom.log`).
  - Run at `7073072` from a clean tree; it takes 0.6 min. The first run, at `4d44e59`, printed
    the same numbers.
  - It reads the §0a run's out-of-fold files. As a guard, it rebuilds both blends for FY4B,
    FY4B_QM and GPM_uncal in all 4,760,619 cells of both schemes, and its scores match
    `metrics_pooled.csv`.

**Attribution: the run's own nested rows, RMSE against ADW, 95% day-block interval.**
- The own-anchor blend `blend_mgwr` uses only the target product.
- The three-product blend `blend_agrenv_mgwr` bands on FY4B_QM, GPM_uncal and GSMaP.

| balanced spatial | own-anchor vs ADW | three-product vs ADW | three-product vs own-anchor |
|---|---|---|---|
| FY4B | +0.81% (−0.18..+1.77) | −0.11% (−0.97..+0.69) | −0.91% (−1.24..−0.63) |
| FY4B_QM | −0.38% (−1.19..+0.40) | **−0.91% (−1.61..−0.25)** | −0.53% (−0.75..−0.35) |
| GPM_uncal | −2.56% (−3.45..−1.71) | −2.84% (−3.63..−2.09) | −0.29% (−0.47..−0.13) |

- **QM against plain FY4B, same blend:** own-anchor −1.18% (−1.77..−0.64); three-product −0.80%
  (−1.19..−0.44).
- **Random split:** FY4B_QM's own-anchor blend is +0.98% (+0.35..+1.61), *worse* than ADW. Its
  three-product blend is +0.28% (−0.21..+0.77).
- **Reading:**
  - With FY4B's own information alone, FY4B_QM's blend is not shown to beat gauges alone: the
    interval spans 0.
  - The significant −0.91% needs the three-product bands, which give −0.53% of it. That part
    mixes three things:
    - the other products' rain calls;
    - 9 weights instead of 1;
    - the extra cells blended where only another product is wet.

    K4 below (FY4B-only bands with 3 weights) suggests the extra weights alone add little. So
    most of it is the other products, but not all of it is shown to be.
  - So never write "corrected FY4B beats gauges alone" on the strength of `blend_agrenv_mgwr`.
  - The QM's own effect is clearest in the own-anchor blend: −1.2% against plain FY4B.

**Where FY4B_QM's blend loses to GPM_uncal's** (balanced, `gap_by_band.csv`, share of the SSE
gap):
- FY4B_QM wet and gauge dry, i.e. false rain: 51%.
- FY4B_QM dry and gauge wet, i.e. misses: 20%.
- Both wet: 16%. Both dry: 13%.
- Band 0, where no product is wet: 15% net, and its gauge-dry cells alone give 26%.
  - That is consistent with FY4B's training-gauge errors reaching cells no product calls wet,
    through the residual surface. It is not shown directly.
  - In band 0 (`attenuation.csv`), FY4B_QM's blend has bias +0.006 and RMSE 0.2717; GPM_uncal's
    has −0.001 and 0.2649.
- Cells where only FY4B_QM is wet (bands 1–3): 17%.

**Headroom, exploratory.** λ is chosen on the other folds' out-of-fold cells, and those
predictions came from models that saw the held-out fold, so this is not leakage-free. Balanced
spatial, against K0:

| key | what it adds | RMSE | bias | CSI_8 |
|---|---|---|---|---|
| K1 | × FY4B_QM dry/wet | −0.00% (−0.03..+0.03) | −0.0013 | −0.0004 |
| K2 | × month mapped / unmapped | −0.02% (−0.10..+0.06) | +0.0000 | +0.0006 |
| K3 | band 0 blended (ID 10) | −0.17% (−0.27..−0.09) | +0.0029 | +0.0003 |
| K4 | FY4B_QM-only bands | +0.84% (+0.61..+1.12) | −0.0064 | −0.0069 (sig) |

- **No key meets the rule** (RMSE −0.3% or better, with heavy-rain detection and bias held).
  - "|bias| does not grow materially" is coded as at most +0.005 mm/h (`BIAS_TOLERANCE`, about
    4% of mean gauge rain). That number is my choice, not the user's. It decides nothing here,
    because K3 already fails the RMSE bar.
  - Every interval here is a day-block bootstrap. It holds the stations fixed, so it leaves out
    the spatial (station and fold) sampling, as the earlier screens do.
- **K3 is the only real gain, and it is small.**
  - Band 0's weight is 0.4–0.6 (balanced) and 0.7–0.8 (random). In the random split K3 gives
    −0.46% (−0.75..−0.16), as context only.
  - It is not FY4B-specific: it blends the all-dry cells of any product.
- **K1 and K2 add nothing**, so FY4B's own wet call and its month quality class don't help the
  weight.
- **K4 shows the bands need GPM and GSMaP.**
- **K0, the protocol control,** scores 0.7356, against the benchmark's own 0.7391. That gap mixes
  the leaky cross-fit, the tuning, and a protocol difference: K0 picks λ on all hours, unweighted,
  where the benchmark uses 336 weighted tuning hours. It is not a tuning-bias estimate.
- **Outcome:** stop, per the rule. Frame the paper on the improvement over ADW, about 0.1% → 0.9%:
  - more than half of it comes with the three-product bands, mostly the other products' rain
    calls (see the attribution reading);
  - FY4B's own part is the own-anchor row, whose interval spans 0.
- **Open for the user:**
  - K3 as a general "blend band 0" option. It is below the 0.3% bar, applies to every product,
    and would need a nested benchmark run before it could be cited.
  - Band 0 toward the anchor instead of ADW (a possible K5). In band 0 the raw anchor beats every
    corrected field: FY4B_QM raw 0.2613 vs mgwr 0.2717, GPM_uncal raw 0.2608 vs mgwr 0.2649
    (`attenuation.csv`). So where all three products are dry, the gauge correction adds error on
    net. This was spotted in held-out numbers, so it needs the exploratory screen and then a
    nested run.

## 0. Earlier (2026-11-01): two-stage counting tables, screened on raw FY4B

- **Branches:**
  - `fy4b/with-fy4b-benchmark` was merged into master as `805bad6` and pushed.
  - This work was on `fy4b/two-stage`: `ad9bc97` (code and screen), `1c58cca` (review bump),
    `0262a36` (first version of this note), `978ef74` (the amount-only ablation), later
    revisions of this note, then merged into master at the user's request.
- **User's decisions:** no machine learning, so counting tables; no 15-minute quarters yet; stop
  after the screen.
- **What it is** (`src/models/SatelliteTwoStage.jl`). Per fold, on training gauges only, each
  station-hour gets a class from FY4B alone, at the fold's shift:
  - intensity: 8 classes;
  - rain in the 8 surrounding shifts: none, some or all;
  - rain at t−1 / t+1: neither, one or both;
  - all within regime × season.

  Per class it counts P(gauge ≥ 0.1) and the mean wet amount. A class with fewer than 200 pairs
  pools to its regime, then all cells, then its intensity class.
  - `FY4B_TS_E` = P × amount.
  - `FY4B_TS_I`, the amount-only ablation: P × amount on the intensity class alone, with no rain
    nearby and no adjacent hours. A side agent asked for it: r rewards pulling FY4B's large, often
    false values toward typical gauge amounts, and that needs no context.
  - `FY4B_TS_O`: the classes with the highest P are called wet until they cover the training
    gauges' wet count. A wet hour keeps `FY4B_QM`'s amount.
- **Screen:** `scripts/screen_fy4b_two_stage.jl` → `output/fy4b_two_stage_screen/` (log
  `output/fy4b_two_stage_screen.log`).
  - Run at `978ef74` from a clean tree; it takes 2 min.
  - It uses the §0a run's grid, folds and per-fold `FY4B_QM`. Its held-out `FY4B_QM` equals that
    run's `oof_raw.csv` in all 4,760,619 cells of both schemes.
  - It scores cells where all seven series are finite, which is not the benchmark's mask.

Held-out raw scores, balanced spatial:

| series | r | RMSE | CSI 0.1 / 2 / 8 | FBI 0.1 / 2 / 8 |
|---|---|---|---|---|
| FY4B | 0.224 | 1.142 | 0.176 / 0.124 / 0.047 | 1.08 / 0.99 / 2.00 |
| FY4B_QM | 0.247 | 0.949 | 0.175 / 0.123 / 0.055 | 0.99 / 0.87 / 0.65 |
| FY4B_TS_I | 0.279 | 0.827 | 0.122 / 0.047 / 0.000 | 3.58 / 0.14 / 0.00 |
| FY4B_TS_E | 0.290 | 0.825 | 0.181 / 0.051 / 0.000 | 1.92 / 0.16 / 0.00 |
| FY4B_TS_O | 0.244 | 0.951 | 0.173 / 0.122 / 0.054 | 0.95 / 0.88 / 0.65 |
| GPM_uncal | 0.437 | 0.800 | 0.338 / 0.266 / 0.075 | 1.38 / 0.77 / 0.47 |
| GSMaP | 0.382 | 1.077 | 0.323 / 0.260 / 0.075 | 1.56 / 1.05 / 1.63 |

- **Paired** (balanced spatial, day-block bootstrap, 95%):
  - TS_E against `FY4B_QM`: r +0.043 (+0.021 to +0.066), RMSE −13.1% (sig). Random split:
    r +0.050, RMSE −13.0%.
  - It splits into two parts:
    - TS_I against `FY4B_QM`, the amount table alone: r +0.032 (+0.011 to +0.054), RMSE −12.8%;
    - TS_E against TS_I, what the context adds: r +0.011 (+0.005 to +0.017), RMSE −0.33%
      (−0.0045 to −0.0012 mm/h). Random split: r +0.013, RMSE −0.41%.
  - TS_O against `FY4B_QM`: r −0.003 and RMSE +0.2%; both CIs span 0.
- **Standing:** TS_E's r is 0.147 below GPM_uncal and 0.092 below GSMaP. Its RMSE is +3.1%
  against GPM_uncal and −23% against GSMaP (all sig).
- **Reading:**
  - **Most of the gain is not context.** Three quarters of TS_E's r gain over `FY4B_QM`
    (+0.032 of +0.043) comes from the amount table alone: it replaces FY4B's values with the
    training gauges' mean rain per intensity class, which pulls FY4B's large, often false values
    down toward typical gauge amounts. The context adds +0.011 in r, significant but small, and
    0.3% in RMSE.
  - Where the context does show is light-rain occurrence: CSI at 0.1 mm/h 0.122 → 0.181, and
    FBI_0.1 3.58 → 1.92 (TS_I against TS_E).
  - TS_E's gain over `FY4B_QM` comes from summer and autumn: +0.02 to +0.05 per regime × season.
    Winter and spring gain nothing.
  - TS_E is a conditional mean, so it flattens the rain. It never reaches 8 mm/h (CSI_8 and
    FBI_8 are 0), CSI at 2 mm/h falls from 0.123 to 0.051, and it spreads drizzle (FBI_0.1 is
    1.92). ROADMAP warned of exactly this: better r and RMSE, heavy rain lost.
  - The occurrence rule adds nothing over `FY4B_QM`, whose threshold already matches the gauges'
    wet frequency.
  - Context sees FY4B's misses better than its false rain. In 133E V1.0 JJA (fold 1):
    - a dry FY4B hour with rain nearby or in both adjacent hours is wet at the gauge 10–12% of
      the time, against 2% when isolated;
    - even FY4B ≥ 16 mm/h with rain all around and in both hours is wet only 42% of the time.
- **Fixed before the recorded run:** TS_O's frequency rule first ran on into the huge "dry,
  nothing nearby" class.
  - Where: three winter cells and one spring cell, where FY4B's wet classes cover fewer hours
    than the gauges are wet.
  - Effect: on training pairs, 8–115 times the gauges' wet count was called wet; held out,
    FBI_0.1 was 4.3.
  - Fix: the rule now stops before a class that would leave the count further from the target.
  - Spotted from the held-out FBI in the trial run; the cause and the fix are both read from the
    training counts in `context_tables.csv`.
- **Decided (2026-11-01): stop.** The user chose option 3 below after seeing the screen.
  `FY4B_QM` stays the FY4B quality step. No `--fy4b-twostage` flag, benchmark run or quarter
  rebuild was made; `SatelliteTwoStage` and the screen stay in the repo as the record.
- **The options that were on the table:**
  1. Add `FY4B_TS_E` to the benchmark (`--fy4b-twostage`, a derived product like `FY4B_QM`) and
     do the ~4 h full run. It would show whether a flatter anchor with a higher r helps MGWR and
     the agreement blend. TS_O is not worth a run.
  2. Add the 15-minute quarter context: a rebuild extension, ~18 min on the T9. The rule was "only
     if context helps". The context we already have adds +0.011 in r, so more of it is unlikely to
     add much.
  3. Or stop: `FY4B_QM` stays the FY4B quality step. The screen then records that a per-class
     conditional mean raises r but loses heavy rain, and that context adds little.
  - Do not write "context re-ranks hours" in the paper: the ablation shows that is a small part
    of the gain.

## 0a. Earlier (2026-10-31): adaptive map and QM-fed fusion scored

- **Commits on `fy4b/with-fy4b-benchmark`** (merged into master as `805bad6`, 2026-11-01), after §0b's:
  - `48f61df` this note; `4896f92` `ROADMAP.md` (the two-stage occurrence/amount model);
  - `e683517` adaptive map, `cdbd860` CLI flags `--fy4b-qm-adaptive`, `--fy4b-qm-fusion`;
  - `2e9f0e4` strength chosen per month, not per season;
  - review bumps `693b2bf`, `b55a3a9`.
- **What the two options do:**
  - `--fy4b-qm-adaptive` (`_fy4bqmadapt`): in each fold the training gauges pick, per regime ×
    year-month, how strongly the season map applies: none, threshold only, or 25/50/75/100%.
    The map itself is still fitted per season × regime. Choices: `fy4b_qm_strengths.csv`.
  - `--fy4b-qm-fusion` (`_qmfusion`): each fold's `FY4B_QM` replaces plain FY4B in
    MERGED_MEAN/MERGED_OLS and in the agreement-envelope bands.
- **Full run:** `..._fusedanchor_fy4bw3h0_fy4bqm_fy4bqmadapt_qmfusion_gpmuncal`, at `b55a3a9`, clean
  tree, the §0b flags plus both options. FY4B's and GPM_uncal's raw and mgwr rows are identical
  to §0b's run (checked cell by cell).
- **Analysis:** `output/fy4b_partB_analysis/qmfusion/` (`paired_qmfusion.csv` against §0b's
  run, `rank_qmfusion.csv` within this run, `monthly_qm.txt`, and their scripts).

Balanced spatial, change in RMSE; "sig" means the 95% day-block CI excludes 0:

| comparison | raw | mgwr | blend_agrenv_mgwr |
|---|---|---|---|
| adaptive `FY4B_QM` vs full-map `FY4B_QM` | −7.6% sig | −4.2% sig | −0.44% sig |
| adaptive `FY4B_QM` vs plain FY4B (w3) | −16.9% sig | −5.4% sig | −0.80% sig |
| MERGED_MEAN fed `FY4B_QM` vs fed FY4B | −4.6% sig | −1.2% sig | −0.42% sig |
| MERGED_OLS fed `FY4B_QM` vs fed FY4B | −0.07%, CI spans 0 | −0.12% sig | −0.11% sig |

- **July–August are repaired and June keeps its gain** (raw RMSE, plain / full map / adaptive):
  - Jun 1.261 / 1.182 / 1.171; Jul 1.551 / 1.707 / 1.528; Aug 1.438 / 1.664 / 1.397;
  - May also improves (1.001 → 0.944). Spring and autumn months are unchanged.
  - r 0.223 (plain) → 0.237 (full) → 0.246 (adaptive).
- **What the gauges chose, consistently across folds:**
  - 2022 summer (133E V1.0): the full map (Aug 2022: 0.75 in two folds of five);
  - summer 2023 and 2024 (May–Aug 2023 under V1.0.0, Jun–Aug 2024): no adjustment at all, not
    even the threshold; also Dec 2023–Feb 2024;
  - spring and autumn: the full map. The 8-hour "2025-01" cell is all dry and takes "none".
  - So the July–August loss came from the later two regimes, where FY4B's convective cells are
    displaced too much for amount mapping to help.
- **Standing** (paired, balanced spatial):
  - mgwr: `FY4B_QM` RMSE +1.6% vs GPM_uncal (95% CI −0.9% to +4.5%) and −6.8% vs GSMaP
    (−10.9% to −2.7%).
  - blend: +2.0% vs GPM_uncal (+1.4% to +2.7%), +0.5% vs GSMaP (−0.25% to +1.3%).
  - A CI that spans 0 can't rule out a tie, but it doesn't show one either. Cite the gap with its
    range, never "level with" or "tied". Against GPM_uncal after MGWR, most of the range is on
    the worse side. Against GSMaP in the blend, "within about 1%" is fair.
  - Every blend on one shared mask: MERGED_OLS 0.7189, MERGED_MEAN 0.7211, GPM_uncal 0.7257,
    GSMaP 0.7361, `FY4B_QM` 0.7400, FY4B 0.7459.
- **Fusion:** MERGED_OLS's FY4B weight rises from 0.022 to 0.063 (GPM_uncal 0.47, GSMaP 0.15).
  QM-fed bands leave the real products' blends where they were: FY4B +0.10% (CI spans 0),
  GPM_uncal 0.00%, GSMaP +0.08%.
- **Caveats:**
  - Random split: past raw, the gain mostly vanishes (mgwr −0.12%, CI spans 0). Neighbouring
    gauges in training already give MGWR the local field.
  - Raw CSI at 2 mm/h dips 0.128 → 0.123. The months left unmapped keep FY4B's false light rain.
  - The strength is chosen on the same month's training gauges, as the map is. This is
    in-period calibration under spatial CV, like the hourly GWR fits, not a forecast.

## 0b. Earlier (2026-10-31): Part B done, both full runs scored

- **Branch `fy4b/with-fy4b-benchmark`**, merged as `805bad6` on 2026-11-01. Since master:
  - rebuild code: `37e27d7`, `14d70d4`, `0924141` (duplicate-copy fix);
  - `FY4B_QM`: `ede88b8`; CLI flags `--fy4b-file`, `--fy4b-qm`: `3c6361d`;
  - review bumps `ef64be0`, `2ae2dbc`, `ca52142`, `d492f42`; this note.
- **Download:** the strict 4-of-4 record went from 13,474 to 20,087 complete hours.
  - Still missing: 2022 Jan–May, Feb 2024, scattered gaps, and the 14:15 UTC file every day. So
    the 22–23 Beijing hour never exists.
  - One corrupt file with no good copy: `...NOM_20220710071500_20220710072959...`.
- **Full runs:** `-t 4`, clean tree; the rebuild was reproduced byte for byte first. Both used
  `full --nested-covariates --satellite-wet-blend --blend-axis agreement_envelope --fused-anchor --with-fy4b`.
  - `..._fusedanchor_fy4bw3h0_fy4bqm_gpmuncal`: `--fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv --fy4b-qm`.
  - `..._fusedanchor_fy4bw1h0_gpmuncal`: `--fy4b-file .../fy4b_w1_h0.csv` (completeness only).
- **Paired analysis** on shared cells, with a day-block bootstrap: `output/fy4b_partB_analysis/`
  (`paired.csv`, `pooled.txt`, and the scripts that made them).

Balanced spatial, change in RMSE (mm/h); "sig" means the 95% CI excludes 0:

| step | raw | mgwr | blend_agrenv_mgwr |
|---|---|---|---|
| 3×3 window vs centre pixel (w3 vs w1) | −0.7% sig | −0.5% sig | −0.1% sig |
| shift + quantile map vs w3 (`FY4B_QM`) | −10.0% sig | −1.2%, CI spans 0 | −0.25%, CI spans 0 |

- **What the quantile map does to raw FY4B:**
  - bias +0.022 → −0.004 mm/h;
  - dry-hour RMSE 0.746 → 0.623;
  - r 0.223 → 0.237;
  - CSI at ≥ 8 mm/h 0.047 → 0.059.
- **Where it helps:** the gains are in Mar, Apr and Sep (RMSE roughly halves in Mar and Sep).
  Jul and Aug get worse (1.55 → 1.71, 1.44 → 1.66): the map restores the gauge tail onto
  misplaced convective cells.
- **The shift is a stable property of the product.** Every fold, in both schemes, picks the
  same or an adjacent shift (dx east+, dy south+):
  - 133E V1.0: dx −2/−3 (west);
  - 133E V1.0.0: (+1, −1);
  - 105E: (+2/+3, −1), i.e. east and north.
- **Standing:**
  - Raw `FY4B_QM` RMSE is now below GSMaP's (1.029 vs 1.079, CI spans 0).
  - Its r stays far lower: 0.24, against 0.38 for GSMaP and 0.44 for GPM_uncal.
  - With the headline blend, `FY4B_QM` stays 2.4% behind GPM_uncal and 1.1% behind GSMaP (both
    sig).
- **FY4B in the merges:**
  - MERGED_OLS gives FY4B a weight of about 0.02. With and without FY4B, the merged raw series is
    unchanged.
  - The blend moves by ≤ 0.5%, and that comparison is confounded by the two runs' different hour
    grids.
- **Reading:** the quality steps clearly raise raw FY4B. The GWR corrections absorb most of the
  same systematic error themselves, so the corrected fields barely move.

## 0c. Earlier (2026-10-31): Part A done, waiting for the download

- **Branch `fy4b/with-fy4b-benchmark`**, merged as `805bad6` on 2026-11-01:
  - `37e27d7`: rebuild code;
  - `ef64be0`: review bump;
  - `14d70d4`: offset grid;
  - `2ae2dbc`: review bump.

  Full test suite green, and `make_paper.jl` reports nothing CHANGED.
- **`scripts/prepare_fy4b_rebuild.jl`** reads `data/FY4B` (the T9 link) and writes
  `output/fy4b_rebuild_2022_2024/`:
  - `fy4b_{w1_h0,w3_h0,w5_h0,w3_h10}.csv` and their `_capfrac.csv`;
  - `fy4b_offset_grid.nc` (3×3 means at every pixel shift up to ±6);
  - the QC table, hour metadata and provenance.
- **Checks on June 2022:** `w1_h0` matches the shipped file in 159,027/159,027 cells, and the
  grid's zero shift matches `w3_h0`. Time is 0.7 min per month, so the full run is about 25 min.
- **Findings:**
  - The fixed-height parallax shift makes FY4B worse (June 2022: r 0.351 → 0.343 at 10 km).
  - Pixel-offset maps show a systematic displacement: 2022 west 2–5 px, 2023 none, 2024 east
    1–4 px.
  - The user chose (2026-10-30) to fix it with a **per-fold, per-regime shift picked on training
    gauges** (Part B).
- **One corrupt raw file:** `...NOM_20240608064500_20240608065959...` (HDF error). Later recovered:
  the download brought a good duplicate, which `0924141` now picks.
- **Tomorrow, once the download is on the T9:**
  1. Run `find /Volumes/T9/FY4B -name '._*' | wc -l`. The reader skips such files anyway.
  2. Rescan coverage with the scratch script `fy4b_missing_list.jl` (the plan has its path).
  3. Run the full rebuild: `julia --project=. scripts/prepare_fy4b_rebuild.jl`.
  4. Build Part B: the per-fold shift plus quantile map in the benchmark. The plan is
     `~/.claude/plans/mutable-dazzling-bengio.md`.

## 1. State at the restart

- `master` is at `1fbc1ea`, pushed, and the tree is clean (`## master...origin/master`).
- No Julia job was running. Only the julia-repl MCP adapter was up, and it held no state worth keeping.
- **Branch `fy4b/with-fy4b-benchmark` is kept for this work.** It has no commits and sits at `3662386`,
  before the GSMaP merge. Bring it up to master before using it:
  ```sh
  git checkout fy4b/with-fy4b-benchmark && git merge --ff-only master
  ```
- **To resume:** in `~/GGG/GGGG`, run `claude --continue` to reopen the same conversation, or start a
  new session and say "read `Projects_Pending_Handover/fy4b_quality_handover.md`".

---

## 2. Why this work exists

The tutor wants FY4B back in the research. The framing is **how to raise FY4B's quality**, not
dismissing or disparaging it. The paper currently leaves FY4B out, and its scope note says FY4B
"scores worst" (`paper/paper_plan.typ:174`). That note has to be rewritten once the improved FY4B
results exist.

---

## 3. What is already done

### 3.1 The with-FY4B benchmark (step 1, finished)

This is the command of the paper's `benchmark` step in `paper/sources.toml`, plus `--with-fy4b`:

```sh
julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates \
  --satellite-wet-blend --blend-axis agreement_envelope --fused-anchor --with-fy4b
```

- Output: `output/interpolation_benchmark_full_joint_covariates_nested_mgwrintercept_only_satwetblend_blendagrenv_fusedanchor_gpmuncal/`
- Log: `output/full_benchmark_fusedanchor_gpmuncal.log` (exit 0, about 2 h)
- Common grid: 13,471 hours, which is FY4B's own coverage. The no-FY4B paper run has 26,304.
- Still valid after the GSMaP fix (§3.2): full mode always read the full-year GSMaP file.

Pooled results, `balanced_spatial` scheme, `overall` group (`metrics_pooled.csv`):

| Product | raw r | MGWR r | blend_agrenv_mgwr r | RMSE raw → MGWR → blend |
|---|---|---|---|---|
| FY4B | 0.228 | 0.429 | **0.495** | 1.317 → 0.964 → 0.869 |
| GPM_uncal | 0.427 | 0.493 | 0.529 | 0.925 → 0.886 → 0.843 |
| GSMaP | 0.369 | 0.476 | 0.521 | 1.219 → 0.959 → 0.852 |
| IDW (gauges only) | | | r 0.482 | 0.877 |

- Correction more than doubles FY4B's raw r, and the blend beats IDW. This is the "raise its
  quality" story in one line.
- Heavy hours (gauge ≥ 8 mm/h): FY4B `residual_gwr` beats IDW by 4.7% (significant). GPM beats
  it by 3.3% and GSMaP by 6.7%.
- Dry hours: FY4B raw RMSE is 0.85 and MGWR's is 0.50. GPM raw is 0.33. False rain in dry hours is
  FY4B's main remaining weakness.
- Fold 4 of `balanced_spatial` was partial ("prediction coverage below minimum") for the joint
  methods, in every product. Coverage is 94.6–96.2%.

### 3.2 GSMaP: now one version (merged as `1fbc1ea`)

- Every script now reads `hubei_gsmap_hourly_2022_2024_full_aligned.csv`, the paper's v8
  `hourlyPrecipRate` file.
- The old `..._2022_2025_JunSep_aligned.csv` is a *different* GSMaP field, not a subset of the
  paper's file. On the shared hours r is 0.72 and no wet value is identical. Nothing reads it any
  more, so it can be deleted.
- The memory note is `gsmap-junsep-old-file.md`.

---

## 4. FY4B diagnosis and the correction chain

Taken from the current processed file,
`data/processed/study_area/hubei_fy4b_hourly_2022_2024_full_strict_navcorrected.csv`
(13,471 rows, 237 stations, values in 0.025 steps):

| Problem | Evidence | Fix |
|---|---|---|
| Download gaps | Months present: 2022 Jun–Sep; 2023 May–Sep plus 19% of Apr; 2024 every month but Feb. The original download was Jun–Sep only | Redownload the missing months (§6) |
| Strict 4-of-4 rule | An hour with any one of its four 15-min files missing becomes NaN. The 14:15 UTC file never exists in 2022–2024, so the 22–23 Beijing hour (label 23:00) is **always** missing | Keep the strict rule (user's decision, 2026-10-30) and close the gaps at the source: download the missing 15-min files (list in `output/fy4b_download/`). **Never** fill from the native one-hour files: they are gauge-corrected (§4.1) |
| Ceiling at 20 mm/h | 1,430 station-hours sit exactly at 20.0 and none are above. The product itself declares `valid_range = [0, 20]` | Treat 20 as censored (≥ 20), not as a true value |
| Location error | Point r is 0.23 but areal r is 0.43; POD 0.34, FAR 0.70 | Parallax correction (subpoints 133°E / 105°E) plus a neighbourhood mean |
| Seasonal bias and timing | SON +132%, DJF −86%; the JJA diurnal peak comes 2 h late | Per-season (and per-regime) rain threshold plus quantile mapping, fitted per fold on the training gauges only, inside the benchmark. Fitting on all 237 gauges would move r by only about 0.001 (measured), but per-fold settles the reviewer question |
| False rain in dry hours | What is left after the fixes above | The main target after re-extraction |

- **Already checked and rejected:** a side agent claimed the 30 mm/h filter in
  `extract_precipitation` cuts heavy rain. It cannot: the product tops out at 20.0, so the filter
  never fires.
- The filter does also drop DQF ≥ 2, and heavy hours show no extra NaN (0.0% vs 1.9% elsewhere).
### 4.1 What the raw files say (read 2026-10-30, step 1 done)

- **15-minute files:** satellite-only rain rate in mm/h, `valid_range = [0, 20]`. Each one is a
  single scan of about 13 minutes.
- **The 1-hour, 3-hour and 24-hour files are a different product:** "corrected with raingauge in
  China", unit mm, ranges 0–50, 0–100 and 0–250. Never mix them into the satellite series: they
  may contain our own validation gauges.
- **Fill values:** the attribute is named `FillValue`, not CF `_FillValue`, so NCDatasets does
  not mask 65533–65535. In practice the `> 30` filter is what removes them; values in (20, 30]
  cannot occur.
- **DQF:** 0 good, 1 conditionally usable, 2 out of range, 3 no value. Dropping DQF ≥ 2 is right.
- **Three regimes, not one series:**
  - Jun 2022–May 2023: old software V1.0 (2017-09-12), satellite at 133°E.
  - Jun 2023–Jan 2024: software V1.0.0 (2022-04-20), still at 133°E.
  - Mar–Dec 2024: V1.0.0 at 105°E. The viewing zenith over Hubei is 37–40°, against 44–45° from
    133°E. For a 10 km cloud top the parallax shift is about 8 km toward N/NNE instead of about
    10 km toward NW: a bit smaller, and in a different direction.

  Bias and parallax corrections should therefore be fitted per regime.
- **Coverage by 15-minute completeness, 2022–2024:** 13,474 hours have all four files, 2,907 have
  three, 202 have two, and 4,435 have one. The one-file hours are the months where only the :00
  file exists: 2022 Oct–Dec and 2023 Oct–Dec.

---

## 5. Getting the raw data onto the Mac with the T9

**Disk facts (2026-10-30):**
- T9: **exFAT**, 1 TB, 931 GB free, mounted at `/Volumes/T9`.
- Mac: 105 GB free. `data/` is 7 GB and `output/` 30 GB.

**Where the reader looks.** `src/sources/FY4BPreprocessing.jl` `walkdir`s `data/FY4B`, so any
folder layout below it works (year folders, batch folders, anything). Files are found by name:

```
FY4B-_AGRI--_N_DISK_1330E_L2-_QPE-_MULT_NOM_<YYYYMMDDhhmmss>_<YYYYMMDDhhmmss>_4000M_V0001.NC
```

- The disk tag is `1330E` or `1050E`; the times are UTC.
- **Bring the native one-hour files too.** Their start is `hh0000` and their end is `hh5959`. The
  "fill incomplete hours" fix needs them.

**⚠ AppleDouble trap (found in the code, not hit yet).**
- When macOS writes to an exFAT disk, it leaves a `._<name>` file next to every file that has
  extended attributes.
- A `._FY4B-…NC` file passes the reader's name filter.
- `find_nc_files` (`src/sources/FY4BPreprocessing.jl:353`) sorts each duplicate group and keeps
  `first(paths)`. `._` sorts before `F`, so **the junk file wins and the real file is dropped**: the
  hour then counts as unreadable or incomplete, without any error.
- So copy without them, then check:

```sh
rsync -a --progress --exclude='._*' /Volumes/T9/FY4B/ ~/GGG/GGGG/data/FY4B/
find ~/GGG/GGGG/data/FY4B -name '._*' | wc -l      # must print 0
```

- Files written onto the T9 from Windows don't get `._` files, but check anyway.
- `dot_clean -m /Volumes/T9/FY4B` cleans the disk side.

**If the raw set is bigger than about 90 GB,** leave it on the T9 and link it instead of copying:
`ln -s /Volumes/T9/FY4B ~/GGG/GGGG/data/FY4B`. `walkdir` reads through a symlinked root. The T9
must then be plugged in for every FY4B step.

**Optional:** back up `data/` to the T9. It is gitignored, so it is not on GitHub:
`rsync -a --exclude='._*' ~/GGG/GGGG/data/ /Volumes/T9/GGGG_data_backup/`.

**Data safety (CLAUDE.md):**
- Never `git clean -x` and never `git stash --all`; both reach `data/`.
- Claude can't Edit/Write inside `data/`. Copying the raw files in is the user's step, or a
  shell step the user approves.

---

## 6. What to tell the next session

1. **Which months your copy covers.** The gaps to fill:
   - 2022: Jan–May (check whether they exist at all) and Oct–Dec
   - 2023: Jan–Apr (April is only 19% present) and Oct–Dec
   - 2024: Feb
   - 2025 Jun–Sep: only needed to rebuild the missing Jun–Sep file (§8), not for the paper's
     2022–2024 window
2. **Whether the native one-hour files are included.**
3. **The total size**, to choose between copying and linking (§5).

---

## 7. Next steps once the data is in

1. ~~Read the NetCDF attributes~~ Done 2026-10-30, see §4.1. The raw data is on the T9 and linked
   as `data/FY4B`.
2. Check coverage: count files per month and disk tag, and list the months that are still missing.
3. ~~Download, re-extract, per-fold fixes~~ Done 2026-10-31, see §0b. (Plan text kept:) Download the missing 15-minute files (`output/fy4b_download/`). Then re-extract with the
   strict 4-of-4 rule, censoring at 20, and parallax plus neighbourhood. The seasonal threshold and
   quantile map are fitted per fold inside the benchmark. The full plan is
   `~/.claude/plans/mutable-dazzling-bengio.md`. The reusable parts go in `src/sources/`; the run goes
   in a `scripts/prepare_fy4b_*.jl`; results go to `output/` or `data/processed/` per CLAUDE.md.
4. ~~Rerun the with-FY4B benchmark (§3.1) on the new file~~ Done 2026-10-31, see §0b; adaptive map and QM-fed fusion in §0a.
5. Paper:
   - Add a `[[step]]` for the with-FY4B run in `paper/sources.toml`, plus `[[table]]` blocks
     for every cited number.
   - Rewrite the scope note at `paper_plan.typ:174`.
   - Run `/paper-check`.

Workflow: plan first, one branch per change, merge and push only when asked.

---

## 8. Known open problem (it existed before this work)

- `data/processed/study_area/hubei_fy4b_hourly_2022_2025_JunSep_strict_navcorrected.csv` is
  missing.
- Without it, these fail: the four `run_*_variable_selection.jl` scripts,
  `audit_mger_inputs.jl` without `--smoke-202206`, `run_mger_full.jl` without `--full-year`, and
  both `verify_mger_*.jl`.
- Rebuilding FY4B fixes this. That needs Jun–Sep 2022–2025; the check counts in the verify scripts
  are for that window.

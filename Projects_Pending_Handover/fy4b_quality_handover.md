# Handover: FY4B quality work (raw data coming via the Samsung T9)

**Date:** 2026-11-01 · **State:** the FY4B work is merged (`805bad6`); two-stage counting tables
screened on raw FY4B and stopped there by the user's decision (§0), so `FY4B_QM` stays the FY4B
quality step; a blend screen (below) found no band key that clears the decision rule · **Next
step:** the user's call on the blend screen, then the paper framing (§0a) and a paper step

---

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
  - **The blend leans on ADW**, because MGWR alone loses to ADW for every product.
  - **One subtlety:** the shift, QM and strength are fitted on all of a fold's training stations,
    including the inner hold-out group that κ and λ are tuned on.
- **Screen:** `scripts/screen_fy4b_blend_headroom.jl` → `output/fy4b_blend_headroom/` (log
  `output/fy4b_blend_headroom.log`).
  - Run at `4d44e59` from a clean tree; it takes 0.6 min.
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
  - The significant −0.91% needs the GPM and GSMaP bands, which supply more than half of it.
  - So never write "corrected FY4B beats gauges alone" on the strength of `blend_agrenv_mgwr`.
  - The QM's own effect is clearest in the own-anchor blend: −1.2% against plain FY4B.

**Where FY4B_QM's blend loses to GPM_uncal's** (balanced, `gap_by_band.csv`, share of the SSE
gap):
- FY4B_QM wet and gauge dry, i.e. false rain: 51%.
- FY4B_QM dry and gauge wet, i.e. misses: 20%.
- Both wet: 16%. Both dry: 13%.
- Band 0, where no product is wet: 15% net. Its gauge-dry cells alone give 26%: through the
  residual surface, FY4B's training-gauge errors reach cells no product calls wet.
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
- **K3 is the only real gain, and it is small.**
  - Band 0's weight is 0.4–0.6 (balanced) and 0.7–0.8 (random). In the random split K3 gives
    −0.46% (−0.75..−0.16), as context only.
  - It is not FY4B-specific: it blends the all-dry cells of any product.
- **K1 and K2 add nothing**, so FY4B's own wet call and its month quality class don't help the
  weight.
- **K4 shows the bands need GPM and GSMaP.**
- **K0, the protocol control,** scores 0.7356, against the benchmark's own 0.7391. That gap mixes
  the leaky cross-fit with the tuning. It is not a tuning-bias estimate.
- **Outcome:** stop, per the rule. Frame the paper on the improvement over ADW, about 0.1% → 0.9%:
  - more than half of it comes from the other products' bands;
  - FY4B's own part is the own-anchor row, whose interval spans 0.
- **Open for the user:** K3 as a general "blend band 0" option. It is below the 0.3% bar, applies
  to every product, and would need a nested benchmark run before it could be cited.

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

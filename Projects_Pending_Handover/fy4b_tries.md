# Register: every FY4B fix tried on the benchmark's gauges and hours

**Why this exists.** Every fix below was judged on the same 237 gauges and the same hours, so each new
try is one more chance for a fix that does nothing to clear its bar by luck. A side agent raised this
on 2026-11-09, before the last try. Fixing the rule before each run stops tuning within one test, but
not this. So:
- the list is kept whole, for the paper;
- the last try's bar shares the usual 5% over the tries counted here (`N_TRIES` in
  `scripts/compare_fy4b_fusion_same_hours.jl`).

**What counts as a try.** One row per test that had its own run or comparison script. A row can hold
more than one variant when one test judged them together. Diagnoses that built no corrected series are
listed separately, below the table, and are not counted.

Details are in `fy4b_quality_handover.md`, under the section named in the last column.

**In the paper:** Draft Table S1 (`paper_plan.typ`, supplementary material) reads `paper/fy4b_tries.csv`,
the words-only version of this table; its step and table in `paper/sources.toml` are both `fy4b_tries`.
When a try is added, add it in both places and re-export with `scripts/export_paper_numbers.py`.

| # | Date | Try | Judged on | Rule fixed first | Verdict | Record |
|---|---|---|---|---|---|---|
| 1 | 2026-10-30/31 | Rebuild choices: a fixed-height parallax shift (10 km tops), and a 3×3 mean against the centre pixel | Parallax: June 2022 raw r. Window: two full runs (w3, w1) | No | Parallax worse (r 0.351 → 0.343), dropped. The 3×3 mean kept (raw −0.7%, blend −0.1%) | §0c, §0b; `37e27d7`, `14d70d4` |
| 2 | 2026-10-31 | Per-fold, per-regime pixel shift plus quantile map (`FY4B_QM`) | Full run `_fy4bqm` | No | Raw −10.0% (sig); MGWR −1.2% and blend −0.25%, both intervals span 0 | §0b; `ede88b8` |
| 3 | 2026-10-31 | Map strength per month (`--fy4b-qm-adaptive`). The same run fed FY4B_QM to the fusion (`--fy4b-qm-fusion`) | Full run `_fy4bqmadapt_qmfusion` | No | Kept as FY4B_QM: raw −16.9%, MGWR −5.4%, blend −0.80% (all sig). MERGED_MEAN blend −0.42% (sig) against the fused-FY4B run | §0a; `e683517`, `2e9f0e4`, `b55a3a9` |
| 4 | 2026-11-01 | Two-stage counting tables (`SatelliteTwoStage`) | Screen on held-out raw FY4B | No; the user chose "stop after a raw screen" | r +0.043, but never reaches 8 mm/h, and context adds +0.011. Stopped | §0; `ad9bc97`, `978ef74` |
| 5 | 2026-11-01 | Blend-weight keys K1–K4 (FY4B wet call, mapped month, band 0, FY4B-only bands) | Screen, cross-fitted out-of-fold, exploratory | Yes (the blend screen's rule) | No key clears. K3 (band 0 blended) −0.17% | "Latest (2026-11-01)"; `4d44e59` |
| 6 | 2026-11-01 | Dry-band blends, band 0 toward ADW or toward the anchor (`--blend-dry-band`) | Full nested run | Yes | Neither clears: toward ADW −0.20%; toward the anchor −0.71% but bias −0.012 mm/h | "Earlier (2026-11-01): nested run"; `8232154`, `fd59d5f` |
| 7 | 2026-11-01 | Whole-hour time shift | Screen | Yes | No: the best shift is 0 in JJA and SON, and held-out Δr is 0 | "two cheap screens"; `1dc0262` |
| 8 | 2026-11-01 | FY4B_QM in the fusion, first check (`compare_fy4b_in_fusion.jl`, against the then paper run) | The two existing runs | Yes | "Confounded" (hours and seeds differ). At face value no agreement blend clears | "does FY4B help the fusion?"; `80287f2` |
| 9 | 2026-11-02 | Cloud-type filter from FY4B's CLT product | July 2023 pilot | Yes | Fails: ΔCSI_0.1 +0.0056 against the +0.01 bar | "the CLT cloud-type filter"; `1a61f2a` |
| 10 | 2026-11-06 | FY4B_QM in against out of the fusion, same hours and seeds (`--fy4b-out-of-fusion`) | New full run against the qmfusion run | Yes | Keep FY4B out. MERGED_MEAN −0.45% but heavy rain worse; MERGED_OLS −0.22%, under the bar. (On 2026-11-08 the user put plain FY4B into the fusion for the three-product comparison.) | "FY4B in the fusion, on the same hours and seeds"; `d273b46`, `e59f367` |
| 11 | 2026-11-06 / 11-08 | Near-surface humidity gate on FY4B_QM (`HumidityGate`, `--fy4b-qm-rh-gate`) | Screen, then a full run | Yes, both | Screen GO for the gate that spares ≥ 8 mm/h. Full run: raw −1.7%, but blend +0.03%. Stop | "humidity gate passes its screen", "humidity gate in the benchmark"; `3f1c3cc`, `d7a8cc7`, `423b2f4` |
| 12 | 2026-11-07 | AGRI L1 channels (BT10.8, WV − IR, split window) to tell false rain from real | Six scans, read inline | Yes (the reading) | Stop, a weak "no" from a small pilot | "a first look at two AGRI L1 scans"; nothing committed |
| 13 | 2026-11-09 | FY4B_QM in place of plain FY4B in the merges and bands | Scratch comparison: qmfusion run against the paper run | Yes (the script's rule) | No agreement blend clears. MERGED_MEAN −0.42%, but CSI_8 and \|FBI_8 − 1\| worse | "why FY4B_QM stays out of the merges and the bands"; scratch in `output/fy4b_qm_vs_plain_in_fusion_scratch/` |
| 14 | 2026-11-09 | FY4B_QM fitted out of time, in place of plain FY4B (`--fy4b-qm-crossed --fy4b-qm-fusion`) | New full run against the paper run | Yes ("no worse") | Fails \|FBI_8 − 1\| on GPM_uncal and MERGED_MEAN; keep plain FY4B | "FY4B_QM fitted out of time"; `7368f66`, `2e9d832` |
| 15 | 2026-11-09 | **Confirmed FY4B** in place of plain FY4B in the merges and bands: FY4B kept only where GPM_uncal or GSMaP also calls the cell wet (`--fy4b-confirm-fusion`). Gauge-free, fits nothing | New full run against the paper run | Yes, stricter (below) | Fails. MERGED_MEAN blend −0.04% (99.67%: −0.22..+0.14); GSMaP's CSI_8 and MERGED_MEAN's \|FBI_8 − 1\| worse; R1 +0.05%; random +0.09% (sig). Keep plain FY4B. Its raw series is better on its own (RMSE −5.1%, FAR 0.71 → 0.57), and the raw merge gains −0.97%, gone after MGWR | "the last try, confirmed FY4B"; `e0b208f`, `b044b3f`, run at `15db848` |

**Try 15's rule**, chosen by the user and fixed before its run. It is judged on balanced_spatial, with
7-day blocks:
1. **Primary:** MERGED_MEAN's agreement blend improves RMSE by 0.3% or more, and the upper end of its
   1 − 0.05/15 ≈ 99.67% interval is below 0.
2. **No worse:** all four agreement blends pass try 14's checks at 95%:
   - the RMSE upper end is below +0.3%;
   - CSI_8 is not lower;
   - |FBI_8 − 1| is not worse;
   - |bias| rises by at most 0.005 mm/h.
3. **Replication:** MERGED_MEAN's agreement blend is lower in each of FY4B's periods R1, R2 and R3.
   In the random split it is not significantly worse. That check was first "lower in random"; it was
   changed before the run, because the random split has never moved more than 0.12%.

**Not counted.** These are diagnoses; none built a corrected series:
- 105°E against 133°E;
- the 2023 software update;
- the error budget;
- the time-held-out check of FY4B_QM;
- the week-block recheck of the paper's intervals.

**Closed (2026-11-10).** Try 15 was agreed with the user as the last, and it fails, so the line of FY4B fixes
for the fusion ends here. Plain FY4B stays in the merges and the bands.

**Parked:** confirmed FY4B in the merges only, the bands keeping plain FY4B. It is in `ROADMAP.md` as
exploratory. It would be try 16, with the bar 1 − 0.05/16, and it is not independent of try 15: its
merged anchor is try 15's, whose gain is gone after MGWR.

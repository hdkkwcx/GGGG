// FY4B's path through the interpolation benchmark, in brief: the stages in order, with the name of
// each method and no formula.
// Compile:  typst compile NOWAY/fy4b_benchmark_workflow.typ   (it reads no data file)
//
// The short companion of the methods page (scripts/plot_fy4b_methods.py ->
// output/fy4b_process_trace/fy4b_methods.html), which gives the same stages with their formulas, each
// checked against the code by test/test-fy4b-methods-formulas.jl. The stage order follows that page.
// Method and mask lists are BENCHMARK_RUNS and MASK_METHODS in
// src/benchmark/InterpolationBenchmarkConfig.jl; the covariate groups are JOINT_GROUP_ORDER in
// src/models/CovariateGroups.jl; the run's settings are in its benchmark_scope.csv; the switches are
// described in CLAUDE.md ("Interpolation benchmark invariants").

#set document(
  title: "FY4B in the benchmark: the workflow in brief",
  date: datetime(year: 2026, month: 11, day: 16),
)
#set page(
  paper: "a4",
  margin: (x: 2.1cm, top: 1.9cm, bottom: 1.9cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 8pt, fill: luma(110))
      FY4B in the benchmark: the workflow in brief
      #h(1fr) 16 November 2026
    ]
  },
)
#set text(font: ("Libertinus Serif", "New Computer Modern"), size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.62em, spacing: 1.05em)
#set heading(numbering: none)
#show heading: set text(font: "Arial", weight: "semibold")
#show heading.where(level: 1): it => { v(0.4em); it; v(0.15em) }
#set figure(numbering: none)
#show figure.caption: set text(size: 9pt)
#show figure.where(kind: table): set figure.caption(position: top)
#show table.cell: set text(size: 9pt)
#show table: set par(justify: false)
#set list(indent: 0.6em)

#let accent = rgb("#2a78d6")
#let rule = table.hline(stroke: 0.8pt)
#let thin = table.hline(stroke: 0.4pt)

/// One stage of the flow chart: its number, its name, the methods used, and which FY4B is in it.
#let stage(number, title, body, fy4b) = block(
  width: 100%, breakable: false, above: 0pt, below: 0pt,
  inset: (x: 9pt, y: 6.5pt), radius: 4pt, fill: luma(249), stroke: 0.6pt + luma(160),
  grid(
    columns: (17pt, 1fr), column-gutter: 9pt,
    align(center + horizon, box(width: 17pt, height: 17pt, radius: 8.5pt, fill: accent,
      align(center + horizon, text(font: "Arial", size: 9pt, weight: "bold", fill: white)[#number]))),
    {
      set par(justify: false, leading: 0.55em, spacing: 0.5em)
      set text(size: 9.5pt)
      text(font: "Arial", size: 10pt, weight: "semibold")[#title]
      parbreak()
      body
      parbreak()
      text(size: 8.5pt, fill: accent)[FY4B here: #fy4b]
    },
  ),
)
/// The arrow between two stages.
#let arrow = align(center, block(above: 4pt, below: 4pt, box(width: 7pt, stack(dir: ttb,
  align(center, rect(width: 0.9pt, height: 5pt, fill: luma(120), stroke: none)),
  polygon(fill: luma(120), (0pt, 0pt), (7pt, 0pt), (3.5pt, 4.5pt))))))

#align(center)[
  #text(font: "Arial", size: 17pt, weight: "semibold")[FY4B in the benchmark: the workflow in brief]
  #v(-0.3em)
  #text(size: 11pt, fill: luma(60))[The stages in order, with the name of each method. No formulas.]
]

#v(0.9em)

#stage(0, [Input, made before the benchmark],
  [*Rebuilt hourly FY4B* at the gauges, from the satellite-only 15-minute files: *strict hour* (all four
  files) · *3 × 3 pixel mean* · *offset grid* of whole-pixel shifts, kept for stage 3.],
  [plain FY4B.])
#arrow
#stage(1, [What the run loads],
  [Six anchors, in this order: *FY4B*, *GPM (IMERG Uncal)*, *GSMaP*, *equal-weight merge*, *OLS merge*,
  *corrected FY4B*. The run uses the hours every product file has, which are FY4B's hours.],
  [plain FY4B is the first anchor. The last three anchors are built in stages 3 and 4.])
#arrow
#stage(2, [Cross-validation],
  [*Balanced spatial 5-fold* (main) · *random station 5-fold* (second). Every choice in stages 3 to 6
  is made on the fold's training gauges only. Tuning uses an *inner 5-group spatial split* of them.],
  [the same folds as every other anchor.])
#arrow
#stage(3, [Corrected FY4B, fitted inside each fold],
  [*Pixel shift* per regime (satellite position and software version) · *wet threshold* · *quantile
  mapping* per season and regime · *adaptive strength*, month by month.],
  [plain FY4B goes in, corrected FY4B (`FY4B_QM`) comes out.])
#arrow
#stage(4, [Merged anchors, fitted inside each fold],
  [*Equal-weight merge* (`MERGED_MEAN`) · *OLS merge* (`MERGED_OLS`) of FY4B, GPM and GSMaP.],
  [plain FY4B, as one of the three products merged.])
#arrow
#stage(5, [The methods, run for every anchor and fold],
  [*Raw* anchor · from the gauges alone: *IDW*, *ADW*, *TPS*, *GWR* · the anchor corrected with the
  gauges: *residual GWR*, *mixed GWR*, *MGWR*.],
  [plain FY4B and corrected FY4B each go through all of them.])
#arrow
#stage(6, [Blends with ADW, for the three correcting methods],
  [*Constant blend* · *agreement blend*, with one weight per *agreement band* (how many of the three
  products report rain, and their largest rate). MGWR with the agreement blend is the headline method.],
  [plain FY4B is one of the three votes, also when corrected FY4B is the anchor.])
#arrow
#stage(7, [Scores at the held-out gauges],
  [*Common evaluation mask* · *RMSE*, *MAE*, *bias*, *r* · *POD*, *FAR*, *CSI*, *FBI* at rain
  thresholds · by *intensity class* and by *distance to the nearest training gauge*.],
  [each of its anchors is scored under every method.])
#arrow
#stage(8, [Tests],
  [*Paired day-block bootstrap* against IDW, ADW and TPS · *Holm's correction*. For the paper the run is
  then rescored on the *matched sample*, with *7-day blocks*.],
  [corrected FY4B is also compared with plain FY4B, GPM and GSMaP, pair by pair.])

#pagebreak()

= Inside stage 5: what the correcting methods use

- *Covariates, nine groups.* From the DEM: elevation, slope, aspect. From ERA5-Land, hourly: 2 m
  temperature, 2 m dew point, 10 m wind (two components), surface pressure. And NDVI.
- *Nested covariate selection*, repeated in every fold for every anchor: *screening* of each group on
  its own (permutation test, Benjamini–Hochberg adjustment) · *collinearity check* (VIF) ·
  *spatial-variability test*, which makes each kept group local or global.
- *Local fit.* Kernels: Gaussian, exponential, bisquare, tricube, boxcar. Bandwidth: adaptive
  (nearest gauges) or fixed (distance).
- *Mixed GWR and MGWR* are fitted by *back-fitting*, MGWR with over-relaxation and one bandwidth per
  covariate group, found by *coordinate descent*.
- *Tuning.* Kernel, bandwidth and *residual shrinkage* are chosen on the *tuning hours* (the wettest
  hours and evenly spaced hours), by the error on the inner split.

= The three places FY4B appears

+ *As an anchor:* plain FY4B, corrected by the gauges in stages 5 and 6.
+ *As a corrected anchor:* `FY4B_QM`, the output of stage 3, run through the same stages 5 and 6. It
  is the last anchor in the list, and it enters neither the merges nor the agreement bands.
+ *Inside the merges and the agreement bands:* always plain FY4B, beside GPM and GSMaP.

= Switches that are off in the paper's run

#figure(
  table(
    columns: (auto, 1fr),
    align: (left, left),
    stroke: none,
    inset: (x: 5pt, y: 4pt),
    rule,
    table.header([*Switch*], [*What it changes*]),
    thin,
    [`--fy4b-qm-fusion`], [The merges and the agreement bands take corrected FY4B in place of plain FY4B.],
    [`--fy4b-out-of-fusion`], [FY4B stays an anchor, but the merges and the bands use GPM and GSMaP only.],
    [`--fy4b-confirm-fusion`], [The merges and the bands take a FY4B kept only where GPM or GSMaP also reports rain.],
    [`--fy4b-qm-crossed`], [Stage 3 is also fitted out of time: each 7-day block takes the correction fitted on the other blocks.],
    [`--fy4b-qm-rh-gate`], [Adds one more anchor: corrected FY4B set to zero where ERA5-Land humidity is too low.],
    [`--fy4b-twostage`], [Adds one more anchor: FY4B's two-stage series (the chance of rain, then the wet amount).],
    rule,
  ),
)

= The run, and where to read more

The paper's benchmark step is this command:

#block(width: 100%, inset: (x: 8pt, y: 6pt), radius: 3pt, fill: luma(246))[
  #set text(size: 8.5pt)
  #set par(justify: false)
  `julia -t 4 --project=. scripts/run_interpolation_benchmark.jl full --nested-covariates --satellite-wet-blend --blend-axis agreement_envelope --fused-anchor --with-fy4b --fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv --fy4b-qm --fy4b-qm-adaptive`
]

#set par(justify: false)
- *The formulas* of every stage: the methods page, drawn by `scripts/plot_fy4b_methods.py`.
- *The counts* after every stage: `scripts/trace_fy4b_process.jl`, drawn by
  `scripts/plot_fy4b_process.py`.
- *The results:* `NOWAY/why_fy4b_scored_worst.typ` and §3.3–§3.4 of `paper/paper_plan.typ`.

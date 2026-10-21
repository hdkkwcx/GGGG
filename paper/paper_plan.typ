// Research plan: gauge, GPM and GSMaP precipitation in a sparse-gauge mountain area.
// Compile from the repository root:  typst compile paper/paper_plan.typ
// Figures are copies of the `--gpm-gsmap` outputs of scripts/plot_heavy_rain_events.py,
// scripts/plot_satellite_temporal_evaluation.py and scripts/plot_gauge_satellite_fusion.py
// (the fusion figures from its FY4B-free IMERG Uncal report, output/gauge_satellite_fusion_nofy4b_gpmuncal).
// Numbers are read from the CSVs beside those figures in output/.

#set document(
  title: "Research plan: satellite precipitation and gauge fusion (GPM, GSMaP)",
  date: datetime(year: 2026, month: 10, day: 19),
)
#set page(
  paper: "a4",
  margin: (x: 2.1cm, top: 2.3cm, bottom: 2.2cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 8pt, fill: luma(110))
      Research plan — GPM and GSMaP
      #h(1fr) 19 October 2026
    ]
  },
)
#set text(font: ("Libertinus Serif", "New Computer Modern"), size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.62em, spacing: 1.05em)
#set heading(numbering: "1.1")
#show heading: set text(font: ("Segoe UI", "Arial"), weight: "semibold")
#show heading.where(level: 1): it => { v(0.6em); it; v(0.25em) }
// Captions carry their planned paper numbers ("Draft Fig. 3"), so Typst adds none of its own.
#set figure(numbering: none)
#show figure.caption: set text(size: 9pt)
#show figure.where(kind: table): set figure.caption(position: top)
#show table.cell: set text(size: 8.5pt)
#show table: set par(justify: false)
#set list(indent: 0.6em)

#let note(body) = block(width: 100%, inset: (x: 10pt, y: 8pt), radius: 3pt, fill: luma(246),
  stroke: (left: 2.5pt + rgb("#2a78d6")), body)
#let rule = table.hline(stroke: 0.8pt)
#let thin = table.hline(stroke: 0.4pt)
// A literature slot still to be filled with a paper actually read.
#let refslot = text(fill: rgb("#c0392b"))[\[ref\]]
// A figure planned for the paper but not shown in this plan.
#let planned(body) = block(inset: (left: 8pt), text(size: 9pt, fill: luma(80))[_Planned:_ #body])

#align(center)[
  #text(size: 17pt, weight: "semibold", font: ("Segoe UI", "Arial"))[Research plan]
  #v(0.2em)
  #text(size: 11.5pt, fill: luma(60))[Can GPM and GSMaP stand in for gauges where the network is sparse, and what does fusing them with the gauges gain?]
]

#note[
  *Main line and status.* Three research questions, one results section each (§3).
  #v(0.3em)
  #table(
    columns: (auto, 1.3fr, 2fr, auto, 1.4fr),
    align: (left, left, left, left, left),
    stroke: none,
    rule,
    table.header([*RQ*], [*Question*], [*Answer so far*], [*Status*], [*Remaining*]),
    thin,
    [1], [Does the satellite place heavy rain correctly in space?],
    [Widespread rain yes (r ≈ 0.55, both products); localized storm cores no. GSMaP far too wet.],
    [Answered], [Rerun on IMERG Uncal],
    [2], [Does it follow the hourly evolution at each intensity?],
    [Rain is detected and timed well; peak magnitude is lost, by 60–85% at ≥ 8 mm h⁻¹.],
    [Answered], [Landform panel in Fig. 4; rerun on IMERG Uncal],
    [3], [What does gauge fusion gain, with which method and product?],
    [MGWR + blend beats raw GPM by 10% and gauges alone by 3%; gain largest on heavy hours and far from gauges. GPM is the better anchor, with or without gauge calibration; fused GSMaP does not beat gauges alone.],
    [Answered, gaps], [Choose headline method],
    rule,
  )
]

= Introduction

Flood forecasting in mountain catchments needs hourly rainfall, and above all the moderate and heavy
rain that drives runoff. Rain gauges measure it directly but are sparse in the mountains, so large
parts of a catchment sit tens of kilometres from the nearest gauge. Satellite products such as GPM
IMERG and GSMaP cover every pixel every hour, but their errors depend on rain type and intensity.

Most evaluations of these products work at daily or monthly scale and report pooled scores.
Pooled scores do not say whether a product puts a heavy-rain event in the right place, whether it
follows the hour-by-hour evolution of rain of a given intensity, or where its error comes from. Gauge–satellite
fusion is widely used to correct the products, but it is usually judged by average skill. What fusion
adds *far from the gauges* and *on heavy hours*, which is exactly the sparse-network flood case, is less
clear.

*Objective.* To establish how far GPM and GSMaP can stand in for gauges in a sparse mountain network,
and what fusing them with the gauges gains. Three questions:

- *RQ1 — Space.* Does the satellite field reproduce *where* heavy rain falls?
- *RQ2 — Time.* Does it reproduce *how rain evolves in time* at each intensity?
- *RQ3 — Fusion.* What does gauge–satellite fusion gain, which fusion method and which anchor product
  are best, and does fusion fix the moderate and heavy rain that floods depend on?

*Scope.* Two products, GPM IMERG and GSMaP. FY4B is left out: on its own hours it scores worst even
when no product is gauge-calibrated (hourly r 0.23 vs 0.37–0.43, areal-mean r 0.43 vs 0.66–0.83,
wet-hour detection 0.34 vs 0.64–0.66). It enters no part of the analysis: the merged product and the
agreement blend in RQ3 use GPM and GSMaP only.

= Data and Methods

== Study area and data

237 hourly gauges in the Shiyan study area, 2022–2024 (values quantised to 0.5 mm). All sit in hills
and mountains: by local relief from the Copernicus GLO-30 DEM, 129 are below 500 m, 97 at 500–1000 m
and 11 at 1000 m or more. GPM IMERG V07 and GSMaP v8 operational (`hourlyPrecipRate`), both hourly
at 0.1°, are sampled at the gauge pixels. The primary GPM product is IMERG Uncal, so that neither
product carries gauge information: IMERG Final is calibrated monthly to GPCC gauge totals, GSMaP
operational has no gauge correction. IMERG Final is kept as a sensitivity check. RQ1 and RQ2 are so
far evaluated on IMERG Final and are still to be rerun on IMERG Uncal. Days run 08–08 Beijing time.

#planned[Fig. 1 — study area: relief, gauges, landform regions, location inset, river basins.]

== Evaluating the satellite products (RQ1, RQ2)

*Space (RQ1).* Heavy-rain days are 08–08 days on which more than three gauges exceed 50 mm. Of 56 such
days, six are selected, one per rain system: four widespread, two localized. At each gauge the satellite
24 h total is compared with the gauge total by spatial correlation, spatial CV ratio, rain-centre
shift, mean bias, RMSE and CSI at 50 mm.

*Time (RQ2).* Only gauge-wet hours are scored, in five intensity classes (0.1–2, 2–4, 4–8, 8–20 and
≥ 20 mm h⁻¹), stratified by season and landform region, with day-block bootstrap intervals. Rain
events are runs of wet hours separated by at least 3 dry hours and are scored for detection, timing and
volume. Areal-mean series are compared at hourly and daily scale.

== Gauge–satellite fusion and its validation (RQ3)

Gauge-only interpolators (IDW, ADW, TPS, GWR) are the baseline. Fusion methods use a satellite product
as the anchor and correct it with the gauges (residual GWR, mixed GWR, MGWR), optionally blended with
ADW where the satellite reports rain. The agreement blend chooses the ADW weight separately for each
band of how many of the two products call the cell wet (≥ 0.1 mm h⁻¹) and the larger of their two
rates (edges 0.1, 2.5 and 8 mm h⁻¹). Anchors are GPM (IMERG Uncal), GSMaP and an OLS merge of the two;
the whole fusion is repeated with IMERG Final in place of IMERG Uncal as a sensitivity run. Every
method is scored at gauges held out of the fit, by balanced spatial 5-fold cross-validation over all
26,304 hours of 2022–2024, using RMSE, correlation and CSI, overall and by intensity and distance to
the nearest training gauge.

= Research content and results

== RQ1 — Does the satellite place heavy rain correctly in space?

*Answer so far.*
- *Placement is similar for both products, and fair.* Mean spatial r over the six events is 0.55 for
  GPM and GSMaP alike, 0.6–0.8 on four events. Both fail on 2022-06-26 (r ≤ 0.12), although each 8 h
  phase of that day scores r = 0.43–0.88: the rain moved north to south, so the gauge total is nearly
  flat and varies mostly below 10 km, while the products keep a north–south gradient the gauges do
  not show.
- *Amounts differ sharply.* GPM is nearly unbiased (MB −2.0 mm, RMSE 20.5 mm); GSMaP is wet
  (MB +18.9 mm, RMSE 39.9 mm), up to +73 mm on 2023-08-26.
- *Localized storm cores are missed by both* (CSI at 50 mm is 0 except GPM on 2024-07-30). Pixel
  averaging explains only 16–43% of the core deficit. The cores fell as short bursts (median peak
  25–27 mm h⁻¹, up to 87% of a gauge's total in three hours), which both products spread over more
  hours at about a fifth of the peak rate. GPM's partial success on 2024-07-30 comes from its monthly
  gauge calibration (×1.7 over the uncalibrated estimate).

#figure(
  image("figures/fig02_event_maps.png", width: 100%),
  caption: [Draft Fig. 2. The six heavy-rain events (columns). Top row: gauge 24 h totals (08–08 BJT).
  Lower rows: satellite minus gauge at each gauge, from the pixel sampled there; red means the satellite is
  drier, blue wetter, grey within ±10 mm. MB is the mean bias (mean of satellite − gauge over the
  event's gauges).],
)

#planned[Fig. 3 — spatial scores per event and product.]

*Still to do.* Rerun on IMERG Uncal, the paper's primary GPM product; the numbers above are IMERG Final.

== RQ2 — Does the satellite follow the hourly evolution at each intensity?

*Answer so far.*
- *Underestimation grows with intensity, in every season.* GPM relative bias runs from −11% (light) to
  −38% (4–8), −68% (8–20) and −85% (≥ 20 mm h⁻¹); GSMaP from +16% to +8%, −17%, −59% and −80%. At
  ≥ 8 mm h⁻¹ at most 13% of hours land in the right class.
- *Detection and timing are not the problem; magnitude is.* Rain events peaking at 2–4 mm h⁻¹ or more
  are detected 88–96% of the time, and rain-centre timing errors stay within ±0.5 h.
- *Winter is where both fail* (bias about −80%, only 14–18% of wet hours detected). Summer light rain
  is overestimated (GPM +39%, GSMaP +52%).
- *Terrain matters little.* Between relief below 500 m and 500–1000 m the class biases differ by at most
  8 percentage points, far less than between classes or seasons.
- *Areal means are much better than points.* Hourly r is 0.83 for GPM and 0.66 for GSMaP (daily 0.91
  and 0.80); overall bias +6% and +33%.

#figure(
  image("figures/fig04a_intensity_season.png", width: 100%),
  caption: [Draft Fig. 4. Hourly skill by gauge intensity class and season.],
)

#planned[Fig. 5 — rain events: detection, timing and volume by class.]

*Still to do.*
- Add the landform-region panel to Fig. 4 (the result is known; only the figure is missing).
- Rerun on IMERG Uncal, the paper's primary GPM product; the numbers above are IMERG Final.

== RQ3 — What does gauge fusion gain, with which method and which product?

*Answer so far.*
- *Fusion beats both of its inputs on GPM, not on GSMaP.* MGWR followed by a blend with ADW where the
  satellites report rain cuts RMSE by 10% from raw GPM and 32% from raw GSMaP. Against gauge-only ADW it
  gains 2.8% on GPM and 3.1% on the merged anchor, both significant, but only 1.3% on GSMaP, which is
  not (p = 0.10). Every blended variant on GPM and on the merged anchor beats ADW significantly; over
  all hours none on GSMaP does.
- *The gain is where floods need it.* On heavy hours (≥ 8 mm h⁻¹) MGWR alone improves on ADW by 5.5%
  (GPM) and 8.0% (GSMaP), the agreement blend by 2.8% and 4.0%, all significant. Unblended regression
  *hurts* dry, light and moderate hours (MGWR on GPM −22%, −25% and −33% against ADW), which is why the
  blend is needed. At 2.5–8 mm h⁻¹ the blend gains nothing significant (+0.6% on GPM).
- *The gain grows with distance from the network.* On GPM the improvement over ADW is 0.6% within
  20 km of a training gauge, 2.7% at 20–50 km and 8.2% at 50–100 km (GSMaP 0.3%, 1.8%, 1.9%).
- *Best method.* Among gauge-only interpolators ADW, with IDW indistinguishable (0.695 vs 0.697 mm h⁻¹);
  TPS and GWR are 5–6% worse. Best fusion: MGWR + blend; MGWR alone only for heavy hours.
- *Best product.* GPM, with neither product gauge-calibrated: lower raw error over all 2022–2024 hours
  (RMSE 0.75 vs 1.02 mm h⁻¹, r 0.43 vs 0.38). GPM's near-zero event bias and better areal timeline
  (RQ1, RQ2) are so far measured on IMERG Final. Once fused, the two anchors differ by about 1.6%
  (0.676 vs 0.687), but only fused GPM beats ADW significantly. Merging gives a better raw field
  (RMSE 0.718) but adds almost nothing after fusion (0.674 vs 0.676); its OLS weights are 0.49 for GPM
  and 0.15 for GSMaP, plus a 0.04 mm h⁻¹ intercept.
- *IMERG Final gives the same fused result.* Gauge-calibrated IMERG Final has higher raw error than
  IMERG Uncal (RMSE 0.781 vs 0.750 mm h⁻¹), but after fusion the two are level (0.676 both; +2.7% vs
  +2.8% over ADW). On Final, MGWR alone loses more to ADW (−6.4%, against a non-significant −2.7% on
  Uncal) and gains more on heavy hours (6.8% vs 5.5%). The OLS merge puts less weight on Final than on
  Uncal (0.41 vs 0.49).

#figure(
  image("figures/fig06b_fusion_methods.png", width: 100%),
  caption: [Draft Fig. 6. RMSE improvement of each fusion method over gauge-only ADW, per anchor, for
  all hours and heavy hours. Rows are grouped by how the satellite enters: direct correction, blend with
  ADW, agreement blend. Hollow cells are not significant.],
)

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right),
    stroke: none,
    rule,
    table.header([*Anchor*], [*Method*], [*RMSE*], [*r*], [*CSI 2*], [*vs raw*], [*vs ADW*]),
    thin,
    [—], [ADW (gauges only)], [0.695], [0.505], [0.318], [—], [0.0%],
    [—], [IDW (gauges only)], [0.697], [0.501], [0.317], [—], [−0.4%],
    [—], [TPS (gauges only)], [0.732], [0.468], [0.314], [—], [−5.4%],
    [—], [GWR (gauges only)], [0.738], [0.438], [0.284], [—], [−6.2%],
    thin,
    [GPM], [raw], [0.750], [0.426], [0.250], [0.0%], [−7.9%],
    [GPM], [MGWR + agreement blend], [0.676], [0.540], [0.347], [+9.9%], [+2.8%],
    [GSMaP], [raw], [1.017], [0.382], [0.257], [0.0%], [−46.0%],
    [GSMaP], [MGWR + agreement blend], [0.687], [0.540], [0.354], [+32.4%], [+1.3%],
    [Merged (OLS)], [raw], [0.718], [0.447], [0.224], [0.0%], [−3.3%],
    [Merged (OLS)], [MGWR + agreement blend], [0.674], [0.544], [0.335], [+6.2%], [+3.1%],
    thin,
    [GPM (Final)], [raw], [0.781], [0.425], [0.264], [0.0%], [−12.5%],
    [GPM (Final)], [MGWR + agreement blend], [0.676], [0.540], [0.353], [+13.5%], [+2.7%],
    rule,
  ),
  caption: [Draft Table 3. Held-out gauges, balanced spatial CV over all 26,304 hours of 2022–2024;
  RMSE in mm h⁻¹; improvements are RMSE reductions; CSI at 2 mm h⁻¹. GPM is IMERG Uncal. No row uses
  FY4B: the merge and the agreement blend are built from GPM and GSMaP only. The GPM (Final) rows are a
  sensitivity run with IMERG Final in IMERG Uncal's place, otherwise identical. Gauge-only rows are on
  GPM's evaluation mask.],
)

#planned[Fig. 7 — where fusion helps: by intensity, by distance to the nearest gauge, and in CSI.
Fig. 8 — raw vs fused skill per anchor, and the merge weights.]

*Still to do.*
- Decide whether the newer stacked combination of MGWR and ADW replaces the agreement blend as the
  headline fusion method.

= Discussion

== Uncertainties

- *Point versus pixel.* A gauge samples a point, the products a 0.1° pixel. Where rain varies below
  10 km (2022-06-26, the localized cores) part of the "error" is representativeness, not product skill;
  pixel averaging explains 16–43% of the core deficit, and the rest comes from the products spreading
  short bursts over too many hours.
- *Gauge resolution.* The 0.5 mm quantisation limits what can be said about the lightest class.
- *Independence of GPM.* The fusion uses IMERG Uncal, which carries no gauge information, and IMERG
  Final gives the same fused result. RQ1 and RQ2 are still on IMERG Final, which is calibrated to GPCC
  gauges; GPM's partial capture of the 2024-07-30 core depends on that calibration.
- *Small samples.* Only 11 gauges lie at relief ≥ 1000 m, too few to interpret; heavy hours and the six
  events are few, so their intervals are wide.
- *Residual heavy-rain deficit.* Fusion narrows but does not close the underestimation at ≥ 20 mm h⁻¹.
- *Generality.* One region and three years; the transfer to other climates and denser or sparser
  networks is untested.

== How others have studied these problems

- Validation of IMERG and GSMaP in complex terrain, and by rain intensity #refslot. Most report pooled or
  daily scores; here the error is decomposed into placement (RQ1) and hourly magnitude by class (RQ2).
- Satellite–gauge merging by residual interpolation (kriging, IDW) #refslot, geographically weighted
  regression #refslot and machine learning #refslot. This study compares them on held-out gauges and asks
  where, not only whether, the satellite adds skill.
- Multiscale GWR, which lets each covariate have its own bandwidth #refslot. Applied here as the
  correction step of the fusion.
- Spatial rather than random cross-validation for interpolation and merging #refslot. Random splits
  leave a held-out gauge next to a training one and overstate skill far from the network; the balanced
  spatial folds used here measure that case directly.
- Satellite rainfall as input to flood and rainfall–runoff models #refslot, which motivates the outlook
  below.

== Outlook: flood forecasting

RQ2 shows that raw satellite rain loses about 60–85% of the volume at 8 mm h⁻¹ and above, the input
flood models are most sensitive to, so raw products should not drive them. Fusion recovers part of
this (RQ3). The planned extension: aggregate raw, gauge-only and fused fields to basin-mean rainfall
for the gauged catchments during the six RQ1 events and the wettest seasons; correlate basin rainfall
(and antecedent totals) with observed discharge, raw against fused; and only if the correlation
improves clearly, drive a simple rainfall–runoff model as a demonstration. This needs discharge data,
which the project does not yet have.

= Work remaining before drafting

- Fig. 1: study-area map.
- Fig. 4: add the landform-region panel (RQ2).
- Rerun RQ1 and RQ2 on IMERG Uncal. `scripts/run_heavy_rain_event_evaluation.jl` and
  `scripts/run_satellite_temporal_evaluation.jl` read only IMERG Final today.
- Choose the headline fusion method: agreement blend or stacked combination (RQ3).
- Fill the #refslot slots in §4.2 with papers actually read.
- Move the diurnal-cycle results to supplementary material.
- Obtain hydrological-station discharge data (outlook only; not needed for RQ1–RQ3).

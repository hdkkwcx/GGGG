// Research plan: gauge, FY4B, GPM and GSMaP precipitation in a sparse-gauge mountain area.
// Compile from the repository root:  typst compile paper/paper_plan.typ
// Every comparison uses one matched sample: the gauge-hours from June 2022 at which the gauge and
// FY4B (the rebuilt series of §2.5, output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv), GPM (IMERG Uncal)
// and GSMaP all have a value (2026-11-08; before that RQ1-RQ3 covered GPM and GSMaP over 2022-2024).
// Figures are copies of the outputs of scripts/plot_heavy_rain_events.py,
// scripts/plot_satellite_temporal_evaluation.py and scripts/plot_gauge_satellite_fusion.py; numbers
// are read from the CSVs beside them in output/. Intervals and p-values resample whole 7-day blocks.
// Fig. 1 is scripts/plot_study_area_map.py's output/study_area_map/figures/fig1_study_area.png, drawn
// straight from data/ (its plot_provenance.csv records the inputs' SHA-256).
// RQ1, RQ2 and Figs. 2-4 come from clean runs with `--fy4b-file ... --from 2022-06`:
// output/heavy_rain_events_fy4bw3h0_gpmuncal_from202206, heavy_core_diagnostics_fy4bw3h0_gpmuncal_from202206
// and satellite_temporal_evaluation_fy4bw3h0_gpmuncal_from202206_block7d (its `all_products` sample).
// Table 3, Figs. 6-7 and §3.4 come from a clean benchmark run with FY4B in the merges and bands and
// corrected FY4B appended (..._fusedanchor_fy4bw3h0_fy4bqm_fy4bqmadapt_gpmuncal, recorded in its
// benchmark_scope.csv with its Manifest.toml and a SHA-256 for every input file), rescored on the
// matched sample by scripts/run_claim_reassessment.jl and run_gauge_satellite_fusion.jl
// (`--matched-products`), compare_fy4b_quality.jl, check_fy4b_qm_time_holdout.jl and the July 2023
// cloud-type pilot screen_fy4b_clt.jl; Table 4 by run_heavy_rain_event_fy4b.jl. Each output's
// provenance.csv records its commit. §3.4's sentence on false-rain episodes rests on
// scripts/explain_fy4b_qm_time_holdout.jl at 3e45ae6 (output/fy4b_qm_time_holdout_explained), on the
// earlier FY4B run's cells; corrected FY4B's raw field is the same in both runs.
// scripts/make_paper.jl rebuilds every figure and number here in dependency order; the [[step]]
// graph in paper/sources.toml records each step's command and what it needs.
// paper/references.bib is the BibTeX export of the Zotero library's items tagged GGGG (local API,
// format=bibtex), without the `file` and `abstract` fields, and with the two book chapters typed
// @inbook and no series `number`: typst's APA style repeats the book title for @incollection.

#set document(
  title: "Research plan: satellite precipitation and gauge fusion (FY4B, GPM, GSMaP)",
  date: datetime(year: 2026, month: 10, day: 23),
)
#set page(
  paper: "a4",
  margin: (x: 2.1cm, top: 2.3cm, bottom: 2.2cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 8pt, fill: luma(110))
      Research plan — FY4B, GPM and GSMaP
      #h(1fr) 23 October 2026
    ]
  },
)
#set text(font: ("Libertinus Serif", "New Computer Modern"), size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.62em, spacing: 1.05em)
#set heading(numbering: "1.1")
#show heading: set text(font: "Arial", weight: "semibold")
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

// Every result the text cites is read from paper/numbers/*.json, never typed by hand:
// scripts/export_paper_numbers.py exports them from the source CSVs named in paper/sources.toml, and
// scripts/check_paper_provenance.py fails if a JSON no longer matches its source. A missing row
// prints a red "?".
#let load(name) = json("numbers/" + name + ".json").rows
#let mu = load("methods_uncal")   // pooled metrics on the matched sample (Table 3)
#let cu = load("vs_adw_uncal")    // bootstrap comparisons with ADW, by stratum
#let su = load("strata_uncal")    // gain over ADW and mean bias by rain intensity and distance to a gauge
#let q3 = load("rq3_sample")      // §2.4: the gauge-hours every anchor and method is scored on
#let r1 = load("rq1_summary")      // RQ1: mean over the six events
#let e1 = load("rq1_events")       // RQ1: scores per event
#let k1 = load("rq1_cores")        // RQ1: storm-core deficit per event and product
#let c2 = load("rq2_classes")      // RQ2: by intensity class, season and relief
#let v2 = load("rq2_events")       // RQ2: rain-event detection and timing
#let a2 = load("rq2_areal")        // RQ2: areal-mean series
#let n0 = load("study_gauges")     // Study area: gauges RQ2 scores, in all and by landform region
#let dh = load("data_hours")       // §2.1/§2.5: the matched sample's hours, days and gauge-hours
#let fe = load("fy4b_events_summary") // §3.4, Table 4: corrected FY4B on the RQ1 events, mean over the events
#let fv = load("fy4b_event_scores")   // §3.4: scores per RQ1 event
#let fq = load("fy4b_quality")     // §3.4: corrected FY4B against plain FY4B, GPM and GSMaP, FY4B run
#let ft = load("fy4b_time_holdout") // §3.4: corrected FY4B's raw change with whole weeks held out as well
#let fc = load("fy4b_clt")         // §3.4: July 2023 cloud-type pilot, hours dominated by ice cloud
#let fcf = load("fy4b_clt_filter") // §3.4: the cloud-type filter against an amount-matched control
#let G = "GPM_uncal"              // the primary GPM product
#let B = "blend_agrenv_mgwr"      // MGWR + agreement blend, the headline fusion method
#let find(rows, ..conditions) = {
  let hits = rows.filter(row => conditions.named().pairs().all(pair => row.at(pair.at(0)) == pair.at(1)))
  if hits.len() == 0 { none } else { hits.at(0) }
}
/// One value as typeset text: ASCII hyphens become minus signs.
#let val(rows, column, ..conditions) = {
  let hit = find(rows, ..conditions)
  if hit == none { text(fill: red)[?] } else { hit.at(column).replace("-", "−") }
}
/// A percentage for prose: a leading plus is dropped (the sentence says "gains"), a minus is kept.
#let pct(rows, column, ..conditions) = {
  let hit = find(rows, ..conditions)
  if hit == none { text(fill: red)[?] } else [#hit.at(column).trim("+", at: start).replace("-", "−")%]
}
/// A percentage with its sign dropped, for "x% worse".
#let mag(rows, column, ..conditions) = {
  let hit = find(rows, ..conditions)
  if hit == none { text(fill: red)[?] } else [#hit.at(column).trim("+", at: start).trim("-", at: start)%]
}
/// A number worked out from several rows (a range, a median, a largest gap), exported under `stats`.
#let stat(name, key) = {
  let stats = json("numbers/" + name + ".json").at("stats", default: (:))
  if key in stats { stats.at(key).replace("-", "−") } else { text(fill: red)[?] }
}
/// One Table 3 row's numbers: RMSE, r, CSI 2, vs raw (a dash on gauge-only rows) and vs ADW.
#let t3(rows, product, method, raw: true) = (
  val(rows, "RMSE", product: product, method: method),
  val(rows, "r", product: product, method: method),
  val(rows, "CSI2", product: product, method: method),
  if raw [#val(rows, "vs_raw", product: product, method: method)%] else [—],
  [#val(rows, "vs_adw", product: product, method: method)%],
)

/// Corrected FY4B's change in RMSE against `baseline`, same method, with its 95% interval over
/// whole 7-day blocks, signed: "−16.9% (−27.6 to −7.0%)".
#let dq(baseline, method, scheme: "balanced_spatial") = {
  let row = (scheme: scheme, product: "FY4B_QM", baseline: baseline, method: method)
  [#val(fq, "change", ..row)% (#val(fq, "lo", ..row) to #val(fq, "hi", ..row)%)]
}

/// Raw corrected FY4B against plain FY4B in the check that holds out whole weeks as well as gauges
/// (balanced spatial), with 95% intervals over whole weeks: its change in RMSE, "−9.5% (−15.4 to
/// −3.5%)" (`dt`), or in r, "−0.004 (−0.021 to +0.015)" (`dtr`).
#let dt(candidate) = {
  let row = (scheme: "balanced_spatial", candidate: candidate, baseline: "FY4B")
  [#val(ft, "change", ..row)% (#val(ft, "lo", ..row) to #val(ft, "hi", ..row)%)]
}
#let dtr(candidate) = {
  let row = (scheme: "balanced_spatial", candidate: candidate, baseline: "FY4B")
  [#val(ft, "dr", ..row) (#val(ft, "dr_lo", ..row) to #val(ft, "dr_hi", ..row))]
}

/// One Table 4 row: a product's means over the RQ1 events, over their shared hours (r, MB, RMSE, CSI 50).
#let t4(product) = ("r", "MB", "RMSE", "CSI50").map(column =>
  val(fe, column, product: product, window: "fy4b_hours"))
/// One product's spatial r on one RQ1 event, over the event's shared hours (RQ1).
#let ev(product, day) = val(e1, "r", product: product, event_day: day)
/// One product's spatial r on one RQ1 event in §3.4's run, where corrected FY4B is scored too.
#let evq(product, day) = val(fv, "r", product: product, event_day: day)

/// Mean bias (mm h⁻¹) of one anchor and method on hours of ≥ 20 mm h⁻¹ at the held-out gauges.
#let bias20(product, method) = val(su, "bias", product: product, method: method,
  group: "rain_intensity", level: "severe_rainstorm")

#align(center)[
  #text(size: 17pt, weight: "semibold", font: "Arial")[Research plan]
  #v(0.2em)
  #text(size: 11.5pt, fill: luma(60))[Correcting FY4B, GPM and GSMaP satellite precipitation against rain gauges with GWR / MGWR]
]

#note[
  *Problem.* Flood forecasting in mountain catchments needs hourly heavy rain, which a sparse gauge
  network cannot resolve. Satellite products cover every pixel, but how their error depends on place,
  intensity and hour is poorly known, and so is what fusing them with the gauges adds far from the
  gauges and on heavy hours.

  *Objective.* To correct FY4B, GPM and GSMaP satellite precipitation against rain gauges with GWR /
  MGWR. RQ1 and RQ2 measure the errors to correct, in space and in time; RQ3 tests the correction.
  Three research contents, one results section each (§3), all on the gauge-hours where the gauge and
  all three products have a value.
  #v(0.3em)
  #table(
    columns: (auto, 1.4fr, 1.8fr, 0.8fr, 1.1fr),
    align: left,
    stroke: none,
    rule,
    table.header([*RQ*], [*Question*], [*Data and method*], [*Status*], [*To finish*]),
    thin,
    [1], [Space: does the satellite put heavy rain in the right place, at the right amount?],
    [Six heavy-rain days; totals over the hours all three products share, at the gauges, scored by
    spatial r, bias, RMSE and CSI at 50 mm],
    [Answered], [—],
    [2], [Time: does it follow the hourly evolution of rain at each intensity?],
    [Gauge-wet hours in five intensity classes, by season and relief; rain events; areal means],
    [Answered], [—],
    [3], [Fusion: what does fusing it with the gauges gain, with which method and anchor?],
    [Gauge-only and fusion methods on FY4B, GPM, GSMaP and their merge, scored at held-out gauges by
    spatial cross-validation],
    [Answered], [—],
    rule,
  )
]

= Introduction

Flood forecasting in mountain catchments needs hourly rainfall, and above all the moderate and heavy
rain that drives runoff. Rain gauges measure it directly but are sparse in the mountains, so large
parts of a catchment sit tens of kilometres from the nearest gauge. Satellite products such as GPM
IMERG and GSMaP cover every pixel every hour, and geostationary imagers such as that of China's
FengYun-4B every 15 minutes, but their errors depend on rain type and intensity.

Most evaluations of these products work at daily or monthly scale and report pooled scores. Pooled
scores do not say whether a product puts a heavy-rain event in the right place, whether it follows the
hour-by-hour evolution of rain of a given intensity, or where its error comes from. Gauge–satellite
fusion is widely used to correct the products, but it is usually judged by average skill. What fusion
adds *far from the gauges* and *on heavy hours*, which is exactly the sparse-network flood case, is less
clear.

*Objective.* To correct FY4B, GPM and GSMaP satellite precipitation against rain gauges with GWR /
MGWR. RQ1 and RQ2 measure the errors to correct, in space and in time; RQ3 tests the correction. §3
takes each as a research content, with the questions it must answer:

- *RQ1 — Space.* Does the satellite field reproduce *where* heavy rain falls?
- *RQ2 — Time.* Does it reproduce *how rain evolves in time* at each intensity?
- *RQ3 — Fusion.* What does gauge–satellite fusion gain, which fusion method and which anchor product
  are best, and does fusion fix the moderate and heavy rain that floods depend on?

*Scope.* RQ1–RQ3 compare three products, none of them adjusted to gauges: GPM IMERG, GSMaP and
FY4B, the rain estimate of China's geostationary FengYun-4B, rebuilt here from its 15-minute files
(§2.5). FY4B's record starts only in June 2022 and has gaps, so every comparison is made on one
matched sample: the gauge-hours from June 2022 to December 2024 at which the gauge and all three
products have a value (§2.1). §3.4 asks in addition how far FY4B's quality can be raised by a
correction fitted to the gauges.

= Data and Methods

== Study area and data

#val(n0, "n", region: "all") hourly gauges in the Shiyan study area, June 2022 – December 2024 (values quantised to 0.5 mm). All sit in hills
and mountains: by local relief from the Copernicus GLO-30 DEM, #val(n0, "n", region: "relief_lt500") are below 500 m, #val(n0, "n", region: "relief_500_1000") at 500–1000 m
and #val(n0, "n", region: "relief_ge1000") at 1000 m or more. GPM IMERG V07 Uncal @huffman_integrated_2020 @huffman_nasa_2023 and GSMaP v8 @kubota_global_2020, both at
0.1°, are sampled at the gauge pixels: IMERG Uncal from NASA GES DISC's half-hourly files, averaged to
hours, and GSMaP as the satellite-only `hourlyPrecipRate` of Google Earth Engine's operational
collection, which for 2022–2024 holds the permanent standard product GSMaP_MVK. Neither is adjusted
to gauges. IMERG Uncal is the multi-satellite estimate before IMERG's monthly calibration to the GPCC
gauge analysis. Its microwave estimates are calibrated to GPM's radar–radiometer product, and V07
applies a climatological ratio to the satellite–gauge GPCP analysis only over high-latitude oceans
@huffman_nasa_2023; in V05 that ratio also covered all land @huffman_integrated_2020. GSMaP_MVK has no
gauge correction @kubota_global_2020. FY4B, the third product, is rebuilt from FengYun-4B's
satellite-only 15-minute files (§2.5). Days run 08–08 Beijing time. FY4B's files start on 1 June 2022
and have gaps (§2.5), so every evaluation uses one matched sample: the gauge-hours at which the gauge
and all three products have a value, #val(dh, "cells", sample: "all_products") of them on
#val(dh, "hours", sample: "all_products") hours and #val(dh, "days", sample: "all_products") days of
June 2022 – December 2024.

#figure(
  image("figures/fig01_study_area.png", width: 80%),
  caption: [Draft Fig. 1. The study area: shaded relief from a 1 km SRTM DEM, and the hourly gauges by
  landform region, classed by local relief in the Copernicus GLO-30 DEM. Hollow markers: the region has
  too few gauges to be interpreted on its own.],
)

== Heavy-rain events and spatial scores (RQ1)

Heavy-rain days are 08–08 days from June 2022 on which more than three gauges exceed 50 mm. Each day
is scored over its shared hours: those at which the gauges and each of the three products have a
value at 90% of the gauges or more. FY4B never has the hour ending 23:00 BJT, so no day is complete. A
day is eligible when its shared hours number at least 20 and leave out at most 20% of its gauge rain,
over the network and at its wettest gauge. Of #stat("rq1_days", "days") heavy-rain days,
#stat("rq1_days", "eligible") are eligible, and six are selected, one per rain system: four
widespread, two localized. At each gauge with a value in every shared hour, the satellite total over
those hours is compared with the gauge total by spatial correlation, spatial CV ratio, rain-centre
shift, mean bias, RMSE and CSI at 50 mm.

== Hourly evaluation by intensity (RQ2)

Only gauge-wet hours are scored, in five intensity classes (0.1–2, 2–4, 4–8, 8–20 and ≥ 20 mm h⁻¹),
stratified by season and landform region, with 95% bootstrap intervals that resample whole 7-day blocks
(rain episodes often span several days, so single days are not independent). Rain events are runs of
wet hours separated by at least 3 dry hours and are scored for detection, timing and volume. Areal-mean
series are compared at hourly and 3-hourly scale; with the hour ending 23:00 BJT always missing,
there is no daily series. All of it is scored on the matched sample (§2.1).

== Gauge–satellite fusion and its validation (RQ3)

Gauge-only interpolators (IDW, ADW, TPS, GWR) are the baseline. Fusion methods use a satellite product
as the anchor and correct it with the gauges (residual GWR, mixed GWR, MGWR), optionally blended with
ADW where the satellite reports rain. The agreement blend chooses the ADW weight separately for each
band of how many of the three products call the cell wet (≥ 0.1 mm h⁻¹) and the largest of their
rates (edges 0.1, 2.5 and 8 mm h⁻¹). MGWR with the agreement blend is the headline fusion method.
Anchors are FY4B, GPM (IMERG Uncal), GSMaP and an OLS merge of the three; corrected FY4B (§2.5) is a
further anchor that enters neither the merge nor the bands.
Every method is scored at gauges held out of the fit, by balanced spatial 5-fold cross-validation over
FY4B's hours of June 2022 – December 2024, using RMSE, correlation and CSI, overall and by intensity
and distance to the nearest training gauge. Every anchor and method is scored on the same
#val(q3, "n", scheme: "balanced_spatial") gauge-hours: those of the matched sample (§2.1) at which
every anchor's prediction by every method exists, on #val(q3, "hours", scheme: "balanced_spatial")
hours. Differences from ADW are tested by a paired bootstrap that resamples whole 7-day blocks
(counted from 1 June 2022) 2,000 times, with Holm's correction across the methods and anchors
compared in each stratum. These intervals, like those of RQ2 and §3.4, resample time only: the
gauges and the five cross-validation folds stay fixed, so they do not cover which gauges the network
has or how they were split.

== Raising FY4B's quality

FY4B's rain estimate comes from the AGRI imager of FengYun-4B. It is rebuilt here from the
satellite-only 15-minute full-disk files (4 km); the hourly and longer files are corrected with
China's rain gauges and are not used. An hour is kept only when all four of its 15-minute files exist,
and its value at a gauge is the mean over the 3 × 3 pixels around it. Its files start on 1 June
2022, and none exist for February 2024, between the satellite's last hours at 133°E and its first at
105°E. The 14:15 UTC file is never produced, so the hour 22–23 BJT is always missing, and the strict
rule drops other scattered hours. FY4B thus has a value at the gauges on
#val(dh, "hours", sample: "all_products") hours on #val(dh, "days", sample: "all_products") days of
June 2022 – December 2024, the hours of the matched sample (§2.1). This rebuilt series is the FY4B of
RQ1–RQ3 and, before the corrections below, the plain FY4B of §3.4. Two corrections
are then fitted inside each cross-validation fold, on the training gauges only, so a held-out gauge never shapes the
anchor it is scored against. First, a pixel shift. FY4B falls into three regimes: two versions of its
processing at 133°E, switched in May 2023, and the satellite at 105°E from March 2024. For each
regime the training gauges choose the shift under which FY4B correlates best with them. Second, a
quantile map per season and regime, with a wet threshold that matches the gauges' wet-hour frequency.
It is applied month by month, only as strongly as the training gauges support: not at all, the
threshold alone, or the threshold with a quarter to all of the map. The corrected product, FY4B_QM,
is scored in the same benchmark run as RQ3, as a further anchor beside the three products and on
the same matched sample (§2.4); it enters neither the merge nor the agreement bands. Corrected FY4B
is also scored on the RQ1 events (§2.2), over the hours of each event that the gauges and all four
products share, which are the hours RQ1 scores; its value at each gauge comes from the cross-validation fold
that held that gauge out.

= Results: research content

Each research content lists the questions it must answer, what has been found so far, and what is
left before it is finished.

== RQ1 — Does the satellite place heavy rain correctly in space?

*Questions to answer.*
- Is the spatial pattern of event rainfall reproduced?
- Are the event amounts right?
- Are localized storm cores resolved, and if not, is the cause pixel averaging or rain spread over too
  many hours?

*Findings so far.*
- *Placement is similar for GPM and GSMaP, and fair; FY4B's is worse on average.* Mean spatial r
  over the six events is #val(r1, "r", product: G) for GPM, #val(r1, "r", product: "GSMaP") for GSMaP and
  #val(r1, "r", product: "FY4B") for FY4B; GPM and GSMaP score #stat("rq1_events", "r_four") on four events. All three fail on 2022-06-26
  (r ≤ #stat("rq1_events", "r_0626_max")), when the rain moved north to south and left a nearly flat gauge total that the
  products' north–south gradient does not match. Event by event FY4B is mixed: on 2022-10-03 its pattern is
  the best of the three (r #ev("FY4B", "2022-10-03") against #ev(G, "2022-10-03") and
  #ev("GSMaP", "2022-10-03")), but it fails on 2023-08-26 (#ev("FY4B", "2023-08-26")) and trails on
  both localized events (#ev("FY4B", "2023-07-31") and #ev("FY4B", "2023-09-10")).
- *Amounts differ sharply.* GPM is somewhat dry (MB #val(r1, "MB", product: G) mm, RMSE #val(r1, "RMSE", product: G) mm), GSMaP wet
  (MB #val(r1, "MB", product: "GSMaP") mm, RMSE #val(r1, "RMSE", product: "GSMaP") mm). FY4B's mean bias is near zero
  (#val(r1, "MB", product: "FY4B") mm) only because one wet day
  (#val(e1, "MB0", product: "FY4B", event_day: "2022-06-26") mm on 2022-06-26) offsets dry ones; its
  RMSE is #val(r1, "RMSE", product: "FY4B") mm.
- *Localized storm cores are missed by all three* (CSI at 50 mm on 2023-07-31: GPM
  #val(e1, "CSI50", product: G, event_day: "2023-07-31"), GSMaP #val(e1, "CSI50", product: "GSMaP", event_day: "2023-07-31"), FY4B #val(e1, "CSI50", product: "FY4B", event_day: "2023-07-31")\;
  on 2023-09-10 #val(e1, "CSI50", product: G, event_day: "2023-09-10") for every product). How much of
  the core deficit pixel averaging explains depends on the storm. On 2023-09-10 it explains only
  #stat("rq1_cores", "repr_0910")%, and the rest is the products spreading short bursts over more hours.
  On 2023-07-31, whose core is narrower than a pixel (the 10 km mean around the core gauges is well
  below their totals), it explains #stat("rq1_cores", "repr_0731")%. On both days the products reach
  #stat("rq1_cores", "peak_ratio_local") of the gauges' peak hourly rate.

Corrected FY4B on these events is in §3.4.

#figure(
  image("figures/fig02_event_maps.png", width: 100%),
  caption: [Draft Fig. 2. The six heavy-rain events (columns). Top row: gauge totals over the hours the
  gauges and all three products share within the 08–08 BJT day (§2.2). Lower rows: satellite minus
  gauge at each gauge, from the pixel sampled there; red means the satellite is
  drier, blue wetter, grey within ±10 mm. MB is the mean bias (mean of satellite − gauge over the
  event's gauges). GPM is IMERG Uncal.],
)

*To finish.* Nothing: the question is answered.

== RQ2 — Does the satellite follow the hourly evolution at each intensity?

*Questions to answer.*
- How do bias and the share of hours in the right class change with intensity, season and relief?
- Are rain events detected, and timed correctly?
- Are areal means better than single points?

*Findings so far.*
- *Sorted by the gauge's class, the shortfall grows with intensity, in every season.* GPM relative bias runs from #val(c2, "RB0", product: G, season: "all", region: "all", gauge_class: "light")% on light
  rain to #val(c2, "RB0", product: G, season: "all", region: "all", gauge_class: "severe_rainstorm")% at ≥ 20 mm h⁻¹, GSMaP's from #val(c2, "RB0", product: "GSMaP", season: "all", region: "all", gauge_class: "light")% to #val(c2, "RB0", product: "GSMaP", season: "all", region: "all", gauge_class: "severe_rainstorm")% and FY4B's from
  #val(c2, "RB0", product: "FY4B", season: "all", region: "all", gauge_class: "light")% to #val(c2, "RB0", product: "FY4B", season: "all", region: "all", gauge_class: "severe_rainstorm")%. GPM is drier than GSMaP
  in every class, and FY4B drier than GPM up to 8 mm h⁻¹. At ≥ 8 mm h⁻¹ at most #stat("rq2_classes", "hit_heavy_max")% of hours land in the right class.
  Sorting hours by the gauge's own reading makes any imperfect estimate look low on its heaviest hours,
  so part of this shortfall comes from the sorting. Counted without it, raw GPM reaches ≥ 20 mm h⁻¹
  #val(mu, "FBI20", product: G, method: "raw") times as often as the gauges, GSMaP #val(mu, "FBI20", product: "GSMaP", method: "raw") times and FY4B
  #val(mu, "FBI20", product: "FY4B", method: "raw") times. A 0.1° pixel averages a short burst over
  about 10 km, so even a perfect product would reach that rate less often than a point gauge. But
  GSMaP, on the same grid, reaches it more often than the gauges, so pixel size alone does not explain
  GPM's shortfall: GPM is short of heavy hours, and GSMaP has them in the wrong place or hour.
- *For GPM and GSMaP, detection and timing are not the problem; magnitude is. FY4B also misses
  rain.* Rain events peaking at 2–4 mm h⁻¹ or more are detected #stat("rq2_events", "pod")% of the time
  by GPM and GSMaP but #stat("rq2_events", "pod_fy4b")% by FY4B, and median rain-centre timing errors
  stay between #stat("rq2_events", "centroid_min") and #stat("rq2_events", "centroid_max") h for all three.
  FY4B calls only #pct(c2, "POD", product: "FY4B", season: "all", region: "all", gauge_class: "all_wet") of gauge-wet hours wet, against
  #pct(c2, "POD", product: G, season: "all", region: "all", gauge_class: "all_wet") for GPM and #pct(c2, "POD", product: "GSMaP", season: "all", region: "all", gauge_class: "all_wet") for GSMaP.
- *Season matters, terrain does not.* All three fail in winter (bias #val(c2, "RB0", product: G, season: "DJF", region: "all", gauge_class: "all_wet")% for GPM,
  #val(c2, "RB0", product: "GSMaP", season: "DJF", region: "all", gauge_class: "all_wet")% for GSMaP and #val(c2, "RB0", product: "FY4B", season: "DJF", region: "all", gauge_class: "all_wet")% for FY4B), while between relief below 500 m and 500–1000 m the class biases differ
  by at most #stat("rq2_classes", "relief_gap_max") percentage points.
- *Areal means are much better than points* for GPM and GSMaP (hourly r #val(a2, "r", product: G, timescale_h: "1") and #val(a2, "r", product: "GSMaP", timescale_h: "1"),
  against #val(a2, "r", product: "FY4B", timescale_h: "1") for FY4B).

#figure(
  image("figures/fig04a_intensity_season.png", width: 100%),
  caption: [Draft Fig. 4. Hourly skill of FY4B, GPM (IMERG Uncal) and GSMaP by gauge intensity class
  and season, on the matched sample.],
)

*To finish.* Nothing: the question is answered.

== RQ3 — What does gauge fusion gain, with which method and which product?

*Questions to answer.*
- Does fusion beat both the raw product and gauges alone, overall and on heavy hours?
- How does the gain change with distance to the nearest gauge?
- Which fusion method and which anchor product are best?

*Findings so far.*
- *Fusion beats both of its inputs on GPM and on the merge, not on GSMaP or FY4B.* MGWR followed by a
  blend with ADW where the satellites report rain cuts RMSE by #pct(mu, "vs_raw", product: G, method: B) from raw GPM, #pct(mu, "vs_raw", product: "GSMaP", method: B) from raw GSMaP
  and #pct(mu, "vs_raw", product: "FY4B", method: B) from raw FY4B. Against gauge-only ADW it
  gains #pct(mu, "vs_adw", product: G, method: B) on GPM and #pct(mu, "vs_adw", product: "MERGED_OLS", method: B) on the merged anchor, both significant, but only
  #pct(mu, "vs_adw", product: "GSMaP", method: B) on GSMaP (p = #val(cu, "p_holm", product: "GSMaP", method: B, stratum: "overall")) and #pct(mu, "vs_adw", product: "FY4B", method: B) on FY4B
  (p = #val(cu, "p_holm", product: "FY4B", method: B, stratum: "overall")), which are not.
- *The gain holds on heavy hours and grows far from the gauges, but the heavy-rain peaks are not
  restored.* On heavy hours
  (≥ 8 mm h⁻¹) MGWR alone improves on ADW by #pct(cu, "change", product: G, method: "mgwr", stratum: "8_plus") (GPM), #pct(cu, "change", product: "GSMaP", method: "mgwr", stratum: "8_plus") (GSMaP) and
  #pct(cu, "change", product: "FY4B", method: "mgwr", stratum: "8_plus") (FY4B), the agreement
  blend by #pct(cu, "change", product: G, method: B, stratum: "8_plus"), #pct(cu, "change", product: "GSMaP", method: B, stratum: "8_plus") and #pct(cu, "change", product: "FY4B", method: B, stratum: "8_plus"), all significant. Unblended MGWR *hurts* dry, light and moderate hours,
  which is why the blend is needed. On GPM the improvement over ADW is #pct(su, "change", product: G, method: B, group: "nearest_train_km", level: "0_20") within 20 km of a training gauge,
  #pct(su, "change", product: G, method: B, group: "nearest_train_km", level: "20_50") at 20–50 km and #pct(su, "change", product: G, method: B, group: "nearest_train_km", level: "50_100") at 50–100 km. #cite(<hu_satellite_2015>, form: "prose") found a matching pattern by thinning a gauge network:
  merging gains little over interpolation where gauges are dense. At ≥ 20 mm h⁻¹ the mean bias of fused GPM is #bias20(G, B) mm h⁻¹, against
  #bias20(G, "raw") for raw GPM and #bias20(G, "adw") for gauge-only ADW. On GSMaP the fused field (#bias20("GSMaP", B)) is drier there than the
  raw product (#bias20("GSMaP", "raw")). Sorting by the gauge's reading pulls every estimate low on these
  hours. Counted without it, fused GPM reaches ≥ 20 mm h⁻¹ #val(mu, "FBI20", product: G, method: B) times as often as the gauges, against
  #val(mu, "FBI20", product: G, method: "raw") for raw GPM and #val(mu, "FBI20", product: G, method: "adw") for gauge-only ADW. Fusion flattens the peaks further, as
  expected: an average of gauges, or a fit for the lowest mean error, pulls extremes toward typical values.
- *Best method and product.* Among gauge-only interpolators ADW is best, with IDW indistinguishable
  (#val(mu, "RMSE", product: G, method: "adw") vs #val(mu, "RMSE", product: G, method: "idw") mm h⁻¹). The best fusion is MGWR + blend. GPM is the
  best of the three products as an anchor, with none gauge-calibrated (raw RMSE #val(mu, "RMSE", product: G, method: "raw") vs #val(mu, "RMSE", product: "GSMaP", method: "raw") for GSMaP and
  #val(mu, "RMSE", product: "FY4B", method: "raw") for FY4B, mm h⁻¹), and of
  the three only fused GPM beats ADW significantly. Merging the three gives a better raw field (RMSE #val(mu, "RMSE", product: "MERGED_OLS", method: "raw")) but adds
  little after fusion (#val(mu, "RMSE", product: "MERGED_OLS", method: B) vs #val(mu, "RMSE", product: G, method: B)).

#figure(
  image("figures/fig06b_fusion_methods.png", width: 100%),
  caption: [Draft Fig. 6. RMSE improvement of each fusion method over gauge-only ADW, per anchor, for
  all hours and heavy hours, on the matched sample. Rows are grouped by how the satellite enters: direct
  correction, blend with ADW, agreement blend. Hollow cells are not significant. The equal and OLS merges
  are of FY4B, GPM and GSMaP; corrected FY4B (§3.4) is the last column.],
)

#figure(
  table(
    columns: (auto, auto, auto, auto, auto, auto, auto),
    align: (left, left, right, right, right, right, right),
    stroke: none,
    rule,
    table.header([*Anchor*], [*Method*], [*RMSE*], [*r*], [*CSI 2*], [*vs raw*], [*vs ADW*]),
    thin,
    [—], [ADW (gauges only)], ..t3(mu, "GPM_uncal", "adw", raw: false),
    [—], [IDW (gauges only)], ..t3(mu, "GPM_uncal", "idw", raw: false),
    [—], [TPS (gauges only)], ..t3(mu, "GPM_uncal", "tps", raw: false),
    [—], [GWR (gauges only)], ..t3(mu, "GPM_uncal", "gwr", raw: false),
    thin,
    [FY4B], [raw], ..t3(mu, "FY4B", "raw"),
    [FY4B], [MGWR + agreement blend], ..t3(mu, "FY4B", "blend_agrenv_mgwr"),
    [GPM], [raw], ..t3(mu, "GPM_uncal", "raw"),
    [GPM], [MGWR + agreement blend], ..t3(mu, "GPM_uncal", "blend_agrenv_mgwr"),
    [GSMaP], [raw], ..t3(mu, "GSMaP", "raw"),
    [GSMaP], [MGWR + agreement blend], ..t3(mu, "GSMaP", "blend_agrenv_mgwr"),
    [Merged (OLS)], [raw], ..t3(mu, "MERGED_OLS", "raw"),
    [Merged (OLS)], [MGWR + agreement blend], ..t3(mu, "MERGED_OLS", "blend_agrenv_mgwr"),
    thin,
    [Corrected FY4B], [raw], ..t3(mu, "FY4B_QM", "raw"),
    [Corrected FY4B], [MGWR + agreement blend], ..t3(mu, "FY4B_QM", "blend_agrenv_mgwr"),
    rule,
  ),
  caption: [Draft Table 3. Held-out gauges, balanced spatial CV on the matched sample of June 2022 –
  December 2024: every row is scored on the same #val(q3, "n", scheme: "balanced_spatial")
  gauge-hours (§2.4). RMSE in mm h⁻¹; improvements are RMSE reductions; CSI at 2 mm h⁻¹. GPM is IMERG
  Uncal. The merge and the agreement blend are built from FY4B, GPM and GSMaP; corrected FY4B (§2.5)
  enters neither. The run's commit is recorded in its benchmark_scope.csv, with no uncommitted changes.],
)

*To finish.* Nothing: the question is answered.

== FY4B — how far can its quality be raised?

Corrected FY4B (§2.5) gains at every stage of the hourly evaluation. Every comparison here is paired,
on the matched sample of Table 3 (§2.4). Against plain FY4B (the rebuilt series before the shift
and the map) on the same cells, its RMSE changes by
#dq("FY4B", "raw") raw, #dq("FY4B", "mgwr") after MGWR and #dq("FY4B", B) after MGWR and the
agreement blend (95% intervals over whole 7-day blocks, balanced spatial CV), and raw r over all gauge-hours rises
from
#val(fq, "r0", scheme: "balanced_spatial", baseline: "FY4B", method: "raw") to
#val(fq, "r", scheme: "balanced_spatial", baseline: "FY4B", method: "raw"). These gains come from a
correction fitted on gauges in the same weeks as the held-out gauges it is scored at: a held-out gauge
never shapes the map, but it shares those weeks' weather. With whole weeks held out as well (7-day
blocks, a day either side, crossed with the gauge folds), the raw change is #dt("FY4B_QM_crossed"),
against #dt("FY4B_QM_gauges") for the same correction, both with 95% intervals over whole weeks, and
raw r no longer rises: it changes by #dtr("FY4B_QM_crossed"). Most of the difference falls on a few false-rain
episodes whose largest values the map learns to shrink. MGWR and the blend were not refitted this way,
so their gains above hold for the same weeks only. Corrected FY4B does not reach GPM:
against GPM's agreement blend its RMSE differs by #dq(G, B), and after MGWR alone by
#dq(G, "mgwr"), an interval that spans zero but lies mostly on the worse side. Against GSMaP it does
better after MGWR alone, #dq("GSMaP", "mgwr"), though its r stays lower
(#val(fq, "r", scheme: "balanced_spatial", baseline: "GSMaP", method: "mgwr") vs
#val(fq, "r0", scheme: "balanced_spatial", baseline: "GSMaP", method: "mgwr")); after the blend the
difference is #dq("GSMaP", B), an interval that spans zero. The gain over plain FY4B after MGWR shows only where held-out gauges are far from
training ones: in the random split, where a held-out gauge usually has a training gauge nearby, it is
#dq("FY4B", "mgwr", scheme: "random"). What
remains is false rain, and FY4B's own cloud-type product cannot pick it out. In a one-month pilot (July
2023), hours dominated by ice cloud held #pct(fc, "false_alarms", name: "ice") of corrected FY4B's
false alarms but also #pct(fc, "hits", name: "ice") of its hits; their false-alarm ratio was
#val(fc, "far", name: "ice"), against #val(fc, "far_gpm", name: "ice") for GPM. Zeroing FY4B where the
training gauges' chosen cloud class dominated raised CSI at 0.1 mm h⁻¹ by only
#val(fcf, "csi", candidate: "FY4B_QM_CLT") (#val(fcf, "csi_lo", candidate: "FY4B_QM_CLT") to
#val(fcf, "csi_hi", candidate: "FY4B_QM_CLT"), over day blocks, as one month holds too few weeks) over
removing as many of its lightest rain hours.

On the RQ1 heavy-rain events the hourly gain does not carry over to storm totals. Draft Table 4 scores
plain and corrected FY4B with GPM and GSMaP on the six events of §3.1, each over its shared hours,
which corrected FY4B does not change. The correction makes most events' totals drier. That removes
part of the wet bias of 2022-06-26 (#val(fv, "MB", product: "FY4B", event_day: "2022-06-26") to
#val(fv, "MB", product: "FY4B_QM", event_day: "2022-06-26") mm), which is why the mean RMSE falls
slightly, from #val(fe, "RMSE", product: "FY4B", window: "fy4b_hours") to
#val(fe, "RMSE", product: "FY4B_QM", window: "fy4b_hours") mm, but it deepens the deficit on dry days:
on 2022-10-03 the bias goes from #val(fv, "MB", product: "FY4B", event_day: "2022-10-03") to
#val(fv, "MB", product: "FY4B_QM", event_day: "2022-10-03") mm (CSI at 50 mm
#val(fv, "CSI50", product: "FY4B", event_day: "2022-10-03") to
#val(fv, "CSI50", product: "FY4B_QM", event_day: "2022-10-03")). Mean CSI at 50 mm falls from
#val(fe, "CSI50", product: "FY4B", window: "fy4b_hours") to
#val(fe, "CSI50", product: "FY4B_QM", window: "fy4b_hours"). Mean r barely moves (from
#val(fe, "r", product: "FY4B", window: "fy4b_hours") to
#val(fe, "r", product: "FY4B_QM", window: "fy4b_hours")): it falls on 2023-09-10
(#evq("FY4B", "2023-09-10") to #evq("FY4B_QM", "2023-09-10")) and rises on 2023-07-31
(#evq("FY4B", "2023-07-31") to #evq("FY4B_QM", "2023-07-31")).

#figure(
  table(
    columns: (auto, auto, auto, auto, auto),
    align: (left, right, right, right, right),
    stroke: none,
    rule,
    table.header([*Product*], [*r*], [*MB*], [*RMSE*], [*CSI 50*]),
    thin,
    [FY4B], ..t4("FY4B"),
    [Corrected FY4B], ..t4("FY4B_QM"),
    [GPM], ..t4(G),
    [GSMaP], ..t4("GSMaP"),
    rule,
  ),
  caption: [Draft Table 4. Plain and corrected FY4B on the RQ1 heavy-rain events: totals at the
  gauges, mean over the #val(fe, "n", product: "FY4B", window: "fy4b_hours") events. Each event is
  scored over the #stat("fy4b_event_windows", "hours") hours of its 08–08 day that the gauges and all
  four products share, the same hours and gauges for every row and the same as in §3.1.
  Corrected FY4B's value at each gauge comes from the cross-validation fold that held that gauge out.
  MB is the mean bias and RMSE the root-mean-square error, both in mm; CSI at 50 mm. GPM is IMERG Uncal.],
)

Corrected FY4B therefore stays out of the merge and the agreement blend.

= Discussion

== Uncertainties

- *Point versus pixel.* A gauge samples a point, the products a 0.1° pixel (FY4B a 3 × 3 box of 4 km
  pixels). Where rain varies below 10 km (2022-06-26, the localized cores) part of the "error" is
  representativeness, not product skill; pixel averaging explains #stat("rq1_cores", "repr_0910")% of
  the core deficit on 2023-09-10, where the rest comes from the products spreading short bursts over
  too many hours, but #stat("rq1_cores", "repr_0731")% on 2023-07-31 (§3.1). The same averaging means a
  pixel reaches heavy rates less often than a gauge (§3.2).
- *Gauge resolution.* The 0.5 mm quantisation limits what can be said about the lightest class.
- *Small samples.* Only #val(n0, "n", region: "relief_ge1000") gauges lie at relief ≥ 1000 m, too few to interpret; heavy hours and the six
  events are few, so their intervals are wide.
- *The matched sample.* Requiring FY4B leaves out January–May 2022, February 2024, the hour ending
  23:00 BJT and scattered hours (§2.5). Winter and spring are therefore thin, and no heavy-rain event is
  scored over its whole day (§2.2). What the three products are compared on is the same, but it is
  not the full record of GPM or GSMaP.
- *Residual heavy-rain deficit.* Fusion lowers the RMSE of heavy hours but does not restore their
  peaks (§3.3). Gauge-only ADW is as low as raw GPM there, so neither input carries the peak to a
  held-out gauge. Part of this is built into the scores and the methods: sorting hours by the gauge's
  reading pulls the bias low, and averaging or least-squares fitting smooths peaks. Event or basin
  totals would measure the lost volume with neither effect. Matching the satellite's intensity
  distribution to the gauges before merging targets this loss: PDF matching, as
  #cite(<xie_conceptual_2011>, form: "prose") do over China, or the quantile mapping that
  #cite(<belay_merging_2022>, form: "prose") recommend. Here it is applied only to FY4B (§2.5). On the
  RQ1 storm days corrected FY4B's totals come out drier (§3.4); the effect on hourly peaks is not
  assessed.
- *Generality.* One region and two and a half years; the transfer to other climates and denser or
  sparser networks is untested.

== How others have studied these problems

- Validation of IMERG and GSMaP in complex terrain, and by rain intensity @tang_comparative_2020
  @tang_evaluation_2023 @wang_comparison_2023 @guo_has_2024. #cite(<zhou_capability_2021>, form: "prose")
  scored IMERG V6 Early, Late and Final on four Hubei storms, one of them in Shiyan in June 2016. These
  studies score IMERG Final or its near-real-time runs, and #cite(<tang_evaluation_2023>, form: "prose") the
  gauge-adjusted GSMaP_Gauge; here neither product is adjusted to gauges, and the error is split
  into placement (RQ1) and hourly magnitude by class (RQ2). In a mountain catchment of south-west China
  the accuracy of IMERG and GSMaP falls with elevation and slope @tang_evaluation_2023, and across
  mainland China IMERG's errors change with elevation @guo_has_2024. Here class biases barely differ
  between relief classes (§3.2), perhaps because local relief is not elevation, and only two relief
  classes have enough gauges.
- Satellite–gauge merging by residual interpolation (kriging, IDW, optimal interpolation)
  @xie_conceptual_2011 @belay_merging_2022, geographically weighted regression
  @brunsdon_geographically_1996, used for merging by #cite(<hu_satellite_2015>, form: "prose"), and machine
  learning @shen_merging_2022. This study
  compares them on held-out gauges and asks where, not only whether, the satellite adds skill.
- Multiscale GWR, which lets each covariate have its own bandwidth @oshan_mgwr_2019.
  #cite(<tu_spatial_2024>, form: "prose") use it to downscale annual and monthly IMERG; here it corrects
  hourly rain, as the correction step of the fusion.
- Spatial rather than random cross-validation: for large-scale mapping models,
  #cite(<ploton_spatial_2020>, form: "prose") show that random splits suggest skill that spatial splits
  remove. Here random splits would leave a held-out gauge next to a training one and overstate skill far
  from the network; the balanced spatial folds used here measure that case directly.
- Satellite rainfall as input to flood and rainfall–runoff models @saouabe_evaluation_2020
  @seo_hydrologic_2023, which motivates the outlook below.

== Outlook: flood forecasting

On hours when the gauge records 8 mm h⁻¹ or more, raw satellite rain is #stat("rq2_classes", "loss_heavy")% below it (RQ2), and
raw GPM reaches those rates only #val(mu, "FBI8", product: G, method: "raw") times as often as the gauges. Flood models are most
sensitive to these hours, so raw products should not drive them. Fusion lowers the error on these
hours but does not restore the peaks (RQ3). A follow-up outside RQ1–RQ3: aggregate raw, gauge-only and fused fields to basin-mean
rainfall for the gauged catchments during the six RQ1 events and the wettest seasons, correlate it (and
antecedent totals) with observed discharge, and only if fused rain correlates clearly better, drive a
simple rainfall–runoff model as a demonstration. Basin means also sidestep the point-versus-pixel
mismatch. This needs discharge data, which the project does not
yet have.

= Research plan: what remains

The study is finished when RQ1–RQ3 have answered the questions listed under each in §3. Progress
reports report against these three research contents.

*Research content.*
- RQ1 (space): complete.
- RQ2 (time): complete.
- RQ3 (fusion): complete.
- FY4B: in RQ1–RQ3 as the third product, on the matched sample; §3.4 complete; corrected FY4B stays
  out of the merge and the agreement blend.
- Literature (§4.2): the citations come from the 18 papers with a PDF in Zotero and the IMERG V07 ATBD.
  Claude checked each cited sentence against the paper's full text; the author has still to read
  them. The other 14 GGGG items without a PDF are unread, among them the MGWR paper of Fotheringham et
  al. (2017), the two machine-learning merging papers and two spatial cross-validation papers.

*Paper preparation* (presentation only; no research question depends on it).
- Fig. 1: add a location inset and the river basins, which need boundary and basin data not yet in
  `data/`.
- Fig. 4: add the landform-region panel; the result is already known (§3.2).
- Figures planned but not yet in this plan: Fig. 3, spatial scores per event and product; Fig. 5, rain
  events by class (detection, timing and volume); Fig. 7, where fusion helps (by intensity, by distance
  to the nearest gauge, and in CSI); Fig. 8, raw vs fused skill per anchor, and the merge weights.
- Move the diurnal-cycle results to supplementary material.

*Outside RQ1–RQ3.* Obtain hydrological-station discharge data for the flood-forecasting outlook (§4.3).

#bibliography("references.bib", title: "References", style: "apa")

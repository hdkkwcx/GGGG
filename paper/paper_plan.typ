// Research plan: gauge, GPM and GSMaP precipitation in a sparse-gauge mountain area.
// Compile from the repository root:  typst compile paper/paper_plan.typ
// Figures are copies of the `--gpm-gsmap` outputs of scripts/plot_heavy_rain_events.py,
// scripts/plot_satellite_temporal_evaluation.py and scripts/plot_gauge_satellite_fusion.py
// (the fusion figures from its FY4B-free IMERG Uncal report, output/gauge_satellite_fusion_nofy4b_gpmuncal).
// Numbers are read from the CSVs beside those figures in output/.
// Fig. 1 is scripts/plot_study_area_map.py's output/study_area_map/figures/fig1_study_area.png, drawn
// straight from data/ (its plot_provenance.csv records the inputs' SHA-256).
// Table 3 and Figs. 6-7 come from a clean benchmark run at commit 94afc78 (git_dirty=false, recorded
// in its benchmark_scope.csv with its Manifest.toml and a SHA-256 for every input file):
// ..._fusedanchor_nofy4b_gpmuncal. It reproduces the earlier runs at ec320ab and on uncommitted code
// (ef3616d) byte for byte.
// RQ1/RQ2 and Figs. 2-4 come from clean IMERG Uncal runs at commit 3bce549, plotted at 072eef6:
// output/heavy_rain_events_gpm_gsmap_gpmuncal, heavy_core_diagnostics_gpmuncal,
// satellite_temporal_evaluation_nofy4b_gpmuncal and, for the Scope paragraph's FY4B numbers,
// satellite_temporal_evaluation_gpmuncal.
// scripts/make_paper.jl rebuilds every figure and number here in dependency order; the [[step]]
// graph in paper/sources.toml records each step's command and what it needs.

#set document(
  title: "Research plan: satellite precipitation and gauge fusion (GPM, GSMaP)",
  date: datetime(year: 2026, month: 10, day: 23),
)
#set page(
  paper: "a4",
  margin: (x: 2.1cm, top: 2.3cm, bottom: 2.2cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 8pt, fill: luma(110))
      Research plan — GPM and GSMaP
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
// A literature slot still to be filled with a paper actually read.
#let refslot = text(fill: rgb("#c0392b"))[\[ref\]]

// Every result the text cites is read from paper/numbers/*.json, never typed by hand:
// scripts/export_paper_numbers.py exports them from the source CSVs named in paper/sources.toml, and
// scripts/check_paper_provenance.py fails if a JSON no longer matches its source. A missing row
// prints a red "?".
#let load(name) = json("numbers/" + name + ".json").rows
#let mu = load("methods_uncal")   // pooled metrics, IMERG Uncal run (Table 3)
#let cu = load("vs_adw_uncal")    // bootstrap comparisons with ADW, by stratum
#let su = load("strata_uncal")    // gain over ADW and mean bias by rain intensity and distance to a gauge
#let r1 = load("rq1_summary")      // RQ1: mean over the six events, IMERG Uncal run
#let e1 = load("rq1_events")       // RQ1: scores per event
#let k1 = load("rq1_cores")        // RQ1: storm-core deficit per event and product
#let c2 = load("rq2_classes")      // RQ2: by intensity class, season and relief
#let v2 = load("rq2_events")       // RQ2: rain-event detection and timing
#let a2 = load("rq2_areal")        // RQ2: areal-mean series
#let n0 = load("study_gauges")     // Study area: gauges RQ2 scores, in all and by landform region
#let fr = load("scope_fy4b_r")     // Scope: FY4B on its own hours, IMERG Uncal run with FY4B
#let fa = load("scope_fy4b_areal")
#let fp = load("scope_fy4b_pod")
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

/// Mean bias (mm h⁻¹) of one anchor and method on hours of ≥ 20 mm h⁻¹ at the held-out gauges.
#let bias20(product, method) = val(su, "bias", product: product, method: method,
  group: "rain_intensity", level: "severe_rainstorm")

#align(center)[
  #text(size: 17pt, weight: "semibold", font: "Arial")[Research plan]
  #v(0.2em)
  #text(size: 11.5pt, fill: luma(60))[Correcting GPM and GSMaP satellite precipitation against rain gauges with GWR / MGWR]
]

#note[
  *Problem.* Flood forecasting in mountain catchments needs hourly heavy rain, which a sparse gauge
  network cannot resolve. Satellite products cover every pixel, but how their error depends on place,
  intensity and hour is poorly known, and so is what fusing them with the gauges adds far from the
  gauges and on heavy hours.

  *Objective.* To correct GPM and GSMaP satellite precipitation against rain gauges with GWR / MGWR.
  RQ1 and RQ2 measure the errors to correct, in space and in time; RQ3 tests the correction. Three
  research contents, one results section each (§3).
  #v(0.3em)
  #table(
    columns: (auto, 1.4fr, 1.8fr, 0.8fr, 1.1fr),
    align: left,
    stroke: none,
    rule,
    table.header([*RQ*], [*Question*], [*Data and method*], [*Status*], [*To finish*]),
    thin,
    [1], [Space: does the satellite put heavy rain in the right place, at the right amount?],
    [Six heavy-rain days; event totals at the gauges, scored by spatial r, bias, RMSE and CSI at 50 mm],
    [Answered], [—],
    [2], [Time: does it follow the hourly evolution of rain at each intensity?],
    [Gauge-wet hours in five intensity classes, by season and relief; rain events; areal means],
    [Answered], [—],
    [3], [Fusion: what does fusing it with the gauges gain, with which method and anchor?],
    [Gauge-only and fusion methods on GPM, GSMaP and their merge, scored at held-out gauges by spatial
    cross-validation],
    [Answered], [—],
    rule,
  )
]

= Introduction

Flood forecasting in mountain catchments needs hourly rainfall, and above all the moderate and heavy
rain that drives runoff. Rain gauges measure it directly but are sparse in the mountains, so large
parts of a catchment sit tens of kilometres from the nearest gauge. Satellite products such as GPM
IMERG and GSMaP cover every pixel every hour, but their errors depend on rain type and intensity.

Most evaluations of these products work at daily or monthly scale and report pooled scores. Pooled
scores do not say whether a product puts a heavy-rain event in the right place, whether it follows the
hour-by-hour evolution of rain of a given intensity, or where its error comes from. Gauge–satellite
fusion is widely used to correct the products, but it is usually judged by average skill. What fusion
adds *far from the gauges* and *on heavy hours*, which is exactly the sparse-network flood case, is less
clear.

*Objective.* To correct GPM and GSMaP satellite precipitation against rain gauges with GWR / MGWR.
RQ1 and RQ2 measure the errors to correct, in space and in time; RQ3 tests the correction. §3 takes
each as a research content, with the questions it must answer:

- *RQ1 — Space.* Does the satellite field reproduce *where* heavy rain falls?
- *RQ2 — Time.* Does it reproduce *how rain evolves in time* at each intensity?
- *RQ3 — Fusion.* What does gauge–satellite fusion gain, which fusion method and which anchor product
  are best, and does fusion fix the moderate and heavy rain that floods depend on?

*Scope.* Two products, GPM IMERG and GSMaP. FY4B is left out: on its own hours it scores worst even
when no product is gauge-calibrated (hourly r #val(fr, "r", product: "FY4B") vs #stat("scope_fy4b_r", "others"), areal-mean r #val(fa, "r", product: "FY4B") vs #stat("scope_fy4b_areal", "others"),
wet-hour detection #val(fp, "POD", product: "FY4B") vs #stat("scope_fy4b_pod", "others")), so neither the merged product nor the agreement blend in
RQ3 uses it.

= Data and Methods

== Study area and data

#val(n0, "n", region: "all") hourly gauges in the Shiyan study area, 2022–2024 (values quantised to 0.5 mm). All sit in hills
and mountains: by local relief from the Copernicus GLO-30 DEM, #val(n0, "n", region: "relief_lt500") are below 500 m, #val(n0, "n", region: "relief_500_1000") at 500–1000 m
and #val(n0, "n", region: "relief_ge1000") at 1000 m or more. GPM IMERG V07 Uncal and GSMaP v8 operational (`hourlyPrecipRate`), both at
0.1°, are sampled at the gauge pixels: IMERG Uncal from NASA GES DISC's half-hourly files, averaged to
hours, and GSMaP through Google Earth Engine. Neither carries gauge information: IMERG Uncal is the
multi-satellite estimate before IMERG's monthly gauge calibration, and GSMaP operational has no gauge
correction. Days run 08–08 Beijing time.

#figure(
  image("figures/fig01_study_area.png", width: 80%),
  caption: [Draft Fig. 1. The study area: shaded relief from a 1 km SRTM DEM, and the hourly gauges by
  landform region, classed by local relief in the Copernicus GLO-30 DEM. Hollow markers: the region has
  too few gauges to be interpreted on its own.],
)

== Heavy-rain events and spatial scores (RQ1)

Heavy-rain days are 08–08 days on which more than three gauges exceed 50 mm. Of 56 such days, six are
selected, one per rain system: four widespread, two localized. At each gauge the satellite 24 h total
is compared with the gauge total by spatial correlation, spatial CV ratio, rain-centre shift, mean
bias, RMSE and CSI at 50 mm.

== Hourly evaluation by intensity (RQ2)

Only gauge-wet hours are scored, in five intensity classes (0.1–2, 2–4, 4–8, 8–20 and ≥ 20 mm h⁻¹),
stratified by season and landform region, with day-block bootstrap intervals. Rain events are runs of
wet hours separated by at least 3 dry hours and are scored for detection, timing and volume. Areal-mean
series are compared at hourly and daily scale.

== Gauge–satellite fusion and its validation (RQ3)

Gauge-only interpolators (IDW, ADW, TPS, GWR) are the baseline. Fusion methods use a satellite product
as the anchor and correct it with the gauges (residual GWR, mixed GWR, MGWR), optionally blended with
ADW where the satellite reports rain. The agreement blend chooses the ADW weight separately for each
band of how many of the two products call the cell wet (≥ 0.1 mm h⁻¹) and the larger of their two
rates (edges 0.1, 2.5 and 8 mm h⁻¹). MGWR with the agreement blend is the headline fusion method.
Anchors are GPM (IMERG Uncal), GSMaP and an OLS merge of the two.
Every method is scored at gauges held out of the fit, by balanced spatial 5-fold cross-validation over all
26,304 hours of 2022–2024, using RMSE, correlation and CSI, overall and by intensity and distance to
the nearest training gauge.

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
- *Placement is similar for both products, and fair.* Mean spatial r over the six events is
  #val(r1, "r", product: G) for GPM and #val(r1, "r", product: "GSMaP") for GSMaP, #stat("rq1_events", "r_four") on four events. Both fail on 2022-06-26
  (r ≤ #stat("rq1_events", "r_0626_max")), when the rain moved north to south and left a nearly flat gauge total that the
  products' north–south gradient does not match.
- *Amounts differ sharply.* GPM is somewhat dry (MB #val(r1, "MB", product: G) mm, RMSE #val(r1, "RMSE", product: G) mm), GSMaP wet
  (MB #val(r1, "MB", product: "GSMaP") mm, RMSE #val(r1, "RMSE", product: "GSMaP") mm).
- *Localized storm cores are missed by both* (CSI at 50 mm on the two localized events: GPM
  #val(e1, "CSI50", product: G, event_day: "2023-09-10") and #val(e1, "CSI50", product: G, event_day: "2024-07-30"), GSMaP #val(e1, "CSI50", product: "GSMaP", event_day: "2023-09-10") and #val(e1, "CSI50", product: "GSMaP", event_day: "2024-07-30")). Pixel averaging explains only
  #stat("rq1_cores", "repr_local")% of the core deficit; the rest is the products spreading short bursts over more hours, at
  #stat("rq1_cores", "peak_ratio_local") of the peak rate.

#figure(
  image("figures/fig02_event_maps.png", width: 100%),
  caption: [Draft Fig. 2. The six heavy-rain events (columns). Top row: gauge 24 h totals (08–08 BJT).
  Lower rows: satellite minus gauge at each gauge, from the pixel sampled there; red means the satellite is
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
- *Underestimation grows with intensity, in every season.* GPM relative bias runs from #val(c2, "RB0", product: G, season: "all", region: "all", gauge_class: "light")% on light
  rain to #val(c2, "RB0", product: G, season: "all", region: "all", gauge_class: "severe_rainstorm")% at ≥ 20 mm h⁻¹, GSMaP's from #val(c2, "RB0", product: "GSMaP", season: "all", region: "all", gauge_class: "light")% to #val(c2, "RB0", product: "GSMaP", season: "all", region: "all", gauge_class: "severe_rainstorm")%; GPM is the drier of the two
  in every class. At ≥ 8 mm h⁻¹ at most #stat("rq2_classes", "hit_heavy_max")% of hours land in the right class.
- *Detection and timing are not the problem; magnitude is.* Rain events peaking at 2–4 mm h⁻¹ or more
  are detected #stat("rq2_events", "pod")% of the time, and median rain-centre timing errors stay between
  #stat("rq2_events", "centroid_min") and #stat("rq2_events", "centroid_max") h.
- *Season matters, terrain does not.* Both products fail in winter (bias #val(c2, "RB0", product: G, season: "DJF", region: "all", gauge_class: "all_wet")% for GPM and
  #val(c2, "RB0", product: "GSMaP", season: "DJF", region: "all", gauge_class: "all_wet")% for GSMaP), while between relief below 500 m and 500–1000 m the class biases differ
  by at most #stat("rq2_classes", "relief_gap_max") percentage points.
- *Areal means are much better than points* (hourly r #val(a2, "r", product: G, timescale_h: "1") for GPM and #val(a2, "r", product: "GSMaP", timescale_h: "1") for GSMaP).

#figure(
  image("figures/fig04a_intensity_season.png", width: 100%),
  caption: [Draft Fig. 4. Hourly skill by gauge intensity class and season.],
)

*To finish.* Nothing: the question is answered.

== RQ3 — What does gauge fusion gain, with which method and which product?

*Questions to answer.*
- Does fusion beat both the raw product and gauges alone, overall and on heavy hours?
- How does the gain change with distance to the nearest gauge?
- Which fusion method and which anchor product are best?

*Findings so far.*
- *Fusion beats both of its inputs on GPM, not on GSMaP.* MGWR followed by a blend with ADW where the
  satellites report rain cuts RMSE by #pct(mu, "vs_raw", product: G, method: B) from raw GPM and #pct(mu, "vs_raw", product: "GSMaP", method: B) from raw GSMaP. Against gauge-only ADW it
  gains #pct(mu, "vs_adw", product: G, method: B) on GPM and #pct(mu, "vs_adw", product: "MERGED_OLS", method: B) on the merged anchor, both significant, but only
  #pct(mu, "vs_adw", product: "GSMaP", method: B) on GSMaP, which is not (p = #val(cu, "p_holm", product: "GSMaP", method: B, stratum: "overall")).
- *The gain holds on heavy hours and grows far from the gauges, but the heavy-rain volume is not
  restored.* On heavy hours
  (≥ 8 mm h⁻¹) MGWR alone improves on ADW by #pct(cu, "change", product: G, method: "mgwr", stratum: "8_plus") (GPM) and #pct(cu, "change", product: "GSMaP", method: "mgwr", stratum: "8_plus") (GSMaP), the agreement
  blend by #pct(cu, "change", product: G, method: B, stratum: "8_plus") and #pct(cu, "change", product: "GSMaP", method: B, stratum: "8_plus"), all significant. Unblended MGWR *hurts* dry, light and moderate hours,
  which is why the blend is needed. On GPM the improvement over ADW is #pct(su, "change", product: G, method: B, group: "nearest_train_km", level: "0_20") within 20 km of a training gauge,
  #pct(su, "change", product: G, method: B, group: "nearest_train_km", level: "20_50") at 20–50 km and #pct(su, "change", product: G, method: B, group: "nearest_train_km", level: "50_100") at 50–100 km. At ≥ 20 mm h⁻¹ the mean bias of fused GPM is #bias20(G, B) mm h⁻¹, against
  #bias20(G, "raw") for raw GPM and #bias20(G, "adw") for gauge-only ADW. On GSMaP the fused field (#bias20("GSMaP", B)) is drier there than the
  raw product (#bias20("GSMaP", "raw")).
- *Best method and product.* Among gauge-only interpolators ADW is best, with IDW indistinguishable
  (#val(mu, "RMSE", product: G, method: "adw") vs #val(mu, "RMSE", product: G, method: "idw") mm h⁻¹). The best fusion is MGWR + blend. GPM is the
  better anchor, with neither product gauge-calibrated (raw RMSE #val(mu, "RMSE", product: G, method: "raw") vs #val(mu, "RMSE", product: "GSMaP", method: "raw") mm h⁻¹), and only
  fused GPM beats ADW significantly. Merging the two gives a better raw field (RMSE #val(mu, "RMSE", product: "MERGED_OLS", method: "raw")) but adds
  almost nothing after fusion (#val(mu, "RMSE", product: "MERGED_OLS", method: B) vs #val(mu, "RMSE", product: G, method: B)).

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
    [—], [ADW (gauges only)], ..t3(mu, "GPM_uncal", "adw", raw: false),
    [—], [IDW (gauges only)], ..t3(mu, "GPM_uncal", "idw", raw: false),
    [—], [TPS (gauges only)], ..t3(mu, "GPM_uncal", "tps", raw: false),
    [—], [GWR (gauges only)], ..t3(mu, "GPM_uncal", "gwr", raw: false),
    thin,
    [GPM], [raw], ..t3(mu, "GPM_uncal", "raw"),
    [GPM], [MGWR + agreement blend], ..t3(mu, "GPM_uncal", "blend_agrenv_mgwr"),
    [GSMaP], [raw], ..t3(mu, "GSMaP", "raw"),
    [GSMaP], [MGWR + agreement blend], ..t3(mu, "GSMaP", "blend_agrenv_mgwr"),
    [Merged (OLS)], [raw], ..t3(mu, "MERGED_OLS", "raw"),
    [Merged (OLS)], [MGWR + agreement blend], ..t3(mu, "MERGED_OLS", "blend_agrenv_mgwr"),
    rule,
  ),
  caption: [Draft Table 3. Held-out gauges, balanced spatial CV over all 26,304 hours of 2022–2024;
  RMSE in mm h⁻¹; improvements are RMSE reductions; CSI at 2 mm h⁻¹. GPM is IMERG Uncal. No row uses
  FY4B: the merge and the agreement blend are built from GPM and GSMaP only. Gauge-only rows are on
  GPM's evaluation mask. The run is from commit `94afc78`, with no uncommitted changes.],
)

*To finish.* Nothing: the question is answered.

= Discussion

== Uncertainties

- *Point versus pixel.* A gauge samples a point, the products a 0.1° pixel. Where rain varies below
  10 km (2022-06-26, the localized cores) part of the "error" is representativeness, not product skill;
  pixel averaging explains #stat("rq1_cores", "repr_local")% of the core deficit, and the rest comes from the products spreading
  short bursts over too many hours.
- *Gauge resolution.* The 0.5 mm quantisation limits what can be said about the lightest class.
- *Small samples.* Only #val(n0, "n", region: "relief_ge1000") gauges lie at relief ≥ 1000 m, too few to interpret; heavy hours and the six
  events are few, so their intervals are wide.
- *Residual heavy-rain deficit.* Fusion lowers the RMSE of heavy hours but restores almost none of
  their missing volume (§3.3). Gauge-only ADW is as low as raw GPM there, so neither input carries the
  peak to a held-out gauge.
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

RQ2 shows that raw satellite rain loses #stat("rq2_classes", "loss_heavy")% of the volume at 8 mm h⁻¹ and above, the input
flood models are most sensitive to, so raw products should not drive them; fusion lowers the error on
these hours but recovers little of the missing volume (RQ3). A follow-up outside RQ1–RQ3: aggregate raw, gauge-only and fused fields to basin-mean
rainfall for the gauged catchments during the six RQ1 events and the wettest seasons, correlate it (and
antecedent totals) with observed discharge, and only if fused rain correlates clearly better, drive a
simple rainfall–runoff model as a demonstration. This needs discharge data, which the project does not
yet have.

= Research plan: what remains

The study is finished when RQ1–RQ3 have answered the questions listed under each in §3. Progress
reports report against these three research contents.

*Research content.*
- RQ1 (space): complete.
- RQ2 (time): complete.
- RQ3 (fusion): complete.
- Literature (§4.2): fill the #refslot slots with papers actually read, so the Discussion can set each
  result against earlier work.

*Paper preparation* (presentation only; no research question depends on it).
- Fig. 1: add a location inset and the river basins, which need boundary and basin data not yet in
  `data/`.
- Fig. 4: add the landform-region panel; the result is already known (§3.2).
- Figures planned but not yet in this plan: Fig. 3, spatial scores per event and product; Fig. 5, rain
  events by class (detection, timing and volume); Fig. 7, where fusion helps (by intensity, by distance
  to the nearest gauge, and in CSI); Fig. 8, raw vs fused skill per anchor, and the merge weights.
- Move the diurnal-cycle results to supplementary material.

*Outside RQ1–RQ3.* Obtain hydrological-station discharge data for the flood-forecasting outlook (§4.3).

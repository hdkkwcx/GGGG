// Research plan: gauge, GPM and GSMaP precipitation in a sparse-gauge mountain area.
// Compile from the repository root:  typst compile paper/paper_plan.typ
// Figures are copies of the `--gpm-gsmap` outputs of scripts/plot_heavy_rain_events.py,
// scripts/plot_satellite_temporal_evaluation.py and scripts/plot_gauge_satellite_fusion.py
// (the fusion figures from its FY4B-free IMERG Uncal report, output/gauge_satellite_fusion_nofy4b_gpmuncal).
// Numbers are read from the CSVs beside those figures in output/.
// Table 3 and Figs. 6-7 come from clean benchmark runs at commit ec320ab (git_dirty=false, recorded
// in each run's benchmark_scope.csv with its Manifest.toml): ..._fusedanchor_nofy4b_gpmuncal and,
// for the GPM (Final) rows, ..._fusedanchor_nofy4b. They reproduce the earlier uncommitted-code
// runs (output/*_uncommitted_ef3616d) byte for byte.

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

// Every result the text cites is read from paper/numbers/*.json, never typed by hand:
// scripts/export_paper_numbers.py exports them from the source CSVs named in paper/sources.toml, and
// scripts/check_paper_provenance.py fails if a JSON no longer matches its source. A missing row
// prints a red "?".
#let load(name) = json("numbers/" + name + ".json").rows
#let mu = load("methods_uncal")   // pooled metrics, IMERG Uncal run (Table 3)
#let mf = load("methods_final")   // pooled metrics, IMERG Final sensitivity run
#let cu = load("vs_adw_uncal")    // bootstrap comparisons with ADW, by stratum
#let cf = load("vs_adw_final")
#let su = load("strata_uncal")    // gain over ADW by rain intensity and distance to a gauge
#let sf = load("strata_final")
#let wu = load("merge_uncal")     // OLS merge weights
#let wf = load("merge_final")
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
/// One Table 3 row's numbers: RMSE, r, CSI 2, vs raw (a dash on gauge-only rows) and vs ADW.
#let t3(rows, product, method, raw: true) = (
  val(rows, "RMSE", product: product, method: method),
  val(rows, "r", product: product, method: method),
  val(rows, "CSI2", product: product, method: method),
  if raw [#val(rows, "vs_raw", product: product, method: method)%] else [—],
  [#val(rows, "vs_adw", product: product, method: method)%],
)

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
    [MGWR + blend beats raw GPM by #pct(mu, "vs_raw0", product: G, method: B) and gauges alone by #pct(mu, "vs_adw0", product: G, method: B); gain largest on heavy hours and far from gauges. GPM is the better anchor, with or without gauge calibration; fused GSMaP does not beat gauges alone.],
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
  satellites report rain cuts RMSE by #pct(mu, "vs_raw", product: G, method: B) from raw GPM and #pct(mu, "vs_raw", product: "GSMaP", method: B) from raw GSMaP. Against gauge-only ADW it
  gains #pct(mu, "vs_adw", product: G, method: B) on GPM and #pct(mu, "vs_adw", product: "MERGED_OLS", method: B) on the merged anchor, both significant, but only
  #pct(mu, "vs_adw", product: "GSMaP", method: B) on GSMaP, which is not (p = #val(cu, "p_holm", product: "GSMaP", method: B, stratum: "overall")). Every blended variant on GPM and on the merged anchor beats ADW significantly; over
  all hours none on GSMaP does.
- *The gain is where floods need it.* On heavy hours (≥ 8 mm h⁻¹) MGWR alone improves on ADW by
  #pct(cu, "change", product: G, method: "mgwr", stratum: "8_plus") (GPM) and #pct(cu, "change", product: "GSMaP", method: "mgwr", stratum: "8_plus") (GSMaP), the agreement blend by
  #pct(cu, "change", product: G, method: B, stratum: "8_plus") and #pct(cu, "change", product: "GSMaP", method: B, stratum: "8_plus"), all significant. Unblended regression
  *hurts* dry, light and moderate hours (MGWR on GPM #pct(su, "change0", product: G, method: "mgwr", group: "rain_intensity", level: "no_rain"), #pct(su, "change0", product: G, method: "mgwr", group: "rain_intensity", level: "light") and #pct(su, "change0", product: G, method: "mgwr", group: "rain_intensity", level: "moderate") against ADW), which is why the
  blend is needed. At 2.5–8 mm h⁻¹ the blend gains nothing significant (#val(cu, "change", product: G, method: B, stratum: "2.5_8")% on GPM).
- *The gain grows with distance from the network.* On GPM the improvement over ADW is
  #pct(su, "change", product: G, method: B, group: "nearest_train_km", level: "0_20") within 20 km of a training gauge, #pct(su, "change", product: G, method: B, group: "nearest_train_km", level: "20_50") at 20–50 km and
  #pct(su, "change", product: G, method: B, group: "nearest_train_km", level: "50_100") at 50–100 km (GSMaP #pct(su, "change", product: "GSMaP", method: B, group: "nearest_train_km", level: "0_20"), #pct(su, "change", product: "GSMaP", method: B, group: "nearest_train_km", level: "20_50"), #pct(su, "change", product: "GSMaP", method: B, group: "nearest_train_km", level: "50_100")).
- *Best method.* Among gauge-only interpolators ADW, with IDW indistinguishable (#val(mu, "RMSE", product: G, method: "adw") vs #val(mu, "RMSE", product: G, method: "idw") mm h⁻¹);
  TPS and GWR are #mag(mu, "vs_adw", product: G, method: "tps") and #mag(mu, "vs_adw", product: G, method: "gwr") worse. Best fusion: MGWR + blend; MGWR alone only for heavy hours.
- *Best product.* GPM, with neither product gauge-calibrated: lower raw error over all 2022–2024 hours
  (RMSE #val(mu, "RMSE", product: G, method: "raw") vs #val(mu, "RMSE", product: "GSMaP", method: "raw") mm h⁻¹,
  r #val(mu, "r", product: G, method: "raw") vs #val(mu, "r", product: "GSMaP", method: "raw")). GPM's near-zero event bias and better areal timeline
  (RQ1, RQ2) are so far measured on IMERG Final. Once fused, the two anchors are close
  (#val(mu, "RMSE", product: G, method: B) vs #val(mu, "RMSE", product: "GSMaP", method: B)), but only fused GPM beats ADW significantly. Merging gives a better raw field
  (RMSE #val(mu, "RMSE", product: "MERGED_OLS", method: "raw")) but adds almost nothing after fusion (#val(mu, "RMSE", product: "MERGED_OLS", method: B) vs #val(mu, "RMSE", product: G, method: B)); its OLS weights are
  #val(wu, "beta_gpm", product: "MERGED_OLS") for GPM and #val(wu, "beta_gsmap", product: "MERGED_OLS") for GSMaP, plus a #val(wu, "intercept", product: "MERGED_OLS") mm h⁻¹ intercept.
- *IMERG Final gives the same fused result.* Gauge-calibrated IMERG Final has higher raw error than
  IMERG Uncal (RMSE #val(mf, "RMSE", product: "GPM", method: "raw") vs #val(mu, "RMSE", product: G, method: "raw") mm h⁻¹), but after fusion the two are level
  (#val(mf, "RMSE", product: "GPM", method: B) vs #val(mu, "RMSE", product: G, method: B); #val(mf, "vs_adw", product: "GPM", method: B)% vs #val(mu, "vs_adw", product: G, method: B)% over ADW). On Final, MGWR alone loses more to ADW
  (#pct(cf, "change", product: "GPM", method: "mgwr", stratum: "overall"), against a non-significant #pct(cu, "change", product: G, method: "mgwr", stratum: "overall") on Uncal) and gains more on heavy hours
  (#pct(cf, "change", product: "GPM", method: "mgwr", stratum: "8_plus") vs #pct(cu, "change", product: G, method: "mgwr", stratum: "8_plus")). The OLS merge puts less weight on Final than on Uncal
  (#val(wf, "beta_gpm", product: "MERGED_OLS") vs #val(wu, "beta_gpm", product: "MERGED_OLS")).

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
    thin,
    [GPM (Final)], [raw], ..t3(mf, "GPM", "raw"),
    [GPM (Final)], [MGWR + agreement blend], ..t3(mf, "GPM", "blend_agrenv_mgwr"),
    rule,
  ),
  caption: [Draft Table 3. Held-out gauges, balanced spatial CV over all 26,304 hours of 2022–2024;
  RMSE in mm h⁻¹; improvements are RMSE reductions; CSI at 2 mm h⁻¹. GPM is IMERG Uncal. No row uses
  FY4B: the merge and the agreement blend are built from GPM and GSMaP only. The GPM (Final) rows are a
  sensitivity run with IMERG Final in IMERG Uncal's place, otherwise identical. Gauge-only rows are on
  GPM's evaluation mask. Both runs are from commit `ec320ab`, with no uncommitted changes.],
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

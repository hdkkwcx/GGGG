// Why FY4B scored worst of the three satellite products: the evidence in this repository, in one note.
// Compile from the repository root:  typst compile --root . NOWAY/why_fy4b_scored_worst.typ
// (`--root .` because the numbers are read from paper/numbers/, outside this folder.)
//
// Numbers. Every number the paper exports is read from paper/numbers/*.json with the paper's own
// helpers (paper/paper_plan.typ), never typed, so this note cannot drift from the paper. A missing row
// prints a red "?". The few numbers the paper does not export are typed, and each paragraph that holds
// one names the file it was read from: Projects_Pending_Handover/fy4b_quality_handover.md,
// Projects_Pending_Handover/fy4b_tries.md, or a CSV under output/ (fy4b_error_budget,
// fy4b_process_trace, fy4b_two_stage_benchmark_comparison).
// Nothing here is a new analysis: no run was made for this note.

#set document(
  title: "Why FY4B scored worst of the three products",
  date: datetime(year: 2026, month: 11, day: 16),
)
#set page(
  paper: "a4",
  margin: (x: 2.1cm, top: 2.3cm, bottom: 2.2cm),
  numbering: "1",
  header: context {
    if counter(page).get().first() > 1 [
      #set text(size: 8pt, fill: luma(110))
      Why FY4B scored worst
      #h(1fr) 16 November 2026
    ]
  },
)
#set text(font: ("Libertinus Serif", "New Computer Modern"), size: 10.5pt, lang: "en")
#set par(justify: true, leading: 0.62em, spacing: 1.05em)
#set heading(numbering: "1.")
#show heading: set text(font: "Arial", weight: "semibold")
#show heading.where(level: 1): it => { v(0.5em); it; v(0.2em) }
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
/// Where a typed number was read from.
#let src(body) = text(size: 8.5pt, fill: luma(90))[(#body)]

// The paper's helpers, with the path made absolute from the repository root.
#let numbers(name) = json("/paper/numbers/" + name + ".json")
#let load(name) = numbers(name).rows
#let find(rows, ..conditions) = {
  let hits = rows.filter(row => conditions.named().pairs().all(pair => row.at(pair.at(0)) == pair.at(1)))
  if hits.len() == 0 { none } else { hits.at(0) }
}
#let val(rows, column, ..conditions) = {
  let hit = find(rows, ..conditions)
  if hit == none { text(fill: red)[?] } else { hit.at(column).replace("-", "−") }
}
#let pct(rows, column, ..conditions) = {
  let hit = find(rows, ..conditions)
  if hit == none { text(fill: red)[?] } else [#hit.at(column).trim("+", at: start).replace("-", "−")%]
}
#let stat(name, key) = {
  let stats = numbers(name).at("stats", default: (:))
  if key in stats { stats.at(key).replace("-", "−") } else { text(fill: red)[?] }
}

#let mu = load("methods_uncal")        // Table 3 of the paper: pooled scores at held-out gauges
#let cu = load("vs_adw_uncal")         // paired comparisons with gauge-only ADW
#let c2 = load("rq2_classes")          // RQ2: gauge-wet hours by class and season
#let a2 = load("rq2_areal")            // RQ2: areal-mean series
#let r1 = load("rq1_summary")          // RQ1: mean over the six heavy-rain events
#let e1 = load("rq1_events")           // RQ1: scores per event
#let dh = load("data_hours")           // the matched sample's hours and days
#let fe = load("fy4b_events_summary")  // corrected FY4B on the events, mean
#let fv = load("fy4b_event_scores")    // corrected FY4B on the events, per event
#let fq = load("fy4b_quality")         // corrected FY4B against plain FY4B, GPM and GSMaP
#let ft = load("fy4b_time_holdout")    // the same correction with whole weeks held out
#let fc = load("fy4b_clt")             // July 2023 cloud-type pilot, ice-cloud hours
#let fcf = load("fy4b_clt_filter")     // the cloud-type filter against an amount-matched control
#let fg = load("fusion_configs")       // plain against corrected FY4B in the merges and bands
#let tr = load("fy4b_tries")           // every FY4B fix tried, in words
#let G = "GPM_uncal"
#let B = "blend_agrenv_mgwr"
#let three = ("FY4B", G, "GSMaP")

/// Corrected FY4B's change in RMSE against `baseline`, same method, with its 95% interval.
#let dq(baseline, method, scheme: "balanced_spatial") = {
  let row = (scheme: scheme, product: "FY4B_QM", baseline: baseline, method: method)
  [#val(fq, "change", ..row)% (#val(fq, "lo", ..row) to #val(fq, "hi", ..row)%)]
}
/// One anchor and method with corrected instead of plain FY4B in the merges and the bands.
#let dg(product, method) = {
  let row = (product: product, method: method)
  [#val(fg, "change", ..row)% (#val(fg, "lo", ..row) to #val(fg, "hi", ..row)%)]
}
#let dt(candidate) = {
  let row = (scheme: "balanced_spatial", candidate: candidate, baseline: "FY4B")
  [#val(ft, "change", ..row)% (#val(ft, "lo", ..row) to #val(ft, "hi", ..row)%)]
}
#let dtr(candidate) = {
  let row = (scheme: "balanced_spatial", candidate: candidate, baseline: "FY4B")
  [#val(ft, "dr", ..row) (#val(ft, "dr_lo", ..row) to #val(ft, "dr_hi", ..row))]
}
/// One scorecard row: the same column for FY4B, GPM and GSMaP.
#let across(rows, column, ..conditions) = three.map(product =>
  val(rows, column, product: product, ..conditions))
#let all-wet = (season: "all", region: "all", gauge_class: "all_wet")

#align(center)[
  #text(font: "Arial", size: 17pt, weight: "semibold")[Why FY4B scored worst of the three products]
  #v(-0.3em)
  #text(size: 11.5pt, fill: luma(60))[FY4B against GPM (IMERG Uncal) and GSMaP at the rain gauges, June 2022 – December 2024]
  #v(-0.2em)
  #text(size: 9.5pt, fill: luma(90))[16 November 2026 · a summary of results already in this repository; no new run]
]

#note[
  *In short.* On the #val(dh, "hours", sample: "all_products") hours that the gauges and all three
  products share, FY4B is last on every score in Table 1, and after fusion with the gauges its field
  is no better than the gauges alone (#val(mu, "vs_adw", product: "FY4B", method: B)%,
  p = #val(cu, "p_holm", product: "FY4B", method: B, stratum: "overall")). Six reasons follow:
  its rain rate holds one infrared channel's information (§2); it reports rain under cold cloud
  that is not raining (§3); it also misses about two thirds of the wet hours (§4); its record is short
  and broken (§5); a correction against gauges helps, but part of the gain does not last (§6); and
  the fusion finds almost nothing in it that GPM and GSMaP do not already give (§7).
  #stat("fy4b_tries", "n") registered fixes were tried, and none changed its place.
  "Worst" means worst on average: §8 lists where FY4B is not last, and what was never tested.
]

= The scores

#figure(
  table(
    columns: (1fr, auto, auto, auto),
    align: (left, right, right, right),
    stroke: none,
    rule,
    table.header([*Score*], [*FY4B*], [*GPM*], [*GSMaP*]),
    thin,
    table.cell(colspan: 4)[_Every hour, raw product at the gauge_],
    [RMSE, mm h⁻¹ (lower is better)], ..across(mu, "RMSE", method: "raw"),
    [Correlation r], ..across(mu, "r", method: "raw"),
    [CSI at 2 mm h⁻¹], ..across(mu, "CSI2", method: "raw"),
    thin,
    table.cell(colspan: 4)[_Hours when the gauge is wet_],
    [Share the product also calls wet, %], ..across(c2, "POD", ..all-wet).map(cell => cell.trim("+")),
    [Rain events detected, % (range over event classes)], stat("rq2_events", "pod_fy4b"),
    table.cell(colspan: 2, align: center, stat("rq2_events", "pod")),
    thin,
    table.cell(colspan: 4)[_Averages and storms_],
    [Hourly r of the mean over all gauges], ..across(a2, "r", timescale_h: "1"),
    [Six heavy-rain events: mean spatial r of the totals], ..across(r1, "r"),
    [Six heavy-rain events: mean CSI at 50 mm], ..across(fe, "CSI50", window: "fy4b_hours"),
    thin,
    table.cell(colspan: 4)[_After MGWR and the agreement blend, at held-out gauges_],
    [RMSE, mm h⁻¹ (lower is better)], ..across(mu, "RMSE", method: B),
    [Gain in RMSE over the gauges alone (ADW), %], ..across(mu, "vs_adw", method: B),
    [p of that gain (Holm)], ..across(cu, "p_holm", method: B, stratum: "overall"),
    rule,
  ),
  caption: [Table 1. FY4B, GPM (IMERG Uncal) and GSMaP on one matched sample: the gauge-hours of
  June 2022 – December 2024 at which the gauge and all three products have a value. Higher is better
  except for RMSE. The last block is scored on the gauge-hours of that sample at which every method
  has a prediction. FY4B is last in every row. All values are read from the paper's exported numbers
  (`paper/numbers/`).],
)

Two of these rows carry most of the story. FY4B calls
#pct(c2, "POD", product: "FY4B", ..all-wet) of the gauge-wet hours wet, about half the share of GPM
(#pct(c2, "POD", product: G, ..all-wet)) and GSMaP (#pct(c2, "POD", product: "GSMaP", ..all-wet)).
And averaging over all gauges, which cancels errors of place, lifts GPM's hourly r to
#val(a2, "r", product: G, timescale_h: "1") but FY4B's only to
#val(a2, "r", product: "FY4B", timescale_h: "1"). So FY4B's main error is not where it puts the rain.
It is whether it reports rain at all.

= One infrared channel's information

FY4B's satellite-only rain rate is, in effect, one number read from one infrared channel. The
provider's own slides describe it as a lookup from one infrared brightness temperature, matched
offline to microwave rain from polar orbiters. The files agree: within any 5° × 5° box, the 10.8 µm
brightness temperature alone explains 93–99.5% of the variance of the rain rate, and at 220 K or
colder 94–100% of pixels are wet #src[two AGRI L1 scans; handover, "QPE is a one-channel lookup per
region"].

An infrared channel sees the temperature of the cloud top. It does not see the rain below it. GPM and
GSMaP are built on passive-microwave estimates, which respond to the rain and ice in the column
itself, with infrared in a supporting part only (general background on the two products, not measured
here). That difference in the input is the root of everything below.

It also explains why the fixes ran out. A pixel shift moves the number, a quantile map rescales it,
and a threshold cuts its low end. None of them adds information about whether it is raining.

= Rain reported under cold cloud that is not raining

- *Most of FY4B's wet calls are false.* At 0.1 mm h⁻¹ its false-alarm ratio is 0.71, against 0.56 for
  GPM and 0.60 for GSMaP. False rain makes up 39.8% of its squared error, against 11.8% for GPM and
  13.6% for GSMaP #src[handover, error budget table; `output/fy4b_error_budget/budget.csv`].
- *The false rain comes in episodes.* In 398 hours (2.0% of all hours) FY4B is wet at half the gauges
  or more while at most 5% of the gauges are wet. Those hours hold 38.7% of its false-rain squared
  error. GPM has 94 such hours, holding 6.6% of its own. In the worst hour (label 2023-03-19 08 BJT)
  FY4B's mean over the network is 16.3 mm h⁻¹, wet at every gauge, while the gauges read 0.00 and GPM
  0.08. The two thresholds were chosen after a first look; under the other choices tried, the share
  runs from 27.6% to 48.4% #src[handover, "False rain comes in episodes";
  `output/fy4b_error_budget/gross_summary.csv`].
- *It sits under ice cloud, mixed with the real rain.* In the July 2023 pilot with FY4B's own
  cloud-type product, hours dominated by ice cloud held #pct(fc, "false_alarms", name: "ice") of
  corrected FY4B's false alarms but also #pct(fc, "hits", name: "ice") of its hits. Corrected FY4B
  called rain in 80% of those hours and GPM in 44% #src[handover, cloud-type pilot]. Its false-alarm
  ratio there was #val(fc, "far", name: "ice"), against #val(fc, "far_gpm", name: "ice") for GPM. So
  no cloud class can be removed without removing the hits: the filter raised CSI at 0.1 mm h⁻¹ by
  #val(fcf, "csi", candidate: "FY4B_QM_CLT") (#val(fcf, "csi_lo", candidate: "FY4B_QM_CLT") to
  #val(fcf, "csi_hi", candidate: "FY4B_QM_CLT")), about half the +0.01 it needed.
- *The false rain looks more like rain than the real rain does.* In six AGRI scans, at gauges where
  FY4B is wet, the false alarms have colder tops than the hits (median 223.6 K against 227.9 K), and
  no single channel feature separates the two (AUC 0.60 at best) #src[handover, channel pilot]. A
  one-channel lookup is most sure of rain exactly there. Six scans are a small pilot, so this is a
  weak result.
- *Its heavy values fall in the wrong hours.* Raw FY4B reaches 8 mm h⁻¹
  #val(mu, "FBI8", product: "FY4B", method: "raw") times as often as the gauges, yet its CSI at
  8 mm h⁻¹ is #val(fg, "csi8_plain", product: "FY4B", method: "raw"), against
  #val(fg, "csi8_plain", product: G, method: "raw") for GPM and
  #val(fg, "csi8_plain", product: "GSMaP", method: "raw") for GSMaP.

= Rain that is missed

- *Misses count nearly as much as false rain.* Corrected FY4B's mean squared error is 0.262 above
  GPM's. Split by cause, false rain adds 0.177 and misses add 0.139, while the amounts on its hits take
  0.054 off. That last term is partly because it has half as many hits, so it does not show that its
  amounts are better #src[handover, error budget; `output/fy4b_error_budget/budget_gap.csv`].
- *Removing the false rain would not be enough.* With every false-rain hour set to zero by an oracle
  (a bound, not a method), corrected FY4B's r is 0.386, still below GPM's 0.436 as it stands
  #src[`output/fy4b_error_budget/oracles.csv`].
- *Rain events are missed too.* Events peaking at 2–4 mm h⁻¹ or more are detected
  #stat("rq2_events", "pod_fy4b")% of the time by FY4B, against #stat("rq2_events", "pod")% by GPM and
  GSMaP.
- *Averaging does not rescue it.* For the daily mean over the network, r is 0.46 for FY4B, against
  0.80 for GSMaP and 0.89 for GPM #src[handover, "It is not mainly misplacement";
  `output/fy4b_error_budget/time_scales.csv`]. An error of place would average out. This one does not.
- *Winter is nearly blank.* In December–February FY4B calls
  #pct(c2, "POD", product: "FY4B", season: "DJF", region: "all", gauge_class: "all_wet") of gauge-wet
  hours wet. All three products fail in winter
  (GPM #pct(c2, "POD", product: G, season: "DJF", region: "all", gauge_class: "all_wet"),
  GSMaP #pct(c2, "POD", product: "GSMaP", season: "DJF", region: "all", gauge_class: "all_wet")),
  FY4B the most.

= A short and broken record

Of the 26,304 hours of 2022–2024, FY4B has a usable value on
#val(dh, "hours", sample: "all_products") #src[`output/fy4b_process_trace/hour_funnel.csv`]:

- 3,624 hours fall before its first file, on 1 June 2022;
- 792 hours fall in a gap with no files, 1 February – 4 March 2024, between the satellite's last hour
  at 133°E and its first at 105°E;
- 912 hours are the hour 22–23 Beijing time, whose 14:15 UTC file is never produced, so FY4B has no
  complete day and no daily total;
- 889 more hours have no file or lack one to three of their four 15-minute files;
- of the 20,087 complete hours, 269 have no value at any gauge.

GPM and GSMaP cover all 26,304. Three more limits come with the files:

- *A ceiling.* The files declare a valid range of 0–20 mm h⁻¹, and no value in the rebuilt series is
  above 20.0. GPM reaches 43.8 and GSMaP 95.7 on the same sample, which has 1,420 gauge-hours at
  20 mm h⁻¹ or more #src[`output/fy4b_error_budget/ceiling.csv`, balanced spatial].
- *The hourly files cannot be used.* FY4B's own 1-, 3- and 24-hour files are corrected with China's
  rain gauges, so they cannot be scored against those gauges. The series had to be rebuilt from the
  15-minute files.
- *Three regimes in 31 months.* Two versions of the processing at 133°E, switched in May 2023, and the
  move to 105°E in March 2024. FY4B's June–August total against the gauges was 1.76 times theirs in
  2022, 0.89 in 2023 and 0.59 in 2024 #src[handover, "The product changes between regimes"]. Each
  regime has one summer, so the software, the position and the weather cannot be told apart. Whatever
  the cause, FY4B's relation to the gauges was not stable from year to year.

= A correction that helps, but only in part

The corrected series (a pixel shift and a quantile map, fitted on the training gauges of each fold)
is the best that the #stat("fy4b_tries", "n") tries produced. Against plain FY4B its RMSE changes by #dq("FY4B", "raw") raw,
#dq("FY4B", "mgwr") after MGWR and #dq("FY4B", B) after MGWR and the agreement blend. Four things
limit that gain.

- *Part of it belongs to the weeks it was fitted in.* With whole weeks held out as well as gauges, the
  raw change is #dt("FY4B_QM_crossed") in place of #dt("FY4B_QM_gauges"), and r no longer rises: it
  changes by #dtr("FY4B_QM_crossed"). Most of the difference is a few false-rain episodes whose
  largest values the map learns to shrink. This is not a leak. It means the map corrects weeks it has
  seen better than weeks it has not.
- *It still does not reach GPM.* Against GPM's agreement blend, corrected FY4B's RMSE differs by
  #dq(G, B).
- *It shows only far from gauges.* In the random split, where a held-out gauge usually has a training
  gauge nearby, the change after MGWR is #dq("FY4B", "mgwr", scheme: "random").
- *It makes storms drier.* On the six heavy-rain events the mean bias goes from
  #val(fe, "MB", product: "FY4B", window: "fy4b_hours") to
  #val(fe, "MB", product: "FY4B_QM", window: "fy4b_hours") mm and mean CSI at 50 mm from
  #val(fe, "CSI50", product: "FY4B", window: "fy4b_hours") to
  #val(fe, "CSI50", product: "FY4B_QM", window: "fy4b_hours"); mean r barely moves
  (#val(fe, "r", product: "FY4B", window: "fy4b_hours") to
  #val(fe, "r", product: "FY4B_QM", window: "fy4b_hours")). On 2022-10-03, the one event where
  FY4B's pattern was the best of the three, CSI at 50 mm falls from
  #val(fv, "CSI50", product: "FY4B", event_day: "2022-10-03") to
  #val(fv, "CSI50", product: "FY4B_QM", event_day: "2022-10-03").

After fusion, corrected FY4B's gain over the gauges alone is
#val(mu, "vs_adw", product: "FY4B_QM", method: B)%
(p = #val(cu, "p_holm", product: "FY4B_QM", method: B, stratum: "overall")), still not significant.

= Almost nothing for the fusion to use

- *The merge gives it almost no weight.* In the least-squares merge of the three products, fitted on
  the training gauges of each fold, FY4B's weight is 0.019–0.023, against 0.46–0.49 for GPM and
  0.14–0.16 for GSMaP #src[`output/fy4b_process_trace/fusion_merge_weights.csv`, balanced spatial].
- *Correcting it first changes little.* With corrected FY4B in the merges and the agreement bands in
  place of plain FY4B, the blend on the least-squares merge changes by #dg("MERGED_OLS", B), on GPM by
  #dg(G, B) and on GSMaP by #dg("GSMaP", B). Only the equal-weight merge, where FY4B is given a third
  of the weight, gains: #dg("MERGED_MEAN", B), and its CSI at 8 mm h⁻¹ changes by
  #val(fg, "dcsi8", product: "MERGED_MEAN", method: B)
  (#val(fg, "dcsi8_lo", product: "MERGED_MEAN", method: B) to
  #val(fg, "dcsi8_hi", product: "MERGED_MEAN", method: B)).
- *Its rain calls that no other product shares are, on balance, wrong.* Keeping FY4B only where GPM or
  GSMaP also reports rain, with no gauge involved, cuts its raw RMSE by 5.1% and its false-alarm ratio
  from 0.71 to 0.57. But the fused fields do not gain from it: the merge's 0.97% raw gain is gone after
  MGWR #src[`fy4b_tries.md`, try 15].

Table 2 lists every fix tried. Each was judged on the same gauges and hours, so the later ones were
held to a stricter level, and the line of fixes for the fusion was closed after the last.

#figure(
  table(
    columns: (auto, 1.05fr, 1fr),
    align: (right, left, left),
    stroke: none,
    inset: (x: 5pt, y: 3.5pt),
    rule,
    table.header([*\#*], [*Fix tried*], [*Outcome*]),
    thin,
    ..tr.map(row => (row.n, row.fix, row.outcome)).flatten(),
    rule,
  ),
  caption: [Table 2. Every fix tried for FY4B, #stat("fy4b_tries", "n") in all, read from
  `paper/fy4b_tries.csv` (the paper's Draft Table S1), so the section numbers in it are the paper's.
  Dates, tests and rules are in `Projects_Pending_Handover/fy4b_tries.md`.],
)

One more run came after the register closed. On 2026-11-15 FY4B's two-stage series (the chance of rain
times the mean wet amount of each hour's class) was scored in the full benchmark. Against corrected
FY4B its MGWR improved by 4.8%, but its agreement blend changed by +0.05% (−0.55 to +0.59), with CSI
at 8 mm h⁻¹ lower by 0.010. All four checks of the rule fixed before the run failed
#src[`output/fy4b_two_stage_benchmark_comparison/rule.csv` and `paired.csv`].

= Where FY4B is not the worst, and what was not tested

"Worst" holds for the averages in Table 1. It does not hold everywhere.

- *Single events.* On 2022-10-03 FY4B's pattern is the best of the three
  (r #val(e1, "r", product: "FY4B", event_day: "2022-10-03") against
  #val(e1, "r", product: G, event_day: "2022-10-03") and
  #val(e1, "r", product: "GSMaP", event_day: "2022-10-03")), and on 2023-06-17 it beats GPM
  (#val(e1, "r", product: "FY4B", event_day: "2023-06-17") against
  #val(e1, "r", product: G, event_day: "2023-06-17")).
- *Storm totals.* GSMaP's error on the event totals is larger: mean RMSE
  #val(r1, "RMSE", product: "GSMaP") mm against #val(r1, "RMSE", product: "FY4B") mm for FY4B
  (GPM #val(r1, "RMSE", product: G) mm), because GSMaP is far too wet
  (#val(r1, "MB", product: "GSMaP") mm).
- *Against GSMaP after correction.* Raw corrected FY4B against raw GSMaP differs by #dq("GSMaP", "raw")
  in RMSE, and after MGWR alone by #dq("GSMaP", "mgwr"), though its r stays lower
  (#val(fq, "r", scheme: "balanced_spatial", baseline: "GSMaP", method: "mgwr") against
  #val(fq, "r0", scheme: "balanced_spatial", baseline: "GSMaP", method: "mgwr")). After the blend the
  difference is #dq("GSMaP", B). These are same-weeks figures (§6). GSMaP's own fused gain over the
  gauges alone is not significant either
  (#val(mu, "vs_adw", product: "GSMaP", method: B)%,
  p = #val(cu, "p_holm", product: "GSMaP", method: B, stratum: "overall")).
- *Heavy hours.* On hours of 8 mm h⁻¹ or more, MGWR on FY4B improves on the gauges alone by
  #pct(cu, "change", product: "FY4B", method: "mgwr", stratum: "8_plus") and the agreement blend by
  #pct(cu, "change", product: "FY4B", method: B, stratum: "8_plus"), close to GPM's
  #pct(cu, "change", product: G, method: "mgwr", stratum: "8_plus") and
  #pct(cu, "change", product: G, method: B, stratum: "8_plus"). FY4B does carry some heavy-rain signal.

Not tested, so not ruled out:

- *More channels.* The channel pilot read six scans and only infrared and water-vapour channels. The
  visible and near-infrared channels are untried, and a few dozen more scans would give a firmer answer.
- *The 2023 software update and the move to 105°E.* Neither is shown to help or to hurt: each regime
  has one year, and 2023 was a wet one.
- *Other regions.* Every result here is from one gauge network of 237 stations, in one study area.

= Sources

- *The paper's numbers:* `paper/numbers/*.json`, exported from the run outputs by
  `scripts/export_paper_numbers.py` and checked by `scripts/check_paper_provenance.py`. Their text is
  in `paper/paper_plan.typ`, §3.1–§3.4.
- *The register of tries:* `Projects_Pending_Handover/fy4b_tries.md` and `paper/fy4b_tries.csv`.
- *The handover:* `Projects_Pending_Handover/fy4b_quality_handover.md`, for the error budget, the
  cloud-type pilot, the channel pilot and the raw-file facts.
- *Outputs read directly:* `output/fy4b_error_budget/`, `output/fy4b_process_trace/` and
  `output/fy4b_two_stage_benchmark_comparison/`.

"""Draw the FY4B process page from the tables of scripts/trace_fy4b_process.jl.

    $MIXEDGWR_PYTHON scripts/plot_fy4b_process.py     # -> output/fy4b_process_trace/fy4b_process.html

One self-contained HTML file (inline CSS, inline SVG, no library): FY4B from the raw 15-minute files
to the numbers the paper cites, on the one path the paper uses. Every count comes from a CSV beside
the page or from paper/numbers/*.json; the page computes shares and nothing else. It stops if the
trace's matched sample differs from the numbers the paper exports.
"""

from __future__ import annotations

import csv
import html
import json
import os
import sys
from pathlib import Path

from plot_provenance import write_plot_provenance

ROOT = Path(__file__).resolve().parents[1]
TRACE = ROOT / "output" / "fy4b_process_trace"
MANIFEST = ROOT / "output" / "fy4b_raw_manifest"
NUMBERS = ROOT / "paper" / "numbers"
PAGE = TRACE / "fy4b_process.html"
# The companion page (scripts/plot_fy4b_methods.py): beside this one on disk; a published copy
# passes its own address.
METHODS_LINK = os.environ.get("FY4B_METHODS_URL", "fy4b_methods.html")

# Loss causes keep one colour each wherever they appear; "complete" is the kept series.
CAUSE_COLOR = {"complete": "var(--s1)", "before_first_file": "var(--s7)", "relocation_gap": "var(--s2)",
               "quarter_1415_utc": "var(--s4)", "no_files": "var(--s5)", "quarters_missing": "var(--s3)"}
CAUSE_SHORT = {"complete": "Kept", "before_first_file": "Before 1 June 2022", "relocation_gap": "133°E → 105°E gap",
               "quarter_1415_utc": "14:15 UTC file", "no_files": "No file", "quarters_missing": "1–3 quarters missing"}
REGIME_NAME = {"133.0E_V1.0": "133°E, software V1.0", "133.0E_V1.0.0": "133°E, software V1.0.0",
               "105.0E_V1.0.0": "105°E, software V1.0.0"}
PRODUCT_NAME = {"FY4B": "FY4B", "FY4B_QM": "Corrected FY4B", "GPM_uncal": "GPM", "GSMaP": "GSMaP",
                "MERGED_MEAN": "Merged (equal weights)", "MERGED_OLS": "Merged (OLS)"}
METHOD_NAME = {"raw": "raw", "mgwr": "MGWR", "blend_agrenv_mgwr": "MGWR + agreement blend", "adw": "ADW (gauges only)"}


def rows(name: str, folder: Path = TRACE) -> list[dict[str, str]]:
    with open(folder / name, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def numbers(name: str) -> dict:
    return json.loads((NUMBERS / f"{name}.json").read_text(encoding="utf-8"))


def esc(text) -> str:
    return html.escape(str(text), quote=True)


def n(value) -> str:
    """An integer with thousands separators."""
    return f"{int(float(value)):,}"


def f(value, digits: int = 2) -> str:
    return f"{float(value):.{digits}f}".replace("-", "−")


def pct(part, whole, digits: int = 1) -> str:
    return f"{100 * float(part) / float(whole):.{digits}f}%"


def minus(text: str) -> str:
    return str(text).replace("-", "−")


def table(head: list[str], body: list[list[str]], num: tuple[int, ...] = (), cls: str = "") -> str:
    """An HTML table in its own scroll box; `num` lists the right-aligned columns. Cells are HTML."""
    th = "".join(f'<th class="{"num" if i in num else ""}">{cell}</th>' for i, cell in enumerate(head))
    trs = "".join("<tr>" + "".join(f'<td class="{"num" if i in num else ""}">{cell}</td>'
                                   for i, cell in enumerate(row)) + "</tr>" for row in body)
    return f'<div class="scroll"><table class="{cls}"><thead><tr>{th}</tr></thead><tbody>{trs}</tbody></table></div>'


def source(*paths: str) -> str:
    return '<p class="src">Source: ' + ", ".join(f"<code>{esc(p)}</code>" for p in paths) + "</p>"


def legend(items: list[tuple[str, str]]) -> str:
    return '<div class="legend">' + "".join(
        f'<span><i style="background:{color}"></i>{esc(label)}</span>' for label, color in items) + "</div>"


def section(anchor: str, number: str, title: str, lead: str, body: str) -> str:
    return (f'<section id="{anchor}"><header class="sec"><span class="no">{number}</span><h2>{esc(title)}</h2></header>'
            f'<p class="lead">{lead}</p>{body}</section>')


# --- Charts -----------------------------------------------------------------------------------------

def stacked_bar(parts: list[tuple[str, str, int]], total: int) -> str:
    """One horizontal bar split into parts (label, colour, value), with a 2px surface gap."""
    width, height, x = 1000, 44, 0.0
    marks = []
    for label, color, value in parts:
        w = width * value / total
        marks.append(f'<rect x="{x + 1:.2f}" y="0" width="{max(w - 2, 0.6):.2f}" height="{height}" rx="3" fill="{color}" '
                     f'data-tip="{esc(label)}: {value:,} hours ({pct(value, total)})"></rect>')
        x += w
    return (f'<svg class="chart" viewBox="0 0 {width} {height}" preserveAspectRatio="none" role="img" '
            f'aria-label="Hours of 2022 to 2024 by what happened to them" style="height:44px">{"".join(marks)}</svg>')


def hbars(items: list[tuple[str, float, str]], color: str = "var(--s1)", label_width: int = 300,
          unit: str = "") -> str:
    """Horizontal bars with the label left and the value at the bar's end. items: (label, value, tip)."""
    width, row, top = 1000, 34, 4
    peak = max(value for _, value, _ in items)
    room = width - label_width - 150
    out = []
    for i, (label, value, tip) in enumerate(items):
        y = top + i * row
        w = max(room * value / peak, 1)
        out.append(f'<text x="{label_width - 10}" y="{y + 17}" text-anchor="end" class="lab">{esc(label)}</text>'
                   f'<rect x="{label_width}" y="{y + 4}" width="{w:.1f}" height="20" rx="3" fill="{color}" '
                   f'data-tip="{esc(tip)}"></rect>'
                   f'<text x="{label_width + w + 8:.1f}" y="{y + 19}" class="val">{value:,.0f}{unit}</text>')
    height = top * 2 + len(items) * row
    return (f'<div class="scroll"><svg class="chart wide" viewBox="0 0 {width} {height}" role="img" '
            f'style="min-width:640px">{"".join(out)}</svg></div>')


def columns(items: list[tuple[str, float, str]], color, height: int = 190, peak: float | None = None,
            tick_every: int = 1, y_label: str = "", fmt="{:,.0f}", ticks: int = 4) -> str:
    """Vertical bars from a zero baseline. items: (x label, value, tip); `color` a string or a function of the value."""
    width, left, bottom, top = 1000, 56, 26, 12
    peak = peak or max(value for _, value, _ in items) or 1
    plot_h, plot_w = height - bottom - top, width - left - 8
    step = plot_w / len(items)
    out = []
    for k in range(ticks + 1):
        value = peak * k / ticks
        y = top + plot_h * (1 - k / ticks)
        out.append(f'<line x1="{left}" x2="{width - 8}" y1="{y:.1f}" y2="{y:.1f}" class="{"axis" if k == 0 else "grid"}"></line>'
                   f'<text x="{left - 8}" y="{y + 4:.1f}" text-anchor="end" class="tick">{fmt.format(value)}</text>')
    for i, (label, value, tip) in enumerate(items):
        h = plot_h * value / peak
        x = left + i * step
        fill = color(value) if callable(color) else color
        if value > 0:
            out.append(f'<rect x="{x + 1:.1f}" y="{top + plot_h - h:.1f}" width="{max(step - 2, 1):.1f}" height="{max(h, 1):.1f}" '
                       f'rx="2" fill="{fill}"></rect>')
        out.append(f'<rect x="{x:.1f}" y="{top}" width="{step:.1f}" height="{plot_h}" fill="transparent" data-tip="{esc(tip)}"></rect>')
        if i % tick_every == 0:
            out.append(f'<text x="{x + step / 2:.1f}" y="{height - 8}" text-anchor="middle" class="tick">{esc(label)}</text>')
    caption = f'<p class="axis-label">{esc(y_label)}</p>' if y_label else ""
    return (f'{caption}<div class="scroll"><svg class="chart wide" viewBox="0 0 {width} {height}" role="img" '
            f'style="min-width:640px">{"".join(out)}</svg></div>')


def diverging_columns(items: list[tuple[str, float, str]], height: int = 230, limit: float = 200) -> str:
    """Bars above and below zero: FY4B wetter than the gauges up (blue), drier down (orange)."""
    width, left, bottom, top = 1000, 56, 26, 12
    plot_h, plot_w = height - bottom - top, width - left - 8
    low = -100.0
    scale = plot_h / (limit - low)
    zero = top + limit * scale
    step = plot_w / len(items)
    out = []
    for value in (-100, -50, 0, 50, 100, 150, 200):
        y = zero - value * scale
        label = "0%" if value == 0 else minus(f"{value:+d}%")
        out.append(f'<line x1="{left}" x2="{width - 8}" y1="{y:.1f}" y2="{y:.1f}" class="{"axis" if value == 0 else "grid"}"></line>'
                   f'<text x="{left - 8}" y="{y + 4:.1f}" text-anchor="end" class="tick">{label}</text>')
    for i, (label, value, tip) in enumerate(items):
        x = left + i * step
        h = abs(min(value, limit)) * scale
        y = zero - h if value >= 0 else zero
        out.append(f'<rect x="{x + 1.5:.1f}" y="{y:.1f}" width="{max(step - 3, 1):.1f}" height="{max(h, 1):.1f}" rx="2" '
                   f'fill="{"var(--s1)" if value >= 0 else "var(--s2)"}"></rect>'
                   f'<rect x="{x:.1f}" y="{top}" width="{step:.1f}" height="{plot_h}" fill="transparent" data-tip="{esc(tip)}"></rect>')
        if label.endswith(("-01", "-06")):
            out.append(f'<text x="{x + step / 2:.1f}" y="{height - 8}" text-anchor="middle" class="tick">{esc(label)}</text>')
    return (f'<div class="scroll"><svg class="chart wide" viewBox="0 0 {width} {height}" role="img" '
            f'style="min-width:640px">{"".join(out)}</svg></div>')


def heat(value: float) -> str:
    """Sequential blue step for a share in 0..1; 0 is the empty cell."""
    if value <= 0:
        return "var(--empty)"
    return f"var(--q{min(int(value * 5), 4) + 1})"


def month_grid(cells: dict[str, tuple[float, str, str]], years: list[int], title_of_row=str) -> str:
    """Rows of twelve month cells. cells: 'yyyy-mm' -> (share 0..1, text in the cell, tip)."""
    head = "".join(f"<th>{m}</th>" for m in "JFMAMJJASOND")
    body = []
    for year in years:
        tds = []
        for month in range(1, 13):
            share, text, tip = cells.get(f"{year}-{month:02d}", (0.0, "", "no hours"))
            ink = "var(--on-high)" if share >= 0.6 else "var(--fg)"
            tds.append(f'<td style="background:{heat(share)};color:{ink}" data-tip="{esc(tip)}">{esc(text)}</td>')
        body.append(f"<tr><th>{esc(title_of_row(year))}</th>{''.join(tds)}</tr>")
    return f'<div class="scroll"><table class="grid"><thead><tr><th></th>{head}</tr></thead><tbody>{"".join(body)}</tbody></table></div>'


def regime_timeline(regimes: list[dict[str, str]]) -> str:
    from datetime import datetime
    start, end = datetime(2022, 1, 1), datetime(2025, 1, 1, 8)
    span = (end - start).total_seconds()
    width = 1000
    colors = ["var(--s1)", "var(--s3)", "var(--s7)"]
    out = []
    for year in (2022, 2023, 2024, 2025):
        x = width * (datetime(year, 1, 1) - start).total_seconds() / span
        out.append(f'<line x1="{x:.1f}" x2="{x:.1f}" y1="0" y2="62" class="grid"></line>'
                   f'<text x="{min(x + 4, width - 34):.1f}" y="78" class="tick">{year}</text>')
    for regime, color in zip(regimes, colors):
        a = datetime.fromisoformat(regime["first_hour"])
        b = datetime.fromisoformat(regime["last_hour"])
        x1, x2 = (width * (t - start).total_seconds() / span for t in (a, b))
        out.append(f'<rect x="{x1 + 1:.1f}" y="14" width="{x2 - x1 - 2:.1f}" height="34" rx="4" fill="{color}" '
                   f'data-tip="{esc(REGIME_NAME[regime["regime"]])}: {a:%d %b %Y} to {b:%d %b %Y}, {n(regime["hours"])} complete hours"></rect>')
    return (f'<div class="scroll"><svg class="chart wide" viewBox="0 0 {width} 84" role="img" style="min-width:640px">'
            f'{"".join(out)}</svg></div>'
            + legend([(REGIME_NAME[r["regime"]], c) for r, c in zip(regimes, colors)]))


# --- Sections ---------------------------------------------------------------------------------------

def build() -> str:
    funnel = {r["cause"]: r for r in rows("hour_funnel.csv")}
    sample = {r["stage"]: r for r in rows("sample_funnel.csv")}
    raw = {r["item"]: r["value"] for r in rows("raw_files.csv")}
    checks = rows("checks.csv")
    regimes = rows("regimes.csv")
    steps = rows("steps.csv")
    scope = {r["key"]: r["value"] for r in rows("benchmark_scope.csv")}
    trace_prov = {r["key"]: r["value"] for r in rows("provenance.csv")}
    has_manifest = (MANIFEST / "folders.csv").is_file()

    # The page and the paper must describe one sample.
    data_hours = numbers("data_hours")["rows"][0]
    rq3 = numbers("rq3_sample")["rows"][0]
    matched, scored = sample["matched_sample"], sample["scored"]
    for label, ours, theirs in (("matched hours", matched["hours"], data_hours["hours"]),
                                ("matched days", matched["days"], data_hours["days"]),
                                ("matched gauge-hours", matched["gauge_hours"], data_hours["cells"]),
                                ("scored gauge-hours", scored["gauge_hours"], rq3["n"]),
                                ("scored hours", scored["hours"], rq3["hours"])):
        if n(ours) != theirs:
            sys.exit(f"{label}: the trace has {n(ours)}, paper/numbers has {theirs}")

    total_hours = sum(int(r["hours"]) for r in funnel.values())
    complete = int(funnel["complete"]["hours"])
    passed = sum(r["pass"] == "true" for r in checks)
    all_pass = passed == len(checks)
    out: list[str] = []

    # --- Header -------------------------------------------------------------------------------
    status = (f'<span class="pill ok">✓ {passed} of {len(checks)} checks pass</span>' if all_pass else
              f'<span class="pill bad">✕ {len(checks) - passed} of {len(checks)} checks fail</span>')
    dirty = trace_prov.get("git_dirty") == "true" or trace_prov.get("git_dirty_at_end") == "true"
    out.append(f'''<header class="top">
<p class="eyebrow">Satellite precipitation correction · FY-4B AGRI QPE · Shiyan study area</p>
<h1>FY4B process trace</h1>
<p class="lead">How FY4B gets from the raw 15-minute files on the drive to the numbers in <code>paper/paper_plan.typ</code>,
stage by stage, with the count left after each stage. Only the path the paper uses is followed: the rebuilt 3 × 3 series
as the third product of RQ1–RQ3, and corrected FY4B from the <code>benchmark</code> step.</p>
<div class="meta">{status}
<span>Tables: <code>scripts/trace_fy4b_process.jl</code> at <code>{esc(trace_prov.get("git_commit", "")[:7])}</code>{" (uncommitted changes)" if dirty else ", clean tree"}</span>
<span>Rebuild: <code>{esc(steps[1]["git_commit"][:7])}</code></span>
<span>Benchmark: <code>{esc(scope["git_commit"][:7])}</code></span>
<span>How each step works, with formulas: <a href="{esc(METHODS_LINK)}">FY4B methods</a></span></div>
<div class="tiles">
<div><b>{n(raw["15-minute files used, after duplicates are dropped"])}</b><span>15-minute files used</span></div>
<div><b>{n(complete)}</b><span>hours with all four files, of {n(total_hours)}</span></div>
<div><b>{n(matched["hours"])}</b><span>hours in the matched sample, on {n(matched["days"])} days</span></div>
<div><b>{n(scored["gauge_hours"])}</b><span>gauge-hours scored in Table 3</span></div>
</div></header>''')

    # --- 0. The path at a glance ----------------------------------------------------------------
    drive_files = raw.get("15-minute files on the drive, every copy (raw drive)")
    path_rows = [
        ("Raw files", "FY-4B AGRI L2 QPE, full disk, 4 km, one file per 15-minute scan. Repeated downloads are dropped: of several copies, the first that opens.",
         (f"{n(drive_files)} files → " if drive_files else "") + f'{n(raw["15-minute files used, after duplicates are dropped"])} scans', "FY4BPreprocessing.find_nc_files"),
        ("Strict hour", "An hour is kept only when all four of its 15-minute files exist and open.",
         f"{n(complete)} of {n(total_hours)} hours", "FY4BPreprocessing.build_hourly_qc"),
        ("Value at a gauge", "Mean of the 3 × 3 pixels around the gauge in each quarter (5 of 9 must be valid), then the mean of the four quarters. No parallax shift.",
         f'{n(sample["fy4b_value"]["hours"])} hours with a value, {int(sample["record"]["gauge_hours"]) // int(sample["record"]["hours"])} gauges', "FY4BPreprocessing.sample_variants"),
        ("Matched sample", "Gauge-hours at which the gauge, FY4B, GPM (IMERG Uncal) and GSMaP all have a value. RQ1–RQ3 are scored on it.",
         f'{n(matched["gauge_hours"])} gauge-hours on {n(matched["hours"])} hours', "SatelliteTemporalEvaluation.sample_mask"),
        ("Pixel shift, per fold", "For each of three regimes, the shift (up to ±6 pixels) under which FY4B correlates best with the fold's training gauges.",
         "3 regimes × 5 folds", "SatelliteCalibration.shift_fold"),
        ("Wet threshold and quantile map, per fold", "Per season and regime: a threshold that gives FY4B the training gauges' wet-hour frequency, and a map of its wet amounts onto theirs.",
         "12 season × regime cells per fold", "SatelliteCalibration.fit_quantile_map"),
        ("Strength, month by month", "The map is applied only as strongly as the training gauges support: not at all, the threshold alone, or a quarter to all of the map.",
         f'{len({r["cell"] + r["regime"] for r in rows("strengths.csv")})} month cells per fold', "SatelliteCalibration.choose_strength"),
        ("Benchmark", "Plain FY4B is an anchor, and enters the merges and the agreement bands. Corrected FY4B (FY4B_QM) is appended last as a further anchor and enters neither.",
         "6 anchors, 5 folds, 2 CV schemes", "scripts/run_interpolation_benchmark.jl"),
        ("Rescoring on one sample", "Every anchor and method is rescored on the gauge-hours where every anchor's prediction by every method exists.",
         f'{n(scored["gauge_hours"])} gauge-hours on {n(scored["hours"])} hours', "BenchmarkDiagnostics.matched_evaluation_mask"),
        ("Paper", "Table 3, §3.4, Tables 4 and 5 read the exported numbers; nothing is typed by hand.",
         "paper/numbers/*.json", "scripts/export_paper_numbers.py"),
    ]
    flow = "".join(f'<li><div class="what"><h3>{esc(name)}</h3><p>{esc(text)}</p><code>{esc(code)}</code></div>'
                   f'<div class="count">{esc(count)}</div></li>' for name, text, count, code in path_rows)
    out.append(section("path", "0", "The path at a glance",
                       "Ten stages, in the order the data passes through them. Stages 1–4 use no gauge; 5–7 are fitted inside each cross-validation fold on its training gauges only.",
                       f'<ol class="flow">{flow}</ol>'))

    # --- 1. Raw files -------------------------------------------------------------------------
    raw_table = table(["What", "Count"], [[esc(k), esc(n(v) if v.isdigit() else v)] for k, v in raw.items()], num=(1,))
    folders = ""
    if has_manifest:
        folder_rows = rows("folders.csv", MANIFEST)
        folders = ("<h3>Where the copies lie on the drive</h3>" + table(
            ["Folder", "15-minute files", "Of which used", "Size, GB", "First scan", "Last scan"],
            [[f'<code>{esc(r["folder"])}</code>', n(r["files"]), n(r["kept_files"]), f(int(r["bytes"]) / 1e9, 1),
              esc(r["first_scan"][:10]), esc(r["last_scan"][:10])] for r in folder_rows], num=(1, 2, 3))
            + source("output/fy4b_raw_manifest/folders.csv"))
    out.append(section("raw", "1", "Raw files", "What is on the drive, what is read, and what the record can and cannot prove.", f'''
<dl class="facts">
<dt>Product</dt><dd>FY-4B AGRI Level 2 QPE (quantitative precipitation estimate), full disk, 4 km at nadir, 2748 × 2748 pixels, one file per 15-minute scan, from NSMC's FengYun data service. No tracked script downloads them.</dd>
<dt>File name</dt><dd><code>FY4B-_AGRI--_N_DISK_1330E_L2-_QPE-_MULT_NOM_&lt;start&gt;_&lt;end&gt;_4000M_V0001.NC</code>. The satellite position (<code>1330E</code> or <code>1050E</code>) and the scan's start and end (UTC) are read from the name; a file counts as a 15-minute scan when it starts at :00, :15, :30 or :45 and ends 14 min 59 s later.</dd>
<dt>Never read</dt><dd>The 1-, 3- and 24-hour QPE files in the same tree. They are corrected with China's rain gauges, so they would leak gauge information into a product compared against gauges. macOS <code>._</code> companion files are skipped too.</dd>
<dt>Repeated downloads</dt><dd>The tree holds several download batches. Of several copies of one scan, the first that opens is kept, in sorted path order, so a damaged download cannot shadow a good copy.</dd>
<dt>Location</dt><dd>External drive, linked as <code>data/FY4B</code>. Only read, never written.</dd>
</dl>
{raw_table}{source("output/fy4b_process_trace/raw_files.csv")}
{folders}
<div class="note"><h3>What is fingerprinted, and what is not</h3>
<p>Every other input of the paper's run has a SHA-256 saved with the run. The raw 15-minute files do not: the rebuild records only the name of each file it used.
{"<b>Names and sizes are now checked</b> by <code>scripts/list_fy4b_raw_files.jl</code>: the drive holds exactly the scans the rebuild used, name for name, and the count of repeated copies equals the rebuild's. Seven scans have copies of different size; in each, the odd copy is a damaged download that does not open, and the rule keeps the copy that opens. The bytes of the files are not hashed, so a file changed without a change of name or size would not be seen." if has_manifest else "The raw drive has not been listed yet: run <code>scripts/list_fy4b_raw_files.jl</code> with the drive mounted."}</p></div>'''))

    # --- 2. Strict hour -------------------------------------------------------------------------
    order = ["complete", "before_first_file", "relocation_gap", "quarter_1415_utc", "quarters_missing", "no_files"]
    bar = stacked_bar([(CAUSE_SHORT[k], CAUSE_COLOR[k], int(funnel[k]["hours"])) for k in order], total_hours)
    cause_table = table(["", "What happened to the hour", "Hours", "Share", "Files in them", "First (UTC)", "Last (UTC)"],
                        [[f'<i class="sw" style="background:{CAUSE_COLOR[k]}"></i>', esc(funnel[k]["what"]), n(funnel[k]["hours"]),
                          pct(funnel[k]["hours"], total_hours), n(funnel[k]["files_in_them"]),
                          esc(funnel[k]["first_hour_utc"].replace("T", " ")[:16]), esc(funnel[k]["last_hour_utc"].replace("T", " ")[:16])]
                         for k in order], num=(2, 3, 4))
    months = {r["month"]: r for r in rows("hours_by_month.csv")}
    grid = month_grid({m: (int(r["hours_complete"]) / int(r["hours_possible"]), f'{100 * int(r["hours_complete"]) / int(r["hours_possible"]):.0f}',
                           f'{m}: {n(r["hours_complete"])} of {n(r["hours_possible"])} hours complete; {n(r["files_found"])} of {n(r["files_expected"])} files')
                       for m, r in months.items()}, [2022, 2023, 2024])
    by_hour = rows("hours_by_utc_hour.csv")
    hour_chart = columns([(f'{int(r["utc_hour"]):02d}', 100 * int(r["hours_complete"]) / int(r["hours_in_record"]),
                           f'{int(r["utc_hour"]):02d}–{(int(r["utc_hour"]) + 1) % 24:02d} UTC (hour ending {int(r["beijing_hour_ending"]):02d}:00 Beijing time): '
                           f'{n(r["hours_complete"])} of {n(r["hours_in_record"])} days complete')
                          for r in by_hour], "var(--s1)", peak=100, fmt="{:.0f}%",
                         y_label="Share of days on which the hour is complete, by UTC hour (1 June 2022 – 31 December 2024, the relocation gap left out)")
    quarters = rows("missing_quarters.csv")
    quarter_table = table(["Quarters missing (minute of the scan start)", "Hours"],
                          [[esc(r["missing_minutes"].replace(";", ", ") or "none listed"), n(r["hours"])] for r in quarters], num=(1,))
    out.append(section("hour", "2", "From 15 minutes to an hour",
                       f"The strict rule: all four files of the hour, or no value. Of the {n(total_hours)} hours of 2022–2024, {n(complete)} pass ({pct(complete, total_hours)}).", f'''
{bar}{legend([(CAUSE_SHORT[k], CAUSE_COLOR[k]) for k in order])}
{cause_table}{source("output/fy4b_process_trace/hour_funnel.csv")}
<p>The causes do not overlap: an hour is given the first that applies, in the order of the table below the kept hours.
The hour 14–15 UTC is never complete because its 14:15 file is never produced, so no 08–08 day has all 24 hours.</p>
<h3>Complete hours by month, % of the month's hours</h3>
{grid}{source("output/fy4b_process_trace/hours_by_month.csv")}
<h3>By hour of the day</h3>
{hour_chart}{source("output/fy4b_process_trace/hours_by_utc_hour.csv")}
<h3>Which quarters are missing in the {n(funnel["quarters_missing"]["hours"])} partly missing hours</h3>
{quarter_table}{source("output/fy4b_process_trace/missing_quarters.csv")}'''))

    # --- 3. Sampling at a gauge -----------------------------------------------------------------
    window = rows("window_choice.csv")
    valueless = rows("valueless_hours.csv")
    v_all = next(r for r in valueless if r["group"] == "all")
    v_hours = [r for r in valueless if r["group"] == "utc_hour"]
    v_months = {r["level"]: r for r in valueless if r["group"] == "month"}
    v_chart = columns([(r["level"], int(r["no_value_at_any_gauge"]),
                        f'{r["level"]}–{(int(r["level"]) + 1) % 24:02d} UTC (hour ending {(int(r["level"]) + 9) % 24:02d}:00 Beijing time): '
                        f'{n(r["no_value_at_any_gauge"])} of {n(r["complete_hours"])} complete hours have no value')
                       for r in v_hours], "var(--s2)", height=170,
                      y_label="Complete hours with no FY4B value at any gauge, by UTC hour")
    v_peak = max(int(r["no_value_at_any_gauge"]) for r in v_months.values()) or 1
    v_grid = month_grid({m: (int(r["no_value_at_any_gauge"]) / v_peak if int(r["no_value_at_any_gauge"]) else 0.0,
                             r["no_value_at_any_gauge"] if int(r["complete_hours"]) else "",
                             f'{m}: {r["no_value_at_any_gauge"]} of {n(r["complete_hours"])} complete hours have no value')
                         for m, r in v_months.items()}, [2022, 2023, 2024])
    out.append(section("sampling", "3", "The value at a gauge",
                       "How one number per gauge and hour is read from four full-disk images. Nothing in this stage looks at a gauge reading.", f'''
<dl class="facts">
<dt>Finding the pixel</dt><dd>The gauge's latitude and longitude are turned into scan angles of the geostationary view and then into a column and row of the 2748 × 2748 grid, using the satellite longitude stored in the file (133°E, later 105°E). No parallax correction: a fixed 10 km cloud-top shift was built and scored worse on June 2022, so it is not used.</dd>
<dt>Accepting a pixel</dt><dd>A pixel is valid when its value is not the fill value, lies in 0–30 mm h⁻¹ (which drops the 65533–65535 sentinels; real values stop at 20), and its quality flag <code>DQF</code> is below 2.</dd>
<dt>3 × 3 mean</dt><dd>Each quarter's value is the mean of the valid pixels in the 3 × 3 box around the gauge's pixel (about 12 km), and needs at least 5 of the 9. The rule was fixed in advance, not tuned.</dd>
<dt>The hour</dt><dd>The mean of the four quarters' rates, which is the hour's accumulation in mm. One quarter without a value voids the gauge-hour. The label is the hour's end in Beijing time (UTC + 8).</dd>
<dt>Ceiling</dt><dd>The 15-minute product stops at 20 mm h⁻¹. The share of pixels at the ceiling is saved beside the series (<code>fy4b_w3_h0_capfrac.csv</code>).</dd>
<dt>Built in the same pass</dt><dd>Four variants: centre pixel (<code>w1_h0</code>), 3 × 3 (<code>w3_h0</code>, used), 5 × 5 (<code>w5_h0</code>) and 3 × 3 with the 10 km parallax shift (<code>w3_h10</code>). Also <code>fy4b_offset_grid.nc</code>: the 3 × 3 mean with the pixel moved by every shift up to ±6 in x and y (169 shifts), from which stage 5 picks.</dd>
</dl>
<h3>The 3 × 3 mean against the centre pixel, raw, on the matched sample</h3>
{table(["Variant", "Gauge-hours", "r", "RMSE, mm h⁻¹", "Wet share", "Mean, mm h⁻¹"],
       [[esc(r["variant"]), n(r["cells"]), f(r["r"], 3), f(r["rmse"], 3), pct(r["wet_share"], 1), f(r["mean"], 3)] for r in window], num=(1, 2, 3, 4, 5))}
{source("output/fy4b_process_trace/window_choice.csv")}
<h3>Complete hours that still have no value</h3>
<p>{n(v_all["no_value_at_any_gauge"])} of the {n(v_all["complete_hours"])} complete hours have no FY4B value at any of the 237 gauges: all four files exist, but the values over the study area are rejected.
The centre-pixel series is empty on the same {n(v_all["centre_pixel_no_value_at_any_gauge"])} hours, so the 5-of-9 rule is not the cause. {n(v_all["value_at_some_gauges_only"])} further hours have a value at some gauges only.
This is the step from {n(v_all["complete_hours"])} to the paper's {n(matched["hours"])} hours.</p>
{v_chart}
<h3>The same hours by month</h3>
{v_grid}{source("output/fy4b_process_trace/valueless_hours.csv")}
<p class="aside">They fall almost entirely at 15–19 UTC (the hours ending 00:00–03:00 Beijing time) in March–April and August–October. That is local midnight in the weeks around the equinoxes, when a geostationary imager looks close to the sun. The timing fits that explanation; the trace only shows that the values are rejected, and the quality flags of those files have not been read to confirm it.</p>'''))

    # --- 4. Regimes -----------------------------------------------------------------------------
    out.append(section("regimes", "4", "Three regimes",
                       "FY4B's record falls into three periods, read from each file's satellite longitude and software version. The shift and the map of stages 5–6 are fitted separately in each.",
                       regime_timeline(regimes) + table(
                           ["Regime", "First hour (Beijing time, hour ending)", "Last hour", "Complete hours", "Hours in the matched sample"],
                           [[esc(REGIME_NAME[r["regime"]]), esc(r["first_hour"].replace("T", " ")[:16]), esc(r["last_hour"].replace("T", " ")[:16]),
                             n(r["hours"]), n(r["matched_hours"])] for r in regimes], num=(3, 4))
                       + source("output/fy4b_process_trace/regimes.csv", "output/fy4b_rebuild_2022_2024/fy4b_hour_meta.csv")
                       + f'<p>The software changes on 16 May 2023. The satellite is at 133°E until the hour ending 08:00 on 1 February 2024 and at 105°E from the hour ending 09:00 on 5 March 2024; '
                         f'the {n(funnel["relocation_gap"]["hours"])} hours between hold no file at all.</p>'))

    # --- 5. Matched sample -----------------------------------------------------------------------
    stage_names = {"record": "Record: every gauge, every hour", "strict_hour": "Hour has all four FY4B files",
                   "fy4b_value": "FY4B has a value at the gauge", "gauge_value": "… and the gauge has a value",
                   "gpm_value": "… and GPM has a value", "matched_sample": "… and GSMaP has a value",
                   "scored": "… and every prediction exists"}
    bars = hbars([(stage_names[k], int(r["gauge_hours"]), f'{r["what"]}: {n(r["gauge_hours"])} gauge-hours on {n(r["hours"])} hours')
                  for k, r in sample.items()])
    previous = None
    sample_rows = []
    for k, r in sample.items():
        removed = "" if previous is None else n(int(previous["gauge_hours"]) - int(r["gauge_hours"]))
        sample_rows.append([esc(r["what"]), n(r["hours"]), n(r["days"]) if int(r["days"]) else "—", n(r["gauge_hours"]), removed])
        previous = r
    folds = [r for r in rows("sample_by_fold.csv") if r["scheme"] == "balanced_spatial"]
    fold_table = table(["Anchor", "Fold", "Gauge-hours scored on its own mask", "On the shared sample", "Share kept", "Hours"],
                       [[esc(PRODUCT_NAME[r["product"]]), esc(r["fold"]), n(r["n_own"]), n(r["n_matched"]),
                         pct(r["n_matched"], r["n_own"], 2), n(r["hours_matched"])] for r in folds], num=(2, 3, 4, 5))
    out.append(section("sample", "5", "The matched sample",
                       "Every comparison in the paper is made on the gauge-hours at which the gauge and all three products have a value.", f'''
{bars}
{table(["Condition added", "Hours", "08–08 days", "Gauge-hours", "Removed by this step"], sample_rows, num=(1, 2, 3, 4))}
{source("output/fy4b_process_trace/sample_funnel.csv")}
<p>GPM and GSMaP remove nothing: both have a value at every gauge and hour of the record, so the sample is set by FY4B and by the gauges alone.
The last step is the benchmark's: {n(int(matched["gauge_hours"]) - int(scored["gauge_hours"]))} gauge-hours are left out because some method has no prediction there for some anchor.</p>
<h3>Held-out gauge-hours per fold, balanced spatial CV</h3>
{fold_table}{source("output/gauge_satellite_fusion_fy4bw3h0_fy4bqm_fy4bqmadapt_gpmuncal_matched_block7d/matched_sample.csv")}
<p>Fold 3 loses the most (2.7%): its shared hours fall to 19,292, because the shared sample drops an hour for every anchor when one anchor has no prediction by some method there.</p>'''))

    # --- 6. The series against the gauges --------------------------------------------------------
    summary = rows("series_summary.csv")
    def summary_row(r):
        name = REGIME_NAME.get(r["level"], r["level"])
        return [esc({"all": "All"}.get(name, name)), n(r["hours"]), n(r["cells"]), pct(r["gauge_wet_share"], 1), pct(r["fy4b_wet_share"], 1),
                f(r["gauge_mean"], 3), f(r["fy4b_mean"], 3), minus(f'{float(r["relative_bias_pct"]):+.0f}%'), f(r["r"], 2), f(r["rmse"], 2),
                f(r["pod"], 2), f(r["far"], 2), n(r["gauge_cells_ge20"]), n(r["fy4b_cells_ge20"])]
    head = ["", "Hours", "Gauge-hours", "Gauge wet", "FY4B wet", "Gauge mean", "FY4B mean", "Bias", "r", "RMSE", "POD", "FAR", "Gauge ≥ 20", "FY4B ≥ 20"]
    blocks = "".join(f"<h3>{title}</h3>" + table(head, [summary_row(r) for r in summary if r["group"] == group], num=tuple(range(1, 14)))
                     for group, title in (("all", "Whole matched sample"), ("regime", "By regime"), ("season", "By season")))
    monthly = rows("series_by_month.csv")
    bias_chart = diverging_columns([(r["level"], float(r["relative_bias_pct"]),
                                     f'{r["level"]}: FY4B {float(r["relative_bias_pct"]):+.0f}% against the gauges '
                                     f'(mean {float(r["fy4b_mean"]):.3f} vs {float(r["gauge_mean"]):.3f} mm/h, r {float(r["r"]):.2f}, {n(r["hours"])} hours)')
                                    for r in monthly])
    out.append(section("series", "6", "The rebuilt series against the gauges",
                       "Plain FY4B, before any correction, on the matched sample. This is the FY4B of RQ1–RQ3 and the error the correction of stage 5–7 works on.", f'''
{blocks}{source("output/fy4b_process_trace/series_summary.csv")}
<p class="aside">Wet is ≥ 0.1 mm h⁻¹. Means in mm h⁻¹; bias is (FY4B total − gauge total) / gauge total. POD is the share of gauge-wet hours FY4B calls wet, FAR the share of FY4B-wet hours the gauge calls dry. The last two columns count gauge-hours at 20 mm h⁻¹ or more: FY4B cannot exceed its ceiling of 20, and reaches it only when the whole box stays there for the whole hour.</p>
<h3>Bias by month</h3>
<p class="axis-label">FY4B total against the gauge total, by month (bars above +200% are cut at the top)</p>
{bias_chart}{legend([("FY4B wetter than the gauges", "var(--s1)"), ("FY4B drier", "var(--s2)")])}
{source("output/fy4b_process_trace/series_by_month.csv")}'''))

    # --- 7. Correction inside each fold ----------------------------------------------------------
    shifts = [r for r in rows("shifts.csv") if r["scheme"] == "balanced_spatial"]
    regime_order = [r["regime"] for r in regimes]
    shift_rows = []
    for fold in sorted({r["fold"] for r in shifts}):
        cells = [esc(fold), n(next(r["n_train_station"] for r in shifts if r["fold"] == fold))]
        for regime in regime_order:
            r = next(x for x in shifts if x["fold"] == fold and x["regime"] == regime)
            cells += [minus(f'{int(r["dx"]):+d}, {int(r["dy"]):+d}'), f'{f(r["r_train_unshifted"], 3)} → {f(r["r_train"], 3)}']
        shift_rows.append(cells)
    shift_head = ["Fold", "Training gauges"] + [x for regime in regime_order for x in (f"{REGIME_NAME[regime]}: shift x, y", "training r, before → after")]
    maps = [r for r in rows("maps.csv") if r["scheme"] == "balanced_spatial"]
    map_rows = []
    for regime in regime_order:
        for season in ("DJF", "MAM", "JJA", "SON"):
            group = [r for r in maps if r["regime"] == regime and r["season"] == season]
            if not group:
                continue
            mean = lambda key: sum(float(r[key]) for r in group) / len(group)
            taus = [float(r["tau"]) for r in group]
            pools = sorted({r["pool"] for r in group})
            map_rows.append([esc(REGIME_NAME[regime]), season, esc({"cell": "own", "regime": "whole regime", "pooled": "all hours"}.get(pools[0], pools[0]) if len(pools) == 1 else "mixed"),
                             f"{min(taus):.2f}–{max(taus):.2f}", f'{mean("q50_sat"):.2f} → {mean("q50_obs"):.2f}', f'{mean("q99_sat"):.1f} → {mean("q99_obs"):.1f}',
                             n(mean("n_wet_obs"))])
    strengths = [r for r in rows("strengths.csv") if r["scheme"] == "balanced_spatial"]
    level_color = {"none": "var(--empty)", "threshold": "var(--q1)", "0.25": "var(--q2)", "0.5": "var(--q3)", "0.75": "var(--q4)", "1.0": "var(--q5)"}
    level_text = {"none": "·", "threshold": "T", "0.25": "¼", "0.5": "½", "0.75": "¾", "1.0": "1"}
    def level(r):
        if r["threshold_applied"] == "false":
            return "none"
        return "threshold" if float(r["strength"]) == 0 else str(float(r["strength"]))
    cells_in_order = []
    for regime in regime_order:
        for cell in sorted({r["cell"] for r in strengths if r["regime"] == regime}):
            cells_in_order.append((regime, cell))
    head_cells = "".join(f'<th title="{esc(REGIME_NAME[g])}">{esc(c[2:4] + "/" + c[5:])}</th>' for g, c in cells_in_order)
    grid_rows = []
    for fold in sorted({r["fold"] for r in strengths}):
        tds = []
        for regime, cell in cells_in_order:
            r = next(x for x in strengths if x["fold"] == fold and x["regime"] == regime and x["cell"] == cell)
            key = level(r)
            ink = "var(--on-high)" if key in ("0.75", "1.0") else "var(--fg)"
            tds.append(f'<td style="background:{level_color[key]};color:{ink}" data-tip="Fold {fold}, {cell}, {esc(REGIME_NAME[regime])}: '
                       f'{ {"none": "no adjustment", "threshold": "wet threshold only"}.get(key, "threshold + " + key + " of the map") }; '
                       f'training RMSE {float(r["train_rmse_none"]):.3f} with none, {float(r["train_rmse_full"]):.3f} with the full map, {float(r["train_rmse_chosen"]):.3f} chosen">{level_text[key]}</td>')
        grid_rows.append(f"<tr><th>Fold {esc(fold)}</th>{''.join(tds)}</tr>")
    strength_grid = (f'<div class="scroll"><table class="grid tight"><thead><tr><th></th>{head_cells}</tr></thead><tbody>{"".join(grid_rows)}</tbody></table></div>'
                     + legend([("no adjustment", level_color["none"]), ("T: wet threshold only", level_color["threshold"]), ("¼ of the map", level_color["0.25"]),
                               ("½", level_color["0.5"]), ("¾", level_color["0.75"]), ("1: full map", level_color["1.0"])]))
    counts = {key: sum(level(r) == key for r in strengths) for key in level_color}
    out.append(section("correction", "7", "The correction inside each fold",
                       "Corrected FY4B (FY4B_QM) is fitted five times, once per cross-validation fold, on that fold's training gauges only, and applied to every gauge. A held-out gauge never shapes the anchor it is scored against.", f'''
<h3>Step 1: a pixel shift per regime</h3>
<p>For each regime, every shift of the 3 × 3 box up to ±6 pixels is scored by the pooled correlation with the training gauges over all hours of the regime; the best is used for that regime's hours. Shifts are tried nearest first, so a tie keeps the smaller one. x is positive to the east, y positive to the south (image rows).</p>
{table(shift_head, shift_rows, num=(1,))}{source("…/fy4b_qm_shift_selection.csv (the benchmark run)", "output/fy4b_process_trace/shifts.csv")}
<p>The shifts agree across folds: about 2 pixels west before the May 2023 software change, one pixel north after it, and 2 pixels east and one north at 105°E. The gain in training r is small each time (third decimal).</p>
<h3>Step 2: a wet threshold and a quantile map per season and regime</h3>
<p>In each season × regime cell, the threshold τ is the FY4B value above which FY4B is wet as often as the training gauges are (≥ 0.1 mm h⁻¹); values at or below τ become 0. Above τ, FY4B's wet values are mapped onto the gauges' wet values at 199 quantiles (0.5% to 99.5%). Above the top quantile a value is scaled by the top quantiles' ratio, so the 20 mm h⁻¹ ceiling can be exceeded. A cell with fewer than 500 wet pairs on either side takes its regime's map.</p>
{table(["Regime", "Season", "Map fitted on", "τ, mm h⁻¹ (range over folds)", "Median wet value: FY4B → gauge", "99th percentile: FY4B → gauge", "Gauge wet pairs (mean)"], map_rows, num=(3, 4, 5, 6))}
{source("…/fy4b_qm_maps.csv (the benchmark run)", "output/fy4b_process_trace/maps.csv")}
<h3>Step 3: how strongly, month by month</h3>
<p>For every month of each regime, six candidates are scored by squared error on the training gauges: no adjustment, the threshold alone, and the threshold with ¼, ½, ¾ or all of the way to the mapped value. The lowest wins; a tie goes to the smaller adjustment. Of the {len(strengths)} fold-months, {counts["1.0"]} take the full map, {counts["0.75"] + counts["0.5"] + counts["0.25"]} part of it, {counts["threshold"]} the threshold alone and {counts["none"]} nothing.</p>
{strength_grid}{source("…/fy4b_qm_strengths.csv (the benchmark run)", "output/fy4b_process_trace/strengths.csv")}
<p>Columns are months (yy/mm) in regime order; May 2023 appears twice because the software changed on the 16th. The summers of 2023 and 2024 get no adjustment in every fold: there the full map raised the training error.</p>'''))

    # --- 8. In the benchmark -----------------------------------------------------------------------
    merge = [r for r in rows("fusion_merge_weights.csv") if r["scheme"] == "balanced_spatial" and r["product"] == "MERGED_OLS"]
    blend = [r for r in rows("fusion_blend_weights.csv") if r["scheme"] == "balanced_spatial"]
    scope_rows = [(k, scope[k]) for k in ("primary_cv", "secondary_cv", "repeated_cv_seeds", "common_evaluation_mask", "covariate_selection_mode",
                                           "fused_anchor", "fused_anchor_leakage", "fy4b_qm", "fy4b_qm_leakage", "fy4b_qm_adaptive",
                                           "satellite_wet_blend_selection", "blend_axes", "git_commit", "git_dirty", "julia_version", "manifest_sha256") if k in scope]
    out.append(section("benchmark", "8", "In the benchmark",
                       "One run scores every anchor and method. Plain FY4B has three roles in it; corrected FY4B has one.", f'''
<pre class="cmd">julia -t 4 --project=. {esc(scope["command"])}</pre>
<dl class="facts">
<dt>Anchors, in run order</dt><dd>FY4B, GPM (IMERG Uncal), GSMaP, the equal-weight merge, the OLS merge, and corrected FY4B last. The position sets each anchor's covariate-selection seed, so appending FY4B_QM last moves no other row.</dd>
<dt>Plain FY4B as an anchor</dt><dd>Raw, and corrected against the training gauges by residual GWR, mixed GWR and MGWR, each also blended with ADW.</dd>
<dt>Plain FY4B in the merges</dt><dd>The equal-weight merge gives it a third. The OLS merge fits an intercept and three weights on the fold's training gauges (table below).</dd>
<dt>Plain FY4B in the agreement bands</dt><dd>The blend's ADW weight is chosen per band of how many of the three products call the cell wet and the largest of their rates. FY4B's wet call is one of the three.</dd>
<dt>Corrected FY4B</dt><dd>An anchor only. It enters neither the merges nor the bands.</dd>
<dt>Validation</dt><dd>Balanced spatial 5-fold CV (primary) and random 5-fold; every choice (covariates, bandwidths, blend weights, merge weights, shift, map, strength) is made inside the training fold.</dd>
</dl>
<h3>OLS merge weights per fold</h3>
{table(["Fold", "Training gauge-hours", "Intercept", "FY4B", "GPM", "GSMaP"],
       [[esc(r["fold"]), n(r["n_train_cell"]), f(r["intercept"], 3), f(r["beta_fy4b"], 3), f(r["beta_gpm_uncal"], 3), f(r["beta_gsmap"], 3)] for r in merge], num=(1, 2, 3, 4, 5))}
{source("…/fused_anchor_selection.csv (the benchmark run)")}
<p>FY4B's weight in the OLS merge is about 0.02, against about 0.48 for GPM and 0.15 for GSMaP: the fit gives plain FY4B almost no say.</p>
<h3>Agreement-blend weights where FY4B or corrected FY4B is the anchor</h3>
{table(["Anchor", "Fold", "Weights per band (9 bands, in the run's order)", "Inner RMSE, blended", "Inner RMSE, unblended"],
       [[esc(PRODUCT_NAME[r["product"]]), esc(r["fold"]), esc(r["lambdas"].replace("|", "  ")), f(r["inner_RMSE"], 3), f(r["unblended_RMSE"], 3)] for r in blend], num=(3, 4))}
{source("…/blend_selection.csv (the benchmark run)")}
<h3>What the run recorded about itself</h3>
{table(["Key", "Value"], [[f"<code>{esc(k)}</code>", esc(v)] for k, v in scope_rows], cls="wrap")}
{source("…/benchmark_scope.csv (the benchmark run)")}'''))

    # --- 9. What the paper cites -------------------------------------------------------------------
    mu = numbers("methods_uncal")["rows"]
    def t3(product, method):
        r = next(x for x in mu if x["product"] == product and x["method"] == method)
        return [esc(PRODUCT_NAME.get(product, product) if method != "adw" else "—"), esc(METHOD_NAME[method]), r["RMSE"], r["r"], r["CSI2"],
                "—" if method == "adw" else minus(r["vs_raw"]) + "%", minus(r["vs_adw"]) + "%"]
    table3 = table(["Anchor", "Method", "RMSE", "r", "CSI 2", "vs raw", "vs ADW"],
                   [t3("GPM_uncal", "adw")] + [t3(p, m) for p in ("FY4B", "FY4B_QM", "GPM_uncal", "GSMaP", "MERGED_OLS")
                                               for m in ("raw", "mgwr", "blend_agrenv_mgwr")], num=(2, 3, 4, 5, 6))
    fq = [r for r in numbers("fy4b_quality")["rows"] if r["scheme"] == "balanced_spatial"]
    quality = table(["Corrected FY4B against", "Stage", "RMSE change", "95% interval", "r, corrected", "r, the other"],
                    [[esc(PRODUCT_NAME[r["baseline"]]), esc(METHOD_NAME[r["method"]]), minus(r["change"]) + "%",
                      f'{minus(r["lo"])} to {minus(r["hi"])}%', r["r"], r["r0"]] for r in fq], num=(2, 3, 4, 5))
    ft = [r for r in numbers("fy4b_time_holdout")["rows"] if r["scheme"] == "balanced_spatial"]
    holdout = table(["Correction fitted with", "RMSE change, raw", "95% interval", "Change in r", "95% interval"],
                    [[{"FY4B_QM_gauges": "gauges held out (as in the benchmark)", "FY4B_QM_crossed": "gauges and whole weeks held out"}.get(r["candidate"], esc(r["candidate"])),
                      minus(r["change"]) + "%", f'{minus(r["lo"])} to {minus(r["hi"])}%', minus(r["dr"]), f'{minus(r["dr_lo"])} to {minus(r["dr_hi"])}']
                     for r in ft], num=(1, 2, 3, 4))
    fe = [r for r in numbers("fy4b_events_summary")["rows"] if r["window"] == "fy4b_hours"]
    events = table(["Product", "r", "MB, mm", "RMSE, mm", "CSI 50"],
                   [[esc(PRODUCT_NAME[r["product"]]), r["r"], minus(r["MB"]), r["RMSE"], r["CSI50"]] for r in fe], num=(1, 2, 3, 4))
    fg = [r for r in numbers("fusion_configs")["rows"] if r["method"] == "blend_agrenv_mgwr"]
    configs = table(["Anchor", "RMSE, plain FY4B in the merges and bands", "RMSE, corrected FY4B there", "Change", "95% interval", "CSI 8, plain", "CSI 8, corrected"],
                    [[esc(PRODUCT_NAME[r["product"]]), r["rmse_plain"], r["rmse_qm"], minus(r["change"]) + "%", f'{minus(r["lo"])} to {minus(r["hi"])}%',
                      r["csi8_plain"], r["csi8_qm"]] for r in fg], num=(1, 2, 3, 4, 5, 6))
    clt = numbers("fy4b_clt")["rows"][0]
    clt_filter = numbers("fy4b_clt_filter")["rows"][0]
    out.append(section("paper", "9", "What the paper cites",
                       f'The numbers below are the paper\'s own, read from <code>paper/numbers/*.json</code>. All held-out scores are on the {n(scored["gauge_hours"])} gauge-hours of stage 5.', f'''
<h3>Draft Table 3: FY4B and corrected FY4B beside the other anchors</h3>
{table3}{source("paper/numbers/methods_uncal.json", "output/gauge_satellite_fusion_…_matched_block7d/method_summary.csv")}
<p class="aside">RMSE in mm h⁻¹; CSI at 2 mm h⁻¹; the last two columns are RMSE reductions (positive is better). The paper's Table 3 prints the raw and blend rows.</p>
<h3>§3.4: corrected FY4B against plain FY4B, GPM and GSMaP</h3>
{quality}{source("paper/numbers/fy4b_quality.json", "scripts/compare_fy4b_quality.jl")}
<p class="aside">Negative is better for corrected FY4B. Intervals resample whole 7-day blocks.</p>
<h3>§3.4: the same raw gain with whole weeks held out as well</h3>
{holdout}{source("paper/numbers/fy4b_time_holdout.json", "scripts/check_fy4b_qm_time_holdout.jl")}
<p>The raw gain against plain FY4B shrinks when the correction never sees the weeks it is scored on, and the gain in r is gone. The benchmark's corrected FY4B is the first row.</p>
<h3>§3.4: the cloud-type pilot, July 2023</h3>
<p>Hours dominated by ice cloud held {esc(clt["false_alarms"].lstrip("+"))}% of corrected FY4B's false alarms and {esc(clt["hits"].lstrip("+"))}% of its hits (false-alarm ratio {esc(clt["far"])}, against {esc(clt["far_gpm"])} for GPM).
A filter on FY4B's own cloud type raised CSI at 0.1 mm h⁻¹ by {minus(clt_filter["csi"])} ({minus(clt_filter["csi_lo"])} to {minus(clt_filter["csi_hi"])}) over removing as many of its lightest rain hours. Not used.</p>
{source("paper/numbers/fy4b_clt.json", "paper/numbers/fy4b_clt_filter.json")}
<h3>Draft Table 4: the six RQ1 heavy-rain events, over their shared hours</h3>
{events}{source("paper/numbers/fy4b_events_summary.json", "scripts/run_heavy_rain_event_fy4b.jl")}
<h3>Draft Table 5: plain or corrected FY4B in the merges and the bands (MGWR + agreement blend)</h3>
{configs}{source("paper/numbers/fusion_configs.json", "scripts/compare_fy4b_fusion_configs.jl")}
<p>Only the equal-weight merge gains, and its CSI at 8 mm h⁻¹ falls. Plain FY4B stays in the merges and the bands.</p>'''))

    # --- 10. Checks --------------------------------------------------------------------------------
    check_rows = [[f'<span class="pill {"ok" if r["pass"] == "true" else "bad"}">{"✓ pass" if r["pass"] == "true" else "✕ fail"}</span>',
                   esc(r["check"]), esc(r["expected"] if len(r["expected"]) < 70 else r["expected"][:12] + "…"),
                   esc(r["found"] if len(r["found"]) < 70 else r["found"][:12] + "…")] for r in checks]
    out.append(section("checks", "10", "Checks",
                       f"{passed} of {len(checks)} cross-checks pass. Each compares two independent records of the same thing; the trace script stops with an error if one fails.",
                       table(["", "Check", "Expected", "Found"], check_rows, cls="wrap checks") + source("output/fy4b_process_trace/checks.csv")))

    # --- 11. Provenance ----------------------------------------------------------------------------
    step_rows = [[f'<code>{esc(r["step"])}</code>', f'<code class="cmdcell">{esc(r["command"]) or "no script (downloaded by hand)"}</code>',
                  f'<code>{esc(r["git_commit"][:7])}</code>' if r["git_commit"] else "—",
                  "—" if not r["git_dirty"] else ("clean" if r["git_dirty"] == "false" and r["git_dirty_at_end"] in ("false", "") else "uncommitted changes"),
                  "same file" if r["fy4b_file_sha256"] else "—"] for r in steps]
    out.append(section("provenance", "11", "Provenance of each step",
                       "The steps of <code>paper/sources.toml</code> that make or read FY4B, with the commit and tree state each run recorded.",
                       table(["Step", "Command", "Commit", "Tree", "fy4b_w3_h0.csv (SHA-256)"], step_rows, cls="wrap")
                       + source("paper/sources.toml", "each step's provenance file", "output/fy4b_process_trace/steps.csv")
                       + '<p class="aside">"Same file" means the run recorded the SHA-256 of the <code>fy4b_w3_h0.csv</code> that is on disk now. Steps that read FY4B through the benchmark run\'s outputs record that run instead.</p>'))

    # --- 12. Tried and not kept ---------------------------------------------------------------------
    tries = numbers("fy4b_tries")["rows"]
    out.append(section("tries", "12", "Tried and not kept",
                       f"The {len(tries)} FY4B fixes tried on the same gauges and hours (Draft Table S1). Tries 1–3 produced the path above; the rest changed nothing in it.",
                       table(["No.", "Date", "Fix tried", "Outcome"], [[esc(r["n"]), esc(r["date"]), esc(r["fix"]), esc(r["outcome"])] for r in tries], cls="wrap")
                       + source("paper/fy4b_tries.csv", "Projects_Pending_Handover/fy4b_tries.md")))

    # --- 13. Points for the paper text ---------------------------------------------------------------
    gap = funnel["relocation_gap"]
    out.append(section("points", "13", "Two points for the paper's text",
                       "Both are wording, not numbers: every count the paper cites is reproduced above. Nothing in the paper was changed.", f'''
<div class="note"><h3>The gap at the move to 105°E is longer than February 2024</h3>
<p>§2.5 says no files exist "for February 2024". The gap runs from {esc(gap["first_hour_utc"][:10])} to {esc(gap["last_hour_utc"][:10])} (UTC), {n(gap["hours"])} hours with no file at all: February and the first four days of March.</p></div>
<div class="note"><h3>{n(v_all["no_value_at_any_gauge"])} complete hours have no value, and the strict rule is not why</h3>
<p>§2.5 says FY4B has a value at the gauges on {n(matched["hours"])} hours, and names the missing 14:15 file and "the strict rule" as the reasons for the losses. The strict rule leaves {n(complete)} hours; the other {n(v_all["no_value_at_any_gauge"])} are lost because the values over the study area are rejected in complete hours, mostly near local midnight around the equinoxes (stage 3).</p></div>'''))

    nav = "".join(f'<a href="#{anchor}">{label}</a>' for anchor, label in (
        ("path", "0 Path"), ("raw", "1 Raw files"), ("hour", "2 Strict hour"), ("sampling", "3 Value at a gauge"), ("regimes", "4 Regimes"),
        ("sample", "5 Matched sample"), ("series", "6 Against gauges"), ("correction", "7 Correction"), ("benchmark", "8 Benchmark"),
        ("paper", "9 Paper"), ("checks", "10 Checks"), ("provenance", "11 Provenance"), ("tries", "12 Tried"), ("points", "13 Points")))
    return PAGE_HEAD + out[0] + f'<nav class="toc">{nav}</nav>' + "".join(out[1:]) + PAGE_TAIL


PAGE_HEAD = """<meta name="claude:viewport-fit-sides" content="cover">
<title>FY4B process trace</title>
<style>
/* Layout: one reading column of numbered stages; wide tables and charts scroll inside their own box. */
:root {
  --bg: #f7f8f8; --surface: #ffffff; --fg: #101416; --fg2: #4a5258; --muted: #7c858b; --line: #dfe3e5; --axis: #b9c0c4;
  --accent: #1c5cab; --note: #eef3f9; --on-dark: #ffffff; --on-high: #ffffff; --empty: #eceeef;
  --s1: #2a78d6; --s2: #eb6834; --s3: #1baf7a; --s4: #eda100; --s5: #e87ba4; --s7: #6250d6;
  --q1: #cde2fb; --q2: #9ec5f4; --q3: #6da7ec; --q4: #2a78d6; --q5: #184f95;
  --ok: #006300; --ok-bg: #e4f3e4; --bad: #a82222; --bad-bg: #fbe6e6;
  --sans: "Inter", ui-sans-serif, system-ui, -apple-system, "Segoe UI", "PingFang SC", "Microsoft YaHei", sans-serif;
  --mono: ui-monospace, "SF Mono", "JetBrains Mono", Menlo, Consolas, monospace;
}
@media (prefers-color-scheme: dark) { :root:not([data-theme="light"]) {
  --bg: #0e1113; --surface: #171b1e; --fg: #eef0f1; --fg2: #c0c6ca; --muted: #8a9399; --line: #2a3035; --axis: #3b4349;
  --accent: #86b6ef; --note: #172230; --on-dark: #ffffff; --on-high: #0e1113; --empty: #22282c;
  --s1: #3987e5; --s2: #d95926; --s3: #199e70; --s4: #c98500; --s5: #d55181; --s7: #9085e9;
  --q1: #17385f; --q2: #1c5cab; --q3: #2a78d6; --q4: #5598e7; --q5: #9ec5f4;
  --ok: #5fd35f; --ok-bg: #12301a; --bad: #f08a8a; --bad-bg: #3a1717; color-scheme: dark } }
:root[data-theme="dark"] {
  --bg: #0e1113; --surface: #171b1e; --fg: #eef0f1; --fg2: #c0c6ca; --muted: #8a9399; --line: #2a3035; --axis: #3b4349;
  --accent: #86b6ef; --note: #172230; --on-dark: #ffffff; --on-high: #0e1113; --empty: #22282c;
  --s1: #3987e5; --s2: #d95926; --s3: #199e70; --s4: #c98500; --s5: #d55181; --s7: #9085e9;
  --q1: #17385f; --q2: #1c5cab; --q3: #2a78d6; --q4: #5598e7; --q5: #9ec5f4;
  --ok: #5fd35f; --ok-bg: #12301a; --bad: #f08a8a; --bad-bg: #3a1717; color-scheme: dark }
* { box-sizing: border-box }
body { background: var(--bg); color: var(--fg); font: 15px/1.55 var(--sans); padding-block: 28px 64px;
  padding-left: max(16px, env(safe-area-inset-left, 0px)); padding-right: max(16px, env(safe-area-inset-right, 0px)) }
.top, .toc, section { max-width: 1040px; margin-inline: auto }
h1 { font-size: 34px; line-height: 1.1; margin: 4px 0 10px; letter-spacing: -0.02em; text-wrap: balance }
h2 { font-size: 22px; margin: 0; letter-spacing: -0.01em; text-wrap: balance }
h3 { font-size: 15px; margin: 26px 0 8px }
p { margin: 10px 0; max-width: 76ch }
code, pre { font-family: var(--mono); font-size: 0.88em }
code { overflow-wrap: anywhere }
.eyebrow { font-size: 12px; letter-spacing: 0.08em; text-transform: uppercase; color: var(--muted); margin: 0 }
.lead { color: var(--fg2); font-size: 16px; max-width: 78ch }
.meta { display: flex; flex-wrap: wrap; gap: 8px 18px; align-items: center; color: var(--fg2); font-size: 13px; margin: 14px 0 }
.pill { display: inline-block; padding: 2px 9px; border-radius: 999px; font-size: 12.5px; font-weight: 600; white-space: nowrap }
.pill.ok { background: var(--ok-bg); color: var(--ok) } .pill.bad { background: var(--bad-bg); color: var(--bad) }
.tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 1px; background: var(--line);
  border: 1px solid var(--line); border-radius: 8px; overflow: hidden; margin-top: 18px }
.tiles div { background: var(--surface); padding: 14px 16px; min-width: 0 }
.tiles b { display: block; font: 600 26px/1.15 var(--mono); font-variant-numeric: tabular-nums; letter-spacing: -0.02em }
.tiles span { color: var(--fg2); font-size: 13px }
.toc { display: flex; flex-wrap: wrap; gap: 4px 14px; margin-top: 22px; padding: 10px 0; border-block: 1px solid var(--line); font-size: 13px }
.toc a { color: var(--accent); text-decoration: none; white-space: nowrap } .toc a:hover, .toc a:focus-visible { text-decoration: underline }
section { margin-top: 52px }
header.sec { display: flex; align-items: baseline; gap: 12px; border-bottom: 2px solid var(--fg); padding-bottom: 6px }
.no { font: 600 14px var(--mono); color: var(--muted); min-width: 1.6em }
.flow { list-style: none; margin: 16px 0 0; padding: 0; counter-reset: stage; border-top: 1px solid var(--line) }
.flow li { counter-increment: stage; display: grid; grid-template-columns: 2.2em minmax(0, 1fr) minmax(0, 16em); gap: 4px 14px;
  padding: 12px 0; border-bottom: 1px solid var(--line); align-items: start }
.flow li::before { content: counter(stage); font: 600 14px var(--mono); color: var(--muted); padding-top: 2px }
.flow h3 { margin: 0 0 2px } .flow p { margin: 0 0 4px; color: var(--fg2); font-size: 14px } .flow code { color: var(--muted); font-size: 12px }
.flow .count { font: 600 14px/1.4 var(--mono); font-variant-numeric: tabular-nums; text-align: right }
@media (max-width: 640px) { .flow li { grid-template-columns: 2.2em minmax(0, 1fr) } .flow .count { grid-column: 2; text-align: left } h1 { font-size: 27px } }
.facts { display: grid; grid-template-columns: minmax(0, 11em) minmax(0, 1fr); gap: 8px 18px; margin: 16px 0 }
.facts dt { font-weight: 600 } .facts dd { margin: 0; color: var(--fg2); max-width: 74ch }
@media (max-width: 640px) { .facts { grid-template-columns: minmax(0, 1fr) } .facts dd { margin-bottom: 8px } }
.scroll { overflow-x: auto; margin: 10px 0 4px }
table { border-collapse: collapse; font-size: 13.5px; font-variant-numeric: tabular-nums; min-width: 100% }
th, td { text-align: left; padding: 6px 12px 6px 0; border-bottom: 1px solid var(--line); vertical-align: top; white-space: nowrap }
thead th { font-size: 12px; color: var(--fg2); font-weight: 600; border-bottom: 1px solid var(--axis); white-space: normal; vertical-align: bottom }
.num { text-align: right; font-family: var(--mono); font-size: 12.5px }
table.wrap td, table.wrap th { white-space: normal } table.wrap td:last-child { min-width: 14em }
table.wrap td:first-child, table.wrap td:first-child code { white-space: nowrap; overflow-wrap: normal } table.checks td:nth-child(2) { min-width: 24em }
.cmdcell { display: block; min-width: 22em; max-width: 46em; color: var(--fg2) }
table.grid { min-width: 0; border-collapse: separate; border-spacing: 2px }
table.grid th, table.grid td { border: 0; padding: 0; text-align: center; font: 12px var(--mono) }
table.grid thead th { color: var(--muted); font-weight: 500; padding-bottom: 2px }
table.grid tbody th { text-align: right; padding-right: 8px; color: var(--fg2); white-space: nowrap }
table.grid td { width: 46px; height: 30px; border-radius: 3px } table.grid.tight td { width: 30px; min-width: 30px }
table.grid.tight thead th { font-size: 10px; writing-mode: vertical-rl; padding: 0 0 4px; height: 44px }
.sw { display: inline-block; width: 11px; height: 11px; border-radius: 2px; vertical-align: -1px }
.legend { display: flex; flex-wrap: wrap; gap: 4px 16px; font-size: 12.5px; color: var(--fg2); margin: 8px 0 }
.legend i { display: inline-block; width: 11px; height: 11px; border-radius: 2px; margin-right: 6px; vertical-align: -1px }
.chart { display: block; width: 100%; height: auto } .chart.wide { height: auto }
.chart .grid { stroke: var(--line); stroke-width: 1 } .chart .axis { stroke: var(--axis); stroke-width: 1 }
.chart .tick { fill: var(--muted); font: 11px var(--mono) } .chart .lab { fill: var(--fg2); font: 13px var(--sans) }
.chart .val { fill: var(--fg); font: 12px var(--mono) }
.axis-label { font-size: 12.5px; color: var(--fg2); margin: 14px 0 2px }
.meta a { color: var(--accent) }
.src { font-size: 12px; color: var(--muted); margin: 2px 0 0 } .aside { font-size: 13.5px; color: var(--fg2) }
.note { background: var(--note); border-radius: 6px; padding: 14px 18px; margin: 18px 0 } .note h3 { margin: 0 0 4px } .note p { margin: 4px 0 0 }
pre.cmd { background: var(--surface); border: 1px solid var(--line); border-radius: 6px; padding: 12px 14px; overflow-x: auto; white-space: pre-wrap; overflow-wrap: anywhere }
[data-tip] { cursor: default }
#tip { position: fixed; z-index: 20; max-width: 320px; background: var(--fg); color: var(--bg); font-size: 12.5px; line-height: 1.4;
  padding: 6px 9px; border-radius: 5px; pointer-events: none }
a:focus-visible { outline: 2px solid var(--accent); outline-offset: 2px }
</style>
"""

PAGE_TAIL = """<div id="tip" hidden></div>
<script>
(function () {
  var tip = document.getElementById("tip");
  function show(event) {
    var target = event.target.closest ? event.target.closest("[data-tip]") : null;
    if (!target) { tip.hidden = true; return; }
    tip.textContent = target.getAttribute("data-tip");
    tip.hidden = false;
    var x = Math.min(event.clientX + 14, window.innerWidth - tip.offsetWidth - 8);
    var y = event.clientY + 16;
    if (y + tip.offsetHeight > window.innerHeight - 8) { y = event.clientY - tip.offsetHeight - 10; }
    tip.style.left = Math.max(8, x) + "px";
    tip.style.top = Math.max(8, y) + "px";
  }
  document.addEventListener("pointermove", show);
  document.addEventListener("pointerdown", show);
  document.addEventListener("scroll", function () { tip.hidden = true; }, true);
})();
</script>
"""


def main() -> None:
    page = build()
    PAGE.write_text(page, encoding="utf-8")
    inputs = tuple(sorted(TRACE.glob("*.csv"))) + tuple(sorted(MANIFEST.glob("folders.csv")))
    inputs = tuple(p for p in inputs if p.name != "plot_provenance.csv")
    write_plot_provenance(TRACE / "page", ROOT, inputs)
    print(f"wrote {PAGE.relative_to(ROOT)} ({len(page) / 1024:.0f} KB)")


if __name__ == "__main__":
    main()

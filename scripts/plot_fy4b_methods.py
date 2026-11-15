"""Draw the FY4B methods page: every step and formula from a raw 15-minute file to the benchmark's scores.

    $MIXEDGWR_PYTHON scripts/plot_fy4b_methods.py     # -> output/fy4b_process_trace/fy4b_methods.html

The companion of scripts/plot_fy4b_process.py (counts) for the question "how": what the code does to
FY4B, written as formulas. The page follows the codebase, not the paper: each step names the function
it describes, and the script stops if a named function is no longer in the named file. A step is
marked "checked by number" only when test/test-fy4b-methods-formulas.jl holds a testset that compares
the formula, written out plainly, with that function; the script stops if the testset is missing.

The text and the formulas live here. The values under each step are read from the trace tables
(scripts/trace_fy4b_process.jl) and from the benchmark run those tables point to.
"""

from __future__ import annotations

import math
import os
import re
import sys
from pathlib import Path

from plot_fy4b_process import PAGE_HEAD, PAGE_TAIL, REGIME_NAME, ROOT, TRACE, esc, f, n, rows, source, table
from plot_provenance import write_plot_provenance

PAGE = TRACE / "fy4b_methods.html"
# The same page without the document shell, for publishing where a shell is added (KaTeX needs a
# doctype, so the local file carries its own).
FRAGMENT = TRACE / "methods_page" / "fy4b_methods.html"
# Where the companion page is: beside this one on disk; a published copy passes its own address.
PROCESS_LINK = os.environ.get("FY4B_PROCESS_URL", "fy4b_process.html")
SHELL = '<!doctype html>\n<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"></head><body>\n{}\n</body></html>\n'
TEST_FILE = ROOT / "test" / "test-fy4b-methods-formulas.jl"
KATEX = "https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.16.9/katex.min.js"
KERNELS = ["Gaussian", "exponential", "bisquare", "tricube", "boxcar"]
PRODUCT_LABEL = {"FY4B": "FY4B", "FY4B_QM": "Corrected FY4B"}

_cited: list[tuple[str, str]] = []
_tests = TEST_FILE.read_text(encoding="utf-8")


def md(text: str) -> str:
    """Paragraphs, `code`, **bold**, $inline$ and $$display$$ math, and '- ' lists, as HTML."""
    out = []
    for block in re.split(r"\n\s*\n", text.strip()):
        block = block.strip()
        if block.startswith("$$") and block.endswith("$$"):
            out.append(f'<div class="scroll"><div class="math">{esc(block[2:-2].strip())}</div></div>')
            continue
        html = esc(block)
        html = re.sub(r"\$(.+?)\$", lambda m: f'<span class="math">{m.group(1)}</span>', html, flags=re.S)
        html = re.sub(r"`(.+?)`", r"<code>\1</code>", html)
        html = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", html)
        if html.startswith("- "):
            items = "".join(f"<li>{item}</li>" for item in re.split(r"\n- ", html[2:]))
            out.append(f"<ul>{items}</ul>")
        else:
            out.append(f"<p>{html}</p>")
    return "".join(out)


def where(*refs: tuple[str, str]) -> str:
    """The functions a step describes, as `file`: `name`; each must exist."""
    for path, name in refs:
        text = (ROOT / path).read_text(encoding="utf-8")
        if name not in text:
            sys.exit(f"{name} is not in {path}")
        _cited.append((path, name))
    return '<p class="where">Code: ' + "; ".join(f"<code>{esc(name)}</code> in <code>{esc(path)}</code>"
                                                  for path, name in refs) + "</p>"


def checked(*testsets: str) -> str:
    if not testsets:
        return '<span class="pill warn">read from the code, not checked by number</span>'
    for name in testsets:
        if f'@testset "formula: {name}"' not in _tests:
            sys.exit(f'no testset "formula: {name}" in {TEST_FILE.name}')
    return '<span class="pill ok">✓ checked by number</span> <span class="src">testset: ' + ", ".join(
        f"“{esc(name)}”" for name in testsets) + "</span>"


def step(number: int, title: str, refs: str, status: str, body: str) -> str:
    return (f'<article class="step" id="s{number}"><h3><span class="no">{number}</span>{esc(title)}</h3>'
            f'<div class="stepmeta">{status}</div>{body}{refs}</article>')


def part(anchor: str, letter: str, title: str, lead: str, body: str) -> str:
    return (f'<section id="{anchor}"><header class="sec"><span class="no">{letter}</span><h2>{esc(title)}</h2></header>'
            f'<p class="lead">{lead}</p>{body}</section>')


def navigation(lat: float, lon: float, sat_lon: float) -> tuple[float, float]:
    """The page's navigation formula (step 2), for the worked example."""
    r_e, r_p, big_h, size, res = 6378.137e3, 6356.75231414e3, 42164.0e3, 2748, 4000.0
    phi, lam = math.radians(lat), math.radians(lon - sat_lon)
    flat = (r_e - r_p) / r_e
    e2 = flat * (2 - flat)
    phi_c = math.atan((1 - e2) * math.tan(phi))
    r_c = r_p / math.sqrt(1 - e2 * math.cos(phi_c) ** 2)
    x, y, z = r_c * math.cos(phi_c) * math.cos(lam), r_c * math.cos(phi_c) * math.sin(lam), r_c * math.sin(phi_c)
    delta = res / (big_h - r_e)
    return ((size - 1) / 2 + math.atan(y / (big_h - x)) / delta,
            (size - 1) / 2 - math.atan(z / math.hypot(big_h - x, y)) / delta)


def bw_text(row: dict[str, str]) -> str:
    if row["bw"] == "Inf":
        return "global"
    return f'{float(row["bw"]):.0f} nearest' if row["adaptive"] == "true" else f'{float(row["bw"]):.0f} km'


def build() -> str:
    steps = {r["step"]: r for r in rows("steps.csv")}
    run_dir = (ROOT / steps["benchmark"]["marker"]).parent
    run = lambda name: rows(name, run_dir)
    scope = {r["key"]: r["value"] for r in rows("benchmark_scope.csv")}
    regimes = rows("regimes.csv")
    regime_order = [r["regime"] for r in regimes]
    out: list[str] = []

    # --- Header and notation --------------------------------------------------------------------
    n_tests = len(re.findall(r'@testset "formula: ', _tests))
    out.append(f'''<header class="top">
<p class="eyebrow">Satellite precipitation correction · FY-4B AGRI QPE · interpolation benchmark</p>
<h1>FY4B methods</h1>
<p class="lead">What the code does to FY4B, from a raw 15-minute file to its scores in the interpolation benchmark:
every step, every method, every formula. The page follows the codebase. Each step names the function it describes,
and a formula is marked as checked only when a test compares it, written out plainly, with that function.</p>
<div class="meta"><span class="pill ok">✓ {n_tests} formula tests pass</span>
<span>Test: <code>test/test-fy4b-methods-formulas.jl</code></span>
<span>Values shown: the run at <code>{esc(scope["git_commit"][:7])}</code></span>
<span>Counts and checks of the data: <a href="{esc(PROCESS_LINK)}">FY4B process trace</a></span></div></header>''')

    notation = [
        ("i, j", "a gauge (station); 237 in the study area"),
        ("t", "an hour, labelled by its end in Beijing time"),
        ("q", "one of the four 15-minute scans of an hour"),
        ("G_{i,t}", "gauge rain, mm in the hour"),
        ("F_{i,t}", "plain FY4B at the gauge, mm h⁻¹ (steps 1–6)"),
        ("\\tilde F_{i,t}", "corrected FY4B, FY4B_QM (steps 10–13)"),
        ("S_{i,t}", "the anchor a method corrects: a satellite product or a merge"),
        ("\\mathcal T_f,\\ \\mathcal V_f", "training and held-out gauges of fold f"),
        ("d_{ij}", "great-circle distance between two gauges, km"),
        ("b,\\ k", "a kernel bandwidth in km, or as a count of nearest gauges"),
        ("\\hat y_{i,t}", "a method's prediction at a held-out gauge"),
        ("\\mathbb 1[\\cdot]", "1 when the condition holds, else 0"),
    ]
    out.append(part("notation", "·", "Notation", "Symbols used throughout. Rain rates are in mm h⁻¹; an hourly mean rate equals the hour's accumulation in mm.",
                    table(["Symbol", "Meaning"], [[f'<span class="math">{esc(sym)}</span>', esc(text)] for sym, text in notation])))

    # === A. Raw file to hourly value ===============================================================
    a: list[str] = []
    a.append(step(1, "Choosing the files", where(("src/sources/FY4BPreprocessing.jl", "find_nc_files"),
                                                   ("src/sources/FY4BPreprocessing.jl", "is_15min_file"),
                                                   ("src/sources/FY4BPreprocessing.jl", "segment_problem")), checked(), md(r"""
The raw tree is walked once. A file is a 15-minute scan when its name carries a start time $t_0$ and an end time $t_1$ with

$$t_0 \in \{{:}00,\ {:}15,\ {:}30,\ {:}45\} \ \text{(seconds 0)}, \qquad t_1 = t_0 + 15\ \text{min} - 1\ \text{s}.$$

The satellite longitude $\lambda_s$ (133° or 105°E) is read from the name (`N_DISK_1330E`). Files are grouped by $(t_0, \lambda_s)$. When a scan has several copies, the paths are sorted and the first copy that opens as a 2748 × 2748 `Precipitation` field is kept; a lone copy is kept as it is. Names starting with `._` (macOS companions) are skipped. The gauge-corrected 1-, 3- and 24-hour files fail the 15-minute test and are never read.

This step has no formula to check; the trace page checks its result against the drive (names and sizes of all 84,630 scans).
""")))

    stations = rows("station_meta.csv", ROOT / "data" / "processed" / "study_area")[:3]
    nav_rows = []
    for s in stations:
        lat, lon = float(s["lat"]), float(s["lon"])
        c1, r1 = navigation(lat, lon, 133.0)
        c2, r2 = navigation(lat, lon, 105.0)
        nav_rows.append([esc(s["station_id"]), f(lon, 3), f(lat, 3), f(c1, 2), f(r1, 2), f"{round(c1) + 1}, {round(r1) + 1}",
                         f(c2, 2), f(r2, 2), f"{round(c2) + 1}, {round(r2) + 1}"])
    a.append(step(2, "Navigation: from latitude and longitude to a pixel",
                  where(("src/sources/FY4BPreprocessing.jl", "latlon_to_scan_angles"),
                        ("src/sources/FY4BPreprocessing.jl", "scan_angles_to_grid")), checked("navigation"), md(r"""
A gauge at geodetic latitude $\varphi$ and longitude $\lambda_g$ is placed in the full-disk image of a satellite at longitude $\lambda_s$. With $\lambda = \lambda_g - \lambda_s$, the WGS84 radii $r_e = 6\,378\,137$ m and $r_p = 6\,356\,752.314$ m, and $e^2 = f(2-f)$, $f = (r_e - r_p)/r_e$:

$$\varphi_c = \arctan\!\big((1-e^2)\tan\varphi\big), \qquad r_c = \frac{r_p}{\sqrt{1 - e^2\cos^2\varphi_c}}$$

$$x = r_c\cos\varphi_c\cos\lambda, \qquad y = r_c\cos\varphi_c\sin\lambda, \qquad z = r_c\sin\varphi_c$$

These are Earth-centred coordinates with the $x$ axis through the sub-satellite point. The satellite sits at distance $H = 42\,164$ km on that axis. The two scan angles are

$$\alpha_x = \arctan\frac{y}{H - x}, \qquad \alpha_y = \arctan\frac{z}{\sqrt{(H-x)^2 + y^2}}$$

and one pixel spans $\Delta = 4000 / (H - r_e)$ rad (about 111.8 µrad), so on the $N = 2748$ grid

$$c = \frac{N-1}{2} + \frac{\alpha_x}{\Delta}, \qquad r = \frac{N-1}{2} - \frac{\alpha_y}{\Delta}.$$

$c$ grows to the east and $r$ to the south; both count from 0. The pixel read is $(\operatorname{round}(c)+1,\ \operatorname{round}(r)+1)$ in 1-based indexing.

**Parallax (built, not used).** A cloud top at height $h$ above the gauge is found by adding $h(\cos\varphi\cos\lambda,\ \cos\varphi\sin\lambda,\ \sin\varphi)$ to $(x, y, z)$ before the scan angles. The series used has $h = 0$.

**A point is not visible** when $x(x-H) + y^2 + z^2 > 0$ and $(H-x)/\sqrt{(H-x)^2+y^2+z^2} < 0.156$; no gauge of the study area is near that limit.

**My reading, not a statement of the code:** this is the scan-angle form used for GOES-R, with a small-angle pixel spacing. It is not a published FY-4 column/line formula with its own scaling constants. Step 11 then measures a displacement of 1 to 3 pixels against the gauges; part of it may come from this approximation, part from the product itself. The code's only test of the navigation is against gauges.
""") + "<h4>The first three gauges of the station file</h4>" + table(
        ["Gauge", "Lon", "Lat", "c at 133°E", "r at 133°E", "Pixel (1-based)", "c at 105°E", "r at 105°E", "Pixel (1-based)"],
        nav_rows, num=(1, 2, 3, 4, 5, 6, 7, 8)) + source("data/processed/study_area/station_meta.csv", "computed here with the formula above")))

    a.append(step(3, "Which pixel values are accepted", where(("src/sources/FY4BPreprocessing.jl", "accepts_value")),
                  checked("pixel rule and 3x3 mean"), md(r"""
A raw value $P$ with quality flag $Q$ (`DQF`) is accepted when

$$P \neq 65534, \qquad 0 \le P \le 30, \qquad Q < 2,$$

and neither is missing. The product's valid range is 0–20 mm h⁻¹; the bound of 30 is what removes the sentinels 65533 and 65535, which the files do not declare as fill values.
""")))

    a.append(step(4, "The 3 × 3 mean of one scan", where(("src/sources/FY4BPreprocessing.jl", "_window_value"),
                                                           ("src/sources/FY4BPreprocessing.jl", "sample_variants")),
                  checked("pixel rule and 3x3 mean"), md(r"""
Let $\mathcal A_{i,q}$ be the accepted pixels among the nine of the 3 × 3 box centred on gauge $i$'s pixel in scan $q$. The scan's value is their mean, and it exists only when more than half are accepted:

$$v_{i,q} = \frac{1}{|\mathcal A_{i,q}|}\sum_{p \in \mathcal A_{i,q}} P_p \quad \text{if } |\mathcal A_{i,q}| \ge 5, \qquad \text{else undefined.}$$

The share of accepted pixels at the product's ceiling is kept beside it:

$$\kappa_{i,q} = \frac{|\{p \in \mathcal A_{i,q} : P_p \ge 20\}|}{|\mathcal A_{i,q}|}.$$

The same function gives the other variants built in the same pass: the centre pixel alone (window 1, needs 1 pixel), a 5 × 5 box (needs 13 of 25), and the 3 × 3 box at the parallax-shifted pixel.
""")))

    a.append(step(5, "The strict hour", where(("src/sources/FY4BPreprocessing.jl", "aggregate_hourly_variants"),
                                                ("src/sources/FY4BPreprocessing.jl", "build_hourly_qc")), checked(), md(r"""
An hour is processed only when all four scans (:00, :15, :30, :45) exist and open. Then

$$F_{i,t} = \frac{1}{4}\sum_{q=1}^{4} v_{i,q}, \qquad \text{undefined if any } v_{i,q} \text{ is undefined.}$$

The scans hold rain rates in mm h⁻¹, so their mean over the hour is the hour's accumulation in mm: $\sum_q v_{i,q}\cdot 0.25\ \text{h}$. The label of the hour is its end in Beijing time, $t = t_0^{\text{UTC}} + 1\ \text{h} + 8\ \text{h}$. The hour's satellite longitude and software version are taken from its first scan; they define the regimes of step 10.

The arithmetic is a mean of four numbers and has no test of its own; the trace page checks that every kept hour has exactly the four quarters.
""")))

    a.append(step(6, "The offset grid kept for the shift", where(("src/sources/FY4BPreprocessing.jl", "sample_variants"),
                                                                   ("scripts/prepare_fy4b_rebuild.jl", "OFFSET_GRID")), checked(), md(r"""
In the same pass the 3 × 3 value is also computed with the box moved by every whole-pixel shift

$$(\delta_x, \delta_y) \in \{-6,\dots,6\}^2 \quad (169\ \text{shifts}), \qquad F_{i,t}^{(\delta_x,\delta_y)} = \tfrac14\sum_q v_{i,q}\big(c_i + \delta_x,\ r_i + \delta_y\big),$$

with $\delta_x$ positive to the east and $\delta_y$ positive to the south. They are written to `fy4b_offset_grid.nc`. No gauge reading is used: this only tabulates the candidates from which step 11 chooses.
""")))
    out.append(part("raw", "A", "From a raw file to an hourly value at a gauge",
                    "Six steps turn four full-disk images into one number per gauge and hour. None of them reads a gauge value, so the series can enter cross-validation without leaking.",
                    "".join(a)))

    # === B. What the benchmark loads ================================================================
    b = step(7, "Products, hours and the two matrices", where(("scripts/run_interpolation_benchmark.jl", "--fy4b-file"),
                                                                ("src/benchmark/InterpolationBenchmarkFold.jl", "run_seed")), checked(), md(r"""
The benchmark is started with FY4B switched on and pointed at the rebuilt series:

`--with-fy4b --fy4b-file output/fy4b_rebuild_2022_2024/fy4b_w3_h0.csv --fy4b-qm --fy4b-qm-adaptive --fused-anchor --satellite-wet-blend --blend-axis agreement_envelope --nested-covariates`

- **Products, in this order:** FY4B, GPM (IMERG Uncal), GSMaP, the equal-weight merge `MERGED_MEAN`, the OLS merge `MERGED_OLS`, and `FY4B_QM` last. A product's position sets the seed of its covariate selection, so a product appended last moves no other product's rows.
- **Hours:** the hours every loaded product file has. With FY4B in the list these are FY4B's hours.
- **Matrices:** gauges $G_{i,t}$ and each anchor $S_{i,t}$, stations × hours, with NaN where a value is missing. Every method works hour by hour on these.

Each anchor is then run through the same seven methods (steps 14–23), fold by fold.
"""))
    out.append(part("load", "B", "What the benchmark loads", "How FY4B enters the run. No formula here, only the switches and the order.", b))

    # === C. Cross-validation ==========================================================================
    split = rows("split_common.csv", run_dir / "balanced_spatial") if (run_dir / "balanced_spatial" / "split_common.csv").is_file() else []
    fold_sizes = ""
    if split:
        key = "fold" if "fold" in split[0] else list(split[0].keys())[-1]
        counts: dict[str, int] = {}
        for r in split:
            counts[r[key]] = counts.get(r[key], 0) + 1
        fold_sizes = "<h4>Gauges per fold in the run</h4>" + table(["Fold", "Held-out gauges", "Training gauges"],
                                                                    [[esc(k), n(v), n(len(split) - v)] for k, v in sorted(counts.items())], num=(1, 2))
    c = step(8, "Balanced spatial folds, and the split inside each fold",
             where(("src/benchmark/InterpolationBenchmarkFolds.jl", "split_stations_balanced_spatial_kfold"),
                   ("src/benchmark/InterpolationBenchmarkFolds.jl", "_hilbert_centers"),
                   ("src/benchmark/InterpolationBenchmarkFolds.jl", "_capacity_assignment"),
                   ("src/benchmark/InterpolationBenchmarkFolds.jl", "selection_folds"),
                   ("src/models/TraditionalInterpolation.jl", "local_km_coordinates"),
                   ("src/core/distance.jl", "haversine_distance_matrix")), checked("distances and local coordinates"), md(r"""
**Coordinates.** Gauges are projected to km about the network's mean position $(\lambda_0, \varphi_0)$, with $R = 6371.0088$ km:

$$X_i = (\lambda_i - \lambda_0)\cos\varphi_0\,\frac{\pi}{180}\,R, \qquad Y_i = (\varphi_i - \varphi_0)\,\frac{\pi}{180}\,R.$$

**Starting centres.** Each axis is scaled to the integers $0 \dots 2^{16}-1$, each gauge gets its index along a Hilbert curve of order 16, and the gauges are sorted by that index (ties by station id). The sorted list is cut into $k = 5$ runs of $\lfloor n/k \rfloor$ or $\lfloor n/k \rfloor + 1$ gauges; the mean position of each run is a starting centre. No random number is used.

**Balanced assignment.** Fold sizes are fixed in advance ($\lfloor n/k\rfloor$, the first $n \bmod k$ folds one larger). With squared distances $D_{ic}$ from gauge $i$ to centre $c$, gauges are taken in order of decreasing certainty

$$u_i = D_{i,(2)} - D_{i,(1)} \quad \text{(second-nearest minus nearest centre)},$$

and each goes to its nearest centre that still has room. Centres are then set to the mean of their gauges, and the assignment is repeated until it stops changing (at most 100 rounds).

**The random scheme** shuffles the station ids with seed 20260627 and deals them into five folds in turn.

**Inside a fold.** Every tuning choice (steps 11–13 and 20–23) is scored on an inner split of the fold's training gauges, built by the same procedure with $k = 5$, so a tuning score is the same kind of prediction as the reported one: at gauges held out as a block.

**Distance to the nearest training gauge** (a reporting stratum) is the great-circle distance

$$d_{ij} = 2R\,\arcsin\sqrt{\sin^2\frac{\varphi_j-\varphi_i}{2} + \cos\varphi_i\cos\varphi_j\,\sin^2\frac{\lambda_j-\lambda_i}{2}}$$

minimised over the training gauges. The same $d_{ij}$ is used by every method below.

The Hilbert ordering and the assignment loop are described from the code and have no test in the formula file; the two coordinate formulas do.
""") + fold_sizes)
    out.append(part("cv", "C", "Cross-validation", "A method is scored at gauges it never saw. Five folds; each fold's gauges are held out together as a compact block.", c))

    # === D. Corrected FY4B ===============================================================================
    d: list[str] = []
    cells = sorted({(r["regime"], r["cell"]) for r in rows("strengths.csv")})
    d.append(step(9, "Why the correction sits inside the fold", where(("src/models/SatelliteCalibration.jl", "calibrate_fold")), checked(), md(r"""
Corrected FY4B, $\tilde F$ (`FY4B_QM`), is plain FY4B after a pixel shift and a quantile map. Both are fitted to gauges, so they are fitted once per fold $f$ on the training gauges $\mathcal T_f$ only and then applied to every gauge:

$$\tilde F^{(f)}_{i,t} = a_{c^\star}\Big(F^{(\delta^\star)}_{i,t}\Big), \qquad \delta^\star,\ \text{map},\ c^\star \ \text{fitted on } \{(i,t): i \in \mathcal T_f\}.$$

A held-out gauge is scored against $\tilde F^{(f)}$ of the fold that held it out, so its own readings never shaped the anchor. Steps 10–13 define the three pieces.
""")))
    d.append(step(10, "Regimes and cells", where(("src/models/SatelliteCalibration.jl", "product_regimes"),
                                                   ("src/models/SatelliteCalibration.jl", "calibration_cells"),
                                                   ("src/models/SatelliteCalibration.jl", "strength_cells")), checked(), md(r"""
Each hour belongs to a **regime** $g(t)$: one per distinct pair (satellite longitude, software version) found in the files. Seasons are $s(t) \in \{\text{DJF}, \text{MAM}, \text{JJA}, \text{SON}\}$ by calendar month. Two kinds of cell are used:

- a **map cell** is a regime and a season: $\text{cell}(t) = 4\,(g(t) - 1) + s(t)$;
- a **month cell** is a regime and a calendar month of one year, for example 2023-07. A month in which the regime changes gives two cells.
""") + table(["Regime", "First hour", "Last hour", "Hours"], [[esc(REGIME_NAME[r["regime"]]), esc(r["first_hour"].replace("T", " ")[:16]),
                                                              esc(r["last_hour"].replace("T", " ")[:16]), n(r["hours"])] for r in regimes], num=(3,))
        + f'<p class="aside">The run has {len(regimes)} regimes, {len({(r["regime"], r["season"]) for r in rows("maps.csv")})} map cells and {len(cells)} month cells.</p>'))

    shifts = [r for r in rows("shifts.csv") if r["scheme"] == "balanced_spatial"]
    shift_rows = []
    for fold in sorted({r["fold"] for r in shifts}):
        cells_row = [esc(fold)]
        for regime in regime_order:
            r = next(x for x in shifts if x["fold"] == fold and x["regime"] == regime)
            cells_row += [f'{int(r["dx"]):+d}, {int(r["dy"]):+d}'.replace("-", "−"), f'{f(r["r_train_unshifted"], 3)} → {f(r["r_train"], 3)}']
        shift_rows.append(cells_row)
    d.append(step(11, "The pixel shift", where(("src/models/SatelliteCalibration.jl", "shift_sums"),
                                                 ("src/models/SatelliteCalibration.jl", "_pooled_r"),
                                                 ("src/models/SatelliteCalibration.jl", "choose_shift")),
                  checked("pooled correlation and the shift"), md(r"""
For each shift $\delta$ of step 6, each gauge and each regime, six sums are accumulated once over the hours where both values exist, with $x = F^{(\delta)}_{i,t}$ and $y = G_{i,t}$:

$$n,\quad \textstyle\sum x,\quad \sum y,\quad \sum x^2,\quad \sum y^2,\quad \sum xy.$$

A fold adds up the rows of its training gauges and gets the pooled correlation of the regime:

$$r_g(\delta) = \frac{\tfrac1n\sum xy - \bar x\,\bar y}{\sqrt{\big(\tfrac1n\sum x^2 - \bar x^2\big)\big(\tfrac1n\sum y^2 - \bar y^2\big)}}, \qquad \bar x = \tfrac1n\textstyle\sum x,\ \ \bar y = \tfrac1n\sum y.$$

The shift used for the regime's hours is

$$\delta^\star_g = \operatorname*{arg\,max}_{\delta}\ r_g(\delta),$$

with the shifts tried in order of increasing $\delta_x^2 + \delta_y^2$ and replaced only by a strictly larger $r$, so a tie keeps the smaller shift and $(0,0)$ wins when nothing beats it.
""") + "<h4>Shift (x, y) chosen in the run, with the training correlation before → after</h4>"
        + table(["Fold"] + [x for g in regime_order for x in (REGIME_NAME[g], "r")], shift_rows)
        + source("fy4b_qm_shift_selection.csv of the run")))

    maps = [r for r in rows("maps.csv") if r["scheme"] == "balanced_spatial" and r["fold"] == "1"]
    d.append(step(12, "The wet threshold and the quantile map", where(("src/models/SatelliteCalibration.jl", "fit_quantile_map"),
                                                                        ("src/models/SatelliteCalibration.jl", "_fit_cell"),
                                                                        ("src/models/SatelliteCalibration.jl", "_map_value")),
                  checked("wet threshold and quantile map"), md(r"""
In one map cell, take all training pairs $(x, y) = (F^{(\delta^\star)}_{i,t},\ G_{i,t})$ with both values present.

**Threshold.** The gauges are wet in a share $w = \operatorname{mean}\,\mathbb 1[y \ge 0.1]$ of the pairs. FY4B is given the same wet frequency:

$$\tau = Q_x(1 - w),$$

the $(1-w)$ quantile of the FY4B values (quantiles are interpolated linearly between order statistics, Julia's default).

**Knots.** With probabilities $p_j = 0.005\,j$, $j = 1,\dots,199$:

$$q^x_j = Q\big(\{x : x > \tau\},\ p_j\big), \qquad q^y_j = Q\big(\{y : y \ge 0.1\},\ p_j\big).$$

**The map.**

$$m(x) = \begin{cases} 0 & x \le \tau \\[2pt] q^y_{199}\,\dfrac{x}{q^x_{199}} & x \ge q^x_{199} \\[6pt] q^y_1 & x \le q^x_1 \\[2pt] q^y_j + \dfrac{x - q^x_j}{q^x_{j+1} - q^x_j}\,\big(q^y_{j+1} - q^y_j\big) & q^x_j \le x < q^x_{j+1} \end{cases}$$

The cases are tested in this order. Above the top knot a value is scaled by the top knots' ratio, so FY4B's 20 mm h⁻¹ ceiling becomes a lower bound of the largest values, not a ceiling.

**Sparse cells.** A cell keeps its own map only when both wet sets hold at least 500 pairs. Otherwise it takes the map fitted on all hours of its regime, and if that is also sparse, the map fitted on all hours.
""") + "<h4>Fold 1 of the run</h4>" + table(
        ["Regime", "Season", "Fitted on", "τ", "Pairs", "FY4B wet", "Gauge wet", "Median: FY4B → gauge", "99%: FY4B → gauge"],
        [[esc(REGIME_NAME[r["regime"]]), r["season"], {"cell": "own cell", "regime": "regime", "pooled": "all hours"}[r["pool"]], f(r["tau"], 3),
          n(r["n_pairs"]), n(r["n_wet_sat"]), n(r["n_wet_obs"]), f'{f(r["q50_sat"], 2)} → {f(r["q50_obs"], 2)}',
          f'{f(r["q99_sat"], 1)} → {f(r["q99_obs"], 1)}'] for r in maps], num=(3, 4, 5, 6, 7, 8))
        + source("fy4b_qm_maps.csv of the run")))

    strengths = [r for r in rows("strengths.csv") if r["scheme"] == "balanced_spatial"]
    def level(r):
        return "no adjustment" if r["threshold_applied"] == "false" else ("threshold only" if float(r["strength"]) == 0 else f'threshold + {float(r["strength"]):g} of the map')
    level_counts: dict[str, int] = {}
    for r in strengths:
        level_counts[level(r)] = level_counts.get(level(r), 0) + 1
    d.append(step(13, "How strongly the map is applied, month by month", where(("src/models/SatelliteCalibration.jl", "choose_strength"),
                                                                                ("src/models/SatelliteCalibration.jl", "_adjusted_value"),
                                                                                ("src/models/SatelliteCalibration.jl", "QM_STRENGTHS")),
                  checked("strength of the map"), md(r"""
With `--fy4b-qm-adaptive`, six candidates are defined from the map $m$ and threshold $\tau$ of the hour's map cell. Candidate 0 leaves the value alone; candidates $c = 1,\dots,5$ use strengths $s_c \in \{0, \tfrac14, \tfrac12, \tfrac34, 1\}$:

$$a_0(x) = x, \qquad a_c(x) = \begin{cases} 0 & x \le \tau \\ (1 - s_c)\,x + s_c\,m(x) & x > \tau \end{cases}$$

So $c = 1$ is the threshold alone and $c = 5$ the full map. In each month cell $M$ the candidate with the smallest squared error on the training gauges is used:

$$c^\star_M = \operatorname*{arg\,min}_{c \in \{0,\dots,5\}}\ \sum_{i \in \mathcal T_f}\ \sum_{t \in M}\big(a_c(x_{i,t}) - y_{i,t}\big)^2,$$

the smaller $c$ winning a tie. The corrected value is $\tilde F_{i,t} = a_{c^\star_{M(t)}}\big(F^{(\delta^\star)}_{i,t}\big)$ at every gauge. Without the flag every hour gets the full map.
""") + "<h4>Choices in the run, over all folds and month cells</h4>"
        + table(["Applied", "Fold-months"], [[esc(k), n(v)] for k, v in sorted(level_counts.items(), key=lambda kv: -kv[1])], num=(1,))
        + source("fy4b_qm_strengths.csv of the run")))
    out.append(part("qm", "D", "Corrected FY4B, inside each fold",
                    "A second anchor made from FY4B and the training gauges: a pixel shift, a wet threshold, a quantile map, and a month-by-month strength.", "".join(d)))

    # === E. Methods ==========================================================================================
    scan = [r for r in run("parameter_scan.csv") if r["scheme"] == "balanced_spatial" and r["product"] in PRODUCT_LABEL and r["selected"] == "true"]
    def chosen(product, mode, method, fold):
        return [r for r in scan if r["product"] == product and r["mode"] == mode and r["method"] == method and r["fold"] == fold]
    folds = sorted({r["fold"] for r in scan})
    e: list[str] = []

    e.append(step(14, "The residual scheme: how an anchor is corrected", where(("src/benchmark/InterpolationBenchmarkPredictors.jl", "predict_selected"),
                                                                                ("src/benchmark/InterpolationBenchmarkConfig.jl", "BENCHMARK_RUNS")), checked(), md(r"""
Seven runs are made for every anchor and fold. Four interpolate the gauges directly and ignore the anchor (`idw`, `adw`, `tps`, `gwr`). Three correct the anchor (`residual_gwr`, `mixed_gwr`, `mgwr`): they interpolate the **residual** at the training gauges,

$$\rho_{i,t} = G_{i,t} - S_{i,t}, \qquad i \in \mathcal T_f,$$

hour by hour, to the held-out gauge, and add it back to the anchor there:

$$\hat y_{0,t} = \max\big(S_{0,t} + \gamma\,\hat\rho_{0,t},\ 0\big) \quad \text{(residual methods)}, \qquad \hat y_{0,t} = \max\big(\hat G_{0,t},\ 0\big) \quad \text{(direct methods)}.$$

$\gamma \in \{0.1, 0.2, \dots, 1\}$ is a shrinkage scale chosen with the bandwidth (step 22). The raw anchor $S$ itself is scored as method `raw`. Every prediction is clipped at zero. The anchor is $F$ (FY4B), $\tilde F$ (FY4B_QM), GPM, GSMaP or a merge (step 24).

This is the bookkeeping around the methods; each method below has its own test.
""")))

    adw = [chosen("FY4B", "direct", m, fold)[0] for m in ("idw", "adw") for fold in folds]
    e.append(step(15, "Gauge-only baselines: IDW and ADW", where(("src/models/TraditionalInterpolation.jl", "_weighted_predict"),
                                                                   ("src/models/TraditionalInterpolation.jl", "_selection_weights")),
                  checked("IDW and ADW"), md(r"""
For a target point and an hour, let $\mathcal N$ be the training gauges that report in that hour (all of them, or the $K$ nearest reporting ones), at distances $d_j$. Inverse-distance weighting with power $p$:

$$w_j = d_j^{-p}, \qquad \hat G_0 = \frac{\sum_{j \in \mathcal N} w_j\,G_j}{\sum_{j \in \mathcal N} w_j}.$$

Angular-distance weighting multiplies each weight by Shepard's directional factor. With $\theta_{jl}$ the angle at the target between gauges $j$ and $l$ (from unit vectors in local km coordinates centred on the target):

$$w'_j = w_j\left[1 + \frac{\sum_{l \ne j} w_l\,(1 - \cos\theta_{jl})}{\sum_{l \ne j} w_l}\right], \qquad \hat G_0 = \frac{\sum_j w'_j\,G_j}{\sum_j w'_j}.$$

A gauge whose neighbours lie in the same direction is down-weighted; an isolated one is up-weighted. If a gauge lies on the target ($d_j \le \sqrt{\varepsilon}$) its value is returned. The neighbour set and the factor are formed from the gauges that actually report in the hour, not from the full network.

Tuned per fold on the inner split: $p \in \{1, 1.5, 2, 2.5, 3\}$ and $K \in \{8, 16, 32, 64, \text{all}\}$.
""") + "<h4>Chosen in the run (these two methods do not use the anchor)</h4>" + table(
        ["Method", "Fold", "Power p", "Neighbours K"], [[esc(r["method"].upper()), esc(r["fold"]), f(r["power"], 1), "all" if r["neighbors"] == "0" else n(r["neighbors"])] for r in adw], num=(2,))
        + source("parameter_scan.csv of the run")))

    tps = [chosen("FY4B", "direct", "tps", fold)[0] for fold in folds]
    e.append(step(16, "Thin-plate spline", where(("src/models/TraditionalInterpolation.jl", "tps_predict"),
                                                   ("src/models/TraditionalInterpolation.jl", "_tps_system"),
                                                   ("src/models/TraditionalInterpolation.jl", "_tps_scale")), checked("thin-plate spline"), md(r"""
In local km coordinates $\mathbf x_j = (X_j, Y_j)$, with the kernel

$$\phi(r) = \tfrac12\,r^2\ln r^2 \quad (\phi(0) = 0),$$

the spline through the $n$ reporting gauges of the hour solves

$$\begin{bmatrix} K + \lambda\,\sigma\,I & P \\ P^\top & 0 \end{bmatrix}\begin{bmatrix} \mathbf a \\ \mathbf b \end{bmatrix} = \begin{bmatrix} \mathbf G \\ \mathbf 0 \end{bmatrix}, \qquad K_{jl} = \phi(\lVert \mathbf x_j - \mathbf x_l\rVert), \quad P = [\,\mathbf 1\ \ \mathbf X\ \ \mathbf Y\,],$$

and predicts

$$\hat G_0 = \sum_{j=1}^{n} a_j\,\phi(\lVert \mathbf x_0 - \mathbf x_j\rVert) + b_1 + b_2 X_0 + b_3 Y_0.$$

$\lambda$ is a dimensionless smoothing parameter and $\sigma$ the median of the non-zero $|K_{jl}|$ over the gauges that report at least once in the call, so the smoothing means the same thing in every hour. At least three non-collinear gauges are needed. Tuned over $\lambda \in \{10^{-4}, 10^{-3}, 10^{-2}, 10^{-1}, 1\}$.
""") + table(["Fold", "λ chosen in the run"], [[esc(r["fold"]), f'{float(r["smooth"]):g}'] for r in tps], num=(1,))))

    e.append(step(17, "Kernels and bandwidths", where(("src/core/kernel.jl", "gw_weight_gaussian"), ("src/core/gw_weight.jl", "gw_weight!")),
                  checked("kernels and the adaptive bandwidth"), md(r"""
Every GWR-family method weights the training gauges by distance $d$ with one of five kernels and a bandwidth $b$:

$$\begin{aligned} \text{Gaussian:}\quad & w = \exp\!\big(-d^2 / 2b^2\big) & \text{exponential:}\quad & w = \exp(-d/b) \\ \text{bisquare:}\quad & w = \big(1 - (d/b)^2\big)^2\ \ (d \le b) & \text{tricube:}\quad & w = \big(1 - (d/b)^3\big)^3\ \ (d \le b) \\ \text{boxcar:}\quad & w = 1\ \ (d \le b) \end{aligned}$$

with $w = 0$ beyond $b$ for the last three. The bandwidth comes from one of two families:

- **fixed:** $b$ in km, the same at every target;
- **adaptive:** $b$ is the distance from the target to its $k$-th nearest reporting gauge. If only $m < k$ gauges report, $b = (k/m)\,\max_j d_j$.

$b = \infty$ gives weight 1 everywhere, which is an ordinary (global) least-squares fit; it is offered as a candidate of the fixed family.
""")))

    e.append(step(18, "The local fit, and direct GWR", where(("src/core/local_hat.jl", "local_hat_matrix"),
                                                               ("src/benchmark/InterpolationBenchmarkPredictors.jl", "_gwr_predict"),
                                                               ("src/mger/pipeline/mger_residual_predict.jl", "st_gwr_predict_nanaware")),
                  checked("local fit (hat matrix)"), md(r"""
Given a design matrix $X$ (one row per training gauge), a response $\mathbf y$ and kernel weights $W_0 = \operatorname{diag}(w_1,\dots,w_n)$ for a target with design row $\mathbf x_0$, the local fit and its prediction are

$$\hat{\boldsymbol\beta}_0 = \big(X^\top W_0 X + \epsilon I\big)^{-1} X^\top W_0\,\mathbf y, \qquad \hat y_0 = \mathbf x_0^\top\hat{\boldsymbol\beta}_0 = \mathbf h_0^\top\mathbf y, \qquad \mathbf h_0^\top = \mathbf x_0^\top\big(X^\top W_0 X + \epsilon I\big)^{-1} X^\top W_0,$$

with a small ridge $\epsilon = 10^{-8}$. $\mathbf h_0$ is the target's row of the **hat matrix** $H$; stacking the rows of the training gauges themselves gives the $n \times n$ matrix that maps a response to its local fit. A target with fewer than $p + 1$ positively weighted gauges ($p$ columns of $X$) is not fitted and its prediction is missing.

**Direct GWR** (`gwr`) uses the response $\mathbf y = \mathbf G_t$ and the design $[\,1,\ \lambda_j - \bar\lambda,\ \varphi_j - \bar\varphi\,]$ (longitude and latitude in degrees, centred on the training gauges), needing at least 3 reporting gauges with weight.

The test covers the hat row and the too-few-gauges rule; the direct-GWR wrapper is read from the code.
""")))

    e.append(step(19, "Covariates: preparation", where(("src/models/JointCovariateModels.jl", "build_joint_fold_context"),
                                                         ("src/models/JointCovariateModels.jl", "_hierarchical_weights"),
                                                         ("src/models/JointCovariateModels.jl", "_spatial_design"),
                                                         ("src/models/CovariateGroups.jl", "JOINT_GROUP_ORDER")),
                  checked("covariate centring and scaling"), md(r"""
The three residual methods may use covariates. Nine groups are candidates: elevation, slope, aspect (sine and cosine) from the DEM; 2 m temperature, 2 m dew point, 10 m wind $u$ and $v$, surface pressure from ERA5-Land, hourly; and NDVI, 16-day.

**Spatial block.** Longitude and latitude are standardised over the training gauges, $z^\lambda_j = (\lambda_j - \bar\lambda)/s_\lambda$ and $z^\varphi_j = (\varphi_j - \bar\varphi)/s_\varphi$, giving the columns $[\,1,\ z^\lambda,\ z^\varphi\,]$.

**A covariate** $x_{j,t}$ is first centred hour by hour on the training gauges, then standardised once:

$$x'_{j,t} = x_{j,t} - \frac{1}{|\mathcal T_f|}\sum_{l \in \mathcal T_f} x_{l,t}, \qquad \tilde x_{j,t} = \frac{x'_{j,t} - \mu}{\sigma},$$

$$\mu = \frac{\sum_{j,t} \omega_{j,t}\,x'_{j,t}}{\sum_{j,t}\omega_{j,t}}, \qquad \sigma^2 = \frac{\sum_{j,t}\omega_{j,t}\,(x'_{j,t} - \mu)^2}{\sum_{j,t}\omega_{j,t}}.$$

The weights $\omega_{j,t}$ are non-zero only on the training gauges' wet hours ($G \ge 0.1$, anchor and covariates present) and give every gauge the same total: $\omega_{j,t} = 1/n_j$ with $n_j$ the gauge's number of such hours. (When NDVI is in the model the gauge's weight is first split equally over its NDVI periods.) Held-out gauges are transformed with the training fold's hour means, $\mu$ and $\sigma$.
""")))

    roles = [r for r in run("joint_fold_roles.csv") if r["scheme"] == "balanced_spatial" and r["product"] in PRODUCT_LABEL and r["final_included"] == "true"]
    role_rows = [[esc(PRODUCT_LABEL[p]), esc(fold), esc(", ".join(f'{r["variable_group"]} ({r["role"]})' for r in roles if r["product"] == p and r["fold"] == fold) or "none")]
                 for p in PRODUCT_LABEL for fold in folds]
    e.append(step(20, "Covariates: selection inside the fold", where(("src/selection/JointVariableSelection.jl", "select_joint_covariates"),
                                                                       ("src/selection/JointVariableSelection.jl", "joint_vif"),
                                                                       ("src/selection/JointVariableSelection.jl", "joint_spatial_variability_test"),
                                                                       ("src/models/dem/screening.jl", "terrain_screen"),
                                                                       ("src/selection/era5/era5_panel_screen.jl", "dynamic_panel_screen")),
                  checked("VIF, permutation p, Benjamini-Hochberg, joint F"), md(r"""
With `--nested-covariates` the covariates are chosen again in every fold, per anchor, from the training gauges only. Three stages.

**1. Each group on its own.** The response is the wet-hour residual $G - S$. For a DEM group it is each gauge's mean wet residual (gauges with at least 100 wet hours); the statistic is $|r|$ (Pearson) for elevation and slope, and for aspect the joint $F$ of its two columns:

$$F = \frac{(\mathrm{RSS}_0 - \mathrm{RSS}_1)/p}{\mathrm{RSS}_1/(n - p - 1)}.$$

Its p-value comes from 999 permutations of the response, and is adjusted within the family by Benjamini–Hochberg:

$$p = \frac{1 + \#\{\text{permuted statistic} \ge \text{observed}\}}{1 + 999}, \qquad q_{(i)} = \min_{j \ge i}\ \min\!\Big(1,\ \frac{m\,p_{(j)}}{j}\Big).$$

A group passes with $q < 0.05$ and a stable direction (Pearson and Spearman of the same sign). ERA5 and NDVI groups are screened the same way on the station-hour panel.

**2. Collinearity across groups.** For each column $j$ of the surviving groups, with the weights $\omega$ of step 19,

$$\mathrm{VIF}_j = \frac{1}{1 - R_j^2},$$

$R_j^2$ from the weighted regression of column $j$ on the others. While the largest group VIF is 5 or more, the group that has it is dropped (ties: the one with the larger $q$) and the VIFs are recomputed.

**3. Local or global.** A GWR of the residual on $[\,1, z^\lambda, z^\varphi,$ covariates$\,]$ is fitted at every training gauge (bisquare kernel, $k$ nearest gauges with $k$ chosen from $\{30, 50, 80, 120, 160\}$ by leave-one-out error). For each group the spread of its coefficient over the gauges,

$$T = \sum_{\text{columns}} \operatorname{Var}_j\big(\hat\beta_j\big),$$

is compared with 999 permutations of the gauges' positions. A group with $q < 0.05$ varies in space and is **local**; otherwise it is **global**.

The four statistics are tested. The ERA5 and NDVI panel screens (their weighting and the block permutation of stations) are described from the code and not rederived here.
""") + "<h4>Covariates kept in the run, with their role</h4>" + table(["Anchor", "Fold", "Groups"], role_rows, cls="wrap")
        + source("joint_fold_roles.csv of the run")))

    e.append(step(21, "Residual GWR, mixed GWR and MGWR", where(("src/models/JointCovariateModels.jl", "_design_at"),
                                                                  ("src/models/JointCovariateModels.jl", "dynamic_covariate_predict"),
                                                                  ("src/core/mixed_gwr.jl", "mixed_gwr_predict"),
                                                                  ("src/core/mixed_gwr.jl", "_backfit_components"),
                                                                  ("src/models/JointCovariateModels.jl", "_multiscale_predict_damped")),
                  checked("mixed GWR back-fit", "MGWR back-fit with over-relaxation"), md(r"""
All three fit the residual $\boldsymbol\rho_t$ of one hour over the training gauges that have it. They differ in the design.

**Residual GWR.** One local fit (step 18) with the design $[\,1, z^\lambda, z^\varphi,$ every selected covariate$\,]$ and one bandwidth: $\hat\rho_0 = \mathbf h_0^\top\boldsymbol\rho$.

**Mixed GWR.** Local covariates go in the local design $X_L = [\,1, z^\lambda, z^\varphi,$ local covariates$\,]$; global covariates $X_G$ get one coefficient for the whole network:

$$\rho_j = \mathbf x_{L,j}^\top\boldsymbol\beta(\mathbf u_j) + \mathbf x_{G,j}^\top\boldsymbol\gamma + \varepsilon_j.$$

With the local hat matrix $H$ and the global projection $P_G = X_G(X_G^\top X_G + \epsilon I)^{-1}X_G^\top$, the two parts are fitted in turn, starting from $\mathbf g = P_G\boldsymbol\rho$:

$$\mathbf l \leftarrow H(\boldsymbol\rho - \mathbf g), \qquad \mathbf g \leftarrow P_G(\boldsymbol\rho - \mathbf l),$$

until the residual sum of squares changes by less than $10^{-5}$ in relative terms. Then

$$\hat{\boldsymbol\gamma} = (X_G^\top X_G + \epsilon I)^{-1}X_G^\top(\boldsymbol\rho - \mathbf l), \qquad \hat\rho_0 = \mathbf h_0^\top(\boldsymbol\rho - X_G\hat{\boldsymbol\gamma}) + \mathbf x_{G,0}^\top\hat{\boldsymbol\gamma}.$$

When no covariate is global the two models are the same fit; the run records that.

**MGWR.** Each local group gets its own bandwidth: a locally varying intercept (group 0), and one group per local covariate. The coordinate columns are dropped (`mgwr_spatial_grouping = intercept_only`). With group hat matrices $H_j$ at bandwidths $b_j$ and components $\mathbf f_j$ starting at zero, one sweep is

$$\mathbf f_j \leftarrow \omega\,H_j\Big(\boldsymbol\rho - \mathbf g - \sum_{l \ne j}\mathbf f_l\Big) + (1 - \omega)\,\mathbf f_j \quad (j = 0, 1, \dots), \qquad \mathbf g \leftarrow \omega\,P_G\Big(\boldsymbol\rho - \sum_l \mathbf f_l\Big) + (1-\omega)\,\mathbf g,$$

with over-relaxation $\omega = 1.5$. Sweeps stop when the fitted vector $\hat{\boldsymbol\rho} = \sum_l\mathbf f_l + \mathbf g$ moves by less than $10^{-5}$ in relative norm, $\lVert\hat{\boldsymbol\rho}^{(s)} - \hat{\boldsymbol\rho}^{(s-1)}\rVert / \lVert\hat{\boldsymbol\rho}^{(s-1)}\rVert$, within 1000 sweeps; an hour that does not converge has no prediction. At a target,

$$\hat\rho_0 = \mathbf x_{G,0}^\top\hat{\boldsymbol\gamma} + \sum_j \mathbf h_{j,0}^\top\Big(\boldsymbol\rho - X_G\hat{\boldsymbol\gamma} - \sum_{l \ne j}\mathbf f_l\Big).$$

An hour is fitted only when the number of usable training gauges exceeds the number of columns plus 2. An adaptive bandwidth is capped at the number of usable gauges minus 1.
""")))

    def chosen_text(product, method, fold):
        got = chosen(product, "residual", method, fold)
        if not got:
            return "—"
        if method != "mgwr":
            r = got[0]
            return f'{KERNELS[int(r["kernel"])]}, {bw_text(r)}, γ {float(r["shrink"]):g}'
        parts = [f'{r["group"]}: {bw_text(r)}' for r in got]
        return f'{KERNELS[int(got[0]["kernel"])]}; ' + "; ".join(parts) + f'; γ {float(got[0]["shrink"]):g}'
    tune_rows = [[esc(PRODUCT_LABEL[p]), esc({"gwr": "residual GWR", "mixed_gwr": "mixed GWR", "mgwr": "MGWR"}[m]), esc(fold), esc(chosen_text(p, m, fold))]
                 for p in PRODUCT_LABEL for m in ("gwr", "mixed_gwr", "mgwr") for fold in folds]
    e.append(step(22, "Tuning: hours, score, grids", where(("src/benchmark/InterpolationBenchmarkTuning.jl", "_tuning_time_sample"),
                                                             ("src/benchmark/InterpolationBenchmarkTuning.jl", "_candidate_metrics"),
                                                             ("src/benchmark/InterpolationBenchmarkTuning.jl", "_select_candidate!"),
                                                             ("src/benchmark/InterpolationBenchmarkJoint.jl", "_select_joint_coordinate_descent!"),
                                                             ("src/benchmark/InterpolationBenchmarkJoint.jl", "_joint_candidate_metrics")),
                  checked("tuning hours"), md(r"""
Every parameter is chosen per fold, anchor and method, on the training gauges only.

**Hours.** At most 336 hours are used: the 168 hours with the largest mean gauge rain over the training gauges, joined with 168 hours spaced evenly through the record,

$$\mathcal H = \text{top}_{168}\big(\bar G_t\big)\ \cup\ \Big\{\operatorname{round}\Big(1 + \tfrac{(T-1)(m-1)}{167}\Big) : m = 1,\dots,168\Big\}.$$

**Score.** A candidate predicts each inner group of training gauges from the other four (step 8). Over the cells with a gauge value, an anchor value and a prediction, its score is

$$\mathrm{RMSE} = \sqrt{\frac{1}{N}\sum_{(i,t)}\big(\hat y_{i,t} - G_{i,t}\big)^2}.$$

A candidate must predict at least 95% of the scorable cells. The smallest RMSE wins, then the smallest MAE, then the larger coverage.

**Grids.** Kernels: all five. Bandwidths: adaptive $k \in \{30, 50, 80, 120\}$ nearest gauges; fixed $b \in \{10, 20, 30, 50\}$ km and $\infty$ (global). Shrinkage $\gamma \in \{0.1, \dots, 1\}$ is rescored for every candidate at no extra fit, on the clipped prediction $\max(S + \gamma\hat\rho, 0)$.

**MGWR** has one bandwidth per group. For each kernel and bandwidth family, all groups start at the widest candidate; then, group by group, each candidate is tried with the other groups held fixed and the best is kept; sweeps over the groups repeat until no bandwidth changes (at most 5). The kernel and family with the lowest final RMSE win.

Only the choice of hours has a test here; the search loops are described from the code.
""") + "<h4>Chosen in the run for the two FY4B anchors</h4>" + table(["Anchor", "Method", "Fold", "Kernel, bandwidth, shrinkage"], tune_rows, cls="wrap")
        + source("parameter_scan.csv of the run")))
    out.append(part("methods", "E", "The methods FY4B is run through",
                    "The gauge-only baselines, the local regression machinery, and the three methods that correct an anchor with it.", "".join(e)))

    # === F. Merges and blends ===================================================================================
    g_parts: list[str] = []
    merge = [r for r in rows("fusion_merge_weights.csv") if r["scheme"] == "balanced_spatial" and r["product"] == "MERGED_OLS"]
    g_parts.append(step(23, "Merged anchors", where(("src/models/SatelliteFusion.jl", "fit_satellite_fusion"),
                                                      ("src/models/SatelliteFusion.jl", "apply_satellite_fusion"),
                                                      ("src/models/SatelliteFusion.jl", "equal_weight_coefficients")), checked("merged anchors"), md(r"""
With `--fused-anchor`, two more anchors are built from the three products $F$ (plain FY4B), $P$ (GPM) and $Q$ (GSMaP):

$$S^{\text{mean}}_{i,t} = \max\Big(\tfrac13\big(F_{i,t} + P_{i,t} + Q_{i,t}\big),\ 0\Big), \qquad S^{\text{ols}}_{i,t} = \max\big(\beta_0 + \beta_1 F_{i,t} + \beta_2 P_{i,t} + \beta_3 Q_{i,t},\ 0\big).$$

The OLS coefficients minimise the squared error against the gauges over the fold's training gauge-hours at which the gauge and all three products have a value,

$$\hat{\boldsymbol\beta} = \operatorname*{arg\,min}_{\boldsymbol\beta}\sum_{i \in \mathcal T_f}\sum_t\big(G_{i,t} - \beta_0 - \beta_1 F_{i,t} - \beta_2 P_{i,t} - \beta_3 Q_{i,t}\big)^2,$$

solved from the normal equations. A cell where any product is missing has no merged value. Plain FY4B is the FY4B of the merges; corrected FY4B is not used here.
""") + "<h4>OLS coefficients in the run</h4>" + table(["Fold", "β₀", "β₁ (FY4B)", "β₂ (GPM)", "β₃ (GSMaP)"],
        [[esc(r["fold"]), f(r["intercept"], 3), f(r["beta_fy4b"], 3), f(r["beta_gpm_uncal"], 3), f(r["beta_gsmap"], 3)] for r in merge], num=(1, 2, 3, 4))
        + source("fused_anchor_selection.csv of the run")))

    blend = [r for r in rows("fusion_blend_weights.csv") if r["scheme"] == "balanced_spatial"]
    band_head = [f"{a} wet, max {e}" for a in ("1", "2", "3") for e in ("< 2.5", "2.5–8", "≥ 8")]
    g_parts.append(step(24, "Blends with ADW, and the agreement bands", where(("src/benchmark/InterpolationBenchmarkPredictors.jl", "satellite_wet_blend_prediction"),
                                                                              ("src/benchmark/InterpolationBenchmarkPredictors.jl", "blend_band_matrix"),
                                                                              ("src/benchmark/InterpolationBenchmarkPredictors.jl", "banded_blend_prediction"),
                                                                              ("src/benchmark/InterpolationBenchmarkPredictors.jl", "select_blend_lambdas"),
                                                                              ("src/benchmark/InterpolationBenchmarkConfig.jl", "BLEND_ENVELOPE_EDGES")),
                        checked("blends and agreement bands"), md(r"""
A residual method passes a satellite false alarm straight into its prediction. With `--satellite-wet-blend`, each of the three residual methods also gets a blended version that leans toward gauge-only ADW where the satellites report rain.

**Constant blend** (`blend_mgwr`, …). One weight $\lambda$ for every cell the anchor calls wet:

$$\hat y^{\,b}_{i,t} = \begin{cases} (1-\lambda)\,\hat y_{i,t} + \lambda\,\hat y^{\text{ADW}}_{i,t} & S_{i,t} \ge 0.1 \\ \hat y_{i,t} & \text{otherwise.} \end{cases}$$

**Agreement blend** (`blend_agrenv_mgwr`, …; `--blend-axis agreement_envelope`). The weight depends on what the three products $F, P, Q$ say about the cell. Let

$$a = \min\big(\#\{F, P, Q \ge 0.1\},\ 3\big), \qquad e = \max(F, P, Q), \qquad \epsilon = \begin{cases} 1 & e < 2.5 \\ 2 & 2.5 \le e < 8 \\ 3 & e \ge 8 \end{cases}$$

(a missing product counts as not wet). The cell's band is

$$B = \begin{cases} 0 & a = 0 \\ 3\,(a - 1) + \epsilon & a \ge 1 \end{cases}$$

and the blend is $\hat y^{\,b} = (1 - \lambda_B)\,\hat y + \lambda_B\,\hat y^{\text{ADW}}$ in bands 1–9, with band 0 left as it is. A cell where ADW has no value keeps $\hat y$. Note that a cell the anchor itself calls dry is blended when another product calls it wet.

**The weights.** On the inner split's predictions at the tuning hours, each band's weight is the grid value with the smallest squared error in that band,

$$\lambda_B = \operatorname*{arg\,min}_{\lambda \in \{0, 0.1, \dots, 1\}}\ \sum_{(i,t) \in B}\Big((1-\lambda)\,\hat y_{i,t} + \lambda\,\hat y^{\text{ADW}}_{i,t} - G_{i,t}\Big)^2,$$

the smaller $\lambda$ winning a tie. The bands separate because a weight only touches its own band's cells.

The bands always use plain FY4B, also when the anchor being corrected is FY4B_QM.
""") + "<h4>MGWR's band weights in the run, where FY4B or corrected FY4B is the anchor (weight on ADW)</h4>" + table(
        ["Anchor", "Fold"] + band_head, [[esc(PRODUCT_LABEL[r["product"]]), esc(r["fold"])] + [f'{float(x):g}' for x in r["lambdas"].split("|")] for r in blend],
        num=tuple(range(2, 11))) + source("blend_selection.csv of the run")))
    out.append(part("fusion", "F", "Merged anchors and blends", "Where plain FY4B meets the other two products: as a third of a merge, and as one of three votes on how much to trust the correction.",
                    "".join(g_parts)))

    # === G. Scores ==================================================================================================
    s_parts: list[str] = []
    s_parts.append(step(25, "The evaluation mask and the scores", where(("src/core/metrics.jl", "metric_continuous"), ("src/core/metrics.jl", "metric_event"),
                                                                        ("src/benchmark/InterpolationBenchmarkMetrics.jl", "_common_method_mask"),
                                                                        ("src/benchmark/InterpolationBenchmarkConfig.jl", "MASK_METHODS")), checked("scores"), md(r"""
Out-of-fold predictions of the five folds are stitched into one stations × hours matrix per anchor and method. For one anchor, a cell is scored when the gauge and all eight mask methods have a value there:

$$M_{i,t} = \mathbb 1[G_{i,t}\ \text{present}]\ \prod_{m}\mathbb 1[\hat y^{(m)}_{i,t}\ \text{present}], \qquad m \in \{\text{raw, idw, adw, tps, gwr, residual\_gwr, mixed\_gwr, mgwr}\}.$$

Blends and reference predictors are scored on this mask but do not define it. With $e_{i,t} = \hat y_{i,t} - G_{i,t}$ and $N = \sum M$:

$$\mathrm{RMSE} = \sqrt{\tfrac1N\textstyle\sum e^2}, \qquad \mathrm{MAE} = \tfrac1N\textstyle\sum |e|, \qquad \mathrm{Bias} = \tfrac1N\textstyle\sum e, \qquad r = \operatorname{corr}(G, \hat y).$$

For a threshold $\theta$ (0.1, 2, 8 mm h⁻¹ and others), with hits $h$ ($G \ge \theta$ and $\hat y \ge \theta$), misses $m$ and false alarms $f$:

$$\mathrm{POD} = \frac{h}{h+m}, \qquad \mathrm{FAR} = \frac{f}{h+f}, \qquad \mathrm{CSI} = \frac{h}{h+m+f}.$$

The same scores are also given by stratum: gauge intensity classes with bounds 0.1, 2, 4, 8, 20 mm h⁻¹, and the distance to the nearest training gauge.
""")))
    s_parts.append(step(26, "The paired day-block bootstrap and Holm's correction", where(("src/benchmark/InterpolationBenchmarkBootstrap.jl", "_daily_bootstrap_delta"),
                                                                                         ("src/benchmark/InterpolationBenchmarkBootstrap.jl", "paired_bootstrap_rows"),
                                                                                         ("src/benchmark/InterpolationBenchmarkBootstrap.jl", "_holm_adjust")),
                        checked("paired day-block bootstrap and Holm"), md(r"""
A method is compared with each gauge-only baseline (IDW, ADW, TPS) on the same cells. For every calendar day $D$ of the hour labels, the day's cell count and squared errors are summed once:

$$n_D = \sum_{t \in D}\sum_i M_{i,t}, \qquad \mathrm{SSE}^{\text{base}}_D,\ \ \mathrm{SSE}^{\text{meth}}_D.$$

One bootstrap replicate draws, with replacement, as many days as there are days with $n_D > 0$, and recomputes both scores from the drawn days:

$$\Delta^\ast = \sqrt{\frac{\sum_{D^\ast}\mathrm{SSE}^{\text{base}}_D}{\sum_{D^\ast} n_D}} - \sqrt{\frac{\sum_{D^\ast}\mathrm{SSE}^{\text{meth}}_D}{\sum_{D^\ast} n_D}}.$$

Both methods see the same drawn days, so the comparison stays paired, and a day's hours and gauges stay together. Over the replicates, the 95% interval is the 2.5% and 97.5% quantiles of $\Delta^\ast$, and the two-sided p-value is

$$p = \min\!\Big(1,\ 2\min\big(\Pr[\Delta^\ast \le 0],\ \Pr[\Delta^\ast \ge 0]\big)\Big).$$

Within a stratum the p-values of the baselines are adjusted by Holm's step-down rule. With $p_{(1)} \le \dots \le p_{(m)}$:

$$\tilde p_{(i)} = \max_{j \le i}\ \min\big(1,\ (m - j + 1)\,p_{(j)}\big).$$

The resampling covers time only: the gauges and the folds stay fixed.
""")))
    out.append(part("scores", "G", "Scores and tests", "What the benchmark writes about each anchor and method: pooled scores on a shared mask, and paired intervals against the gauge-only baselines.",
                    "".join(s_parts)))

    # === H. Other switches ==============================================================================================
    h = step(27, "Other FY4B switches in the code", where(("src/models/SatelliteFusion.jl", "confirmed_by_others"),
                                                            ("src/models/SatelliteCalibration.jl", "calibrate_crossed"),
                                                            ("src/models/HumidityGate.jl", "HumidityGate")), checked(), md(r"""
These options exist in the benchmark and are off in the run shown above. Each changes which FY4B a later step receives; none changes a formula above.

- `--fy4b-qm-fusion`: the merges (step 23) and the agreement bands (step 24) take each fold's corrected FY4B, $\tilde F^{(f)}$, in place of plain $F$.
- `--fy4b-out-of-fusion`: FY4B stays in the product list but the merges and the bands use GPM and GSMaP only.
- `--fy4b-confirm-fusion`: the merges and the bands take a FY4B kept only where another product also reports rain, $F^{c}_{i,t} = F_{i,t}\cdot\mathbb 1\big[\max(P_{i,t}, Q_{i,t}) \ge 0.1\big]$ (a missing other product counts as dry; a missing $F$ stays missing). It uses no gauge.
- `--fy4b-qm-crossed`: steps 11–13 are also fitted out of time. The record is cut into 7-day blocks dealt into time folds; the hours of one time fold take the shift, map and strength fitted on the training gauges in the other time folds, leaving out one day on either side.
- `--fy4b-qm-rh-gate`: adds a further anchor, corrected FY4B set to 0 where ERA5-Land near-surface humidity is below a cut that is fitted across gauge folds and time folds.

Described from the code's own documentation; none of these has a test in the formula file.
"""))
    out.append(part("switches", "H", "Other FY4B switches", "Options the code offers for FY4B that the run shown here does not use.", h))

    # --- Index of code ---------------------------------------------------------------------------------------------------
    seen: dict[str, list[str]] = {}
    for path, name in _cited:
        seen.setdefault(path, [])
        name in seen[path] or seen[path].append(name)
    out.append(part("index", "·", "Where each step lives in the code", "Every function named on this page, by file. The page is rebuilt only if all of them still exist.",
                    table(["File", "Functions and constants"], [[f"<code>{esc(path)}</code>", ", ".join(f"<code>{esc(x)}</code>" for x in names)]
                                                               for path, names in sorted(seen.items())], cls="wrap")))

    nav = "".join(f'<a href="#{anchor}">{label}</a>' for anchor, label in (
        ("notation", "Notation"), ("raw", "A Raw file → hourly value"), ("load", "B What the benchmark loads"), ("cv", "C Cross-validation"),
        ("qm", "D Corrected FY4B"), ("methods", "E Methods"), ("fusion", "F Merges and blends"), ("scores", "G Scores"),
        ("switches", "H Other switches"), ("index", "Code index")))
    head = PAGE_HEAD.replace("<title>FY4B process trace</title>", "<title>FY4B methods</title>").replace("</style>", EXTRA_CSS + "</style>")
    return head + out[0] + f'<nav class="toc">{nav}</nav>' + "".join(out[1:]) + PAGE_TAIL + MATH_SCRIPT


EXTRA_CSS = """
.step { padding: 18px 0 20px; border-bottom: 1px solid var(--line) } .step:last-child { border-bottom: 0 }
.step h3 { font-size: 17px; margin: 0 0 6px; display: flex; gap: 12px; align-items: baseline }
.step h4 { font-size: 13.5px; margin: 18px 0 4px; color: var(--fg2) }
.stepmeta { margin: 0 0 8px } .pill.warn { background: var(--empty); color: var(--fg2); font-weight: 500 }
.where { font-size: 12.5px; color: var(--muted); margin-top: 12px; max-width: none }
.step ul { margin: 8px 0; padding-left: 20px; max-width: 76ch } .step li { margin: 4px 0 }
div.math { padding: 6px 0; font-size: 1.08em; white-space: nowrap; min-width: max-content }
.math { font-family: var(--mono) } .math math { font-family: "STIX Two Math", "Cambria Math", "Latin Modern Math", serif; font-size: 1.12em }
td .math { white-space: nowrap }
.meta a { color: var(--accent) }
"""

MATH_SCRIPT = f"""<script src="{KATEX}"></script>
<script>
(function () {{
  if (!window.katex) {{ return; }}
  document.querySelectorAll(".math").forEach(function (el) {{
    var tex = el.textContent;
    try {{
      katex.render(tex, el, {{ displayMode: el.tagName === "DIV", output: "mathml", throwOnError: false, strict: false }});
    }} catch (error) {{ el.textContent = tex; }}
  }});
}})();
</script>
"""


def main() -> None:
    page = build()
    write_plot_provenance(TRACE / "methods_page", ROOT, (TEST_FILE,))
    PAGE.write_text(SHELL.format(page), encoding="utf-8")
    FRAGMENT.write_text(page, encoding="utf-8")
    print(f"wrote {PAGE.relative_to(ROOT)} ({len(page) / 1024:.0f} KB), {len(set(_cited))} functions cited")


if __name__ == "__main__":
    main()

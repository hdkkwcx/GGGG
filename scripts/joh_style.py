"""Shared figure style for the plot_*.py scripts, sized to Journal of Hydrology (Elsevier) artwork rules.

Elsevier artwork guidance (checked 2026-10-13): widths 90 mm (single column), 140 mm (1.5 column) and
190 mm (full width); lettering 7 pt at the printed size, never below 6 pt; Arial/Helvetica, Times,
Courier or Symbol, embedded; vector EPS/PDF preferred, raster at 300 dpi (halftone), 500 dpi
(combination) or 1000 dpi (line art). A figure drawn at its printed width needs no scaling, so the
font sizes set here are the sizes that get printed.

    import joh_style
    from joh_style import INK, INK_SECONDARY, MUTED, GRID, AXIS, SURFACE
    fig = plt.figure(figsize=joh_style.mm(joh_style.DOUBLE, 120))
    joh_style.save(fig, FIG_DIR, "fig1", submission=True)
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

SINGLE = 90.0      # mm
ONE_HALF = 140.0
DOUBLE = 190.0
MAX_HEIGHT_MM = 240.0   # provisional: the journal's own guide was not reachable to confirm the page depth
MIN_FONT_PT = 6.0

SURFACE = "#ffffff"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"      # reference lines and outlines; spines and ticks use the darker MUTED to survive print

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.edgecolor": MUTED,
    "axes.linewidth": 0.6,
    "axes.labelcolor": INK_SECONDARY,
    "axes.titlesize": 7,
    "axes.titlecolor": INK,
    "text.color": INK,
    "xtick.color": MUTED,
    "ytick.color": MUTED,
    "xtick.labelcolor": INK_SECONDARY,
    "ytick.labelcolor": INK_SECONDARY,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "pdf.fonttype": 42,   # TrueType, embedded: no Type 3 fonts
    # Arial has no superscript minus (U+207B), so units are written "mm h$^{-1}$"; mathtext then sets them
    # in Arial too rather than in DejaVu.
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.default": "regular",
})


def mm(width: float, height: float) -> tuple[float, float]:
    """A figsize in inches from a width and height in millimetres."""
    return width / 25.4, height / 25.4


def check(fig: plt.Figure, name: str) -> list[str]:
    """Problems with `fig` at its printed size: too wide, too tall, or text below MIN_FONT_PT."""
    width, height = (v * 25.4 for v in fig.get_size_inches())
    problems = []
    if width > DOUBLE + 0.5:
        problems.append(f"width {width:.0f} mm > {DOUBLE:g} mm")
    if height > MAX_HEIGHT_MM + 0.5:
        problems.append(f"height {height:.0f} mm > {MAX_HEIGHT_MM:g} mm")
    small = sorted({round(t.get_fontsize(), 1) for t in fig.findobj(mpl.text.Text)
                    if t.get_visible() and t.get_text().strip() and t.get_fontsize() < MIN_FONT_PT})
    if small:
        problems.append(f"text at {', '.join(f'{s:g}' for s in small)} pt < {MIN_FONT_PT:g} pt")
    for problem in problems:
        print(f"  WARNING {name}: {problem}")
    return problems


def save(fig: plt.Figure, fig_dir: Path, name: str, submission: bool = False) -> None:
    """Vector PDF (rasterized layers at 600 dpi) and a 300 dpi PNG preview; `submission` adds a 600 dpi TIFF
    and checks the figure against the journal limits."""
    fig_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(fig_dir / f"{name}.png", dpi=300)
    fig.savefig(fig_dir / f"{name}.pdf", dpi=600)
    written = "png / .pdf"
    if submission:
        fig.savefig(fig_dir / f"{name}.tif", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
        written += " / .tif"
        check(fig, name)
    plt.close(fig)
    print(f"  wrote {name}.{written}")

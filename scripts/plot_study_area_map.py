"""Draw the study-area map (paper Fig. 1): relief, and the gauges by landform region.

Reads the 1 km SRTM DEM clipped to the study box (`data/dem_ShiYan_1km.tif`, written by
`data/clip_dem.R`), the gauge list (`data/processed/study_area/station_meta.csv`) and each gauge's
landform region (`data/processed/covariates/station_landform.csv`, from
`scripts/prepare_station_landform.jl`), and writes PNG (300 dpi) and PDF to
`output/study_area_map/figures/`:

    py -3.13 scripts/plot_study_area_map.py

The plotting environment is pinned in `scripts/requirements-plot.txt`; `plot_provenance.csv` beside the
figure records it and a SHA-256 for each of the three inputs, since no Julia run stands between them
and the figure. The GeoTIFF is read with Pillow, from its GeoKeys, so no GDAL is needed.

The location inset and the river basins the paper plans for Fig. 1 are not drawn: there is no boundary
or basin dataset in `data/` yet.

Gauges follow `plot_satellite_temporal_evaluation.py`'s landform figure: one marker shape per relief
region, and hollow markers for a region with fewer than 20 gauges, which is never interpreted.
"""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LightSource, LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from PIL import Image

from plot_provenance import write_plot_provenance

ROOT = Path(__file__).resolve().parents[1]
DEM = ROOT / "data" / "dem_ShiYan_1km.tif"
STATIONS = ROOT / "data" / "processed" / "study_area" / "station_meta.csv"
LANDFORM = ROOT / "data" / "processed" / "covariates" / "station_landform.csv"
FIG_DIR = ROOT / "output" / "study_area_map" / "figures"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
AXIS = "#c3c2b7"
REGIONS = ["relief_lt500", "relief_500_1000", "relief_ge1000"]
REGION_LABELS = {"relief_lt500": "Relief < 500 m", "relief_500_1000": "Relief 500–1000 m",
                 "relief_ge1000": "Relief ≥ 1000 m"}
REGION_MARKERS = {"relief_lt500": "o", "relief_500_1000": "s", "relief_ge1000": "^"}
# Elevation: a light-to-dark earth ramp, kept off blue so the relief never reads as water or rain.
ELEVATION = LinearSegmentedColormap.from_list(
    "elevation", ["#eef1e4", "#d6dfb8", "#e7d8a4", "#d2b07c", "#a8805c", "#6f5444"])
EARTH_RADIUS_KM = 6371.0

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "Arial", "DejaVu Sans"],
    "font.size": 7,
    "axes.edgecolor": AXIS,
    "axes.linewidth": 0.6,
    "axes.labelcolor": INK_SECONDARY,
    "text.color": INK,
    "xtick.color": AXIS,
    "ytick.color": AXIS,
    "xtick.labelcolor": INK_SECONDARY,
    "ytick.labelcolor": INK_SECONDARY,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "figure.facecolor": SURFACE,
    "axes.facecolor": SURFACE,
    "savefig.facecolor": SURFACE,
    "pdf.fonttype": 42,
})


def read_dem(path: Path) -> tuple[np.ndarray, list[float]]:
    """Elevation (m, north row first) and its [west, east, south, north] extent, from the GeoTIFF tags."""
    image = Image.open(path)
    scale_x, scale_y, _ = image.tag_v2[33550]          # ModelPixelScaleTag
    _, _, _, west, north, _ = image.tag_v2[33922]      # ModelTiepointTag: the upper-left corner
    keys = image.tag_v2[34735]                         # GeoKeyDirectoryTag
    geokeys = {keys[i]: keys[i + 3] for i in range(4, len(keys), 4)}
    if geokeys.get(2048) != 4326 or geokeys.get(1025) != 1:
        raise ValueError(f"{path.name}: expected WGS 84 (EPSG:4326) and pixel-is-area cells")
    elevation = np.array(image, dtype=float)
    nodata = image.tag_v2.get(42113)
    if nodata is not None:
        elevation[elevation == float(nodata)] = np.nan
    rows, cols = elevation.shape
    return elevation, [west, west + cols * scale_x, north - rows * scale_y, north]


def read_gauges() -> list[dict]:
    """Every gauge of the study, with its landform region; refuses a gauge without one."""
    with open(LANDFORM, newline="", encoding="utf-8") as handle:
        landform = {row["station_id"]: row for row in csv.DictReader(handle)}
    with open(STATIONS, newline="", encoding="utf-8") as handle:
        stations = list(csv.DictReader(handle))
    missing = [row["station_id"] for row in stations if row["station_id"] not in landform]
    if missing:
        raise ValueError(f"{len(missing)} gauges have no landform region: {', '.join(missing[:5])}")
    return [{"lon": float(row["lon"]), "lat": float(row["lat"]),
             "region": landform[row["station_id"]]["relief_region"],
             "insufficient": landform[row["station_id"]]["insufficient_gauges"] == "true",
             "window_km": float(landform[row["station_id"]]["window_side_km"])} for row in stations]


def scale_bar(ax: plt.Axes, lon: float, lat: float, km: float) -> None:
    """A `km`-long bar starting at (lon, lat), its length in degrees of longitude at that latitude."""
    degrees = np.degrees(km / (EARTH_RADIUS_KM * np.cos(np.radians(lat))))
    ax.plot([lon, lon + degrees], [lat, lat], color=INK, linewidth=1.6, solid_capstyle="butt", zorder=4)
    for x, label in ((lon, "0"), (lon + degrees, f"{km:g} km")):
        ax.text(x, lat + 0.035, label, ha="center", va="bottom", fontsize=6, color=INK, zorder=4)


def fig1_study_area(elevation: np.ndarray, extent: list[float], gauges: list[dict]) -> None:
    west, east, south, north = extent
    mid_lat = (south + north) / 2
    rows, cols = elevation.shape
    # Cell size in metres, so the hillshade slopes are true slopes.
    dy = np.radians((north - south) / rows) * EARTH_RADIUS_KM * 1000
    dx = dy * np.cos(np.radians(mid_lat))
    top = 500 * np.ceil(np.nanmax(elevation) / 500)
    norm = Normalize(0, top)
    shaded = LightSource(azdeg=315, altdeg=45).shade(np.nan_to_num(elevation, nan=0.0), cmap=ELEVATION, norm=norm,
                                                       blend_mode="soft", vert_exag=2, dx=dx, dy=dy)

    fig = plt.figure(figsize=(5.6, 6.2))
    ax = fig.add_axes([0.1, 0.1, 0.72, 0.77])
    ax.imshow(shaded, extent=extent, origin="upper", interpolation="bilinear", zorder=1)
    counts = {region: sum(g["region"] == region for g in gauges) for region in REGIONS}
    handles = []
    for region in REGIONS:
        members = [g for g in gauges if g["region"] == region]
        hollow = any(g["insufficient"] for g in members)
        ax.scatter([g["lon"] for g in members], [g["lat"] for g in members], marker=REGION_MARKERS[region], s=13,
                   facecolors=SURFACE if hollow else INK, edgecolors=INK if hollow else SURFACE, linewidths=0.6,
                   zorder=3)
        handles.append(Line2D([], [], linestyle="none", marker=REGION_MARKERS[region], markersize=4.5,
                              markerfacecolor=SURFACE if hollow else INK, markeredgecolor=INK if hollow else SURFACE,
                              label=f"{REGION_LABELS[region]} · {counts[region]} gauges"
                                    + (" (too few to interpret)" if hollow else "")))
    scale_bar(ax, west + 0.08, south + 0.07, 50)
    ax.set_xlim(west, east)
    ax.set_ylim(south, north)
    ax.set_aspect(1 / np.cos(np.radians(mid_lat)), adjustable="box")
    xticks = np.arange(np.ceil(west * 2) / 2, east + 1e-9, 0.5)
    yticks = np.arange(np.ceil(south * 2) / 2, north + 1e-9, 0.5)
    ax.set_xticks(xticks, [f"{x:g}°E" for x in xticks])
    ax.set_yticks(yticks, [f"{y:g}°N" for y in yticks])
    ax.tick_params(length=2, pad=1.5, labelsize=6)
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.045), ncol=3, frameon=False,
              fontsize=6, handletextpad=0.3, columnspacing=1.2)

    bar = fig.add_axes([0.85, 0.3, 0.02, 0.4])
    colorbar = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=ELEVATION), cax=bar,
                            ticks=np.arange(0, top + 1, 500))
    colorbar.outline.set_edgecolor(AXIS)
    colorbar.outline.set_linewidth(0.6)
    colorbar.ax.tick_params(length=2, labelsize=6)
    colorbar.set_label("Elevation (m)", fontsize=6.5, color=INK_SECONDARY)

    windows = sorted({g["window_km"] for g in gauges})
    fig.suptitle(f"Study area and the {len(gauges)} hourly gauges", x=0.1, ha="left", fontsize=9, y=0.975)
    fig.text(0.1, 0.94, "Shaded elevation from a 1 km SRTM DEM. Gauge shapes: landform region by local relief\n"
             f"(highest minus lowest point) in a {'/'.join(f'{w:g}' for w in windows)} km window of the "
             "Copernicus GLO-30 DEM.", ha="left", va="top", fontsize=6.3, color=INK_SECONDARY, linespacing=1.4)
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG_DIR / "fig1_study_area.png", dpi=300)
    fig.savefig(FIG_DIR / "fig1_study_area.pdf")
    plt.close(fig)
    print(f"wrote {FIG_DIR / 'fig1_study_area.png'}")


def main() -> None:
    elevation, extent = read_dem(DEM)
    gauges = read_gauges()
    outside = [g for g in gauges if not (extent[0] <= g["lon"] <= extent[1] and extent[2] <= g["lat"] <= extent[3])]
    if outside:
        raise ValueError(f"{len(outside)} gauges lie outside the DEM's extent {extent}")
    fig1_study_area(elevation, extent, gauges)
    write_plot_provenance(FIG_DIR, ROOT, inputs=(DEM, STATIONS, LANDFORM))


if __name__ == "__main__":
    main()

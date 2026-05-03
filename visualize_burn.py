#!/usr/bin/env python3
"""
visualize_burn.py
Reads burn_times.csv and renders a topographic arrival-time map.

Usage:
    python3 visualize_burn.py [burn_times.csv] [output.png]

Defaults to burn_times.csv and burn_map.png in the current directory.
"""

import sys
import csv
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.ticker as ticker
from matplotlib.patches import FancyArrowPatch
matplotlib.use("Agg")

# ---------------------------------------------------------------------------
# Config — edit these to match your simulation parameters
# ---------------------------------------------------------------------------
CELL_SIZE_FT   = 100.0   # ft per cell (must match main.cpp)
WIND_DIR_DEG   = 270.0   # wind FROM direction (met convention, 0=N CW)
WIND_SPEED_MPH = 10.0    # mi/hr
ASPECT_DEG     = 180.0   # slope aspect (direction slope faces, 0=N CW)
SLOPE_PCT      = 20.0    # slope percent (tan * 100)
FUEL_MODEL     = "FM2 — Timber Grass"
# ---------------------------------------------------------------------------

def load_csv(path: str):
    rows, cols, times = [], [], []
    with open(path) as f:
        reader = csv.DictReader(f)
        for rec in reader:
            rows.append(int(rec["row"]))
            cols.append(int(rec["col"]))
            times.append(float(rec["burn_time_min"]))
    return rows, cols, times


def build_grid(rows, cols, times):
    n_rows = max(rows) + 1
    n_cols = max(cols) + 1
    grid = np.full((n_rows, n_cols), np.nan)
    for r, c, t in zip(rows, cols, times):
        if t >= 0.0:
            grid[r, c] = t
    return grid


def wind_arrow_direction(wind_from_deg):
    """Return (dx, dy) unit vector the wind is blowing TOWARD (display coords: y down)."""
    head_deg = (wind_from_deg + 180.0) % 360.0
    rad = math.radians(90.0 - head_deg)   # math convention
    dx =  math.cos(rad)
    dy = -math.sin(rad)                    # flip y for image coords (row 0 = top)
    return dx, dy


def slope_arrow_direction(aspect_deg):
    """Return (dx, dy) unit vector pointing upslope (toward aspect, display coords)."""
    rad = math.radians(90.0 - aspect_deg)
    dx =  math.cos(rad)
    dy = -math.sin(rad)
    return dx, dy


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else "burn_times.csv"
    out_path = sys.argv[2] if len(sys.argv) > 2 else "burn_map.png"

    print(f"Loading {csv_path} ...")
    rows, cols, times = load_csv(csv_path)
    grid = build_grid(rows, cols, times)
    n_rows, n_cols = grid.shape

    burned_times = grid[~np.isnan(grid)]
    t_max = burned_times.max() if len(burned_times) > 0 else 1.0
    t_min = burned_times.min() if len(burned_times) > 0 else 0.0
    burned_cells = np.sum(~np.isnan(grid))
    burned_acres = burned_cells * CELL_SIZE_FT**2 / 43560.0
    total_acres  = n_rows * n_cols * CELL_SIZE_FT**2 / 43560.0

    # Physical extents in feet for axis labels
    extent_ft = [0, n_cols * CELL_SIZE_FT, n_rows * CELL_SIZE_FT, 0]  # left right bottom top

    # ---------------------------------------------------------------------------
    # Figure layout
    # ---------------------------------------------------------------------------
    fig = plt.figure(figsize=(12, 9), facecolor="#0f0f0f")
    ax  = fig.add_axes([0.07, 0.10, 0.78, 0.82])
    cax = fig.add_axes([0.87, 0.10, 0.025, 0.82])

    ax.set_facecolor("#111111")

    # ---------------------------------------------------------------------------
    # Unburned background
    # ---------------------------------------------------------------------------
    unburned_bg = np.ones((n_rows, n_cols, 4))
    unburned_bg[:, :, :3] = np.array([0.13, 0.13, 0.13])  # dark grey
    unburned_bg[:, :, 3]  = 1.0
    ax.imshow(unburned_bg, extent=extent_ft, origin="upper", zorder=0)

    # ---------------------------------------------------------------------------
    # Arrival time colormap — perceptually uniform, fire-themed
    # Uses a custom colormap: deep purple → orange → bright yellow (like embers)
    # ---------------------------------------------------------------------------
    cmap_colors = [
        (0.20, 0.02, 0.35),   # deep purple  (ignition)
        (0.55, 0.05, 0.15),   # dark crimson
        (0.85, 0.22, 0.02),   # burnt orange
        (0.98, 0.60, 0.02),   # amber
        (0.99, 0.95, 0.40),   # pale yellow  (latest arrival)
    ]
    cmap = mcolors.LinearSegmentedColormap.from_list("fire_topo", cmap_colors, N=512)

    masked = np.ma.masked_invalid(grid)
    img = ax.imshow(
        masked,
        extent=extent_ft,
        origin="upper",
        cmap=cmap,
        vmin=t_min,
        vmax=t_max,
        interpolation="bilinear",
        zorder=1,
        alpha=0.92,
    )

    # ---------------------------------------------------------------------------
    # Topographic contour lines (isochrones) — every 5 minutes
    #
    # Build coordinate arrays in ft for contourf/contour.
    # imshow origin=upper means row 0 is at y=n_rows*cell (top in ft coords),
    # so we need to flip the grid vertically for contour (which uses y-up).
    # ---------------------------------------------------------------------------
    x_ft = np.linspace(CELL_SIZE_FT / 2, n_cols * CELL_SIZE_FT - CELL_SIZE_FT / 2, n_cols)
    y_ft = np.linspace(CELL_SIZE_FT / 2, n_rows * CELL_SIZE_FT - CELL_SIZE_FT / 2, n_rows)
    X, Y = np.meshgrid(x_ft, y_ft)

    # For contour, flip Y so row 0 = top = largest y value
    Y_flip = n_rows * CELL_SIZE_FT - Y
    grid_flip = np.flipud(grid)

    contour_interval = max(1.0, round(t_max / 10.0 / 5.0) * 5.0)  # ~10 isochrones, snapped to 5
    levels = np.arange(
        math.ceil(t_min / contour_interval) * contour_interval,
        t_max + contour_interval,
        contour_interval,
    )

    if len(levels) > 1:
        cs = ax.contour(
            X, Y_flip, grid_flip,
            levels=levels,
            colors="white",
            linewidths=0.6,
            alpha=0.35,
            zorder=2,
        )
        ax.clabel(cs, fmt=lambda v: f"{v:.0f} min", fontsize=6.5,
                  colors="white", inline=True, inline_spacing=4)

    # ---------------------------------------------------------------------------
    # Ignition point marker
    # ---------------------------------------------------------------------------
    ign_mask = (grid == t_min)
    ign_rows, ign_cols = np.where(ign_mask)
    if len(ign_rows) > 0:
        ign_r, ign_c = ign_rows[0], ign_cols[0]
        ign_x = (ign_c + 0.5) * CELL_SIZE_FT
        ign_y_flip = (n_rows - ign_r - 0.5) * CELL_SIZE_FT
        ax.plot(ign_x, ign_y_flip, "w+", markersize=12, markeredgewidth=1.5, zorder=5)
        ax.plot(ign_x, ign_y_flip, "wo", markersize=6, fillstyle="none",
                markeredgewidth=1.0, zorder=5)

    # ---------------------------------------------------------------------------
    # Wind and slope arrows (upper-right corner legend)
    # ---------------------------------------------------------------------------
    arrow_x = n_cols * CELL_SIZE_FT * 0.88
    arrow_y = n_rows * CELL_SIZE_FT * 0.12
    arrow_len = n_cols * CELL_SIZE_FT * 0.07

    wdx, wdy = wind_arrow_direction(WIND_DIR_DEG)
    sdx, sdy = slope_arrow_direction(ASPECT_DEG)

    # Wind arrow — blue
    ax.annotate(
        "", xy=(arrow_x + wdx * arrow_len, arrow_y + wdy * arrow_len),
        xytext=(arrow_x, arrow_y),
        arrowprops=dict(arrowstyle="-|>", color="#4fc3f7", lw=1.8),
        zorder=6,
    )
    ax.text(arrow_x + wdx * arrow_len * 1.25, arrow_y + wdy * arrow_len * 1.25,
            f"Wind\n{WIND_SPEED_MPH} mph", color="#4fc3f7", fontsize=7.5,
            ha="center", va="center", zorder=6)

    # Slope arrow — tan/khaki
    slope_ox = arrow_x - arrow_len * 0.3
    ax.annotate(
        "", xy=(slope_ox + sdx * arrow_len, arrow_y + sdy * arrow_len),
        xytext=(slope_ox, arrow_y),
        arrowprops=dict(arrowstyle="-|>", color="#d4b483", lw=1.8),
        zorder=6,
    )
    ax.text(slope_ox + sdx * arrow_len * 1.25, arrow_y + sdy * arrow_len * 1.25,
            f"Upslope\n{SLOPE_PCT:.0f}%", color="#d4b483", fontsize=7.5,
            ha="center", va="center", zorder=6)

    # ---------------------------------------------------------------------------
    # Axes formatting
    # ---------------------------------------------------------------------------
    def ft_to_label(val, _):
        return f"{val / 5280:.2f} mi" if val >= 5280 else f"{val:.0f} ft"

    ax.xaxis.set_major_formatter(ticker.FuncFormatter(ft_to_label))
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(ft_to_label))
    ax.tick_params(colors="#aaaaaa", labelsize=8)
    for spine in ax.spines.values():
        spine.set_edgecolor("#444444")

    ax.set_xlabel("Easting", color="#aaaaaa", fontsize=9)
    ax.set_ylabel("Northing", color="#aaaaaa", fontsize=9)

    # ---------------------------------------------------------------------------
    # Colorbar
    # ---------------------------------------------------------------------------
    cb = fig.colorbar(img, cax=cax)
    cb.set_label("Arrival time [min]", color="#cccccc", fontsize=9)
    cb.ax.yaxis.set_tick_params(color="#aaaaaa", labelsize=8)
    plt.setp(cb.ax.yaxis.get_ticklabels(), color="#aaaaaa")
    cb.outline.set_edgecolor("#444444")

    # ---------------------------------------------------------------------------
    # Title + stats
    # ---------------------------------------------------------------------------
    ax.set_title(
        f"Fire Arrival Time — {FUEL_MODEL}",
        color="white", fontsize=13, pad=10, fontweight="bold",
    )

    stats_str = (
        f"Burned: {burned_acres:.1f} ac / {total_acres:.1f} ac total  |  "
        f"Max arrival: {t_max:.1f} min  |  "
        f"Cell: {CELL_SIZE_FT:.0f} ft"
    )
    fig.text(0.07, 0.035, stats_str, color="#888888", fontsize=8)

    # ---------------------------------------------------------------------------
    # Save
    # ---------------------------------------------------------------------------
    fig.savefig(out_path, dpi=180, bbox_inches="tight", facecolor=fig.get_facecolor())
    print(f"Saved → {out_path}")
    print(f"  Grid : {n_rows} × {n_cols} cells")
    print(f"  Burned: {burned_cells} cells ({burned_acres:.2f} ac)")
    print(f"  Arrival time range: {t_min:.2f} – {t_max:.2f} min")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Basic visualization for generated topography NetCDF.

- Elevation (m)
- Land mask (0/1)
- Base albedo (unitless)
- Surface friction (s^-1)

Usage:
  python -m scripts.plot_topography                   # auto-pick latest *.nc in ./data
  python -m scripts.plot_topography --file data/xxx.nc

Outputs:
  - Writes PNG(s) to the same directory as the input .nc:
      <stem>_overview.png
"""

import argparse
import glob
import os
from datetime import datetime
from typing import Optional

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset


def find_latest_nc(data_dir: str = "data") -> Optional[str]:
    paths = sorted(glob.glob(os.path.join(data_dir, "*.nc")))
    return paths[-1] if paths else None


def load_fields(nc_path: str):
    with Dataset(nc_path) as ds:
        lat = ds["lat"][:]
        lon = ds["lon"][:]
        elevation = ds["elevation"][:]
        land_mask = ds["land_mask"][:]
        base_albedo = ds["base_albedo"][:]
        friction = ds["friction"][:]
        # Optional attrs
        sea_level_m = getattr(ds, "sea_level_m", np.nan)
        target_land_fraction = getattr(ds, "target_land_fraction", np.nan)
    return lat, lon, elevation, land_mask, base_albedo, friction, sea_level_m, target_land_fraction


def plot_overview(lat, lon, elevation, land_mask, base_albedo, friction, out_png: str, sea_level_m: float, target_land_frac: float):
    LON, LAT = np.meshgrid(lon, lat)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)

    # Elevation
    ax = axes[0, 0]
    elev_plot = ax.pcolormesh(LON, LAT, elevation, shading="auto", cmap="terrain")
    cb = fig.colorbar(elev_plot, ax=ax, orientation="vertical", pad=0.02)
    cb.set_label("m")
    title = "Elevation"
    if np.isfinite(sea_level_m):
        title += f" (sea level ≈ {sea_level_m:.0f} m)"
    ax.set_title(title)
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")

    # Land mask
    ax = axes[0, 1]
    mask_plot = ax.pcolormesh(LON, LAT, land_mask, shading="auto", cmap="Greys", vmin=0, vmax=1)
    cb = fig.colorbar(mask_plot, ax=ax, orientation="vertical", pad=0.02)
    cb.set_label("0=ocean, 1=land")
    title = "Land-Sea Mask"
    if np.isfinite(target_land_frac):
        # Compute achieved land fraction with cos(lat) weights
        area_w = np.cos(np.deg2rad(LAT))
        achieved = float((area_w * (land_mask == 1)).sum() / (area_w.sum() + 1e-15))
        title += f" (target {target_land_frac:.2f}, achieved {achieved:.2f})"
    ax.set_title(title)
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")

    # Base albedo
    ax = axes[1, 0]
    alb_plot = ax.pcolormesh(LON, LAT, base_albedo, shading="auto", cmap="cividis", vmin=0, vmax=0.85)
    cb = fig.colorbar(alb_plot, ax=ax, orientation="vertical", pad=0.02)
    cb.set_label("unitless")
    ax.set_title("Base Albedo (ice-free)")
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")

    # Friction
    ax = axes[1, 1]
    fr_plot = ax.pcolormesh(LON, LAT, friction, shading="auto", cmap="magma")
    cb = fig.colorbar(fr_plot, ax=ax, orientation="vertical", pad=0.02)
    cb.set_label("s$^{-1}$")
    ax.set_title("Surface Friction Coefficient")
    ax.set_xlabel("Longitude (deg)")
    ax.set_ylabel("Latitude (deg)")

    fig.suptitle("Qingdai Topography & Surface Properties (Project 004)", fontsize=14)
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot basic visualizations for Qingdai topography NetCDF.")
    parser.add_argument("--file", "-f", type=str, default=None, help="Path to NetCDF file under ./data")
    args = parser.parse_args()

    nc_path = args.file or find_latest_nc("data")
    if not nc_path or not os.path.exists(nc_path):
        raise SystemExit("No NetCDF found. Run `python -m scripts.generate_topography` first, or pass --file data/xxx.nc.")

    print(f"[Plot] Using file: {nc_path}")
    lat, lon, elevation, land_mask, base_albedo, friction, sea_level_m, target_land_frac = load_fields(nc_path)

    out_dir = os.path.dirname(nc_path)
    stem = os.path.splitext(os.path.basename(nc_path))[0]
    out_png = os.path.join(out_dir, f"{stem}_overview.png")

    print(f"[Plot] Writing: {out_png}")
    plot_overview(lat, lon, elevation, land_mask, base_albedo, friction, out_png, sea_level_m, target_land_frac)
    print("[Plot] Done.")


if __name__ == "__main__":
    main()

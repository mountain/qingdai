#!/usr/bin/env python3
"""
Generate procedural topography, land-sea mask, and base surface properties,
then export a multi-field NetCDF file under ./data/ for GCM and other projects.

Implements Project 004 (Milestone 1: L1+L3 + sea-level adaptation).

Environment variables (optional):
  QD_N_LAT                int   default 181
  QD_N_LON                int   default 360
  QD_SEED                 int   default 42
  QD_TARGET_LAND_FRAC     float default 0.40

  # Optional parameter overrides passed to generate_elevation_map:
  QD_N_CONTINENTS         int   default 3
  QD_CONT_SIGMA_DEG       float default 30.0
  QD_CONT_SHAPE_P         float default 2.0
  QD_CONT_MIN_DIST_DEG    float default 40.0
  QD_W_VLF                float default 0.35
  QD_FBM_OCTAVES          int   default 5
  QD_HURST_H              float default 0.8
  QD_W1                   float default 1.0
  QD_W3                   float default 0.6
  QD_SCALE_M              float default 4500.0

Usage:
  python scripts/generate_topography.py
"""
import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from datetime import datetime

from pygcm.grid import SphericalGrid
from pygcm.topography import (
    generate_elevation_map,
    create_land_sea_mask_from_elevation,
    generate_base_properties,
    export_topography_to_netcdf,
)


def getenv_int(name: str, default: int) -> int:
    try:
        return int(os.getenv(name, default))
    except Exception:
        return default


def getenv_float(name: str, default: float) -> float:
    try:
        return float(os.getenv(name, default))
    except Exception:
        return default


def main():
    # Grid/config
    n_lat = getenv_int("QD_N_LAT", 181)
    n_lon = getenv_int("QD_N_LON", 360)
    seed = getenv_int("QD_SEED", 42)
    target_land_frac = getenv_float("QD_TARGET_LAND_FRAC", 0.40)

    # Elevation param overrides
    params = {
        "N_CONTINENTS": getenv_int("QD_N_CONTINENTS", 3),
        "CONTINENT_SIGMA_DEG": getenv_float("QD_CONT_SIGMA_DEG", 30.0),
        "CONTINENT_SHAPE_P": getenv_float("QD_CONT_SHAPE_P", 2.0),
        "CONT_MIN_DIST_DEG": getenv_float("QD_CONT_MIN_DIST_DEG", 40.0),
        "W_VLF": getenv_float("QD_W_VLF", 0.35),
        "FBM_OCTAVES": getenv_int("QD_FBM_OCTAVES", 5),
        "HURST_H": getenv_float("QD_HURST_H", 0.8),
        "W1": getenv_float("QD_W1", 1.0),
        "W3": getenv_float("QD_W3", 0.6),
        "SCALE_M": getenv_float("QD_SCALE_M", 4500.0),
    }

    print(f"[Topo] Grid {n_lat}x{n_lon}, seed={seed}, target_land_frac={target_land_frac}")
    print(f"[Topo] Params: {params}")

    # Build grid
    grid = SphericalGrid(n_lat=n_lat, n_lon=n_lon)

    # Generate fields
    elevation = generate_elevation_map(grid, seed=seed, params=params)
    land_mask, sea_level_m = create_land_sea_mask_from_elevation(
        elevation, grid, target_land_frac=target_land_frac
    )
    base_albedo, friction = generate_base_properties(land_mask, elevation=elevation, grid=grid)

    # Output path
    os.makedirs("data", exist_ok=True)
    timestamp = datetime.utcnow().strftime("%Y%m%dT%H%M%SZ")
    out_name = f"topography_qingdai_{n_lat}x{n_lon}_seed{seed}_{timestamp}.nc"
    out_path = os.path.join("data", out_name)

    print(f"[Topo] Exporting to NetCDF: {out_path}")
    export_topography_to_netcdf(
        grid=grid,
        elevation=elevation,
        land_mask=land_mask,
        base_albedo=base_albedo,
        friction=friction,
        sea_level_m=sea_level_m,
        out_path=out_path,
    )
    print("[Topo] Done.")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"[Topo] ERROR: {e}", file=sys.stderr)
        sys.exit(1)

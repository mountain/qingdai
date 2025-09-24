# -*- coding: utf-8 -*-
"""
Analyze topography NetCDF to report:
- Maximum elevation (m)
- Minimum depth (m, absolute of most negative elevation)
- Mean slope (area-weighted) in degrees for global/land/ocean
Also saves results to output/topography_stats.json and .txt
"""
import json
import glob
import os
import sys
from typing import Dict, Any

import numpy as np

try:
    import netCDF4 as nc
except Exception as e:
    print("ERROR: netCDF4 is required. Please install dependencies from requirements.txt")
    raise

# Planet radius for metric conversion (meters)
def get_planet_radius() -> float:
    try:
        from pygcm.constants import PLANET_RADIUS as R
        if not np.isfinite(R) or R <= 0:
            raise ValueError("Invalid PLANET_RADIUS in pygcm.constants")
        return float(R)
    except Exception:
        # Fallback to Earth radius if not available
        return 6.371e6


def find_topography_file() -> str:
    candidates = sorted(
        glob.glob("data/topography_*.nc") + glob.glob("data/topography_qingdai_*.nc")
    )
    if not candidates:
        raise FileNotFoundError("No topography NetCDF found under data/.")
    return candidates[-1]


def area_weights(lat: np.ndarray, shape2d) -> np.ndarray:
    """
    Construct area weights proportional to cos(lat) for each grid cell.
    Assumes regular lat-lon grid.
    """
    w = np.cos(np.deg2rad(lat))
    w = np.clip(w, 0.0, None)  # guard numerical tiny negatives
    return np.broadcast_to(w[:, None], shape2d)


def compute_slope_metrics(H: np.ndarray, lat: np.ndarray, lon: np.ndarray, R: float) -> Dict[str, np.ndarray]:
    """
    Compute central-difference slope on a sphere and return:
      - slope angle in degrees
      - dimensionless slope magnitude |∇H|
    Uses periodic boundary in longitude and one-sided at the lat edges.
    """
    lat_rad = np.deg2rad(lat)
    lon_rad = np.deg2rad(lon)

    # Use mean spacing (regular grid expected)
    dphi = float(np.mean(np.diff(lat_rad)))
    dlmb = float(np.mean(np.diff(lon_rad)))

    # Metric factors
    phi2d = np.deg2rad(np.broadcast_to(lat[:, None], H.shape))
    cosphi = np.clip(np.cos(phi2d), 1e-6, None)  # avoid singularity at poles

    # Neighbors
    H_e = np.roll(H, -1, axis=1)
    H_w = np.roll(H,  1, axis=1)
    H_n = np.empty_like(H)
    H_s = np.empty_like(H)
    # North (toward +lat), South (toward -lat)
    H_n[:-1, :] = H[1:, :]
    H_n[-1,  :] = H[-1, :]
    H_s[1:,  :] = H[:-1, :]
    H_s[0,   :] = H[0,  :]

    # Gradients on sphere (m/m -> dimensionless)
    dHdx = (H_e - H_w) / (2.0 * R * dlmb * cosphi)
    dHdy = (H_n - H_s) / (2.0 * R * dphi)
    S = np.sqrt(dHdx**2 + dHdy**2)
    theta_deg = np.degrees(np.arctan(S))

    return {"theta_deg": theta_deg, "S": S}


def wmean(field: np.ndarray, weights: np.ndarray, mask: np.ndarray = None) -> float:
    if mask is not None:
        weights = weights * mask.astype(np.float64)
    num = np.nansum(field * weights)
    den = np.nansum(weights)
    return float(num / den) if den > 0 else float("nan")


def analyze() -> Dict[str, Any]:
    R = get_planet_radius()
    path = find_topography_file()

    with nc.Dataset(path) as ds:
        H = ds.variables["elevation"][:].astype(np.float64)
        lat = ds.variables["lat"][:].astype(np.float64)
        lon = ds.variables["lon"][:].astype(np.float64)
        if "land_mask" in ds.variables:
            land_mask = ds.variables["land_mask"][:].astype(np.int8)
        else:
            # Fallback: elevation >= 0 as land
            land_mask = (H >= 0).astype(np.int8)

        # Optional attribute for sea level
        sea_level = getattr(ds, "sea_level_m", 0.0)

    h_min = float(np.nanmin(H))
    h_max = float(np.nanmax(H))

    # Depth is absolute of the most negative elevation relative to sea level
    # If sea_level attribute is present and meaningful, compute depth accordingly
    # Otherwise use 0 m as reference.
    ref_sea = float(sea_level)
    min_depth = float(max(0.0, ref_sea - h_min))

    slope = compute_slope_metrics(H, lat, lon, R)
    theta_deg = slope["theta_deg"]
    S = slope["S"]

    weights = area_weights(lat, H.shape)
    ocean_mask = 1 - land_mask

    # Area-weighted means
    avg_slope_deg_global = wmean(theta_deg, weights)
    avg_slope_deg_land   = wmean(theta_deg, weights, land_mask)
    avg_slope_deg_ocean  = wmean(theta_deg, weights, ocean_mask)

    avg_S_global = wmean(S, weights)
    avg_S_land   = wmean(S, weights, land_mask)
    avg_S_ocean  = wmean(S, weights, ocean_mask)

    result = {
        "file": path,
        "planet_radius_m": R,
        "sea_level_m": ref_sea,
        "elevation_max_m": h_max,
        "elevation_min_m": h_min,
        "min_depth_m": min_depth,  # absolute depth of deepest point below sea level
        "mean_slope_angle_deg": {
            "global": avg_slope_deg_global,
            "land":   avg_slope_deg_land,
            "ocean":  avg_slope_deg_ocean,
        },
        "mean_dimensionless_slope": {
            "global": avg_S_global,
            "land":   avg_S_land,
            "ocean":  avg_S_ocean,
        },
        "earth_reference": {
            "elevation_max_m": 8849.0,   # Everest
            "min_depth_m":      10994.0, # Mariana Trench (approx)
            # Global mean slope (cell-scale) depends on resolution/definition; typically small (<~1–3 deg)
        },
    }
    return result


def main():
    try:
        res = analyze()
    except Exception as e:
        print(f"ERROR: {e}")
        sys.exit(1)

    os.makedirs("output", exist_ok=True)

    # Write JSON
    json_path = "output/topography_stats.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=2, ensure_ascii=False)

    # Write TXT summary
    txt_path = "output/topography_stats.txt"
    with open(txt_path, "w", encoding="utf-8") as f:
        f.write(f"Using: {res['file']}\n")
        f.write(f"Planet radius (m): {res['planet_radius_m']:.2f}\n")
        f.write(f"Sea level (m): {res['sea_level_m']:.2f}\n")
        f.write(f"Elevation max (m): {res['elevation_max_m']:.2f}\n")
        f.write(f"Elevation min (m): {res['elevation_min_m']:.2f}\n")
        f.write(f"Min depth (m): {res['min_depth_m']:.2f}\n")
        f.write(f"Mean slope angle (deg) - global: {res['mean_slope_angle_deg']['global']:.4f}\n")
        f.write(f"Mean slope angle (deg) - land:   {res['mean_slope_angle_deg']['land']:.4f}\n")
        f.write(f"Mean slope angle (deg) - ocean:  {res['mean_slope_angle_deg']['ocean']:.4f}\n")
        f.write(f"Mean |grad H| (dimensionless) - global: {res['mean_dimensionless_slope']['global']:.6f}\n")
        f.write(f"Mean |grad H| (dimensionless) - land:   {res['mean_dimensionless_slope']['land']:.6f}\n")
        f.write(f"Mean |grad H| (dimensionless) - ocean:  {res['mean_dimensionless_slope']['ocean']:.6f}\n")
        f.write("\nEarth reference:\n")
        f.write(f"  Elevation max (m): {res['earth_reference']['elevation_max_m']:.0f}\n")
        f.write(f"  Min depth (m):     {res['earth_reference']['min_depth_m']:.0f}\n")

    # Also print JSON to stdout for convenience
    print(json.dumps(res, indent=2, ensure_ascii=False))
    print(f"\nWrote:\n  {json_path}\n  {txt_path}")


if __name__ == "__main__":
    main()

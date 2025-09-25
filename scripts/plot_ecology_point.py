#!/usr/bin/env python3
"""
Plot ecology distributions at a given location:

- Select top-3 species by LAI (Σ_k LAI_s,k) at the target cell
- For each species, draw 4 subplots (total 12):
  (1) Canopy height distribution (in a neighborhood)
  (2) Leaf area increment distribution (samples) [m^2/day]  – requires running IndividualPool with E_day
  (3) Root development proxy (alloc_root · E_day)           – requires running IndividualPool with E_day
  (4) Lifespan (days) as a vertical line

If IndividualPool is not available (offline use), (2)(3) will be empty histograms, while (1)(4) are still shown.

Usage:
  python3 -m scripts.plot_ecology_point --lat 35.0 --lon 110.0 --out output/top3_ecology.png

Optional:
  --eco-npz data/eco_autosave.npz         # to load LAI and species weights
  --topo data/topography_*.nc             # to load grid lat/lon and land mask
  --nbins 24 --neigh 1                    # histogram bins and neighborhood radius
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from typing import Optional, Tuple

import numpy as np

# Matplotlib is only needed for plotting; defer import error to runtime
try:
    import matplotlib.pyplot as plt  # type: ignore
except Exception as _e:
    plt = None

# Optional NetCDF4 for reading topography dims
try:
    from netCDF4 import Dataset  # type: ignore
except Exception:
    Dataset = None

# Project imports
try:
    from pygcm.topography import load_topography_from_netcdf, create_land_sea_mask
except Exception:
    load_topography_from_netcdf = None
    create_land_sea_mask = None

try:
    from pygcm.grid import SphericalGrid
    from pygcm.ecology import EcologyAdapter
except Exception as e:
    print(f"[Error] Failed to import core modules: {e}", file=sys.stderr)
    sys.exit(1)

try:
    from pygcm.ploter import plot_top3_species_distributions
except Exception as e:
    print(f"[Error] Failed to import plotting utility from pygcm.ploter: {e}", file=sys.stderr)
    sys.exit(1)


def _require_matplotlib():
    if plt is None:
        raise RuntimeError("matplotlib is required. Please install 'matplotlib'.")


def _latest_file(pattern: str) -> Optional[str]:
    files = glob.glob(pattern)
    if not files:
        return None
    files.sort(key=lambda p: os.path.getmtime(p), reverse=True)
    return files[0]


def _peek_eco_npz_dims(path_npz: str) -> Optional[Tuple[int, int]]:
    """
    Peek LAI dims from eco_autosave npz without building adapter.
    Returns (H,W) or None.
    """
    try:
        data = np.load(path_npz)
        LAI = np.asarray(data.get("LAI"))
        if LAI is None or LAI.ndim != 2:
            return None
        return int(LAI.shape[0]), int(LAI.shape[1])
    except Exception:
        return None


def _read_nc_dims(path_nc: str) -> Optional[Tuple[int, int]]:
    """
    Read (H,W) dims from NetCDF by 'lat' and 'lon' variables if available.
    """
    if Dataset is None:
        return None
    try:
        with Dataset(path_nc) as ds:
            lat = np.asarray(ds.variables["lat"][:], dtype=float)
            lon = np.asarray(ds.variables["lon"][:], dtype=float)
        return int(lat.size), int(lon.size)
    except Exception:
        return None


def _load_grid_and_mask(topo_path: Optional[str],
                        target_hw: Optional[Tuple[int, int]] = None,
                        H_default: int = 181,
                        W_default: int = 360) -> Tuple[SphericalGrid, np.ndarray]:
    """
    Build SphericalGrid and land_mask.

    Preference:
      1) If topo_path provided and dims readable, use those dims only if they match target_hw (if specified).
      2) Else if target_hw provided, build grid with target_hw and procedural land mask.
      3) Else fallback to (H_default,W_default) grid and procedural land mask.
    """
    topo_dims = _read_nc_dims(topo_path) if topo_path else None
    if topo_path and topo_dims and (target_hw is None or topo_dims == target_hw):
        H, W = topo_dims
        grid = SphericalGrid(H, W)
        if load_topography_from_netcdf is not None:
            try:
                elevation, land_mask, base_albedo, friction = load_topography_from_netcdf(topo_path, grid)  # type: ignore
                return grid, np.asarray(land_mask, dtype=int)
            except Exception as e:
                print(f"[Topo] Failed to load '{topo_path}': {e}; falling back to procedural mask.", file=sys.stderr)
        # procedural mask fallback
        if create_land_sea_mask is not None:
            try:
                land_mask = create_land_sea_mask(grid)
                return grid, np.asarray(land_mask, dtype=int)
            except Exception:
                pass
        return grid, np.ones((H, W), dtype=int)

    # No topo or mismatched topo → use target_hw or defaults
    if target_hw is not None:
        Ht, Wt = target_hw
        grid = SphericalGrid(int(Ht), int(Wt))
    else:
        grid = SphericalGrid(H_default, W_default)
    if create_land_sea_mask is not None:
        try:
            land_mask = create_land_sea_mask(grid)
            return grid, np.asarray(land_mask, dtype=int)
        except Exception:
            pass
    return grid, np.ones(np.asarray(grid.lat_mesh).shape, dtype=int)


def _load_eco_autosave_into_adapter(eco: EcologyAdapter, path_npz: str) -> bool:
    """
    Lightweight recreation of run_simulation.load_eco_autosave logic:
      - Load LAI (H,W) and species_weights (S,)
      - Rebuild LAI_layers_SK [S,K,H,W] as weights×(LAI/K)
    """
    try:
        data = np.load(path_npz)
        LAI = np.asarray(data.get("LAI"))
        w = np.asarray(data.get("species_weights"))
        if LAI is None or w is None or LAI.ndim != 2 or w.ndim != 1:
            print(f"[Ecology] Autosave malformed: LAI shape or species_weights missing.", file=sys.stderr)
            return False
        pop = eco.pop
        H, W = LAI.shape
        S = int(w.size)
        K = int(getattr(pop, "K", 1))
        if K < 1:
            K = 1
            pop.K = 1
        # Normalize weights
        ssum = float(np.sum(w))
        ww = (w / ssum) if ssum > 0 else np.full((S,), 1.0 / float(S), dtype=float)
        # Build LAI layers tensor
        LAI_layers_SK = np.zeros((S, K, H, W), dtype=float)
        for s_idx in range(S):
            frac_s = float(ww[s_idx])
            for k in range(K):
                LAI_layers_SK[s_idx, k, :, :] = frac_s * (LAI / float(K))
        pop.LAI_layers_SK = np.clip(LAI_layers_SK, 0.0, float(getattr(pop.params, "lai_max", 5.0)))
        pop.LAI_layers = np.sum(pop.LAI_layers_SK, axis=0)
        pop.LAI = np.sum(pop.LAI_layers, axis=0)
        pop.species_weights = ww
        # Optional diag
        if int(os.getenv("QD_ECO_DIAG", "1")) == 1:
            sdiag = pop.summary()
            print(f"[Ecology] Autosave loaded: LAI(min/mean/max)={sdiag['LAI_min']:.2f}/{sdiag['LAI_mean']:.2f}/{sdiag['LAI_max']:.2f}; S={S}, K={K}")
        return True
    except Exception as e:
        print(f"[Ecology] Autosave load failed: {e}", file=sys.stderr)
        return False


def main():
    _require_matplotlib()

    ap = argparse.ArgumentParser(description="Plot top-3 ecology distributions at a given location.")
    ap.add_argument("--lat", type=float, required=True, help="latitude in degrees")
    ap.add_argument("--lon", type=float, required=True, help="longitude in degrees")
    ap.add_argument("--out", type=str, default="output/top3_ecology.png", help="output PNG path")
    ap.add_argument("--eco-npz", type=str, default="data/eco_autosave.npz", help="path to ecology autosave npz")
    ap.add_argument("--topo", type=str, default=None, help="path or glob to topography NetCDF (optional)")
    ap.add_argument("--nbins", type=int, default=24, help="histogram bins")
    ap.add_argument("--neigh", type=int, default=1, help="neighborhood radius for canopy height (1 -> 3x3)")
    args = ap.parse_args()

    lat_deg = float(args.lat)
    lon_deg = float(args.lon)
    out_png = str(args.out)
    eco_npz = str(args.eco_npz)
    topo_arg = args.topo

    # Resolve topo (support glob)
    topo_path = None
    if topo_arg:
        if any(ch in topo_arg for ch in ["*", "?", "["]):
            topo_path = _latest_file(topo_arg)
            if topo_path is None:
                print(f"[Topo] No files matched pattern '{topo_arg}', proceeding with fallback grid.", file=sys.stderr)
        else:
            topo_path = topo_arg if os.path.exists(topo_arg) else None
            if topo_path is None:
                print(f"[Topo] File '{topo_arg}' not found, proceeding with fallback grid.", file=sys.stderr)
    else:
        topo_path = _latest_file("data/topography_*.nc")

    # Peek dims from eco autosave (if present) to align grid with LAI arrays
    target_hw = _peek_eco_npz_dims(eco_npz) if (eco_npz and os.path.exists(eco_npz)) else None
    if target_hw:
        print(f"[Ecology] Using LAI dims from autosave: H={target_hw[0]} W={target_hw[1]}")
    # Build grid/mask with alignment preference
    grid, land_mask = _load_grid_and_mask(topo_path, target_hw=target_hw)

    # Initialize EcologyAdapter (creates PopulationManager and bands/genes)
    eco = EcologyAdapter(grid, land_mask)

    # Optionally load autosave LAI/species weights (now grid is aligned to LAI if provided)
    if eco_npz and os.path.exists(eco_npz):
        ok = _load_eco_autosave_into_adapter(eco, eco_npz)
        if not ok:
            print("[Ecology] Proceeding without autosave (using adapter defaults).", file=sys.stderr)
    else:
        print(f"[Ecology] Autosave '{eco_npz}' not found. Using adapter defaults.", file=sys.stderr)

    # IndividualPool not reconstructed offline
    indiv = None

    # Plot
    title = f"Top-3 species at ({lat_deg:.2f}°, {lon_deg:.2f}°)"
    try:
        fig, _ = plot_top3_species_distributions(
            eco, grid, indiv,
            lat_deg=lat_deg,
            lon_deg=lon_deg,
            nbins=int(args.nbins),
            neigh_radius=int(args.neigh),
            save_path=out_png,
            title=title
        )
        print(f"[Plot] Saved 12-panel figure to '{out_png}'")
    except Exception as e:
        print(f"[Error] Plotting failed: {e}", file=sys.stderr)
        sys.exit(2)


if __name__ == "__main__":
    main()

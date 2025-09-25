# pygcm/topography.py
"""
Topography and surface properties generation for the Qingdai GCM.

Implements Project 004 (Milestone 1): L1 (large-scale continents) + L3 (fractal detail)
with adaptive sea level to achieve a target land fraction. Also provides base
surface properties and NetCDF export.

Key public functions:
- generate_elevation_map(grid, seed=42, params=None) -> elevation_m
- create_land_sea_mask_from_elevation(elevation_m, grid, target_land_frac=0.29) -> (mask, H_sea_m)
- create_land_sea_mask(grid, target_land_frac=0.29, seed=42, params=None) -> mask      (back-compat shim)
- generate_base_properties(mask, elevation=None, grid=None) -> (base_albedo_map, friction_map)
- export_topography_to_netcdf(grid, elevation, land_mask, base_albedo, friction, sea_level_m, out_path)

Notes:
- All computations assume lat in degrees and lon in degrees on a regular lon-lat grid.
- Land-sea area fractions are computed with cosine(latitude) area weights.
"""

from __future__ import annotations

import os
from typing import Dict, Tuple, Optional

import numpy as np
from scipy.ndimage import gaussian_filter

from . import constants


# ----------------------------
# Utilities
# ----------------------------

def _great_circle_distance_rad(lat_deg: np.ndarray, lon_deg: np.ndarray,
                               lat0_deg: float, lon0_deg: float) -> np.ndarray:
    """
    Great-circle angular distance (in radians) between grid (lat, lon) and a point (lat0, lon0).
    Uses the spherical law of cosines with numeric safety.
    """
    lat = np.deg2rad(lat_deg)
    lon = np.deg2rad(lon_deg)
    lat0 = np.deg2rad(lat0_deg)
    lon0 = np.deg2rad(lon0_deg)

    # spherical law of cosines
    cos_d = np.sin(lat) * np.sin(lat0) + np.cos(lat) * np.cos(lat0) * np.cos(lon - lon0)
    cos_d = np.clip(cos_d, -1.0, 1.0)
    d = np.arccos(cos_d)
    return d


def _rng(seed: Optional[int]) -> np.random.Generator:
    return np.random.default_rng(None if seed is None else int(seed))


def _weighted_quantile(values: np.ndarray, weights: np.ndarray, q: float) -> float:
    """
    Weighted quantile in [0,1], robust to NaNs (ignores them).
    Returns value v such that weighted fraction >= q is above v.
    """
    v = values.ravel()
    w = weights.ravel()
    mask = np.isfinite(v) & np.isfinite(w)
    v = v[mask]
    w = w[mask]
    if v.size == 0:
        return np.nan
    sorter = np.argsort(v)
    v = v[sorter]
    w = w[sorter]
    cum_w = np.cumsum(w)
    if cum_w[-1] <= 0:
        return np.nan
    # normalize 0..1
    cum_w /= cum_w[-1]
    # We want threshold t so that fraction of values >= t equals target_land_frac.
    # If land is defined as v >= t, then weighted fraction of (v >= t) = target_land_frac.
    # Equivalent to choose t as the (1 - target_land_frac) quantile from the bottom.
    idx = np.searchsorted(cum_w, q, side="left")
    idx = np.clip(idx, 0, v.size - 1)
    return float(v[idx])


# ----------------------------
# L1 + L3 Elevation Generation
# ----------------------------

def _generate_L1_continents(grid, seed: int, params: Dict) -> np.ndarray:
    """
    Large-scale continents via sum of Gaussian hills centered at random seeds on the sphere.
    Number of continents controlled by N_CONTINENTS. Width controlled by CONTINENT_SIGMA_DEG.
    Also blends with a very low-frequency band-limited random field for natural irregularity.
    """
    lat_mesh = grid.lat_mesh
    lon_mesh = grid.lon_mesh
    n_lat, n_lon = lat_mesh.shape
    rng = _rng(seed)

    N_CONT = int(params.get("N_CONTINENTS", 3))
    SIGMA_DEG = float(params.get("CONTINENT_SIGMA_DEG", 30.0))  # ~ 3300 km half-width
    SHAPE_P = float(params.get("CONTINENT_SHAPE_P", 2.0))       # generalized Gaussian exponent
    A_MIN, A_MAX = params.get("CONTINENT_AMP_RANGE", (0.8, 1.2))

    # Sample seed centers with optional minimum great-circle spacing.
    # Default: lon ~ U(0, 360), sin(lat) ~ U(-1, 1) (area-uniform).
    # If CONT_MIN_DIST_DEG > 0, enforce a Poisson-disk-like minimum separation.
    MIN_DIST_DEG = float(params.get("CONT_MIN_DIST_DEG", 0.0))
    if MIN_DIST_DEG <= 0.0:
        cont_lats = np.rad2deg(np.arcsin(rng.uniform(-1.0, 1.0, size=N_CONT)))
        cont_lons = rng.uniform(0.0, 360.0, size=N_CONT)
    else:
        cont_lats_list = []
        cont_lons_list = []
        max_tries = 10000
        tries = 0
        while len(cont_lats_list) < N_CONT and tries < max_tries:
            # Candidate sampled area-uniform in latitude
            lat_cand = np.rad2deg(np.arcsin(rng.uniform(-1.0, 1.0)))
            lon_cand = rng.uniform(0.0, 360.0)
            ok = True
            for la, lo in zip(cont_lats_list, cont_lons_list):
                # Great-circle distance between (lat_cand, lon_cand) and (la, lo)
                lat1 = np.deg2rad(lat_cand); lon1 = np.deg2rad(lon_cand)
                lat2 = np.deg2rad(la);       lon2 = np.deg2rad(lo)
                cos_d = np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon1 - lon2)
                cos_d = np.clip(cos_d, -1.0, 1.0)
                d_deg = np.rad2deg(np.arccos(cos_d))
                if d_deg < MIN_DIST_DEG:
                    ok = False
                    break
            if ok:
                cont_lats_list.append(lat_cand)
                cont_lons_list.append(lon_cand)
            tries += 1
        if len(cont_lats_list) < N_CONT:
            # Fallback: fill the remainder without spacing if not enough found
            n_rem = N_CONT - len(cont_lats_list)
            extra_lats = np.rad2deg(np.arcsin(rng.uniform(-1.0, 1.0, size=n_rem)))
            extra_lons = rng.uniform(0.0, 360.0, size=n_rem)
            cont_lats_list.extend(list(extra_lats))
            cont_lons_list.extend(list(extra_lons))
        cont_lats = np.asarray(cont_lats_list[:N_CONT])
        cont_lons = np.asarray(cont_lons_list[:N_CONT])
    cont_amps = rng.uniform(A_MIN, A_MAX, size=N_CONT)

    H_l1 = np.zeros_like(lat_mesh, dtype=float)
    for lat0, lon0, A in zip(cont_lats, cont_lons, cont_amps):
        d = _great_circle_distance_rad(lat_mesh, lon_mesh, lat0, lon0)  # radians
        sigma_rad = np.deg2rad(SIGMA_DEG)
        # generalized Gaussian bump on the sphere
        bump = A * np.exp(- (d / sigma_rad) ** SHAPE_P)
        H_l1 += bump

    # Normalize to zero mean, unit std (avoid degenerate)
    H_l1 = (H_l1 - np.mean(H_l1)) / (np.std(H_l1) + 1e-8)

    # Blend with very-low-frequency noise for irregularity
    VLF_SIGMA_LAT = float(params.get("VLF_SIGMA_LAT", max(4, n_lat // 12)))
    VLF_SIGMA_LON = float(params.get("VLF_SIGMA_LON", max(8, n_lon // 12)))
    noise = rng.standard_normal(size=(n_lat, n_lon))
    # Apply anisotropic smoothing; emulate via two passes
    vlf = gaussian_filter(noise, sigma=(VLF_SIGMA_LAT, VLF_SIGMA_LON), mode=("nearest", "wrap"))
    vlf = (vlf - vlf.mean()) / (vlf.std() + 1e-8)

    W_VLF = float(params.get("W_VLF", 0.35))
    H_l1 = (1 - W_VLF) * H_l1 + W_VLF * vlf
    # Final normalize
    H_l1 = (H_l1 - np.mean(H_l1)) / (np.std(H_l1) + 1e-8)
    return H_l1


def _generate_L3_fbm(grid, seed: int, params: Dict) -> np.ndarray:
    """
    Fractal Brownian Motion (fBm)-like texture using Gaussian-filtered octaves
    of white noise at decreasing smoothing scales (increasing spatial frequency).
    """
    lat_mesh = grid.lat_mesh
    n_lat, n_lon = lat_mesh.shape
    rng = _rng(seed)

    OCT = int(params.get("FBM_OCTAVES", 5))
    HURST = float(params.get("HURST_H", 0.8))  # amplitude decay per octave ~ 2^{-H}
    BASE_SIGMA_LAT = float(params.get("FBM_BASE_SIGMA_LAT", max(1, n_lat // 20)))
    BASE_SIGMA_LON = float(params.get("FBM_BASE_SIGMA_LON", max(1, n_lon // 20)))

    fbm = np.zeros((n_lat, n_lon), dtype=float)
    amp = 1.0
    sigma_lat = BASE_SIGMA_LAT
    sigma_lon = BASE_SIGMA_LON
    for _ in range(OCT):
        noise = rng.standard_normal(size=(n_lat, n_lon))
        layer = gaussian_filter(noise, sigma=(sigma_lat, sigma_lon), mode=("nearest", "wrap"))
        layer = (layer - layer.mean()) / (layer.std() + 1e-8)
        fbm += amp * layer
        amp *= 2 ** (-HURST)
        # next octave: smaller smoothing (higher freq)
        sigma_lat = max(0.5, sigma_lat / 2.0)
        sigma_lon = max(0.5, sigma_lon / 2.0)

    fbm = (fbm - fbm.mean()) / (fbm.std() + 1e-8)
    return fbm


def generate_elevation_map(grid, seed: int = 42, params: Optional[Dict] = None) -> np.ndarray:
    """
    Generate elevation (in meters) combining L1 continents and L3 fractal detail.

    Args:
        grid (SphericalGrid): Grid with lat_mesh, lon_mesh (deg).
        seed (int): RNG seed for reproducibility.
        params (dict): Optional parameter overrides:
            - N_CONTINENTS (default 3)
            - CONTINENT_SIGMA_DEG (default 30)
            - CONTINENT_SHAPE_P (default 2)
            - CONTINENT_AMP_RANGE (default (0.8,1.2))
            - W_VLF (default 0.35)
            - FBM_OCTAVES (default 5)
            - HURST_H (default 0.8)
            - FBM_BASE_SIGMA_LAT/LON (defaults from grid size)
            - W1, W3 weights (defaults 1.0 and 0.6)
            - SCALE_M (default 4500) total rough vertical scale (pre-sea-level)

    Returns:
        elevation_m (np.ndarray): 2D array, meters, positive above reference, negative below.
    """
    if params is None:
        params = {}
    seed = int(seed)

    H_l1 = _generate_L1_continents(grid, seed=seed, params=params)
    H_l3 = _generate_L3_fbm(grid, seed=seed + 1, params=params)

    W1 = float(params.get("W1", 1.0))
    W3 = float(params.get("W3", 0.6))
    combined = W1 * H_l1 + W3 * H_l3

    # Normalize and scale to meters (pre sea-level threshold)
    combined = (combined - combined.mean()) / (combined.std() + 1e-8)
    SCALE_M = float(params.get("SCALE_M", 4500.0))
    elevation_m = combined * SCALE_M

    # Gentle smoothing to remove pixel noise (does not remove structures)
    elevation_m = gaussian_filter(elevation_m, sigma=(0.5, 0.5), mode=("nearest", "wrap"))
    return elevation_m


# ----------------------------
# Sea level and land-sea mask
# ----------------------------

def create_land_sea_mask_from_elevation(elevation_m: np.ndarray,
                                        grid,
                                        target_land_frac: float = 0.29) -> Tuple[np.ndarray, float]:
    """
    Compute sea level by weighted quantile so that area fraction of land (H >= H_sea)
    matches target_land_frac. Returns binary mask (1=land, 0=ocean) and sea level (m).

    Uses cosine(latitude) area weights.
    """
    lat_rad = np.deg2rad(grid.lat_mesh)
    area_w = np.cos(lat_rad)
    area_w = np.maximum(area_w, 0.0)

    # Choose threshold so that weighted fraction of values >= threshold equals target_land_frac
    # That means threshold is the (1 - target_land_frac) quantile from the bottom.
    q = 1.0 - float(target_land_frac)
    H_sea = _weighted_quantile(elevation_m, area_w, q=q)

    mask = (elevation_m >= H_sea).astype(np.uint8)

    # Report achieved fraction
    land_frac = float((area_w * (mask == 1)).sum() / (area_w.sum() + 1e-15))
    print(f"[Topography] Target land fraction={target_land_frac:.3f}, achieved={land_frac:.3f}, sea_level={H_sea:.1f} m")
    return mask, float(H_sea)


def create_land_sea_mask(grid, target_land_frac: float = 0.29, seed: int = 42, params: Optional[Dict] = None) -> np.ndarray:
    """
    Back-compat shim: generate a natural land-sea mask by creating elevation first,
    then applying adaptive sea level to hit the target land fraction.

    This replaces the previous rectangular continents.
    """
    elevation = generate_elevation_map(grid, seed=seed, params=params)
    mask, _ = create_land_sea_mask_from_elevation(elevation, grid, target_land_frac=target_land_frac)
    return mask


# ----------------------------
# Base surface properties
# ----------------------------

def generate_base_properties(mask: np.ndarray,
                             elevation: Optional[np.ndarray] = None,
                             grid=None) -> Tuple[np.ndarray, np.ndarray]:
    """
    Generate ice-free base albedo and surface friction maps.

    Albedo:
      - Ocean baseline low (~0.08-0.10)
      - Land baseline higher (~0.28), slightly increases with |lat| and elevation.

    Friction:
      - Ocean weak friction (1e-6)
      - Land stronger (1e-5), increases with elevation.

    Args:
        mask: 2D uint8 array (1=land, 0=ocean).
        elevation: Optional elevation (m); if None, treated as zeros for modulation.
        grid: Optional SphericalGrid for latitude-dependent terms.

    Returns:
        (albedo_map, friction_map)
    """
    mask = mask.astype(np.uint8)
    if elevation is None:
        elevation = np.zeros_like(mask, dtype=float)

    if grid is not None:
        lat_abs = np.abs(grid.lat_mesh)
        lat_factor = (lat_abs / 90.0) ** 2  # 0..1, higher at poles
    else:
        lat_factor = np.zeros_like(mask, dtype=float)

    # Base albedo
    albedo_ocean = 0.08
    albedo_land = 0.28

    # Elevation effect (only land)
    elev_norm = np.clip(np.maximum(elevation, 0.0) / 4000.0, 0.0, 1.0)
    albedo = np.where(mask == 1, albedo_land, albedo_ocean)
    albedo += 0.08 * lat_factor  # both land & ocean slightly brighter toward poles
    albedo += 0.05 * elev_norm * (mask == 1)  # land elevation brightening
    albedo = np.clip(albedo, 0.05, 0.85)

    # Friction (Rayleigh-type coefficient for lower boundary)
    friction_ocean = 1.0e-6
    friction_land = 1.0e-5
    friction = np.where(mask == 1, friction_land, friction_ocean)
    friction += 6.0e-6 * elev_norm * (mask == 1)  # mountain drag
    # Clip to safe bounds
    friction = np.clip(friction, 5e-7, 3e-5)

    return albedo, friction


# ----------------------------
# NetCDF Export
# ----------------------------

def export_topography_to_netcdf(grid,
                                elevation: np.ndarray,
                                land_mask: np.ndarray,
                                base_albedo: np.ndarray,
                                friction: np.ndarray,
                                sea_level_m: float,
                                out_path: str) -> None:
    """
    Export multiple fields to NetCDF:
      - coords: lat (deg), lon (deg)
      - variables: elevation (m), land_mask (0/1), base_albedo (1), friction (s^-1 approx)
      - global attributes include planet constants and sea level
    """
    # Lazy import to keep core runnable without netCDF4 present
    try:
        from netCDF4 import Dataset
    except Exception as e:
        raise RuntimeError("netCDF4 is required for export. Please install 'netCDF4' and retry.") from e

    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    n_lat = grid.n_lat
    n_lon = grid.n_lon

    with Dataset(out_path, "w") as ds:
        # Dimensions
        ds.createDimension("lat", n_lat)
        ds.createDimension("lon", n_lon)

        # Coordinates
        lat_var = ds.createVariable("lat", "f4", ("lat",))
        lon_var = ds.createVariable("lon", "f4", ("lon",))
        lat_var.units = "degrees_north"
        lon_var.units = "degrees_east"
        lat_var[:] = grid.lat
        lon_var[:] = grid.lon

        # Variables
        elev_var = ds.createVariable("elevation", "f4", ("lat", "lon"), zlib=True, complevel=4)
        mask_var = ds.createVariable("land_mask", "i1", ("lat", "lon"), zlib=True, complevel=4)
        alb_var = ds.createVariable("base_albedo", "f4", ("lat", "lon"), zlib=True, complevel=4)
        fric_var = ds.createVariable("friction", "f4", ("lat", "lon"), zlib=True, complevel=4)

        elev_var.units = "m"
        elev_var.long_name = "surface_elevation_above_reference_sea_level"
        mask_var.long_name = "land_sea_mask"
        mask_var.flag_values = "0, 1"
        mask_var.flag_meanings = "ocean land"
        alb_var.units = "1"
        alb_var.long_name = "base_surface_albedo_ice_free"
        fric_var.units = "s-1"
        fric_var.long_name = "surface_friction_coefficient"

        elev_var[:, :] = elevation.astype(np.float32)
        mask_var[:, :] = land_mask.astype(np.int8)
        alb_var[:, :] = base_albedo.astype(np.float32)
        fric_var[:, :] = friction.astype(np.float32)

        # Global attributes
        ds.title = "Qingdai Topography and Surface Properties"
        ds.institution = "PyGCM for Qingdai"
        ds.source = "Procedural generation (Project 004 Milestone 1)"
        ds.history = "Created by pygcm.topography.export_topography_to_netcdf"
        ds.references = "docs/projects/004-topography-generation.md"
        ds.sea_level_m = float(sea_level_m)
        ds.target_land_fraction = 0.29

        # Planet constants
        ds.planet_radius_m = constants.PLANET_RADIUS
        ds.planet_omega_rad_s = constants.PLANET_OMEGA
        ds.planet_axial_tilt_deg = constants.PLANET_AXIAL_TILT

# ----------------------------
# NetCDF Load & Regrid
# ----------------------------
def load_topography_from_netcdf(path: str, grid, *, regrid: str = "auto"):
    """
    Load topography and surface properties from a NetCDF file and (optionally) regrid to the model grid.

    Args:
        path (str): Path to NetCDF file produced by export_topography_to_netcdf or compatible layout.
        grid (SphericalGrid): Target grid to map fields onto.
        regrid (str): "auto" (default) attempts bilinear (nearest for mask) if grids differ; "never" to require exact match.

    Returns:
        (elevation, land_mask, base_albedo, friction): 2D arrays aligned to 'grid' resolution.
    """
    try:
        from netCDF4 import Dataset
    except Exception as e:
        raise RuntimeError("netCDF4 is required to load topography. Please install 'netCDF4'.") from e

    import numpy as np
    try:
        from scipy.interpolate import RegularGridInterpolator
    except Exception as e:
        raise RuntimeError("SciPy is required for regridding. Please install 'scipy'.") from e

    def _to_0360(lon_arr):
        lon = np.asarray(lon_arr, dtype=float).copy()
        lon = np.mod(lon, 360.0)
        lon[lon < 0] += 360.0
        return lon

    def _prepare_src_coords(ds):
        lat = np.asarray(ds["lat"][:], dtype=float)
        lon = np.asarray(ds["lon"][:], dtype=float)

        # Normalize longitude to [0,360)
        if np.nanmin(lon) < 0.0 or np.nanmax(lon) <= 180.0:
            lon = _to_0360(lon)

        # Ensure strictly increasing coordinates for interpolator
        lat_increasing = np.all(np.diff(lat) > 0)
        if not lat_increasing:
            lat = lat[::-1]
        # Sort longitude ascending and roll fields accordingly later
        lon_sort_idx = np.argsort(lon)
        lon = lon[lon_sort_idx]
        return lat, lon, lat_increasing, lon_sort_idx

    def _read_field(ds, name, lat_increasing, lon_sort_idx):
        arr = np.asarray(ds[name][:])
        # Reorder dims if necessary to (lat, lon)
        # Convention here is already (lat, lon) as written by our exporter.
        # Fix latitude descending
        if not lat_increasing:
            arr = arr[::-1, :]
        # Sort longitude
        arr = arr[:, lon_sort_idx]
        return arr

    def _interp_field(src_lat, src_lon, src_field, tgt_lat_mesh, tgt_lon_mesh, method="linear", is_mask=False):
        # Build cyclic extension in longitude to avoid seam artifacts
        lon_ext = np.concatenate([src_lon - 360.0, src_lon, src_lon + 360.0])
        field_ext = np.concatenate([src_field, src_field, src_field], axis=1)

        # Interpolator requires strictly increasing coords
        interp = RegularGridInterpolator(
            (src_lat, lon_ext), field_ext,
            bounds_error=False,
            fill_value=None,
            method=("nearest" if is_mask else method)
        )

        # Query points
        pts_lat = tgt_lat_mesh.ravel()
        pts_lon = tgt_lon_mesh.ravel()
        # Clip latitude to source bounds to avoid None
        pts_lat = np.clip(pts_lat, src_lat.min(), src_lat.max())

        vals = interp(np.stack([pts_lat, pts_lon], axis=-1)).reshape(tgt_lat_mesh.shape)

        if is_mask:
            # Force binary after nearest
            vals = np.where(vals >= 0.5, 1, 0).astype(np.uint8)
        else:
            # If any NaNs occurred (e.g., outside bounds), fall back to nearest for those pixels
            if np.any(~np.isfinite(vals)):
                interp_nn = RegularGridInterpolator(
                    (src_lat, lon_ext), field_ext,
                    bounds_error=False,
                    fill_value=None,
                    method="nearest"
                )
                nn_vals = interp_nn(np.stack([pts_lat, pts_lon], axis=-1)).reshape(tgt_lat_mesh.shape)
                vals = np.where(np.isfinite(vals), vals, nn_vals)

        return vals

    with Dataset(path, "r") as ds:
        src_lat, src_lon, lat_increasing, lon_sort_idx = _prepare_src_coords(ds)

        elev = _read_field(ds, "elevation", lat_increasing, lon_sort_idx)
        mask = _read_field(ds, "land_mask", lat_increasing, lon_sort_idx)
        base = _read_field(ds, "base_albedo", lat_increasing, lon_sort_idx)
        fric = _read_field(ds, "friction", lat_increasing, lon_sort_idx)

        # Remove duplicate seam if present (e.g., lon includes both 0 and 360)
        if src_lon.size >= 2 and (np.isclose(src_lon[0], 0.0) and np.isclose(src_lon[-1], 360.0)):
            src_lon = src_lon[:-1]
            elev = elev[:, :-1]
            mask = mask[:, :-1]
            base = base[:, :-1]
            fric = fric[:, :-1]

        # Quick exact-shape match path
        same_shape = (elev.shape == (grid.n_lat, grid.n_lon))
        if same_shape:
            # Also check coords close
            if regrid == "never" or (
                np.allclose(src_lat, grid.lat, atol=1e-6) and np.allclose(src_lon, grid.lon, atol=1e-6)
            ):
                elevation = elev.astype(float)
                land_mask = mask.astype(np.uint8)
                base_albedo = base.astype(float)
                friction = fric.astype(float)
            else:
                # Shapes match but coords differ slightly; treat as regrid
                elevation = _interp_field(src_lat, src_lon, elev, grid.lat_mesh, grid.lon_mesh, method="linear", is_mask=False)
                land_mask = _interp_field(src_lat, src_lon, mask, grid.lat_mesh, grid.lon_mesh, method="nearest", is_mask=True)
                base_albedo = _interp_field(src_lat, src_lon, base, grid.lat_mesh, grid.lon_mesh, method="linear", is_mask=False)
                friction = _interp_field(src_lat, src_lon, fric, grid.lat_mesh, grid.lon_mesh, method="linear", is_mask=False)
        else:
            if regrid == "never":
                raise ValueError(f"Topography grid mismatch: source {elev.shape} vs target {(grid.n_lat, grid.n_lon)} and regrid='never'.")
            elevation = _interp_field(src_lat, src_lon, elev, grid.lat_mesh, grid.lon_mesh, method="linear", is_mask=False)
            land_mask = _interp_field(src_lat, src_lon, mask, grid.lat_mesh, grid.lon_mesh, method="nearest", is_mask=True)
            base_albedo = _interp_field(src_lat, src_lon, base, grid.lat_mesh, grid.lon_mesh, method="linear", is_mask=False)
            friction = _interp_field(src_lat, src_lon, fric, grid.lat_mesh, grid.lon_mesh, method="linear", is_mask=False)

    # Basic sanity/logging
    lat_rad = np.deg2rad(grid.lat_mesh)
    area_w = np.cos(lat_rad)
    achieved = float((area_w * (land_mask == 1)).sum() / (area_w.sum() + 1e-15))
    print(f"[Topo] Loaded: {path}")
    print(f"[Topo] Land fraction (achieved): {achieved:.3f}")
    print(f"[Topo] Albedo stats (min/mean/max): {np.nanmin(base_albedo):.3f}/{np.nanmean(base_albedo):.3f}/{np.nanmax(base_albedo):.3f}")
    print(f"[Topo] Friction stats (min/mean/max): {np.nanmin(friction):.2e}/{np.nanmean(friction):.2e}/{np.nanmax(friction):.2e}")
    if np.isfinite(elevation).any():
        print(f"[Topo] Elevation stats (m): {np.nanmin(elevation):.1f}/{np.nanmean(elevation):.1f}/{np.nanmax(elevation):.1f}")

    return elevation, land_mask, base_albedo, friction

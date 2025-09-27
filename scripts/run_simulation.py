# scripts/run_simulation.py

"""
Main simulation script for the Qingdai GCM.
"""

import numpy as np
import sys
import os
import matplotlib.pyplot as plt
import signal
import atexit
import json

# Add the project root to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcm.constants as constants

from pygcm.grid import SphericalGrid
from pygcm.orbital import OrbitalSystem
from pygcm.forcing import ThermalForcing
from pygcm.dynamics import SpectralModel
from pygcm.topography import create_land_sea_mask, generate_base_properties, load_topography_from_netcdf
from pygcm.physics import diagnose_precipitation, diagnose_precipitation_hybrid, parameterize_cloud_cover, calculate_dynamic_albedo, cloud_from_precip, compute_orographic_factor
from pygcm.hydrology import (
    get_hydrology_params_from_env,
    partition_precip_phase,
    partition_precip_phase_smooth,
    snow_step,
    snowpack_step,
    update_land_bucket,
    diagnose_water_closure,
)
from pygcm.routing import RiverRouting
from pygcm import energy as energy
from pygcm.ocean import WindDrivenSlabOcean
from pygcm.jax_compat import is_enabled as JAX_IS_ENABLED
# P020 Phase 0: OO façade (optional)
try:
    from pygcm.world import QingdaiWorld  # Phase 0 skeleton
except Exception:
    QingdaiWorld = None
# Ecology (P015 M1) adapter + Phyto (P017)
try:
    from pygcm.ecology import EcologyAdapter, PhytoManager
except Exception:
    EcologyAdapter = None
    PhytoManager = None
# Diversity diagnostics (alpha/beta)
try:
    from pygcm.ecology import diversity as eco_diversity
except Exception:
    eco_diversity = None
# Vectorized individual pool (subdaily spectral adaptation, sampled cells)
try:
    from pygcm.ecology.individuals import IndividualPool
except Exception:
    IndividualPool = None

# --- Restart I/O (NetCDF) and Initialization Helpers ---

def save_restart(path, grid, gcm, ocean, land_mask, W_land=None, S_snow=None, C_snow=None, t_seconds: float | None = None):
    """
    Save minimal prognostic state to a NetCDF restart file.
    Includes: lat/lon coords, u/v/h, T_s, cloud_cover, q (if exists), h_ice (if exists),
              ocean state (uo/vo/eta/Ts if ocean provided),
              hydrology reservoirs (W_land, S_snow) if provided,
              snow optical coverage C_snow (optional; derived from SWE, persisted for viz continuity).
    """
    from netCDF4 import Dataset
    import numpy as np
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with Dataset(path, "w") as ds:
        nlat, nlon = grid.n_lat, grid.n_lon
        ds.createDimension("lat", nlat)
        ds.createDimension("lon", nlon)
        vlat = ds.createVariable("lat", "f4", ("lat",))
        vlon = ds.createVariable("lon", "f4", ("lon",))
        vlat[:] = grid.lat
        vlon[:] = grid.lon

        def wvar(name, data):
            if data is None:
                return
            var = ds.createVariable(name, "f4", ("lat", "lon"))
            var[:] = np.asarray(data, dtype=np.float32)

        # Atmospheric / surface
        wvar("u", gcm.u)
        wvar("v", gcm.v)
        wvar("h", gcm.h)
        wvar("T_s", gcm.T_s)
        wvar("cloud_cover", getattr(gcm, "cloud_cover", None))
        wvar("q", getattr(gcm, "q", None))
        wvar("h_ice", getattr(gcm, "h_ice", None))

        # Ocean
        if ocean is not None:
            wvar("uo", getattr(ocean, "uo", None))
            wvar("vo", getattr(ocean, "vo", None))
            wvar("eta", getattr(ocean, "eta", None))
            wvar("Ts", getattr(ocean, "Ts", None))

        # Hydrology / Cryo
        wvar("W_land", W_land)
        wvar("S_snow", S_snow)
        wvar("C_snow", C_snow)

        # Masks for reference
        wvar("land_mask", land_mask)

        # Optional: astronomical epoch (simulation time in seconds)
        try:
            vts = ds.createVariable("t_seconds", "f8")
            vts[...] = float(t_seconds) if (t_seconds is not None) else 0.0
        except Exception:
            pass

        # Minimal metadata
        ds.setncattr("title", "Qingdai GCM Restart")
        ds.setncattr("creator", "PyGCM for Qingdai")
        ds.setncattr("note", "Contains minimal prognostic fields for warm restart (incl. t_seconds).")
        ds.setncattr("format", "v1")

def save_topography(path, grid, land_mask, base_albedo_map, friction_map, elevation=None):
    """
    Write standardized topography file:
      - path: data/topography.nc
      - variables: lat, lon, land_mask (u1), base_albedo (f4), friction (f4), optional elevation (f4)
    """
    from netCDF4 import Dataset
    import numpy as _np
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with Dataset(path, "w") as ds:
        nlat, nlon = grid.n_lat, grid.n_lon
        ds.createDimension("lat", nlat)
        ds.createDimension("lon", nlon)
        vlat = ds.createVariable("lat", "f4", ("lat",))
        vlon = ds.createVariable("lon", "f4", ("lon",))
        vlat[:] = grid.lat.astype(_np.float32)
        vlon[:] = grid.lon.astype(_np.float32)

        vmask = ds.createVariable("land_mask", "u1", ("lat", "lon"))
        vmask[:] = _np.asarray(land_mask, dtype=_np.uint8)

        vba = ds.createVariable("base_albedo", "f4", ("lat", "lon"))
        vba[:] = _np.asarray(base_albedo_map, dtype=_np.float32)

        vfr = ds.createVariable("friction", "f4", ("lat", "lon"))
        vfr[:] = _np.asarray(friction_map, dtype=_np.float32)

        if elevation is not None:
            vel = ds.createVariable("elevation", "f4", ("lat", "lon"))
            vel[:] = _np.asarray(elevation, dtype=_np.float32)

        ds.setncattr("title", "Qingdai Topography")
        ds.setncattr("source", "scripts/run_simulation.py")
        ds.setncattr("format", "v1")

def load_restart(path):
    """
    Load restart file and return a dict of arrays. Missing variables are returned as None.
    """
    from netCDF4 import Dataset
    out = {}
    with Dataset(path, "r") as ds:
        def rvar(name):
            try:
                return ds.variables[name][:].data
            except Exception:
                return None
        out["lat"] = ds.variables["lat"][:].data
        out["lon"] = ds.variables["lon"][:].data
        for name in ["u", "v", "h", "T_s", "cloud_cover", "q", "h_ice",
                     "uo", "vo", "eta", "Ts", "W_land", "S_snow", "C_snow", "land_mask"]:
            out[name] = rvar(name)
        # Optional astronomical epoch
        try:
            out["t_seconds"] = float(ds.variables["t_seconds"][...])
        except Exception:
            out["t_seconds"] = None
    return out

def save_ocean(path: str, grid, ocean, day_value: float | None = None) -> bool:
    """
    Write standardized ocean physical state to data/ocean.nc:
      - lat, lon
      - uo, vo (m/s), eta (m), Ts (K)
      - attributes: day (planetary), title/source
    """
    try:
        from netCDF4 import Dataset
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with Dataset(path, "w") as ds:
            nlat, nlon = grid.n_lat, grid.n_lon
            ds.createDimension("lat", nlat)
            ds.createDimension("lon", nlon)
            vlat = ds.createVariable("lat", "f4", ("lat",))
            vlon = ds.createVariable("lon", "f4", ("lon",))
            vlat[:] = grid.lat
            vlon[:] = grid.lon

            def wvar(name, data):
                var = ds.createVariable(name, "f4", ("lat", "lon"))
                var[:] = np.asarray(data, dtype=np.float32)

            if getattr(ocean, "uo", None) is not None: wvar("uo", ocean.uo)
            if getattr(ocean, "vo", None) is not None: wvar("vo", ocean.vo)
            if getattr(ocean, "eta", None) is not None: wvar("eta", ocean.eta)
            if getattr(ocean, "Ts", None) is not None: wvar("Ts", ocean.Ts)

            ds.setncattr("title", "Qingdai Ocean State")
            ds.setncattr("source", "scripts/run_simulation.py")
            if day_value is not None:
                ds.setncattr("day", float(day_value))
        return True
    except Exception as e:
        print(f"[Ocean] Save failed: {e}")
        return False

def load_ocean(path: str) -> dict:
    """
    Load ocean physical state from ocean.nc. Returns dict with keys uo, vo, eta, Ts, day(optional).
    Missing variables are None.
    """
    out = {"uo": None, "vo": None, "eta": None, "Ts": None, "day": None}
    try:
        from netCDF4 import Dataset
        with Dataset(path, "r") as ds:
            def r(name):
                try:
                    return ds.variables[name][:].data
                except Exception:
                    return None
            out["uo"] = r("uo")
            out["vo"] = r("vo")
            out["eta"] = r("eta")
            out["Ts"] = r("Ts")
            try:
                out["day"] = float(ds.getncattr("day"))
            except Exception:
                out["day"] = None
    except Exception as e:
        print(f"[Ocean] Load failed '{path}': {e}")
    return out

def save_autosave(data_dir: str, grid, gcm, ocean, land_mask, W_land, S_snow, eco, day_value: float) -> None:
    """
    Save autosave checkpoint into data/:
      - data/atmosphere.nc   (core model state via NetCDF, includes t_seconds epoch)
      - data/ecology.nc      (ecology extended state via NetCDF: LAI/species_weights/species reflectance/bands)
      - data/genes.json      (species genes snapshot with band weights)
    """
    try:
        os.makedirs(data_dir, exist_ok=True)
    except Exception:
        pass
    # Core model state
    try:
        autosave_nc = os.path.join(data_dir, "atmosphere.nc")
        # Convert current planetary day to seconds for astronomical epoch persistence
        try:
            t_sec = float(day_value) * (2 * np.pi / constants.PLANET_OMEGA)
        except Exception:
            t_sec = 0.0
        save_restart(
            autosave_nc, grid, gcm, ocean, land_mask,
            W_land=W_land, S_snow=S_snow, C_snow=getattr(gcm, "C_snow_map_last", None), t_seconds=t_sec
        )
        print(f"[Autosave] Core state saved to '{autosave_nc}' (standardized atmosphere.nc)")
    except Exception as e:
        print(f"[Autosave] NetCDF save failed: {e}")
    # Ecology (extended if possible, fallback to legacy)
    try:
        if eco is not None and getattr(eco, "pop", None) is not None:
            # Resolve autosave path (allow override by env, enforce NetCDF target)
            path_env = os.getenv("QD_ECO_AUTOSAVE_PATH")
            path_ec = path_env if (path_env and path_env.lower().endswith(".nc")) else os.path.join(data_dir, "ecology.nc")
            # Ensure directory exists for custom path
            try:
                os.makedirs(os.path.dirname(path_ec) or ".", exist_ok=True)
            except Exception:
                pass
            # Write extended ecology state via adapter (NetCDF)
            try:
                if hasattr(eco, "save_autosave"):
                    ok = bool(eco.save_autosave(path_ec, day_value=float(day_value)))
                    if not ok:
                        print(f"[Autosave] Ecology extended save returned False for '{path_ec}'")
                else:
                    print("[Autosave] Ecology adapter has no save_autosave; skipping ecology.nc")
                    ok = False
            except Exception as _esa:
                print(f"[Autosave] Ecology extended save failed: {_esa}")
                ok = False
            # Always attempt to write genes.json alongside ecology state (best-effort)
            try:
                if hasattr(eco, "save_genes_json"):
                    eco.save_genes_json(os.path.join(data_dir, "genes.json"), day_value=float(day_value))
            except Exception:
                pass
    except Exception as e:
        print(f"[Autosave] Ecology save failed: {e}")


# (legacy NPZ loader removed; ecology persistence now uses NetCDF via EcologyAdapter.load_autosave and genes.json)


def apply_banded_initial_ts(grid, gcm, ocean, land_mask):
    """
    Apply latitudinally banded initial surface temperature:
      T(φ) = T_pole + (T_eq - T_pole) * cos^2(φ)
    Controlled by env: QD_INIT_BANDED=1, QD_INIT_T_EQ (K), QD_INIT_T_POLE (K).
    """
    if int(os.getenv("QD_INIT_BANDED", "0")) != 1:
        return
    T_eq = float(os.getenv("QD_INIT_T_EQ", "295.0"))
    T_pole = float(os.getenv("QD_INIT_T_POLE", "265.0"))
    phi = np.deg2rad(grid.lat_mesh)
    Ts0 = T_pole + (T_eq - T_pole) * (np.cos(phi) ** 2)
    # Apply to atmospheric surface temperature
    gcm.T_s = Ts0.copy()
    # If dynamic ocean is enabled, set SST over ocean (preserve land)
    if ocean is not None:
        ocean_mask = (land_mask == 0)
        ocean.Ts = np.where(ocean_mask, Ts0, ocean.Ts)
    print(f"[Init] Applied banded initial Ts: T_eq={T_eq} K, T_pole={T_pole} K")

def plot_state(grid, gcm, land_mask, precip, cloud_cover, albedo, t_days, output_dir, ocean=None, routing=None):
    """
    Generates and saves a 3-column x 5-row diagnostic plot of the current model state.
    Panels (left→right, top→bottom):
      1) Ts (°C), 2) Ta (°C), 3) Sea-level Pressure (hPa)
      4) SST (°C), 5) Precip (1-day, mm/day), 6) Cloud Cover
      7) Wind (streamlines, m/s), 8) Ocean currents (m/s) or h anomaly, 9) Vorticity (1/s)
      10) Incoming Shortwave (W/m²), 11) Dynamic Albedo, 12) OLR (W/m²)
      13) Specific Humidity q (g/kg), 14) Evaporation E (mm/day), 15) Condensation P_cond (mm/day)
    """
    fig, axes = plt.subplots(5, 3, figsize=(22, 28), constrained_layout=True)
    fig.suptitle(f"Qingdai GCM State at Day {t_days:.2f}", fontsize=16)

    g_const = 9.81
    # Temperature diagnostics
    T_a = 288.0 + (g_const / 1004.0) * gcm.h
    ta_c = T_a - 273.15
    ts_c = np.nan_to_num(gcm.T_s - 273.15)
    sst_c = np.nan_to_num((ocean.Ts if ocean is not None else gcm.T_s) - 273.15)
    tmin = float(np.nanmin([ts_c.min(), ta_c.min(), sst_c.min()]))
    tmax = float(np.nanmax([ts_c.max(), ta_c.max(), sst_c.max()]))
    t_levels = np.linspace(tmin, tmax, 20)

    # 1) Ts
    ax = axes[0, 0]
    cs = ax.contourf(grid.lon, grid.lat, ts_c, levels=t_levels, cmap="coolwarm")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)
    ax.set_title("Surface Temperature (°C)")
    fig.colorbar(cs, ax=ax, label="°C")

    # 2) Ta
    ax = axes[0, 1]
    cs = ax.contourf(grid.lon, grid.lat, ta_c, levels=t_levels, cmap="coolwarm")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)
    ax.set_title("Atmospheric Temperature (°C)")
    fig.colorbar(cs, ax=ax, label="°C")

    # 3) Sea-level pressure (hPa) from shallow-water mass (diagnostic)
    p0 = float(getattr(getattr(gcm, "hum_params", None), "p0", 1.0e5))
    rho_air = float(getattr(getattr(gcm, "hum_params", None), "rho_a", 1.2))
    # Interpret h as thickness perturbation. Provide two plotting modes:
    # - anom: pressure anomaly (hPa) = rho*g*h / 100
    # - abs:  absolute pressure (hPa) = (p0 + rho*g*h) / 100
    ps_mode = os.getenv("QD_PLOT_PS_MODE", "anom").lower()
    if ps_mode == "abs":
        ps_field = (p0 + rho_air * g_const * gcm.h) * 1e-2
        title_ps = "Sea-level Pressure (hPa, diag)"
    else:
        ps_field = (rho_air * g_const * gcm.h) * 1e-2
        title_ps = "Sea-level Pressure Anomaly (hPa, diag)"
    ax = axes[0, 2]
    cs = ax.contourf(grid.lon, grid.lat, ps_field, levels=20, cmap="viridis")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)
    ax.set_title(title_ps)
    fig.colorbar(cs, ax=ax, label="hPa")

    # 4) SST
    ax = axes[1, 0]
    cs = ax.contourf(grid.lon, grid.lat, sst_c, levels=t_levels, cmap="coolwarm")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)
    ax.set_title("SST (°C)")
    fig.colorbar(cs, ax=ax, label="°C")

    # 5) Precip (instantaneous, mm/day)
    ax = axes[1, 1]
    precip_mmday = np.nan_to_num(precip) * 86400.0  # kg m^-2 s^-1 → mm/day
    cs = ax.contourf(grid.lon, grid.lat, precip_mmday, levels=np.linspace(0, 30, 11), cmap="Blues", extend="max")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)
    ax.set_title("Precipitation (instant, mm/day)")
    fig.colorbar(cs, ax=ax, label="mm/day")

    # 6) Cloud cover
    ax = axes[1, 2]
    cs = ax.contourf(grid.lon, grid.lat, cloud_cover, levels=np.linspace(0, 1, 11), cmap="Greys")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)
    ax.set_title("Cloud Cover Fraction")
    fig.colorbar(cs, ax=ax, label="Fraction")

    # 7) Wind field (streamlines)
    ax = axes[2, 0]
    speed_w = np.sqrt(np.nan_to_num(gcm.u)**2 + np.nan_to_num(gcm.v)**2)
    strm_w = ax.streamplot(grid.lon, grid.lat, gcm.u, gcm.v, color=speed_w, cmap="viridis", density=1.5)
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)
    ax.set_title("Wind Field (m/s)")
    fig.colorbar(strm_w.lines, ax=ax, label="m/s")

    # 8) Ocean currents or height anomaly
    ax = axes[2, 1]
    if ocean is not None:
        uo = np.nan_to_num(ocean.uo); vo = np.nan_to_num(ocean.vo)
        sp_o = np.sqrt(uo**2 + vo**2)
        strm_o = ax.streamplot(grid.lon, grid.lat, uo, vo, color=sp_o, cmap="viridis", density=1.2)
        ax.set_title("Ocean Currents (m/s)")
        fig.colorbar(strm_o.lines, ax=ax, label="m/s")
    else:
        cs = ax.contourf(grid.lon, grid.lat, gcm.h - float(getattr(gcm, "H", 8000.0)), levels=20, cmap="RdBu_r")
        ax.set_title("Geopotential Height Anomaly (m)")
        fig.colorbar(cs, ax=ax, label="m")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)

    # 9) Vorticity (1/s)
    ax = axes[2, 2]
    vort = grid.vorticity(gcm.u, gcm.v)
    vmax = np.nanmax(np.abs(vort))
    levels = np.linspace(-vmax, vmax, 21) if np.isfinite(vmax) and vmax > 0 else 20
    cs = ax.contourf(grid.lon, grid.lat, vort, levels=levels, cmap="PuOr")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)
    ax.set_title("Relative Vorticity (1/s)")
    fig.colorbar(cs, ax=ax, label="1/s")

    # 10) Incoming Shortwave
    ax = axes[3, 0]
    cs = ax.contourf(grid.lon, grid.lat, gcm.isr, levels=20, cmap="magma")
    try:
        idxA = np.unravel_index(np.argmax(gcm.isr_A), gcm.isr_A.shape)
        idxB = np.unravel_index(np.argmax(gcm.isr_B), gcm.isr_B.shape)
        lonA, latA = grid.lon[idxA[1]], grid.lat[idxA[0]]
        lonB, latB = grid.lon[idxB[1]], grid.lat[idxB[0]]
        ax.scatter([lonA], [latA], c="cyan", s=30, marker="x", label="Star A center")
        ax.scatter([lonB], [latB], c="yellow", s=30, marker="+", label="Star B center")
        ax.legend(loc="upper right", fontsize=8)
    except Exception:
        pass
    ax.set_title("Incoming Shortwave (W/m²)")
    fig.colorbar(cs, ax=ax, label="W/m²")

    # 11) Dynamic Albedo
    ax = axes[3, 1]
    cs = ax.contourf(grid.lon, grid.lat, albedo, levels=np.linspace(0, 0.8, 17), cmap="cividis")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)
    ax.set_title("Dynamic Albedo")
    fig.colorbar(cs, ax=ax, label="Albedo")

    # 12) OLR
    ax = axes[3, 2]
    cs = ax.contourf(grid.lon, grid.lat, gcm.olr, levels=20, cmap="plasma")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="white", linewidths=0.5)
    ax.set_title("Outgoing Longwave (W/m²)")
    fig.colorbar(cs, ax=ax, label="W/m²")

    # 13) Specific humidity q (g/kg)
    ax = axes[4, 0]
    if hasattr(gcm, "q"):
        q_gkg = 1e3 * np.nan_to_num(gcm.q)
        cs = ax.contourf(grid.lon, grid.lat, q_gkg, levels=20, cmap="GnBu")
        ax.set_title("Specific Humidity q (g/kg)")
        fig.colorbar(cs, ax=ax, label="g/kg")
    else:
        ax.text(0.5, 0.5, "q not enabled", ha="center", va="center")
        ax.set_title("Specific Humidity q")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)

    # 14) Evaporation E (mm/day)
    ax = axes[4, 1]
    E = getattr(gcm, "E_flux_last", 0.0)
    if np.isscalar(E):
        E = np.full_like(gcm.T_s, float(E))
    E_mmday = np.nan_to_num(E) * 86400.0
    cs = ax.contourf(grid.lon, grid.lat, E_mmday, levels=20, cmap="YlGn")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)
    ax.set_title("Evaporation (mm/day)")
    fig.colorbar(cs, ax=ax, label="mm/day")

    # 15) Condensation/Precip source P_cond (mm/day)
    ax = axes[4, 2]
    Pcond = getattr(gcm, "P_cond_flux_last", 0.0)
    if np.isscalar(Pcond):
        Pcond = np.full_like(gcm.T_s, float(Pcond))
    Pcond_mmday = np.nan_to_num(Pcond) * 86400.0
    cs = ax.contourf(grid.lon, grid.lat, Pcond_mmday, levels=20, cmap="BuPu")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)
    ax.set_title("Condensation P_cond (mm/day)")
    fig.colorbar(cs, ax=ax, label="mm/day")

    # Cosmetics
    for ax in axes.flatten():
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_xlim(0, 360)
        ax.set_ylim(-90, 90)

    # Overlay rivers and lakes (P014) on selected panels
    try:
        if routing is not None and int(os.getenv("QD_PLOT_RIVERS", "1")) == 1:
            import numpy as _np
            rd = routing.diagnostics()
            flow = _np.asarray(rd.get("flow_accum_kgps", _np.zeros_like(grid.lat_mesh)))
            river_min = float(os.getenv("QD_RIVER_MIN_KGPS", "1e6"))
            river_alpha = float(os.getenv("QD_RIVER_ALPHA", "0.35"))
            if _np.any(flow >= river_min):
                # Binary mask of major rivers (land-only), draw as contour lines
                river_mask = ((flow >= river_min) & (land_mask == 1)).astype(float)
                for _ax in (axes[0, 0], axes[2, 1]):  # Ts panel and Ocean panel
                    _ax.contour(grid.lon, grid.lat, river_mask, levels=[0.5],
                                colors="deepskyblue", linewidths=1.0, alpha=river_alpha)
            lake_mask = getattr(routing, "lake_mask", None)
            if lake_mask is not None and _np.any(lake_mask):
                lake_alpha = float(os.getenv("QD_LAKE_ALPHA", "0.40"))
                for _ax in (axes[0, 0], axes[2, 1]):
                    _ax.contour(grid.lon, grid.lat, lake_mask.astype(float), levels=[0.5],
                                colors="dodgerblue", linewidths=0.8, alpha=lake_alpha)
    except Exception:
        pass

    # Save
    filename = os.path.join(output_dir, f"state_day_{t_days:05.1f}.png")
    plt.savefig(filename, dpi=140)
    plt.close(fig)

def plot_true_color(grid, gcm, land_mask, t_days, output_dir, routing=None, eco=None, phyto=None):
    """
    Generates and saves a pseudo-true-color plot of the planet.

    Notes:
    - Sea-ice is now rendered from thickness (h_ice → optical ice_frac), NOT by T_s threshold,
      to be consistent with diagnostics与反照率。
    - Clouds are blended with configurable opacity，避免“整片纯白误判为冰”。
    """
    # Base colors
    ocean_color = np.array([0.10, 0.20, 0.50])
    land_color  = np.array([0.40, 0.30, 0.20])
    ice_color   = np.array([0.90, 0.90, 0.95])

    # Initialize RGB map
    rgb_map = np.zeros((grid.n_lat, grid.n_lon, 3), dtype=float)
    rgb_map[land_mask == 0] = ocean_color
    rgb_map[land_mask == 1] = land_color

    # Sea-ice from thickness (optical fraction)
    H_ice_ref = float(os.getenv("QD_HICE_REF", "0.5"))  # m
    ice_frac = 1.0 - np.exp(-np.maximum(gcm.h_ice, 0.0) / max(1e-6, H_ice_ref))
    # Render as "ice" only when optical coverage exceeds a small threshold
    ice_frac_thresh = float(os.getenv("QD_TRUECOLOR_ICE_FRAC", "0.15"))
    sea_ice_mask = (land_mask == 0) & (ice_frac >= ice_frac_thresh)
    rgb_map[sea_ice_mask] = ice_color

    # Land snow (from SWE/C_snow) overlay — render snow cover fraction on land
    try:
        if int(os.getenv("QD_TRUECOLOR_SNOW_BY_SWE", "1")) == 1 and hasattr(gcm, "C_snow_map_last"):
            C = np.nan_to_num(getattr(gcm, "C_snow_map_last"), nan=0.0)
            frac_thr = float(os.getenv("QD_SNOW_COVER_FRAC", "0.20"))  # coverage threshold
            vis_alpha = float(os.getenv("QD_SNOW_VIS_ALPHA", "0.60"))   # max blend strength
            land_snow_mask = (land_mask == 1) & (C >= frac_thr)
            # Scale alpha by coverage for smoother look
            alpha_map = vis_alpha * np.clip(C, 0.0, 1.0)
            alpha3 = alpha_map[..., None]
            lm3 = land_snow_mask[..., None]
            rgb_map = np.where(lm3, rgb_map * (1.0 - alpha3) + ice_color * alpha3, rgb_map)
    except Exception:
        pass

    # Vegetation coloring (optional TrueColor vegetation overlay)
    try:
        if int(os.getenv("QD_ECO_TRUECOLOR_VEG", "1")) == 1 and eco is not None:
            # Canopy factor f(LAI) ∈ [0,1] for land; fallback to 1 where pop is missing
            if getattr(eco, "pop", None) is not None:
                f_canopy = np.nan_to_num(eco.pop.canopy_reflectance_factor(), nan=0.0)
            else:
                f_canopy = np.zeros_like(gcm.T_s)
                f_canopy[land_mask == 1] = 1.0

            # Get banded surface albedo A_b^surface (NB×lat×lon); may be from daily cache
            Abands, _w_b = eco.get_surface_albedo_bands()
            if Abands is not None:
                NB = Abands.shape[0]
                lam = getattr(getattr(eco, "bands", None), "lambda_centers", None)
                if lam is None or len(lam) != NB:
                    lam = np.linspace(420.0, 680.0, NB)  # coarse fallback

                # Channel weights (simple Gaussians around canonical wavelengths)
                def _norm_gauss(x, mu, sigma):
                    w = np.exp(-((x - mu) ** 2) / (2.0 * sigma ** 2))
                    s = float(np.sum(w)) + 1e-12
                    return w / s

                wr = _norm_gauss(lam, 610.0, 50.0)
                wg = _norm_gauss(lam, 550.0, 40.0)
                wb = _norm_gauss(lam, 460.0, 40.0)

                # Dynamic per-band irradiance from dual stars → day/night & color modulation
                try:
                    from pygcm.ecology.spectral import dual_star_insolation_to_bands
                    I_b = dual_star_insolation_to_bands(getattr(gcm, "isr_A", np.zeros_like(gcm.T_s)),
                                                        getattr(gcm, "isr_B", np.zeros_like(gcm.T_s)),
                                                        eco.bands)  # [NB,lat,lon]
                    I_tot = np.maximum(getattr(gcm, "isr", np.zeros_like(gcm.T_s)), 0.0)
                    eps = 1e-12
                    w_rel = np.zeros_like(I_b)
                    mask = (I_tot > eps)
                    if np.any(mask):
                        w_rel[:, mask] = I_b[:, mask] / (I_tot[mask][None, ...] + eps)
                except Exception:
                    # Fallback: flat weights
                    w_rel = np.ones((NB, *gcm.T_s.shape), dtype=float) / float(NB)

                # Compute per-channel reflected intensity maps by band summation with I_b-relative weights
                Rr = np.nansum(Abands * (wr[:, None, None] * w_rel), axis=0)
                Rg = np.nansum(Abands * (wg[:, None, None] * w_rel), axis=0)
                Rb = np.nansum(Abands * (wb[:, None, None] * w_rel), axis=0)
                veg_rgb = np.stack([Rr, Rg, Rb], axis=-1)
                veg_rgb = np.clip(veg_rgb, 0.0, 1.0)

                # Optional saturation/gamma shaping for vegetation appearance
                try:
                    gamma = float(os.getenv("QD_ECO_TRUECOLOR_GAMMA", "1.8"))
                except Exception:
                    gamma = 1.8
                if gamma > 0:
                    veg_rgb = np.clip(veg_rgb, 0.0, 1.0) ** (1.0 / gamma)
                # Saturation boost: scale deviation from per-pixel mean
                try:
                    sat = float(os.getenv("QD_ECO_TRUECOLOR_SAT", "1.35"))
                except Exception:
                    sat = 1.35
                if sat != 1.0:
                    m = np.mean(veg_rgb, axis=-1, keepdims=True)
                    veg_rgb = np.clip(m + sat * (veg_rgb - m), 0.0, 1.0)

                # Mix vegetation with soil on land using canopy factor
                f = np.clip(f_canopy, 0.0, 1.0)[..., None]
                land3 = (land_mask == 1)[..., None]
                rgb_map = np.where(land3, rgb_map * (1.0 - f) + veg_rgb * f, rgb_map)
    except Exception:
        pass

    # Ocean color overlay from phytoplankton (P017)
    try:
        if int(os.getenv("QD_PLOT_OCEANCOLOR", "1")) == 1 and phyto is not None:
            alpha_bands, _alpha_scalar = phyto.get_alpha_maps()
            if alpha_bands is not None:
                NB = alpha_bands.shape[0]
                lam = getattr(getattr(phyto, "bands", None), "lambda_centers", None)
                if lam is None or len(lam) != NB:
                    lam = np.linspace(420.0, 680.0, NB)

                # Channel weights similar to vegetation block
                def _norm_gauss(x, mu, sigma):
                    w = np.exp(-((x - mu) ** 2) / (2.0 * sigma ** 2))
                    s = float(np.sum(w)) + 1e-12
                    return w / s

                wr = _norm_gauss(lam, 610.0, 50.0)
                wg = _norm_gauss(lam, 550.0, 40.0)
                wb = _norm_gauss(lam, 460.0, 40.0)

                # Dynamic per-band irradiance weights where possible
                try:
                    from pygcm.ecology.spectral import dual_star_insolation_to_bands
                    bands = getattr(phyto, "bands", None)
                    I_b = dual_star_insolation_to_bands(getattr(gcm, "isr_A", np.zeros_like(gcm.T_s)),
                                                        getattr(gcm, "isr_B", np.zeros_like(gcm.T_s)),
                                                        bands) if bands is not None else None
                    I_tot = np.maximum(getattr(gcm, "isr", np.zeros_like(gcm.T_s)), 0.0)
                    eps = 1e-12
                    if I_b is not None:
                        w_rel = np.zeros_like(I_b)
                        mask = (I_tot > eps)
                        if np.any(mask):
                            w_rel[:, mask] = I_b[:, mask] / (I_tot[mask][None, ...] + eps)
                    else:
                        w_rel = np.ones((NB, *gcm.T_s.shape), dtype=float) / float(NB)
                except Exception:
                    w_rel = np.ones((NB, *gcm.T_s.shape), dtype=float) / float(NB)

                # Per-channel water reflectance by bands
                Rr = np.nansum(alpha_bands * (wr[:, None, None] * w_rel), axis=0)
                Rg = np.nansum(alpha_bands * (wg[:, None, None] * w_rel), axis=0)
                Rb = np.nansum(alpha_bands * (wb[:, None, None] * w_rel), axis=0)
                water_rgb = np.stack([Rr, Rg, Rb], axis=-1)
                water_rgb = np.clip(water_rgb, 0.0, 1.0)

                # Optional shaping
                try:
                    gamma = float(os.getenv("QD_OC_GAMMA", os.getenv("QD_ECO_TRUECOLOR_GAMMA", "2.2")))
                except Exception:
                    gamma = 2.2
                if gamma > 0:
                    water_rgb = np.clip(water_rgb, 0.0, 1.0) ** (1.0 / gamma)

                # Blend factor for overlay
                try:
                    blend = float(os.getenv("QD_OC_BLEND", "0.85"))
                except Exception:
                    blend = 0.85

                # Apply to open ocean (avoid overwriting sea-ice tiles)
                ocean_open_mask = (land_mask == 0) & (~sea_ice_mask)
                rgb_map[ocean_open_mask] = (
                    rgb_map[ocean_open_mask] * (1.0 - blend) + water_rgb[ocean_open_mask] * blend
                )
    except Exception:
        pass

    if int(os.getenv("QD_TRUECOLOR_SNOW_BY_TS", "0")) == 1:
        snow_thresh = float(os.getenv("QD_SNOW_THRESH", "273.15"))
        land_snow_mask = (land_mask == 1) & (gcm.T_s <= snow_thresh)
        # 轻微偏白，避免与云混淆
        rgb_map[land_snow_mask] = 0.97 * ice_color

    # Cloud overlay (semi-transparent white)
    cloud_alpha = float(os.getenv("QD_TRUECOLOR_CLOUD_ALPHA", "0.60"))  # 0..1
    cloud_white = float(os.getenv("QD_TRUECOLOR_CLOUD_WHITE", "0.95"))  # 0..1
    cloud_layer = np.stack([gcm.cloud_cover, gcm.cloud_cover, gcm.cloud_cover], axis=-1)
    rgb_map = rgb_map * (1.0 - cloud_alpha * cloud_layer) + (cloud_alpha * cloud_layer) * cloud_white

    # Rivers/Lakes overlay (blend into RGB map)
    try:
        if routing is not None and int(os.getenv("QD_PLOT_RIVERS", "1")) == 1:
            import numpy as _np
            rd = routing.diagnostics()
            flow = _np.asarray(rd.get("flow_accum_kgps", _np.zeros_like(gcm.T_s)))
            river_min = float(os.getenv("QD_RIVER_MIN_KGPS", "1e6"))
            river_color = np.array([0.05, 0.35, 0.90])
            river_alpha = float(os.getenv("QD_RIVER_ALPHA", "0.45"))
            land3 = (land_mask == 1).astype(float)[..., None]
            river_mask = ((flow >= river_min).astype(float)[..., None]) * land3
            rgb_map = rgb_map * (1.0 - river_alpha * river_mask) + river_color * (river_alpha * river_mask)
        lake_mask = getattr(routing, "lake_mask", None)
        if lake_mask is not None and np.any(lake_mask):
            lake_color = np.array([0.15, 0.55, 0.95])
            lake_alpha = float(os.getenv("QD_LAKE_ALPHA", "0.40"))
            lake_mask3 = (lake_mask.astype(float) * (land_mask == 1).astype(float))[..., None]
            rgb_map = rgb_map * (1.0 - lake_alpha * lake_mask3) + lake_color * (lake_alpha * lake_mask3)
    except Exception:
        pass

    # Clamp
    rgb_map = np.clip(rgb_map, 0.0, 1.0)

    # Plotting
    fig, ax = plt.subplots(1, 1, figsize=(12, 6), constrained_layout=True)
    ax.imshow(rgb_map, extent=[0, 360, -90, 90], origin='lower')
    ax.set_title(f"Qingdai 'True Color' at Day {t_days:.2f}")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")

    # Save
    filename = os.path.join(output_dir, f"true_color_day_{t_days:05.1f}.png")
    plt.savefig(filename)
    plt.close(fig)

    # Console diagnostics for consistency with SeaIce logs
    try:
        w = np.maximum(np.cos(np.deg2rad(grid.lat_mesh)), 0.0)
        sea_ice_area = float((w * sea_ice_mask).sum() / (w.sum() + 1e-15))
        mean_h_ice = float(gcm.h_ice[sea_ice_mask].mean()) if np.any(sea_ice_mask) else 0.0
        print(f"[TrueColor] sea_ice_area≈{sea_ice_area:.3f}, mean_h_ice={mean_h_ice:.3f} m (thr={ice_frac_thresh}, alpha={cloud_alpha})")
    except Exception:
        pass

def plot_ocean(grid, ocean, land_mask, t_days, output_dir):
    """
    Plot ocean diagnostics:
    - SST (°C)
    - Ocean surface currents (quiver, sub-sampled)
    """
    import numpy as np
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6), constrained_layout=True)

    # 1) SST (°C)
    sst_c = np.nan_to_num(ocean.Ts - 273.15)
    sst_plot = ax1.contourf(grid.lon, grid.lat, sst_c, levels=20, cmap="coolwarm")
    ax1.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)
    ax1.set_title(f"SST (°C) at Day {t_days:.2f}")
    ax1.set_xlabel("Longitude")
    ax1.set_ylabel("Latitude")
    ax1.set_xlim(0, 360)
    ax1.set_ylim(-90, 90)
    fig.colorbar(sst_plot, ax=ax1, label="°C")

    # 2) Ocean currents (quiver, sub-sampled)
    # Sub-sample for readability
    step_lat = max(1, grid.n_lat // 30)
    step_lon = max(1, grid.n_lon // 30)
    lon_q = grid.lon_mesh[::step_lat, ::step_lon]
    lat_q = grid.lat_mesh[::step_lat, ::step_lon]
    uo_q = np.nan_to_num(ocean.uo[::step_lat, ::step_lon])
    vo_q = np.nan_to_num(ocean.vo[::step_lat, ::step_lon])

    speed = np.sqrt(ocean.uo**2 + ocean.vo**2)
    sp_plot = ax2.contourf(grid.lon, grid.lat, speed, levels=20, cmap="viridis")
    ax2.quiver(lon_q, lat_q, uo_q, vo_q, color="white", scale=400, width=0.002)
    ax2.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.7)
    ax2.set_title(f"Ocean Currents (m/s) at Day {t_days:.2f}")
    ax2.set_xlabel("Longitude")
    ax2.set_ylabel("Latitude")
    ax2.set_xlim(0, 360)
    ax2.set_ylim(-90, 90)
    fig.colorbar(sp_plot, ax=ax2, label="m/s")

    # Save figure
    fname = os.path.join(output_dir, f"ocean_day_{t_days:05.1f}.png")
    plt.savefig(fname)
    plt.close(fig)


def plot_plankton_species(grid, phyto, land_mask, t_days, output_dir):
    """
    Save plankton species 0/1 density (mg Chl m^-3) maps to output/plankton/.
    Note:
    - This plot uses raw C_phyto_s (mg Chl m^-3) without any day/night (irradiance) weighting
      or cloud/true-color overlays. Only land is masked to NaN for clarity.
    - Optional env QD_PHYTO_VMAX can fix the upper color limit; otherwise uses 99th percentile.
    """
    import os as _os
    import numpy as _np
    import matplotlib.pyplot as _plt
    try:
        C = getattr(phyto, "C_phyto_s", None)  # shape [S, n_lat, n_lon]
        if C is None or C.ndim != 3:
            return
        S = int(C.shape[0])
        if S <= 0:
            return
        _os.makedirs(_os.path.join(output_dir, "plankton"), exist_ok=True)

        land_bool = (land_mask == 1)

        def _plot_one(spec_idx: int):
            # Raw concentration; do NOT modulate by irradiance, clouds, or day/night.
            field = C[spec_idx, :, :].astype(float).copy()
            # Mask land to NaN so color scale reflects ocean values only.
            field[land_bool] = _np.nan

            # Determine color limits (vmin=0, vmax from env or 99th percentile of ocean values)
            try:
                vmax_env = _os.getenv("QD_PHYTO_VMAX")
                if vmax_env is not None and vmax_env.strip() != "":
                    vmax = float(vmax_env)
                else:
                    vmax = float(_np.nanpercentile(field, 99.0))
            except Exception:
                vmax = float(_np.nanmax(field))
            if not _np.isfinite(vmax) or vmax <= 0.0:
                vmax = 1.0e-3
            vmin = 0.0
            levels = _np.linspace(vmin, vmax, 21)

            fig, ax = _plt.subplots(1, 1, figsize=(12, 5), constrained_layout=True)
            cs = ax.contourf(grid.lon, grid.lat, field, levels=levels, cmap="viridis")
            # Coastline
            ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.6)
            ax.set_title(f"Plankton species {spec_idx} (mg Chl m$^{{-3}}$) Day {t_days:.2f}")
            ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
            ax.set_xlim(0, 360); ax.set_ylim(-90, 90)
            cbar = fig.colorbar(cs, ax=ax, label="mg Chl m$^{-3}$")

            fname = _os.path.join(output_dir, "plankton", f"species{spec_idx}_day_{t_days:05.1f}.png")
            _plt.savefig(fname, dpi=140)
            _plt.close(fig)

        # Always plot species 0；若存在则再画 species 1
        _plot_one(0)
        if S >= 2:
            _plot_one(1)
    except Exception:
        # Non-fatal: just skip visualization
        pass

def plot_isr_components(grid, gcm, t_days, output_dir):
    """
    Save a diagnostic figure showing per-star incoming shortwave (ISR) components
    to verify the expected double centers (subsolar points) from the two stars.
    """
    import numpy as np
    import matplotlib.pyplot as plt

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)

    # Choose common levels for easier visual comparison
    vmin = 0.0
    vmax = max(np.max(gcm.isr_A), np.max(gcm.isr_B))
    levels = np.linspace(vmin, vmax, 21)

    csA = axA.contourf(grid.lon, grid.lat, gcm.isr_A, levels=levels, cmap='magma')
    axA.set_title(f"ISR - Star A (Day {t_days:.2f})")
    axA.set_xlabel("Longitude")
    axA.set_ylabel("Latitude")
    axA.set_xlim(0, 360)
    axA.set_ylim(-90, 90)
    fig.colorbar(csA, ax=axA, label="W/m^2")

    csB = axB.contourf(grid.lon, grid.lat, gcm.isr_B, levels=levels, cmap='magma')
    axB.set_title(f"ISR - Star B (Day {t_days:.2f})")
    axB.set_xlabel("Longitude")
    axB.set_ylabel("Latitude")
    axB.set_xlim(0, 360)
    axB.set_ylim(-90, 90)
    fig.colorbar(csB, ax=axB, label="W/m^2")

    # Mark subsolar points (maxima) for each component and report their great-circle separation
    try:
        idxA = np.unravel_index(np.argmax(gcm.isr_A), gcm.isr_A.shape)
        idxB = np.unravel_index(np.argmax(gcm.isr_B), gcm.isr_B.shape)
        lonA, latA = grid.lon[idxA[1]], grid.lat[idxA[0]]
        lonB, latB = grid.lon[idxB[1]], grid.lat[idxB[0]]
        axA.scatter([lonA], [latA], c='cyan', s=40, marker='x', label='A center')
        axB.scatter([lonB], [latB], c='yellow', s=40, marker='+', label='B center')
        axA.legend(loc='upper right', fontsize=8)
        axB.legend(loc='upper right', fontsize=8)

        # Great-circle separation
        import math
        phi1 = math.radians(latA); lam1 = math.radians(lonA)
        phi2 = math.radians(latB); lam2 = math.radians(lonB)
        dlam = lam2 - lam1
        # Haversine
        d_sigma = 2 * math.asin(math.sqrt(
            math.sin((phi2 - phi1)/2)**2 +
            math.cos(phi1)*math.cos(phi2)*math.sin(dlam/2)**2
        ))
        separation_deg = math.degrees(d_sigma)
        print(f"[Diagnostics] Day {t_days:.2f}: Subsolar separation ≈ {separation_deg:.2f}° "
              f"(A: lon={lonA:.1f}°, lat={latA:.1f}°; B: lon={lonB:.1f}°, lat={latB:.1f}°)")
    except Exception as e:
        print(f"[Diagnostics] Could not compute subsolar separation: {e}")

    fname = os.path.join(output_dir, f"isr_components_day_{t_days:05.1f}.png")
    plt.savefig(fname)
    plt.close(fig)


def plot_ecology(grid, land_mask, t_days, output_dir, *, lai=None, alpha_ecology=None, alpha_banded=None, canopy_height=None, species_density=None):
    """
    Save ecology diagnostics:
    - LAI map
    - Ecology alpha (scalar) map
    - Optional banded alpha map
    - Optional canopy height proxy (m)
    - Optional species density summary (per species)
    """
    import numpy as _np
    import matplotlib.pyplot as _plt

    have_bands = alpha_banded is not None
    have_h = canopy_height is not None
    nspecies = 0 if (species_density is None) else len([m for m in species_density if m is not None])
    ncols = 3 + (1 if have_h else 0) + (1 if nspecies > 0 else 0)
    fig, axes = _plt.subplots(1, ncols, figsize=(6.5 * ncols, 5), constrained_layout=True)
    col = 0

    # 1) LAI
    ax = axes[col]; col += 1
    if lai is None:
        ax.text(0.5, 0.5, "LAI not available", ha="center", va="center")
    else:
        cs = ax.contourf(grid.lon, grid.lat, _np.nan_to_num(lai), levels=20, cmap="YlGn")
        fig.colorbar(cs, ax=ax, label="LAI")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.6)
    ax.set_title(f"LAI (Day {t_days:.2f})")
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude"); ax.set_xlim(0, 360); ax.set_ylim(-90, 90)

    # 2) Ecology alpha (scalar)
    ax = axes[col]; col += 1
    if alpha_ecology is None:
        ax.text(0.5, 0.5, "alpha_ecology not available", ha="center", va="center")
    else:
        cs = ax.contourf(grid.lon, grid.lat, _np.nan_to_num(alpha_ecology), levels=_np.linspace(0, 0.8, 17), cmap="cividis")
        fig.colorbar(cs, ax=ax, label="alpha")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.6)
    ax.set_title("Ecology alpha (scalar)")
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude"); ax.set_xlim(0, 360); ax.set_ylim(-90, 90)

    # 3) Banded alpha
    ax = axes[col]; col += 1
    if not have_bands:
        ax.text(0.5, 0.5, "alpha_banded not available", ha="center", va="center")
    else:
        cs = ax.contourf(grid.lon, grid.lat, _np.nan_to_num(alpha_banded), levels=_np.linspace(0, 0.8, 17), cmap="cividis")
        fig.colorbar(cs, ax=ax, label="alpha")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.6)
    ax.set_title("Ecology alpha (banded)")
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude"); ax.set_xlim(0, 360); ax.set_ylim(-90, 90)

    # 4) Canopy height
    if have_h:
        ax = axes[col]; col += 1
        cs = ax.contourf(grid.lon, grid.lat, _np.nan_to_num(canopy_height), levels=20, cmap="Greens")
        fig.colorbar(cs, ax=ax, label="m")
        ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.6)
        ax.set_title("Canopy height (m)")
        ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude"); ax.set_xlim(0, 360); ax.set_ylim(-90, 90)

    # 5) Species density (first species)
    if nspecies > 0:
        ax = axes[col]; col += 1
        m = species_density[0]
        cs = ax.contourf(grid.lon, grid.lat, _np.nan_to_num(m), levels=20, cmap="viridis")
        fig.colorbar(cs, ax=ax, label="arb. units")
        ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.6)
        ax.set_title("Species 0 density (proxy)")
        ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude"); ax.set_xlim(0, 360); ax.set_ylim(-90, 90)

    # --- Beta diversity (Whittaker) annotation on ecology panel ---
    try:
        if eco is not None and getattr(eco, "pop", None) is not None and getattr(eco.pop, "LAI_layers_SK", None) is not None:
            # L_s: per-species canopy mass proxy (Σ_k LAI_s,k) over land
            L_s = _np.sum(_np.maximum(eco.pop.LAI_layers_SK, 0.0), axis=1)  # [S,lat,lon]
            land = (land_mask == 1)
            # Area weights
            w = _np.maximum(_np.cos(_np.deg2rad(grid.lat_mesh)), 0.0)
            w_sum_land = float(_np.sum(w[land])) + 1e-15
            w_norm = w / w_sum_land
            # Per-pixel alpha diversity (effective species number) α_eff = exp(H), H=−Σ p_i ln p_i
            tot = _np.sum(L_s, axis=0)  # [lat,lon]
            alpha_eff = _np.full_like(tot, _np.nan, dtype=float)
            mask = land & (tot > 0)
            if _np.any(mask):
                P = L_s[:, mask] / (tot[mask][None, ...] + 1e-15)  # [S, nmask]
                H = -_np.sum(P * _np.log(P + 1e-15), axis=0)
                alpha_eff[mask] = _np.exp(H)
            alpha_mean = float(_np.nansum(alpha_eff[land] * w_norm[land]))
            # Gamma diversity: from global composition (area-weighted across land)
            T_s = _np.array([float(_np.nansum(L_s[s, :, :][land] * w_norm[land])) for s in range(L_s.shape[0])], dtype=float)
            T_sum = float(_np.sum(T_s)) + 1e-15
            p_gamma = T_s / T_sum
            H_gamma = float(-_np.sum(p_gamma * _np.log(p_gamma + 1e-15)))
            gamma_eff = float(_np.exp(H_gamma))
            beta_whittaker = float(gamma_eff / max(alpha_mean, 1e-12))
            # Annotate figure
            txt = f"Beta diversity (Whittaker): β≈{beta_whittaker:.2f} (ᾱ≈{alpha_mean:.2f}, γ≈{gamma_eff:.2f})"
            fig.suptitle(txt, fontsize=12)
    except Exception:
        # Non-fatal: simply skip annotation
        pass

    fname = os.path.join(output_dir, f"ecology_day_{t_days:05.1f}.png")
    _plt.savefig(fname, dpi=140)
    _plt.close(fig)


def _try_autogen_hydro_network(grid, land_mask, elevation, topo_nc, out_path) -> bool:
    """
    Attempt to auto-generate a hydrology routing network NetCDF when missing.
    Uses the same core routines as scripts/generate_hydrology_maps.py.
    Returns True on success, False otherwise.
    """
    try:
        from netCDF4 import Dataset
        import numpy as _np
        # Reuse helpers from the generator script
        from scripts.generate_hydrology_maps import (
            pit_fill,
            compute_flow_to_index,
            identify_lakes,
            compute_lake_outlets,
            topo_sort_flow_order,
        )
        elev0 = elevation if elevation is not None else _np.zeros_like(grid.lat_mesh, dtype=float)
        land = (land_mask.astype(_np.uint8))
        src = "procedural" if (not topo_nc or not os.path.exists(str(topo_nc))) else os.path.basename(str(topo_nc))
        print(f"[HydroRouting] Auto-generating network to '{out_path}' (source={src})...")
        elev_filled = pit_fill(elev0.copy(), land, max_iters=200, eps=1e-3)
        flow_to = compute_flow_to_index(grid, elev_filled, land)
        lake_mask, lake_id, n_lakes = identify_lakes(flow_to, land)
        lake_outlet_index = compute_lake_outlets(grid, elev_filled, lake_mask, lake_id, land) if int(_np.max(lake_id)) > 0 else None
        flow_order = topo_sort_flow_order(flow_to, land)

        os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
        with Dataset(out_path, "w") as ds:
            nlat, nlon = grid.n_lat, grid.n_lon
            ds.createDimension("lat", nlat)
            ds.createDimension("lon", nlon)
            n_land = int((land == 1).sum())
            ds.createDimension("n_land", n_land)
            n_lakes = int(_np.max(lake_id)) if lake_id is not None else 0
            if n_lakes > 0:
                ds.createDimension("n_lakes", n_lakes)

            vlat = ds.createVariable("lat", "f4", ("lat",))
            vlon = ds.createVariable("lon", "f4", ("lon",))
            vlat[:] = grid.lat.astype(_np.float32)
            vlon[:] = grid.lon.astype(_np.float32)

            def wvar(name, dtype, dims, data):
                var = ds.createVariable(name, dtype, dims)
                var[:] = data

            wvar("land_mask", "u1", ("lat", "lon"), land.astype(_np.uint8))
            wvar("elevation_filled", "f4", ("lat", "lon"), elev_filled.astype(_np.float32))
            wvar("flow_to_index", "i4", ("lat", "lon"), flow_to.astype(_np.int32))
            wvar("flow_order", "i4", ("n_land",), flow_order.astype(_np.int32))
            wvar("lake_mask", "u1", ("lat", "lon"), lake_mask.astype(_np.uint8))
            wvar("lake_id", "i4", ("lat", "lon"), lake_id.astype(_np.int32))
            if n_lakes > 0 and lake_outlet_index is not None:
                wvar("lake_outlet_index", "i4", ("n_lakes",), lake_outlet_index.astype(_np.int32))

            ds.setncattr("title", "Qingdai Hydrology Network (auto-generated)")
            ds.setncattr("indexing", "row-major (i=lon index, j=lat index), idx=j*n_lon+i")
            ds.setncattr("projection", "latlon")
            ds.setncattr("created_by", "scripts/run_simulation.py (auto)")
        print("[HydroRouting] Network auto-generation complete.")
        return True
    except Exception as e:
        print(f"[HydroRouting] Auto-generation failed: {e}")
        return False


# --- Optional: scalar advection helper for periodic lat-lon (cloud/diagnostics) ---
def _advect_scalar_periodic(field, u, v, dt, grid):
    """
    Semi-Lagrangian advection on regular lat-lon with longitudinal periodicity.
    Uses scipy.ndimage.map_coordinates; wraps in longitude (periodic). Latitude is approximated.
    """
    import numpy as _np
    try:
        from scipy.ndimage import map_coordinates, gaussian_filter  # noqa: F401 (gaussian optional)
    except Exception:
        from scipy.ndimage import map_coordinates  # type: ignore

    a = getattr(grid, "a", 6.371e6)
    dlat = grid.dlat_rad
    dlon = grid.dlon_rad
    coslat = _np.maximum(_np.cos(_np.deg2rad(grid.lat_mesh)), 0.5)

    # Convert velocities to grid displacement in index space
    dlam = u * dt / (a * coslat)
    dphi = v * dt / a
    dx = dlam / dlon
    dy = dphi / dlat

    JJ, II = _np.meshgrid(_np.arange(grid.n_lat), _np.arange(grid.n_lon), indexing="ij")
    dep_J = JJ - dy
    dep_I = II - dx

    adv = map_coordinates(field, [dep_J, dep_I], order=1, mode="wrap", prefilter=False)
    return adv


def main():
    """
    Main function to run the simulation.
    """
    print("--- Initializing Qingdai GCM ---")
    try:
        print(f"[JAX] Acceleration enabled: {JAX_IS_ENABLED()} (toggle via QD_USE_JAX=1; platform via QD_JAX_PLATFORM=cpu|gpu|tpu)")
    except Exception:
        pass

    # P020 Phase 0 switch (non-intrusive): when QD_USE_OO=1, instantiate façade and run stub.
    # Legacy path continues unless QD_USE_OO_STRICT=1.
    try:
        USE_OO = int(os.getenv("QD_USE_OO", "0")) == 1
        OO_STRICT = int(os.getenv("QD_USE_OO_STRICT", "0")) == 1
    except Exception:
        USE_OO, OO_STRICT = False, False
    if USE_OO:
        if QingdaiWorld is not None:
            try:
                world = QingdaiWorld.create_default()
                print("[P020] QD_USE_OO=1 → QingdaiWorld façade active (Phase 0).")
                # Phase 0 stub (does not alter legacy behavior)
                world.run()
            except Exception as _wo:
                print(f"[P020] world façade run stub raised: {_wo}")
            if OO_STRICT:
                print("[P020] QD_USE_OO_STRICT=1 → exiting legacy engine after façade run.")
                return
        else:
            print("[P020] QD_USE_OO=1 but pygcm.world is unavailable; continuing with legacy engine.")

    # 1. Initialization
    print("Creating grid...")
    grid = SphericalGrid(n_lat=121, n_lon=240)

    print("Creating topography...")
    topo_nc = os.getenv("QD_TOPO_NC")
    elevation = None
    if topo_nc and os.path.exists(topo_nc):
        try:
            elevation, land_mask, base_albedo_map, friction_map = load_topography_from_netcdf(topo_nc, grid)
        except Exception as e:
            print(f"[Topo] Failed to load '{topo_nc}': {e}\nFalling back to procedural generation.")
            land_mask = create_land_sea_mask(grid)
            base_albedo_map, friction_map = generate_base_properties(land_mask)
            elevation = None
        else:
            # Loader already prints stats
            pass
    else:
        land_mask = create_land_sea_mask(grid)
        base_albedo_map, friction_map = generate_base_properties(land_mask)
        # Log fallback stats
        LAT = grid.lat_mesh
        area_w = np.cos(np.deg2rad(LAT))
        achieved = float((area_w * (land_mask == 1)).sum() / (area_w.sum() + 1e-15))
        print(f"[Topo] Procedural topography (no external NetCDF). Land fraction: {achieved:.3f}")
        print(f"[Topo] Albedo stats (min/mean/max): {np.min(base_albedo_map):.3f}/{np.mean(base_albedo_map):.3f}/{np.max(base_albedo_map):.3f}")
        print(f"[Topo] Friction stats (min/mean/max): {np.min(friction_map):.2e}/{np.mean(friction_map):.2e}/{np.max(friction_map):.2e}")
    # Write standardized topography.nc (always)
    try:
        save_topography(os.path.join("data", "topography.nc"), grid, land_mask, base_albedo_map, friction_map, elevation=elevation)
        print("[Topo] Wrote standardized topography.nc")
    except Exception as _tw:
        print(f"[Topo] topography.nc write skipped: {_tw}")

    # --- Slab Ocean (P007 M1): construct per-grid surface heat capacity map ---
    # C_s_ocean = rho_w * c_p_w * H_mld; land uses smaller constant C_s_land
    rho_w = float(os.getenv("QD_RHO_W", "1000"))      # kg/m^3
    cp_w = float(os.getenv("QD_CP_W", "4200"))        # J/(kg K)
    H_mld = float(os.getenv("QD_MLD_M", "50"))        # m
    Cs_ocean = rho_w * cp_w * H_mld                   # J/m^2/K
    Cs_land = float(os.getenv("QD_CS_LAND", "3e6"))   # J/m^2/K
    Cs_ice = float(os.getenv("QD_CS_ICE", "5e6"))     # J/m^2/K (thin ice/snow effective capacity)
    C_s_map = np.where(land_mask == 1, Cs_land, Cs_ocean).astype(float)

    # Diagnostics
    LAT = grid.lat_mesh
    area_w = np.cos(np.deg2rad(LAT))
    area_w = np.maximum(area_w, 0.0)
    land_area = float((area_w * (land_mask == 1)).sum() / (area_w.sum() + 1e-15))
    print(f"[SlabOcean] H_mld={H_mld:.1f} m -> C_s_ocean={Cs_ocean:.2e} J/m^2/K, C_s_land={Cs_land:.2e}")
    print(f"[SlabOcean] Land fraction={land_area:.3f}; C_s stats (min/mean/max): {np.min(C_s_map):.2e}/{np.mean(C_s_map):.2e}/{np.max(C_s_map):.2e}")

    print("Initializing orbital mechanics...")
    orbital_sys = OrbitalSystem()

    print("Initializing thermal forcing...")
    forcing = ThermalForcing(grid, orbital_sys)

    # Initialize energy parameters and optional greenhouse autotuning
    print("Initializing energy parameters...")
    eparams = energy.get_energy_params_from_env()
    GH_LOCK = int(os.getenv("QD_GH_LOCK", "1")) == 1
    AUTOTUNE = (not GH_LOCK) and (int(os.getenv("QD_ENERGY_AUTOTUNE", "0")) == 1)
    TUNE_EVERY = int(os.getenv("QD_ENERGY_TUNE_EVERY", "50"))
    if GH_LOCK:
        try:
            g_fixed = float(os.getenv("QD_GH_FACTOR", "0.40"))
        except Exception:
            g_fixed = 0.40
        print(f"[Greenhouse] Lock enabled: fixed g={g_fixed:.2f}; autotune disabled.")

    print("Initializing dynamics core with surface friction and greenhouse effect...")
    gcm = SpectralModel(
        grid, friction_map, H=8000, tau_rad=10 * 24 * 3600, greenhouse_factor=float(os.getenv("QD_GH_FACTOR", "0.40")),
        C_s_map=C_s_map, land_mask=land_mask, Cs_ocean=Cs_ocean, Cs_land=Cs_land, Cs_ice=Cs_ice
    )

    # --- Ocean model (P011): optional dynamic slab ocean (M1+M2+M3) ---
    # Default: enabled (no configuration needed). Set QD_USE_OCEAN=0 to disable explicitly.
    USE_OCEAN = int(os.getenv("QD_USE_OCEAN", "1")) == 1
    ocean = None
    if USE_OCEAN:
        try:
            H_ocean = float(os.getenv("QD_OCEAN_H_M", str(H_mld)))
        except Exception:
            H_ocean = H_mld
        # Initialize ocean SST from current surface temperature over ocean; fill land with 288 K placeholder
        init_Ts = np.where(land_mask == 0, gcm.T_s, 288.0)
        ocean = WindDrivenSlabOcean(grid, land_mask, H_ocean, init_Ts=init_Ts)
        print(f"[Ocean] Dynamic slab ocean enabled: H={H_ocean:.1f} m, CD={float(os.getenv('QD_CD', '1.5e-3'))}, R_bot={float(os.getenv('QD_R_BOT', '1.0e-6'))}")
    else:
        print("[Ocean] Dynamic slab ocean disabled (QD_USE_OCEAN=0).")

    # --- Hydrology (P009): reservoirs and parameters ---
    hydro_params = get_hydrology_params_from_env()
    W_land = np.zeros_like(grid.lat_mesh, dtype=float)   # land bucket water (kg m^-2)
    S_snow = np.zeros_like(grid.lat_mesh, dtype=float)   # land snow water-equivalent (kg m^-2)
    _hydro_prev_total = None
    _hydro_prev_time = None

    # Routing (P014): optional river routing and lakes via offline network
    routing = None
    try:
        HYDRO_ENABLED = (int(os.getenv("QD_HYDRO_ENABLE", "1")) == 1)
        hydro_net = os.getenv("QD_HYDRO_NETCDF", "data/hydrology.nc")
        if HYDRO_ENABLED:
            # If network file missing, attempt auto-generation once
            if not (hydro_net and os.path.exists(hydro_net)):
                ok_gen = _try_autogen_hydro_network(grid, land_mask, elevation, topo_nc, hydro_net)
                if not ok_gen or not os.path.exists(hydro_net):
                    print(f"[HydroRouting] Enabled but network not available; running WITHOUT routing "
                          f"(QD_HYDRO_NETCDF='{hydro_net}').")
            # Initialize routing if file now exists
            if hydro_net and os.path.exists(hydro_net):
                routing = RiverRouting(
                    grid,
                    hydro_net,
                    dt_hydro_hours=float(os.getenv("QD_HYDRO_DT_HOURS", "6")),
                    treat_lake_as_water=(int(os.getenv("QD_TREAT_LAKE_AS_WATER", "1")) == 1),
                    alpha_lake=(float(os.getenv("QD_ALPHA_LAKE")) if os.getenv("QD_ALPHA_LAKE") else None),
                    diag=(int(os.getenv("QD_HYDRO_DIAG", "1")) == 1),
                )
                print(f"[HydroRouting] Enabled with network '{hydro_net}'.")
        else:
            print("[HydroRouting] Disabled by QD_HYDRO_ENABLE=0.")
    except Exception as e:
        print(f"[HydroRouting] Initialization skipped due to error: {e}")
        routing = None

    # --- Ecology (P015 M1): hourly substep adapter (land-only scalar alpha) ---
    ECO_ENABLED = int(os.getenv("QD_ECO_ENABLE", "1")) == 1
    SUBDAILY_ENABLED = int(os.getenv("QD_ECO_SUBDAILY_ENABLE", "1")) == 1
    eco = None
    # Basic env echo to help users verify switches
    try:
        _nb = os.getenv("QD_ECO_SPECTRAL_BANDS", "default")
        _mode = os.getenv("QD_ECO_TOA_TO_SURF_MODE", "simple")
        _use_lai = os.getenv("QD_ECO_USE_LAI", "1")
        print(f"[Ecology] env: ENABLE={ECO_ENABLED} SUBDAILY={SUBDAILY_ENABLED} NB={_nb} MODE={_mode} USE_LAI={_use_lai}")
    except Exception:
        pass
    if ECO_ENABLED and EcologyAdapter is not None:
        try:
            eco = EcologyAdapter(grid, land_mask)
            if eco is not None and int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                print("[Ecology] Adapter initialized successfully.")
        except Exception as e:
            print(f"[Ecology] Adapter init failed: {e}")
            eco = None
    elif ECO_ENABLED:
        print("[Ecology] Adapter not available (import failed).")

    # --- Phytoplankton (P017): daily mixed-layer ocean color coupling ---
    PHYTO_ENABLED = int(os.getenv("QD_PHYTO_ENABLE", "1")) == 1
    PHYTO_COUPLE = int(os.getenv("QD_PHYTO_ALBEDO_COUPLE", "1")) == 1
    PHYTO_FEEDBACK_MODE = os.getenv("QD_PHYTO_FEEDBACK_MODE", "daily").strip().lower()
    # Enable ocean-current advection for phytoplankton (per-physics-step), default on
    PHYTO_ADVECTION = int(os.getenv("QD_PHYTO_ADVECTION", "1")) == 1
    phyto = None
    if PHYTO_ENABLED and 'PhytoManager' in globals() and PhytoManager is not None:
        try:
            try:
                H_ocean_for_phyto = float(os.getenv("QD_OCEAN_H_M", os.getenv("QD_MLD_M", "50")))
            except Exception:
                H_ocean_for_phyto = 50.0
            phyto = PhytoManager(grid, land_mask, H_mld_m=H_ocean_for_phyto, diag=(int(os.getenv("QD_PHYTO_DIAG", "1")) == 1))
            if int(os.getenv("QD_PHYTO_DIAG", "1")) == 1:
                print("[Phyto] Manager initialized.")
        except Exception as e:
            print(f"[Phyto] Init failed: {e}")
            phyto = None
    elif PHYTO_ENABLED:
        print("[Phyto] Manager not available (import failed).")

    # Phyto NPZ autosave deprecated:
    # Initialization now relies on standardized files in data/:
    #   - data/plankton.nc     (gridded distributions: C_phyto_s/alpha/Kd490/N)
    #   - data/plankton.json   (bio/optics & bands parameters)
    # See the section "Optional: load plankton bio/distribution on startup" below.
    # To force randomized/default init, use:
    #   QD_PHYTO_INIT_RANDOM=1  → randomize_state()
    # Otherwise falls back to reset_default_state() via that block when no files are present.

    # Optional: load plankton bio/distribution on startup (controlled by QD_LOAD_PLANKTON=1)
    try:
        if phyto is not None and int(os.getenv("QD_LOAD_PLANKTON", "1")) == 1:
            # Bio/optics (JSON)
            try:
                pj = os.path.join("data", "plankton.json")
                if os.path.exists(pj):
                    ok_bio = phyto.load_bio_json(pj, on_mismatch=os.getenv("QD_PLANKTON_BIO_ON_MISMATCH", "keep"))
                    if int(os.getenv("QD_PHYTO_DIAG", "1")) == 1:
                        print(f"[Phyto] plankton.json load {'OK' if ok_bio else 'skipped/failed'}.")
            except Exception as _plj:
                if int(os.getenv("QD_PHYTO_DIAG", "1")) == 1:
                    print(f"[Phyto] plankton.json load skipped: {_plj}")
            # Distributions (NetCDF)
            try:
                pnc = os.path.join("data", "plankton.nc")
                if os.path.exists(pnc):
                    ok_nc = phyto.load_distribution_nc(pnc, on_mismatch=os.getenv("QD_PLANKTON_DIST_ON_MISMATCH", "keep"))
                    if int(os.getenv("QD_PHYTO_DIAG", "1")) == 1:
                        print(f"[Phyto] plankton.nc load {'OK' if ok_nc else 'skipped/failed'}.")
            except Exception as _pln:
                if int(os.getenv("QD_PHYTO_DIAG", "1")) == 1:
                    print(f"[Phyto] plankton.nc load skipped: {_pln}")
    except Exception:
        pass

    # --- Individual pool settings (vectorized sampled individuals) ---
    INDIV_ENABLED = int(os.getenv("QD_ECO_INDIV_ENABLE", "1")) == 1
    indiv = None
    if INDIV_ENABLED and (eco is not None) and (IndividualPool is not None):
        try:
            indiv = IndividualPool(
                grid,
                land_mask,
                eco,
                # defaults are read inside pool from env as well; passing for clarity
                sample_frac=float(os.getenv("QD_ECO_INDIV_SAMPLE_FRAC", "0.02")),
                per_cell=int(os.getenv("QD_ECO_INDIV_PER_CELL", "150")),
                substeps_per_day=int(os.getenv("QD_ECO_INDIV_SUBSTEPS_PER_DAY", "10")),
                diag=(int(os.getenv("QD_ECO_DIAG", "1")) == 1),
            )
        except Exception as _ie:
            print(f"[EcoIndiv] init failed: {_ie}")
            indiv = None
    elif INDIV_ENABLED and IndividualPool is None:
        print("[EcoIndiv] Module not available (import failed).")

    # --- Diversity diagnostics settings ---
    DIVERSITY_ENABLED = int(os.getenv("QD_ECO_DIVERSITY_ENABLE", "0")) == 1
    try:
        DIVERSITY_EVERY_DAYS = float(os.getenv("QD_ECO_DIVERSITY_EVERY_DAYS", "10"))
    except Exception:
        DIVERSITY_EVERY_DAYS = 10.0
    if DIVERSITY_ENABLED and eco_diversity is None:
        print("[Diversity] Module not available (import failed).")

    # --- Restart load or banded initialization ---
    t0_seconds = 0.0
    restart_in = os.getenv("QD_RESTART_IN")
    if restart_in and os.path.exists(restart_in):
        try:
            rst = load_restart(restart_in)
            # Basic shape checks (optional): assume matching grid
            # Atmospheric / surface
            if rst.get("u") is not None: gcm.u = rst["u"]
            if rst.get("v") is not None: gcm.v = rst["v"]
            if rst.get("h") is not None: gcm.h = rst["h"]
            if rst.get("T_s") is not None: gcm.T_s = rst["T_s"]
            if rst.get("cloud_cover") is not None: gcm.cloud_cover = np.clip(rst["cloud_cover"], 0.0, 1.0)
            if rst.get("q") is not None and hasattr(gcm, "q"): gcm.q = rst["q"]
            if rst.get("h_ice") is not None and hasattr(gcm, "h_ice"): gcm.h_ice = np.maximum(rst["h_ice"], 0.0)
            # Ocean
            if ocean is not None:
                if rst.get("uo") is not None: ocean.uo = rst["uo"]
                if rst.get("vo") is not None: ocean.vo = rst["vo"]
                if rst.get("eta") is not None: ocean.eta = rst["eta"]
                if rst.get("Ts") is not None: ocean.Ts = rst["Ts"]
            # Hydrology / Cryo
            if rst.get("W_land") is not None: W_land = rst["W_land"]
            if rst.get("S_snow") is not None: S_snow = rst["S_snow"]
            # Persisted snow optical coverage for immediate visualization continuity
            try:
                C_snow_rst = rst.get("C_snow", None)
                if C_snow_rst is not None:
                    gcm.C_snow_map_last = C_snow_rst
            except Exception:
                pass
            # Load genes autosave JSON (if present) before LAI/species weights
            if int(os.getenv("QD_AUTOSAVE_LOAD", "1")) == 1 and eco is not None:
                try:
                    genes_path = os.getenv("QD_ECO_GENES_JSON_PATH") or os.path.join("data", "genes.json")
                    if os.path.exists(genes_path) and hasattr(eco, "load_genes_json"):
                        eco.load_genes_json(genes_path)
                except Exception as _egl:
                    if int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                        print(f"[Ecology] genes autosave load skipped: {_egl}")
            # Try to load ecology autosave (optional, persists LAI/species weights & species bands)
            if int(os.getenv("QD_AUTOSAVE_LOAD", "1")) == 1 and eco is not None:
                try:
                    eco_path = os.getenv("QD_ECO_AUTOSAVE_PATH") or os.path.join("data", "ecology.nc")
                    if os.path.exists(eco_path):
                        ok_eco = False
                        if hasattr(eco, "load_autosave"):
                            ok_eco = bool(eco.load_autosave(eco_path, on_mismatch=os.getenv("QD_ECO_ON_MISMATCH", "fallback")))
                        else:
                            if int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                                print("[Ecology] adapter has no load_autosave; cannot load ecology.nc")
                        if int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                            print(f"[Ecology] autosave load {'OK' if ok_eco else 'failed'} from '{eco_path}'")
                except Exception as _ecl:
                    if int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                        print(f"[Ecology] autosave load skipped: {_ecl}")
            # Load astronomical epoch if present
            try:
                t_loaded = rst.get("t_seconds", None)
                if t_loaded is not None:
                    t0_seconds = float(t_loaded)
            except Exception:
                pass
            print(f"[Restart] Loaded state from '{restart_in}'.")
            # Optionally load standardized ocean.nc to override ocean fields
            try:
                if ocean is not None and int(os.getenv("QD_LOAD_OCEAN", "1")) == 1:
                    ocean_nc = os.path.join("data", "ocean.nc")
                    if os.path.exists(ocean_nc):
                        rst_o = load_ocean(ocean_nc)
                        if rst_o.get("uo") is not None: ocean.uo = rst_o["uo"]
                        if rst_o.get("vo") is not None: ocean.vo = rst_o["vo"]
                        if rst_o.get("eta") is not None: ocean.eta = rst_o["eta"]
                        if rst_o.get("Ts") is not None: ocean.Ts = rst_o["Ts"]
                        print("[Restart] Ocean state overridden from 'data/ocean.nc'.")
            except Exception as _loe:
                print(f"[Restart] ocean.nc load skipped: {_loe}")
        except Exception as e:
            print(f"[Restart] Failed to load '{restart_in}': {e}\nContinuing with fresh init.")
            apply_banded_initial_ts(grid, gcm, ocean, land_mask)
    else:
        # Fallback: try autosave checkpoint if enabled
        autosave_nc = os.path.join("data", "atmosphere.nc")
        if int(os.getenv("QD_AUTOSAVE_LOAD", "1")) == 1 and os.path.exists(autosave_nc):
            try:
                rst = load_restart(autosave_nc)
                if rst.get("u") is not None: gcm.u = rst["u"]
                if rst.get("v") is not None: gcm.v = rst["v"]
                if rst.get("h") is not None: gcm.h = rst["h"]
                if rst.get("T_s") is not None: gcm.T_s = rst["T_s"]
                if rst.get("cloud_cover") is not None: gcm.cloud_cover = np.clip(rst["cloud_cover"], 0.0, 1.0)
                if rst.get("q") is not None and hasattr(gcm, "q"): gcm.q = rst["q"]
                if rst.get("h_ice") is not None and hasattr(gcm, "h_ice"): gcm.h_ice = np.maximum(rst["h_ice"], 0.0)
                if ocean is not None:
                    if rst.get("uo") is not None: ocean.uo = rst["uo"]
                    if rst.get("vo") is not None: ocean.vo = rst["vo"]
                    if rst.get("eta") is not None: ocean.eta = rst["eta"]
                    if rst.get("Ts") is not None: ocean.Ts = rst["Ts"]
                if rst.get("W_land") is not None: W_land = rst["W_land"]
                if rst.get("S_snow") is not None: S_snow = rst["S_snow"]
                # Persisted snow optical coverage for immediate visualization continuity
                try:
                    C_snow_rst = rst.get("C_snow", None)
                    if C_snow_rst is not None:
                        gcm.C_snow_map_last = C_snow_rst
                except Exception:
                    pass
                print(f"[Autosave] Loaded checkpoint from '{autosave_nc}'.")
                # Optionally load standardized ocean.nc to override ocean fields
                try:
                    if ocean is not None and int(os.getenv("QD_LOAD_OCEAN", "1")) == 1:
                        ocean_nc = os.path.join("data", "ocean.nc")
                        if os.path.exists(ocean_nc):
                            rst_o = load_ocean(ocean_nc)
                            if rst_o.get("uo") is not None: ocean.uo = rst_o["uo"]
                            if rst_o.get("vo") is not None: ocean.vo = rst_o["vo"]
                            if rst_o.get("eta") is not None: ocean.eta = rst_o["eta"]
                            if rst_o.get("Ts") is not None: ocean.Ts = rst_o["Ts"]
                            print("[Autosave] Ocean state overridden from 'data/ocean.nc'.")
                except Exception as _loe:
                    print(f"[Autosave] ocean.nc load skipped: {_loe}")
                # Load astronomical epoch if present
                try:
                    t_loaded = rst.get("t_seconds", None)
                    if t_loaded is not None:
                        t0_seconds = float(t_loaded)
                except Exception:
                    pass
            except Exception as e:
                print(f"[Autosave] Failed to load '{autosave_nc}': {e}\nApplying banded initialization.")
                apply_banded_initial_ts(grid, gcm, ocean, land_mask)
        else:
            apply_banded_initial_ts(grid, gcm, ocean, land_mask)
        # Load genes autosave JSON (if present) before LAI/species weights
        if int(os.getenv("QD_AUTOSAVE_LOAD", "1")) == 1 and eco is not None:
            try:
                genes_path = os.getenv("QD_ECO_GENES_JSON_PATH") or os.path.join("data", "genes.json")
                if os.path.exists(genes_path) and hasattr(eco, "load_genes_json"):
                    eco.load_genes_json(genes_path)
            except Exception as _egl:
                if int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                    print(f"[Ecology] genes autosave load skipped: {_egl}")
        # Try to load ecology autosave (optional; NetCDF preferred, JSON genes loaded earlier)
        if int(os.getenv("QD_AUTOSAVE_LOAD", "1")) == 1 and eco is not None:
            try:
                eco_path = os.getenv("QD_ECO_AUTOSAVE_PATH") or os.path.join("data", "ecology.nc")
                if os.path.exists(eco_path):
                    ok_eco = False
                    if hasattr(eco, "load_autosave"):
                        ok_eco = bool(eco.load_autosave(eco_path, on_mismatch=os.getenv("QD_ECO_ON_MISMATCH", "fallback")))
                    else:
                        if int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                            print("[Ecology] adapter has no load_autosave; cannot load ecology.nc")
                    if int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                        print(f"[Ecology] autosave load {'OK' if ok_eco else 'failed'} from '{eco_path}'")
            except Exception as _ecl:
                if int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                    print(f"[Ecology] autosave load skipped: {_ecl}")

    # --- Simulation Parameters ---
    dt = int(os.getenv("QD_DT_SECONDS", "300"))  # default 300s
    day_in_seconds = 2 * np.pi / constants.PLANET_OMEGA
    # Duration priority: QD_TOTAL_YEARS -> QD_SIM_DAYS -> default (5 planetary years)
    if os.getenv("QD_TOTAL_YEARS"):
        sim_duration_seconds = float(os.getenv("QD_TOTAL_YEARS")) * orbital_sys.T_planet
    elif os.getenv("QD_SIM_DAYS"):
        sim_duration_seconds = float(os.getenv("QD_SIM_DAYS")) * day_in_seconds
    else:
        sim_duration_seconds = 5 * orbital_sys.T_planet
    
    # --- Physics Parameters (Cloud-Precipitation-Albedo Feedback) ---
    print("Setting physics parameters...")
    D_crit = -1e-7  # Critical convergence for precipitation (s^-1).
    k_precip = 1e5    # Precipitation efficiency.
    alpha_water = 0.1   # Albedo of water/land (fallback when not using base map)
    alpha_ice = 0.6     # Albedo of ice
    alpha_cloud = 0.5   # Average cloud albedo
    # Topography-driven options
    USE_TOPO_ALBEDO = int(os.getenv("QD_USE_TOPO_ALBEDO", "1")) == 1
    OROG_ENABLED = int(os.getenv("QD_OROG", "0")) == 1
    K_OROG = float(os.getenv("QD_OROG_K", "7e-4"))
    # Diagnostics toggle: per-star ISR components (disabled by default)
    PLOT_ISR = int(os.getenv("QD_PLOT_ISR", "0")) == 1

    # --- P019: Lapse & geometry constraints (docs/18, projects/019) ---
    LAPSE_ENABLE = int(os.getenv("QD_LAPSE_ENABLE", "1")) == 1
    GAMMA_KPM = float(os.getenv("QD_LAPSE_K_KPM", "6.5"))     # K/km for air
    GAMMA_S_KPM = float(os.getenv("QD_LAPSE_KS_KPM", os.getenv("QD_LAPSE_K_KPM", "6.5")))  # K/km for surface gate
    LAND_ELEV_MAX_M = float(os.getenv("QD_LAND_ELEV_MAX_M", "10000"))
    POLAR_ICE_THICK_MAX_M = float(os.getenv("QD_POLAR_ICE_THICK_MAX_M", "4500"))
    POLAR_LAT_THRESH = float(os.getenv("QD_POLAR_LAT_THRESH", "60"))
    RHO_SNOW = float(os.getenv("QD_RHO_SNOW", "300"))  # kg/m^3, geometric conversion for snow
    # Glacier mask thresholds（陆地冰相掩膜阈值）
    GLACIER_FRAC = float(os.getenv("QD_GLACIER_FRAC", "0.60"))     # C_snow ≥ 0.60 → 冰盖
    GLACIER_SWE_MM = float(os.getenv("QD_GLACIER_SWE_MM", "50.0")) # 或 SWE ≥ 50 mm → 冰盖
    
    # Optional epoch override (only if no restart/autosave time was loaded)
    if t0_seconds == 0.0:
        try:
            if os.getenv("QD_ORBIT_EPOCH_SECONDS"):
                t0_seconds = float(os.getenv("QD_ORBIT_EPOCH_SECONDS"))
            elif os.getenv("QD_ORBIT_EPOCH_DAYS"):
                t0_seconds = float(os.getenv("QD_ORBIT_EPOCH_DAYS")) * day_in_seconds
        except Exception:
            pass

    time_steps = np.arange(t0_seconds, t0_seconds + sim_duration_seconds, dt)

    # --- Visualization Parameters ---
    output_dir = "output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    # Ensure data/ exists for autosave/restart files
    data_dir = "data"
    os.makedirs(data_dir, exist_ok=True)
    
    plot_every_days = float(os.getenv("QD_PLOT_EVERY_DAYS", "10"))
    plot_interval_seconds = plot_every_days * 24 * 3600
    plot_interval_steps = max(1, int(plot_interval_seconds / dt))
    print(f"Generating a plot every {plot_interval_steps} steps.")
    
    print(f"\n--- Starting Simulation ---")
    print(f"Grid resolution: {grid.n_lat} lat x {grid.n_lon} lon")
    print(f"Time step (dt): {dt} s")
    print(f"Simulation duration: {sim_duration_seconds / day_in_seconds:.1f} planetary days")
    print(f"Total time steps: {len(time_steps)}")

    # --- Autosave hooks (Ctrl-C safe) ---
    AUTOSAVE_ENABLED = int(os.getenv("QD_AUTOSAVE_ENABLE", "1")) == 1
    current_day = [0.0]
    # Initialize autosave day marker from epoch
    try:
        current_day[0] = float(time_steps[0]) / day_in_seconds
    except Exception:
        current_day[0] = 0.0

    def _autosave_hook():
        try:
            # Core and ecology
            save_autosave("data", grid, gcm, ocean, land_mask, W_land, S_snow, eco, day_value=current_day[0])
            # Ocean physical state
            if ocean is not None:
                save_ocean(os.path.join("data", "ocean.nc"), grid, ocean, current_day[0])
            # Phytoplankton: bio/optical JSON and distributions
            if phyto is not None:
                try:
                    phyto.save_bio_json(os.path.join("data", "plankton.json"), current_day[0])
                except Exception as _pb:
                    print(f"[Plankton] bio json save failed: {_pb}")
                try:
                    phyto.save_distribution_nc(os.path.join("data", "plankton.nc"), current_day[0])
                except Exception as _pn:
                    print(f"[Plankton] distribution save failed: {_pn}")
        except Exception as _ae:
            print(f"[Autosave] Save failed: {_ae}")

    def _signal_handler(signum, frame):
        if AUTOSAVE_ENABLED:
            print(f"[Autosave] Caught signal {signum}, saving checkpoint...")
            _autosave_hook()
        # Exit with conventional codes: 130 for SIGINT, 143 for SIGTERM
        try:
            import sys as _sys
            _sys.exit(130 if signum == signal.SIGINT else 143)
        except SystemExit:
            raise

    if AUTOSAVE_ENABLED:
        atexit.register(_autosave_hook)
        try:
            signal.signal(signal.SIGINT, _signal_handler)
            signal.signal(signal.SIGTERM, _signal_handler)
        except Exception:
            pass

    # --- Precipitation accumulation over one planetary day (kg m^-2 ≡ mm/day) ---
    precip_acc_day = np.zeros_like(grid.lat_mesh, dtype=float)
    accum_t_day = 0.0
    precip_day_last = None
    # Ecology viz caches
    last_alpha_ecology_map = None
    last_alpha_banded = None
    alpha_banded_daily = None
    # Bootstrap ecology alpha at t=0 so first panel is not 'not available'
    try:
        if eco is not None:
            insA0, insB0 = forcing.calculate_insolation_components(0.0)
            gcm.isr_A, gcm.isr_B = insA0, insB0
            gcm.isr = insA0 + insB0
            if SUBDAILY_ENABLED and int(os.getenv("QD_ECO_ALBEDO_COUPLE", "1")) == 1:
                alpha0 = eco.step_subdaily(gcm.isr, getattr(gcm, "cloud_cover", 0.0), dt)
                if alpha0 is not None:
                    last_alpha_ecology_map = alpha0.copy()
            # Precompute banded alpha once for initial panel if needed
            try:
                Ab0, w0 = eco.get_surface_albedo_bands()
                if Ab0 is not None and w0 is not None:
                    last_alpha_banded = np.clip(np.nansum(Ab0 * w0[:, None, None], axis=0), 0.0, 1.0)
                    alpha_banded_daily = last_alpha_banded.copy()
            except Exception:
                pass
    except Exception:
        pass
    # Phytoplankton coupling state (P017)
    ocean_mask = (land_mask == 0)
    phyto_next_time = 0.0
    last_alpha_water_scalar = None
    # Diversity scheduling (in planetary days)
    diversity_next_day = [0.0]

    # 2. Time Integration Loop
    try:
        from tqdm import tqdm
        iterator = tqdm(time_steps)
    except ImportError:
        print("tqdm not found, using simple print statements for progress.")
        iterator = time_steps

    # Periodic autosave interval (planetary hours)
    try:
        _hours = float(os.getenv("QD_ECO_AUTOSAVE_EVERY_HOURS", "6"))
        AUTOSAVE_INTERVAL_SECONDS = (_hours * (day_in_seconds / 24.0))
        next_autosave_t = time_steps[0] + AUTOSAVE_INTERVAL_SECONDS
    except Exception:
        AUTOSAVE_INTERVAL_SECONDS = None
        next_autosave_t = None

    for i, t in enumerate(iterator):
        # Periodic autosave tick
        if AUTOSAVE_ENABLED and (AUTOSAVE_INTERVAL_SECONDS is not None) and (next_autosave_t is not None) and (t >= next_autosave_t):
            _autosave_hook()
            next_autosave_t += AUTOSAVE_INTERVAL_SECONDS

        # 1. Physics Step
        # 1) Humidity-aware precipitation（混合方案，使用 P_cond + 动力再分配 + 可选地形强化）
        #    先计算地形增强因子（若有外部 elevation）
        orog_factor = None
        if OROG_ENABLED and (elevation is not None):
            try:
                orog_factor = compute_orographic_factor(grid, elevation, gcm.u, gcm.v, k_orog=K_OROG)
            except Exception as e:
                if i == 0:
                    print(f"[Orog] Disabled due to error: {e}")
        #    使用湿度模块的 P_cond（若存在）作为总量基准；用辐合和地形做空间再分配；全局重标定保持 ⟨P⟩=⟨P_cond⟩
        beta_div = float(os.getenv("QD_P_BETADIV", "0.4"))
        precip = diagnose_precipitation_hybrid(
            gcm, grid, D_crit=D_crit, k_precip=k_precip,
            orog_factor=orog_factor, smooth_sigma=1.0, beta_div=beta_div, renorm=True
        )

        # Accumulate precipitation over one planetary day (kg m^-2 over last day window)
        precip_acc_day += np.nan_to_num(precip) * dt
        accum_t_day += dt
        while accum_t_day >= day_in_seconds:
            precip_day_last = precip_acc_day.copy()
            precip_acc_day[:] = 0.0
            accum_t_day -= day_in_seconds

            # --- Ecology M2 daily step: update LAI from soil water proxy ∈[0,1] ---
            if eco is not None:
                try:
                    soil_cap_env = os.getenv("QD_ECO_SOIL_WATER_CAP")
                    if soil_cap_env is not None:
                        try:
                            soil_cap = float(soil_cap_env)
                        except Exception:
                            soil_cap = 50.0  # mm ~ kg/m^2
                    else:
                        soil_cap = 50.0
                    # W_land 单位≈kg/m^2 ≡ mm；归一化后截断到 [0,1]
                    soil_idx = np.clip(W_land / max(1e-6, soil_cap), 0.0, 1.0)
                    # 冰盖处植物无法生存：将冰盖处土壤指数置 0，抑制 LAI 增长
                    try:
                        if hasattr(gcm, "glacier_mask_last"):
                            soil_idx = soil_idx * (~gcm.glacier_mask_last)
                    except Exception:
                        pass
                    eco.step_daily(soil_idx)
                    # 保险起见，若生态内核没有处理，强制将冰盖处 LAI 置 0
                    try:
                        if getattr(eco, "pop", None) is not None and hasattr(eco.pop, "LAI") and hasattr(gcm, "glacier_mask_last"):
                            mgl = gcm.glacier_mask_last
                            eco.pop.LAI = np.where(mgl, 0.0, eco.pop.LAI)
                    except Exception:
                        pass
                    # Individual pool daily aggregation -> adjust species splits for sampled cells
                    if 'indiv' in locals() and indiv is not None:
                        try:
                            # Exclude glacier pixels from individual sampling to save compute
                            try:
                                glacier_mask = getattr(gcm, "glacier_mask_last", None)
                                if glacier_mask is not None:
                                    active_mask = ((land_mask == 1) & (~glacier_mask))
                                    if hasattr(indiv, "set_active_mask"):
                                        indiv.set_active_mask(active_mask)
                                    elif hasattr(indiv, "land_mask"):
                                        indiv.land_mask = active_mask.astype(int)
                            except Exception:
                                pass
                            indiv.step_daily(eco, soil_idx, Ts_map=gcm.T_s, day_length_hours=24.0)
                        except Exception as _ied:
                            if int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                                print(f"[EcoIndiv] daily step skipped: {_ied}")
                    # Daily banded albedo update (move heavy call to daily boundary)
                    if int(os.getenv("QD_ECO_BANDS_COUPLE", "0")) == 1:
                        try:
                            Abands, w_b = eco.get_surface_albedo_bands()
                        except Exception:
                            Abands, w_b = (None, None)
                        if Abands is not None and w_b is not None:
                            alpha_banded_daily = np.nansum(Abands * w_b[:, None, None], axis=0)
                            last_alpha_banded = np.clip(alpha_banded_daily, 0.0, 1.0).copy()
                            if int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                                print(f"[Ecology] Daily banded alpha updated from {Abands.shape[0]} bands.")
                            # Optional: export genes snapshot daily (M3)
                            try:
                                if int(os.getenv("QD_ECO_GENES_EXPORT", "0")) == 1 and hasattr(eco, "export_genes"):
                                    # Estimate day value at the completed daily boundary
                                    day_value_est = float((t - accum_t_day) / day_in_seconds)
                                    eco.export_genes(output_dir, day_value_est)
                            except Exception:
                                pass
                    # Fallback/export genes even if band coupling disabled or failed
                    if int(os.getenv("QD_ECO_GENES_EXPORT", "0")) == 1 and hasattr(eco, "export_genes"):
                        try:
                            day_value_est = float((t - accum_t_day) / day_in_seconds)
                            eco.export_genes(output_dir, day_value_est)
                        except Exception:
                            pass
                except Exception as _ed:
                    if int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                        print(f"[Ecology] daily step skipped: {_ed}")
        
        # 1b) Convert precipitation to cloud fraction via smooth saturating relation
        if np.any(precip > 0):
            P_ref_env = os.getenv("QD_PREF")
            if P_ref_env:
                P_ref = float(P_ref_env)
            else:
                P_pos = precip[precip > 0]
                P_ref = float(np.median(P_pos)) if P_pos.size > 0 else 1e-6
        else:
            P_ref = 1e-6
        C_from_P = cloud_from_precip(
            precip,
            C_max=float(os.getenv("QD_CMAX", "0.95")),
            P_ref=P_ref,
            smooth_sigma=1.0
        )

        # 1c) Additional prognostic cloud source from thermodynamic/dynamic proxies
        cloud_source = parameterize_cloud_cover(gcm, grid, land_mask)

        # 1d) Integrate cloud cover with blending:
        #     - retain memory (previous cloud_cover)
        #     - enforce strong coupling to precipitation-driven cloud
        #     - include additional source term
        tendency = cloud_source * (dt / (6 * 3600))
        # Blending weights (tunable via env: QD_W_MEM, QD_W_P, QD_W_SRC), normalized to sum to 1
        W_MEM = float(os.getenv("QD_W_MEM", "0.4"))
        W_P   = float(os.getenv("QD_W_P", "0.4"))
        W_SRC = float(os.getenv("QD_W_SRC", "0.2"))
        W_sum = W_MEM + W_P + W_SRC
        if W_sum <= 0:
            W_MEM, W_P, W_SRC, W_sum = 0.5, 0.4, 0.1, 1.0
        W_MEM /= W_sum; W_P /= W_sum; W_SRC /= W_sum

        gcm.cloud_cover = (
            W_MEM * gcm.cloud_cover +
            W_P   * C_from_P +
            W_SRC * np.clip(gcm.cloud_cover + tendency, 0.0, 1.0)
        )
        if i == 0:
            print(f"[CloudBlend] Weights: W_MEM={W_MEM:.2f}, W_P={W_P:.2f}, W_SRC={W_SRC:.2f}; P_ref={P_ref:.3e}, C_max={float(os.getenv('QD_CMAX', '0.95')):.2f}")
        # Physical consistency: where it rains, clouds must be present (event-floor from precip)
        C_floor = float(os.getenv("QD_CLOUD_FROM_P_FLOOR", "0.8"))  # 0..1
        if C_floor > 0.0:
            cloud_floor = np.clip(C_floor * C_from_P, 0.0, 1.0)
            gcm.cloud_cover = np.maximum(gcm.cloud_cover, cloud_floor)

        gcm.cloud_cover = np.clip(gcm.cloud_cover, 0.0, 1.0)

        # Optional: advect cloud as a tracer to reduce stickiness (experimental)
        if int(os.getenv("QD_CLOUD_ADVECT", "1")) == 1:
            try:
                adv_alpha = float(os.getenv("QD_CLOUD_ADV_ALPHA", "0.7"))
            except Exception:
                adv_alpha = 0.7
            try:
                cloud_adv = _advect_scalar_periodic(gcm.cloud_cover, gcm.u, gcm.v, dt, grid)
                # Optional smoothing after advection
                try:
                    sig = float(os.getenv("QD_CLOUD_SMOOTH_SIGMA", "0.2"))
                except Exception:
                    sig = 0.2
                if sig > 0.0:
                    try:
                        from scipy.ndimage import gaussian_filter
                        cloud_adv = gaussian_filter(cloud_adv, sigma=sig, mode="wrap")
                    except Exception:
                        pass
                gcm.cloud_cover = np.clip((1.0 - adv_alpha) * gcm.cloud_cover + adv_alpha * cloud_adv, 0.0, 1.0)
                if i == 0:
                    print(f"[CloudAdvect] enabled: ADV_ALPHA={adv_alpha:.2f}, SMOOTH_SIGMA={sig}")
            except Exception as _cae:
                if i == 0:
                    print(f"[CloudAdvect] skipped: {_cae}")

        # 2) Forcing Step (moved earlier for ecology coupling): compute insolation components
        insA, insB = forcing.calculate_insolation_components(t)
        gcm.isr_A, gcm.isr_B = insA, insB
        gcm.isr = insA + insB

        # --- P019: Lapse-adjusted temperatures and smooth phase split + provisional snowpack update ---
        # Compute atmospheric temperature proxy (consistent with dynamics core simplification)
        g_const = 9.81
        T_a_proxy = 288.0 + (g_const / 1004.0) * gcm.h

        # Geometry: bedrock elevation and snow geometric thickness (land only)
        H_bedrock = elevation if (('elevation' in locals()) and (elevation is not None)) else np.zeros_like(gcm.T_s)
        # SWE unit here is kg/m^2 ≡ mm of water; convert to geometric snow thickness by rho_snow
        h_snow_geom = np.where(land_mask == 1, np.maximum(S_snow, 0.0) / max(RHO_SNOW, 1e-6), 0.0)
        # Polar cap thickness limit
        polar_mask = (np.abs(grid.lat_mesh) >= POLAR_LAT_THRESH)
        h_ice_eff = np.where(polar_mask, np.minimum(h_snow_geom, POLAR_ICE_THICK_MAX_M), h_snow_geom)
        # Effective elevation (cap to LAND_ELEV_MAX_M)
        H_eff = np.minimum(H_bedrock + h_ice_eff, LAND_ELEV_MAX_M)

        # Lapse-adjusted temps
        if LAPSE_ENABLE:
            T_hat_a = T_a_proxy - GAMMA_KPM * (H_eff / 1000.0)
            T_hat_s = gcm.T_s - GAMMA_S_KPM * (H_eff / 1000.0)
        else:
            T_hat_a = T_a_proxy
            T_hat_s = gcm.T_s

        # Smooth phase split using T_hat_a (Sigmoid)
        # Use hydrology params loaded earlier for thresholds/ΔT
        try:
            dT_half = float(getattr(hydro_params, "snow_t_band_K", 1.5))
        except Exception:
            dT_half = 1.5
        P_rain_p019, P_snow_p019, f_snow_p019 = partition_precip_phase_smooth(
            P_flux=gcm.isr*0.0 + (0.0 if 'precip' not in locals() else precip),  # safe shape; replaced below if precip exists
            T_hat_a=T_hat_a,
            T_thresh=hydro_params.snow_thresh_K,
            dT_half_K=dT_half
        )
        # Replace with actual precip field (ensure arrays)
        P_flux_arr = precip if not np.isscalar(precip) else np.full_like(gcm.T_s, float(precip))
        P_rain_p019, P_snow_p019, f_snow_p019 = partition_precip_phase_smooth(
            P_flux=P_flux_arr, T_hat_a=T_hat_a, T_thresh=hydro_params.snow_thresh_K, dT_half_K=dT_half
        )

        # Provisional snowpack update (do not commit S_snow yet; reuse later to keep single-step update)
        if getattr(hydro_params, "swe_enable", True):
            land = (land_mask == 1)
            P_snow_land_p019 = P_snow_p019 * land
            S_snow_next_p019, melt_flux_land_p019, C_snow_map, alpha_snow_map_p019 = snowpack_step(
                S_snow=S_snow, P_snow_land=P_snow_land_p019, T_hat_a=T_hat_a, params=hydro_params, dt=dt
            )
            # --- Glacier mask（陆地冰相掩膜）---
            glacier_mask = (land_mask == 1) & ((C_snow_map >= GLACIER_FRAC) | (S_snow_next_p019 >= GLACIER_SWE_MM))
            # 雨落在冰盖上：按“沉积”处理，转入 SWE（冻结沉积），不进入地表桶
            try:
                P_rain_land_glacier = (P_rain_p019 * (land_mask == 1)) * glacier_mask
                if np.any(P_rain_land_glacier):
                    S_snow_next_p019 = S_snow_next_p019 + P_rain_land_glacier * dt
            except Exception:
                pass
            # Expose for TrueColor & downstream modules
            try:
                gcm.C_snow_map_last = C_snow_map
                gcm.glacier_mask_last = glacier_mask
            except Exception:
                pass
        else:
            C_snow_map = np.zeros_like(gcm.T_s, dtype=float)
            # Expose for TrueColor even when SWE disabled (renders nothing by default)
            try:
                gcm.C_snow_map_last = C_snow_map
                gcm.glacier_mask_last = (land_mask == 1) & (C_snow_map >= GLACIER_FRAC)
            except Exception:
                pass
            alpha_snow_map_p019 = np.full_like(gcm.T_s, float(os.getenv("QD_SNOW_ALBEDO_FRESH", "0.70")), dtype=float)
            S_snow_next_p019 = S_snow.copy()
            melt_flux_land_p019 = np.zeros_like(gcm.T_s, dtype=float)

        # 2a) Eco individual pool subdaily step (vectorized) - accumulate per-individual energy/day
        if 'indiv' in locals() and indiv is not None:
            try:
                # Build subdaily soil index proxy from W_land (kg/m^2 ≈ mm) with cap
                soil_cap_env = os.getenv("QD_ECO_SOIL_WATER_CAP")
                if soil_cap_env is not None:
                    try:
                        soil_cap_sub = float(soil_cap_env)
                    except Exception:
                        soil_cap_sub = 50.0
                else:
                    soil_cap_sub = 50.0
                soil_idx_sub = np.clip(W_land / max(1e-6, soil_cap_sub), 0.0, 1.0)
                # Exclude glacier pixels from individual sampling to save compute
                try:
                    glacier_mask = getattr(gcm, "glacier_mask_last", None)
                    if glacier_mask is not None:
                        active_mask = ((land_mask == 1) & (~glacier_mask))
                        if hasattr(indiv, "set_active_mask"):
                            indiv.set_active_mask(active_mask)
                        elif hasattr(indiv, "land_mask"):
                            # Some implementations store mask as int map
                            indiv.land_mask = active_mask.astype(int)
                except Exception:
                    pass
                indiv.try_substep(gcm.isr_A, gcm.isr_B, eco, soil_idx_sub, dt, day_in_seconds)
            except Exception as _ies:
                if int(os.getenv("QD_ECO_DIAG", "1")) == 1 and i == 0:
                    print(f"[EcoIndiv] substep skipped: {_ies}")

        # P017: Daily phytoplankton step (after ISR available)
        if phyto is not None and PHYTO_ENABLED:
            if t >= phyto_next_time:
                try:
                    T_w = ocean.Ts if ocean is not None else gcm.T_s
                    _ab, alpha_scalar = phyto.step_daily(gcm.isr_A, gcm.isr_B, T_w, dt_days=1.0)
                    last_alpha_water_scalar = alpha_scalar
                except Exception as _pe:
                    if int(os.getenv("QD_PHYTO_DIAG", "1")) == 1:
                        print(f"[Phyto] daily step skipped: {_pe}")
                phyto_next_time = t + day_in_seconds

        # 2a) Radiative inputs and sea-ice fraction for albedo synthesis
        H_ice_ref = float(os.getenv("QD_HICE_REF", "0.5"))  # m, e-folding thickness for ice optical effect
        ice_frac = 1.0 - np.exp(-np.maximum(gcm.h_ice, 0.0) / max(1e-6, H_ice_ref))
        cloud_for_rad = getattr(gcm, "cloud_eff_last", gcm.cloud_cover)

        # 2b) Ecology M1: land-only scalar surface alpha (optional hourly coupling)
        base_input = None
        if USE_TOPO_ALBEDO:
            base_input = base_albedo_map.copy()
        else:
            base_input = np.full_like(gcm.T_s, float(alpha_water))

        if eco is not None and SUBDAILY_ENABLED and int(os.getenv("QD_ECO_ALBEDO_COUPLE", "1")) == 1:
            try:
                alpha_map = eco.step_subdaily(gcm.isr, cloud_for_rad, dt)
            except Exception as _ee:
                if i == 0:
                    print(f"[Ecology] subdaily step skipped due to error: {_ee}")
                alpha_map = None
            # Decide which alpha to apply this step: new one or last cached
            if alpha_map is None and last_alpha_ecology_map is not None:
                alpha_apply = last_alpha_ecology_map
            else:
                alpha_apply = alpha_map
            if alpha_apply is not None:
                try:
                    W_LAI = float(os.getenv("QD_ECO_LAI_ALBEDO_WEIGHT", "1.0"))
                except Exception:
                    W_LAI = 1.0
                # Blend ecology alpha into land base albedo
                land = (land_mask == 1)
                # 冰盖区域不混入生态反照率（保持雪/冰主导）
                try:
                    glacier_mask = getattr(gcm, "glacier_mask_last", np.zeros_like(land_mask, dtype=bool))
                except Exception:
                    glacier_mask = np.zeros_like(land_mask, dtype=bool)
                m = land & (~glacier_mask) & np.isfinite(alpha_apply)
                base_input[m] = (1.0 - W_LAI) * base_input[m] + W_LAI * alpha_apply[m]
                # cache for viz if new map arrived
                if alpha_map is not None:
                    last_alpha_ecology_map = alpha_map.copy()
                # First-step confirmation

        # 2b.1) M3b（可选）：若启用带化耦合，则用 A_b^surface 做一次加权降维混入 base_input
        if eco is not None and int(os.getenv("QD_ECO_BANDS_COUPLE", "0")) == 1:
            # Use daily-cached banded alpha to avoid heavy call every physics step
            alpha_banded_step = alpha_banded_daily if alpha_banded_daily is not None else last_alpha_banded
            if alpha_banded_step is not None:
                land = (land_mask == 1)
                m2 = land & np.isfinite(alpha_banded_step)
                base_input[m2] = np.clip(alpha_banded_step[m2], 0.0, 1.0)
                if i == 0 and int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                    try:
                        ab = alpha_banded_step[land]
                        print(f"[Ecology] bands-couple: alpha_banded(min/mean/max)={np.nanmin(ab):.3f}/{np.nanmean(ab):.3f}/{np.nanmax(ab):.3f}")
                    except Exception:
                        pass

        # Apply phytoplankton ocean-color feedback (override ocean base albedo)
        if PHYTO_ENABLED and PHYTO_COUPLE and (last_alpha_water_scalar is not None):
            try:
                m_o = ocean_mask & np.isfinite(last_alpha_water_scalar)
                base_input[m_o] = np.clip(last_alpha_water_scalar[m_o], 0.0, 1.0)
            except Exception as _pc:
                if int(os.getenv("QD_PHYTO_DIAG", "1")) == 1 and i == 0:
                    print(f"[Phyto] coupling skipped: {_pc}")

        # --- P019: Blend snow cover into land surface base albedo ---
        try:
            if getattr(hydro_params, "swe_enable", True):
                land = (land_mask == 1)
                # α_surface_eff = α_base·(1 − C_snow) + α_snow·C_snow
                base_input[land] = np.clip(
                    (1.0 - C_snow_map[land]) * base_input[land] + C_snow_map[land] * alpha_snow_map_p019[land],
                    0.0, 1.0
                )
        except Exception as _sb:
            if i == 0:
                print(f"[P019] snow-albedo blend skipped: {_sb}")

        # 2c) Final dynamic albedo (surface base + cloud + ice)
        albedo = calculate_dynamic_albedo(
            cloud_for_rad, gcm.T_s, base_input, alpha_ice, alpha_cloud, land_mask=land_mask, ice_frac=ice_frac
        )

        # Global energy budget diagnostic (TOA/SFC/ATM), independent of autotune
        try:
            if int(os.getenv("QD_ENERGY_DIAG", "1")) == 1 and (i % 200 == 0):
                # Shortwave split and reflection with current albedo/cloud
                SW_atm_dbg, SW_sfc_dbg, R_dbg = energy.shortwave_radiation(gcm.isr, albedo, cloud_for_rad, eparams)
                # Longwave using chosen scheme
                g_const = 9.81
                T_a_dbg = 288.0 + (g_const / 1004.0) * gcm.h
                use_lw_v2 = int(os.getenv("QD_LW_V2", "1")) == 1
                if use_lw_v2:
                    eps_sfc_dbg = energy.surface_emissivity_map(land_mask, ice_frac)
                    _LW_atm_dbg, LW_sfc_dbg, OLR_dbg, DLR_dbg, _eps_dbg = energy.longwave_radiation_v2(
                        gcm.T_s, T_a_dbg, cloud_for_rad, eps_sfc_dbg, eparams
                    )
                else:
                    _LW_atm_dbg, LW_sfc_dbg, OLR_dbg, DLR_dbg, _eps_dbg = energy.longwave_radiation(
                        gcm.T_s, T_a_dbg, cloud_for_rad, eparams
                    )
                # SH/LH
                C_H = float(os.getenv("QD_CH", "1.5e-3"))
                cp_air = float(os.getenv("QD_CP_A", "1004.0"))
                rho_air = float(getattr(getattr(gcm, "hum_params", None), "rho_a", 1.2))
                B_land = float(os.getenv("QD_BOWEN_LAND", "0.7"))
                B_ocean = float(os.getenv("QD_BOWEN_OCEAN", "0.3"))
                SH_dbg, _ = energy.boundary_layer_fluxes(
                    gcm.T_s, T_a_dbg, gcm.u, gcm.v, land_mask,
                    C_H=C_H, rho=rho_air, c_p=cp_air, B_land=B_land, B_ocean=B_ocean
                )
                LH_dbg = getattr(gcm, "LH_last", 0.0)
                if np.isscalar(LH_dbg):
                    LH_dbg = np.full_like(gcm.T_s, float(LH_dbg))
                # Compute and print diagnostics
                diagE = energy.compute_energy_diagnostics(
                    grid.lat_mesh, gcm.isr, R_dbg, OLR_dbg, SW_sfc_dbg, LW_sfc_dbg, SH_dbg, LH_dbg
                )
                print(f"[EnergyDiag] TOA_net={diagE['TOA_net']:.2f} W/m^2 | "
                      f"SFC_net={diagE['SFC_net']:.2f} | ATM_net={diagE['ATM_net']:.2f} | "
                      f"<Ts>={diagE.get('Ts_mean', float(np.nanmean(gcm.T_s))):.2f} K")
        except Exception as _edbg:
            if i == 0:
                print(f"[EnergyDiag] skipped: {_edbg}")

        # 2d) Equilibrium temp with updated albedo
        Teq = forcing.calculate_equilibrium_temp(t, albedo)

        # 3. Dynamics Step
        gcm.time_step(Teq, dt)

        # 3a. Ocean step (P011): wind-driven currents + SST advection + optional Q_net coupling
        if ocean is not None:
            try:
                # Prepare inputs
                # Sea-ice mask: treat any positive thickness as ice-covered
                ice_mask = (getattr(gcm, "h_ice", np.zeros_like(gcm.T_s)) > 0.0)
                # Cloud optical field for radiation consistency if available
                cloud_eff = getattr(gcm, "cloud_eff_last", gcm.cloud_cover)
                # Energy params for radiation calculation (persistent; may be auto-tuned)
                # eparams is initialized once before the loop and optionally auto-tuned
                # Shortwave components at surface (use same albedo as current step)
                SW_atm, SW_sfc, _R = energy.shortwave_radiation(gcm.isr, albedo, cloud_eff, eparams)
                # Atmospheric temperature proxy (consistent with dynamics core simplification)
                g_const = 9.81
                T_a = 288.0 + (g_const / 1004.0) * gcm.h
                # Longwave at surface
                use_lw_v2 = int(os.getenv("QD_LW_V2", "1")) == 1
                H_ice_ref = float(os.getenv("QD_HICE_REF", "0.5"))
                ice_frac = 1.0 - np.exp(-np.maximum(getattr(gcm, "h_ice", np.zeros_like(gcm.T_s)), 0.0) / max(1e-6, H_ice_ref))
                if use_lw_v2:
                    eps_sfc_map = energy.surface_emissivity_map(land_mask, ice_frac)
                    _LW_atm, LW_sfc, _OLR, _DLR, _eps = energy.longwave_radiation_v2(
                        gcm.T_s, T_a, cloud_eff, eps_sfc_map, eparams
                    )
                else:
                    _LW_atm, LW_sfc, _OLR, _DLR, _eps = energy.longwave_radiation(
                        gcm.T_s, T_a, cloud_eff, eparams
                    )
                # Sensible heat flux (SH)
                C_H = float(os.getenv("QD_CH", "1.5e-3"))
                cp_air = float(os.getenv("QD_CP_A", "1004.0"))
                rho_air = float(getattr(getattr(gcm, "hum_params", None), "rho_a", 1.2))
                B_land = float(os.getenv("QD_BOWEN_LAND", "0.7"))
                B_ocean = float(os.getenv("QD_BOWEN_OCEAN", "0.3"))
                SH_arr, _LH_bowen = energy.boundary_layer_fluxes(
                    gcm.T_s, T_a, gcm.u, gcm.v, land_mask,
                    C_H=C_H, rho=rho_air, c_p=cp_air, B_land=B_land, B_ocean=B_ocean
                )
                # Latent heat (LH) from humidity module diagnostics (already computed in dynamics)
                LH_arr = getattr(gcm, "LH_last", 0.0)
                if np.isscalar(LH_arr):
                    LH_arr = np.full_like(gcm.T_s, float(LH_arr))
                # Net heat into surface (W/m^2)
                Q_net = SW_sfc - LW_sfc - SH_arr - LH_arr

                # Optional greenhouse autotuning using global energy diagnostics
                if AUTOTUNE and (i % TUNE_EVERY == 0):
                    diagE = energy.compute_energy_diagnostics(
                        grid.lat_mesh, gcm.isr, _R, _OLR, SW_sfc, LW_sfc, SH_arr, LH_arr
                    )
                    eparams = energy.autotune_greenhouse_params(eparams, diagE)

                # Advance ocean
                ocean.step(dt, gcm.u, gcm.v, Q_net=Q_net, ice_mask=ice_mask)

                # Inject SST back into atmospheric surface temperature over open ocean (no ice)
                ocean_open = (land_mask == 0) & (~ice_mask)
                gcm.T_s = np.where(ocean_open, ocean.Ts, gcm.T_s)

                # P017 M3: Advect phytoplankton by updated ocean currents (optional)
                if phyto is not None and PHYTO_ENABLED and PHYTO_ADVECTION:
                    try:
                        phyto.advect_diffuse(ocean.uo, ocean.vo, dt)
                    except Exception as _pae:
                        if int(os.getenv("QD_PHYTO_DIAG", "1")) == 1 and i == 0:
                            print(f"[Phyto] advection skipped: {_pae}")

                # Optional diagnostics
                if int(os.getenv("QD_OCEAN_DIAG", "1")) == 1 and (i % 200 == 0):
                    od = ocean.diagnostics()
                    print(f"[OceanDiag] KE_mean={od['KE_mean']:.3e} m2/s2 | Umax={od['U_max']:.2f} m/s | "
                          f"eta[{od['eta_min']:.3f},{od['eta_max']:.3f}] m | cfl/sqrt(gH)/dx={od['cfl_per_s']:.3e} s^-1")
            except Exception as e:
                if i == 0:
                    print(f"[Ocean] step skipped due to error: {e}")

        # 3b. Humidity diagnostics (P008): global means of E, P_cond, LH, LH_release
        try:
            if getattr(gcm, "hum_params", None) is not None and getattr(gcm.hum_params, "diag", False):
                if i % 200 == 0:
                    w = np.maximum(np.cos(np.deg2rad(grid.lat_mesh)), 0.0)
                    wsum = np.sum(w) + 1e-15
                    def wmean(x):
                        return float(np.sum(x * w) / wsum)
                    E_mean = wmean(getattr(gcm, "E_flux_last", 0.0))
                    Pcond_mean = wmean(getattr(gcm, "P_cond_flux_last", 0.0))
                    LH_mean = wmean(getattr(gcm, "LH_last", 0.0))
                    LHrel_mean = wmean(getattr(gcm, "LH_release_last", 0.0))
                    print(f"[HumidityDiag] ⟨E⟩={E_mean:.3e} kg/m^2/s | ⟨P_cond⟩={Pcond_mean:.3e} kg/m^2/s | "
                          f"⟨LH⟩={LH_mean:.2f} W/m^2 | ⟨LH_release⟩={LHrel_mean:.2f} W/m^2")
        except Exception:
            pass

        # 3c. Hydrology step (P009): E–P–R, snow and water-closure diagnostics
        try:
            # Fluxes for hydrology (kg m^-2 s^-1)
            # Use diagnosed precipitation (hybrid precip) for hydrology/routing instead of P_cond,
            # otherwise P=0 will prevent buckets from filling and R=0 forever.
            E_flux = getattr(gcm, "E_flux_last", 0.0)
            P_flux = precip

            # Ensure array shape
            if np.isscalar(E_flux):
                E_flux = np.full_like(gcm.T_s, float(E_flux))
            if np.isscalar(P_flux):
                P_flux = np.full_like(gcm.T_s, float(P_flux))

            # Phase partition & snowpack update (P019) — reuse provisional values computed before albedo
            land = (land_mask == 1)
            if 'P_rain_p019' in locals():
                P_rain = P_rain_p019
                P_snow = P_snow_p019
                P_rain_land = P_rain * land
                P_snow_land = P_snow * land
                E_land = E_flux * land
                # Commit the provisional snowpack update
                S_snow = S_snow_next_p019
                melt_flux_land = melt_flux_land_p019
            else:
                # Fallback legacy (no lapse/sigmoid), threshold on T_s
                P_rain, P_snow = partition_precip_phase(P_flux, gcm.T_s, T_thresh=hydro_params.snow_thresh_K)
                P_rain_land = P_rain * land
                P_snow_land = P_snow * land
                E_land = E_flux * land
                S_snow, melt_flux_land = snow_step(S_snow, P_snow_land, gcm.T_s, hydro_params, dt)

            # Land bucket update: 冰盖处不入桶；冰下管网仅接收融水
            try:
                glacier_mask = getattr(gcm, "glacier_mask_last", np.zeros_like(land_mask, dtype=bool))
            except Exception:
                glacier_mask = np.zeros_like(land_mask, dtype=bool)
            non_glacier = (land_mask == 1) & (~glacier_mask)

            # 非冰盖：入桶（雨 + 融水），蒸发按原值；冰盖：不入桶（雨已在上文沉积进 SWE）
            P_in_land_non_glacier = (P_rain_land + melt_flux_land) * non_glacier
            E_land_non_glacier = E_land * non_glacier

            W_land, R_flux_land_bucket = update_land_bucket(W_land, P_in_land_non_glacier, E_land_non_glacier, hydro_params, dt)

            # 冰盖融水（kg m^-2 s^-1）直接作为“下游流量源”进入路由（不通过桶）
            R_flux_glacier_melt = melt_flux_land * glacier_mask

            # 合成对路由的总陆面出流
            R_flux_land_total = R_flux_land_bucket + R_flux_glacier_melt

            # P014: route runoff along offline network (if enabled)
            if 'routing' in locals() and routing is not None:
                try:
                    routing.step(R_land_flux=R_flux_land_total, dt_seconds=dt, precip_flux=P_flux, evap_flux=E_flux)
                except Exception as _re:
                    if i == 0:
                        print(f"[HydroRouting] step skipped due to error: {_re}")

            # Diagnostics: global water closure (area-weighted)
            if getattr(hydro_params, "diag", True) and (i % 200 == 0):
                # Time since previous diagnostic
                t_now = (i * dt)
                dt_since_prev = None if _hydro_prev_time is None else (t_now - _hydro_prev_time)

                # Required densities/heights from modules
                rho_a = float(getattr(gcm.hum_params, "rho_a", 1.2))
                h_mbl = float(getattr(gcm.hum_params, "h_mbl", 800.0))
                rho_i = float(getattr(gcm, "rho_i", 917.0))

                diag_h2o = diagnose_water_closure(
                    lat_mesh=grid.lat_mesh,
                    q=getattr(gcm, "q", np.zeros_like(gcm.T_s)),
                    rho_a=rho_a,
                    h_mbl=h_mbl,
                    h_ice=getattr(gcm, "h_ice", np.zeros_like(gcm.T_s)),
                    rho_i=rho_i,
                    W_land=W_land,
                    S_snow=S_snow,
                    E_flux=E_flux,
                    P_flux=P_flux,        # use diagnosed precip (hybrid), not P_cond
                    R_flux=R_flux_land_total,   # runoff only from land (bucket + glacier melt)
                    dt_since_prev=dt_since_prev,
                    prev_total=_hydro_prev_total
                )

                # Print concise diagnostics
                msg = (f"[WaterDiag] ⟨E⟩={diag_h2o['E_mean']:.3e} kg/m^2/s | "
                       f"⟨P⟩={diag_h2o['P_mean']:.3e} | ⟨R⟩={diag_h2o['R_mean']:.3e} | "
                       f"⟨CWV⟩={diag_h2o['CWV_mean']:.3e} kg/m^2 | ⟨ICE⟩={diag_h2o['ICE_mean']:.3e} | "
                       f"⟨W_land⟩={diag_h2o['W_land_mean']:.3e} | ⟨S_snow⟩={diag_h2o['S_snow_mean']:.3e}")
                if "closure_residual" in diag_h2o and "d/dt_total_mean" in diag_h2o:
                    msg += (f" | d/dt Σ={diag_h2o['d/dt_total_mean']:.3e} vs (E−P−R) -> "
                            f"residual={diag_h2o['closure_residual']:.3e}")
                print(msg)
                # Optional routing diagnostics
                if 'routing' in locals() and routing is not None:
                    rd = routing.diagnostics()
                    try:
                        max_flow = float(np.nanmax(rd["flow_accum_kgps"]))
                    except Exception:
                        max_flow = 0.0
                    print(f"[HydroRoutingDiag] ocean_inflow={rd['ocean_inflow_kgps']:.3e} kg/s | "
                          f"mass_error={rd['mass_closure_error_kg']:.3e} kg | "
                          f"max_flow={max_flow:.3e} kg/s")

                # Update prev totals/time
                _hydro_prev_total = diag_h2o["total_reservoir_mean"]
                _hydro_prev_time = t_now
        except Exception as _e:
            if i == 0:
                print(f"[Hydrology] step skipped due to error: {_e}")

        # 4. (Optional) Print diagnostics and generate plots
        t_days = t / day_in_seconds

        # Diversity diagnostics (every DIVERSITY_EVERY_DAYS)
        if 'DIVERSITY_ENABLED' in locals() and DIVERSITY_ENABLED and eco_diversity is not None and eco is not None and getattr(eco, "pop", None) is not None and getattr(eco.pop, "LAI_layers_SK", None) is not None:
            if t_days >= diversity_next_day[0]:
                try:
                    eco_diversity.compute_and_plot(grid, eco, land_mask, t_days, output_dir)
                    diversity_next_day[0] = t_days + DIVERSITY_EVERY_DAYS
                except Exception as _de:
                    if int(os.getenv("QD_ECO_DIAG", "1")) == 1 and i == 0:
                        print(f"[Diversity] diagnostics skipped: {_de}")

        # Update autosave day marker
        current_day[0] = t_days
        if i % 100 == 0 and i > 0:
            if not isinstance(iterator, type(time_steps)):
                iterator.set_description(
                    f"t={t_days:.1f}d | "
                    f"max|u|={np.max(np.abs(gcm.u)):.2f} m/s | "
                    f"max|v|={np.max(np.abs(gcm.v)):.2f} m/s"
                )
        
        if i % plot_interval_steps == 0:
            # Plot instantaneous precipitation rate (kg m^-2 s^-1 ⇒ mm/day)
            plot_state(grid, gcm, land_mask, precip, gcm.cloud_cover, albedo, t_days, output_dir, ocean=ocean, routing=routing)
            plot_true_color(grid, gcm, land_mask, t_days, output_dir, routing=routing, eco=eco, phyto=phyto)
            # Optional: plankton species maps (species 0/1) when enabled
            if phyto is not None and int(os.getenv("QD_PLOT_PHYTO", "1")) == 1:
                try:
                    plot_plankton_species(grid, phyto, land_mask, t_days, output_dir)
                except Exception as _pp:
                    if int(os.getenv("QD_PHYTO_DIAG", "1")) == 1:
                        print(f"[PlanktonViz] skipped: {_pp}")
            # Ecology panel: always output (with placeholders if eco is None)
            try:
                lai_map = None
                ch_map = None
                sd_maps = None
                if eco is not None and getattr(eco, "pop", None) is not None:
                    try:
                        lai_map = eco.pop.LAI
                    except Exception:
                        lai_map = None
                    # M2.5: canopy height and species density maps (optional diagnostics)
                    try:
                        ch_map = eco.pop.canopy_height_map()
                    except Exception:
                        ch_map = None
                    try:
                        sd_maps = eco.pop.species_density_maps()
                    except Exception:
                        sd_maps = None

                # Ensure banded alpha available for panel (compute on the fly if not cached)
                if last_alpha_banded is None and eco is not None:
                    try:
                        Abands_now, w_b_now = eco.get_surface_albedo_bands()
                    except Exception:
                        Abands_now, w_b_now = (None, None)
                    if Abands_now is not None and w_b_now is not None:
                        last_alpha_banded = np.clip(np.nansum(Abands_now * w_b_now[:, None, None], axis=0), 0.0, 1.0)
                        if int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                            print(f"[Ecology] Panel fallback: banded alpha computed from {Abands_now.shape[0]} bands.")
                band_ok = (last_alpha_banded is not None)
                # Ecology panel: gated by QD_ECO_PLOT (default off)
                if int(os.getenv("QD_ECO_PLOT", "1")) == 1:
                    plot_ecology(
                        grid, land_mask, t_days, output_dir,
                        lai=lai_map,
                        alpha_ecology=last_alpha_ecology_map,
                        alpha_banded=(last_alpha_banded if band_ok else None),
                        canopy_height=ch_map,
                        species_density=sd_maps
                    )
                # Optional: auto-open ecology panel on macOS at first plot
                try:
                    if sys.platform == "darwin" and int(os.getenv("QD_ECO_PLOT", "1")) == 1 and int(os.getenv("QD_ECO_OPEN", "0")) == 1 and i == 0:
                        panel_path = os.path.join(output_dir, f"ecology_day_{t_days:05.1f}.png")
                        os.system(f"open '{panel_path}'")
                except Exception:
                    pass
            except Exception as _ev:
                if int(os.getenv("QD_ECO_DIAG", "1")) == 1:
                    print(f"[EcologyViz] skipped: {_ev}")
            # Diagnostics: per-star ISR components (disabled by default; enable with QD_PLOT_ISR=1)
            if PLOT_ISR:
                plot_isr_components(grid, gcm, t_days, output_dir)


    # --- Optional: Save restart at end ---
    restart_out = os.getenv("QD_RESTART_OUT")
    if restart_out:
        try:
            t_sec_final = float(current_day[0]) * day_in_seconds
        except Exception:
            t_sec_final = 0.0
        try:
            save_restart(
                restart_out, grid, gcm, ocean, land_mask,
                W_land=W_land, S_snow=S_snow, C_snow=getattr(gcm, "C_snow_map_last", None), t_seconds=t_sec_final
            )
            print(f"[Restart] Saved final state to '{restart_out}'.")
            # Also export ocean.nc alongside restart_out (standardized ocean state)
            try:
                out_dir = os.path.dirname(restart_out) or "."
                oce_path = os.path.join(out_dir, "ocean.nc")
                if ocean is not None:
                    _ok_o = save_ocean(oce_path, grid, ocean, day_value=current_day[0])
                    if _ok_o:
                        print(f"[Restart] Ocean state saved to '{oce_path}'.")
            except Exception as _roe:
                print(f"[Restart] Ocean state save skipped: {_roe}")
        except Exception as e:
            print(f"[Restart] Failed to save '{restart_out}': {e}")

    print("\n--- Simulation Finished ---")
    print("Final state diagnostics:")
    print(f"  Max absolute zonal wind (u): {np.max(np.abs(gcm.u)):.2f} m/s")
    print(f"  Max absolute meridional wind (v): {np.max(np.abs(gcm.v)):.2f} m/s")
    print(f"  Max absolute height anomaly (h): {np.max(np.abs(gcm.h)):.1f} m")

if __name__ == "__main__":
    main()

# scripts/run_simulation.py

"""
Main simulation script for the Qingdai GCM.
"""

import numpy as np
import sys
import os
import matplotlib.pyplot as plt

# Add the project root to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pygcm.constants as constants

from pygcm.grid import SphericalGrid
from pygcm.orbital import OrbitalSystem
from pygcm.forcing import ThermalForcing
from pygcm.dynamics import SpectralModel
from pygcm.topography import create_land_sea_mask, generate_base_properties, load_topography_from_netcdf
from pygcm.physics import diagnose_precipitation, diagnose_precipitation_hybrid, parameterize_cloud_cover, calculate_dynamic_albedo, cloud_from_precip, compute_orographic_factor
from pygcm.hydrology import get_hydrology_params_from_env, partition_precip_phase, snow_step, update_land_bucket, diagnose_water_closure
from pygcm.routing import RiverRouting
from pygcm import energy as energy
from pygcm.ocean import WindDrivenSlabOcean
from pygcm.jax_compat import is_enabled as JAX_IS_ENABLED

# --- Restart I/O (NetCDF) and Initialization Helpers ---

def save_restart(path, grid, gcm, ocean, land_mask, W_land=None, S_snow=None):
    """
    Save minimal prognostic state to a NetCDF restart file.
    Includes: lat/lon coords, u/v/h, T_s, cloud_cover, q (if exists), h_ice (if exists),
              ocean state (uo/vo/eta/Ts if ocean provided),
              hydrology reservoirs (W_land, S_snow) if provided.
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

        # Hydrology
        wvar("W_land", W_land)
        wvar("S_snow", S_snow)

        # Masks for reference
        wvar("land_mask", land_mask)

        # Minimal metadata
        ds.setncattr("title", "Qingdai GCM Restart")
        ds.setncattr("creator", "PyGCM for Qingdai")
        ds.setncattr("note", "Contains minimal prognostic fields for warm restart.")
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
                     "uo", "vo", "eta", "Ts", "W_land", "S_snow", "land_mask"]:
            out[name] = rvar(name)
    return out

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
                # Binary mask of major rivers, draw as contour lines
                river_mask = (flow >= river_min).astype(float)
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

def plot_true_color(grid, gcm, land_mask, t_days, output_dir, routing=None):
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

    # Optional land snow rendering（默认关闭；如需可由环境变量开启）
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
            river_mask = (flow >= river_min).astype(float)[..., None]
            rgb_map = rgb_map * (1.0 - river_alpha * river_mask) + river_color * (river_alpha * river_mask)
        lake_mask = getattr(routing, "lake_mask", None)
        if lake_mask is not None and np.any(lake_mask):
            lake_color = np.array([0.15, 0.55, 0.95])
            lake_alpha = float(os.getenv("QD_LAKE_ALPHA", "0.40"))
            lake_mask3 = lake_mask.astype(float)[..., None]
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


def main():
    """
    Main function to run the simulation.
    """
    print("--- Initializing Qingdai GCM ---")
    try:
        print(f"[JAX] Acceleration enabled: {JAX_IS_ENABLED()} (toggle via QD_USE_JAX=1; platform via QD_JAX_PLATFORM=cpu|gpu|tpu)")
    except Exception:
        pass

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
        hydro_net = os.getenv("QD_HYDRO_NETCDF", "data/hydrology_network.nc")
        if int(os.getenv("QD_HYDRO_ENABLE", "1")) == 1 and hydro_net and os.path.exists(hydro_net):
            routing = RiverRouting(
                grid,
                hydro_net,
                dt_hydro_hours=float(os.getenv("QD_HYDRO_DT_HOURS", "6")),
                treat_lake_as_water=(int(os.getenv("QD_TREAT_LAKE_AS_WATER", "1")) == 1),
                alpha_lake=(float(os.getenv("QD_ALPHA_LAKE")) if os.getenv("QD_ALPHA_LAKE") else None),
                diag=(int(os.getenv("QD_HYDRO_DIAG", "1")) == 1),
            )
        else:
            print(f"[HydroRouting] Disabled or network file not found (QD_HYDRO_NETCDF='{hydro_net}').")
    except Exception as e:
        print(f"[HydroRouting] Initialization skipped due to error: {e}")
        routing = None

    # --- Restart load or banded initialization ---
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
            # Hydrology
            if rst.get("W_land") is not None: W_land = rst["W_land"]
            if rst.get("S_snow") is not None: S_snow = rst["S_snow"]
            print(f"[Restart] Loaded state from '{restart_in}'.")
        except Exception as e:
            print(f"[Restart] Failed to load '{restart_in}': {e}\nContinuing with fresh init.")
            apply_banded_initial_ts(grid, gcm, ocean, land_mask)
    else:
        apply_banded_initial_ts(grid, gcm, ocean, land_mask)

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
    
    time_steps = np.arange(0, sim_duration_seconds, dt)

    # --- Visualization Parameters ---
    output_dir = "output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    plot_every_days = float(os.getenv("QD_PLOT_EVERY_DAYS", "10"))
    plot_interval_seconds = plot_every_days * 24 * 3600
    plot_interval_steps = max(1, int(plot_interval_seconds / dt))
    print(f"Generating a plot every {plot_interval_steps} steps.")
    
    print(f"\n--- Starting Simulation ---")
    print(f"Grid resolution: {grid.n_lat} lat x {grid.n_lon} lon")
    print(f"Time step (dt): {dt} s")
    print(f"Simulation duration: {sim_duration_seconds / day_in_seconds:.1f} planetary days")
    print(f"Total time steps: {len(time_steps)}")

    # --- Precipitation accumulation over one planetary day (kg m^-2 ≡ mm/day) ---
    precip_acc_day = np.zeros_like(grid.lat_mesh, dtype=float)
    accum_t_day = 0.0
    precip_day_last = None

    # 2. Time Integration Loop
    try:
        from tqdm import tqdm
        iterator = tqdm(time_steps)
    except ImportError:
        print("tqdm not found, using simple print statements for progress.")
        iterator = time_steps

    for i, t in enumerate(iterator):
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

        # 2) Radiative albedo using updated cloud cover
        #    Use sea-ice fraction derived from thickness (saturating), consistent with M2 thermodynamics.
        H_ice_ref = float(os.getenv("QD_HICE_REF", "0.5"))  # m, e-folding thickness for ice optical effect
        ice_frac = 1.0 - np.exp(-np.maximum(gcm.h_ice, 0.0) / max(1e-6, H_ice_ref))
        # M4: Use humidity/precipitation-coupled cloud for radiation/albedo if available
        cloud_for_rad = getattr(gcm, "cloud_eff_last", gcm.cloud_cover)

        if USE_TOPO_ALBEDO:
            albedo = calculate_dynamic_albedo(
                cloud_for_rad, gcm.T_s, base_albedo_map, alpha_ice, alpha_cloud, land_mask=land_mask, ice_frac=ice_frac
            )
        else:
            albedo = calculate_dynamic_albedo(
                cloud_for_rad, gcm.T_s, alpha_water, alpha_ice, alpha_cloud, land_mask=land_mask, ice_frac=ice_frac
            )

        # 2. Forcing Step
        # Compute two-star insolation components and diagnostics
        insA, insB = forcing.calculate_insolation_components(t)
        gcm.isr_A, gcm.isr_B = insA, insB
        gcm.isr = insA + insB
        Teq = forcing.calculate_equilibrium_temp(t, albedo)

        # 3. Dynamics Step
        gcm.time_step(Teq, dt, albedo=albedo)

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
            # Fluxes from humidity module (kg m^-2 s^-1)
            E_flux = getattr(gcm, "E_flux_last", 0.0)
            P_flux = getattr(gcm, "P_cond_flux_last", 0.0)

            # Ensure array shape
            if np.isscalar(E_flux):
                E_flux = np.full_like(gcm.T_s, float(E_flux))
            if np.isscalar(P_flux):
                P_flux = np.full_like(gcm.T_s, float(P_flux))

            # Phase partition (rain/snow) by surface temperature
            P_rain, P_snow = partition_precip_phase(P_flux, gcm.T_s, T_thresh=hydro_params.snow_thresh_K)

            land = (land_mask == 1)
            # Land components
            P_rain_land = P_rain * land
            P_snow_land = P_snow * land
            E_land = E_flux * land

            # Update land snow reservoir and compute melt flux (kg m^-2 s^-1)
            S_snow, melt_flux_land = snow_step(S_snow, P_snow_land, gcm.T_s, hydro_params, dt)

            # Land bucket update: inputs = rain + snowmelt; evaporation removes; runoff returned to ocean
            P_in_land = P_rain_land + melt_flux_land
            W_land, R_flux_land = update_land_bucket(W_land, P_in_land, E_land, hydro_params, dt)
            # P014: route runoff along offline network (if enabled)
            if 'routing' in locals() and routing is not None:
                try:
                    routing.step(R_land_flux=R_flux_land, dt_seconds=dt, precip_flux=P_flux, evap_flux=E_flux)
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
                    P_flux=P_flux,
                    R_flux=R_flux_land,   # runoff only from land
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
            plot_true_color(grid, gcm, land_mask, t_days, output_dir, routing=routing)
            # Diagnostics: per-star ISR components (disabled by default; enable with QD_PLOT_ISR=1)
            if PLOT_ISR:
                plot_isr_components(grid, gcm, t_days, output_dir)


    # --- Optional: Save restart at end ---
    restart_out = os.getenv("QD_RESTART_OUT")
    if restart_out:
        try:
            save_restart(restart_out, grid, gcm, ocean, land_mask, W_land=W_land, S_snow=S_snow)
            print(f"[Restart] Saved final state to '{restart_out}'.")
        except Exception as e:
            print(f"[Restart] Failed to save '{restart_out}': {e}")

    print("\n--- Simulation Finished ---")
    print("Final state diagnostics:")
    print(f"  Max absolute zonal wind (u): {np.max(np.abs(gcm.u)):.2f} m/s")
    print(f"  Max absolute meridional wind (v): {np.max(np.abs(gcm.v)):.2f} m/s")
    print(f"  Max absolute height anomaly (h): {np.max(np.abs(gcm.h)):.1f} m")

if __name__ == "__main__":
    main()

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
from pygcm.physics import diagnose_precipitation, parameterize_cloud_cover, calculate_dynamic_albedo, cloud_from_precip, compute_orographic_factor
from pygcm.hydrology import get_hydrology_params_from_env, partition_precip_phase, snow_step, update_land_bucket, diagnose_water_closure

def plot_state(grid, gcm, land_mask, precip, cloud_cover, albedo, t_days, output_dir):
    """
    Generates and saves a plot of the current model state.
    """
    fig, axes = plt.subplots(4, 2, figsize=(18, 20), constrained_layout=True)
    fig.suptitle(f"Qingdai GCM State at Day {t_days:.2f}", fontsize=16)
    ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8 = axes.flatten()

    # 1. Plot Surface Temperature and Land Mask
    ts_plot = ax1.contourf(grid.lon, grid.lat, gcm.T_s - 273.15, levels=20, cmap='coolwarm')
    ax1.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors='black', linewidths=1)
    ax1.set_title("Surface Temperature (°C)")
    fig.colorbar(ts_plot, ax=ax1, label="°C")

    # 2. Plot Geopotential Height Anomaly
    h_plot = ax2.contourf(grid.lon, grid.lat, gcm.h, levels=20, cmap='RdBu_r')
    ax2.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors='black', linewidths=1)
    ax2.set_title("Geopotential Height Anomaly (m)")
    fig.colorbar(h_plot, ax=ax2, label="m")

    # 3. Plot Wind Field (Streamline plot)
    speed = np.sqrt(gcm.u**2 + gcm.v**2)
    ax3.streamplot(grid.lon, grid.lat, gcm.u, gcm.v, color=speed, cmap='viridis', density=1.5)
    ax3.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors='black', linewidths=1)
    ax3.set_title("Wind Field (m/s)")
    
    # 4. Plot Precipitation
    precip_plot = ax4.contourf(grid.lon, grid.lat, precip * 1e5, levels=np.linspace(0.1, 5, 10), cmap='Blues', extend='max')
    ax4.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors='black', linewidths=1)
    ax4.set_title("Precipitation Rate (arbitrary units)")
    fig.colorbar(precip_plot, ax=ax4, label="Rate")

    # 5. Plot Cloud Cover
    cloud_plot = ax5.contourf(grid.lon, grid.lat, cloud_cover, levels=np.linspace(0, 1, 11), cmap='Greys')
    ax5.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors='black', linewidths=1)
    ax5.set_title("Cloud Cover Fraction")
    fig.colorbar(cloud_plot, ax=ax5, label="Fraction")

    # 6. Plot Albedo
    albedo_plot = ax6.contourf(grid.lon, grid.lat, albedo, levels=np.linspace(0, 0.8, 11), cmap='cividis')
    ax6.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors='black', linewidths=1)
    ax6.set_title("Dynamic Albedo")
    fig.colorbar(albedo_plot, ax=ax6, label="Albedo")

    # 7. Incoming Shortwave Radiation (Total) with star centers overlaid
    isr_plot = ax7.contourf(grid.lon, grid.lat, gcm.isr, levels=20, cmap='magma')
    ax7.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors='white', linewidths=0.5)

    # Overlay star A/B subsolar points (maxima of each component)
    try:
        idxA = np.unravel_index(np.argmax(gcm.isr_A), gcm.isr_A.shape)
        idxB = np.unravel_index(np.argmax(gcm.isr_B), gcm.isr_B.shape)
        lonA, latA = grid.lon[idxA[1]], grid.lat[idxA[0]]
        lonB, latB = grid.lon[idxB[1]], grid.lat[idxB[0]]
        ax7.scatter([lonA], [latA], c='cyan', s=30, marker='x', label='Star A center')
        ax7.scatter([lonB], [latB], c='yellow', s=30, marker='+', label='Star B center')
        ax7.legend(loc='upper right', fontsize=8)
    except Exception:
        pass

    ax7.set_title("Incoming Shortwave (W/m^2, total)")
    fig.colorbar(isr_plot, ax=ax7, label="W/m^2")

    # 8. Outgoing Longwave Radiation
    olr_plot = ax8.contourf(grid.lon, grid.lat, gcm.olr, levels=20, cmap='plasma')
    ax8.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors='white', linewidths=0.5)
    ax8.set_title("Outgoing Longwave (W/m^2)")
    fig.colorbar(olr_plot, ax=ax8, label="W/m^2")

    for ax in axes.flatten():
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_xlim(0, 360)
        ax.set_ylim(-90, 90)

    # Save the figure
    filename = os.path.join(output_dir, f"state_day_{t_days:05.1f}.png")
    plt.savefig(filename)
    plt.close(fig)

def plot_true_color(grid, gcm, land_mask, t_days, output_dir):
    """
    Generates and saves a pseudo-true-color plot of the planet.
    """
    # Define colors
    ocean_color = np.array([0.1, 0.2, 0.5])
    land_color = np.array([0.4, 0.3, 0.2])
    ice_color = np.array([0.9, 0.9, 0.95])
    
    # Create base map (land/ocean/ice)
    rgb_map = np.zeros((grid.n_lat, grid.n_lon, 3))
    
    # Ocean
    rgb_map[land_mask == 0] = ocean_color
    # Land
    rgb_map[land_mask == 1] = land_color
    # Ice
    rgb_map[gcm.T_s < 273.15] = ice_color
    
    # Add clouds (as a white layer)
    # The cloud layer is semi-transparent, so we blend it
    cloud_layer = np.stack([gcm.cloud_cover, gcm.cloud_cover, gcm.cloud_cover], axis=-1)
    rgb_map = rgb_map * (1 - cloud_layer) + cloud_layer * 0.9 # 0.9 makes clouds slightly off-white
    
    # Clamp values to be safe
    rgb_map = np.clip(rgb_map, 0, 1)
    
    # Plotting
    fig, ax = plt.subplots(1, 1, figsize=(12, 6), constrained_layout=True)
    ax.imshow(rgb_map, extent=[0, 360, -90, 90], origin='lower')
    ax.set_title(f"Qingdai 'True Color' at Day {t_days:.2f}")
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    
    # Save the figure
    filename = os.path.join(output_dir, f"true_color_day_{t_days:05.1f}.png")
    plt.savefig(filename)
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

    print("Initializing dynamics core with surface friction and greenhouse effect...")
    gcm = SpectralModel(
        grid, friction_map, H=8000, tau_rad=10 * 24 * 3600, greenhouse_factor=0.3,
        C_s_map=C_s_map, land_mask=land_mask, Cs_ocean=Cs_ocean, Cs_land=Cs_land, Cs_ice=Cs_ice
    )

    # --- Hydrology (P009): reservoirs and parameters ---
    hydro_params = get_hydrology_params_from_env()
    W_land = np.zeros_like(grid.lat_mesh, dtype=float)   # land bucket water (kg m^-2)
    S_snow = np.zeros_like(grid.lat_mesh, dtype=float)   # land snow water-equivalent (kg m^-2)
    _hydro_prev_total = None
    _hydro_prev_time = None

    # --- Simulation Parameters ---
    dt = int(os.getenv("QD_DT_SECONDS", "300"))  # default 300s
    day_in_seconds = 2 * np.pi / constants.PLANET_OMEGA
    if os.getenv("QD_SIM_DAYS"):
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

    # 2. Time Integration Loop
    try:
        from tqdm import tqdm
        iterator = tqdm(time_steps)
    except ImportError:
        print("tqdm not found, using simple print statements for progress.")
        iterator = time_steps

    for i, t in enumerate(iterator):
        # 1. Physics Step
        # 1) Diagnose precipitation from dynamics WITHOUT cloud gating (avoid circularity)
        precip = diagnose_precipitation(gcm, grid, D_crit, k_precip, cloud_threshold=None)
        # Optional orographic enhancement (uses upslope wind and elevation)
        if OROG_ENABLED and (elevation is not None):
            try:
                factor = compute_orographic_factor(grid, elevation, gcm.u, gcm.v, k_orog=K_OROG)
                precip = precip * factor
            except Exception as e:
                if i == 0:
                    print(f"[Orog] Disabled due to error: {e}")
        
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
            plot_state(grid, gcm, land_mask, precip, gcm.cloud_cover, albedo, t_days, output_dir)
            plot_true_color(grid, gcm, land_mask, t_days, output_dir)
            # Diagnostics: per-star ISR components (disabled by default; enable with QD_PLOT_ISR=1)
            if PLOT_ISR:
                plot_isr_components(grid, gcm, t_days, output_dir)


    print("\n--- Simulation Finished ---")
    print("Final state diagnostics:")
    print(f"  Max absolute zonal wind (u): {np.max(np.abs(gcm.u)):.2f} m/s")
    print(f"  Max absolute meridional wind (v): {np.max(np.abs(gcm.v)):.2f} m/s")
    print(f"  Max absolute height anomaly (h): {np.max(np.abs(gcm.h)):.1f} m")

if __name__ == "__main__":
    main()

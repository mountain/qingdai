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
from pygcm.topography import create_land_sea_mask, generate_base_properties
from pygcm.physics import diagnose_precipitation, parameterize_cloud_cover, calculate_dynamic_albedo, cloud_from_precip

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
    land_mask = create_land_sea_mask(grid)

    print("Initializing orbital mechanics...")
    orbital_sys = OrbitalSystem()

    print("Initializing thermal forcing...")
    forcing = ThermalForcing(grid, orbital_sys)

    print("Initializing dynamics core with surface friction and greenhouse effect...")
    base_albedo_map, friction_map = generate_base_properties(land_mask)
    gcm = SpectralModel(grid, friction_map, H=8000, tau_rad=10 * 24 * 3600, greenhouse_factor=0.3)

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
    alpha_water = 0.1   # Albedo of water/land
    alpha_ice = 0.6     # Albedo of ice
    alpha_cloud = 0.5   # Average cloud albedo
    
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
        albedo = calculate_dynamic_albedo(gcm.cloud_cover, gcm.T_s, alpha_water, alpha_ice, alpha_cloud)

        # 2. Forcing Step
        # Compute two-star insolation components and diagnostics
        insA, insB = forcing.calculate_insolation_components(t)
        gcm.isr_A, gcm.isr_B = insA, insB
        gcm.isr = insA + insB
        Teq = forcing.calculate_equilibrium_temp(t, albedo)

        # 3. Dynamics Step
        gcm.time_step(Teq, dt)

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
            # Diagnostics: save per-star ISR components to verify double centers
            plot_isr_components(grid, gcm, t_days, output_dir)


    print("\n--- Simulation Finished ---")
    print("Final state diagnostics:")
    print(f"  Max absolute zonal wind (u): {np.max(np.abs(gcm.u)):.2f} m/s")
    print(f"  Max absolute meridional wind (v): {np.max(np.abs(gcm.v)):.2f} m/s")
    print(f"  Max absolute height anomaly (h): {np.max(np.abs(gcm.h)):.1f} m")

if __name__ == "__main__":
    main()

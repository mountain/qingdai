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

def plot_state(grid, gcm, t_days, output_dir):
    """
    Generates and saves a plot of the current model state.
    """
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10), constrained_layout=True)
    fig.suptitle(f"Qingdai GCM State at Day {t_days:.2f}", fontsize=16)

    # 1. Plot Geopotential Height Anomaly
    h_plot = ax1.contourf(grid.lon, grid.lat, gcm.h, levels=20, cmap='RdBu_r')
    ax1.set_title("Geopotential Height Anomaly (m)")
    ax1.set_xlabel("Longitude")
    ax1.set_ylabel("Latitude")
    fig.colorbar(h_plot, ax=ax1, label="Height (m)")

    # 2. Plot Wind Field (Streamline plot)
    speed = np.sqrt(gcm.u**2 + gcm.v**2)
    ax2.streamplot(grid.lon, grid.lat, gcm.u, gcm.v, color=speed, cmap='viridis', density=2)
    ax2.set_title("Wind Field Streamlines (m/s)")
    ax2.set_xlabel("Longitude")
    ax2.set_ylabel("Latitude")
    ax2.set_xlim(0, 360)
    ax2.set_ylim(-90, 90)

    # Save the figure
    filename = os.path.join(output_dir, f"state_day_{t_days:05.1f}.png")
    plt.savefig(filename)
    plt.close(fig)

def main():
    """
    Main function to run the simulation.
    """
    print("--- Initializing Qingdai GCM ---")

    # 1. Initialization
    print("Creating grid...")
    grid = SphericalGrid(n_lat=121, n_lon=240)

    print("Initializing orbital mechanics...")
    orbital_sys = OrbitalSystem()

    print("Initializing thermal forcing...")
    forcing = ThermalForcing(grid, orbital_sys)

    print("Initializing dynamics core (Spectral Model)...")
    gcm = SpectralModel(grid, H=8000, tau_rad=10 * 24 * 3600) # 10 day relaxation

    # --- Simulation Parameters ---
    dt = 300  # Time step in seconds (5 minutes)
    
    # Simulate for one full planetary year
    day_in_seconds = 2 * np.pi / constants.PLANET_OMEGA
    sim_duration_seconds = orbital_sys.T_planet
    
    time_steps = np.arange(0, sim_duration_seconds, dt)

    # --- Visualization Parameters ---
    output_dir = "output"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # A "week" is arbitrary, let's define it as ~10 Earth days
    plot_interval_seconds = 10 * 24 * 3600
    plot_interval_steps = int(plot_interval_seconds / dt)
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
        # 1. Calculate thermal forcing for the current time
        Teq = forcing.calculate_equilibrium_temp(t)

        # 2. Advance the GCM state by one time step
        gcm.time_step(Teq, dt)

        # 3. (Optional) Print diagnostics and generate plots
        t_days = t / day_in_seconds
        if i % 100 == 0 and i > 0:
            if not isinstance(iterator, type(time_steps)): # if tqdm is used
                iterator.set_description(
                    f"t={t_days:.1f}d | "
                    f"max|u|={np.max(np.abs(gcm.u)):.2f} m/s | "
                    f"max|v|={np.max(np.abs(gcm.v)):.2f} m/s"
                )
        
        # Generate plot at specified interval
        if i % plot_interval_steps == 0:
            plot_state(grid, gcm, t_days, output_dir)


    print("\n--- Simulation Finished ---")
    print("Final state diagnostics:")
    print(f"  Max absolute zonal wind (u): {np.max(np.abs(gcm.u)):.2f} m/s")
    print(f"  Max absolute meridional wind (v): {np.max(np.abs(gcm.v)):.2f} m/s")
    print(f"  Max absolute height anomaly (h): {np.max(np.abs(gcm.h)):.1f} m")

if __name__ == "__main__":
    main()

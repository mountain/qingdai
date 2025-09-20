#!/usr/bin/env python3
"""
Diagnostic: visualize per-star incoming shortwave (ISR) to verify the expected
'double centers' (subsolar points) from the two stars.

This runs quickly without the full GCM time integration.
"""

import os
import math
import numpy as np
import matplotlib.pyplot as plt
import sys

# Add the project root to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pygcm.grid import SphericalGrid
from pygcm.orbital import OrbitalSystem
from pygcm.forcing import ThermalForcing
import pygcm.constants as const


def plot_isr_components(grid, isr_A, isr_B, t_days, out_path):
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(16, 6), constrained_layout=True)

    vmin = 0.0
    vmax = max(np.max(isr_A), np.max(isr_B))
    levels = np.linspace(vmin, vmax, 21)

    csA = axA.contourf(grid.lon, grid.lat, isr_A, levels=levels, cmap='magma')
    axA.set_title(f"ISR - Star A (Day {t_days:.2f})")
    axA.set_xlabel("Longitude")
    axA.set_ylabel("Latitude")
    axA.set_xlim(0, 360)
    axA.set_ylim(-90, 90)
    fig.colorbar(csA, ax=axA, label="W/m^2")

    csB = axB.contourf(grid.lon, grid.lat, isr_B, levels=levels, cmap='magma')
    axB.set_title(f"ISR - Star B (Day {t_days:.2f})")
    axB.set_xlabel("Longitude")
    axB.set_ylabel("Latitude")
    axB.set_xlim(0, 360)
    axB.set_ylim(-90, 90)
    fig.colorbar(csB, ax=axB, label="W/m^2")

    # Compute and print subsolar points and great-circle separation
    try:
        idxA = np.unravel_index(np.argmax(isr_A), isr_A.shape)
        idxB = np.unravel_index(np.argmax(isr_B), isr_B.shape)
        lonA, latA = grid.lon[idxA[1]], grid.lat[idxA[0]]
        lonB, latB = grid.lon[idxB[1]], grid.lat[idxB[0]]
        axA.scatter([lonA], [latA], c='cyan', s=40, marker='x', label='A center')
        axB.scatter([lonB], [latB], c='yellow', s=40, marker='+', label='B center')
        axA.legend(loc='upper right', fontsize=8)
        axB.legend(loc='upper right', fontsize=8)

        phi1 = math.radians(latA); lam1 = math.radians(lonA)
        phi2 = math.radians(latB); lam2 = math.radians(lonB)
        dlam = lam2 - lam1
        d_sigma = 2 * math.asin(math.sqrt(
            math.sin((phi2 - phi1)/2)**2 +
            math.cos(phi1)*math.cos(phi2)*math.sin(dlam/2)**2
        ))
        separation_deg = math.degrees(d_sigma)
        print(f"[Diag] Day {t_days:.2f}: Subsolar separation ≈ {separation_deg:.2f}° "
              f"(A: lon={lonA:.1f}°, lat={latA:.1f}°; B: lon={lonB:.1f}°, lat={latB:.1f}°)")
    except Exception as e:
        print(f"[Diag] Could not compute subsolar separation: {e}")

    plt.savefig(out_path)
    plt.close(fig)


def main():
    grid = SphericalGrid(n_lat=121, n_lon=240)
    orbital_sys = OrbitalSystem()
    forcing = ThermalForcing(grid, orbital_sys)

    output_dir = "output"
    os.makedirs(output_dir, exist_ok=True)

    # Planetary day in seconds for labeling
    day_in_seconds = 2 * np.pi / const.PLANET_OMEGA

    # Choose a few times: t=0, quarter and half of the binary period
    times = [
        0.0,
        0.25 * orbital_sys.T_binary,
        0.50 * orbital_sys.T_binary,
        0.75 * orbital_sys.T_binary,
    ]

    for t in times:
        isr_A, isr_B = forcing.calculate_insolation_components(t)
        t_days = t / day_in_seconds
        out_path = os.path.join(output_dir, f"diag_isr_components_day_{t_days:05.1f}.png")
        plot_isr_components(grid, isr_A, isr_B, t_days, out_path)
        # Also save the total field for reference
        fig, ax = plt.subplots(1, 1, figsize=(10, 5), constrained_layout=True)
        total = isr_A + isr_B
        cs = ax.contourf(grid.lon, grid.lat, total, levels=20, cmap='magma')
        ax.set_title(f"ISR - Total (Day {t_days:.2f})")
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_xlim(0, 360)
        ax.set_ylim(-90, 90)
        fig.colorbar(cs, ax=ax, label="W/m^2")
        out_total = os.path.join(output_dir, f"diag_isr_total_day_{t_days:05.1f}.png")
        plt.savefig(out_total)
        plt.close(fig)

    print("Diagnostics complete. See PNGs under ./output/.")


if __name__ == "__main__":
    main()

# scripts/test_orbital_module.py

"""
Test script to verify the functionality of the pygcm.orbital module.
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Add the project root to the Python path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from pygcm.orbital import OrbitalSystem

def main():
    """
    Main function to run the test.
    """
    print("Testing pygcm.orbital module...")

    # 1. Initialize the orbital system
    orbital_system = OrbitalSystem()

    # 2. Set up a time array for one planetary year
    sim_duration_days = orbital_system.T_planet / (3600 * 24)
    t_days = np.linspace(0, sim_duration_days, 1000)
    t_seconds = t_days * 3600 * 24

    # 3. Calculate the flux at each time step
    flux = orbital_system.calculate_total_flux(t_seconds)

    # 4. Print some key results
    print(f"Binary Period (Pulse Season): {orbital_system.T_binary / (3600 * 24):.2f} Earth days")
    print(f"Planet Period (Year): {orbital_system.T_planet / (3600 * 24):.2f} Earth days")
    print(f"Mean Flux: {np.mean(flux):.2f} W/m^2")
    print(f"Max Flux: {np.max(flux):.2f} W/m^2")
    print(f"Min Flux: {np.min(flux):.2f} W/m^2")

    # 5. Plot the results for visual verification
    plt.figure(figsize=(12, 6))
    plt.plot(t_days, flux)
    plt.title("Energy Flux on Qingdai (Test Run)")
    plt.xlabel("Time (Earth Days)")
    plt.ylabel("Flux (W/m^2)")
    plt.grid(True)
    plt.show()

    print("\nTest complete. Check the plot to verify results.")

if __name__ == "__main__":
    main()

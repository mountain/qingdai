# pygcm/constants.py

"""
Central repository for all physical and astronomical constants and parameters
for the Qingdai simulation.
"""

# --- Physical Constants (SI units) ---
G = 6.67430e-11      # Gravitational constant (m^3 kg^-1 s^-2)
SIGMA = 5.670374e-8  # Stefan-Boltzmann constant (W m^-2 K^-4)

# --- Astronomical Units ---
M_SUN = 1.989e30     # Mass of the Sun (kg)
L_SUN = 3.828e26     # Luminosity of the Sun (W)
AU = 1.496e11        # Astronomical Unit (m)

# --- Harmony Star System Parameters (from docs/01-astronomical-setting.md) ---
# Star A (G6V)
M_A = 0.914 * M_SUN
L_A = 0.7 * L_SUN

# Star B (K1V)
M_B = 0.8 * M_SUN
L_B = 0.410 * L_SUN

# Binary System
M_TOTAL_STARS = M_A + M_B
A_BINARY = 0.5 * AU  # Semi-major axis of the binary pair

# --- Qingdai Planet Parameters ---
A_PLANET = 1.32 * AU     # Semi-major axis of the planet's orbit (final, final tuning)
PLANET_RADIUS = 6.371e6  # Planet radius (m, Earth-like for now)
PLANET_ALBEDO = 0.3      # Planetary albedo (Earth-like for now)
PLANET_OMEGA = 8.726646259971648e-5  # Angular velocity (rad/s, for a 20-hour day)
PLANET_AXIAL_TILT = 27.0 # Axial tilt in degrees

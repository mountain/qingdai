# pygcm/orbital.py

"""
Calculates the orbital mechanics of the Harmony-Qingdai system.
"""

import numpy as np
from . import constants as const

class OrbitalSystem:
    """
    Handles the calculation of stellar positions and total energy flux
    received by the planet Qingdai.
    """
    def __init__(self):
        """
        Initializes the orbital system with parameters from the constants module.
        """
        # Binary star orbital period
        self.T_binary = 2 * np.pi * np.sqrt(const.A_BINARY**3 / (const.G * const.M_TOTAL_STARS))
        
        # Planet orbital period
        self.T_planet = 2 * np.pi * np.sqrt(const.A_PLANET**3 / (const.G * const.M_TOTAL_STARS))

        # Angular velocities
        self.omega_binary = 2 * np.pi / self.T_binary
        self.omega_planet = 2 * np.pi / self.T_planet

        # Individual stellar orbital radii around the barycenter
        self.r_A = const.A_BINARY * (const.M_B / const.M_TOTAL_STARS)
        self.r_B = const.A_BINARY * (const.M_A / const.M_TOTAL_STARS)

    def calculate_stellar_positions(self, t):
        """
        Calculates the Cartesian coordinates of Star A and Star B at a given time t.

        Args:
            t (float or np.ndarray): Time in seconds.

        Returns:
            tuple: A tuple containing the (x, y) coordinates of Star A and Star B.
                   (x_A, y_A, x_B, y_B)
        """
        x_A = self.r_A * np.cos(self.omega_binary * t)
        y_A = self.r_A * np.sin(self.omega_binary * t)
        x_B = -self.r_B * np.cos(self.omega_binary * t)
        y_B = -self.r_B * np.sin(self.omega_binary * t)
        return x_A, y_A, x_B, y_B

    def calculate_total_flux(self, t):
        """
        Calculates the total stellar energy flux (W/m^2) received by Qingdai
        at a given time t.

        Args:
            t (float or np.ndarray): Time in seconds.

        Returns:
            float or np.ndarray: The total energy flux.
        """
        # Get stellar positions
        x_A, y_A, x_B, y_B = self.calculate_stellar_positions(t)

        # Calculate planet's position
        x_p = const.A_PLANET * np.cos(self.omega_planet * t)
        y_p = const.A_PLANET * np.sin(self.omega_planet * t)

        # Calculate distances from the planet to each star
        d_A = np.sqrt((x_p - x_A)**2 + (y_p - y_A)**2)
        d_B = np.sqrt((x_p - x_B)**2 + (y_p - y_B)**2)

        # Calculate flux from each star
        S_A = const.L_A / (4 * np.pi * d_A**2)
        S_B = const.L_B / (4 * np.pi * d_B**2)

        # Total flux is the sum of the two
        return S_A + S_B

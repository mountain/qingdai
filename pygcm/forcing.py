# pygcm/forcing.py

"""
Calculates the thermal forcing for the GCM.
"""

import numpy as np
from . import constants as const
from .grid import SphericalGrid
from .orbital import OrbitalSystem

class ThermalForcing:
    """
    Calculates the dynamic equilibrium temperature field for the planet.
    """
    def __init__(self, grid: SphericalGrid, orbital_system: OrbitalSystem):
        """
        Initializes the thermal forcing module.

        Args:
            grid (SphericalGrid): The model's grid system.
            orbital_system (OrbitalSystem): The orbital system providing energy flux.
        """
        self.grid = grid
        self.orbital_system = orbital_system
        self.planet_params = {
            'axial_tilt': const.PLANET_AXIAL_TILT,
            'omega': const.PLANET_OMEGA,
            'T_planet': self.orbital_system.T_planet
        }
        # Precompute planet-fixed equatorial frame (in inertial coords)
        # Rotation axis (tilted from orbital normal toward +x by axial tilt)
        tilt_rad = np.deg2rad(self.planet_params['axial_tilt'])
        self.n_hat = np.array([np.sin(tilt_rad), 0.0, np.cos(tilt_rad)])  # unit rotation axis
        # x_eq: projection of inertial +x onto equatorial plane
        x_inertial = np.array([1.0, 0.0, 0.0])
        self.x_eq = x_inertial - np.dot(x_inertial, self.n_hat) * self.n_hat
        self.x_eq /= np.linalg.norm(self.x_eq)
        # y_eq completes right-handed system
        self.y_eq = np.cross(self.n_hat, self.x_eq)

    def calculate_insolation(self, t):
        """
        Calculates the instantaneous solar radiation (insolation) at each
        grid point, accounting for the two separate stars.

        Args:
            t (float): Time in seconds.

        Returns:
            np.ndarray: A 2D array of total insolation values (W/m^2).
        """
        # Planet's orbital angle
        planet_orbital_angle = self.orbital_system.omega_planet * t
        
        # --- Calculate for Star A ---
        x_A, y_A, _, _ = self.orbital_system.calculate_stellar_positions(t)
        x_p = const.A_PLANET * np.cos(planet_orbital_angle)
        y_p = const.A_PLANET * np.sin(planet_orbital_angle)
        
        # Vector from planet to star A
        vec_p_A = np.array([x_A - x_p, y_A - y_p, 0])
        dist_A = np.linalg.norm(vec_p_A)
        flux_A = const.L_A / (4 * np.pi * dist_A**2)
        
        # --- Calculate for Star B ---
        _, _, x_B, y_B = self.orbital_system.calculate_stellar_positions(t)
        vec_p_B = np.array([x_B - x_p, y_B - y_p, 0])
        dist_B = np.linalg.norm(vec_p_B)
        flux_B = const.L_B / (4 * np.pi * dist_B**2)

        # --- Calculate Insolation from Each Star ---
        insolation_A = self._get_single_star_insolation(t, flux_A, vec_p_A)
        insolation_B = self._get_single_star_insolation(t, flux_B, vec_p_B)

        return insolation_A + insolation_B

    def calculate_insolation_components(self, t):
        """
        Same geometry as calculate_insolation, but returns the two-star
        contributions separately for diagnostics.
        Returns:
            (insolation_A, insolation_B): Tuple of 2D arrays (W/m^2)
        """
        planet_orbital_angle = self.orbital_system.omega_planet * t

        x_A, y_A, x_B, y_B = self.orbital_system.calculate_stellar_positions(t)
        x_p = const.A_PLANET * np.cos(planet_orbital_angle)
        y_p = const.A_PLANET * np.sin(planet_orbital_angle)

        vec_p_A = np.array([x_A - x_p, y_A - y_p, 0.0])
        vec_p_B = np.array([x_B - x_p, y_B - y_p, 0.0])

        dist_A = np.linalg.norm(vec_p_A)
        dist_B = np.linalg.norm(vec_p_B)

        flux_A = const.L_A / (4 * np.pi * (dist_A**2))
        flux_B = const.L_B / (4 * np.pi * (dist_B**2))

        insolation_A = self._get_single_star_insolation(t, flux_A, vec_p_A)
        insolation_B = self._get_single_star_insolation(t, flux_B, vec_p_B)

        return insolation_A, insolation_B

    def _get_single_star_insolation(self, t, flux, star_vector):
        """
        Calculate insolation from one star using proper per-star geometry:
        - Declination delta_s = asin( s_hat · n_hat )
        - Right ascension alpha_s from projection on (x_eq, y_eq)
        - Hour angle h = theta + lon - alpha_s, with theta = omega * t
        """
        # Unit vector from planet to star (in inertial frame)
        s_hat = star_vector / (np.linalg.norm(star_vector) + 1e-15)

        # Planet equatorial frame (precomputed in __init__)
        n_hat = self.n_hat
        x_eq = self.x_eq
        y_eq = self.y_eq

        # Per-star declination and right ascension
        delta = np.arcsin(np.clip(np.dot(s_hat, n_hat), -1.0, 1.0))
        alpha = np.arctan2(np.dot(s_hat, y_eq), np.dot(s_hat, x_eq))

        # Local hour angle field
        theta = (t * self.planet_params['omega']) % (2 * np.pi)
        lon_rad = np.deg2rad(self.grid.lon_mesh)
        h = theta + lon_rad - alpha

        # Cosine of solar zenith angle
        lat_rad = np.deg2rad(self.grid.lat_mesh)
        cos_z = np.sin(lat_rad) * np.sin(delta) + np.cos(lat_rad) * np.cos(delta) * np.cos(h)

        # Night side clamp
        cos_z = np.maximum(0.0, cos_z)

        return flux * cos_z

    def calculate_equilibrium_temp(self, t, albedo):
        """
        Calculates the global radiative equilibrium temperature field T_eq
        for a given time t, considering a given albedo field.

        Args:
            t (float): Time in seconds.
            albedo (np.ndarray): Current albedo field.

        Returns:
            np.ndarray: A 2D array of equilibrium temperatures (K) on the grid.
        """
        insolation = self.calculate_insolation(t)
        
        # Stefan-Boltzmann Law: I * (1 - albedo) = sigma * T^4
        # T = (I * (1 - albedo) / sigma)^(1/4)
        
        numerator = insolation * (1 - albedo)
        # Avoid division by zero or negative roots for the night side
        # Where insolation is zero, temperature should be zero (or a minimum background temp)
        # For now, we'll let it be zero.
        
        # To prevent taking the root of a negative number if numerator is somehow negative
        numerator[numerator < 0] = 0
        
        temp_eq_field = (numerator / const.SIGMA)**0.25
        
        return temp_eq_field

if __name__ == '__main__':
    # Example usage
    from .physics import calculate_dynamic_albedo

    grid = SphericalGrid(n_lat=73, n_lon=144)
    orbital_sys = OrbitalSystem()
    forcing = ThermalForcing(grid, orbital_sys)

    # Calculate for a specific time (e.g., one quarter into the orbit)
    time_t = orbital_sys.T_planet / 4.0
    
    # Dummy cloud cover and albedo for testing
    dummy_cloud_cover = np.full(grid.lat_mesh.shape, 0.4)
    dummy_albedo = calculate_dynamic_albedo(dummy_cloud_cover, 288.0, 0.1, 0.6, 0.5)

    insolation_field = forcing.calculate_insolation(time_t)
    temp_eq_field = forcing.calculate_equilibrium_temp(time_t, dummy_albedo)

    print(f"Calculating for time t = {time_t / (3600*24):.1f} days")
    print("Insolation field shape:", insolation_field.shape)
    print(f"Max insolation: {np.max(insolation_field):.2f} W/m^2")
    print("Equilibrium temperature field shape:", temp_eq_field.shape)
    print(f"Max equilibrium temperature: {np.max(temp_eq_field):.2f} K")
    print(f"Min equilibrium temperature: {np.min(temp_eq_field):.2f} K")

    # Check a point on the equator at "noon"
    # Noon is where hour angle h is 0. h = time_of_day_angle + lon_rad
    # So lon_rad = -time_of_day_angle.
    # Let's find the longitude index closest to that.
    time_of_day_angle = (time_t * const.PLANET_OMEGA) % (2 * np.pi)
    noon_lon_rad = -time_of_day_angle
    noon_lon_deg = np.rad2deg(noon_lon_rad) % 360
    
    equator_idx = grid.n_lat // 2
    noon_lon_idx = np.argmin(np.abs(grid.lon - noon_lon_deg))

    print(f"Equator noon temperature: {temp_eq_field[equator_idx, noon_lon_idx]:.2f} K")

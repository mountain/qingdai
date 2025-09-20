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
    def __init__(self, grid: SphericalGrid, orbital_system: OrbitalSystem, base_albedo_map: np.ndarray):
        """
        Initializes the thermal forcing module.

        Args:
            grid (SphericalGrid): The model's grid system.
            orbital_system (OrbitalSystem): The orbital system providing energy flux.
            base_albedo_map (np.ndarray): A 2D array of ice-free surface albedo.
        """
        self.grid = grid
        self.orbital_system = orbital_system
        self.base_albedo_map = base_albedo_map
        self.planet_params = {
            'axial_tilt': const.PLANET_AXIAL_TILT,
            'omega': const.PLANET_OMEGA,
            'T_planet': self.orbital_system.T_planet
        }

    def calculate_insolation(self, t):
        """
        Calculates the instantaneous solar radiation (insolation) at each
        grid point (lat, lon) for a given time t.

        Args:
            t (float): Time in seconds.

        Returns:
            np.ndarray: A 2D array of insolation values (W/m^2) on the grid.
        """
        # 1. Get total flux from the orbital system
        s_total = self.orbital_system.calculate_total_flux(t)

        # 2. Calculate solar declination (delta)
        # This approximates the seasonal change due to axial tilt
        axial_tilt_rad = np.deg2rad(self.planet_params['axial_tilt'])
        # The angle in the orbit (0 to 2*pi)
        orbital_angle = 2 * np.pi * (t / self.planet_params['T_planet'])
        delta = axial_tilt_rad * np.sin(orbital_angle)

        # 3. Calculate hour angle (h) for each longitude
        # This represents the time of day at each longitude
        # We need to define a "day" based on rotation.
        # Let's define time_of_day as fraction of a full rotation
        time_of_day_angle = (t * self.planet_params['omega']) % (2 * np.pi)
        lon_rad = np.deg2rad(self.grid.lon_mesh)
        h = time_of_day_angle + lon_rad

        # 4. Calculate cosine of solar zenith angle (cos(Z))
        lat_rad = np.deg2rad(self.grid.lat_mesh)
        sin_lat = np.sin(lat_rad)
        cos_lat = np.cos(lat_rad)
        sin_delta = np.sin(delta)
        cos_delta = np.cos(delta)

        cos_z = sin_lat * sin_delta + cos_lat * cos_delta * np.cos(h)

        # 5. Clamp to zero on the night side
        cos_z = np.maximum(0, cos_z)

        # 6. Calculate insolation
        insolation = s_total * cos_z
        return insolation

    def _calculate_dynamic_albedo(self, T_s):
        """Calculates the albedo, accounting for ice formation."""
        albedo_ice = 0.6
        freezing_point = 273.15
        # Where surface temperature is below freezing, albedo increases to ice_albedo
        ice_mask = T_s < freezing_point
        dynamic_albedo = np.where(ice_mask, albedo_ice, self.base_albedo_map)
        return dynamic_albedo

    def calculate_equilibrium_temp(self, t, T_s):
        """
        Calculates the global radiative equilibrium temperature field T_eq
        for a given time t, considering the current surface temperature for albedo.

        Args:
            t (float): Time in seconds.
            T_s (np.ndarray): Current surface temperature field (K).

        Returns:
            np.ndarray: A 2D array of equilibrium temperatures (K) on the grid.
        """
        insolation = self.calculate_insolation(t)
        dynamic_albedo = self._calculate_dynamic_albedo(T_s)
        
        # Stefan-Boltzmann Law: I * (1 - albedo) = sigma * T^4
        # T = (I * (1 - albedo) / sigma)^(1/4)
        
        numerator = insolation * (1 - dynamic_albedo)
        # Avoid division by zero or negative roots for the night side
        # Where insolation is zero, temperature should be zero (or a minimum background temp)
        # For now, we'll let it be zero.
        
        # To prevent taking the root of a negative number if numerator is somehow negative
        numerator[numerator < 0] = 0
        
        temp_eq_field = (numerator / const.SIGMA)**0.25
        
        return temp_eq_field

if __name__ == '__main__':
    # Example usage
    grid = SphericalGrid(n_lat=73, n_lon=144)
    orbital_sys = OrbitalSystem()
    # Create a dummy uniform albedo map for testing
    base_albedo_map = np.full(grid.lat_mesh.shape, const.PLANET_ALBEDO)
    forcing = ThermalForcing(grid, orbital_sys, base_albedo_map)

    # Calculate for a specific time (e.g., one quarter into the orbit)
    time_t = orbital_sys.T_planet / 4.0
    
    # Dummy surface temperature for testing
    dummy_T_s = np.full(grid.lat_mesh.shape, 288.0)
    
    insolation_field = forcing.calculate_insolation(time_t)
    temp_eq_field = forcing.calculate_equilibrium_temp(time_t, dummy_T_s)

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

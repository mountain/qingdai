# pygcm/grid.py

"""
Defines the spherical grid system for the planet Qingdai.
"""

import numpy as np
from . import constants

class SphericalGrid:
    """
    Represents a global spherical grid.
    """
    def __init__(self, n_lat, n_lon):
        """
        Initializes a grid with a specified resolution.

        Args:
            n_lat (int): Number of latitude points.
            n_lon (int): Number of longitude points.
        """
        self.n_lat = n_lat
        self.n_lon = n_lon

        # Create 1D arrays for latitude and longitude
        # Latitude from -90 to +90 degrees
        self.lat = np.linspace(-90, 90, n_lat)
        # Longitude from 0 to 360 degrees
        self.lon = np.linspace(0, 360, n_lon)

        # Create a 2D meshgrid for calculations
        self.lon_mesh, self.lat_mesh = np.meshgrid(self.lon, self.lat)

        # Calculate the Coriolis parameter f = 2 * Omega * sin(lat)
        self.coriolis_param = self._calculate_coriolis()

    def _calculate_coriolis(self):
        """
        Calculates the Coriolis parameter for each latitude point.
        """
        lat_rad = np.deg2rad(self.lat_mesh)
        f = 2 * constants.PLANET_OMEGA * np.sin(lat_rad)
        return f

if __name__ == '__main__':
    # Example usage: Create a grid and print some properties
    grid = SphericalGrid(n_lat=73, n_lon=144)
    print(f"Grid created with {grid.n_lat} latitude points and {grid.n_lon} longitude points.")
    print("Latitude points:", grid.lat)
    print("Longitude points:", grid.lon)
    print("Shape of Coriolis parameter array:", grid.coriolis_param.shape)
    print("Coriolis parameter at the equator:", grid.coriolis_param[grid.n_lat // 2, 0])
    print("Coriolis parameter at the north pole:", grid.coriolis_param[-1, 0])

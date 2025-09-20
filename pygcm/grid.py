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

        # Grid spacing in radians
        self.dlat_rad = np.deg2rad(self.lat[1] - self.lat[0])
        self.dlon_rad = np.deg2rad(self.lon[1] - self.lon[0])

    def divergence(self, u, v):
        """
        Calculates the horizontal divergence of a vector field (u, v) on the spherical grid
        using a finite difference method.
        div(V) = (1 / (a * cos(phi))) * [d(u)/d(lambda) + d(v * cos(phi))/d(phi)]
        """
        a = constants.PLANET_RADIUS
        lat_rad = np.deg2rad(self.lat_mesh)
        cos_lat = np.cos(lat_rad)
        
        # Avoid division by zero at the poles
        cos_lat_capped = np.maximum(cos_lat, 1e-6)

        # Use np.roll to calculate finite differences, which handles periodic boundaries correctly
        # d(u)/d(lambda)
        du_dlon = (np.roll(u, -1, axis=1) - np.roll(u, 1, axis=1)) / (2 * self.dlon_rad)

        # d(v * cos(phi))/d(phi)
        v_cos_phi = v * cos_lat
        dv_cos_phi_dlat = (np.roll(v_cos_phi, -1, axis=0) - np.roll(v_cos_phi, 1, axis=0)) / (2 * self.dlat_rad)
        
        # At the poles, the gradient calculation is singular. We can set it to 0.
        dv_cos_phi_dlat[0, :] = 0
        dv_cos_phi_dlat[-1, :] = 0

        div = (1 / (a * cos_lat_capped)) * (du_dlon + dv_cos_phi_dlat)
        
        return div

    def vorticity(self, u, v):
        """
        Calculates the vertical component of vorticity on the spherical grid.
        vort_z = (1 / (a * cos(phi))) * [d(v)/d(lambda) - d(u * cos(phi))/d(phi)]
        """
        a = constants.PLANET_RADIUS
        lat_rad = np.deg2rad(self.lat_mesh)
        cos_lat = np.cos(lat_rad)
        cos_lat_capped = np.maximum(cos_lat, 1e-6)

        dv_dlon = (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1)) / (2 * self.dlon_rad)
        
        u_cos_phi = u * cos_lat
        du_cos_phi_dlat = (np.roll(u_cos_phi, -1, axis=0) - np.roll(u_cos_phi, 1, axis=0)) / (2 * self.dlat_rad)
        du_cos_phi_dlat[0, :] = 0
        du_cos_phi_dlat[-1, :] = 0

        vort = (1 / (a * cos_lat_capped)) * (dv_dlon - du_cos_phi_dlat)
        return vort

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

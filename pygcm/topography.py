# pygcm/topography.py

"""
Defines the topography and surface properties of the planet.
"""

import numpy as np

def create_land_sea_mask(grid):
    """
    Creates a simple land-sea mask with two continents.

    Args:
        grid (SphericalGrid): The model grid.

    Returns:
        np.ndarray: A 2D array where 1 represents land and 0 represents ocean.
    """
    land_mask = np.zeros_like(grid.lat_mesh)

    # Continent 1 (Northern Hemisphere)
    lat_min1, lat_max1 = 10, 60
    lon_min1, lon_max1 = 180, 240
    continent1 = (grid.lat_mesh > lat_min1) & (grid.lat_mesh < lat_max1) & \
                 (grid.lon_mesh > lon_min1) & (grid.lon_mesh < lon_max1)
    land_mask[continent1] = 1

    # Continent 2 (Southern Hemisphere)
    lat_min2, lat_max2 = -50, -15
    lon_min2, lon_max2 = 60, 150
    continent2 = (grid.lat_mesh > lat_min2) & (grid.lat_mesh < lat_max2) & \
                 (grid.lon_mesh > lon_min2) & (grid.lon_mesh < lon_max2)
    land_mask[continent2] = 1
    
    # Adjust to get approximately 40% land coverage
    # This is a rough approximation, a more precise method would integrate cell areas.
    current_land_frac = np.mean(land_mask)
    print(f"Initial land fraction: {current_land_frac:.2%}")

    return land_mask

def generate_base_properties(mask):
    """
    Generates base albedo and friction maps based on the land-sea mask.
    These are the ice-free properties.

    Args:
        mask (np.ndarray): The land-sea mask.

    Returns:
        tuple: A tuple containing (base_albedo_map, friction_map).
    """
    # Albedo: higher for land, lower for ocean
    albedo_land = 0.3
    albedo_ocean = 0.1
    albedo_map = np.where(mask == 1, albedo_land, albedo_ocean)

    # Friction coefficient: higher for land, lower for ocean
    friction_land = 1e-5
    friction_ocean = 1e-6
    friction_map = np.where(mask == 1, friction_land, friction_ocean)
    
    return albedo_map, friction_map

if __name__ == '__main__':
    from grid import SphericalGrid
    import matplotlib.pyplot as plt

    grid = SphericalGrid(n_lat=121, n_lon=240)
    mask = create_land_sea_mask(grid)
    albedo, friction = generate_base_properties(mask)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12))
    ax1.contourf(grid.lon, grid.lat, mask)
    ax1.set_title("Land-Sea Mask")
    ax2.contourf(grid.lon, grid.lat, albedo)
    ax2.set_title("Albedo")
    ax3.contourf(grid.lon, grid.lat, friction)
    ax3.set_title("Friction Coefficient")
    plt.tight_layout()
    plt.show()

"""
physics.py

This module contains parameterizations for physical processes in the Qingdai GCM,
such as precipitation, cloud formation, and their feedback on radiation.
"""
import numpy as np
from scipy.ndimage import gaussian_filter
from . import constants

def diagnose_precipitation(gcm, grid, D_crit, k_precip, cloud_threshold=0.05, smooth_sigma=1.0):
    """
    Diagnoses precipitation based on wind field convergence with a smooth ramp.
    Precipitation forms where convergence is sufficiently negative; optionally
    modulated by cloud cover using a soft (logistic) mask to avoid sharp edges.

    Args:
        gcm (SpectralModel): The GCM object containing state variables.
        grid (SphericalGrid): The model grid.
        D_crit (float): Critical convergence for precipitation (s^-1).
        k_precip (float): Precipitation efficiency coefficient.
        cloud_threshold (float): Cloud fraction at which precipitation begins to be favored.
        smooth_sigma (float): Gaussian smoothing sigma (in grid points) applied to the result.

    Returns:
        np.ndarray: Precipitation rate field.
    """
    div = grid.divergence(gcm.u, gcm.v)

    # Smooth ramp instead of hard threshold: P ~ max(0, -(D - D_crit))
    precip_raw = np.maximum(0.0, -(div - D_crit))
    precip = k_precip * precip_raw

    # Soft cloud gating (logistic), avoids rectangular artifacts from hard masks
    if cloud_threshold is not None and cloud_threshold > 0:
        cc = np.clip(gcm.cloud_cover, 0.0, 1.0)
        sharpness = 10.0  # larger = sharper transition
        mask = 1.0 / (1.0 + np.exp(-sharpness * (cc - cloud_threshold)))
        precip *= mask

    # Gentle smoothing to remove pixel-level/blocky artifacts
    if smooth_sigma and smooth_sigma > 0:
        precip = gaussian_filter(precip, sigma=smooth_sigma)

    return precip

def cloud_from_precip(precip, C_max=0.95, P_ref=2e-5, smooth_sigma=1.0):
    """
    Map precipitation rate to cloud cover using a smooth saturating relation:
        C = C_max * tanh(precip / P_ref)
    where:
      - C_max is the maximum achievable cloud fraction (e.g., 0.9–0.95)
      - P_ref sets the sensitivity: precip ~ P_ref gives C ~ 0.76*C_max
    An optional Gaussian smoothing is applied to avoid pixel/blocky artifacts.

    Args:
        precip (np.ndarray): Precipitation rate field (arbitrary internal units).
        C_max (float): Maximum cloud fraction cap.
        P_ref (float): Reference precipitation scale controlling sensitivity.
        smooth_sigma (float): Gaussian smoothing sigma (grid points).

    Returns:
        np.ndarray: Cloud fraction field in [0, 1].
    """
    eps = 1e-12
    C = C_max * np.tanh(precip / (P_ref + eps))
    if smooth_sigma and smooth_sigma > 0:
        C = gaussian_filter(C, sigma=smooth_sigma)
    return np.clip(C, 0.0, 1.0)

def parameterize_cloud_cover(gcm, grid, land_mask):
    """
    Parameterizes cloud cover from local thermodynamic and dynamic proxies:
    1) Evaporation/condensation proxy from surface temperature (applied everywhere to avoid rectangular masks).
    2) Lifting in cyclonic regions (relative vorticity).
    3) Frontal generation from |temperature advection|.
    Returns a source term in [0, 1] that gets time-integrated externally.
    """
    cloud_source = np.zeros_like(gcm.T_s)

    # --- 1) Evaporation/condensation proxy (thermodynamic) ---
    T_threshold = 285.0
    temp_diff = (gcm.T_s - T_threshold) / 12.0
    evap_source = 0.5 * np.clip(np.tanh(temp_diff), 0.0, 1.0)
    cloud_source += evap_source

    # --- 2) Vorticity Source (Lifting in Cyclones) ---
    vort = grid.vorticity(gcm.u, gcm.v)
    f_safe = grid.coriolis_param + 1e-12
    rel_vort = vort / f_safe
    vort_threshold = 0.5
    vsrc = 0.4 * np.clip(np.tanh((rel_vort - vort_threshold) / 2.0), 0.0, 1.0)
    cloud_source += vsrc

    # --- 3) Frontal Source (Temperature Advection) ---
    T_s = gcm.T_s
    u, v = gcm.u, gcm.v
    dx = grid.dlon_rad * constants.PLANET_RADIUS * np.maximum(1e-6, np.cos(np.deg2rad(grid.lat_mesh)))
    dy = grid.dlat_rad * constants.PLANET_RADIUS

    grad_T_x = (np.roll(T_s, -1, axis=1) - np.roll(T_s, 1, axis=1)) / (2 * dx)
    grad_T_y = (np.roll(T_s, -1, axis=0) - np.roll(T_s, 1, axis=0)) / (2 * dy)
    temp_advection = - (u * grad_T_x + v * grad_T_y)

    frontal_threshold = 2e-5  # K/s
    fsrc = 0.3 * np.clip(np.tanh(np.abs(temp_advection) / frontal_threshold), 0.0, 1.0)
    cloud_source += fsrc

    # Gentle spatial smoothing to avoid blocky/rectangular artifacts
    cloud_source = gaussian_filter(cloud_source, sigma=1.0)

    # Clamp the final source term
    return np.clip(cloud_source, 0.0, 1.0)

def calculate_dynamic_albedo(cloud_cover, T_s, alpha_water, alpha_ice, alpha_cloud):
    """
    Calculates the local albedo based on cloud cover and surface ice.

    Args:
        cloud_cover (np.ndarray): Cloud cover fraction field.
        T_s (np.ndarray): Surface temperature field (K).
        alpha_water (float): Albedo of water/land.
        alpha_ice (float): Albedo of ice.
        alpha_cloud (float): Albedo of clouds.

    Returns:
        np.ndarray: Dynamic albedo field.
    """
    # Determine surface albedo based on temperature (ice-albedo feedback)
    freezing_point = 273.15
    surface_albedo = np.where(T_s < freezing_point, alpha_ice, alpha_water)
    
    # Combine surface albedo with cloud albedo
    albedo = surface_albedo * (1 - cloud_cover) + alpha_cloud * cloud_cover
    return albedo

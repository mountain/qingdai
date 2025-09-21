"""
physics.py

This module contains parameterizations for physical processes in the Qingdai GCM,
such as precipitation, cloud formation, and their feedback on radiation.
"""
import os
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

def compute_orographic_factor(grid, elevation, u, v, k_orog=7e-4, cap=2.0, smooth_sigma=1.0):
    """
    Compute multiplicative orographic precipitation enhancement factor based on upslope wind.
      - n_hat = grad(H)/|grad(H)| (slope unit vector); if |grad(H)|~0, factor=1
      - uplift = max(0, U · n_hat), where U=(u, v)
      - factor = clip(1 + k_orog * uplift, 1, cap), optionally smoothed

    Args:
        grid (SphericalGrid): Grid providing spacings (dlat_rad, dlon_rad) and lat_mesh.
        elevation (np.ndarray): Elevation in meters (lat x lon).
        u, v (np.ndarray): Wind components (m/s).
        k_orog (float): Scaling coefficient for uplift strength.
        cap (float): Upper limit for enhancement factor (>=1).
        smooth_sigma (float): Gaussian sigma (grid points) to smooth factor.

    Returns:
        np.ndarray: Orographic factor >= 1 with same shape as elevation.
    """
    a = constants.PLANET_RADIUS
    lat_rad = np.deg2rad(grid.lat_mesh)
    cos_lat = np.maximum(np.cos(lat_rad), 1e-6)

    dx = a * cos_lat * grid.dlon_rad  # meters per lon step
    dy = a * grid.dlat_rad            # meters per lat step

    # Central differences for surface gradient (m/m)
    dHdx = (np.roll(elevation, -1, axis=1) - np.roll(elevation, 1, axis=1)) / (2.0 * dx)
    dHdy = (np.roll(elevation, -1, axis=0) - np.roll(elevation, 1, axis=0)) / (2.0 * dy)

    # Regularize poles
    dHdy[0, :] = 0.0
    dHdy[-1, :] = 0.0

    grad_norm = np.sqrt(dHdx**2 + dHdy**2)
    eps = 1e-12
    n_x = np.where(grad_norm > eps, dHdx / (grad_norm + eps), 0.0)
    n_y = np.where(grad_norm > eps, dHdy / (grad_norm + eps), 0.0)

    uplift = np.maximum(0.0, u * n_x + v * n_y)  # m/s projected upslope wind
    factor = 1.0 + k_orog * uplift
    factor = np.clip(factor, 1.0, cap)

    if smooth_sigma and smooth_sigma > 0:
        factor = gaussian_filter(factor, sigma=smooth_sigma)

    return factor


def calculate_dynamic_albedo(cloud_cover,
                             T_s,
                             base_albedo,
                             alpha_ice,
                             alpha_cloud,
                             land_mask=None,
                             t_freeze: float = 271.35,
                             delta_T: float = 5.0,
                             ice_only_over_ocean: bool = True,
                             ocean_albedo_threshold: float = 0.15,
                             ice_frac: np.ndarray | None = None,
                             h_ice: np.ndarray | None = None,
                             H_ref: float = 0.5,
                             h0: float = 0.05,
                             gamma: float = 1.0):
    """
    Calculates dynamic albedo. Supports either:
      - A smooth temperature-based sea-ice transition (default), or
      - An externally provided ice fraction (ice_frac in [0,1], e.g., from h_ice).

    Optionally restricts sea-ice formation to ocean points.

    Args:
        cloud_cover (np.ndarray): Cloud cover fraction field in [0,1].
        T_s (np.ndarray or float): Surface temperature field (K), broadcastable.
        base_albedo (float or np.ndarray): Baseline surface albedo (water/land or 2D map).
        alpha_ice (float): Albedo of ice (sea-ice or snow/ice).
        alpha_cloud (float): Albedo of clouds.
        land_mask (np.ndarray, optional): 1=land, 0=ocean. If provided and ice_only_over_ocean=True,
                                          ice is only allowed over ocean points.
        t_freeze (float): Freezing temperature (K) for temperature-based transition center.
        delta_T (float): Transition half-width (K) for temperature-based transition.
        ice_only_over_ocean (bool): If True, ice transition applies only to ocean points.
        ocean_albedo_threshold (float): When land_mask is None but base_albedo is a 2D map,
                                        use base_albedo < threshold to approximate ocean.
        ice_frac (np.ndarray|None): If provided, use this [0,1] field as ice fraction instead of
                                    computing from temperature.

    Returns:
        np.ndarray: Dynamic albedo field in [0,1].
    """
    import numpy as np

    # Ensure arrays
    T_s_arr = np.asarray(T_s, dtype=float)
    C = np.clip(np.asarray(cloud_cover, dtype=float), 0.0, 1.0)

    # Prepare base albedo as array
    if isinstance(base_albedo, np.ndarray):
        base = base_albedo.astype(float)
    else:
        base = np.full_like(T_s_arr, float(base_albedo), dtype=float)

    # Determine ice fraction priority:
    # 1) externally provided ice_frac in [0,1]
    # 2) from thickness h_ice via saturating law with threshold h0 and e-folding H_ref
    # 3) fallback: temperature-based smooth transition around t_freeze
    if ice_frac is not None:
        ice_frac_local = np.clip(np.asarray(ice_frac, dtype=float), 0.0, 1.0)
    elif h_ice is not None:
        h = np.maximum(np.asarray(h_ice, dtype=float) - float(h0), 0.0)
        # Thin ice has weak optical effect; saturates for thick ice
        eff = 1.0 - np.exp(-h / max(1e-6, float(H_ref)))
        # Optional nonlinearity control
        eff = np.clip(eff, 0.0, 1.0) ** float(gamma)
        ice_frac_local = eff
    else:
        eps = max(1e-6, float(delta_T))
        ice_frac_local = 0.5 * (1.0 + np.tanh((t_freeze - T_s_arr) / eps))

    # Optionally limit ice to ocean
    if ice_only_over_ocean:
        if land_mask is not None:
            ocean_mask = (land_mask == 0)
        else:
            if isinstance(base_albedo, np.ndarray):
                ocean_mask = (base < float(ocean_albedo_threshold))
            else:
                ocean_mask = np.ones_like(T_s_arr, dtype=bool)
        ice_frac_local = ice_frac_local * ocean_mask

    # Combine base and ice albedos (pre-cloud)
    surface_albedo = base * (1.0 - ice_frac_local) + float(alpha_ice) * ice_frac_local

    # Mix with cloud albedo
    albedo = surface_albedo * (1.0 - C) + float(alpha_cloud) * C
    return np.clip(albedo, 0.0, 1.0)


def diagnose_precipitation_hybrid(gcm,
                                  grid,
                                  D_crit: float = -1e-7,
                                  k_precip: float = 1.0,
                                  *,
                                  orog_factor: np.ndarray | None = None,
                                  smooth_sigma: float = 1.0,
                                  beta_div: float = 0.4,
                                  renorm: bool = True):
    """
    Humidity-aware precipitation diagnosis (hybrid):
      - Base magnitude from humidity module's condensation P_cond (kg m^-2 s^-1)
      - Spatial redistribution by dynamic convergence and optional orographic factor
      - Optional global renormalization to conserve the total P (= total P_cond)

    Rationale:
      - With q available (P008), P_cond provides physically consistent column condensation.
      - Convergence / upslope provides spatial structure without creating mass from nothing.
      - Renorm keeps ⟨P⟩ consistent with humidity/energy closure (LH_release).

    Args:
        gcm: model state (expects attribute `P_cond_flux_last` from humidity module).
        grid: SphericalGrid for divergence/area weights.
        D_crit: critical convergence (s^-1) to start enhancement.
        k_precip: kept for API compatibility (not used directly; scaling is driven by P_cond).
        orog_factor: multiplicative enhancement factor (>=1) from `compute_orographic_factor`.
        smooth_sigma: Gaussian sigma (grid points) to smooth the final field.
        beta_div: weight of dynamic convergence redistribution (0..1 typical 0.3–0.6).
        renorm: if True, rescale to keep global area-weighted ⟨P⟩ equal to ⟨P_cond⟩.

    Returns:
        precip (np.ndarray): precipitation rate (kg m^-2 s^-1), same shape as grid.
    """
    # 1) Base field from humidity condensation (fallback to old diagnostic if missing)
    P_cond = getattr(gcm, "P_cond_flux_last", None)
    if P_cond is None:
        # Fallback: use legacy convergence-based diagnostic (no cloud gating here)
        return diagnose_precipitation(gcm, grid, D_crit, k_precip, cloud_threshold=None, smooth_sigma=smooth_sigma)

    Pq = np.maximum(0.0, np.asarray(P_cond, dtype=float))

    # 2) Dynamic convergence redistribution factor (normalized)
    div = grid.divergence(gcm.u, gcm.v)
    pos = np.maximum(0.0, -(div - float(D_crit)))  # positive where convergence exceeds threshold
    # Normalize by a robust scale (median of positives) and cap to avoid spikes
    if np.any(pos > 0):
        ppos = pos[pos > 0]
        scale = float(np.median(ppos))
        scale = max(scale, 1e-12)
        F_div = np.clip(pos / scale, 0.0, 5.0)  # cap at 5 for stability
    else:
        F_div = np.zeros_like(Pq)

    # 3) Orographic factor (>=1), optional
    if orog_factor is None:
        F_orog = 1.0
    else:
        F_orog = np.clip(np.asarray(orog_factor, dtype=float), 1.0, 3.0)

    # 4) Compose multiplicative redistribution factor
    F = (1.0 + float(beta_div) * F_div) * F_orog

    # 5) Apply and optionally renormalize to conserve total precipitation
    P_raw = Pq * F

    if renorm:
        # Area-weighted global mean conservation
        w = np.maximum(np.cos(np.deg2rad(grid.lat_mesh)), 0.0)
        num = float(np.sum(Pq * w))
        den = float(np.sum(P_raw * w)) + 1e-20
        s = num / den if den > 0 else 1.0
        P = P_raw * s
    else:
        P = P_raw

    # 6) Gentle smoothing to avoid pixel/blocky artifacts
    if smooth_sigma and smooth_sigma > 0:
        P = gaussian_filter(P, sigma=float(smooth_sigma))

    # 7) Optional fallback/blend to dynamic-convergence precipitation when humidity (P_cond) is too weak.
    #    This ensures equatorial convergence (ITCZ-like) can still produce clouds/rain visually,
    #    without altering the water-closure path (which uses P_cond in hydrology).
    try:
        use_fb = int(os.getenv("QD_P_HYBRID_FALLBACK", "1")) == 1
        PQ_MIN = float(os.getenv("QD_PQ_MIN", "1e-8"))          # kg m^-2 s^-1 (~0.86 mm/day)
        ALPHA_LEG = float(os.getenv("QD_P_BLEND", "0.6"))       # blend weight of legacy dyn-precip
    except Exception:
        use_fb = True
        PQ_MIN = 1e-8
        ALPHA_LEG = 0.6

    if use_fb:
        w = np.maximum(np.cos(np.deg2rad(grid.lat_mesh)), 0.0)
        wsum = float(np.sum(w) + 1e-15)
        Pq_mean = float(np.sum(Pq * w) / wsum)
        if Pq_mean < PQ_MIN:
            # Legacy dynamic precipitation purely from convergence (no cloud gating) for structure
            P_dyn = diagnose_precipitation(gcm, grid, D_crit, k_precip, cloud_threshold=None, smooth_sigma=smooth_sigma)
            # Blend fields to inject ITCZ-like signal when moisture supply is weak
            P = (1.0 - ALPHA_LEG) * P + ALPHA_LEG * P_dyn

    return np.clip(P, 0.0, None)

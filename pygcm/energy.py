"""
energy.py

Explicit energy budget modules for Qingdai GCM (Project 006, Milestone 1).
Implements shortwave/longwave radiation partitioning and a surface energy
tendency integrator. Boundary-layer fluxes are stubbed (SH/LH=0) for M1.

Formulas (single-layer gray atmosphere, simplified):
- Shortwave:
    I = incoming SW at TOA (per grid point)
    R = I * alpha_total                      (planetary reflection)
    SW_atm = I * (A_sw0 + k_sw_cloud * C)    (atmospheric SW absorption)
    SW_sfc = I - R - SW_atm                  (surface SW absorption)

- Longwave:
    eps = clip(eps0 + k_lw_cloud * C, 0, 1)  (effective atmospheric emissivity)
    OLR = eps*sigma*Ta^4 + (1 - eps)*sigma*Ts^4  (to space)
    DLR = eps*sigma*Ta^4                         (downward to surface)
    LW_sfc = DLR - sigma*Ts^4                    (net on surface; typically negative)
    LW_atm = eps*(sigma*Ts^4 - 2*sigma*Ta^4)     (net on atmosphere)

- Surface energy tendency:
    C_s dTs/dt = SW_sfc - LW_sfc - SH - LH
    For M1, SH=LH=0. Optional temperature floor avoids night-side collapse.

Environment parameters (defaults tuned conservatively):
    QD_SW_A0=0.06, QD_SW_KC=0.20
    QD_LW_EPS0=0.70, QD_LW_KC=0.20
    QD_T_FLOOR=150.0   (K)
    QD_CS=2e7          (J/m^2/K)
    QD_ENERGY_DIAG=1   (print diagnostics)
"""

from __future__ import annotations

import os
from dataclasses import dataclass
import numpy as np
from typing import Optional, Tuple, Dict

from . import constants as const


@dataclass
class EnergyParams:
    sw_a0: float = 0.06         # base atmospheric SW absorption
    sw_kc: float = 0.20         # cloud shortwave absorption gain
    lw_eps0: float = 0.70       # base atmospheric emissivity (no-cloud)
    lw_kc: float = 0.20         # cloud longwave enhancement
    t_floor: float = 150.0      # night-side temperature floor (K)
    c_sfc: float = 2.0e7        # surface heat capacity (J/m^2/K)
    diag: bool = True           # enable diagnostics printing


def get_energy_params_from_env() -> EnergyParams:
    def _f(env, default):
        try:
            return float(os.getenv(env, str(default)))
        except Exception:
            return default
    def _i(env, default):
        try:
            return int(os.getenv(env, str(default)))
        except Exception:
            return default
    return EnergyParams(
        sw_a0=_f("QD_SW_A0", 0.06),
        sw_kc=_f("QD_SW_KC", 0.20),
        lw_eps0=_f("QD_LW_EPS0", 0.70),
        lw_kc=_f("QD_LW_KC", 0.20),
        t_floor=_f("QD_T_FLOOR", 150.0),
        c_sfc=_f("QD_CS", 2.0e7),
        diag=(_i("QD_ENERGY_DIAG", 1) == 1),
    )


def shortwave_radiation(I: np.ndarray, albedo: np.ndarray, cloud: np.ndarray, params: EnergyParams):
    """
    Partition TOA shortwave I into reflected (R), atmospheric absorption (SW_atm),
    and surface absorption (SW_sfc).
    Returns: (SW_atm, SW_sfc, R)
    """
    alpha = np.clip(albedo, 0.0, 1.0)
    I_clamped = np.maximum(0.0, I)

    # Reflection
    R = I_clamped * alpha

    # Atmospheric SW absorption (bounded)
    A_sw = params.sw_a0 + params.sw_kc * np.clip(cloud, 0.0, 1.0)
    A_sw = np.clip(A_sw, 0.0, 0.95)
    SW_atm = I_clamped * A_sw

    # Surface SW absorption is residual
    SW_sfc = I_clamped - R - SW_atm
    SW_sfc = np.maximum(0.0, SW_sfc)

    return SW_atm, SW_sfc, R


def longwave_radiation(Ts: np.ndarray, Ta: np.ndarray, cloud: np.ndarray, params: EnergyParams):
    """
    Gray single-layer atmosphere longwave partitioning.
    Returns: (LW_atm, LW_sfc, OLR, DLR, eps_eff)
    """
    sigma = const.SIGMA
    Ts4 = np.maximum(0.0, Ts) ** 4
    Ta4 = np.maximum(0.0, Ta) ** 4

    eps = params.lw_eps0 + params.lw_kc * np.clip(cloud, 0.0, 1.0)
    eps = np.clip(eps, 0.0, 1.0)

    OLR = eps * sigma * Ta4 + (1.0 - eps) * sigma * Ts4
    DLR = eps * sigma * Ta4
    LW_sfc = DLR - sigma * Ts4
    LW_atm = eps * (sigma * Ts4 - 2.0 * sigma * Ta4)

    # Enforce fixed greenhouse factor g = 1 - OLR/(σ Ts^4) if enabled.
    # Definition from user: g := 1 - OLR/(σ T^4). On Earth, g ≈ 0.40.
    # We set OLR_target = (1 - g)*σ Ts^4 and DLR_target = g*σ Ts^4 and adjust LW_sfc accordingly.
    try:
        gh_lock = int(os.getenv("QD_GH_LOCK", "1")) == 1
    except Exception:
        gh_lock = True
    if gh_lock:
        try:
            g_target = float(os.getenv("QD_GH_FACTOR", "0.582"))
        except Exception:
            g_target = 0.582
        Ts4_raw = np.maximum(0.0, Ts) ** 4
        OLR_target = (1.0 - g_target) * sigma * Ts4_raw
        DLR_target = g_target * sigma * Ts4_raw
        OLR = OLR_target
        DLR = DLR_target
        LW_sfc = DLR_target - sigma * Ts4

    return LW_atm, LW_sfc, OLR, DLR, eps


# ---------------- P006 extensions: cloud-optical-aware LW and surface emissivity (optional) ----------------
def surface_emissivity_map(land_mask: np.ndarray,
                           ice_frac: np.ndarray | float) -> np.ndarray:
    """
    Build a per-grid surface emissivity map ε_sfc dependent on surface type.
    Defaults are physically plausible and can be tuned via env:
      QD_EPS_OCEAN (default 0.98), QD_EPS_LAND (0.96), QD_EPS_ICE (0.99)
    """
    import numpy as _np
    eps_ocean = float(os.getenv("QD_EPS_OCEAN", "0.98"))
    eps_land  = float(os.getenv("QD_EPS_LAND",  "0.96"))
    eps_ice   = float(os.getenv("QD_EPS_ICE",   "0.99"))
    land = (land_mask == 1)
    ocean = ~land
    eps = _np.full_like(ice_frac, eps_land, dtype=float)
    eps[ocean] = eps_ocean
    # Where ocean has sea-ice, blend towards ε_ice by optical ice_frac
    eps[ocean] = (1.0 - _np.clip(ice_frac[ocean], 0.0, 1.0)) * eps_ocean + _np.clip(ice_frac[ocean], 0.0, 1.0) * eps_ice
    return _np.nan_to_num(eps, copy=False)


def longwave_radiation_v2(Ts: np.ndarray,
                          Ta: np.ndarray,
                          cloud_eff: np.ndarray,
                          eps_sfc: np.ndarray | float,
                          params: EnergyParams,
                          *,
                          tau0: float | None = None,
                          k_tau: float | None = None):
    """
    Cloud-optical-aware single-layer LW with surface emissivity.
    Idea:
      - clear-sky emissivity: eps_clear = clip(lw_eps0,0,1)
      - cloud emissivity:     eps_cloud = 1 - exp(-k_tau * tau_cloud), tau_cloud = tau0 * cloud_eff
      - effective eps:        eps_eff = 1 - (1 - eps_clear) * (1 - eps_cloud)
      - surface emissivity map ε_sfc (ocean/land/ice) enters σ ε_sfc Ts^4

    Returns: (LW_atm, LW_sfc, OLR, DLR, eps_eff)
    """
    sigma = const.SIGMA
    Ts = np.maximum(0.0, Ts)
    Ta = np.maximum(0.0, Ta)
    Ts4 = Ts ** 4
    Ta4 = Ta ** 4

    # Params
    eps_clear = np.clip(float(params.lw_eps0), 0.0, 1.0)
    tau0 = float(os.getenv("QD_LW_TAU0", "6.0")) if tau0 is None else float(tau0)
    k_tau = float(os.getenv("QD_LW_KTAU", "1.0")) if k_tau is None else float(k_tau)

    # Cloud optical depth and corresponding LW emissivity
    cloud_eff = np.clip(cloud_eff, 0.0, 1.0)
    tau_cloud = tau0 * cloud_eff
    eps_cloud = 1.0 - np.exp(-k_tau * tau_cloud)
    eps_cloud = np.clip(eps_cloud, 0.0, 1.0)

    # Combine clear and cloud contributions to effective atmospheric LW emissivity
    eps_eff = 1.0 - (1.0 - eps_clear) * (1.0 - eps_cloud)

    # Surface emissivity map (can be scalar)
    if np.isscalar(eps_sfc):
        eps_sfc_arr = np.full_like(Ts, float(eps_sfc))
    else:
        eps_sfc_arr = np.clip(np.nan_to_num(eps_sfc, copy=False), 0.0, 1.0)

    # Flux partitions
    # TOA OLR: part from atmosphere + transmitted surface emission
    OLR = eps_eff * sigma * Ta4 + (1.0 - eps_eff) * sigma * eps_sfc_arr * Ts4
    # Downward longwave to surface
    DLR = eps_eff * sigma * Ta4
    # Net on surface
    LW_sfc = DLR - sigma * eps_sfc_arr * Ts4
    # Net on atmosphere: absorbed from surface minus its up+down emission
    LW_atm = eps_eff * (sigma * eps_sfc_arr * Ts4 - 2.0 * sigma * Ta4)

    # Enforce fixed greenhouse factor g = 1 - OLR/(σ Ts^4) if enabled (default g≈0.40).
    # For v2 with surface emissivity, use OLR_target = (1 - g)*σ Ts^4 and
    # set DLR_target = g*σ Ts^4, while keeping the surface emissivity in LW_sfc.
    try:
        gh_lock = int(os.getenv("QD_GH_LOCK", "1")) == 1
    except Exception:
        gh_lock = True
    if gh_lock:
        try:
            g_target = float(os.getenv("QD_GH_FACTOR", "0.582"))
        except Exception:
            g_target = 0.582
        Ts4_raw = np.maximum(0.0, Ts) ** 4
        OLR_target = (1.0 - g_target) * sigma * Ts4_raw
        DLR_target = g_target * sigma * Ts4_raw
        OLR = OLR_target
        DLR = DLR_target
        LW_sfc = DLR_target - sigma * eps_sfc_arr * Ts4

    return LW_atm, LW_sfc, OLR, DLR, eps_eff


def integrate_surface_energy(Ts: np.ndarray,
                             SW_sfc: np.ndarray,
                             LW_sfc: np.ndarray,
                             SH: np.ndarray | float,
                             LH: np.ndarray | float,
                             dt: float,
                             params: EnergyParams) -> np.ndarray:
    """
    One-step explicit update of surface temperature from net surface energy.
    Ts_next = max(T_floor, Ts + dt/C_s * (SW_sfc - LW_sfc - SH - LH))
    """
    # Broadcast SH/LH if scalars
    if np.isscalar(SH):
        SH = np.full_like(Ts, float(SH))
    if np.isscalar(LH):
        LH = np.full_like(Ts, float(LH))

    net = SW_sfc - LW_sfc - SH - LH
    dT = (net / max(1e-12, params.c_sfc)) * dt
    Ts_next = Ts + dT
    # Apply temperature floor to avoid runaway night-side collapse
    Ts_next = np.maximum(params.t_floor, Ts_next)
    # Clean NaNs/Infs
    return np.nan_to_num(Ts_next, copy=False)


def integrate_surface_energy_map(Ts: np.ndarray,
                                 SW_sfc: np.ndarray,
                                 LW_sfc: np.ndarray,
                                 SH: np.ndarray | float,
                                 LH: np.ndarray | float,
                                 dt: float,
                                 C_s_map: np.ndarray,
                                 t_floor: float = 150.0) -> np.ndarray:
    """
    Per-grid heat capacity update of surface temperature:
    Ts_next = max(t_floor, Ts + dt/C_s_map * (SW_sfc - LW_sfc - SH - LH))
    C_s_map is in J m^-2 K^-1, broadcastable to Ts. SH/LH can be scalars.
    """
    # Broadcast SH/LH if scalars
    if np.isscalar(SH):
        SH = np.full_like(Ts, float(SH))
    if np.isscalar(LH):
        LH = np.full_like(Ts, float(LH))

    net = SW_sfc - LW_sfc - SH - LH
    # Avoid division by ~0
    C_s_safe = np.where(np.isfinite(C_s_map) & (C_s_map > 1e3), C_s_map, 1e3)
    dT = (net / C_s_safe) * dt
    Ts_next = Ts + dT
    Ts_next = np.maximum(t_floor, Ts_next)
    return np.nan_to_num(Ts_next, copy=False)


def integrate_surface_energy_with_seaice(Ts: np.ndarray,
                                         SW_sfc: np.ndarray,
                                         LW_sfc: np.ndarray,
                                         SH: np.ndarray | float,
                                         LH: np.ndarray | float,
                                         dt: float,
                                         land_mask: np.ndarray,
                                         h_ice: np.ndarray,
                                         C_s_ocean: float | np.ndarray,
                                         C_s_land: float | np.ndarray,
                                         C_s_ice: float | np.ndarray,
                                         t_freeze: float = 271.35,
                                         rho_i: float = 917.0,
                                         L_f: float = 3.34e5,
                                         t_floor: float = 150.0) -> tuple[np.ndarray, np.ndarray]:
    """
    Minimal sea-ice thermodynamics:
      - Compute Q_net = SW_sfc - LW_sfc - SH - LH.
      - Over ocean:
          * If ice present and Q_net>0: melt first (reduce h_ice), residual heats surface.
          * If Q_net<0 and Ts<=t_freeze+0.5: use energy deficit to freeze (increase h_ice), keep Ts near t_freeze.
      - Over land: standard energy integration.
      - Effective heat capacity: land->C_s_land; ocean ice->C_s_ice; open ocean->C_s_ocean.

    Args and Returns:
      Inputs are 2D arrays (lat x lon), scalars broadcast. Returns (Ts_next, h_ice_next).
    """
    import numpy as np

    # Broadcast SH/LH if scalars
    if np.isscalar(SH):
        SH_arr = np.full_like(Ts, float(SH))
    else:
        SH_arr = SH
    if np.isscalar(LH):
        LH_arr = np.full_like(Ts, float(LH))
    else:
        LH_arr = LH

    Q_net = SW_sfc - LW_sfc - SH_arr - LH_arr  # W/m^2

    # Prepare masks
    land = (land_mask == 1)
    ocean = ~land

    Ts_next = Ts.astype(float).copy()
    h_ice_next = h_ice.astype(float).copy()

    # Ensure capacities as arrays
    def _to_arr(v):
        return v if isinstance(v, np.ndarray) else np.full_like(Ts, float(v))
    Cs_ocean_arr = _to_arr(C_s_ocean)
    Cs_land_arr = _to_arr(C_s_land)
    Cs_ice_arr = _to_arr(C_s_ice)

    # Melt first where ice present and heating available
    ice_present = (h_ice_next > 0.0) & ocean
    melt_mask = ice_present & (Q_net > 0.0)
    if np.any(melt_mask):
        melt_energy = Q_net[melt_mask] * dt  # J/m^2 available for melt
        dh_melt = melt_energy / (rho_i * L_f)
        # Cap by available thickness
        dh_cap = np.minimum(dh_melt, h_ice_next[melt_mask])
        h_ice_next[melt_mask] -= dh_cap
        # Remove used energy from Q_net
        Q_net[melt_mask] = Q_net[melt_mask] - (dh_cap * rho_i * L_f) / dt

    # Freeze where cooling and near/below freezing (ocean)
    freeze_tol = 0.5  # K buffer to allow near-freezing transitions
    freeze_mask = ocean & (Q_net < 0.0) & (Ts_next <= (t_freeze + freeze_tol))
    if np.any(freeze_mask):
        freeze_energy = -Q_net[freeze_mask] * dt  # J/m^2 to be converted to ice
        dh_freeze = freeze_energy / (rho_i * L_f)
        h_ice_next[freeze_mask] += dh_freeze
        # Energy fully consumed by phase change
        Q_net[freeze_mask] = 0.0
        # Keep surface near freezing when forming/growing ice
        Ts_next[freeze_mask] = np.minimum(Ts_next[freeze_mask], t_freeze)

    # Effective heat capacity for residual energy update
    Cs_eff = np.where(land, Cs_land_arr, np.where(h_ice_next > 0.0, Cs_ice_arr, Cs_ocean_arr))
    Cs_eff = np.where(np.isfinite(Cs_eff) & (Cs_eff > 1e3), Cs_eff, 1e3)

    # Apply residual energy to temperature
    Ts_next = Ts_next + (Q_net / Cs_eff) * dt

    # South/North-pole hard constraint to suppress polar artifact on regular lat-lon grids:
    # If at the polar ring an ocean point is net-cooling (Q_net<0) but Ts remains above freezing,
    # force Ts to the freezing point. This mitigates numerical warm-pole amplification by ring-averaging.
    try:
        _fix_s = int(os.getenv("QD_POLAR_FREEZE_FIX", "1")) == 1
    except Exception:
        _fix_s = True
    try:
        _fix_n = int(os.getenv("QD_POLAR_FREEZE_FIX_N", "1")) == 1
    except Exception:
        _fix_n = True
    # South pole (row 0)
    if _fix_s:
        try:
            _j = 0
            _ocean = ocean[_j, :]
            _cool = Q_net[_j, :] < 0.0
            _above = Ts_next[_j, :] > t_freeze
            _mask = _ocean & _cool & _above
            if np.any(_mask):
                Ts_next[_j, _mask] = t_freeze
        except Exception:
            # Fail-safe: never crash the timestep on fix failure
            pass
    # North pole (last row)
    if _fix_n:
        try:
            _j = -1
            _ocean = ocean[_j, :]
            _cool = Q_net[_j, :] < 0.0
            _above = Ts_next[_j, :] > t_freeze
            _mask = _ocean & _cool & _above
            if np.any(_mask):
                Ts_next[_j, _mask] = t_freeze
        except Exception:
            # Fail-safe: never crash the timestep on fix failure
            pass

    # Clamp: ice surface not exceeding freezing; global temperature floor
    Ts_next = np.where((h_ice_next > 0.0) & ocean, np.minimum(Ts_next, t_freeze), Ts_next)
    Ts_next = np.maximum(t_floor, Ts_next)
    Ts_next = np.nan_to_num(Ts_next, copy=False)
    h_ice_next = np.nan_to_num(h_ice_next, copy=False)
    return Ts_next, h_ice_next


def boundary_layer_fluxes(Ts: np.ndarray,
                          Ta: np.ndarray,
                          u: np.ndarray,
                          v: np.ndarray,
                          land_mask: np.ndarray,
                          C_H: float = 1.5e-3,
                          rho: float = 1.2,
                          c_p: float = 1004.0,
                          B_land: float = 0.7,
                          B_ocean: float = 0.3):
    """
    Bulk formula for sensible heat; latent via Bowen ratio.
    For M1, this helper is provided but not invoked (SH=LH=0).
    Returns: (SH, LH)
    """
    # Wind speed magnitude
    V = np.sqrt(u**2 + v**2)

    # Sensible heat flux (positive upward from surface)
    SH = rho * c_p * C_H * V * (Ts - Ta)

    # Bowen ratio by surface type
    B = np.where(land_mask == 1, B_land, B_ocean)
    # Avoid division by zero: clamp B
    B = np.maximum(B, 1e-3)
    LH = SH / B
    return SH, LH


def compute_atmos_height_tendency(SW_atm: np.ndarray,
                                  LW_atm: np.ndarray,
                                  SH_from_sfc: np.ndarray | float,
                                  LH_release: np.ndarray | float,
                                  rho_air: float,
                                  H_atm: float,
                                  g: float = 9.81) -> np.ndarray:
    """
    Convert atmospheric net energy flux (W/m^2) into tendency of geopotential height h (m/s),
    using single-layer column mass rho_air*H_atm and hydrostatic scaling by g.

    dh/dt = (SW_atm + LW_atm + SH_from_sfc + LH_release) / (rho_air * H_atm * g)
    """
    # Broadcast SH/LH if scalars
    if np.isscalar(SH_from_sfc):
        SH_from_sfc = np.full_like(SW_atm, float(SH_from_sfc))
    if np.isscalar(LH_release):
        LH_release = np.full_like(SW_atm, float(LH_release))
    F_atm = SW_atm + LW_atm + SH_from_sfc + LH_release
    denom = max(1e-6, float(rho_air)) * max(1.0, float(H_atm)) * float(g)
    return F_atm / denom


def integrate_atmos_energy_height(h: np.ndarray,
                                  SW_atm: np.ndarray,
                                  LW_atm: np.ndarray,
                                  SH_from_sfc: np.ndarray | float,
                                  LH_release: np.ndarray | float,
                                  dt: float,
                                  rho_air: float,
                                  H_atm: float,
                                  g: float = 9.81,
                                  weight: float = 1.0) -> np.ndarray:
    """
    Advance geopotential height h by atmospheric energy tendency over dt.
    The optional 'weight' allows partial coupling (e.g., QD_ENERGY_W).
    """
    dh_dt = compute_atmos_height_tendency(SW_atm, LW_atm, SH_from_sfc, LH_release, rho_air, H_atm, g)
    h_next = h + float(weight) * dh_dt * dt
    return np.nan_to_num(h_next, copy=False)


def compute_energy_diagnostics(lat_mesh: np.ndarray,
                               I: np.ndarray,
                               R: np.ndarray,
                               OLR: np.ndarray,
                               SW_sfc: np.ndarray,
                               LW_sfc: np.ndarray,
                               SH: np.ndarray | float = 0.0,
                               LH: np.ndarray | float = 0.0):
    """
    Area-weighted global means for energy budgets:
      TOA_net = I - R - OLR
      SFC_net = SW_sfc - LW_sfc - SH - LH
      ATM_net = TOA_net - SFC_net
    Returns dict of global means.
    """
    # Broadcast scalars
    if np.isscalar(SH):
        SH = np.full_like(SW_sfc, float(SH))
    if np.isscalar(LH):
        LH = np.full_like(SW_sfc, float(LH))

    TOA_net = I - R - OLR
    SFC_net = SW_sfc - LW_sfc - SH - LH
    ATM_net = TOA_net - SFC_net

    # Area weights ~ cos(lat)
    w = np.cos(np.deg2rad(lat_mesh))
    w = np.maximum(w, 0.0)
    w_sum = np.sum(w)

    def wmean(x):
        return float(np.sum(x * w) / (w_sum + 1e-15))

    return {
        "TOA_net": wmean(TOA_net),
        "SFC_net": wmean(SFC_net),
        "ATM_net": wmean(ATM_net),
        "I_mean": wmean(I),
        "R_mean": wmean(R),
        "OLR_mean": wmean(OLR),
        "SW_sfc_mean": wmean(SW_sfc),
        "LW_sfc_mean": wmean(LW_sfc),
        "SH_mean": wmean(SH),
        "LH_mean": wmean(LH),
    }


# ---------------- Greenhouse autotuning (diagnostic-driven) ----------------


def autotune_greenhouse_params(params, diag, rate_eps=None, rate_kc=None,
                               bounds_eps=(0.30, 0.98), bounds_kc=(0.0, 0.80),
                               verbose=None):
    """
    Nudge greenhouse coefficients based on TOA_net to approach global energy balance.
    Heuristic controller using the sign that dOLR/deps < 0 in gray one-layer model
    (increasing eps reduces OLR if Ts^4 > Ta^4). Thus:
      - If TOA_net > 0 (planet gaining energy), increase OLR -> decrease eps/kc.
      - If TOA_net < 0, decrease OLR -> increase eps/kc.

    Args:
      params: EnergyParams instance (modified in place and also returned).
      diag: dict from compute_energy_diagnostics (must contain "TOA_net" in W/m^2).
      rate_eps: step size for lw_eps0 (default from env QD_TUNE_RATE_EPS or 5e-5 per step).
      rate_kc: step size for lw_kc   (default from env QD_TUNE_RATE_KC  or 2e-5 per step).
      bounds_eps/kc: clipping ranges to keep parameters in plausible intervals.
      verbose: if True prints tuning line (default follows QD_ENERGY_AUTOTUNE_DIAG).

    Returns:
      params (same object) for convenience.
    """
    import numpy as _np
    err = float(diag.get("TOA_net", 0.0))  # W/m^2; target 0
    rate_eps = float(os.getenv("QD_TUNE_RATE_EPS", "5e-5")) if rate_eps is None else float(rate_eps)
    rate_kc  = float(os.getenv("QD_TUNE_RATE_KC",  "2e-5")) if rate_kc  is None else float(rate_kc)
    verbose = (int(os.getenv("QD_ENERGY_AUTOTUNE_DIAG", "1")) == 1) if verbose is None else bool(verbose)

    # Controller: eps_new = eps_old - k * err
    # Positive err (gain energy) -> decrease eps/kc to raise OLR
    params.lw_eps0 = float(_np.clip(params.lw_eps0 - rate_eps * err, bounds_eps[0], bounds_eps[1]))
    params.lw_kc   = float(_np.clip(params.lw_kc   - rate_kc  * err, bounds_kc[0], bounds_kc[1]))

    if verbose:
        print(f"[EnergyTune] TOA_net={err:+.3f} W/m^2 -> eps0={params.lw_eps0:.3f}, kc={params.lw_kc:.3f}")

    return params

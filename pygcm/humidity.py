"""
humidity.py

Project 008: Single-layer atmospheric humidity (q) and E–P–LH coupling.

This module provides:
- HumidityParams: configuration loaded from environment variables.
- q_sat(T, p): saturation specific humidity (kg/kg) using Tetens formula.
- q_init(Ts, RH0, p0): initialize near-surface q from relative humidity and surface temperature.
- surface_evaporation_factor(land_mask, h_ice, params): per-grid surface factor for evaporation.
- evaporation_flux(Ts, q, u, v, factor, params): E (kg m^-2 s^-1) via bulk aerodynamic formula.
- condensation(q, T_a, dt, params): supersaturation sink; returns (P_cond_flux, q_next).

Conventions:
- q: specific humidity (kg/kg), grid 2D.
- E_flux: upward mass flux of water vapor (kg/m^2/s) from surface to atmosphere.
- P_cond_flux: condensation/precip mass flux from air to surface (kg/m^2/s).
- LH (surface) = L_v * E_flux (W/m^2). Positive upward (surface energy sink).
- LH_release (atmos) = L_v * P_cond_flux (W/m^2). Positive heating of atmosphere.

Units:
- rho_a: kg/m^3, h_mbl: m → column mass M_col = rho_a * h_mbl (kg/m^2).
- dq/dt from surface source/sink:
    dq/dt |_evap = + E_flux / M_col
    dq/dt |_cond = - P_cond_flux / M_col
"""

from __future__ import annotations

from dataclasses import dataclass
import os
import numpy as np

EPSILON = 0.622  # ratio of molecular weights Mw/Md for moist/dry air


@dataclass
class HumidityParams:
    # Bulk aero and column properties
    C_E: float = 1.3e-3
    rho_a: float = 1.2       # kg/m^3 (near-surface)
    h_mbl: float = 800.0     # m (effective mixed-layer height)
    L_v: float = 2.5e6       # J/kg (latent heat of vaporization)
    p0: float = 1.0e5        # Pa (reference pressure for q_sat)

    # Surface evaporation scaling by type
    ocean_evap_scale: float = 1.0
    land_evap_scale: float = 0.5
    ice_evap_scale: float = 0.05

    # Microphysics/relaxation
    tau_cond: float = 1800.0  # s; timescale for condensation relaxation

    # Numerical/diagnostics
    diag: bool = True


def get_humidity_params_from_env() -> HumidityParams:
    def _f(env: str, default: float) -> float:
        try:
            return float(os.getenv(env, str(default)))
        except Exception:
            return default

    def _i(env: str, default: int) -> int:
        try:
            return int(os.getenv(env, str(default)))
        except Exception:
            return default

    return HumidityParams(
        C_E=_f("QD_CE", 1.3e-3),
        rho_a=_f("QD_RHO_A", 1.2),
        h_mbl=_f("QD_MBL_H", 800.0),
        L_v=_f("QD_LV", 2.5e6),
        p0=_f("QD_P0", 1.0e5),
        ocean_evap_scale=_f("QD_OCEAN_EVAP_SCALE", 1.0),
        land_evap_scale=_f("QD_LAND_EVAP_SCALE", 0.5),
        ice_evap_scale=_f("QD_ICE_EVAP_SCALE", 0.05),
        tau_cond=_f("QD_TAU_COND", 1800.0),
        diag=(_i("QD_HUMIDITY_DIAG", 1) == 1),
    )


def q_sat(T: np.ndarray | float, p: float = 1.0e5) -> np.ndarray:
    """
    Saturation specific humidity over liquid water using Tetens formula.
    Args:
        T: temperature in K (array or scalar).
        p: ambient pressure in Pa (assumed constant near surface).
    Returns:
        q_sat in kg/kg, same shape as T.
    """
    T_arr = np.asarray(T, dtype=float)
    T_c = np.clip(T_arr - 273.15, -80.0, 60.0)  # Celsius for formula stability
    # Tetens (over water), e_s in Pa
    e_s = 610.94 * np.exp(17.625 * T_c / (T_c + 243.04))
    # Specific humidity from vapor pressure
    denom = np.maximum(p - (1.0 - EPSILON) * e_s, 1.0)  # avoid division by ~0
    qsat = EPSILON * e_s / denom
    return np.clip(qsat, 0.0, 0.5)  # physical upper bound


def q_init(Ts: np.ndarray, RH0: float = 0.5, p0: float = 1.0e5) -> np.ndarray:
    """
    Initialize near-surface specific humidity by relative humidity against surface temperature.
    Args:
        Ts: surface temperature (K), 2D array.
        RH0: initial relative humidity (0..1).
        p0: surface pressure (Pa).
    """
    RH = float(np.clip(RH0, 0.0, 1.0))
    return RH * q_sat(Ts, p=p0)


def surface_evaporation_factor(land_mask: np.ndarray | None,
                               h_ice: np.ndarray | None,
                               params: HumidityParams,
                               ice_threshold: float = 1e-6) -> np.ndarray:
    """
    Construct per-grid surface factor S_type for evaporation:
      - open ocean: ocean_evap_scale
      - sea ice (ocean & h_ice > threshold): ice_evap_scale
      - land: land_evap_scale
    If land_mask is None, assume all ocean (factor=ocean_evap_scale, no ice differentiation).
    """
    if land_mask is None:
        base = np.ones_like(h_ice if h_ice is not None else 0.0) if isinstance(h_ice, np.ndarray) else 1.0
        return np.full_like(base, params.ocean_evap_scale, dtype=float)

    land = (land_mask == 1)
    ocean = ~land
    factor = np.zeros_like(land_mask, dtype=float)
    if h_ice is not None:
        ice = (h_ice > float(ice_threshold)) & ocean
        open_ocean = ocean & (~ice)
        factor[ice] = float(params.ice_evap_scale)
        factor[open_ocean] = float(params.ocean_evap_scale)
    else:
        factor[ocean] = float(params.ocean_evap_scale)
    factor[land] = float(params.land_evap_scale)
    return factor


def evaporation_flux(Ts: np.ndarray,
                     q: np.ndarray,
                     u: np.ndarray,
                     v: np.ndarray,
                     surface_factor: np.ndarray,
                     params: HumidityParams) -> np.ndarray:
    """
    Bulk aerodynamic evaporation mass flux (kg/m^2/s):
        E = rho_a * C_E * |V| * (q_sat(Ts) - q)_+ * S_type
    """
    V = np.sqrt(u**2 + v**2)
    qsat_sfc = q_sat(Ts, p=params.p0)
    deficit = np.maximum(0.0, qsat_sfc - q)
    E_flux = params.rho_a * params.C_E * V * deficit * surface_factor
    return np.nan_to_num(E_flux, copy=False)


def condensation(q: np.ndarray,
                 T_a: np.ndarray,
                 dt: float,
                 params: HumidityParams) -> tuple[np.ndarray, np.ndarray]:
    """
    Supersaturation relaxation to saturation over timescale tau_cond:
      excess = max(0, q - q_sat(T_a))
      P_cond_flux = (excess / tau_cond) * (rho_a * h_mbl)   [kg/m^2/s]
      q_next = q - (P_cond_flux / (rho_a * h_mbl)) * dt     [kg/kg]
    Returns:
      (P_cond_flux, q_next)
    """
    qsat_air = q_sat(T_a, p=params.p0)
    excess = np.maximum(0.0, q - qsat_air)
    # Convert to mass flux using column mass
    M_col = max(1e-6, float(params.rho_a * params.h_mbl))
    P_cond_flux = (excess / max(1e-6, float(params.tau_cond))) * M_col
    q_next = q - (P_cond_flux / M_col) * dt
    # Numerical hygiene
    q_next = np.clip(np.nan_to_num(q_next, copy=False), 0.0, 0.5)
    P_cond_flux = np.nan_to_num(P_cond_flux, copy=False)
    return P_cond_flux, q_next

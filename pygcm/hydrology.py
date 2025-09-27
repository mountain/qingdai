"""
hydrology.py

Project 009: Planetary hydrology closure (E–P–R and surface reservoirs).

This module provides:
- HydrologyParams loaded from environment variables
- partition_precip_phase: split precip into rain/snow by surface temperature
- snow_step: update land snow reservoir and compute melt flux
- update_land_bucket: update land water reservoir (bucket) and compute runoff
- diagnose_water_closure: global diagnostics for water mass closure

Conventions and units:
- All fluxes E/P/R/melt are kg m^-2 s^-1 (1 kg m^-2 s^-1 = 1 mm s^-1 of water)
- Reservoirs W_land and S_snow are kg m^-2 (equivalent to mm of liquid water)
- Ice reservoir mass is rho_i * h_ice (kg m^-2)
- CWV mass is (rho_a * h_mbl * q) (kg m^-2)
"""

from __future__ import annotations

from dataclasses import dataclass
import os
import numpy as np


@dataclass
class HydrologyParams:
    runoff_tau_days: float = 10.0         # days, linear runoff timescale
    wland_cap_mm: float | None = None     # mm; if set, capacity of land bucket (excess -> fast runoff)
    snow_thresh_K: float = 273.15         # K; T threshold for snow vs rain (used with T_hat_a if lapse enabled)
    snow_melt_rate_mm_day: float = 5.0    # mm/day equivalent melt rate (constant-mode fallback)
    rho_w: float = 1000.0                 # kg/m^3, reference fresh water density
    # P019 additions
    snow_t_band_K: float = 1.5            # K; sigmoid half-width for smooth rain/snow split
    snow_melt_mode: str = "degree_day"    # 'degree_day' | 'constant'
    snow_ddf_mm_per_k_day: float = 3.0    # degree-day factor (mm/K/day)
    snow_melt_tref_K: float = 273.15      # K; melt reference temperature for DDF
    swe_enable: bool = True               # enable snowpack reservoir (SWE)
    swe_ref_mm: float = 15.0              # mm; reference SWE for optical coverage C_snow
    swe_max_mm: float | None = None       # optional cap for SWE (mm)
    diag: bool = True                     # enable diagnostics printing

def get_hydrology_params_from_env() -> HydrologyParams:
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
    cap = os.getenv("QD_WLAND_CAP", "")
    try:
        wland_cap_mm = float(cap) if cap not in ("", "None", "none", "null") else None
    except Exception:
        wland_cap_mm = None
    swe_max_env = os.getenv("QD_SWE_MAX_MM", "")
    try:
        swe_max_mm = float(swe_max_env) if swe_max_env not in ("", "None", "none", "null") else None
    except Exception:
        swe_max_mm = None
    return HydrologyParams(
        runoff_tau_days=_f("QD_RUNOFF_TAU_DAYS", 10.0),
        wland_cap_mm=wland_cap_mm,
        snow_thresh_K=_f("QD_SNOW_THRESH", 273.15),
        snow_melt_rate_mm_day=_f("QD_SNOW_MELT_RATE", 5.0),
        rho_w=_f("QD_RHO_W", 1000.0),
        # P019
        snow_t_band_K=_f("QD_SNOW_T_BAND", 1.5),
        snow_melt_mode=os.getenv("QD_SNOW_MELT_MODE", "degree_day").strip().lower(),
        snow_ddf_mm_per_k_day=_f("QD_SNOW_DDF_MM_PER_K_DAY", 3.0),
        snow_melt_tref_K=_f("QD_SNOW_MELT_TREF", 273.15),
        swe_enable=(_i("QD_SWE_ENABLE", 1) == 1),
        swe_ref_mm=_f("QD_SWE_REF_MM", 15.0),
        swe_max_mm=swe_max_mm,
        diag=(_i("QD_WATER_DIAG", 1) == 1),
    )


def partition_precip_phase(P_flux: np.ndarray, T_s: np.ndarray, T_thresh: float = 273.15) -> tuple[np.ndarray, np.ndarray]:
    """
    Split total precipitation mass flux into rain and snow by a temperature threshold.
    Args:
      P_flux: total precip flux (kg m^-2 s^-1)
      T_s: surface temperature (K)
      T_thresh: threshold (K): T < T_thresh -> snow, else rain
    Returns:
      (P_rain_flux, P_snow_flux) in kg m^-2 s^-1
    """
    P_flux = np.asarray(P_flux, dtype=float)
    T_s = np.asarray(T_s, dtype=float)
    snow_mask = (T_s < float(T_thresh))
    P_snow = np.where(snow_mask, P_flux, 0.0)
    P_rain = np.where(snow_mask, 0.0, P_flux)
    return np.nan_to_num(P_rain, copy=False), np.nan_to_num(P_snow, copy=False)

def partition_precip_phase_smooth(P_flux: np.ndarray,
                                  T_hat_a: np.ndarray,
                                  T_thresh: float = 273.15,
                                  dT_half_K: float = 1.5) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Smoothly split total precipitation into rain/snow using a sigmoid in temperature space.
    Args:
      P_flux: total precip flux (kg m^-2 s^-1)
      T_hat_a: lapse-adjusted air temperature (K)
      T_thresh: central threshold (K)
      dT_half_K: sigmoid half-width ΔT (K)
    Returns:
      (P_rain_flux, P_snow_flux, f_snow) where f_snow ∈ [0,1]
    """
    P_flux = np.asarray(P_flux, dtype=float)
    T_hat_a = np.asarray(T_hat_a, dtype=float)
    # f_snow = 1 / (1 + exp((T_hat_a - T_thresh)/ΔT))
    f_snow = 1.0 / (1.0 + np.exp((T_hat_a - float(T_thresh)) / max(1e-6, float(dT_half_K))))
    f_snow = np.clip(f_snow, 0.0, 1.0)
    P_snow = f_snow * P_flux
    P_rain = (1.0 - f_snow) * P_flux
    return np.nan_to_num(P_rain, copy=False), np.nan_to_num(P_snow, copy=False), f_snow

def snowpack_step(S_snow: np.ndarray,
                  P_snow_land: np.ndarray,
                  T_hat_a: np.ndarray,
                  params: HydrologyParams,
                  dt: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Update snowpack SWE and compute melt flux using degree-day or constant mode.
    Also returns C_snow (optical coverage) and a simple snow albedo (fresh/old not tracked here).
    Args:
      S_snow: current SWE (kg m^-2 ≡ mm)
      P_snow_land: snowfall mass flux on land (kg m^-2 s^-1)
      T_hat_a: lapse-adjusted air temperature (K)
      params: HydrologyParams
      dt: step (s)
    Returns:
      (S_next, melt_flux, C_snow, alpha_snow_map)
    """
    S = np.asarray(S_snow, dtype=float).copy()
    Ps = np.asarray(P_snow_land, dtype=float)
    Ta = np.asarray(T_hat_a, dtype=float)

    # Melt rate (kg m^-2 s^-1)
    if str(params.snow_melt_mode).lower() == "degree_day":
        # DDF [mm/K/day] → kg m^-2 s^-1 per K
        ddf_per_s = float(params.snow_ddf_mm_per_k_day) / 86400.0
        positive = np.maximum(Ta - float(params.snow_melt_tref_K), 0.0)
        melt_flux = ddf_per_s * positive  # kg m^-2 s^-1
    else:
        # constant mode (backward-compat)
        melt_flux = float(params.snow_melt_rate_mm_day) / 86400.0
        melt_flux = np.where(Ta >= float(params.snow_thresh_K), melt_flux, 0.0)

    # Potential melt amount over dt, cap by available SWE
    pot_melt_amt = melt_flux * dt
    actual_melt_amt = np.minimum(np.maximum(S, 0.0), pot_melt_amt)
    S_next = S + Ps * dt - actual_melt_amt
    if params.swe_max_mm is not None and params.swe_max_mm > 0:
        S_next = np.minimum(S_next, float(params.swe_max_mm))
    S_next = np.maximum(0.0, S_next)

    # Flux back (kg m^-2 s^-1)
    melt_flux_out = np.where(dt > 0, actual_melt_amt / dt, 0.0)

    # Optical snow cover fraction
    C_snow = 1.0 - np.exp(-np.maximum(S_next, 0.0) / max(1e-6, float(params.swe_ref_mm)))
    C_snow = np.clip(C_snow, 0.0, 1.0)

    # Simple snow albedo map (fresh snow value; aging not tracked here)
    try:
        alpha_snow_val = float(os.getenv("QD_SNOW_ALBEDO_FRESH", "0.70"))
    except Exception:
        alpha_snow_val = 0.70
    alpha_snow_map = np.full_like(S_next, alpha_snow_val, dtype=float)

    return np.nan_to_num(S_next, copy=False), np.nan_to_num(melt_flux_out, copy=False), C_snow, alpha_snow_map


def snow_step(S_snow: np.ndarray,
              P_snow_land: np.ndarray,
              T_s: np.ndarray,
              params: HydrologyParams,
              dt: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Update land snow reservoir and compute melt flux (kg m^-2 s^-1).
    Melt occurs when T_s >= snow_thresh at a fixed rate (mm/day -> kg m^-2 s^-1),
    capped by available snow.

    Args:
      S_snow: current land snow water-equivalent (kg m^-2)
      P_snow_land: snowfall mass flux on land (kg m^-2 s^-1)
      T_s: surface temperature (K)
      params: HydrologyParams
      dt: time step (s)
    Returns:
      (S_next, melt_flux) where melt_flux is kg m^-2 s^-1 on land (else 0)
    """
    S = np.asarray(S_snow, dtype=float).copy()
    P_snow_land = np.asarray(P_snow_land, dtype=float)
    T_s = np.asarray(T_s, dtype=float)

    # Convert mm/day -> kg m^-2 s^-1
    melt_rate = (float(params.snow_melt_rate_mm_day) / 86400.0)  # mm/s
    melt_rate_kg = melt_rate  # 1 mm water over 1 m^2 == 1 kg, so numeric same

    melt_mask = (T_s >= float(params.snow_thresh_K))
    # Potential melt amount over dt
    potential_melt = np.where(melt_mask, melt_rate_kg, 0.0) * dt  # kg m^-2 over dt
    # Cap by available snow
    actual_melt_amt = np.minimum(np.maximum(S, 0.0), potential_melt)
    S_next = S + P_snow_land * dt - actual_melt_amt
    S_next = np.maximum(0.0, S_next)
    # Convert actual melt back to flux (kg m^-2 s^-1)
    melt_flux = np.where(dt > 0, actual_melt_amt / dt, 0.0)
    return np.nan_to_num(S_next, copy=False), np.nan_to_num(melt_flux, copy=False)


def update_land_bucket(W_land: np.ndarray,
                       P_in: np.ndarray,
                       E_land: np.ndarray,
                       params: HydrologyParams,
                       dt: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Update land bucket water storage and compute runoff.
    Linear runoff R = W / tau. Optional capacity: overflow goes to immediate runoff.

    Args:
      W_land: current land water storage (kg m^-2)
      P_in: input water flux to bucket (rain + snowmelt on land), kg m^-2 s^-1
      E_land: land evaporation mass flux (kg m^-2 s^-1)
      params: hydrology params
      dt: time step (s)
    Returns:
      (W_next, R_flux) where R_flux is kg m^-2 s^-1
    """
    W = np.asarray(W_land, dtype=float).copy()
    P_in = np.asarray(P_in, dtype=float)
    E_land = np.asarray(E_land, dtype=float)

    tau_s = max(1.0, float(params.runoff_tau_days) * 86400.0)

    # Base runoff (linear reservoir) computed on previous storage
    R_base = W / tau_s

    # Tend; apply explicit Euler: dW/dt = P_in - E_land - R_base
    W_next = W + (P_in - E_land - R_base) * dt
    W_next = np.maximum(0.0, W_next)

    # Optional capacity overflow as fast runoff this step
    if params.wland_cap_mm is not None and params.wland_cap_mm > 0:
        cap = float(params.wland_cap_mm)
        overflow = np.maximum(0.0, W_next - cap)
        W_next = W_next - overflow
        R_fast = np.where(dt > 0, overflow / dt, 0.0)
    else:
        R_fast = 0.0

    R_flux = R_base + R_fast
    return np.nan_to_num(W_next, copy=False), np.nan_to_num(R_flux, copy=False)


def _area_weights(lat_mesh: np.ndarray) -> np.ndarray:
    w = np.maximum(np.cos(np.deg2rad(lat_mesh)), 0.0)
    return w

def _wmean(x: np.ndarray, w: np.ndarray) -> float:
    return float(np.sum(x * w) / (np.sum(w) + 1e-15))

def diagnose_water_closure(lat_mesh: np.ndarray,
                           q: np.ndarray,
                           rho_a: float,
                           h_mbl: float,
                           h_ice: np.ndarray,
                           rho_i: float,
                           W_land: np.ndarray,
                           S_snow: np.ndarray,
                           E_flux: np.ndarray,
                           P_flux: np.ndarray,
                           R_flux: np.ndarray,
                           dt_since_prev: float | None,
                           prev_total: float | None) -> dict:
    """
    Compute area-weighted global means and closure residual:
      Let reservoirs (kg m^-2):
        CWV = rho_a * h_mbl * q
        ICE = rho_i * h_ice
        W_land, S_snow as given
      Then closure (global mean):
        d/dt <CWV + ICE + W_land + S_snow> ?= <E> - <P> - <R>

    Args:
      lat_mesh: grid lat mesh (deg)
      q: specific humidity (kg/kg)
      rho_a, h_mbl: for CWV mass
      h_ice, rho_i: for ice mass
      W_land, S_snow: land water/snow reservoirs (kg m^-2)
      E_flux, P_flux, R_flux: fluxes (kg m^-2 s^-1), sign: E upward (adds CWV), P downward (removes CWV)
      dt_since_prev: seconds since last diagnostic step (for d/dt estimate); if None, derivative omitted
      prev_total: previous global mean total reservoir; if None, derivative omitted
    Returns:
      dict with means and optional derivative/residual.
    """
    w = _area_weights(lat_mesh)

    # Reservoir masses (kg m^-2)
    CWV = float(rho_a) * float(h_mbl) * np.asarray(q, dtype=float)
    ICE = float(rho_i) * np.asarray(h_ice, dtype=float)
    W = np.asarray(W_land, dtype=float)
    S = np.asarray(S_snow, dtype=float)

    # Means
    CWV_mean = _wmean(CWV, w)
    ICE_mean = _wmean(ICE, w)
    W_mean = _wmean(W, w)
    S_mean = _wmean(S, w)
    E_mean = _wmean(np.asarray(E_flux, dtype=float), w)
    P_mean = _wmean(np.asarray(P_flux, dtype=float), w)
    R_mean = _wmean(np.asarray(R_flux, dtype=float), w)

    total_now = CWV_mean + ICE_mean + W_mean + S_mean

    result = {
        "CWV_mean": CWV_mean,
        "ICE_mean": ICE_mean,
        "W_land_mean": W_mean,
        "S_snow_mean": S_mean,
        "E_mean": E_mean,
        "P_mean": P_mean,
        "R_mean": R_mean,
        "total_reservoir_mean": total_now,
    }

    if (dt_since_prev is not None) and (prev_total is not None) and dt_since_prev > 0:
        ddt_total = (total_now - prev_total) / float(dt_since_prev)
        residual = ddt_total - (E_mean - P_mean - R_mean)
        result["d/dt_total_mean"] = ddt_total
        result["closure_residual"] = residual

    return result

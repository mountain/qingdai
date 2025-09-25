from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass
class WeatherInstant:
    """
    Instantaneous ecological micro-environment for sub-daily updates.
    Minimal M1 fields.
    """
    Ts: float | np.ndarray
    Ta: float | np.ndarray
    wind10: float | np.ndarray
    soil_water_index: float | np.ndarray
    I_bands: np.ndarray  # shape [NB, n_lat, n_lon] or [NB] if global
    cloud_eff: float | np.ndarray = 0.0


@dataclass
class WeatherDaily:
    """
    Daily aggregated environment (M1 optional).
    """
    Ts_mean: float | np.ndarray
    Ta_mean: float | np.ndarray
    wind10_mean: float | np.ndarray
    soil_water_index: float | np.ndarray
    I_bands_mean: np.ndarray  # shape [NB, n_lat, n_lon] or [NB]
    precip_daily: float | np.ndarray = 0.0

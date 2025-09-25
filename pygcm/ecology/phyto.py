from __future__ import annotations

import os
from dataclasses import dataclass
import numpy as np

from .spectral import (
    make_bands,
    dual_star_insolation_to_bands,
    band_weights_from_mode,
    SpectralBands,
)


@dataclass
class PhytoParams:
    # Growth parameters (daily model)
    mu_max: float = 1.5              # d^-1 maximum potential growth
    alpha_P: float = 0.04            # 1/(W m^-2) light utilization coefficient
    Q10: float = 2.0                 # temperature sensitivity
    T_ref: float = 293.15            # K (20°C) reference
    m0: float = 0.05                 # d^-1 background loss (respiration/mortality)
    lambda_sink_m_per_day: float = 0.0  # m/day equivalent sinking (0 for M1)

    # Optics
    kd_exp_m: float = 0.5            # exponent in Kd ~ Chl^m

    # Initialization
    chl0: float = 0.05               # mg/m^3 initial mixed-layer chlorophyll


def _read_env_float(name: str, default: float) -> float:
    v = os.getenv(name)
    if v is None:
        return default
    try:
        return float(v)
    except Exception:
        return default


def _read_env_bool(name: str, default: bool) -> bool:
    v = os.getenv(name)
    if v is None:
        return default
    try:
        return bool(int(v))
    except Exception:
        return default


def _read_env_list(name: str) -> list[float] | None:
    v = os.getenv(name)
    if not v:
        return None
    try:
        parts = [p.strip() for p in v.split(",")]
        out = [float(p) for p in parts if p != ""]
        return out if len(out) > 0 else None
    except Exception:
        return None


def _nearest_index(arr: np.ndarray, value: float) -> int:
    return int(np.argmin(np.abs(arr - value)))


class PhytoManager:
    """
    Minimal mixed-layer phytoplankton (M1) with daily time step:
      - Prognostic chlorophyll concentration C_phyto (mg Chl m^-3) in mixed layer
      - Light limitation via tanh; temperature via Q10; constant background loss
      - Optics across NB shortwave bands:
          Kd_b(Chl) = Kd0_b + k_chl_b * Chl^m
          I_b,MLD   = I_b(surf) * (1 - exp(-Kd_b H)) / (Kd_b H)
          A_b^water = A_b^pure + c_reflect * S_b * Chl^p_reflect, clipped [0,1]
        Reduced to scalar α_water_eff by band weights (Rayleigh/simple).

    Inputs expected for daily stepping:
      - insA, insB: per-star surface shortwave (W m^-2)
      - T_w: mixed-layer water temperature (K) (use SST)
      - dt_days: length of daily step in days (default 1)
    """
    def __init__(
        self,
        grid,
        land_mask: np.ndarray,
        bands: SpectralBands | None = None,
        H_mld_m: float | None = None,
        diag: bool = True,
    ) -> None:
        self.grid = grid
        self.land_mask = (land_mask.astype(int))
        self.ocean_mask = (self.land_mask == 0)
        self.NL, self.NM = self.grid.n_lat, self.grid.n_lon

        # Bands
        self.bands: SpectralBands = bands or make_bands()
        NB = self.bands.nbands

        # Params
        self.params = PhytoParams(
            mu_max=_read_env_float("QD_PHYTO_MU_MAX", 1.5),
            alpha_P=_read_env_float("QD_PHYTO_ALPHA_P", 0.04),
            Q10=_read_env_float("QD_PHYTO_Q10", 2.0),
            T_ref=_read_env_float("QD_PHYTO_T_REF", 293.15),
            m0=_read_env_float("QD_PHYTO_M_LOSS", 0.05),
            lambda_sink_m_per_day=_read_env_float("QD_PHYTO_LAMBDA_SINK", 0.0),
            kd_exp_m=_read_env_float("QD_PHYTO_KD_EXP_M", 0.5),
            chl0=_read_env_float("QD_PHYTO_CHL0", 0.05),
        )
        self.diag = diag

        # Mixed-layer depth (m)
        if H_mld_m is None:
            try:
                H_mld_m = float(os.getenv("QD_OCEAN_H_M", os.getenv("QD_MLD_M", "50")))
            except Exception:
                H_mld_m = 50.0
        self.H_mld = float(max(0.1, H_mld_m))

        # Default/Env optical parameters per band
        # - Base diffuse attenuation Kd0 (m^-1) (clear water)
        # - Chlorophyll contribution factor k_chl (m^-1 per mg^m m^-3)
        # - Pure water base reflectance per band A_b^pure (dimensionless)
        kd0_list = _read_env_list("QD_PHYTO_KD0")          # optional CSV length NB
        kchl_list = _read_env_list("QD_PHYTO_KD_CHL")      # optional CSV length NB
        Apure_list = _read_env_list("QD_PHYTO_APURE")      # optional CSV length NB

        kd0_default = _read_env_float("QD_PHYTO_KD0_DEFAULT", 0.04)
        kchl_default = _read_env_float("QD_PHYTO_KD_CHL_DEFAULT", 0.02)
        Apure_default = _read_env_float("QD_PHYTO_APURE_DEFAULT", 0.06)

        self.Kd0_b = np.full((NB,), kd0_default, dtype=float)
        self.kchl_b = np.full((NB,), kchl_default, dtype=float)
        self.Apure_b = np.full((NB,), Apure_default, dtype=float)

        if kd0_list is not None:
            for i, val in enumerate(kd0_list[:NB]):
                self.Kd0_b[i] = float(val)
        if kchl_list is not None:
            for i, val in enumerate(kchl_list[:NB]):
                self.kchl_b[i] = float(val)
        if Apure_list is not None:
            for i, val in enumerate(Apure_list[:NB]):
                self.Apure_b[i] = float(val)

        # Spectral shape for reflectance (green-enhanced), normalized
        lam = self.bands.lambda_centers
        mu_g = _read_env_float("QD_PHYTO_SHAPE_MU_NM", 550.0)
        sig_g = _read_env_float("QD_PHYTO_SHAPE_SIGMA_NM", 70.0)
        shape = np.exp(-((lam - mu_g) ** 2) / (2.0 * sig_g ** 2))
        ssum = float(np.sum(shape)) + 1e-12
        self.shape_b = shape / ssum

        self.c_reflect = _read_env_float("QD_PHYTO_REFLECT_C", 0.02)
        self.p_reflect = _read_env_float("QD_PHYTO_REFLECT_P", 0.5)
        self.alpha_clip_min = _read_env_float("QD_PHYTO_ALPHA_MIN", 0.0)
        self.alpha_clip_max = _read_env_float("QD_PHYTO_ALPHA_MAX", 1.0)

        # Band weights for reducing to scalar alpha (Rayleigh/simple)
        self.w_b = band_weights_from_mode(self.bands)

        # Prognostic fields
        self.C_phyto = np.zeros((self.NL, self.NM), dtype=float)
        self.C_phyto[self.ocean_mask] = self.params.chl0
        # Diagnostics
        self.alpha_water_scalar = np.zeros((self.NL, self.NM), dtype=float)
        self.alpha_water_bands: np.ndarray | None = None
        self.Kd_490 = np.zeros((self.NL, self.NM), dtype=float)
        self._idx_490 = _nearest_index(self.bands.lambda_centers, 490.0)

        if self.diag:
            print(f"[Phyto] NB={NB} bands, H_mld={self.H_mld:.1f} m | mu_max={self.params.mu_max:.2f} d^-1 | "
                  f"alpha_P={self.params.alpha_P:.3f} | m0={self.params.m0:.3f} d^-1 | Q10={self.params.Q10:.2f}")

    # ---------- Core optics helpers ----------

    def _kd_bands(self, chl_map: np.ndarray) -> np.ndarray:
        """
        Compute band diffuse attenuation Kd_b[NB, NL, NM] from chlorophyll field.
        Kd_b = Kd0_b + k_chl_b * chl^m
        """
        NB = self.bands.nbands
        Kd = np.zeros((NB, self.NL, self.NM), dtype=float)
        # chl exponent m
        m = float(self.params.kd_exp_m)
        chl_pow = np.power(np.maximum(chl_map, 0.0), m)
        for b in range(NB):
            Kd[b, :, :] = self.Kd0_b[b] + self.kchl_b[b] * chl_pow
        # Avoid zero/negatives
        Kd = np.clip(Kd, 1e-6, np.inf)
        return Kd

    def _Ibar_bands_in_mld(self, I_b_surf: np.ndarray, Kd_b: np.ndarray) -> np.ndarray:
        """
        Vertically averaged band irradiance in mixed layer:
          Ī_b = I_b * (1 - exp(-Kd_b H)) / (Kd_b H)
        """
        NB = self.bands.nbands
        H = self.H_mld
        # Safe factor (1 - e^-x)/x with series fallback for small x
        x = Kd_b * H
        out = np.zeros_like(I_b_surf)
        small = (x < 1e-6)
        if np.any(small):
            out[:, small] = I_b_surf[:, small] * (1.0 - 0.5 * x[small] + (x[small] ** 2) / 6.0)
        big = ~small
        if np.any(big):
            out[:, big] = I_b_surf[:, big] * ((1.0 - np.exp(-x[big])) / x[big])
        # Non-negative
        return np.clip(out, 0.0, np.inf)

    def _alpha_bands_from_chl(self, chl_map: np.ndarray) -> np.ndarray:
        """
        Compute water band reflectance A_b^water[NB, NL, NM] from chlorophyll.
        A_b = A_pure + c_reflect * S_b * chl^p
        """
        NB = self.bands.nbands
        chl_pow = np.power(np.maximum(chl_map, 0.0), self.p_reflect)
        A = np.zeros((NB, self.NL, self.NM), dtype=float)
        for b in range(NB):
            A[b, :, :] = self.Apure_b[b] + self.c_reflect * self.shape_b[b] * chl_pow
        return np.clip(A, self.alpha_clip_min, self.alpha_clip_max)

    # ---------- Public interface ----------

    def step_daily(
        self,
        insA: np.ndarray,
        insB: np.ndarray,
        T_w: np.ndarray,
        dt_days: float = 1.0,
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Advance phytoplankton by one daily step and update alpha maps.
        Returns (alpha_bands, alpha_scalar).
        """
        # Only operate on ocean cells; land remains unchanged
        # 1) Bands at surface from dual-star shortwave
        I_b_surf = dual_star_insolation_to_bands(insA, insB, self.bands)  # [NB, NL, NM]

        # 2) Mixed-layer average irradiance per band
        Kd_b = self._kd_bands(self.C_phyto)
        Ibar_b = self._Ibar_bands_in_mld(I_b_surf, Kd_b)

        # 3) PAR proxy: sum over bands (bands already in visible/near-visible)
        I_PAR = np.sum(Ibar_b, axis=0)  # [NL, NM]

        # 4) Growth modifiers
        muL = np.tanh(self.params.alpha_P * I_PAR / max(self.params.mu_max, 1e-6))
        fT = np.power(self.params.Q10, (np.asarray(T_w, dtype=float) - self.params.T_ref) / 10.0)

        # 5) Net specific rate (d^-1)
        sink_term = 0.0
        if self.params.lambda_sink_m_per_day > 0.0:
            sink_term = float(self.params.lambda_sink_m_per_day) / max(1e-6, self.H_mld)
        mu = self.params.mu_max * muL * fT - self.params.m0 - sink_term

        # 6) Update chlorophyll
        dC = mu * self.C_phyto * float(dt_days)
        self.C_phyto = np.clip(self.C_phyto + dC, 0.0, np.inf)
        # Keep land cells at zero
        self.C_phyto[~self.ocean_mask] = 0.0

        # 7) Update optics → albedo maps
        alpha_b = self._alpha_bands_from_chl(self.C_phyto)  # [NB, NL, NM]
        alpha_scalar = np.sum(alpha_b * self.w_b[:, None, None], axis=0)
        alpha_scalar = np.clip(alpha_scalar, self.alpha_clip_min, self.alpha_clip_max)

        self.alpha_water_bands = alpha_b
        self.alpha_water_scalar = alpha_scalar

        # 8) Diagnostics (Kd490)
        self.Kd_490 = Kd_b[self._idx_490, :, :]

        if self.diag:
            try:
                w = np.maximum(np.cos(np.deg2rad(self.grid.lat_mesh)), 0.0)
                wsum = float(np.sum(w)) + 1e-15

                def wmean(x: np.ndarray) -> float:
                    return float(np.sum(np.nan_to_num(x) * w) / wsum)

                print(
                    f"[PhytoDiag] ⟨Chl⟩={wmean(self.C_phyto):.3f} mg/m^3 | "
                    f"⟨Kd490⟩={wmean(self.Kd_490):.3f} m^-1 | "
                    f"⟨α_water⟩={wmean(self.alpha_water_scalar):.3f}"
                )
            except Exception:
                pass

        return self.alpha_water_bands, self.alpha_water_scalar

    def get_alpha_maps(self) -> tuple[np.ndarray | None, np.ndarray]:
        """
        Return (alpha_water_bands, alpha_water_scalar).
        Bands may be None if step_daily has not been called yet.
        """
        return self.alpha_water_bands, self.alpha_water_scalar

    def get_kd490(self) -> np.ndarray:
        """
        Return current Kd(490) field (m^-1).
        """
        return self.Kd_490

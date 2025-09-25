from __future__ import annotations

import os
from dataclasses import dataclass
import numpy as np
from ..jax_compat import (
    is_enabled as _jax_enabled,
    to_numpy as _to_np,
    laplacian_sphere as _j_lap,
    advect_semilag as _j_adv,
)

from .spectral import (
    make_bands,
    dual_star_insolation_to_bands,
    band_weights_from_mode,
    SpectralBands,
)


@dataclass
class PhytoParams:
    # Growth parameters (daily model; can be overridden per-species via env arrays)
    mu_max: float = 1.5              # d^-1 maximum potential growth
    alpha_P: float = 0.04            # 1/(W m^-2) light utilization coefficient
    Q10: float = 2.0                 # temperature sensitivity
    T_ref: float = 293.15            # K (20°C) reference
    m0: float = 0.05                 # d^-1 background loss (respiration/mortality)
    lambda_sink_m_per_day: float = 0.0  # m/day equivalent sinking (0 for M1)

    # Optics
    kd_exp_m: float = 0.5            # exponent in Kd ~ Chl^m

    # Initialization
    chl0: float = 0.05               # mg/m^3 initial mixed-layer chlorophyll (total)


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
    Mixed-layer phytoplankton (M1) with daily time step and multi-species support.

    State:
      - C_phyto_s[S, NL, NM]: per-species chlorophyll (mg Chl m^-3) in the mixed layer
      - alpha_water_bands[NB, NL, NM], alpha_water_scalar[NL, NM]
      - Kd_490[NL, NM]

    Growth (daily):
      - Light limitation via tanh; temperature via Q10; constant background loss
      - Single mixed-layer Kd_b computed from TOTAL chlorophyll C_tot = sum_s C_s
      - Optional species-specific growth rates via env arrays; defaults shared

    Optics:
      - For each band b, water reflectance:
          A_b^water = A_pure_b + sum_s c_reflect_s * Shape_s[b] * (Chl_s)^{p_reflect_s}
        where Shape_s[b] is a Gaussian (center mu_s, width sigma_s) normalized over bands.
      - Scalar α_water_eff by band weights (Rayleigh/simple)
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
        # Grid metrics for transport (spherical geometry)
        self.a = getattr(self.grid, "a", None)
        if self.a is None:
            # Planet radius from constants if grid doesn't carry it
            try:
                from .. import constants as _const
                self.a = _const.PLANET_RADIUS
            except Exception:
                self.a = 6.371e6
        self.dlat = self.grid.dlat_rad
        self.dlon = self.grid.dlon_rad
        self.lat_rad = np.deg2rad(self.grid.lat_mesh)
        # Guard cosφ at high latitudes to avoid metric blow-up
        self.coslat = np.maximum(np.cos(self.lat_rad), 0.5)

        # Horizontal mixing (m^2/s) for tracer; default to ocean K_h if unspecified
        try:
            self.K_h = float(os.getenv("QD_PHYTO_KH", os.getenv("QD_KH_OCEAN", "5.0e3")))
        except Exception:
            self.K_h = 5.0e3

        # Bands
        self.bands: SpectralBands = bands or make_bands()
        NB = self.bands.nbands

        # Params (shared defaults)
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

        # Number of species (default 1; set to 10 to enable ten types)
        S_default = 1
        try:
            S_default = int(os.getenv("QD_PHYTO_NSPECIES", "1"))
        except Exception:
            S_default = 1
        self.S = max(1, S_default)

        # Optical base parameters per band (Kd and pure water reflectance)
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

        # Per-species spectral shapes (Gaussian) and reflectance coefficients
        # Defaults: centers linearly spaced from 460..680 nm; sigma=70 nm; c_reflect=0.02; p_reflect=0.5
        lam = self.bands.lambda_centers  # [NB]
        mu_arr = _read_env_list("QD_PHYTO_SPEC_MU_NM") or []
        sigma_arr = _read_env_list("QD_PHYTO_SPEC_SIGMA_NM") or []
        c_reflect_arr = _read_env_list("QD_PHYTO_SPEC_C_REFLECT") or []
        p_reflect_arr = _read_env_list("QD_PHYTO_SPEC_P_REFLECT") or []
        mu_default_start = 460.0
        mu_default_end = 680.0
        if self.S > 1:
            mu_defaults = np.linspace(mu_default_start, mu_default_end, self.S)
        else:
            mu_defaults = np.array([_read_env_float("QD_PHYTO_SHAPE_MU_NM", 550.0)])
        sigma_default = _read_env_float("QD_PHYTO_SHAPE_SIGMA_NM", 70.0)
        c_reflect_default = _read_env_float("QD_PHYTO_REFLECT_C", 0.02)
        p_reflect_default = _read_env_float("QD_PHYTO_REFLECT_P", 0.5)

        # Initialize per-species arrays
        self.shape_sb = np.zeros((self.S, NB), dtype=float)     # [S, NB]
        self.c_reflect_s = np.zeros((self.S,), dtype=float)      # [S]
        self.p_reflect_s = np.zeros((self.S,), dtype=float)      # [S]
        for s in range(self.S):
            mu_s = mu_arr[s] if s < len(mu_arr) else float(mu_defaults[min(s, len(mu_defaults)-1)])
            sigma_s = sigma_arr[s] if s < len(sigma_arr) else sigma_default
            # Gaussian shape and normalize over bands
            g = np.exp(-((lam - mu_s) ** 2) / (2.0 * sigma_s ** 2))
            gsum = float(np.sum(g)) + 1e-12
            self.shape_sb[s, :] = g / gsum
            # Reflectance controls
            self.c_reflect_s[s] = c_reflect_arr[s] if s < len(c_reflect_arr) else c_reflect_default
            self.p_reflect_s[s] = p_reflect_arr[s] if s < len(p_reflect_arr) else p_reflect_default

        # Clip bounds for alpha
        self.alpha_clip_min = _read_env_float("QD_PHYTO_ALPHA_MIN", 0.0)
        self.alpha_clip_max = _read_env_float("QD_PHYTO_ALPHA_MAX", 1.0)

        # Band weights for reducing to scalar alpha (Rayleigh/simple)
        self.w_b = band_weights_from_mode(self.bands)  # [NB]

        # Species-level growth overrides (optional)
        mu_max_arr = _read_env_list("QD_PHYTO_SPEC_MU_MAX") or []
        m0_arr = _read_env_list("QD_PHYTO_SPEC_M0") or []
        self.mu_max_s = np.array(
            [(mu_max_arr[s] if s < len(mu_max_arr) else self.params.mu_max) for s in range(self.S)],
            dtype=float
        )
        self.m0_s = np.array(
            [(m0_arr[s] if s < len(m0_arr) else self.params.m0) for s in range(self.S)],
            dtype=float
        )

        # Species initial fractions (sum to 1); if not provided, equal split
        frac_arr = _read_env_list("QD_PHYTO_INIT_FRAC") or []
        if len(frac_arr) >= self.S:
            frac = np.clip(np.array(frac_arr[:self.S], dtype=float), 0.0, None)
            s = float(np.sum(frac))
            frac = frac / s if s > 0 else np.full((self.S,), 1.0 / self.S, dtype=float)
        else:
            frac = np.full((self.S,), 1.0 / self.S, dtype=float)
        self.init_frac_s = frac  # [S]

        # Prognostic fields
        # C_phyto_s: [S, NL, NM]; initialize over ocean only with fractions*chl0
        self.C_phyto_s = np.zeros((self.S, self.NL, self.NM), dtype=float)
        for s in range(self.S):
            self.C_phyto_s[s, :, :] = self.init_frac_s[s] * self.params.chl0
        # Land cells zeroed
        for s in range(self.S):
            self.C_phyto_s[s, ~self.ocean_mask] = 0.0

        # Diagnostics
        self.alpha_water_scalar = np.zeros((self.NL, self.NM), dtype=float)
        self.alpha_water_bands: np.ndarray | None = None
        self.Kd_490 = np.zeros((self.NL, self.NM), dtype=float)
        self._idx_490 = _nearest_index(self.bands.lambda_centers, 490.0)

        if self.diag:
            spec_info = f"S={self.S}, mu={self.mu_max_s.min():.2f}..{self.mu_max_s.max():.2f}/d"
            print(f"[Phyto] NB={NB} bands, H_mld={self.H_mld:.1f} m | {spec_info} | alpha_P={self.params.alpha_P:.3f} | m0={self.params.m0:.3f} d^-1 | Q10={self.params.Q10:.2f}")

    # ---------- Core optics helpers ----------

    def _kd_bands(self, chl_total: np.ndarray) -> np.ndarray:
        """
        Compute band diffuse attenuation Kd_b[NB, NL, NM] from TOTAL chlorophyll.
        Kd_b = Kd0_b + k_chl_b * chl_total^m
        """
        NB = self.bands.nbands
        Kd = np.zeros((NB, self.NL, self.NM), dtype=float)
        m = float(self.params.kd_exp_m)
        chl_pow = np.power(np.maximum(chl_total, 0.0), m)  # [NL, NM]
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
        H = self.H_mld
        x = Kd_b * H
        # Safe factor (1 - e^-x)/x with series fallback for small x
        small = x < 1e-6
        factor_small = 1.0 - 0.5 * x + (x ** 2) / 6.0
        factor_big = (1.0 - np.exp(-x)) / np.clip(x, 1e-12, None)
        factor = np.where(small, factor_small, factor_big)
        # Non-negative
        return np.clip(I_b_surf * factor, 0.0, np.inf)

    def _alpha_bands_from_species(self, C_phyto_s: np.ndarray) -> np.ndarray:
        """
        Compute water band reflectance A_b^water[NB, NL, NM] from per-species chlorophyll.
        A_b = A_pure_b + sum_s c_reflect_s * Shape_s[b] * (Chl_s)^p_s
        """
        S = self.S
        NB = self.bands.nbands
        # Start from pure water reflectance broadcast
        A = np.broadcast_to(self.Apure_b[:, None, None], (NB, self.NL, self.NM)).astype(float)
        # Add species contributions
        for s in range(S):
            chl_s = np.maximum(C_phyto_s[s, :, :], 0.0)  # [NL, NM]
            p = float(self.p_reflect_s[s])
            if p == 1.0:
                term_map = chl_s
            else:
                term_map = np.power(chl_s, p)
            coeff = float(self.c_reflect_s[s])
            shape_b = self.shape_sb[s, :]  # [NB]
            # Add for all bands with broadcasting: [NB,1,1] * [1,NL,NM]
            A += (coeff * shape_b[:, None, None]) * term_map[None, :, :]
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
        # 1) Bands at surface from dual-star shortwave
        I_b_surf = dual_star_insolation_to_bands(insA, insB, self.bands)  # [NB, NL, NM]

        # 2) Mixed-layer average irradiance per band (Kd from total chlorophyll)
        C_tot = np.sum(self.C_phyto_s, axis=0)  # [NL, NM]
        Kd_b = self._kd_bands(C_tot)
        Ibar_b = self._Ibar_bands_in_mld(I_b_surf, Kd_b)

        # 3) PAR proxy: sum over bands
        I_PAR = np.sum(Ibar_b, axis=0)  # [NL, NM]

        # 4) Growth modifiers
        muL = np.tanh(self.params.alpha_P * I_PAR / max(1e-6, float(np.max(self.mu_max_s))))
        fT = np.power(self.params.Q10, (np.asarray(T_w, dtype=float) - self.params.T_ref) / 10.0)

        # 5) Net specific rate per SPECIES (d^-1)
        sink_term = 0.0
        if self.params.lambda_sink_m_per_day > 0.0:
            sink_term = float(self.params.lambda_sink_m_per_day) / max(1e-6, self.H_mld)
        # Broadcast to [S, NL, NM]
        mu_s = (self.mu_max_s[:, None, None] * muL[None, :, :] * fT[None, :, :]) - (self.m0_s[:, None, None] + sink_term)

        # 6) Update chlorophyll per species
        dC_s = mu_s * self.C_phyto_s * float(dt_days)
        self.C_phyto_s = np.clip(self.C_phyto_s + dC_s, 0.0, np.inf)
        # Keep land cells at zero
        for s in range(self.S):
            self.C_phyto_s[s, ~self.ocean_mask] = 0.0

        # 7) Update optics → water albedo maps from species mixture
        alpha_b = self._alpha_bands_from_species(self.C_phyto_s)  # [NB, NL, NM]
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

                C_tot_now = np.sum(self.C_phyto_s, axis=0)
                print(
                    f"[PhytoDiag] S={self.S} | ⟨Chl_tot⟩={wmean(C_tot_now):.3f} mg/m^3 | "
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

    # ---------- Ocean current advection & lateral diffusion (M3 transport) ----------

    def _laplacian_sphere(self, F: np.ndarray) -> np.ndarray:
        """
        ∇²F on a regular lat-lon grid using divergence form with cosφ metric.
        Mirrors ocean implementation with optional JAX kernel.
        """
        try:
            if _jax_enabled():
                return _to_np(_j_lap(F, self.dlat, self.dlon, self.coslat, self.a))
        except Exception:
            pass

        F = np.nan_to_num(F, copy=False)
        dF_dphi = np.gradient(F, self.dlat, axis=0)
        term_phi = (1.0 / self.coslat) * np.gradient(self.coslat * dF_dphi, self.dlat, axis=0)
        d2F_dlam2 = (np.roll(F, -1, axis=1) - 2.0 * F + np.roll(F, 1, axis=1)) / (self.dlon ** 2)
        term_lam = d2F_dlam2 / (self.coslat ** 2)
        return (term_phi + term_lam) / (self.a ** 2)

    def _advect_scalar(self, field: np.ndarray, u: np.ndarray, v: np.ndarray, dt: float) -> np.ndarray:
        """
        Semi-Lagrangian advection (bilinear interpolation, lon periodic).
        Uses JAX kernel if available, falls back to scipy.ndimage.map_coordinates.
        """
        try:
            if _jax_enabled():
                adv = _j_adv(field, u, v, dt, self.a, self.dlat, self.dlon, self.coslat)
                return _to_np(adv)
        except Exception:
            pass

        from scipy.ndimage import map_coordinates
        dlam = u * dt / (self.a * self.coslat)   # radians of longitude
        dphi = v * dt / self.a                   # radians of latitude

        dx = dlam / self.dlon
        dy = dphi / self.dlat

        JJ, II = np.meshgrid(np.arange(self.NL), np.arange(self.NM), indexing="ij")
        dep_J = JJ - dy
        dep_I = II - dx

        adv = map_coordinates(field, [dep_J, dep_I], order=1, mode="wrap", prefilter=False)
        return adv

    def advect_diffuse(self, uo: np.ndarray, vo: np.ndarray, dt_seconds: float) -> None:
        """
        Advect each species chlorophyll by ocean surface currents and apply lateral diffusion.

        Args:
            uo, vo: ocean currents (m/s) on the same grid (east/north components).
            dt_seconds: physics step (s).
        Notes:
            - Land cells remain zero (enforced post-update).
            - Uses semi-Lagrangian advection (stable for large dt); lateral diffusion is explicit.
        """
        if dt_seconds <= 0.0:
            return
        S = int(self.S)
        ocean_mask = (self.land_mask == 0)

        for s in range(S):
            C = np.asarray(self.C_phyto_s[s, :, :], dtype=float)
            # Advection
            C_adv = self._advect_scalar(C, uo, vo, float(dt_seconds))
            # Gentle blend toward advected field for stability (similar to SST path)
            adv_alpha = float(os.getenv("QD_PHYTO_ADV_ALPHA", "0.7"))
            C_new = (1.0 - adv_alpha) * C + adv_alpha * C_adv

            # Lateral diffusion (explicit)
            if self.K_h > 0.0:
                C_new = np.nan_to_num(C_new)
                C_new += float(dt_seconds) * self.K_h * self._laplacian_sphere(C_new)

            # Enforce non-negativity and land mask
            C_new = np.clip(C_new, 0.0, np.inf)
            C_new[~ocean_mask] = 0.0

            self.C_phyto_s[s, :, :] = C_new

        # Optional polar ring scalar averaging to avoid multi-longitude single-point inconsistency
        try:
            j_s, j_n = 0, -1
            ocean_row_s = ocean_mask[j_s, :]
            if np.any(ocean_row_s):
                for s in range(S):
                    row = self.C_phyto_s[s, j_s, :]
                    mean_s = float(np.mean(row[ocean_row_s]))
                    self.C_phyto_s[s, j_s, ocean_row_s] = mean_s
            ocean_row_n = ocean_mask[j_n, :]
            if np.any(ocean_row_n):
                for s in range(S):
                    row = self.C_phyto_s[s, j_n, :]
                    mean_n = float(np.mean(row[ocean_row_n]))
                    self.C_phyto_s[s, j_n, ocean_row_n] = mean_n
        except Exception:
            pass

    # ---------- Autosave I/O & Random init ----------

    def save_autosave(self, path_npz: str) -> bool:
        """
        Save minimal prognostic state for PhytoManager.

        Contents:
          - S (int), NL, NM
          - C_phyto_s [S,NL,NM]  (mg/m^3)
          - bands_lambda [NB]
          - params snapshot (subset)

        Returns True if saved.
        """
        try:
            os.makedirs(os.path.dirname(path_npz) or ".", exist_ok=True)
        except Exception:
            pass
        try:
            np.savez_compressed(
                path_npz,
                S=int(self.S),
                NL=int(self.NL),
                NM=int(self.NM),
                C_phyto_s=np.asarray(self.C_phyto_s, dtype=np.float32),
                bands_lambda=np.asarray(self.bands.lambda_centers, dtype=np.float32),
                H_mld=float(self.H_mld),
                mu_max_s=np.asarray(self.mu_max_s, dtype=np.float32),
                m0_s=np.asarray(self.m0_s, dtype=np.float32),
                init_frac_s=np.asarray(self.init_frac_s, dtype=np.float32),
            )
            if self.diag:
                print(f"[Phyto] Autosave written: '{path_npz}' (S={self.S}, NL={self.NL}, NM={self.NM})")
            return True
        except Exception as e:
            if self.diag:
                print(f"[Phyto] Autosave failed: {e}")
            return False

    def load_autosave(self, path_npz: str, *, on_mismatch: str = "random") -> bool:
        """
        Load minimal prognostic state. If shapes mismatch:
          - on_mismatch='random': randomize using current init fractions and chl0
          - on_mismatch='default': reset to defaults (uniform ocean, land=0)

        Returns True on success (or randomized/defaulted), False if hard failure.
        """
        try:
            data = np.load(path_npz)
        except Exception as e:
            if self.diag:
                print(f"[Phyto] Load autosave failed: {e}")
            return False

        try:
            S = int(data.get("S"))
            NL = int(data.get("NL"))
            NM = int(data.get("NM"))
            C_saved = np.asarray(data.get("C_phyto_s"))
        except Exception as e:
            if self.diag:
                print(f"[Phyto] Autosave malformed: {e}")
            return False

        if (S != self.S) or (NL != self.NL) or (NM != self.NM) or C_saved is None:
            if self.diag:
                print(f"[Phyto] Autosave shape mismatch (saved S={S},NL={NL},NM={NM}; current S={self.S},NL={self.NL},NM={self.NM}).")
            if on_mismatch == "random":
                self.randomize_state(seed=None)
                return True
            elif on_mismatch == "default":
                self.reset_default_state()
                return True
            return False

        try:
            self.C_phyto_s = np.clip(np.asarray(C_saved, dtype=float), 0.0, np.inf)
            # Respect ocean mask (land=0)
            for s in range(self.S):
                self.C_phyto_s[s, ~self.ocean_mask] = 0.0
            if self.diag:
                print(f"[Phyto] Autosave loaded: '{path_npz}' (S={S}, NL={NL}, NM={NM})")
            return True
        except Exception as e:
            if self.diag:
                print(f"[Phyto] Applying autosave failed: {e}")
            if on_mismatch == "random":
                self.randomize_state(seed=None)
                return True
            elif on_mismatch == "default":
                self.reset_default_state()
                return True
            return False

    def randomize_state(self, seed: int | None = None, noise_frac: float = 0.3) -> None:
        """
        Randomize phytoplankton state over ocean:
          C_s ≈ init_frac_s[s]*chl0 * (1 + noise in [-noise_frac, +noise_frac])
        Land cells remain 0.
        """
        rng = np.random.default_rng(seed)
        for s in range(self.S):
            base = self.init_frac_s[s] * self.params.chl0
            noise = (rng.random((self.NL, self.NM)) * 2.0 - 1.0) * noise_frac
            field = np.clip(base * (1.0 + noise), 0.0, np.inf)
            # apply ocean mask
            field[~self.ocean_mask] = 0.0
            self.C_phyto_s[s, :, :] = field
        if self.diag:
            print(f"[Phyto] State randomized (seed={seed}, noise_frac={noise_frac}).")

    def reset_default_state(self) -> None:
        """
        Reset to deterministic default initial condition:
          C_s = init_frac_s[s]*chl0 over ocean; 0 over land.
        """
        for s in range(self.S):
            field = np.full((self.NL, self.NM), self.init_frac_s[s] * self.params.chl0, dtype=float)
            field[~self.ocean_mask] = 0.0
            self.C_phyto_s[s, :, :] = field
        if self.diag:
            print("[Phyto] State reset to defaults.")

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
        S_default = 10
        try:
            S_default = int(os.getenv("QD_PHYTO_NSPECIES", "10"))
        except Exception:
            S_default = 10
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

        # --- Nutrient competition (optional, single N pool) ---
        self.enable_N = int(os.getenv("QD_PHYTO_ENABLE_N", "1")) == 1
        # Per-species half-saturation (mmol m^-3) and yield (mg Chl per mmol N)
        KN_list = _read_env_list("QD_PHYTO_KN") or []
        Y_list  = _read_env_list("QD_PHYTO_YIELD") or []
        self.KN_s = np.array(
            [(KN_list[s] if s < len(KN_list) else 0.5) for s in range(self.S)],
            dtype=float
        )
        self.Y_s = np.array(
            [(Y_list[s] if s < len(Y_list) else 1.0) for s in range(self.S)],
            dtype=float
        )
        # Remineralization source (mmol m^-3 d^-1) and initial N
        self.R_remin = _read_env_float("QD_PHYTO_REMIN", 0.01)
        N_init = _read_env_float("QD_PHYTO_N_INIT", 1.0)
        self.N = np.full((self.NL, self.NM), N_init, dtype=float)
        self.N[~self.ocean_mask] = 0.0

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

        # 3) Species‑specific light proxy via band integration:
        #    E_s = Σ_b Ī_b · Shape_s[b] · Δλ_b (if available)  → drives spectral light limitation
        try:
            dlam = np.asarray(self.bands.delta_lambda, dtype=float)
        except Exception:
            dlam = None
        if dlam is not None and dlam.size == self.bands.nbands:
            E_s = np.tensordot(self.shape_sb, Ibar_b * dlam[:, None, None], axes=(1, 0))  # [S, NL, NM]
        else:
            E_s = np.tensordot(self.shape_sb, Ibar_b, axes=(1, 0))  # [S, NL, NM]

        # 4) Growth modifiers: species light limitation (tanh) and shared temperature factor (Q10)
        #    Note: Divide by μ_max_s per species to keep tanh argument scale consistent.
        muL_s = np.tanh(self.params.alpha_P * E_s / np.maximum(self.mu_max_s[:, None, None], 1e-6))
        fT = np.power(self.params.Q10, (np.asarray(T_w, dtype=float) - self.params.T_ref) / 10.0)

        # 5) Net specific rate per SPECIES (d^-1): μ_s = μ_max_s · muL_s · fT − (m0_s + sink/H)
        sink_term = 0.0
        if self.params.lambda_sink_m_per_day > 0.0:
            sink_term = float(self.params.lambda_sink_m_per_day) / max(1e-6, self.H_mld)

        # Nutrient limitation and growth split
        if self.enable_N:
            KN = np.maximum(self.KN_s[:, None, None], 1e-12)  # [S,1,1]
            Nmap = np.asarray(self.N, dtype=float)[None, :, :]  # [1,NL,NM]
            fN_s = Nmap / (KN + Nmap)  # [S,NL,NM] in [0,1]
            mu_grow_s = self.mu_max_s[:, None, None] * muL_s * fT[None, :, :] * np.clip(fN_s, 0.0, 1.0)
        else:
            mu_grow_s = self.mu_max_s[:, None, None] * muL_s * fT[None, :, :]

        mu_s = mu_grow_s - (self.m0_s[:, None, None] + sink_term)

        # 6) Update chlorophyll per species
        dC_s = mu_s * self.C_phyto_s * float(dt_days)
        self.C_phyto_s = np.clip(self.C_phyto_s + dC_s, 0.0, np.inf)
        # Keep land cells at zero
        for s in range(self.S):
            self.C_phyto_s[s, ~self.ocean_mask] = 0.0

        # Nutrient update (mmol m^-3)
        if self.enable_N:
            Y = np.maximum(self.Y_s[:, None, None], 1e-12)  # [S,1,1]
            # Uptake = μ_grow · C / Y  (mmol m^-3 d^-1)
            uptake = (mu_grow_s * self.C_phyto_s) / Y
            total_uptake = np.sum(uptake, axis=0)  # [NL,NM]
            dN = (- total_uptake + float(self.R_remin)) * float(dt_days)
            self.N = np.clip(np.asarray(self.N, dtype=float) + dN, 0.0, np.inf)
            self.N[~self.ocean_mask] = 0.0

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
                N=np.asarray(self.N, dtype=np.float32),
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
            # Optional: load nutrient pool if present
            try:
                N_saved = data.get("N")
                if N_saved is not None:
                    N_arr = np.asarray(N_saved, dtype=float)
                    if N_arr.shape == (self.NL, self.NM):
                        self.N = np.clip(N_arr, 0.0, np.inf)
                        self.N[~self.ocean_mask] = 0.0
            except Exception:
                pass
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

    # --- Standardized IO: plankton.json (bio/optics) and plankton.nc (distributions) ---

    def save_bio_json(self, path: str, day_value: float | None = None) -> bool:
        """
        Save phytoplankton biology & optics configuration to a human-readable JSON:
          - schema_version, source, day
          - bands: nbands, lambda_centers_nm[], delta_lambda_nm[]
          - params (shared): alpha_P, Q10, T_ref, lambda_sink_m_per_day
          - species arrays: mu_max_s[], m0_s[], c_reflect_s[], p_reflect_s[], shape_sb[S][NB]
          - kd optics: Kd0_b[NB], kchl_b[NB], Apure_b[NB]
        This mirrors 'genes.json' style for terrestrial ecology.
        """
        try:
            import json as _json
            os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
            doc = {
                "schema_version": 1,
                "source": "PyGCM.PhytoManager.save_bio_json",
                "day": float(day_value) if day_value is not None else None,
                "bands": {
                    "nbands": int(self.bands.nbands),
                    "lambda_centers_nm": [float(x) for x in np.asarray(self.bands.lambda_centers, dtype=float).tolist()],
                    "delta_lambda_nm": [float(x) for x in np.asarray(self.bands.delta_lambda, dtype=float).tolist()],
                },
                "params": {
                    "alpha_P": float(self.params.alpha_P),
                    "Q10": float(self.params.Q10),
                    "T_ref": float(self.params.T_ref),
                    "lambda_sink_m_per_day": float(self.params.lambda_sink_m_per_day),
                },
                "species": {
                    "mu_max_s": [float(x) for x in np.asarray(self.mu_max_s, dtype=float).tolist()],
                    "m0_s": [float(x) for x in np.asarray(self.m0_s, dtype=float).tolist()],
                    "c_reflect_s": [float(x) for x in np.asarray(self.c_reflect_s, dtype=float).tolist()],
                    "p_reflect_s": [float(x) for x in np.asarray(self.p_reflect_s, dtype=float).tolist()],
                    # Store per-species spectral shapes (Gaussian weights per band, normalized)
                    "shape_sb": np.asarray(self.shape_sb, dtype=float).tolist(),
                },
                "optics": {
                    "Kd0_b": [float(x) for x in np.asarray(self.Kd0_b, dtype=float).tolist()],
                    "kchl_b": [float(x) for x in np.asarray(self.kchl_b, dtype=float).tolist()],
                    "Apure_b": [float(x) for x in np.asarray(self.Apure_b, dtype=float).tolist()],
                },
            }
            with open(path, "w", encoding="utf-8") as f:
                _json.dump(doc, f, ensure_ascii=False, indent=2)
            if self.diag:
                print(f"[Phyto] Bio/optics JSON written: '{path}' (S={self.S}, NB={self.bands.nbands})")
            return True
        except Exception as e:
            if self.diag:
                print(f"[Phyto] save_bio_json failed: {e}")
            return False

    def save_distribution_nc(self, path: str, day_value: float | None = None) -> bool:
        """
        Save gridded distributions to NetCDF 'plankton.nc':
          - dims: species, band, lat, lon
          - lat, lon coordinates
          - C_phyto_s[S,lat,lon] (mg/m^3)
          - alpha_water_bands[band,lat,lon] (if available)
          - alpha_water_scalar[lat,lon]
          - Kd_490[lat,lon]
          - attributes: H_mld, S, NB, day
        """
        try:
            from netCDF4 import Dataset
            os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
            with Dataset(path, "w") as ds:
                NL, NM = self.NL, self.NM
                ds.createDimension("lat", NL)
                ds.createDimension("lon", NM)
                ds.createDimension("species", int(self.S))
                ds.createDimension("band", int(self.bands.nbands))

                vlat = ds.createVariable("lat", "f4", ("lat",))
                vlon = ds.createVariable("lon", "f4", ("lon",))
                vlat[:] = self.grid.lat
                vlon[:] = self.grid.lon

                # Chlorophyll distributions
                vC = ds.createVariable("C_phyto_s", "f4", ("species", "lat", "lon"))
                vC[:] = np.asarray(self.C_phyto_s, dtype=np.float32)

                # Albedo maps
                if self.alpha_water_bands is not None:
                    vab = ds.createVariable("alpha_water_bands", "f4", ("band", "lat", "lon"))
                    vab[:] = np.asarray(self.alpha_water_bands, dtype=np.float32)
                vas = ds.createVariable("alpha_water_scalar", "f4", ("lat", "lon"))
                vas[:] = np.asarray(self.alpha_water_scalar, dtype=np.float32)

                # Kd(490)
                vk = ds.createVariable("Kd_490", "f4", ("lat", "lon"))
                vk[:] = np.asarray(self.Kd_490, dtype=np.float32)

                # Nutrient pool (optional)
                try:
                    vN = ds.createVariable("N", "f4", ("lat", "lon"))
                    vN[:] = np.asarray(self.N, dtype=np.float32)
                except Exception:
                    pass

                # Band centers (for reference)
                vb = ds.createVariable("bands_lambda_centers", "f4", ("band",))
                vb[:] = np.asarray(self.bands.lambda_centers, dtype=np.float32)

                # Attributes
                ds.setncattr("title", "Qingdai Phytoplankton Distributions")
                ds.setncattr("H_mld_m", float(self.H_mld))
                ds.setncattr("S", int(self.S))
                ds.setncattr("NB", int(self.bands.nbands))
                if day_value is not None:
                    ds.setncattr("day", float(day_value))
            if self.diag:
                print(f"[Phyto] Distribution NetCDF written: '{path}'")
            return True
        except Exception as e:
            if self.diag:
                print(f"[Phyto] save_distribution_nc failed: {e}")
            return False

    def load_bio_json(self, path: str, *, on_mismatch: str = "keep") -> bool:
        """
        Load phytoplankton biology & optics configuration from plankton.json.
        on_mismatch:
          - 'keep': keep current bands if JSON bands mismatch; rebuild shapes to current NB
          - 'replace': attempt to replace bands with JSON bands (if compatible with the rest of model)
        Returns True on successful application.
        """
        try:
            import json as _json
            with open(path, "r", encoding="utf-8") as f:
                doc = _json.load(f)
        except Exception as e:
            if self.diag:
                print(f"[Phyto] load_bio_json failed: {e}")
            return False

        try:
            bands_json = doc.get("bands", {}) or {}
            nb_json = int(bands_json.get("nbands", self.bands.nbands))
            lam_cent = np.asarray(bands_json.get("lambda_centers_nm", []), dtype=float)
            dlam = np.asarray(bands_json.get("delta_lambda_nm", []), dtype=float)

            # Decide bands handling
            replace_bands = (on_mismatch == "replace" and nb_json > 0 and lam_cent.size == nb_json)
            if replace_bands:
                # Rebuild SpectralBands using JSON band definition
                try:
                    from .spectral import SpectralBands
                    lam_edges = None
                    if dlam.size == nb_json:
                        # Derive edges from centers and delta if possible (approx)
                        half = 0.5 * dlam
                        lam_edges = np.zeros((nb_json + 1,), dtype=float)
                        lam_edges[0] = lam_cent[0] - half[0]
                        for i in range(1, nb_json):
                            lam_edges[i] = 0.5 * (lam_cent[i-1] + lam_cent[i])
                        lam_edges[-1] = lam_cent[-1] + half[-1]
                    self.bands = SpectralBands(
                        nbands=nb_json,
                        lambda_centers=lam_cent if lam_cent.size == nb_json else None,
                        delta_lambda=dlam if dlam.size == nb_json else None,
                        lambda_edges=lam_edges
                    )
                    if self.diag:
                        print(f"[Phyto] Bands replaced from JSON: NB={self.bands.nbands}")
                except Exception as _be:
                    if self.diag:
                        print(f"[Phyto] bands replace failed, keeping current: {_be}")

            # Update shared params
            p = doc.get("params", {}) or {}
            self.params.alpha_P = float(p.get("alpha_P", self.params.alpha_P))
            self.params.Q10 = float(p.get("Q10", self.params.Q10))
            self.params.T_ref = float(p.get("T_ref", self.params.T_ref))
            self.params.lambda_sink_m_per_day = float(p.get("lambda_sink_m_per_day", self.params.lambda_sink_m_per_day))

            # Update species arrays
            sp = doc.get("species", {}) or {}
            mu_max_s = np.asarray(sp.get("mu_max_s", []), dtype=float)
            m0_s = np.asarray(sp.get("m0_s", []), dtype=float)
            c_reflect_s = np.asarray(sp.get("c_reflect_s", []), dtype=float)
            p_reflect_s = np.asarray(sp.get("p_reflect_s", []), dtype=float)
            shape_sb = np.asarray(sp.get("shape_sb", []), dtype=float)

            # Resize S if JSON provides arrays of a different length
            def _apply_arr(name, arr_json, fallback):
                nonlocal_changed = False
                if arr_json is not None and arr_json.size > 0:
                    return arr_json.astype(float), True
                return fallback, False

            changed_any = False
            if mu_max_s.size > 0:
                self.mu_max_s = mu_max_s.astype(float); changed_any = True
            if m0_s.size > 0:
                self.m0_s = m0_s.astype(float); changed_any = True
            if c_reflect_s.size > 0:
                self.c_reflect_s = c_reflect_s.astype(float); changed_any = True
            if p_reflect_s.size > 0:
                self.p_reflect_s = p_reflect_s.astype(float); changed_any = True
            # Update S from arrays length if consistent
            S_new = max(self.S,
                        self.mu_max_s.size,
                        self.m0_s.size,
                        self.c_reflect_s.size,
                        self.p_reflect_s.size)
            if S_new != self.S:
                self.S = int(S_new)
                # Ensure arrays length match S
                def _ensure_len(a, val=1.0):
                    if a.size == self.S:
                        return a
                    if a.size == 0:
                        return np.full((self.S,), val, dtype=float)
                    if a.size < self.S:
                        pad = np.full((self.S - a.size,), a[-1] if a.size > 0 else val, dtype=float)
                        return np.concatenate([a, pad], axis=0)
                    return a[:self.S]
                self.mu_max_s = _ensure_len(self.mu_max_s, self.params.mu_max)
                self.m0_s = _ensure_len(self.m0_s, self.params.m0)
                self.c_reflect_s = _ensure_len(self.c_reflect_s, 0.02)
                self.p_reflect_s = _ensure_len(self.p_reflect_s, 0.5)
                # Resize shapes to [S, NB]
                NB = self.bands.nbands
                shp = np.zeros((self.S, NB), dtype=float)
                if shape_sb.ndim == 2:
                    for s in range(min(self.S, shape_sb.shape[0])):
                        v = shape_sb[s, :]
                        if v.size == NB:
                            shp[s, :] = v / (float(np.sum(v)) + 1e-12)
                self.shape_sb = np.where(shp > 0, shp, self.shape_sb if getattr(self, "shape_sb", None) is not None else shp)

            # Optics Kd and Apure per band
            opt = doc.get("optics", {}) or {}
            Kd0_b = np.asarray(opt.get("Kd0_b", []), dtype=float)
            kchl_b = np.asarray(opt.get("kchl_b", []), dtype=float)
            Apure_b = np.asarray(opt.get("Apure_b", []), dtype=float)
            NB = self.bands.nbands
            if Kd0_b.size == NB:
                self.Kd0_b = Kd0_b.astype(float)
            if kchl_b.size == NB:
                self.kchl_b = kchl_b.astype(float)
            if Apure_b.size == NB:
                self.Apure_b = Apure_b.astype(float)

            # Recompute band weights with current bands
            from .spectral import band_weights_from_mode
            self.w_b = band_weights_from_mode(self.bands)

            if self.diag:
                print(f"[Phyto] Bio/optics JSON loaded: S={self.S}, NB={self.bands.nbands} (bands {'replaced' if replace_bands else 'kept'})")
            return True
        except Exception as e:
            if self.diag:
                print(f"[Phyto] load_bio_json apply failed: {e}")
            return False

    def load_distribution_nc(self, path: str, *, on_mismatch: str = "keep") -> bool:
        """
        Load gridded distributions from plankton.nc into current manager.
        Applies only if shapes match (S, NB, NL, NM). Returns True on success.
        on_mismatch: 'keep' → silently keep current state if dims mismatch; 'reset' → reset defaults.
        """
        try:
            from netCDF4 import Dataset
            with Dataset(path, "r") as ds:
                C = np.asarray(ds.variables["C_phyto_s"]) if "C_phyto_s" in ds.variables else None
                ab = np.asarray(ds.variables["alpha_water_bands"]) if "alpha_water_bands" in ds.variables else None
                aS = np.asarray(ds.variables["alpha_water_scalar"]) if "alpha_water_scalar" in ds.variables else None
                kd = np.asarray(ds.variables["Kd_490"]) if "Kd_490" in ds.variables else None
                # Optional: bands check
                try:
                    lam_nc = np.asarray(ds.variables["bands_lambda_centers"])
                    if lam_nc.size == self.bands.nbands:
                        pass  # OK; otherwise ignore silently (bands might differ)
                except Exception:
                    pass
        except Exception as e:
            if self.diag:
                print(f"[Phyto] load_distribution_nc failed: {e}")
            return False

        # Validate shapes
        try:
            NL, NM = self.NL, self.NM
            NB = self.bands.nbands
            okC = (C is not None and C.ndim == 3 and C.shape[1] == NL and C.shape[2] == NM)
            okAb = (ab is None) or (ab.ndim == 3 and ab.shape[1] == NL and ab.shape[2] == NM and ab.shape[0] == NB)
            okAS = (aS is None) or (aS.shape == (NL, NM))
            okKd = (kd is None) or (kd.shape == (NL, NM))
            if not okC or not okAb or not okAS or not okKd:
                if self.diag:
                    print(f"[Phyto] plankton.nc dims mismatch; keep={on_mismatch=='keep'}")
                if on_mismatch == "reset":
                    self.reset_default_state()
                return False

            # Apply loaded fields
            self.C_phyto_s = np.clip(C.astype(float), 0.0, np.inf)
            # enforce ocean mask
            for s in range(self.S):
                self.C_phyto_s[s, ~self.ocean_mask] = 0.0
            if ab is not None:
                self.alpha_water_bands = np.clip(ab.astype(float), self.alpha_clip_min, self.alpha_clip_max)
            if aS is not None:
                self.alpha_water_scalar = np.clip(aS.astype(float), self.alpha_clip_min, self.alpha_clip_max)
            if kd is not None:
                self.Kd_490 = np.clip(kd.astype(float), 0.0, np.inf)
            if self.diag:
                print(f"[Phyto] plankton.nc loaded: C_phyto_s[{self.C_phyto_s.shape}], "
                      f"alpha_bands={'OK' if ab is not None else 'none'}, alpha_scalar={'OK' if aS is not None else 'none'}")
            return True
        except Exception as e:
            if self.diag:
                print(f"[Phyto] apply distribution failed: {e}")
            return False

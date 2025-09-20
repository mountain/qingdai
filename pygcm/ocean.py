# pygcm/ocean.py
"""
Wind-driven single-layer (slab) ocean on a spherical grid.

Implements a minimal shallow-water-like barotropic ocean with:
- Prognostic currents (uo, vo) and sea surface height anomaly (eta)
- Wind stress forcing from near-surface atmospheric wind
- Linear bottom drag
- Scale-selective dissipation (∇⁴ hyperdiffusion) + optional Shapiro smoothing
- SST (surface temperature) advection by ocean currents
- Optional vertical heating by Q_net/(rho_w * c_p,w * H) (disabled by default to avoid
  double counting with P006 energy integration; can be enabled via QD_OCEAN_USE_QNET=1)

This module targets Projects/011 M1, M2, M3 (minimal viable dynamics + coupling).
"""

from __future__ import annotations

import os
import numpy as np

from .grid import SphericalGrid
from . import constants as const


class WindDrivenSlabOcean:
    def __init__(self,
                 grid: SphericalGrid,
                 land_mask: np.ndarray,
                 H_m: float,
                 init_Ts: np.ndarray | None = None,
                 rho_w: float | None = None,
                 cp_w: float | None = None):
        """
        Args:
            grid: global spherical grid (same as atmosphere to avoid regridding)
            land_mask: 1 for land, 0 for ocean
            H_m: mixed layer depth (m)
            init_Ts: initial SST (K); if None, uniform 288K
            rho_w: seawater density, default from env QD_RHO_W or 1000
            cp_w: seawater specific heat, default from env QD_CP_W or 4200
        """
        self.grid = grid
        self.land_mask = np.asarray(land_mask, dtype=int)
        self.H = float(H_m)

        # Physical parameters
        self.rho_w = float(os.getenv("QD_RHO_W", str(rho_w if rho_w is not None else 1000.0)))
        self.cp_w = float(os.getenv("QD_CP_W", str(cp_w if cp_w is not None else 4200.0)))
        self.g = 9.81

        # Drag / wind stress
        self.CD = float(os.getenv("QD_CD", "1.5e-3"))          # air-ocean drag coefficient
        self.r_bot = float(os.getenv("QD_R_BOT", "1.0e-6"))     # bottom drag (s^-1)
        self.rho_a = float(os.getenv("QD_RHO_A", "1.2"))        # air density for wind stress

        # Mixing / dissipation
        self.K_h = float(os.getenv("QD_KH_OCEAN", "5.0e3"))     # lateral mixing of SST (m^2/s)
        self.sigma4 = float(os.getenv("QD_SIGMA4_OCEAN", "0.02"))
        self.k4_nsub = int(os.getenv("QD_OCEAN_K4_NSUB", "1"))
        self.diff_every = int(os.getenv("QD_OCEAN_DIFF_EVERY", "1"))
        self.shapiro_n = int(os.getenv("QD_OCEAN_SHAPIRO_N", "0"))           # 0 disables
        self.shapiro_every = int(os.getenv("QD_OCEAN_SHAPIRO_EVERY", "8"))   # default cadence

        # CFL guard (diagnostic)
        self.cfl_target = float(os.getenv("QD_OCEAN_CFL", "0.5"))
        self.max_u_cap = float(os.getenv("QD_OCEAN_MAX_U", "3.0e2"))  # sanity cap

        # Grid metrics
        self.a = const.PLANET_RADIUS
        self.dlat = self.grid.dlat_rad
        self.dlon = self.grid.dlon_rad
        self.lat_rad = np.deg2rad(self.grid.lat_mesh)
        self.coslat = np.maximum(np.cos(self.lat_rad), 0.2)  # guard near poles (strong cap to avoid metric blow-up)
        self.f = self.grid.coriolis_param

        # Prognostic fields
        self.uo = np.zeros_like(self.grid.lat_mesh, dtype=float)
        self.vo = np.zeros_like(self.grid.lat_mesh, dtype=float)
        self.eta = np.zeros_like(self.grid.lat_mesh, dtype=float)

        # SST
        if init_Ts is None:
            self.Ts = np.full_like(self.grid.lat_mesh, 288.0, dtype=float)
        else:
            self.Ts = np.array(init_Ts, dtype=float, copy=True)

        # Internal counter for cadence controls
        self._step = 0

    # ----------------- Numerical utilities -----------------
    def _laplacian_sphere(self, F: np.ndarray) -> np.ndarray:
        """
        ∇²F on a regular lat-lon grid using divergence form with cosφ metric.
        """
        F = np.nan_to_num(F, copy=False)
        dF_dphi = np.gradient(F, self.dlat, axis=0)
        term_phi = (1.0 / self.coslat) * np.gradient(self.coslat * dF_dphi, self.dlat, axis=0)
        d2F_dlam2 = (np.roll(F, -1, axis=1) - 2.0 * F + np.roll(F, 1, axis=1)) / (self.dlon ** 2)
        term_lam = d2F_dlam2 / (self.coslat ** 2)
        return (term_phi + term_lam) / (self.a ** 2)

    def _hyperdiffuse(self, F: np.ndarray, dt: float, k4: float, n_substeps: int = 1) -> np.ndarray:
        """
        Explicit ∇⁴ hyperdiffusion with optional substeps for stability.
        """
        if k4 <= 0.0 or dt <= 0.0:
            return F
        n = max(1, int(n_substeps))
        sub_dt = dt / n
        out = np.nan_to_num(F, copy=True)
        for _ in range(n):
            L = self._laplacian_sphere(out)
            L2 = self._laplacian_sphere(L)
            out = out - k4 * L2 * sub_dt
        return out

    def _shapiro_filter(self, F: np.ndarray, n: int = 2) -> np.ndarray:
        """
        Separable 1-2-1 smoothing applied n times.
        """
        from scipy.ndimage import convolve
        k1 = np.array([1.0, 2.0, 1.0], dtype=float); k1 /= k1.sum()
        out = np.nan_to_num(F, copy=True)
        for _ in range(max(1, int(n))):
            out = convolve(out, k1[np.newaxis, :], mode="wrap")
            out = convolve(out, k1[:, np.newaxis], mode="nearest")
        return out

    def _advect_scalar(self, field: np.ndarray, u: np.ndarray, v: np.ndarray, dt: float) -> np.ndarray:
        """
        Semi-Lagrangian advection (bilinear interpolation, lon periodic).
        """
        from scipy.ndimage import map_coordinates

        # Convert velocities to index displacements
        dlam = u * dt / (self.a * self.coslat)   # radians of longitude
        dphi = v * dt / self.a                   # radians of latitude

        dx = dlam / self.dlon
        dy = dphi / self.dlat

        lats = np.arange(self.grid.n_lat)
        lons = np.arange(self.grid.n_lon)
        JJ, II = np.meshgrid(lats, lons, indexing="ij")

        dep_J = JJ - dy
        dep_I = II - dx

        adv = map_coordinates(field, [dep_J, dep_I], order=1, mode="wrap", prefilter=False)
        return adv

    # ----------------- Main step -----------------
    def step(self,
             dt: float,
             u_atm: np.ndarray,
             v_atm: np.ndarray,
             Q_net: np.ndarray | None = None,
             ice_mask: np.ndarray | None = None) -> None:
        """
        Advance ocean state by one step.

        Args:
            dt: time step (s)
            u_atm, v_atm: near-surface atmospheric wind (m/s) to derive wind stress
            Q_net: optional net surface heat flux into ocean (W/m^2). If provided and
                   QD_OCEAN_USE_QNET=1, used to heat/cool SST.
            ice_mask: optional boolean mask where sea-ice present; used to suppress SST update
        """
        self._step += 1

        # Precompute wind stress for this step (kept constant within substeps)
        Va = np.sqrt(u_atm**2 + v_atm**2)
        tau_x = self.rho_a * self.CD * Va * u_atm
        tau_y = self.rho_a * self.CD * Va * v_atm

        # Determine stable substeps based on gravity wave and advective CFL
        dx_lat = self.a * self.dlat
        min_cos = float(np.min(self.coslat))
        dx_lon_min = self.a * self.dlon * max(1e-3, min_cos)
        dx_min = min(dx_lat, dx_lon_min)
        c = np.sqrt(self.g * self.H)
        uadv = float(np.max(np.sqrt(self.uo**2 + self.vo**2)))
        uadv = max(uadv, float(np.max(Va)))
        target = max(1e-3, self.cfl_target)
        n_sub = int(np.ceil(max(c, uadv) * (dt / max(1e-12, dx_min)) / target))
        n_sub = int(max(1, min(500, n_sub)))
        sub_dt = dt / n_sub

        for _sub in range(n_sub):
            # 2) Pressure gradient from current eta
            deta_dlam = (np.roll(self.eta, -1, axis=1) - np.roll(self.eta, 1, axis=1)) / (2.0 * self.dlon)
            deta_dphi = (np.roll(self.eta, -1, axis=0) - np.roll(self.eta, 1, axis=0)) / (2.0 * self.dlat)
            grad_eta_x = deta_dlam / (self.a * self.coslat)   # ∂η/∂x
            grad_eta_y = deta_dphi / self.a                   # ∂η/∂y

            # 3) Momentum
            du = ( self.f * self.vo
                   - self.g * grad_eta_x
                   + tau_x / (self.rho_w * self.H)
                   - self.r_bot * self.uo )
            dv = ( - self.f * self.uo
                   - self.g * grad_eta_y
                   + tau_y / (self.rho_w * self.H)
                   - self.r_bot * self.vo )

            self.uo += sub_dt * du
            self.vo += sub_dt * dv

            # Land boundary
            on_land = (self.land_mask == 1)
            self.uo[on_land] = 0.0
            self.vo[on_land] = 0.0

            # 4) Scale-selective dissipation (cadence tied to outer step)
            if (self.diff_every > 0) and (self._step % self.diff_every == 0):
                k4_base = self.sigma4 * (dx_min**4) / max(1e-12, sub_dt)
                self.uo = self._hyperdiffuse(self.uo, sub_dt, k4_base, n_substeps=self.k4_nsub)
                self.vo = self._hyperdiffuse(self.vo, sub_dt, k4_base, n_substeps=self.k4_nsub)
                self.eta = self._hyperdiffuse(self.eta, sub_dt, 0.5 * k4_base, n_substeps=self.k4_nsub)

            # Optional Shapiro smoothing (off by default)
            if (self.shapiro_n > 0) and (self.shapiro_every > 0) and (self._step % self.shapiro_every == 0):
                self.uo = self._shapiro_filter(self.uo, n=self.shapiro_n)
                self.vo = self._shapiro_filter(self.vo, n=self.shapiro_n)
                self.eta = self._shapiro_filter(self.eta, n=self.shapiro_n)

            # 5) Continuity
            div = self.grid.divergence(self.uo, self.vo)
            self.eta += - sub_dt * self.H * div
            self.eta[on_land] = 0.0

            # 6) SST advection by ocean currents (semi-Lagrangian, gentle blend)
            adv_alpha = float(os.getenv("QD_OCEAN_ADV_ALPHA", "0.7"))  # 0..1
            Ts_adv = self._advect_scalar(self.Ts, self.uo, self.vo, sub_dt)
            self.Ts = (1.0 - adv_alpha) * self.Ts + adv_alpha * Ts_adv

            # Lateral diffusion of SST (explicit)
            if self.K_h > 0.0:
                self.Ts += sub_dt * self.K_h * self._laplacian_sphere(self.Ts)

            # 7) Optional vertical heat flux (disabled by default to avoid double counting)
            use_qnet = int(os.getenv("QD_OCEAN_USE_QNET", "0")) == 1
            if use_qnet and (Q_net is not None):
                heat_tendency = Q_net / (self.rho_w * self.cp_w * self.H)  # K/s
                if ice_mask is not None:
                    mask = (self.land_mask == 0) & (~ice_mask)
                else:
                    mask = (self.land_mask == 0)
                self.Ts = np.where(mask, self.Ts + sub_dt * heat_tendency, self.Ts)

            # 8) Sanity caps per substep
            self.uo = np.clip(np.nan_to_num(self.uo), -self.max_u_cap, self.max_u_cap)
            self.vo = np.clip(np.nan_to_num(self.vo), -self.max_u_cap, self.max_u_cap)
            self.eta = np.nan_to_num(self.eta)
            # Clamp eta to a safe anomaly range (meters) to prevent runaway
            _eta_cap = float(os.getenv("QD_ETA_CAP", "200.0"))
            self.eta = np.clip(self.eta, -_eta_cap, _eta_cap)
            self.Ts = np.nan_to_num(self.Ts)

        # Final clamp on Ts to safe physical bounds
        ts_min = float(os.getenv("QD_TS_MIN", "150.0"))
        ts_max = float(os.getenv("QD_TS_MAX", "340.0"))
        self.Ts = np.clip(self.Ts, ts_min, ts_max)

    def diagnostics(self) -> dict:
        """
        Return minimal diagnostics: KE, max|u|, eta stats, CFL indicator.
        """
        w = np.maximum(np.cos(self.lat_rad), 0.0)
        wsum = np.sum(w) + 1e-15
        KE = 0.5 * (self.uo**2 + self.vo**2)
        KE_mean = float(np.sum(KE * w) / wsum)
        Umax = float(np.max(np.sqrt(self.uo**2 + self.vo**2)))
        eta_min = float(np.min(self.eta))
        eta_max = float(np.max(self.eta))

        # Estimate most conservative gravity wave speed sqrt(gH) CFL
        c = np.sqrt(self.g * self.H)
        dx_lat = self.a * self.dlat
        min_cos = float(np.min(self.coslat))
        dx_lon_min = self.a * self.dlon * max(1e-3, min_cos)
        dx_min = min(dx_lat, dx_lon_min)
        cfl = float(c / max(1e-12, dx_min))  # per second

        return {
            "KE_mean": KE_mean,
            "U_max": Umax,
            "eta_min": eta_min,
            "eta_max": eta_max,
            "cfl_per_s": cfl
        }

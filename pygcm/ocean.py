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
from .jax_compat import is_enabled as _jax_enabled, to_numpy as _to_np, laplacian_sphere as _j_lap, hyperdiffuse as _j_hyp, advect_semilag as _j_adv


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
        self.r_bot = float(os.getenv("QD_R_BOT", "2.0e-5"))     # bottom drag (s^-1) stronger default for slab ocean
        self.rho_a = float(os.getenv("QD_RHO_A", "1.2"))        # air density for wind stress
        self.vcap = float(os.getenv("QD_WIND_STRESS_VCAP", "15.0"))  # cap for wind speed in stress (m/s)
        self.tau_scale = float(os.getenv("QD_TAU_SCALE", "0.2"))     # momentum transfer efficiency to slab
        # Polar sponge (extra drag near poles)
        self.polar_lat0 = float(os.getenv("QD_POLAR_SPONGE_LAT", "70.0"))   # deg
        self.polar_gain = float(os.getenv("QD_POLAR_SPONGE_GAIN", "5.0e-5")) # s^-1 at pole

        # Mixing / dissipation
        self.K_h = float(os.getenv("QD_KH_OCEAN", "5.0e3"))     # lateral mixing of SST (m^2/s)
        self.sigma4 = float(os.getenv("QD_SIGMA4_OCEAN", "0.02"))
        self.k4_nsub = int(os.getenv("QD_OCEAN_K4_NSUB", "1"))
        self.diff_every = int(os.getenv("QD_OCEAN_DIFF_EVERY", "1"))
        self.shapiro_n = int(os.getenv("QD_OCEAN_SHAPIRO_N", "0"))           # 0 disables
        self.shapiro_every = int(os.getenv("QD_OCEAN_SHAPIRO_EVERY", "8"))   # default cadence

        # CFL guard (diagnostic)
        self.cfl_target = float(os.getenv("QD_OCEAN_CFL", "0.5"))
        self.max_u_cap = float(os.getenv("QD_OCEAN_MAX_U", "3.0"))  # sanity cap (ocean currents, vector cap)
        # Outlier handling for unrealistically fast currents: "mean4" (default) or "clamp"
        self.outlier_method = os.getenv("QD_OCEAN_OUTLIER", "mean4").strip().lower()

        # Grid metrics
        self.a = const.PLANET_RADIUS
        self.dlat = self.grid.dlat_rad
        self.dlon = self.grid.dlon_rad
        self.lat_rad = np.deg2rad(self.grid.lat_mesh)
        self.coslat = np.maximum(np.cos(self.lat_rad), 0.5)  # stronger cap to avoid polar metric blow-up
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
        # Try JAX-jitted kernel
        try:
            if _jax_enabled():
                return _to_np(_j_lap(F, self.dlat, self.dlon, self.coslat, self.a))
        except Exception:
            pass

        # NumPy fallback
        F = np.nan_to_num(F, copy=False)
        dF_dphi = np.gradient(F, self.dlat, axis=0)
        term_phi = (1.0 / self.coslat) * np.gradient(self.coslat * dF_dphi, self.dlat, axis=0)
        d2F_dlam2 = (np.roll(F, -1, axis=1) - 2.0 * F + np.roll(F, 1, axis=1)) / (self.dlon ** 2)
        term_lam = d2F_dlam2 / (self.coslat ** 2)
        return (term_phi + term_lam) / (self.a ** 2)

    def _hyperdiffuse(self, F: np.ndarray, dt: float, k4, n_substeps: int = 1) -> np.ndarray:
        """
        Explicit ∇⁴ hyperdiffusion with optional substeps for stability.
        k4 can be a scalar (float) or a 2D map broadcastable to F.
        """
        # Try JAX-jitted kernel
        try:
            if _jax_enabled():
                return _to_np(_j_hyp(F, k4, dt, n_substeps, self.dlat, self.dlon, self.coslat, self.a))
        except Exception:
            pass

        if dt <= 0.0:
            return F
        # NumPy fallback
        try:
            if np.isscalar(k4):
                k4_arr = float(k4)
                if k4_arr <= 0.0:
                    return F
            else:
                k4_arr = np.nan_to_num(k4, copy=False)
                if np.all(k4_arr <= 0.0):
                    return F
        except Exception:
            return F
        n = max(1, int(n_substeps))
        sub_dt = dt / n
        out = np.nan_to_num(F, copy=True)
        for _ in range(n):
            L = self._laplacian_sphere(out)
            L2 = self._laplacian_sphere(L)
            out = out - k4_arr * L2 * sub_dt
        return np.nan_to_num(out, copy=False)

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
        # Try JAX path
        try:
            if _jax_enabled():
                adv = _j_adv(field, u, v, dt, self.a, self.dlat, self.dlon, self.coslat)
                return _to_np(adv)
        except Exception:
            pass

        # Fallback SciPy path
        from scipy.ndimage import map_coordinates
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

    # ----------------- Polar corrections (poles averaging) -----------------
    def _polar_scalar_average_fill(self, F: np.ndarray, ocean_mask: np.ndarray) -> None:
        """
        Average scalar field F along the polar rings (lat = -90, +90) over ocean longitudes
        and fill the entire polar ocean ring with that mean. Leave land points unchanged.
        Operates in-place.
        """
        # South pole row (index 0)
        j_s = 0
        ocean_j = ocean_mask[j_s, :]
        if np.any(ocean_j):
            mean_s = float(np.mean(F[j_s, ocean_j]))
            F[j_s, ocean_j] = mean_s

        # North pole row (last index)
        j_n = -1
        ocean_j = ocean_mask[j_n, :]
        if np.any(ocean_j):
            mean_n = float(np.mean(F[j_n, ocean_j]))
            F[j_n, ocean_j] = mean_n

    def _polar_vector_average_fill(self, u: np.ndarray, v: np.ndarray, ocean_mask: np.ndarray) -> None:
        """
        For vector field (u: east, v: north), transform each polar-ring vector into a common
        local tangent plane basis at the pole via 3D components, average, then transform back
        to each longitude's local (east, north) and fill across ocean longitudes. Land left unchanged.
        Operates in-place.
        """
        lons_rad = np.deg2rad(self.grid.lon)

        def basis_vectors(lam: np.ndarray, pole: str) -> tuple[np.ndarray, np.ndarray]:
            # East unit vector (independent of latitude)
            e_east = np.stack([-np.sin(lam), np.cos(lam), np.zeros_like(lam)], axis=1)
            if pole == "north":
                # e_north at +90°: (-cosλ, -sinλ, 0)
                e_north = np.stack([-np.cos(lam), -np.sin(lam), np.zeros_like(lam)], axis=1)
            else:
                # e_north at -90°: (cosλ, sinλ, 0)
                e_north = np.stack([np.cos(lam), np.sin(lam), np.zeros_like(lam)], axis=1)
            return e_east, e_north

        # Helper to process one pole
        def process_pole(j_idx: int, pole: str) -> None:
            mask_row = ocean_mask[j_idx, :]
            if not np.any(mask_row):
                return
            idx = np.where(mask_row)[0]
            lam_sel = lons_rad[idx]
            ee_sel, en_sel = basis_vectors(lam_sel, pole)
            u_sel = u[j_idx, idx]
            v_sel = v[j_idx, idx]
            # Assemble 3D vectors on tangent plane
            v3 = ee_sel * u_sel[:, None] + en_sel * v_sel[:, None]  # (N,3)
            v3_mean = np.mean(v3, axis=0)  # (3,)

            # Build basis for all longitudes to refill
            ee_all, en_all = basis_vectors(lons_rad, pole)
            u_fill = ee_all @ v3_mean
            v_fill = en_all @ v3_mean

            # Assign back only on ocean longitudes at this row
            u[j_idx, mask_row] = u_fill[mask_row]
            v[j_idx, mask_row] = v_fill[mask_row]

        # South pole (index 0) and North pole (index -1)
        process_pole(0, "south")
        process_pole(-1, "north")

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
        # Use relative wind (atmosphere minus surface current) for momentum exchange
        u_rel = u_atm - self.uo
        v_rel = v_atm - self.vo
        Va = np.sqrt(u_rel**2 + v_rel**2)
        Va_eff = np.minimum(Va, self.vcap)
        tau_x = self.tau_scale * (self.rho_a * self.CD * Va_eff * u_rel)
        tau_y = self.tau_scale * (self.rho_a * self.CD * Va_eff * v_rel)

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

            # Polar sponge: extra drag toward poles to suppress unrealistically fast polar currents
            try:
                lat_deg = np.abs(np.rad2deg(self.lat_rad))
                s = np.clip((lat_deg - self.polar_lat0) / max(1e-6, 90.0 - self.polar_lat0), 0.0, 1.0)
                r_extra = self.polar_gain * (s ** 2)  # smooth increase to pole
                self.uo -= sub_dt * r_extra * self.uo
                self.vo -= sub_dt * r_extra * self.vo
            except Exception:
                pass

            # 4) Scale-selective dissipation (cadence tied to outer step)
            if (self.diff_every > 0) and (self._step % self.diff_every == 0):
                # Latitude-adaptive K4 map so that (K4*dt/Δx(φ)^4) ≈ sigma4
                cos_map = self.coslat
                dx_lat = self.a * self.dlat                          # scalar
                dx_lon_map = self.a * self.dlon * cos_map            # 2D map
                dx_min_map = np.minimum(dx_lat, dx_lon_map)
                k4_map = self.sigma4 * (dx_min_map**4) / max(1e-12, sub_dt)

                # Allow optional scalar overrides via env; otherwise use maps
                k4_u = float(os.getenv("QD_OCEAN_K4_U")) if ("QD_OCEAN_K4_U" in os.environ) else k4_map
                k4_v = float(os.getenv("QD_OCEAN_K4_V")) if ("QD_OCEAN_K4_V" in os.environ) else k4_map
                k4_eta = float(os.getenv("QD_OCEAN_K4_ETA")) if ("QD_OCEAN_K4_ETA" in os.environ) else (0.5 * k4_map)

                self.uo = self._hyperdiffuse(self.uo, sub_dt, k4_u, n_substeps=self.k4_nsub)
                self.vo = self._hyperdiffuse(self.vo, sub_dt, k4_v, n_substeps=self.k4_nsub)
                self.eta = self._hyperdiffuse(self.eta, sub_dt, k4_eta, n_substeps=self.k4_nsub)

            # Optional Shapiro smoothing (off by default)
            if (self.shapiro_n > 0) and (self.shapiro_every > 0) and (self._step % self.shapiro_every == 0):
                self.uo = self._shapiro_filter(self.uo, n=self.shapiro_n)
                self.vo = self._shapiro_filter(self.vo, n=self.shapiro_n)
                self.eta = self._shapiro_filter(self.eta, n=self.shapiro_n)

            # 5) Continuity
            div = self.grid.divergence(self.uo, self.vo)
            self.eta += - sub_dt * self.H * div
            self.eta[on_land] = 0.0
            # Remove area-weighted ocean mean to avoid drift
            try:
                ocean_mask = (self.land_mask == 0)
                if np.any(ocean_mask):
                    w = np.maximum(np.cos(self.lat_rad), 0.0)
                    w_o = w * ocean_mask
                    eta_mean = float(np.sum(self.eta * w_o) / (np.sum(w_o) + 1e-15))
                    self.eta -= eta_mean
            except Exception:
                pass

            # 6) SST advection by ocean currents (semi-Lagrangian, gentle blend)
            adv_alpha = float(os.getenv("QD_OCEAN_ADV_ALPHA", "0.7"))  # 0..1
            Ts_adv = self._advect_scalar(self.Ts, self.uo, self.vo, sub_dt)
            self.Ts = (1.0 - adv_alpha) * self.Ts + adv_alpha * Ts_adv

            # Lateral diffusion of SST (explicit)
            if self.K_h > 0.0:
                self.Ts += sub_dt * self.K_h * self._laplacian_sphere(self.Ts)

            # 7) Optional vertical heat flux (enabled by default for coupled exchange)
            use_qnet = int(os.getenv("QD_OCEAN_USE_QNET", "1")) == 1
            if use_qnet and (Q_net is not None):
                heat_tendency = Q_net / (self.rho_w * self.cp_w * self.H)  # K/s
                ocean_all = (self.land_mask == 0)
                if ice_mask is not None:
                    open_mask = ocean_all & (~ice_mask)
                    ice_mask_local = ocean_all & (ice_mask)
                    # Allow reduced vertical exchange under ice rather than fully suppressing it,
                    # otherwise polar SST may not cool realistically. Tunable via env:
                    #   QD_OCEAN_ICE_QFAC in [0,1], default 0.2 (20% of open-ocean coupling).
                    ice_qfac = float(os.getenv("QD_OCEAN_ICE_QFAC", "0.2"))
                    Ts_new = self.Ts
                    Ts_new = np.where(open_mask, Ts_new + sub_dt * heat_tendency, Ts_new)
                    if ice_qfac > 0.0:
                        Ts_new = np.where(ice_mask_local, Ts_new + sub_dt * ice_qfac * heat_tendency, Ts_new)
                    self.Ts = Ts_new
                else:
                    self.Ts = np.where(ocean_all, self.Ts + sub_dt * heat_tendency, self.Ts)

            # 8) Sanity caps per substep (outlier handling)
            self.uo = np.nan_to_num(self.uo)
            self.vo = np.nan_to_num(self.vo)
            speed_vec = np.sqrt(self.uo**2 + self.vo**2)
            cap = float(self.max_u_cap)

            if self.outlier_method == "mean4":
                # Replace outliers (> cap) with the mean of 4 neighbors (N,S,E,W).
                u_n = np.roll(self.uo, -1, axis=0); u_s = np.roll(self.uo, 1, axis=0)
                u_e = np.roll(self.uo, -1, axis=1); u_w = np.roll(self.uo, 1, axis=1)
                v_n = np.roll(self.vo, -1, axis=0); v_s = np.roll(self.vo, 1, axis=0)
                v_e = np.roll(self.vo, -1, axis=1); v_w = np.roll(self.vo, 1, axis=1)
                u_mean4 = 0.25 * (u_n + u_s + u_e + u_w)
                v_mean4 = 0.25 * (v_n + v_s + v_e + v_w)
                mask_fast = speed_vec > cap
                self.uo = np.where(mask_fast, u_mean4, self.uo)
                self.vo = np.where(mask_fast, v_mean4, self.vo)
                # Gentle safety clamp in case mean4 still exceeds cap slightly
                speed_vec2 = np.sqrt(self.uo**2 + self.vo**2)
                scale2 = np.where(speed_vec2 > cap, cap / (speed_vec2 + 1e-12), 1.0)
                self.uo *= scale2
                self.vo *= scale2
            else:
                # Fallback: vector clamp (previous behavior)
                scale = np.where(speed_vec > cap, cap / (speed_vec + 1e-12), 1.0)
                self.uo *= scale
                self.vo *= scale

            self.eta = np.nan_to_num(self.eta)
            # Clamp eta to a safe anomaly range (meters) to prevent runaway
            _eta_cap = float(os.getenv("QD_ETA_CAP", "5.0"))
            self.eta = np.clip(self.eta, -_eta_cap, _eta_cap)
            self.Ts = np.nan_to_num(self.Ts)
            # Clamp eta to a safe anomaly range (meters) to prevent runaway
            _eta_cap = float(os.getenv("QD_ETA_CAP", "5.0"))
            self.eta = np.clip(self.eta, -_eta_cap, _eta_cap)
            self.Ts = np.nan_to_num(self.Ts)

        # Ocean energy diagnostics (global and polar), compare implied d<Ts>/dt with effective Q_net
        try:
            if int(os.getenv("QD_OCEAN_ENERGY_DIAG", "1")) == 1:
                every = int(os.getenv("QD_OCEAN_DIAG_EVERY", "200"))
                if every <= 0:
                    every = 200
                if (self._step % every) == 0:
                    # Area weights
                    w = np.maximum(np.cos(self.lat_rad), 0.0)
                    wsum_ocean = float(np.sum(w * (self.land_mask == 0)) + 1e-15)

                    # Effective surface heat flux over ocean used this step (W/m^2):
                    #   open ocean: Q_net
                    #   under-ice:  ice_qfac * Q_net  (if ice mask provided)
                    if Q_net is not None:
                        if ice_mask is not None:
                            ice_qfac = float(os.getenv("QD_OCEAN_ICE_QFAC", "0.2"))
                            eff_Q = np.where((self.land_mask == 0) & (~ice_mask), Q_net, 0.0)
                            if ice_qfac > 0.0:
                                eff_Q += np.where((self.land_mask == 0) & (ice_mask), ice_qfac * Q_net, 0.0)
                        else:
                            eff_Q = np.where((self.land_mask == 0), Q_net, 0.0)
                        Q_mean = float(np.sum(eff_Q * w) / wsum_ocean)
                    else:
                        Q_mean = 0.0

                    # Implied net heat flux from prognosed Ts tendency over this outer step
                    # (includes advection/diffusion/any caps): ρ c_p H d<Ts>/dt
                    # For robustness, compare Ts at start-of-step vs end-of-step
                    # Cache last Ts if available
                    if not hasattr(self, "_Ts_prev_for_diag"):
                        self._Ts_prev_for_diag = self.Ts.copy()
                        implied = 0.0
                        resid = 0.0
                    else:
                        dT = (self.Ts - self._Ts_prev_for_diag) / max(1e-12, dt)
                        dT_mean = float(np.sum(dT * w * (self.land_mask == 0)) / wsum_ocean)
                        implied = float(self.rho_w * self.cp_w * self.H * dT_mean)
                        resid = implied - Q_mean
                        self._Ts_prev_for_diag = self.Ts.copy()

                    # Polar band diagnostics (|lat| >= 60°)
                    lat_deg = np.abs(np.rad2deg(self.lat_rad))
                    polar_mask = (lat_deg >= float(os.getenv("QD_OCEAN_POLAR_LAT", "60.0"))) & (self.land_mask == 0)
                    if np.any(polar_mask):
                        w_p = w * polar_mask
                        wsum_p = float(np.sum(w_p) + 1e-15)
                        if Q_net is not None:
                            if ice_mask is not None:
                                ice_qfac = float(os.getenv("QD_OCEAN_ICE_QFAC", "0.2"))
                                eff_Qp = np.where(polar_mask & (~ice_mask), Q_net, 0.0)
                                if ice_qfac > 0.0:
                                    eff_Qp += np.where(polar_mask & (ice_mask), ice_qfac * Q_net, 0.0)
                            else:
                                eff_Qp = np.where(polar_mask, Q_net, 0.0)
                            Qp_mean = float(np.sum(eff_Qp * w) / wsum_p)
                        else:
                            Qp_mean = 0.0
                        dTp = (self.Ts - getattr(self, "_Ts_prev_for_diag_p", self.Ts)) / max(1e-12, dt)
                        dTp_mean = float(np.sum(dTp * w * polar_mask) / wsum_p)
                        implied_p = float(self.rho_w * self.cp_w * self.H * dTp_mean)
                        resid_p = implied_p - Qp_mean
                        self._Ts_prev_for_diag_p = self.Ts.copy()
                    else:
                        Qp_mean = implied_p = resid_p = 0.0

                    print(f"[OceanE] ⟨Q_net⟩={Q_mean:+.2f} W/m^2 | implied={implied:+.2f} | resid={resid:+.2f}  "
                          f"|| Polar(|lat|>={int(float(os.getenv('QD_OCEAN_POLAR_LAT','60')))}°): "
                          f"⟨Q⟩={Qp_mean:+.2f}, implied={implied_p:+.2f}, resid={resid_p:+.2f}")
        except Exception:
            pass

        # Polar corrections at ±90°: average along longitude ring and refill
        if int(os.getenv("QD_OCEAN_POLAR_FIX", "1")) == 1:
            ocean_mask = (self.land_mask == 0)
            try:
                # Scalars: SST
                self._polar_scalar_average_fill(self.Ts, ocean_mask)
                # Vectors: ocean currents (uo, vo)
                self._polar_vector_average_fill(self.uo, self.vo, ocean_mask)
            except Exception:
                # Be fail-safe: skip polar correction on any unexpected numerical issue
                pass

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

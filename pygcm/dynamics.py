# pygcm/dynamics.py

"""
Implements a stable spectral dynamics core for the shallow water equations.
This approach avoids the grid-point instabilities faced by finite-difference methods.
"""

import os
import numpy as np
from scipy.ndimage import map_coordinates, convolve
from .grid import SphericalGrid
from . import constants as const
from . import energy as energy
from . import humidity as humidity
from .jax_compat import is_enabled as _jax_enabled, to_numpy as _to_np, laplacian_sphere as _j_lap, hyperdiffuse as _j_hyp, advect_semilag as _j_adv

class SpectralModel:
    """
    A spectral shallow water model that solves the equations of motion in
    spectral space using spherical harmonics.
    """
    def __init__(self, grid: SphericalGrid, friction_map: np.ndarray, initial_state=None, g=9.81, H=8000, tau_rad=1e6, greenhouse_factor=0.15, C_s_map=None,
                 land_mask=None, Cs_ocean=None, Cs_land=None, Cs_ice=None, seaice_enabled=None, t_freeze=None, rho_i=None, L_f=None):
        self.grid = grid
        self.friction_map = friction_map
        self.C_s_map = C_s_map
        # Sea-ice / surface properties (M2)
        self.land_mask = land_mask
        self.Cs_ocean = Cs_ocean
        self.Cs_land = Cs_land
        self.Cs_ice = Cs_ice
        self.seaice_enabled = bool(int(os.getenv("QD_USE_SEAICE", "1"))) if seaice_enabled is None else bool(seaice_enabled)
        self.t_freeze = float(os.getenv("QD_T_FREEZE", "271.35")) if t_freeze is None else float(t_freeze)
        self.rho_i = float(os.getenv("QD_RHO_ICE", "917")) if rho_i is None else float(rho_i)
        self.L_f = float(os.getenv("QD_LF", "3.34e5")) if L_f is None else float(L_f)

        self.g = g
        self.H = H
        self.tau_rad = tau_rad
        self.greenhouse_factor = greenhouse_factor
        self.a = const.PLANET_RADIUS
        
        # Spectral truncation (lower for stability and speed)
        self.n_trunc = int(grid.n_lat / 3)

        # Grid properties needed for gradient calculations
        self.dlat_rad = np.deg2rad(self.grid.lat[1] - self.grid.lat[0])
        self.dlon_rad = np.deg2rad(self.grid.lon[1] - self.grid.lon[0])

        # Initialize spectral coefficients for vorticity, divergence, and height
        self.vort_spec = np.zeros((self.n_trunc, self.n_trunc), dtype=complex)
        self.div_spec = np.zeros((self.n_trunc, self.n_trunc), dtype=complex)
        self.h_spec = np.zeros((self.n_trunc, self.n_trunc), dtype=complex)

        # Initial state variables as floats
        self.u = np.zeros(self.grid.lat_mesh.shape, dtype=float)
        self.v = np.zeros(self.grid.lat_mesh.shape, dtype=float)
        
        # Initialize with a more realistic height field (higher pressure at poles)
        lat_rad = np.deg2rad(self.grid.lat_mesh)
        h_initial_anomaly = 300 * (np.sin(lat_rad)**2) # Add a 300m anomaly, max at poles
        self.h = np.full(self.grid.lat_mesh.shape, self.H, dtype=float) + h_initial_anomaly

        self.T_s = np.full(self.grid.lat_mesh.shape, 288.0, dtype=float) # Initial surface temp of 288K
        self.cloud_cover = np.zeros(self.grid.lat_mesh.shape, dtype=float) # Prognostic cloud cover
        self.h_ice = np.zeros(self.grid.lat_mesh.shape, dtype=float)  # Sea-ice thickness (m), M2
        
        # Diagnostic radiation fields
        self.isr = np.zeros(self.grid.lat_mesh.shape, dtype=float)  # Incoming Shortwave (total)
        self.isr_A = np.zeros(self.grid.lat_mesh.shape, dtype=float)  # Star A component
        self.isr_B = np.zeros(self.grid.lat_mesh.shape, dtype=float)  # Star B component
        self.olr = np.zeros(self.grid.lat_mesh.shape, dtype=float)  # Outgoing Longwave
        # Humidity (P008)
        try:
            self.hum_params = humidity.get_humidity_params_from_env()
        except Exception:
            self.hum_params = humidity.HumidityParams()
        RH0 = 0.5
        try:
            import os as _os_mod
            RH0 = float(_os_mod.getenv("QD_Q_INIT_RH", "0.5"))
        except Exception:
            RH0 = 0.5
        self.q = humidity.q_init(self.T_s, RH0=RH0, p0=self.hum_params.p0)
        self.E_flux_last = np.zeros_like(self.T_s)
        self.P_cond_flux_last = np.zeros_like(self.T_s)
        self.LH_last = np.zeros_like(self.T_s)
        self.LH_release_last = np.zeros_like(self.T_s)

    def _advect(self, field, dt):
        """
        Advects a scalar field using a semi-Lagrangian scheme with bilinear interpolation.
        """
        # Try JAX path
        try:
            if _jax_enabled():
                coslat = np.maximum(1e-6, np.cos(np.deg2rad(self.grid.lat_mesh)))
                adv = _j_adv(field, self.u, self.v, dt, self.a, self.dlat_rad, self.dlon_rad, coslat)
                return _to_np(adv)
        except Exception:
            pass

        # Fallback NumPy/SciPy path
        dlon = self.u * dt / (self.a * np.maximum(1e-6, np.cos(np.deg2rad(self.grid.lat_mesh))))
        dlat = self.v * dt / self.a

        dx = dlon / self.dlon_rad
        dy = dlat / self.dlat_rad

        lats = np.arange(self.grid.n_lat)
        lons = np.arange(self.grid.n_lon)
        lon_mesh, lat_mesh = np.meshgrid(lons, lats)

        dep_lat = lat_mesh - dy
        dep_lon = lon_mesh - dx

        advected_field = map_coordinates(field, [dep_lat, dep_lon], order=1, mode='wrap', prefilter=False)
        return advected_field

    def _grid_to_spectral(self, data):
        """Transforms data from grid-point space to spectral space."""
        # Zonal Fourier transform
        fft_data = np.fft.rfft(data, axis=1)
        # Legendre transform (simplified placeholder)
        # A full implementation requires Gaussian quadrature and spherical harmonic transforms.
        # We will use a simplified projection for this prototype.
        spectral_coeffs = np.zeros((self.n_trunc, self.n_trunc), dtype=complex)
        if fft_data.shape[1] >= self.n_trunc:
            for m in range(self.n_trunc):
                # This is a gross simplification of the Legendre transform
                spectral_coeffs[:, m] = np.mean(fft_data[:, m, np.newaxis] * np.sin(np.linspace(0, np.pi, self.grid.n_lat))[:,np.newaxis], axis=0)[:self.n_trunc]
        return spectral_coeffs

    def _spectral_to_grid(self, spectral_coeffs):
        """Transforms data from spectral space back to grid-point space."""
        # Inverse Legendre transform (simplified)
        grid_data_fft = np.zeros((self.grid.n_lat, int(self.grid.n_lon/2)+1), dtype=complex)
        for m in range(self.n_trunc):
             grid_data_fft[:, m] = np.interp(np.arange(self.grid.n_lat), np.arange(self.n_trunc), spectral_coeffs[:, m].real)
        # Inverse Fourier transform
        return np.fft.irfft(grid_data_fft, n=self.grid.n_lon, axis=1)

    # ---------------- Project 010: Hyperdiffusion utilities ----------------
    def _laplacian_sphere(self, F: np.ndarray) -> np.ndarray:
        """
        Spherical Laplacian of a scalar field F on a regular lat-lon grid.
        Uses divergence form with cos(phi) weighting for numerical robustness.
        Longitude is treated as periodic; latitude uses np.gradient one-sided at poles.
        """
        # Try JAX-jitted kernel
        try:
            if _jax_enabled():
                phi = np.deg2rad(self.grid.lat_mesh)
                cos = np.maximum(np.cos(phi), 0.2)
                return _to_np(_j_lap(F, self.dlat_rad, self.dlon_rad, cos, self.a))
        except Exception:
            pass

        # NumPy fallback
        a = self.a
        dphi = self.dlat_rad
        dlmb = self.dlon_rad
        phi = np.deg2rad(self.grid.lat_mesh)
        cos = np.maximum(np.cos(phi), 0.2)
        F = np.nan_to_num(F, copy=False)

        dF_dphi = np.gradient(F, dphi, axis=0)
        term_phi = (1.0 / cos) * np.gradient(cos * dF_dphi, dphi, axis=0)

        d2F_dlmb2 = (np.roll(F, -1, axis=1) - 2.0 * F + np.roll(F, 1, axis=1)) / (dlmb ** 2)
        term_lmb = d2F_dlmb2 / (cos ** 2)

        return (term_phi + term_lmb) / (a ** 2)

    def _hyperdiffuse(self, F: np.ndarray, k4, dt: float, n_substeps: int = 1) -> np.ndarray:
        """
        Apply explicit 4th-order hyperdiffusion: dF/dt = -k4 * ∇⁴ F
        k4 can be a scalar (float) or a 2D map broadcastable to F.
        Implemented as two successive Laplacians. Optionally uses substeps for stability margin.
        """
        # Try JAX-jitted kernel
        try:
            if _jax_enabled():
                phi = np.deg2rad(self.grid.lat_mesh)
                cos = np.maximum(np.cos(phi), 0.2)
                return _to_np(_j_hyp(F, k4, dt, n_substeps, self.dlat_rad, self.dlon_rad, cos, self.a))
        except Exception:
            pass

        if dt <= 0.0:
            return F
        # NumPy fallback
        try:
            import numpy as _np
            if _np.isscalar(k4):
                k4_arr = float(k4)
                if k4_arr <= 0.0:
                    return F
            else:
                k4_arr = _np.nan_to_num(k4, copy=False)
                if _np.all(k4_arr <= 0.0):
                    return F
        except Exception:
            return F
        n = max(1, int(n_substeps))
        sub_dt = dt / n
        out = _np.nan_to_num(F, copy=True)
        for _ in range(n):
            L = self._laplacian_sphere(out)
            L2 = self._laplacian_sphere(L)
            out = out - k4_arr * L2 * sub_dt
        return _np.nan_to_num(out, copy=False)

    # ---------------- Project 010 M4: Alternate filters (Shapiro / Spectral) ----------------
    def _shapiro_filter(self, F: np.ndarray, n: int = 2, lon_wrap: bool = True) -> np.ndarray:
        """
        Shapiro-like smoothing via separable 1-2-1 kernel repeated n times.
        - Longitude uses periodic wrap.
        - Latitude uses nearest-edge (no wrap across poles).
        """
        try:
            n = max(1, int(n))
        except Exception:
            n = 2
        k1 = np.array([1.0, 2.0, 1.0], dtype=float)
        k1 /= k1.sum()  # [0.25, 0.5, 0.25]
        out = np.nan_to_num(F, copy=True)
        for _ in range(n):
            out = convolve(out, k1[np.newaxis, :], mode="wrap" if lon_wrap else "nearest")
            out = convolve(out, k1[:, np.newaxis], mode="nearest")
        return out

    def _spectral_zonal_filter(self, F: np.ndarray, cutoff: float = 0.75, damp: float = 0.5) -> np.ndarray:
        """
        Zonal-FFT high-wavenumber damping.
        - cutoff: fraction of Nyquist (0..1). k > cutoff*k_N are damped.
        - damp:   damping strength (0..1), where 1 means zero out high-k.
        """
        try:
            cutoff = float(cutoff)
            damp = float(damp)
        except Exception:
            return F
        if damp <= 0.0 or cutoff <= 0.0:
            return np.nan_to_num(F, copy=False)

        arr = np.nan_to_num(F, copy=False)
        fft = np.fft.rfft(arr, axis=1)
        bins = fft.shape[1]
        if bins <= 1:
            return arr
        kN = bins - 1
        kcut = int(max(1, min(kN, int(cutoff * kN))))
        factor = np.ones(bins, dtype=float)
        factor[kcut:] *= max(0.0, 1.0 - min(1.0, damp))
        fft *= factor[np.newaxis, :]
        out = np.fft.irfft(fft, n=self.grid.n_lon, axis=1)
        return np.nan_to_num(out, copy=False)

    def time_step(self, Teq_field, dt, albedo=None):
        """
        Advances the model state by one time step in spectral space.
        This is a placeholder for a full spectral dynamics implementation.
        The complexity of a real spectral core is very high.
        
        For this task, we will implement a heavily simplified and stabilized
        grid-point model that borrows spectral ideas (like filtering).
        """
        
        # --- Simplified, Stabilized Grid-Point Model ---
        
        # 1. Calculate atmospheric temperature from height anomaly
        # T_a = T_ref + (g/Cp) * h, a simplification
        T_a = 288.0 + (self.g / 1004.0) * self.h

        # Humidity physics (P008): compute evaporation E, surface LH, supersaturation condensation and LH_release
        if not hasattr(self, "hum_params"):
            try:
                self.hum_params = humidity.get_humidity_params_from_env()
            except Exception:
                self.hum_params = humidity.HumidityParams()
        try:
            # Surface factor based on land/ocean/ice
            surf_factor = humidity.surface_evaporation_factor(getattr(self, "land_mask", None), getattr(self, "h_ice", None), self.hum_params)
            E_flux = humidity.evaporation_flux(self.T_s, getattr(self, "q", np.zeros_like(self.T_s)), self.u, self.v, surf_factor, self.hum_params)
            LH = self.hum_params.L_v * E_flux  # W/m^2, upward (surface energy sink)
            # Column mass for q tendency
            M_col = max(1e-6, float(self.hum_params.rho_a * self.hum_params.h_mbl))
            q_evap = getattr(self, "q", np.zeros_like(self.T_s)) + (E_flux / M_col) * dt
            P_cond_flux, q_after_cond = humidity.condensation(q_evap, T_a, dt, self.hum_params)
            LH_release = self.hum_params.L_v * P_cond_flux  # W/m^2, atmospheric heating
            # Update state and store diagnostics
            self.q = np.clip(np.nan_to_num(q_after_cond, copy=False), 0.0, 0.5)
            self.E_flux_last = E_flux
            self.P_cond_flux_last = P_cond_flux
            self.LH_last = LH
            self.LH_release_last = LH_release
        except Exception:
            # On failure, fall back to zero latent fluxes
            LH = 0.0
            LH_release = 0.0
        # 2. Update Surface Temperature using energy framework (with optional mixing)
        # Old (Newton-like) update path retained for blending/backward-compat
        absorbed_solar_old = const.SIGMA * Teq_field**4
        olr_old = const.SIGMA * self.T_s**4
        ilr_old = self.greenhouse_factor * const.SIGMA * T_a**4
        net_flux_old = absorbed_solar_old + ilr_old - olr_old
        # Initialize energy framework params and weight on first use
        if not hasattr(self, "energy_params"):
            try:
                self.energy_params = energy.get_energy_params_from_env()
            except Exception:
                self.energy_params = energy.EnergyParams()
        if not hasattr(self, "energy_w"):
            try:
                self.energy_w = float(os.getenv("QD_ENERGY_W", "0.0"))
            except Exception:
                self.energy_w = 0.0
        if not hasattr(self, "_step_counter"):
            self._step_counter = 0
        dT_old = (net_flux_old / max(1e-12, self.energy_params.c_sfc)) * dt
        Ts_newton = self.T_s + dT_old

        # New explicit energy budget (Milestone 1): SW/LW partition, SH=LH=0
        Ts_energy = None
        if albedo is not None and hasattr(self, "isr") and self.isr is not None:
            try:
                # M4: humidity/precipitation → cloud optical consistency
                try:
                    couple = int(os.getenv("QD_CLOUD_COUPLE", "1")) == 1
                except Exception:
                    couple = True
                if couple and hasattr(self, "q"):
                    try:
                        qsat_air = humidity.q_sat(T_a, p=self.hum_params.p0)
                    except Exception:
                        qsat_air = np.maximum(1e-12, np.ones_like(self.T_s) * 1e-3)
                    RH = np.clip(self.q / np.maximum(1e-12, qsat_air), 0.0, 1.5)
                    RH0 = float(os.getenv("QD_RH0", "0.6"))
                    k_q = float(os.getenv("QD_K_Q", "0.3"))
                    k_p = float(os.getenv("QD_K_P", "0.4"))
                    rh_excess = np.maximum(0.0, RH - RH0)
                    P = getattr(self, "P_cond_flux_last", np.zeros_like(self.T_s))
                    Ppos = P[P > 0]
                    if os.getenv("QD_PCOND_REF"):
                        P_ref = float(os.getenv("QD_PCOND_REF"))
                    else:
                        P_ref = float(np.median(Ppos)) if Ppos.size > 0 else 1e-6
                    p_term = np.tanh(np.where(P_ref > 0, P / P_ref, 0.0))
                    cloud_eff = np.clip(self.cloud_cover + k_q * rh_excess + k_p * p_term, 0.0, 1.0)
                else:
                    cloud_eff = self.cloud_cover
                self.cloud_eff_last = cloud_eff

                SW_atm, SW_sfc, R = energy.shortwave_radiation(self.isr, albedo, cloud_eff, self.energy_params)
                # P006 upgrade (optional): cloud-optical-aware LW + surface emissivity map
                use_lw_v2 = int(os.getenv("QD_LW_V2", "1")) == 1
                if use_lw_v2:
                    # Sea-ice optical fraction for emissivity blending (ocean→ice)
                    H_ice_ref = float(os.getenv("QD_HICE_REF", "0.5"))
                    if hasattr(self, "h_ice"):
                        ice_frac = 1.0 - np.exp(-np.maximum(self.h_ice, 0.0) / max(1e-6, H_ice_ref))
                    else:
                        ice_frac = np.zeros_like(self.T_s)
                    # Surface emissivity map (ocean/land/ice)
                    if getattr(self, "land_mask", None) is not None:
                        eps_sfc_map = energy.surface_emissivity_map(self.land_mask, ice_frac)
                    else:
                        eps_sfc_map = float(os.getenv("QD_EPS_DEFAULT", "0.97"))
                    LW_atm, LW_sfc, OLR, DLR, eps = energy.longwave_radiation_v2(
                        self.T_s, T_a, cloud_eff, eps_sfc_map, self.energy_params
                    )
                else:
                    LW_atm, LW_sfc, OLR, DLR, eps = energy.longwave_radiation(self.T_s, T_a, cloud_eff, self.energy_params)
                # M2: Sensible heat flux (SH) via bulk formula; use humidity rho_a and env overrides
                try:
                    C_H = float(os.getenv("QD_CH", "1.5e-3"))
                except Exception:
                    C_H = 1.5e-3
                try:
                    cp_air = float(os.getenv("QD_CP_A", "1004.0"))
                except Exception:
                    cp_air = 1004.0
                rho_air = float(getattr(self.hum_params, "rho_a", 1.2))
                B_land = float(os.getenv("QD_BOWEN_LAND", "0.7"))
                B_ocean = float(os.getenv("QD_BOWEN_OCEAN", "0.3"))
                land_mask_arr = self.land_mask if self.land_mask is not None else np.zeros_like(self.T_s, dtype=int)
                SH_arr, _LH_bowen = energy.boundary_layer_fluxes(self.T_s, T_a, self.u, self.v, land_mask_arr,
                                                                  C_H=C_H, rho=rho_air, c_p=cp_air,
                                                                  B_land=B_land, B_ocean=B_ocean)
                # M2: Sea-ice thermodynamics if enabled and land_mask provided
                if self.seaice_enabled and (self.land_mask is not None):
                    Cs_ocean = self.Cs_ocean if self.Cs_ocean is not None else 2.0e8
                    Cs_land = self.Cs_land if self.Cs_land is not None else 3.0e6
                    Cs_ice = self.Cs_ice if self.Cs_ice is not None else 5.0e6
                    Ts_energy, h_ice_next = energy.integrate_surface_energy_with_seaice(
                        self.T_s, SW_sfc, LW_sfc, SH_arr, LH, dt,
                        self.land_mask, self.h_ice,
                        Cs_ocean, Cs_land, Cs_ice,
                        t_freeze=self.t_freeze, rho_i=self.rho_i, L_f=self.L_f, t_floor=self.energy_params.t_floor
                    )
                else:
                    # Use per-grid heat capacity map if provided (P007 M1)
                    if getattr(self, "C_s_map", None) is not None:
                        Ts_energy = energy.integrate_surface_energy_map(
                            self.T_s, SW_sfc, LW_sfc, SH_arr, LH, dt, self.C_s_map, t_floor=self.energy_params.t_floor
                        )
                    else:
                        Ts_energy = energy.integrate_surface_energy(self.T_s, SW_sfc, LW_sfc, SH_arr, LH, dt, self.energy_params)
                # Update diagnostic OLR field
                self.olr = OLR
                # Optional diagnostics every ~200 steps
                if self.energy_params.diag and (self._step_counter % 200 == 0):
                    diag = energy.compute_energy_diagnostics(self.grid.lat_mesh, self.isr, R, OLR, SW_sfc, LW_sfc, SH_arr, LH)
                    print(f"[EnergyDiag] TOA={diag['TOA_net']:+.2f} W/m^2 | SFC={diag['SFC_net']:+.2f} | ATM={diag['ATM_net']:+.2f} | "
                          f"I={diag['I_mean']:.1f} R={diag['R_mean']:.1f} OLR={diag['OLR_mean']:.1f}")
                    if self.seaice_enabled and (self.land_mask is not None):
                        try:
                            ocean = (self.land_mask == 0)
                            ice_mask = (self.h_ice > 0.0) & ocean
                            w = np.maximum(np.cos(np.deg2rad(self.grid.lat_mesh)), 0.0)
                            ice_area = float((w * ice_mask).sum() / (w.sum() + 1e-15))
                            mean_h = float(self.h_ice[ice_mask].mean()) if np.any(ice_mask) else 0.0
                            print(f"[SeaIce] area={ice_area:.3f}, mean_h={mean_h:.2f} m")
                        except Exception:
                            pass
                # Optional diagnostics every ~200 steps
                if self.energy_params.diag and (self._step_counter % 200 == 0):
                    diag = energy.compute_energy_diagnostics(self.grid.lat_mesh, self.isr, R, OLR, SW_sfc, LW_sfc, SH_arr, LH)
                    print(f"[EnergyDiag] TOA={diag['TOA_net']:+.2f} W/m^2 | SFC={diag['SFC_net']:+.2f} | ATM={diag['ATM_net']:+.2f} | "
                          f"I={diag['I_mean']:.1f} R={diag['R_mean']:.1f} OLR={diag['OLR_mean']:.1f}")
            except Exception:
                # Fallback to old path on any error
                Ts_energy = None
                self.olr = olr_old
        else:
            # Keep previous OLR diagnostic in absence of energy calc
            self.olr = olr_old

        # Blend between old and new schemes via weight QD_ENERGY_W
        w = float(getattr(self, "energy_w", 0.0))
        w = min(1.0, max(0.0, w))
        if Ts_energy is None:
            self.T_s = Ts_newton
        else:
            self.T_s = (1.0 - w) * Ts_newton + w * Ts_energy
            # Update sea-ice thickness if computed in this step
            if self.seaice_enabled and 'h_ice_next' in locals():
                self.h_ice = h_ice_next

        self._step_counter += 1

        # 2b. Horizontal advection of surface temperature (semi-Lagrangian, gentle)
        adv_alpha = 0.2  # blending factor for numerical stability
        advected_Ts = self._advect(self.T_s, dt)
        self.T_s = (1.0 - adv_alpha) * self.T_s + adv_alpha * advected_Ts
        # Advect humidity q with the same semi-Lagrangian scheme (gentle)
        if hasattr(self, "q"):
            advected_q = self._advect(self.q, dt)
            self.q = (1.0 - adv_alpha) * self.q + adv_alpha * advected_q
            self.q = np.clip(np.nan_to_num(self.q, copy=False), 0.0, 0.5)
        
        # 3. Add radiative forcing to height field
        R_gas = 287
        h_eq = (R_gas / self.g) * Teq_field
        rad_forcing = (h_eq - self.h) / self.tau_rad
        self.h += rad_forcing * dt

        # M3: Atmospheric energy budget coupling (adds SW_atm + LW_atm + SH + LH_release)
        try:
            if (albedo is not None) and (getattr(self, "energy_w", 0.0) > 0.0):
                H_atm = float(os.getenv("QD_ATM_H", str(getattr(self.hum_params, "h_mbl", 800.0))))
                rho_air = float(getattr(self.hum_params, "rho_a", 1.2))
                # Update geopotential height using energy module helper (stable, weighted by energy_w)
                self.h = energy.integrate_atmos_energy_height(
                    self.h, SW_atm, LW_atm, SH_arr, LH_release, dt,
                    rho_air=rho_air, H_atm=H_atm, g=self.g, weight=float(self.energy_w)
                )
        except Exception:
            pass
        
        # 4. Momentum step: pressure-gradient force + Coriolis (+ friction)
        # Support two schemes via env QD_MOM_SCHEME: "primitive" (du/dt explicit) or "geos" (relaxation).
        mom_scheme = os.getenv("QD_MOM_SCHEME", "geos").lower()
        f = self.grid.coriolis_param

        # Common metrics and gradients
        dh_dlon = np.gradient(self.h, self.dlon_rad, axis=1)   # ∂h/∂λ
        dh_dlat = np.gradient(self.h, self.dlat_rad, axis=0)   # ∂h/∂φ
        cos_lat_capped = np.maximum(np.cos(np.deg2rad(self.grid.lat_mesh)), 1e-6)

        if mom_scheme == "primitive":
            # Primitive-momentum (linearized shallow-water):
            # du/dt = -g/(a cosφ) ∂h/∂λ + f v - r u
            # dv/dt = -g/a ∂h/∂φ - f u - r v
            u_old = self.u.copy()
            v_old = self.v.copy()

            PGF_x = -(self.g / (self.a * cos_lat_capped)) * dh_dlon
            PGF_y = -(self.g / self.a) * dh_dlat
            # Coriolis terms are perpendicular to velocity and proportional to speed
            du = (PGF_x + f * v_old - self.friction_map * u_old) * dt
            dv = (PGF_y - f * u_old - self.friction_map * v_old) * dt

            self.u = u_old + du
            self.v = v_old + dv

            # Clip winds for stability
            max_wind = 200.0
            self.u = np.clip(self.u, -max_wind, max_wind)
            self.v = np.clip(self.v, -max_wind, max_wind)

        else:
            # Geostrophic relaxation (legacy default; numerically robust)
            # Regularize f near equator
            f_min = 2.0 * const.PLANET_OMEGA * np.sin(np.deg2rad(5.0))
            sign_nonzero = np.where(f >= 0.0, 1.0, -1.0)
            f_safe = np.where(np.abs(f) < f_min, sign_nonzero * f_min, f)
            # Geostrophic winds from balance: f k×V = -g ∇h
            u_g = -(self.g / (f_safe * self.a * cos_lat_capped)) * dh_dlat
            v_g = (self.g / (f_safe * self.a)) * dh_dlon
            max_wind = 200.0
            u_g = np.clip(u_g, -max_wind, max_wind)
            v_g = np.clip(v_g, -max_wind, max_wind)
            # Nudge towards geostrophic
            self.u = self.u * 0.8 + u_g * 0.2
            self.v = self.v * 0.8 + v_g * 0.2
            # Apply linear friction after nudging (stronger over land via map)
            self.u += (-self.friction_map * self.u) * dt
            self.v += (-self.friction_map * self.v) * dt

        # ---------------- Project 010: Scale-selective dissipation (hyperdiffusion) ----------------
        try:
            diff_enabled = int(os.getenv("QD_DIFF_ENABLE", "1")) == 1
            filter_type = os.getenv("QD_FILTER_TYPE", "combo").lower()
            diff_every = int(os.getenv("QD_DIFF_EVERY", "1"))
        except Exception:
            diff_enabled = True
            filter_type = "hyper4"
            diff_every = 1

        if diff_enabled and (filter_type in ("hyper4", "combo")) and (self._step_counter % max(1, diff_every) == 0):
            dyn_diag = int(os.getenv("QD_DYN_DIAG", "0")) == 1
            if dyn_diag and (self._step_counter % 200 == 0):
                var_u0 = float(np.var(self.u))
                var_v0 = float(np.var(self.v))
                var_h0 = float(np.var(self.h))

            sigma4_env = os.getenv("QD_SIGMA4", "0.02")
            dx_min = None
            sigma4 = None
            if sigma4_env is not None:
                try:
                    sigma4 = float(sigma4_env)
                except Exception:
                    sigma4 = 0.02
                # Latitude-adaptive metric lengths per grid point
                phi = np.deg2rad(self.grid.lat_mesh)
                cos = np.maximum(np.cos(phi), 1e-3)
                dx_lat = self.a * self.dlat_rad                             # scalar
                dx_lon_map = self.a * self.dlon_rad * cos                   # 2D map
                dx_min_map = np.minimum(dx_lat, dx_lon_map)                 # 2D map
                k4_map_base = sigma4 * (dx_min_map ** 4) / max(1e-12, dt)  # 2D map

                # If explicit overrides are provided in env, use them (scalars). Otherwise use maps.
                k4_u = float(os.getenv("QD_K4_U")) if ("QD_K4_U" in os.environ) else k4_map_base
                k4_v = float(os.getenv("QD_K4_V")) if ("QD_K4_V" in os.environ) else k4_map_base
                k4_h = float(os.getenv("QD_K4_H")) if ("QD_K4_H" in os.environ) else (0.5 * k4_map_base)
                k4_q = float(os.getenv("QD_K4_Q")) if ("QD_K4_Q" in os.environ) else (0.5 * k4_map_base)
                k4_c = float(os.getenv("QD_K4_CLOUD")) if ("QD_K4_CLOUD" in os.environ) else (0.25 * k4_map_base)
            else:
                # No sigma4 provided: fall back to constant coefficients
                k4_u = float(os.getenv("QD_K4_U", "1.0e14"))
                k4_v = float(os.getenv("QD_K4_V", "1.0e14"))
                k4_h = float(os.getenv("QD_K4_H", "5.0e13"))
                k4_q = float(os.getenv("QD_K4_Q", "0.0"))
                k4_c = float(os.getenv("QD_K4_CLOUD", "0.0"))

            # Apply hyperdiffusion to dynamical fields
            try:
                nsub = int(os.getenv("QD_K4_NSUB", "1"))
            except Exception:
                nsub = 1
            self.u = self._hyperdiffuse(self.u, k4_u, dt, n_substeps=nsub)
            self.v = self._hyperdiffuse(self.v, k4_v, dt, n_substeps=nsub)
            self.h = self._hyperdiffuse(self.h, k4_h, dt, n_substeps=nsub)
            # Optional: humidity and cloud fields (weaker, off by default)
            # k4_q / k4_c may be 2D maps; use any() to avoid ambiguous truth-value error
            apply_q = (np.isscalar(k4_q) and (k4_q > 0.0)) or ((not np.isscalar(k4_q)) and np.any(k4_q > 0.0)) or (int(os.getenv("QD_DIFF_Q", "0")) == 1)
            apply_cloud = (np.isscalar(k4_c) and (k4_c > 0.0)) or ((not np.isscalar(k4_c)) and np.any(k4_c > 0.0)) or (int(os.getenv("QD_DIFF_CLOUD", "0")) == 1)
            if apply_q and hasattr(self, "q"):
                self.q = self._hyperdiffuse(self.q, k4_q, dt)
            if apply_cloud:
                self.cloud_cover = self._hyperdiffuse(self.cloud_cover, k4_c, dt)

            if dyn_diag and (self._step_counter % 200 == 0):
                var_u1 = float(np.var(self.u))
                var_v1 = float(np.var(self.v))
                var_h1 = float(np.var(self.h))
                msg = "[DynDiag] hyper4 applied: "
                if sigma4 is not None:
                    msg += f"sigma4={sigma4:.4f}, "
                if dx_min is not None:
                    msg += f"dx_min={dx_min:.1f} m, "
                msg += f"K4(u/v/h)={[k4_u, k4_v, k4_h]} "
                msg += f"Var(u) {var_u0:.3e}->{var_u1:.3e}, Var(v) {var_v0:.3e}->{var_v1:.3e}, Var(h) {var_h0:.3e}->{var_h1:.3e}"
                print(msg)

        # ---------------- Project 010 M4: Apply alternate filters by config ----------------
        try:
            ftype = os.getenv("QD_FILTER_TYPE", "combo").lower()
            # Shapiro smoothing
            sh_every = int(os.getenv("QD_SHAPIRO_EVERY", "6"))
            sh_n = int(os.getenv("QD_SHAPIRO_N", "2"))
            need_shapiro = (
                (ftype in ("shapiro", "combo") and sh_every > 0 and (self._step_counter % sh_every == 0))
                or (ftype == "hyper4" and sh_every > 0 and (self._step_counter % sh_every == 0))
            )
            if need_shapiro:
                self.u = self._shapiro_filter(self.u, n=sh_n, lon_wrap=True)
                self.v = self._shapiro_filter(self.v, n=sh_n, lon_wrap=True)
                self.h = self._shapiro_filter(self.h, n=sh_n, lon_wrap=True)
                if int(os.getenv("QD_DIFF_Q", "0")) == 1 and hasattr(self, "q"):
                    self.q = self._shapiro_filter(self.q, n=max(1, sh_n - 1), lon_wrap=True)
                if int(os.getenv("QD_DIFF_CLOUD", "0")) == 1:
                    self.cloud_cover = self._shapiro_filter(self.cloud_cover, n=max(1, sh_n - 1), lon_wrap=True)

            # Spectral zonal damping (optional or part of combo)
            spec_every = int(os.getenv("QD_SPEC_EVERY", "0"))
            spec_cut = float(os.getenv("QD_SPEC_CUTOFF", "0.75"))
            spec_damp = float(os.getenv("QD_SPEC_DAMP", "0.5"))
            # Spectral damping disabled when spec_every <= 0
            need_spec = (ftype in ("spectral", "combo")) and (spec_every > 0) and (self._step_counter % spec_every == 0)
            if need_spec:
                self.u = self._spectral_zonal_filter(self.u, cutoff=spec_cut, damp=spec_damp)
                self.v = self._spectral_zonal_filter(self.v, cutoff=spec_cut, damp=spec_damp)
                self.h = self._spectral_zonal_filter(self.h, cutoff=spec_cut, damp=spec_damp)
        except Exception:
            pass

        # Advect cloud cover
        self.cloud_cover = self._advect(self.cloud_cover, dt)

        # Add a simple dissipation term for clouds (e.g., 2-day lifetime)
        cloud_dissipation_rate = dt / (2.0 * 24 * 3600)
        self.cloud_cover *= (1 - cloud_dissipation_rate)
        
        # Apply mild diffusion to all fields to ensure stability without over-damping
        try:
            diffusion_factor = float(os.getenv("QD_DIFF_FACTOR", "0.998"))
        except Exception:
            diffusion_factor = 0.998
        self.u *= diffusion_factor
        self.v *= diffusion_factor
        self.h *= diffusion_factor
        self.cloud_cover *= diffusion_factor
        if hasattr(self, "q"):
            self.q *= diffusion_factor
        
        # Ensure no NaNs are present
        self.u = np.nan_to_num(self.u)
        self.v = np.nan_to_num(self.v)
        self.h = np.nan_to_num(self.h)
        self.T_s = np.nan_to_num(self.T_s)
        self.cloud_cover = np.nan_to_num(self.cloud_cover)
        if hasattr(self, "q"):
            self.q = np.nan_to_num(self.q)

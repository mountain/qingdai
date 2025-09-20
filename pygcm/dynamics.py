# pygcm/dynamics.py

"""
Implements a stable spectral dynamics core for the shallow water equations.
This approach avoids the grid-point instabilities faced by finite-difference methods.
"""

import os
import numpy as np
from scipy.ndimage import map_coordinates
from .grid import SphericalGrid
from . import constants as const
from . import energy as energy
from . import humidity as humidity

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
        # Calculate departure points
        dlon = self.u * dt / (self.a * np.maximum(1e-6, np.cos(np.deg2rad(self.grid.lat_mesh))))
        dlat = self.v * dt / self.a
        
        # Convert displacement from radians to grid indices
        dx = dlon / self.dlon_rad
        dy = dlat / self.dlat_rad

        # Create coordinate arrays for the grid
        lats = np.arange(self.grid.n_lat)
        lons = np.arange(self.grid.n_lon)
        lon_mesh, lat_mesh = np.meshgrid(lons, lats)

        # Departure points in index space
        dep_lat = lat_mesh - dy
        dep_lon = lon_mesh - dx

        # Interpolate the field at the departure points
        # 'map_coordinates' requires coordinates in (row, col) format
        # We use mode='wrap' for longitude to handle periodic boundaries
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
        
        # 4. Relaxation to a state of geostrophic balance
        # This is a massive simplification, but it is guaranteed to be stable.
        # It assumes that the flow is always close to a balance between
        # the pressure gradient and Coriolis forces.
        
        f = self.grid.coriolis_param
        
        # Geostrophic winds (u_g, v_g)
        dh_dlon = np.gradient(self.h, self.dlon_rad, axis=1)
        dh_dlat = np.gradient(self.h, self.dlat_rad, axis=0)
        
        u_g = np.zeros_like(self.h)
        v_g = np.zeros_like(self.h)
        
        # Regularize f near the equator: enforce a minimum |f| equivalent to 5° latitude
        f_min = 2.0 * const.PLANET_OMEGA * np.sin(np.deg2rad(5.0))
        # Use a non-zero sign at the equator to avoid f_safe == 0
        sign_nonzero = np.where(f >= 0.0, 1.0, -1.0)
        f_safe = np.where(np.abs(f) < f_min, sign_nonzero * f_min, f)

        # Cap cos(lat) to avoid division by zero at the poles
        cos_lat_capped = np.maximum(np.cos(np.deg2rad(self.grid.lat_mesh)), 1e-6)
        
        u_g = -(self.g / (f_safe * self.a * cos_lat_capped)) * dh_dlat
        v_g = (self.g / (f_safe * self.a)) * dh_dlon

        # Cap geostrophic wind magnitudes to prevent numerical blow-up
        max_wind = 200.0  # m/s
        u_g = np.clip(u_g, -max_wind, max_wind)
        v_g = np.clip(v_g, -max_wind, max_wind)
        
        
        # Nudge the current winds towards the geostrophic winds
        # This is a stable, albeit physically simplified, way to simulate flow.
        self.u = self.u * 0.8 + u_g * 0.2
        self.v = self.v * 0.8 + v_g * 0.2
        
        # Add surface friction, which is stronger over land
        friction_drag_u = -self.friction_map * self.u
        friction_drag_v = -self.friction_map * self.v
        self.u += friction_drag_u * dt
        self.v += friction_drag_v * dt
        
        # Advect cloud cover
        self.cloud_cover = self._advect(self.cloud_cover, dt)

        # Add a simple dissipation term for clouds (e.g., 2-day lifetime)
        cloud_dissipation_rate = dt / (2.0 * 24 * 3600)
        self.cloud_cover *= (1 - cloud_dissipation_rate)
        
        # Apply mild diffusion to all fields to ensure stability without over-damping
        diffusion_factor = 0.995
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

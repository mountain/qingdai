"""
routing.py

Project 014: Surface Hydrology & Runoff Routing

This module provides an offline-network-driven, asynchronous runoff routing with optional
lake handling. It routes the land runoff R (kg m^-2 s^-1) produced by the hydrology
closure (P009) along a precomputed D8 flow network at a coarser "hydrology step" cadence.

Key features:
- Accumulate land runoff mass between routing steps
- On each routing step, route mass along flow_to_index in flow_order (topological order)
- Lakes: optionally pass-through to a designated outlet (storage minimal in v1)
- Diagnostics: flow accumulation map (kg s^-1), total ocean inflow rate, mass closure error

Inputs (from data/hydrology_network.nc; see projects/014 for spec):
- land_mask(lat,lon): 1=land, 0=ocean
- flow_to_index(lat,lon): linear index of downstream cell; -1 = ocean sink
- flow_order(n_land): linear indices of land cells in upstream->downstream order
- lake_mask(lat,lon) (optional)
- lake_id(lat,lon) (optional; 1..n_lakes; 0=non-lake)
- lake_outlet_i(n_lakes), lake_outlet_j(n_lakes) OR lake_outlet_index(n_lakes) (optional)

Units:
- Fluxes R/E/P are kg m^-2 s^-1
- Mass buffers are kg
- flow_accum_kgps is kg s^-1 (event mass / hydrology step seconds)

Notes:
- v1 implements "pass-through lakes": inflow to lake cells is immediately sent to the
  lake's outlet (if provided). Future versions may include explicit storage with overflow.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional, Dict

import numpy as np

try:
    from netCDF4 import Dataset
except Exception:  # pragma: no cover
    Dataset = None  # Will raise at runtime with helpful message

from . import constants


def _assert_has_netcdf():
    if Dataset is None:
        raise RuntimeError("netCDF4 is required to use RiverRouting. Install via pip install netCDF4")


def _ravel_index(i: int, j: int, n_lon: int) -> int:
    # row-major: idx = j * n_lon + i
    return j * n_lon + i


def _unravel_index(idx: int, n_lon: int) -> tuple[int, int]:
    j = idx // n_lon
    i = idx % n_lon
    return i, j


@dataclass
class RoutingDiagnostics:
    flow_accum_kgps: np.ndarray  # (n_lat, n_lon)
    ocean_inflow_kgps: float
    mass_closure_error_kg: float
    lake_volume_kg: Optional[np.ndarray] = None


class RiverRouting:
    """
    Offline-network-based runoff routing with asynchronous execution.
    """

    def __init__(
        self,
        grid,
        network_nc_path: str,
        dt_hydro_hours: float = 6.0,
        treat_lake_as_water: bool = True,
        alpha_lake: Optional[float] = None,
        diag: bool = True,
    ) -> None:
        _assert_has_netcdf()
        if not os.path.exists(network_nc_path):
            raise FileNotFoundError(f"Hydrology network file not found: {network_nc_path}")

        self.grid = grid
        self.dt_hydro_seconds = float(dt_hydro_hours) * 3600.0
        self.treat_lake_as_water = bool(treat_lake_as_water)
        self.alpha_lake = alpha_lake
        self.diag_enabled = bool(diag)

        # Basic shapes
        self.n_lat = int(getattr(grid, "n_lat", grid.lat_mesh.shape[0]))
        self.n_lon = int(getattr(grid, "n_lon", grid.lat_mesh.shape[1]))
        self.shape = (self.n_lat, self.n_lon)
        self.n_cells = self.n_lat * self.n_lon

        # Load network
        with Dataset(network_nc_path, "r") as ds:
            def rvar(name, default=None):
                try:
                    return np.array(ds.variables[name][:])
                except Exception:
                    return default

            self.land_mask = (rvar("land_mask") > 0).astype(np.uint8)
            if self.land_mask is None:
                raise RuntimeError("hydrology_network.nc missing 'land_mask' variable")
            # Flat land mask for quick membership tests (row-major)
            self.land_flat = (self.land_mask.flatten(order="C") == 1)

            flow_to = rvar("flow_to_index")
            if flow_to is None:
                raise RuntimeError("hydrology_network.nc missing 'flow_to_index' variable")
            if flow_to.shape != self.shape:
                raise RuntimeError(f"flow_to_index shape {flow_to.shape} != grid shape {self.shape}")
            self.flow_to_index = flow_to.astype(np.int64)

            flow_order = rvar("flow_order")
            if flow_order is None:
                # Fallback: derive a naive order (not guaranteed DAG safe)
                # Users should generate a proper network file with flow_order.
                land_idxs = np.where(self.land_mask.flatten(order="C") == 1)[0]
                self.flow_order = land_idxs.astype(np.int64)
            else:
                self.flow_order = flow_order.astype(np.int64)

            self.lake_mask = rvar("lake_mask", None)
            self.lake_id = rvar("lake_id", None)

            # Lake outlet indices (either linear or (i,j))
            lake_out_idx = rvar("lake_outlet_index", None)
            if lake_out_idx is None:
                lake_out_i = rvar("lake_outlet_i", None)
                lake_out_j = rvar("lake_outlet_j", None)
                if lake_out_i is not None and lake_out_j is not None:
                    lake_out_idx = (lake_out_j.astype(np.int64) * self.n_lon + lake_out_i.astype(np.int64))
            self.lake_outlet_index = lake_out_idx

            # Optional lake meta
            self.n_lakes = 0
            if self.lake_id is not None:
                self.n_lakes = int(np.max(self.lake_id))
            if self.n_lakes > 0 and self.lake_outlet_index is not None:
                if self.lake_outlet_index.shape[0] != self.n_lakes:
                    # Attempt to coerce length
                    self.n_lakes = min(self.n_lakes, self.lake_outlet_index.shape[0])
                    self.lake_outlet_index = self.lake_outlet_index[: self.n_lakes]

        # Precompute spherical cell areas (m^2)
        self.cell_area = self._compute_cell_areas()

        # Buffers and diagnostics
        self.buffer_kg = np.zeros(self.n_cells, dtype=np.float64)  # accumulated land mass
        self.t_accum = 0.0

        self._flow_accum_kg = np.zeros(self.n_cells, dtype=np.float64)  # event mass per cell
        self._ocean_inflow_kg = 0.0
        self._diag_cache: Optional[RoutingDiagnostics] = None

        # Lake storage (v1: not actively used in routing; placeholder for future)
        self.lake_volume_kg = None
        if self.n_lakes > 0:
            self.lake_volume_kg = np.zeros(self.n_lakes, dtype=np.float64)

        if self.diag_enabled:
            print(f"[Routing] Loaded network: land={int(self.land_mask.sum())} cells, "
                  f"n_lakes={self.n_lakes}, dt_hydro={self.dt_hydro_seconds/3600.0:.1f} h")

    def _compute_cell_areas(self) -> np.ndarray:
        """
        Compute spherical cell areas using lat/lon centers and uniform spacings.
        A = R^2 * dλ * (sin φ+ - sin φ-)
        """
        R = float(getattr(constants, "PLANET_RADIUS", 6.371e6))
        lats = np.asarray(self.grid.lat, dtype=float)  # degrees
        lons = np.asarray(self.grid.lon, dtype=float)  # degrees
        if lats.size != self.n_lat or lons.size != self.n_lon:
            # Fall back to mesh-based spacing
            dphi = np.deg2rad(abs(self.grid.lat[1] - self.grid.lat[0])) if self.n_lat > 1 else np.deg2rad(1.5)
            dlam = np.deg2rad(abs(self.grid.lon[1] - self.grid.lon[0])) if self.n_lon > 1 else np.deg2rad(1.5)
        else:
            dphi = np.deg2rad(abs(lats[1] - lats[0])) if self.n_lat > 1 else np.deg2rad(1.5)
            dlam = np.deg2rad(abs(lons[1] - lons[0])) if self.n_lon > 1 else np.deg2rad(1.5)

        # Latitude edges from centers
        phi_cent = np.deg2rad(self.grid.lat_mesh[:, 0])  # 1D lat per row
        phi_plus = np.clip(phi_cent + 0.5 * dphi, -0.5 * np.pi, 0.5 * np.pi)
        phi_minus = np.clip(phi_cent - 0.5 * dphi, -0.5 * np.pi, 0.5 * np.pi)
        band_area = (np.sin(phi_plus) - np.sin(phi_minus))  # shape (n_lat,)

        A_row = (R * R) * dlam * band_area  # area per cell in row
        A = np.repeat(A_row[:, None], self.n_lon, axis=1)
        return A

    def reset(self) -> None:
        self.buffer_kg.fill(0.0)
        self.t_accum = 0.0
        self._flow_accum_kg.fill(0.0)
        self._ocean_inflow_kg = 0.0
        if self.lake_volume_kg is not None:
            self.lake_volume_kg.fill(0.0)
        self._diag_cache = None

    def step(
        self,
        R_land_flux: np.ndarray,
        dt_seconds: float,
        precip_flux: Optional[np.ndarray] = None,
        evap_flux: Optional[np.ndarray] = None,
    ) -> None:
        """
        Accumulate and, when hydrology step is reached, route mass along the network.

        Args:
          R_land_flux: kg m^-2 s^-1 on land (same grid as model)
          dt_seconds: model time step (s)
          precip_flux: optional kg m^-2 s^-1 over all cells (for lake P)
          evap_flux: optional kg m^-2 s^-1 over all cells (for lake E)
        """
        # Accumulate runoff mass on land cells only
        R = np.asarray(R_land_flux, dtype=float)
        if R.shape != self.shape:
            raise ValueError(f"R_land_flux shape {R.shape} != grid shape {self.shape}")

        land = (self.land_mask == 1)
        mass_incr = R * self.cell_area * float(dt_seconds)
        mass_incr = np.where(land, mass_incr, 0.0)
        self.buffer_kg += mass_incr.flatten(order="C")
        self.t_accum += float(dt_seconds)

        if self.t_accum + 1e-9 < self.dt_hydro_seconds:
            return  # not time yet

        # Route once for the accumulated window (event_dt = t_accum)
        event_dt = self.t_accum
        self.t_accum = 0.0

        # Copy and clear buffer for next accumulation
        acc = self.buffer_kg.copy()
        self.buffer_kg.fill(0.0)

        # Reset diagnostics for this event
        self._flow_accum_kg.fill(0.0)
        self._ocean_inflow_kg = 0.0
        mass_input = float(np.sum(acc))

        # Precompute lake mappings
        has_lakes = self.lake_mask is not None and self.lake_id is not None and self.n_lakes > 0
        lake_cell_is_lake = self.lake_mask.flatten(order="C") > 0 if has_lakes else None
        lake_ids_flat = self.lake_id.flatten(order="C") if has_lakes else None  # 1..n_lakes or 0
        lake_outlet = self.lake_outlet_index if has_lakes else None

        # Route along topological flow order
        for idx in self.flow_order:
            m = acc[idx]
            if m <= 0.0:
                continue

            # Record mass passing through this cell (for flow accumulation map)
            self._flow_accum_kg[idx] += m

            # Lakes: pass-through to outlet; if outlet==-1, treat as ocean sink; else store if no outlet
            if has_lakes and lake_cell_is_lake[idx]:
                lid = int(lake_ids_flat[idx])
                if lid > 0 and lake_outlet is not None and lid <= lake_outlet.shape[0]:
                    outlet_idx = int(lake_outlet[lid - 1])
                    if outlet_idx < 0:
                        # Direct ocean sink for this lake
                        self._ocean_inflow_kg += m
                    elif 0 <= outlet_idx < self.n_cells and self.land_flat[outlet_idx]:
                        acc[outlet_idx] += m
                    else:
                        # Outlet points outside land domain (treat as ocean)
                        self._ocean_inflow_kg += m
                else:
                    # No outlet known -> store internally
                    if self.lake_volume_kg is not None and lid > 0:
                        self.lake_volume_kg[lid - 1] += m

                acc[idx] = 0.0
                continue

            # Normal land cell: route to downstream or ocean
            dn = int(self.flow_to_index.flatten(order="C")[idx])
            # Treat any non-land downstream (including ocean) as ocean sink
            if dn < 0 or not self.land_flat[dn]:
                self._ocean_inflow_kg += m
                acc[idx] = 0.0
            else:
                acc[dn] += m
                acc[idx] = 0.0

        # Any leftover mass in 'acc' corresponds to sinks (should be near 0 except internal lakes)
        residual_cells_mass = float(np.sum(acc))

        # Optionally update lake storage from (P - E) over lake area for the event window
        lake_delta_mass = 0.0
        if has_lakes and self.lake_volume_kg is not None and precip_flux is not None and evap_flux is not None:
            P = np.asarray(precip_flux, dtype=float)
            E = np.asarray(evap_flux, dtype=float)
            lake_mask_bool = self.lake_mask.astype(bool)
            net = (P - E) * self.cell_area * event_dt
            lake_add = float(np.sum(np.where(lake_mask_bool, net, 0.0)))
            # Distribute by lake id proportional to lake area (simplified)
            if lake_add != 0.0 and self.n_lakes > 0:
                ids = self.lake_id
                for k in range(1, self.n_lakes + 1):
                    lake_area = np.sum(np.where((ids == k), self.cell_area, 0.0))
                    frac = 0.0 if lake_area <= 0 else lake_area / np.sum(np.where(lake_mask_bool, self.cell_area, 0.0))
                    self.lake_volume_kg[k - 1] += frac * lake_add
                lake_delta_mass = lake_add

        # Closure: input ≈ ocean_out + lake_delta + residual
        mass_out = self._ocean_inflow_kg + lake_delta_mass + residual_cells_mass
        closure_err = mass_input - mass_out

        # Cache diagnostics
        flow_accum_kgps = (self._flow_accum_kg / max(event_dt, 1e-9)).reshape(self.shape, order="C")
        self._diag_cache = RoutingDiagnostics(
            flow_accum_kgps=flow_accum_kgps,
            ocean_inflow_kgps=float(self._ocean_inflow_kg / max(event_dt, 1e-9)),
            mass_closure_error_kg=float(closure_err),
            lake_volume_kg=(self.lake_volume_kg.copy() if self.lake_volume_kg is not None else None),
        )

        if self.diag_enabled:
            print(f"[HydroRouting] ocean_inflow={self._diag_cache.ocean_inflow_kgps:.3e} kg/s | "
                  f"mass_error={self._diag_cache.mass_closure_error_kg:.3e} kg")

    def diagnostics(self) -> Dict[str, object]:
        """
        Return the last routing diagnostics. If no routing event has occurred yet, returns zeros.
        """
        if self._diag_cache is None:
            zeros = np.zeros(self.shape, dtype=float)
            return {
                "flow_accum_kgps": zeros,
                "ocean_inflow_kgps": 0.0,
                "mass_closure_error_kg": 0.0,
                "lake_volume_kg": (np.zeros(self.n_lakes, dtype=float) if self.n_lakes > 0 else None),
            }
        return {
            "flow_accum_kgps": self._diag_cache.flow_accum_kgps,
            "ocean_inflow_kgps": self._diag_cache.ocean_inflow_kgps,
            "mass_closure_error_kg": self._diag_cache.mass_closure_error_kg,
            "lake_volume_kg": self._diag_cache.lake_volume_kg,
        }

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validate hydrology routing:
- Checks that land cells that flow to ocean are encoded as -1 (ocean sink), not pointing into ocean grid
- Runs RiverRouting with a synthetic constant land runoff to ensure:
  * ocean_inflow_kgps > 0
  * mass_closure_error_kg is small relative to input (|err|/input << 1)
Usage:
  python3 -m scripts.validate_hydro_routing --net data/hydrology_network.nc --hours 6 --runoff 1e-6
"""
import argparse
import os
import sys
import numpy as np
from netCDF4 import Dataset

# Ensure project root in sys.path for local imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pygcm.grid import SphericalGrid
from pygcm.routing import RiverRouting


def validate_network(path: str) -> dict:
    with Dataset(path, "r") as ds:
        land = (np.array(ds["land_mask"][:]) > 0)
        flow = np.array(ds["flow_to_index"][:], dtype=np.int64)
        nlat = ds.dimensions["lat"].size
        nlon = ds.dimensions["lon"].size

    flat_land = land.flatten(order="C")
    flat_flow = flow.flatten(order="C")

    n_land = int(flat_land.sum())
    ocean_sink_cnt = int(np.sum(flat_land & (flat_flow == -1)))
    bad_ptr_cnt = 0
    for idx, dn in enumerate(flat_flow):
        if not flat_land[idx]:
            continue
        if dn >= 0 and not flat_land[dn]:
            bad_ptr_cnt += 1

    return {
        "n_lat": nlat,
        "n_lon": nlon,
        "n_land": n_land,
        "ocean_sink_count": ocean_sink_cnt,
        "bad_ptr_count": bad_ptr_cnt,
    }


def run_routing_test(net_path: str, dt_hours: float, runoff_kg_m2_s: float) -> dict:
    # Read grid size from network and build grid
    with Dataset(net_path, "r") as ds:
        nlat = ds.dimensions["lat"].size
        nlon = ds.dimensions["lon"].size
        land_mask = (np.array(ds["land_mask"][:]) > 0)

    grid = SphericalGrid(n_lat=nlat, n_lon=nlon)

    router = RiverRouting(
        grid=grid,
        network_nc_path=net_path,
        dt_hydro_hours=dt_hours,
        treat_lake_as_water=True,
        diag=True,
    )

    # Synthetic constant land runoff for a single hydro step
    R = np.zeros((nlat, nlon), dtype=float)
    R[land_mask] = float(runoff_kg_m2_s)

    # Route over exactly one hydrology window in two model steps
    dt_total = dt_hours * 3600.0
    router.step(R_land_flux=R, dt_seconds=dt_total / 2.0)
    router.step(R_land_flux=R, dt_seconds=dt_total / 2.0)

    diag = router.diagnostics()

    # Compute expected mass input to compare closure (land area * runoff * dt_total)
    # Recompute cell areas similarly to router._compute_cell_areas (A = R^2 * dλ * (sin φ+ - sin φ-))
    from pygcm import constants as _const
    Rplanet = float(getattr(_const, "PLANET_RADIUS", 6.371e6))
    lats = np.asarray(grid.lat, dtype=float)
    lons = np.asarray(grid.lon, dtype=float)
    dphi = np.deg2rad(abs(lats[1] - lats[0])) if nlat > 1 else np.deg2rad(1.5)
    dlam = np.deg2rad(abs(lons[1] - lons[0])) if nlon > 1 else np.deg2rad(1.5)
    phi_cent = np.deg2rad(grid.lat_mesh[:, 0])
    phi_plus = np.clip(phi_cent + 0.5 * dphi, -0.5 * np.pi, 0.5 * np.pi)
    phi_minus = np.clip(phi_cent - 0.5 * dphi, -0.5 * np.pi, 0.5 * np.pi)
    band_area = (np.sin(phi_plus) - np.sin(phi_minus))
    A_row = (Rplanet * Rplanet) * dlam * band_area
    A = np.repeat(A_row[:, None], nlon, axis=1)
    mass_input = float(np.sum(np.where(land_mask, A, 0.0)) * runoff_kg_m2_s * dt_total)

    return {
        "ocean_inflow_kgps": float(diag["ocean_inflow_kgps"]),
        "mass_error_kg": float(diag["mass_closure_error_kg"]),
        "mass_input_kg": mass_input,
        "mass_error_rel": (float(diag["mass_closure_error_kg"]) / mass_input) if mass_input > 0 else np.nan,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", type=str, default="data/hydrology_network.nc", help="Hydrology network NetCDF path")
    ap.add_argument("--hours", type=float, default=6.0, help="Hydrology routing step hours")
    ap.add_argument("--runoff", type=float, default=1e-6, help="Synthetic constant land runoff (kg m^-2 s^-1)")
    args = ap.parse_args()

    if not os.path.exists(args.net):
        print(f"ERROR: Network file not found: {args.net}", file=sys.stderr)
        sys.exit(2)

    net_info = validate_network(args.net)
    print("[ValidateNet] n_lat={n_lat} n_lon={n_lon} n_land={n_land} ocean_sink={ocean_sink_count} bad_ptr={bad_ptr_count}".format(**net_info))

    res = run_routing_test(args.net, dt_hours=args.hours, runoff_kg_m2_s=args.runoff)
    print("[ValidateRoute] ocean_inflow={o:.3e} kg/s | mass_error={e:.3e} kg | mass_input={m:.3e} kg | rel_err={r:.3e}".format(
        o=res["ocean_inflow_kgps"], e=res["mass_error_kg"], m=res["mass_input_kg"], r=res["mass_error_rel"]
    ))

    # Simple pass/fail heuristic
    ok_inflow = res["ocean_inflow_kgps"] > 0.0
    ok_closure = abs(res["mass_error_rel"]) < 1e-6  # very strict; adjust if needed
    if ok_inflow and ok_closure:
        print("[OK] Routing sends mass to ocean and closes mass within tolerance.")
        sys.exit(0)
    else:
        print("[WARN] Check routing: ok_inflow=%s ok_closure=%s" % (ok_inflow, ok_closure))
        sys.exit(1)


if __name__ == "__main__":
    main()

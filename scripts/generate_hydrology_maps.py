#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
P014: Generate offline hydrology routing network (D8 flow, pit filling, lakes).

Outputs a NetCDF at --out (default: data/hydrology_network.nc) with:
- lat(n_lat), lon(n_lon)
- land_mask(n_lat, n_lon) uint8 (1=land, 0=ocean)
- elevation_filled(n_lat, n_lon) float32 (pit-filled elevation)
- flow_to_index(n_lat, n_lon) int32 (row-major linear index; -1=ocean or terminal sink)
- flow_order(n_land) int32 (row-major linear indices of land cells, upstream→downstream)
- lake_mask(n_lat, n_lon) uint8 (1=lake, else 0)
- lake_id(n_lat, n_lon) int32 (0=non-lake, 1..n_lakes for lakes)
- lake_outlet_index(n_lakes) int32 (row-major index of outlet; -1 if unknown)

Usage:
  python3 -m scripts.generate_hydrology_maps --topo data/topography_*.nc --out data/hydrology_network.nc
"""

from __future__ import annotations

import argparse
import os
import sys
from typing import Tuple, List

import numpy as np

try:
    from netCDF4 import Dataset
except Exception:  # pragma: no cover
    Dataset = None

# Ensure project root in sys.path
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pygcm.grid import SphericalGrid
from pygcm.topography import load_topography_from_netcdf, create_land_sea_mask, generate_base_properties
from pygcm import constants


def _require_netcdf():
    if Dataset is None:
        raise RuntimeError("netCDF4 is required. Install via: python3 -m pip install netCDF4")


def row_major_index(i: int, j: int, n_lon: int) -> int:
    return j * n_lon + i


def neighbors_d8(i: int, j: int, n_lon: int, n_lat: int) -> List[Tuple[int, int]]:
    # Longitude wraps periodically; latitude clamped at edges
    res = []
    for dj in (-1, 0, 1):
        for di in (-1, 0, 1):
            if di == 0 and dj == 0:
                continue
            jj = j + dj
            if jj < 0 or jj >= n_lat:
                continue
            ii = (i + di) % n_lon
            res.append((ii, jj))
    return res


def spherical_distance(grid: SphericalGrid, i1: int, j1: int, i2: int, j2: int) -> float:
    # Approximate physical distance between cell centers using local metric
    R = float(getattr(constants, "PLANET_RADIUS", 6.371e6))
    lat1 = np.deg2rad(grid.lat[j1])
    lon1 = np.deg2rad(grid.lon[i1])
    lat2 = np.deg2rad(grid.lat[j2])
    lon2 = np.deg2rad(grid.lon[i2])
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    # Wrap dlon to [-pi, pi] for periodicity
    if dlon > np.pi:
        dlon -= 2 * np.pi
    elif dlon < -np.pi:
        dlon += 2 * np.pi
    # Equirectangular approximation
    x = dlon * np.cos(0.5 * (lat1 + lat2))
    y = dlat
    return R * np.sqrt(x * x + y * y)


def pit_fill(elev: np.ndarray, land_mask: np.ndarray, max_iters: int = 200, eps: float = 1e-3) -> np.ndarray:
    """
    Simple iterative pit-filling for land cells to remove local sinks.
    For each land cell, if all D8 neighbors have elevation >= cell,
    raise cell to min(neighbor elevation) + eps. Repeat up to max_iters.
    """
    e = elev.copy()
    n_lat, n_lon = e.shape
    changed = True
    it = 0
    while changed and it < max_iters:
        changed = False
        it += 1
        for j in range(n_lat):
            for i in range(n_lon):
                if land_mask[j, i] != 1:
                    continue
                neigh_vals = [e[jj, ii] for (ii, jj) in neighbors_d8(i, j, n_lon, n_lat)]
                if not neigh_vals:
                    continue
                mn = min(neigh_vals)
                if e[j, i] <= mn:
                    new_val = mn + eps
                    if new_val > e[j, i]:
                        e[j, i] = new_val
                        changed = True
    return e


def compute_flow_to_index(grid: SphericalGrid, elev: np.ndarray, land_mask: np.ndarray) -> np.ndarray:
    """
    Compute D8 downstream indices using steepest descent with spherical distances.
    Returns array int64 of shape (n_lat, n_lon), with -1 for ocean sinks and terminal sinks.
    Note: If the steepest neighbor is ocean, encode -1 (ocean sink) instead of pointing to an ocean cell.
    """
    n_lat, n_lon = elev.shape
    flow_to = np.full((n_lat, n_lon), -1, dtype=np.int64)

    for j in range(n_lat):
        for i in range(n_lon):
            if land_mask[j, i] != 1:
                flow_to[j, i] = -1
                continue
            z0 = elev[j, i]
            best_slope = -np.inf
            best_ii = -1
            best_jj = -1
            for ii, jj in neighbors_d8(i, j, n_lon, n_lat):
                dist = spherical_distance(grid, i, j, ii, jj)
                if dist <= 0:
                    continue
                z1 = elev[jj, ii]
                slope = (z0 - z1) / dist  # positive = downhill
                if slope > best_slope:
                    best_slope = slope
                    best_ii, best_jj = ii, jj

            if best_slope > 0 and best_ii >= 0:
                # If downstream is ocean, mark as ocean sink (-1), do not point into ocean grid
                if land_mask[best_jj, best_ii] == 1:
                    flow_to[j, i] = row_major_index(best_ii, best_jj, n_lon)
                else:
                    flow_to[j, i] = -1
            else:
                flow_to[j, i] = -1
    return flow_to


def topo_sort_flow_order(flow_to: np.ndarray, land_mask: np.ndarray) -> np.ndarray:
    """
    Topological order for D8 graph defined on land cells only.
    Kahn's algorithm using indegree; nodes with flow_to==-1 have no outgoing edges.
    """
    n_lat, n_lon = land_mask.shape
    n_cells = n_lat * n_lon
    land_flat = (land_mask.flatten(order="C") == 1)
    flow_flat = flow_to.flatten(order="C")

    indeg = np.zeros(n_cells, dtype=np.int64)
    for idx in range(n_cells):
        if not land_flat[idx]:
            continue
        dn = int(flow_flat[idx])
        if dn >= 0 and land_flat[dn]:
            indeg[dn] += 1

    from collections import deque
    q = deque([idx for idx in range(n_cells) if land_flat[idx] and indeg[idx] == 0])
    order: List[int] = []

    while q:
        u = q.popleft()
        order.append(u)
        dn = int(flow_flat[u])
        if dn >= 0 and land_flat[dn]:
            indeg[dn] -= 1
            if indeg[dn] == 0:
                q.append(dn)

    if len(order) < int(land_flat.sum()):
        seen = set(order)
        rem = [idx for idx in range(n_cells) if land_flat[idx] and idx not in seen]
        order.extend(rem)

    return np.array(order, dtype=np.int64)


def identify_lakes(flow_to: np.ndarray, land_mask: np.ndarray) -> Tuple[np.ndarray, np.ndarray, int]:
    """
    Identify lake cells (terminal sinks with flow_to==-1) and label connected lakes.
    Returns (lake_mask uint8, lake_id int32, n_lakes).
    """
    n_lat, n_lon = land_mask.shape
    lake_mask = np.zeros((n_lat, n_lon), dtype=np.uint8)
    lake_id = np.zeros((n_lat, n_lon), dtype=np.int32)

    term = (land_mask == 1) & (flow_to == -1)
    if not np.any(term):
        return lake_mask, lake_id, 0

    visited = np.zeros_like(term, dtype=bool)
    lake_count = 0

    for j in range(n_lat):
        for i in range(n_lon):
            if not term[j, i] or visited[j, i]:
                continue
            lake_count += 1
            stack = [(i, j)]
            visited[j, i] = True
            while stack:
                ii, jj = stack.pop()
                lake_mask[jj, ii] = 1
                lake_id[jj, ii] = lake_count
                for ni, nj in neighbors_d8(ii, jj, n_lon, n_lat):
                    if term[nj, ni] and not visited[nj, ni]:
                        visited[nj, ni] = True
                        stack.append((ni, nj))

    return lake_mask, lake_id, lake_count


def compute_lake_outlets(grid: SphericalGrid,
                         elev_filled: np.ndarray,
                         lake_mask: np.ndarray,
                         lake_id: np.ndarray,
                         land_mask: np.ndarray) -> np.ndarray:
    """
    For each lake (connected component where lake_mask==1), determine an outlet:
      - If any lake cell touches ocean (land_mask==0) in its D8 neighborhood, mark outlet as -1 (direct ocean sink).
      - Else choose the neighboring non-lake land cell with the lowest filled elevation as the outlet (row-major index).
      - If neither is found (should not happen with proper pit-filling), set -1.
    Returns:
      lake_outlet_index: (n_lakes,) int32 array (row-major indices or -1 for direct ocean sink).
    """
    n_lat, n_lon = lake_mask.shape
    n_lakes = int(np.max(lake_id)) if lake_id is not None else 0
    out = np.full((max(n_lakes, 0),), -1, dtype=np.int32)
    if n_lakes == 0:
        return out

    for k in range(1, n_lakes + 1):
        best_idx = -1
        best_z = np.inf
        touches_ocean = False

        # iterate over cells of this lake
        js, is_ = np.where(lake_id == k)
        for j, i in zip(js, is_):
            for ii, jj in neighbors_d8(i, j, n_lon, n_lat):
                if lake_mask[jj, ii] == 1:
                    continue  # still lake
                if land_mask[jj, ii] == 0:
                    touches_ocean = True
                    break
                # neighbor is land and not lake: candidate outlet via non-lake land cell
                z = float(elev_filled[jj, ii])
                if z < best_z:
                    best_z = z
                    best_idx = row_major_index(ii, jj, n_lon)
            if touches_ocean:
                break

        if touches_ocean:
            out[k - 1] = -1  # direct ocean sink
        else:
            out[k - 1] = int(best_idx) if best_idx >= 0 else -1

    return out


def main():
    _require_netcdf()

    parser = argparse.ArgumentParser(description="Generate P014 hydrology routing network NetCDF.")
    parser.add_argument("--topo", type=str, default=os.getenv("QD_TOPO_NC", ""),
                        help="Path to topography NetCDF (contains land_mask and optional elevation). If empty, fallback is used.")
    parser.add_argument("--out", type=str, default="data/hydrology_network.nc",
                        help="Output NetCDF path")
    parser.add_argument("--nlat", type=int, default=121, help="Grid latitude count if fallback topography is used")
    parser.add_argument("--nlon", type=int, default=240, help="Grid longitude count if fallback topography is used")
    parser.add_argument("--pit-eps", type=float, default=1e-3, help="Pit filling epsilon")
    parser.add_argument("--pit-iters", type=int, default=200, help="Max iterations for pit filling")
    args = parser.parse_args()

    topo_path = args.topo
    out_path = args.out

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)

    grid = SphericalGrid(n_lat=args.nlat, n_lon=args.nlon)

    elevation = None
    if topo_path and os.path.exists(topo_path):
        try:
            elevation, land_mask, _base_albedo, _friction = load_topography_from_netcdf(topo_path, grid)
            print(f"[HydroNet] Loaded topography from '{topo_path}'.")
        except Exception as e:
            print(f"[HydroNet] Failed to load '{topo_path}': {e}\nFalling back to procedural mask.")
            land_mask = create_land_sea_mask(grid)
            _base_albedo, _friction = generate_base_properties(land_mask)
            elevation = np.zeros_like(grid.lat_mesh, dtype=float)
    else:
        print(f"[HydroNet] No topography specified or file missing. Using fallback.")
        land_mask = create_land_sea_mask(grid)
        _base_albedo, _friction = generate_base_properties(land_mask)
        elevation = np.zeros_like(grid.lat_mesh, dtype=float)

    land_mask = land_mask.astype(np.uint8)
    elevation = elevation.astype(float)

    print(f"[HydroNet] Pit filling elevation over land (iters={args.pit_iters}, eps={args.pit_eps})...")
    elev_filled = pit_fill(elevation.copy(), land_mask, max_iters=args.pit_iters, eps=args.pit_eps)

    print("[HydroNet] Computing D8 flow directions...")
    flow_to = compute_flow_to_index(grid, elev_filled, land_mask)

    print("[HydroNet] Identifying lakes...")
    lake_mask, lake_id, n_lakes = identify_lakes(flow_to, land_mask)
    lake_outlet_index = compute_lake_outlets(grid, elev_filled, lake_mask, lake_id, land_mask) if n_lakes > 0 else None

    print("[HydroNet] Computing topological flow order...")
    flow_order = topo_sort_flow_order(flow_to, land_mask)

    print(f"[HydroNet] Writing network to '{out_path}'...")
    with Dataset(out_path, "w") as ds:
        nlat, nlon = grid.n_lat, grid.n_lon
        ds.createDimension("lat", nlat)
        ds.createDimension("lon", nlon)
        # For land list (flow_order)
        n_land = int((land_mask == 1).sum())
        ds.createDimension("n_land", n_land)
        if n_lakes > 0:
            ds.createDimension("n_lakes", int(n_lakes))

        vlat = ds.createVariable("lat", "f4", ("lat",))
        vlon = ds.createVariable("lon", "f4", ("lon",))
        vlat[:] = grid.lat.astype(np.float32)
        vlon[:] = grid.lon.astype(np.float32)

        def wvar(name, dtype, dims, data):
            var = ds.createVariable(name, dtype, dims)
            var[:] = data

        wvar("land_mask", "u1", ("lat", "lon"), land_mask.astype(np.uint8))
        wvar("elevation_filled", "f4", ("lat", "lon"), elev_filled.astype(np.float32))
        wvar("flow_to_index", "i4", ("lat", "lon"), flow_to.astype(np.int32))
        wvar("flow_order", "i4", ("n_land",), flow_order.astype(np.int32))
        wvar("lake_mask", "u1", ("lat", "lon"), lake_mask.astype(np.uint8))
        wvar("lake_id", "i4", ("lat", "lon"), lake_id.astype(np.int32))
        if n_lakes > 0:
            wvar("lake_outlet_index", "i4", ("n_lakes",), lake_outlet_index.astype(np.int32))

        ds.setncattr("title", "Qingdai Hydrology Network")
        ds.setncattr("indexing", "row-major (i=lon index, j=lat index), idx=j*n_lon+i")
        ds.setncattr("projection", "latlon")
        ds.setncattr("created_by", "scripts/generate_hydrology_maps.py")
        ds.setncattr("notes", "D8 routing; simple pit filling; lakes are terminal sinks; outlets not discovered in v1")

    print("[HydroNet] Done.")


if __name__ == "__main__":
    main()

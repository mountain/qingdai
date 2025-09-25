from __future__ import annotations

import os
import math
from typing import Tuple, Optional

import numpy as np

try:
    import matplotlib.pyplot as plt  # type: ignore
except Exception as e:
    plt = None


def _require_matplotlib():
    if plt is None:
        raise RuntimeError("matplotlib is required for plotting. Please install 'matplotlib'.")


def _nearest_ij(grid, lat_deg: float, lon_deg: float) -> Tuple[int, int]:
    """
    Find nearest grid index (j,i) to given lat/lon in degrees.
    Handles periodic longitude.
    """
    j = int(np.argmin(np.abs(np.asarray(grid.lat) - float(lat_deg))))
    lon = np.asarray(grid.lon)
    L0, L1 = float(lon.min()), float(lon.max())
    x = float(lon_deg)
    # normalize to [L0, L1] assuming 360-degree periodicity
    span = 360.0
    if L1 - L0 >= 300.0:  # crude check for full longitudes
        while x < L0:
            x += span
        while x > L1:
            x -= span
    i = int(np.argmin(np.abs(lon - x)))
    return j, i


def _cell_area_m2(grid, j: int, R: float = 6_371_000.0) -> float:
    """
    Regular lat-lon cell area at latitude index j. Uses:
      A = R^2 * dlon * (sin(phi+0.5*dphi) - sin(phi-0.5*dphi))
    """
    lat = np.deg2rad(np.asarray(grid.lat))
    lon = np.deg2rad(np.asarray(grid.lon))
    if lat.size < 2 or lon.size < 2:
        return math.pi * R * R  # fallback: arbitrary non-zero
    dlat = float(np.abs(lat[1] - lat[0]))
    dlon = float(np.abs(lon[1] - lon[0]))
    phi = float(lat[j])
    phi_p = phi + 0.5 * dlat if j < lat.size - 1 else phi
    phi_m = phi - 0.5 * dlat if j > 0 else phi
    return R * R * dlon * (math.sin(phi_p) - math.sin(phi_m))


def _species_height_map(pop, s: int, H_scale: float = 10.0) -> Optional[np.ndarray]:
    """
    Species-resolved canopy height map:
      H_s = H_scale * sum_k (h_k * LAI_s,k) / sum_k LAI_s,k
    """
    LAI_layers_SK = getattr(pop, "LAI_layers_SK", None)
    if LAI_layers_SK is None:
        return None
    S, K, H, W = LAI_layers_SK.shape
    s = int(np.clip(s, 0, S - 1))
    LAI_s_k = np.maximum(LAI_layers_SK[s], 0.0)  # [K,H,W]
    if K <= 0:
        return None
    # Use k = 0..K-1 → h_k = (k+1)/K to match LAI_layers_SK's leading K dimension
    idx = (np.arange(K, dtype=float) + 1.0) / float(K)  # [K]
    num = np.tensordot(idx, LAI_s_k, axes=(0, 0))          # [H,W]
    den = np.sum(LAI_s_k, axis=0) + 1e-12                  # [H,W]
    return float(H_scale) * (num / den)


def _top3_species_at_cell(pop, j: int, i: int) -> np.ndarray:
    """
    Return species indices (size 3) of top-3 species at (j,i) by Σ_k LAI_s,k.
    """
    LAI_layers_SK = getattr(pop, "LAI_layers_SK", None)
    if LAI_layers_SK is None:
        # fallback to species_weights * total_LAI
        w = np.asarray(getattr(pop, "species_weights", np.ones((1,), dtype=float)), dtype=float)
        total_lai = np.asarray(
            getattr(
                pop,
                "LAI",
                pop.total_LAI() if hasattr(pop, "total_LAI") else np.zeros_like(
                    getattr(pop, "species_weights", np.ones((1,)))
                )
            )
        )
        total_lai = total_lai if isinstance(total_lai, np.ndarray) else np.zeros((1,), dtype=float)
        S = int(w.size)
        lai_s = w * float(total_lai[j, i]) if total_lai.ndim == 2 else w * float(total_lai)
    else:
        S = LAI_layers_SK.shape[0]
        lai_s = np.sum(np.maximum(LAI_layers_SK[:, :, j, i], 0.0), axis=1)  # [S]
    order = np.argsort(lai_s)[::-1]
    if order.size >= 3:
        return order[:3]
    # pad if less than 3 species
    pad = [order[-1]] * (3 - order.size) if order.size > 0 else [0, 0, 0]
    return np.concatenate([order, np.asarray(pad, dtype=int)], axis=0)[:3]


def _gather_cell_individuals(indiv, cell_index: int):
    """
    Vectorized gather of individual indices for a given sampled cell (index within pool).
    """
    return np.where(indiv.indiv_cell_index == int(cell_index))[0]


def _nearest_sampled_cell(indiv, j: int, i: int) -> int:
    """
    Find the nearest sampled cell index (in pool) to a given (j,i).
    """
    # L1 distance in index space
    di = np.abs(indiv.sample_i - int(i)) + np.abs(indiv.sample_j - int(j))
    return int(np.argmin(di))


def plot_top3_species_distributions(eco,
                                    grid,
                                    indiv=None,
                                    *,
                                    lat_deg: float,
                                    lon_deg: float,
                                    nbins: int = 24,
                                    neigh_radius: int = 1,
                                    save_path: Optional[str] = None,
                                    title: Optional[str] = None):
    """
    Plot a 3x4 panel (12 subplots) for the top-3 species (by LAI at target cell):
      per species (row):
        - Column 1: Canopy height distribution (in a (2*neigh_radius+1)^2 neighborhood)
        - Column 2: Leaf area increment distribution (samples) [m^2/day]
        - Column 3: Root development proxy distribution (alloc_root · E_day) [J-equiv]
        - Column 4: Lifespan (days) shown as a vertical line

    Args:
      eco: EcologyAdapter (must have .pop and .genes_list)
      grid: SphericalGrid (with 1D lat/lon, 2D lat_mesh)
      indiv: IndividualPool or None (sample-based distributions need it)
      lat_deg/lon_deg: target location for the query
      nbins: histogram bins
      neigh_radius: neighborhood radius for canopy height distribution (1 -> 3x3)
      save_path: optional path to save the figure
      title: optional figure title override

    Returns:
      (fig, axes) matplotlib Figure and Axes for further customization.
    """
    _require_matplotlib()
    pop = getattr(eco, "pop", None)
    if pop is None:
        raise RuntimeError("EcologyAdapter.pop (PopulationManager) is required.")

    H, W = np.asarray(grid.lat_mesh).shape
    j0, i0 = _nearest_ij(grid, lat_deg, lon_deg)
    genes_list = getattr(eco, "genes_list", [])
    if not genes_list or len(genes_list) <= 0:
        raise RuntimeError("EcologyAdapter.genes_list is required with species gene parameters.")

    # Top-3 species by LAI at cell
    top3 = _top3_species_at_cell(pop, j0, i0)

    # Neighborhood indices for canopy height distribution
    jj = np.clip(np.arange(j0 - neigh_radius, j0 + neigh_radius + 1), 0, H - 1)
    ii = np.clip(np.arange(i0 - neigh_radius, i0 + neigh_radius + 1), 0, W - 1)
    JJ, II = np.meshgrid(jj, ii, indexing="ij")  # [Nh, Nw]

    # Sample-based prefetch (if available)
    has_samples = (indiv is not None) and hasattr(indiv, "indiv_cell_index") and hasattr(indiv, "indiv_species_id") and hasattr(indiv, "indiv_E_day")
    if has_samples:
        cidx = _nearest_sampled_cell(indiv, j0, i0)
        cell_mask = (indiv.indiv_cell_index == cidx)
        sp_id = indiv.indiv_species_id
        E_day = indiv.indiv_E_day  # note: reset to 0 at day-end
    else:
        cell_mask = None
        sp_id = None
        E_day = None

    # Species-level gene vectors (alloc_root, leaf_area_per_energy, lifespan_days)
    S_guess = len(genes_list)
    alloc_root_vec = np.array([float(getattr(g, "alloc_root", 0.3)) for g in genes_list], dtype=float)
    leaf_per_E_vec = np.array([float(getattr(g, "leaf_area_per_energy", 1.0e-6)) for g in genes_list], dtype=float)
    lifespan_vec = np.array([float(getattr(g, "lifespan_days", 365.0)) for g in genes_list], dtype=float)

    # Cell area (for absolute leaf area if needed)
    A_cell = _cell_area_m2(grid, j0)

    fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(16, 10), constrained_layout=True)
    for row, s in enumerate(top3):
        s = int(s)

        # 1) Canopy height distribution (neighborhood) using species-resolved height map
        try:
            Hs_map = _species_height_map(pop, s, H_scale=float(os.getenv("QD_ECO_HEIGHT_SCALE_M", "10.0")))
        except Exception:
            Hs_map = _species_height_map(pop, s, H_scale=10.0)
        H_neigh = Hs_map[JJ, II].ravel() if Hs_map is not None else np.array([], dtype=float)
        ax = axes[row, 0]
        if H_neigh.size > 0:
            ax.hist(H_neigh[~np.isnan(H_neigh)], bins=nbins, color="#6699cc")
        ax.set_title(f"Species {s}: Height (m)")
        ax.set_xlabel("m")
        ax.set_ylabel("count")

        # Sample-based per-individual distributions (leaf/roots)
        if has_samples:
            sel = cell_mask & (sp_id == s)
            E_s = E_day[sel]
            dA = E_s * leaf_per_E_vec[s]  # per-individual leaf area increment (m^2/day)
            root_proxy = alloc_root_vec[s] * E_s  # root energy proxy (J-equiv proxy)
        else:
            dA = np.array([], dtype=float)
            root_proxy = np.array([], dtype=float)

        # 2) Leaf area increment distribution (samples)
        ax = axes[row, 1]
        if dA.size > 0:
            ax.hist(dA, bins=nbins, color="#55aa55")
        ax.set_title(f"Species {s}: ΔLeaf area (m²/day, indiv)")
        ax.set_xlabel("m²/day")
        ax.set_ylabel("count")

        # 3) Root development proxy distribution (samples)
        ax = axes[row, 2]
        if root_proxy.size > 0:
            ax.hist(root_proxy, bins=nbins, color="#cc8866")
        ax.set_title(f"Species {s}: Root proxy (alloc_root·E)")
        ax.set_xlabel("J-equiv")
        ax.set_ylabel("count")

        # 4) Lifespan (vertical line)
        ax = axes[row, 3]
        val = float(lifespan_vec[s] if s < lifespan_vec.size else 365.0)
        ax.axvline(val, color="k", lw=2)
        # show +-20% window to make the line visible with some context
        ax.set_xlim(val * 0.8, val * 1.2)
        ax.set_title(f"Species {s}: Lifespan (days)")
        ax.set_xlabel("days")
        ax.set_yticks([])

    if title:
        fig.suptitle(title, fontsize=14)
    else:
        fig.suptitle(f"Top-3 species distributions at ({lat_deg:.2f}°, {lon_deg:.2f}°)", fontsize=14)

    if save_path is not None:
        try:
            os.makedirs(os.path.dirname(save_path) or ".", exist_ok=True)
        except Exception:
            pass
        fig.savefig(save_path, dpi=140)

    return fig, axes


__all__ = [
    "plot_top3_species_distributions",
]

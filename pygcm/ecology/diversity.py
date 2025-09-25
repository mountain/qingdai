from __future__ import annotations

import os
import numpy as np
import matplotlib.pyplot as plt


def _get_species_lai_SK(eco) -> tuple[np.ndarray | None, int, int]:
    """
    Return (L_s, H, W) where L_s is [S, H, W] = sum_k max(LAI_layers_SK, 0) per species.
    If unavailable, returns (None, 0, 0).
    """
    try:
        pop = getattr(eco, "pop", None)
        if pop is None:
            return None, 0, 0
        LAI_layers_SK = getattr(pop, "LAI_layers_SK", None)
        if LAI_layers_SK is None:
            return None, 0, 0
        # [S, K, H, W] -> [S, H, W]
        L_s = np.sum(np.maximum(LAI_layers_SK, 0.0), axis=1)
        H, W = L_s.shape[-2], L_s.shape[-1]
        return L_s, H, W
    except Exception:
        return None, 0, 0


def _area_weights(lat_mesh: np.ndarray) -> np.ndarray:
    """cos(lat) weights, clipped to non-negative."""
    w = np.maximum(np.cos(np.deg2rad(lat_mesh)), 0.0)
    return w


def compute_alpha_eff_map(L_s: np.ndarray, land_mask: np.ndarray) -> np.ndarray:
    """
    Compute per-pixel alpha diversity effective number of species:
      alpha_eff = exp(H), H = -sum_i p_i ln p_i, p_i = L_i / sum_j L_j
    Returns a map with NaN over ocean.
    """
    S, H, W = L_s.shape
    L_tot = np.sum(L_s, axis=0)  # [H, W]
    alpha_map = np.full((H, W), np.nan, dtype=float)
    land = (land_mask == 1)
    mask = land & (L_tot > 0)
    if np.any(mask):
        # Build p over masked indices efficiently
        P = np.zeros((S, np.count_nonzero(mask)), dtype=float)
        idx_flat = np.flatnonzero(mask)
        # For each species, gather abundances at masked cells
        for s in range(S):
            Ps = L_s[s].ravel()[idx_flat]
            P[s, :] = Ps
        denom = np.sum(P, axis=0) + 1e-15
        P = P / denom[None, :]
        H_shannon = -np.sum(P * np.log(P + 1e-15), axis=0)
        alpha_vals = np.exp(H_shannon)
        alpha_map.ravel()[idx_flat] = alpha_vals
    return alpha_map


def compute_whittaker_beta(L_s: np.ndarray, land_mask: np.ndarray, lat_mesh: np.ndarray) -> dict:
    """
    Compute Whittaker beta = gamma_eff / alpha_mean.
    - alpha_mean: area-weighted mean of per-pixel alpha_eff.
    - gamma_eff: exp(H_gamma) where p_gamma is area-weighted global composition over land.
    Returns dict with alpha_mean, gamma_eff, beta_whittaker.
    """
    alpha_map = compute_alpha_eff_map(L_s, land_mask)
    land = (land_mask == 1)
    w = _area_weights(lat_mesh)
    w_sum_land = float(np.sum(w[land])) + 1e-15
    w_norm = w / w_sum_land

    # alpha_mean
    alpha_mean = float(np.nansum(alpha_map[land] * w_norm[land]))

    # gamma from area-weighted species totals
    S = L_s.shape[0]
    T_s = np.zeros((S,), dtype=float)
    for s in range(S):
        T_s[s] = float(np.nansum(L_s[s][land] * w_norm[land]))
    T_sum = float(np.sum(T_s)) + 1e-15
    p_gamma = T_s / T_sum
    H_gamma = float(-np.sum(p_gamma * np.log(p_gamma + 1e-15)))
    gamma_eff = float(np.exp(H_gamma))

    beta_w = float(gamma_eff / max(alpha_mean, 1e-12))
    return {"alpha_mean": alpha_mean, "gamma_eff": gamma_eff, "beta_whittaker": beta_w, "alpha_map": alpha_map}


def compute_local_bray_curtis(L_s: np.ndarray, land_mask: np.ndarray) -> np.ndarray:
    """
    Compute a local beta diversity map using mean Bray–Curtis dissimilarity to 4-neighbors.
    For each land cell c with abundance vector a (S,), and each neighbor n with b:
      BC(a,b) = 1 - 2 * sum_i min(a_i, b_i) / (sum_i a_i + sum_i b_i)
    Returns map with NaN over ocean.
    """
    S, H, W = L_s.shape
    land = (land_mask == 1)
    # Precompute sums per cell
    sum_a = np.sum(L_s, axis=0)  # [H, W]
    # Prepare neighbor shifts: up, down, left, right (with periodic longitude)
    shifts = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    bc_accum = np.zeros((H, W), dtype=float)
    count = np.zeros((H, W), dtype=float)

    # For performance, compute BC for each shift vectorized
    for dj, di in shifts:
        # neighbor indices with periodic wrap in longitude, clipped in latitude
        j_src = np.arange(H)
        i_src = np.arange(W)
        j_nbr = np.clip(j_src[:, None] + dj, 0, H - 1)  # [H,1]
        i_nbr = (i_src[None, :] + di) % W               # [1,W]

        # Broadcast neighbor fields [S,H,W] -> [S,H,W] with shifted indices
        L_n = L_s[:, j_nbr, i_nbr]  # [S,H,W]
        sum_b = np.sum(L_n, axis=0)  # [H,W]

        # Numerator and denominator
        min_sum = np.sum(np.minimum(L_s, L_n), axis=0)  # [H,W]
        denom = (sum_a + sum_b) + 1e-15
        bc = 1.0 - 2.0 * (min_sum / denom)

        # Only count where both center and neighbor are land
        nbr_land = land[j_nbr, i_nbr]
        valid = land & nbr_land
        bc_accum[valid] += bc[valid]
        count[valid] += 1.0

    # Mean over available neighbors
    with np.errstate(invalid="ignore", divide="ignore"):
        bc_mean = np.where(count > 0, bc_accum / count, np.nan)
    # Keep NaN over ocean
    bc_mean[~land] = np.nan
    return bc_mean


def plot_diversity_maps(grid, land_mask, t_days: float, outdir: str,
                        alpha_map: np.ndarray,
                        bc_local: np.ndarray,
                        beta_whittaker: float,
                        alpha_mean: float,
                        gamma_eff: float) -> None:
    os.makedirs(outdir, exist_ok=True)
    # 1) Alpha (effective species) map
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.5), constrained_layout=True)
    levels = np.linspace(0, np.nanmax(alpha_map[land_mask == 1]) if np.any(np.isfinite(alpha_map)) else 1.0, 16)
    cs = ax.contourf(grid.lon, grid.lat, alpha_map, levels=levels, cmap="viridis")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.6)
    ax.set_title(f"Alpha diversity (effective S) — Day {t_days:.2f}")
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    fig.colorbar(cs, ax=ax, label="Effective species number")
    plt.savefig(os.path.join(outdir, f"alpha_effective_day_{t_days:05.1f}.png"), dpi=140)
    plt.close(fig)

    # 2) Local beta (Bray–Curtis to neighbors)
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.5), constrained_layout=True)
    levels = np.linspace(0, 1, 21)
    cs = ax.contourf(grid.lon, grid.lat, bc_local, levels=levels, cmap="magma")
    ax.contour(grid.lon, grid.lat, land_mask, levels=[0.5], colors="black", linewidths=0.6)
    ax.set_title(f"Local beta (Bray–Curtis to neighbors) — Day {t_days:.2f}")
    ax.set_xlabel("Longitude"); ax.set_ylabel("Latitude")
    fig.colorbar(cs, ax=ax, label="Bray–Curtis (0..1)")
    plt.savefig(os.path.join(outdir, f"beta_local_braycurtis_day_{t_days:05.1f}.png"), dpi=140)
    plt.close(fig)

    # 3) Text summary
    with open(os.path.join(outdir, f"diversity_summary_day_{t_days:05.1f}.txt"), "w", encoding="utf-8") as f:
        f.write(f"Day: {t_days:.2f}\n")
        f.write(f"Whittaker beta (β = γ/ᾱ): {beta_whittaker:.4f}\n")
        f.write(f"  alpha_mean (ᾱ): {alpha_mean:.4f}\n")
        f.write(f"  gamma_eff  (γ ): {gamma_eff:.4f}\n")


def save_community_npz(L_s: np.ndarray, land_mask: np.ndarray, t_days: float, outdir: str) -> None:
    """
    Save the community species-level abundance tensor (per cell) for offline analysis.
    File: community_day_{t}.npz with arrays: L_s [S,H,W], land_mask [H,W]
    """
    try:
        os.makedirs(outdir, exist_ok=True)
        path = os.path.join(outdir, f"community_day_{t_days:05.1f}.npz")
        np.savez(path, L_s=L_s.astype(np.float32), land_mask=land_mask.astype(np.int8))
    except Exception:
        pass


def compute_and_plot(grid, eco, land_mask: np.ndarray, t_days: float, base_output_dir: str) -> None:
    """
    Entry point: compute alpha/beta diversity metrics from current ecology population state
    and write plots/files under {base_output_dir}/ecology/.
    """
    try:
        L_s, H, W = _get_species_lai_SK(eco)
        if L_s is None or H == 0 or W == 0:
            return
        outdir = os.path.join(base_output_dir, "ecology")
        # Alpha map
        alpha_map = compute_alpha_eff_map(L_s, land_mask)
        # Local Bray–Curtis
        bc_local = compute_local_bray_curtis(L_s, land_mask)
        # Whittaker summary
        wh = compute_whittaker_beta(L_s, land_mask, grid.lat_mesh)
        # Plots + files
        plot_diversity_maps(
            grid, land_mask, t_days, outdir,
            alpha_map=alpha_map,
            bc_local=bc_local,
            beta_whittaker=wh["beta_whittaker"],
            alpha_mean=wh["alpha_mean"],
            gamma_eff=wh["gamma_eff"],
        )
        save_community_npz(L_s, land_mask, t_days, outdir)
    except Exception:
        # Non-fatal: skip diagnostics
        pass

from __future__ import annotations

import os
import json
import numpy as np
from dataclasses import dataclass

from .spectral import make_bands, band_weights_from_mode, default_leaf_reflectance
from .population import PopulationManager
from .genes import Genes, Peak, reflectance_from_genes


@dataclass
class AdapterConfig:
    substep_every_nphys: int = 1
    lai_albedo_weight: float = 1.0
    feedback_mode: str = "instant"       # instant|daily
    couple_freq: str = "subdaily"        # subdaily|daily


class EcologyAdapter:
    """
    M1 adapter: compute a land-only ecological surface albedo (scalar 2D) each sub-daily call,
    based on a fixed leaf reflectance template and spectral band weights. This is a minimal
    implementation to wire hourly feedback without full population dynamics.

    Returns:
        alpha_surface_ecology_map (2D, same shape as I_total), with ocean cells set to NaN.
    """
    def __init__(self, grid, land_mask: np.ndarray):
        self.grid = grid
        self.land_mask = (land_mask == 1)
        self.cfg = AdapterConfig(
            substep_every_nphys=int(os.getenv("QD_ECO_SUBSTEP_EVERY_NPHYS", "1")),
            lai_albedo_weight=float(os.getenv("QD_ECO_LAI_ALBEDO_WEIGHT", "1.0")),
            feedback_mode=os.getenv("QD_ECO_FEEDBACK_MODE", "instant").strip().lower(),
            couple_freq=os.getenv("QD_ECO_ALBEDO_COUPLE_FREQ", "subdaily").strip().lower(),
        )
        # Evolution/mutation controls
        try:
            self.mut_rate = float(os.getenv("QD_ECO_MUT_RATE", "0.0"))
        except Exception:
            self.mut_rate = 0.0
        try:
            self.mut_eps = float(os.getenv("QD_ECO_MUT_EPS", "0.02"))
        except Exception:
            self.mut_eps = 0.02
        try:
            self.species_max = int(os.getenv("QD_ECO_SPECIES_MAX", "8"))
        except Exception:
            self.species_max = 8

        # Spectral setup
        self.bands = make_bands()
        self.w_b = band_weights_from_mode(self.bands)  # normalized weights sum=1
        self.R_leaf = default_leaf_reflectance(self.bands)  # [NB] in [0,1]
        # Precompute leaf scalar albedo under current weighting
        self.alpha_leaf_scalar = float(np.sum(self.R_leaf * self.w_b))
        # Counters
        self._step_count = 0

        # Band cache for M3b (optional)
        self._last_A_bands = None
        self._last_w_b = None

        # Diagnostics toggle
        self._diag = int(os.getenv("QD_ECO_DIAG", "1")) == 1
        if self._diag:
            print(f"[Ecology] M1 adapter init: NB={self.bands.nbands}, alpha_leaf≈{self.alpha_leaf_scalar:.3f}, "
                  f"substep_every={self.cfg.substep_every_nphys}, W_LAI={self.cfg.lai_albedo_weight:.2f}, "
                  f"feedback={self.cfg.feedback_mode}, couple_freq={self.cfg.couple_freq}")

        # M2: optional prognostic LAI manager
        use_lai = int(os.getenv("QD_ECO_USE_LAI", "1")) == 1
        self.pop = PopulationManager(self.land_mask.astype(int), diag=self._diag) if use_lai else None
        if self._diag and self.pop is not None:
            s = self.pop.summary()
            print(f"[Ecology] LAI init: min/mean/max={s['LAI_min']:.2f}/{s['LAI_mean']:.2f}/{s['LAI_max']:.2f} (use_lai={use_lai})")
        # Genes list cache
        self.genes_list: list[Genes] = []

        # M4 (species bands): build per-species leaf reflectance spectra R_leaf[b] from Genes
        # Species count derived from population species_weights length
        try:
            Ns = int(getattr(self.pop, "species_weights", np.asarray([1.0])).shape[0])
        except Exception:
            Ns = 1
        if Ns <= 0:
            Ns = 1
        R_species = []
        for i in range(Ns):
            # Allow per-species overrides via QD_ECO_SPECIES_{i}_*; fallback to default gene env (QD_ECO_GENE_*)
            try:
                g = Genes.from_env(prefix=f"QD_ECO_SPECIES_{i}_")
                R_i = reflectance_from_genes(self.bands, g)  # [NB]
                if not np.all(np.isfinite(R_i)):
                    raise ValueError("non-finite R_leaf")
                self.genes_list.append(g)
            except Exception:
                # Fallback: default template gene + template reflectance
                self.genes_list.append(Genes.from_env(prefix="QD_ECO_GENE_"))
                R_i = self.R_leaf.copy()
            R_species.append(np.clip(R_i, 0.0, 1.0))
        try:
            R_species_nb = np.stack(R_species, axis=0)  # [Ns, NB]
            if self.pop is not None:
                self.pop.set_species_reflectance_bands(R_species_nb)
            if self._diag:
                print(f"[Ecology] Species bands set: Ns={R_species_nb.shape[0]}, NB={R_species_nb.shape[1]}")
        except Exception as _es:
            if self._diag:
                print(f"[Ecology] Species bands not set (fallback to single template): {_es}")

    def step_subdaily(self,
                      I_total: np.ndarray,
                      cloud_eff: np.ndarray | float,
                      dt_seconds: float) -> np.ndarray | None:
        """
        Compute a land-only surface albedo map for ecology and return it
        for immediate coupling (if enabled). Returns None when not at substep boundary.
        Also accumulate daily energy into PopulationManager when启用 LAI 预报。
        """
        self._step_count += 1
        if self.pop is not None:
            try:
                self.pop.step_subdaily(I_total, dt_seconds)
            except Exception:
                pass

        if (self._step_count % max(1, int(self.cfg.substep_every_nphys))) != 0:
            return None

        # Base leaf scalar under current band weights
        # If LAI manager exists, scale leaf reflectance → canopy reflectance via f(LAI)
        alpha_map = np.full_like(I_total, np.nan, dtype=float)
        if self.pop is None:
            alpha_scalar = float(np.clip(self.alpha_leaf_scalar, 0.0, 1.0))
            alpha_map[self.land_mask] = alpha_scalar
            if self._diag and (self._step_count % 200 == 0):
                print(f"[Ecology] subdaily(M1): alpha_land={alpha_scalar:.3f}")
        else:
            try:
                f_canopy = self.pop.canopy_reflectance_factor()  # [lat,lon] in [0,1]
                # 合成陆面生态反照率：alpha = alpha_leaf_scalar * f(LAI) * W_LAI
                soil_ref = float(os.getenv("QD_ECO_SOIL_REFLECT", "0.20"))
                leaf_s = self.alpha_leaf_scalar  # Σ_b R_leaf[b]·w_b (w_b 归一)
                alpha_map[self.land_mask] = np.clip(
                    leaf_s * f_canopy[self.land_mask] + (1.0 - f_canopy[self.land_mask]) * soil_ref, 0.0, 1.0
                )
                if self._diag and (self._step_count % 200 == 0):
                    s = self.pop.summary()
                    am = alpha_map[self.land_mask]
                    print(f"[Ecology] subdaily(M2): LAI(min/mean/max)={s['LAI_min']:.2f}/{s['LAI_mean']:.2f}/{s['LAI_max']:.2f} | "
                          f"alpha_land(min/mean/max)={np.nanmin(am):.3f}/{np.nanmean(am):.3f}/{np.nanmax(am):.3f}")
            except Exception:
                # fallback to M1 scalar
                alpha_scalar = float(np.clip(self.alpha_leaf_scalar, 0.0, 1.0))
                alpha_map[self.land_mask] = alpha_scalar

        return alpha_map

    def export_genes(self, out_dir: str, day_value: float) -> None:
        """
        Export current species gene table (with weights) to JSON for audit/visualization.
        File: {out_dir}/genes_day_{day:.1f}.json
        """
        try:
            os.makedirs(out_dir, exist_ok=True)
            path = os.path.join(out_dir, f"genes_day_{day_value:05.1f}.json")
            table = []
            weights = None
            try:
                weights = [float(x) for x in np.asarray(getattr(self.pop, "species_weights", []), dtype=float).tolist()]
            except Exception:
                weights = None
            for i, g in enumerate(self.genes_list):
                entry = {
                    "index": i,
                    "identity": getattr(g, "identity", f"sp{i}"),
                    "alloc_root": float(getattr(g, "alloc_root", 0.0)),
                    "alloc_stem": float(getattr(g, "alloc_stem", 0.0)),
                    "alloc_leaf": float(getattr(g, "alloc_leaf", 0.0)),
                    "leaf_area_per_energy": float(getattr(g, "leaf_area_per_energy", 0.0)),
                    "drought_tolerance": float(getattr(g, "drought_tolerance", 0.0)),
                    "gdd_germinate": float(getattr(g, "gdd_germinate", 0.0)),
                    "lifespan_days": int(getattr(g, "lifespan_days", 0)),
                    "peaks": [
                        {"center_nm": float(pk.center_nm), "width_nm": float(pk.width_nm), "height": float(pk.height)}
                        for pk in getattr(g, "absorption_peaks", []) or []
                    ],
                }
                if weights is not None and i < len(weights):
                    entry["weight"] = weights[i]
                table.append(entry)
            with open(path, "w", encoding="utf-8") as f:
                json.dump({"day": day_value, "nbands": int(self.bands.nbands), "genes": table}, f, ensure_ascii=False, indent=2)
            if self._diag:
                print(f"[Ecology] Genes exported: {path} (Ns={len(table)})")
        except Exception as e:
            if self._diag:
                print(f"[Ecology] Genes export failed: {e}")

    # M2: daily LAI update from soil water proxy ∈[0,1]
    def step_daily(self, soil_water_index: np.ndarray | float | None) -> None:
        if self.pop is None:
            return
        try:
            self.pop.step_daily(soil_water_index)
            if self._diag:
                s = self.pop.summary()
                print(f"[Ecology] daily: LAI(min/mean/max)={s['LAI_min']:.2f}/{s['LAI_mean']:.2f}/{s['LAI_max']:.2f}")
            # Evolution: stochastic mutation/new species (minimal M4)
            if self.mut_rate > 0.0 and np.random.rand() < self.mut_rate:
                S_now = int(getattr(self.pop, "Ns", len(self.genes_list) or 1))
                if S_now < self.species_max:
                    # choose parent by current species_weights if available
                    try:
                        w = np.asarray(self.pop.species_weights, dtype=float)
                        w = w / (np.sum(w) + 1e-12)
                        parent = int(np.random.choice(np.arange(S_now), p=w))
                    except Exception:
                        parent = int(np.random.randint(0, max(1, S_now)))
                    idx_new = self.pop.add_species_from_parent(parent, frac=self.mut_eps)
                    # mutate genes from parent
                    try:
                        g_parent = self.genes_list[parent] if parent < len(self.genes_list) else Genes.from_env()
                        g_new = self._mutate_genes(g_parent)
                    except Exception:
                        g_new = Genes.from_env()
                    # append and rebuild species reflectance table
                    if idx_new >= len(self.genes_list):
                        self.genes_list.append(g_new)
                    else:
                        # safety: extend list
                        self.genes_list = (self.genes_list + [g_new])[:idx_new+1]
                    R_species_nb = np.stack([reflectance_from_genes(self.bands, g) for g in self.genes_list], axis=0)
                    self.pop.set_species_reflectance_bands(R_species_nb)
                    if self._diag:
                        print(f"[Ecology] mutation: parent={parent} → new species idx={idx_new}; Ns={len(self.genes_list)}")
                elif self._diag:
                    print(f"[Ecology] mutation skipped: species_max={self.species_max} reached.")
        except Exception as e:
            if self._diag:
                print(f"[Ecology] daily step failed: {e}")

    def _mutate_genes(self, g: Genes) -> Genes:
        """
        Create a slightly perturbed copy of Genes to represent a mutation.
        Jitters spectral peaks and a few physiological/allocation parameters with bounds.
        """
        g2 = Genes(
            identity=(g.identity + "_mut"),
            alloc_root=g.alloc_root,
            alloc_stem=g.alloc_stem,
            alloc_leaf=g.alloc_leaf,
            leaf_area_per_energy=g.leaf_area_per_energy,
            absorption_peaks=[Peak(pk.center_nm, pk.width_nm, pk.height) for pk in g.absorption_peaks],
            drought_tolerance=g.drought_tolerance,
            gdd_germinate=g.gdd_germinate,
            lifespan_days=g.lifespan_days,
        )
        # Allocation jitter then renormalize
        jit = 0.05
        g2.alloc_root = float(np.clip(g2.alloc_root + np.random.uniform(-jit, jit), 0.05, 0.90))
        g2.alloc_stem = float(np.clip(g2.alloc_stem + np.random.uniform(-jit, jit), 0.05, 0.90))
        g2.alloc_leaf = float(np.clip(g2.alloc_leaf + np.random.uniform(-jit, jit), 0.05, 0.90))
        s = g2.alloc_root + g2.alloc_stem + g2.alloc_leaf
        g2.alloc_root /= s; g2.alloc_stem /= s; g2.alloc_leaf /= s
        # Spectral peaks jitter
        for pk in g2.absorption_peaks:
            pk.center_nm = float(np.clip(pk.center_nm + np.random.normal(0.0, 8.0), 380.0, 780.0))
            pk.width_nm  = float(np.clip(pk.width_nm  + np.random.normal(0.0, 5.0), 10.0, 120.0))
            pk.height    = float(np.clip(pk.height    + np.random.normal(0.0, 0.05), 0.05, 0.98))
        # Physiology jitter
        g2.drought_tolerance = float(np.clip(g2.drought_tolerance + np.random.normal(0.0, 0.03), 0.05, 0.95))
        g2.gdd_germinate     = float(np.clip(g2.gdd_germinate + np.random.normal(0.0, 5.0), 10.0, 500.0))
        g2.lifespan_days     = int(np.clip(g2.lifespan_days + np.random.normal(0.0, 30.0), 30, 365*5))
        g2.leaf_area_per_energy = float(np.clip(g2.leaf_area_per_energy * (1.0 + np.random.normal(0.0, 0.1)),
                                                1e-5, 5e-2))
        # Environment-biased spectral drift: nudge centers toward current weighted band center
        try:
            lam = np.asarray(self.bands.lambda_centers, dtype=float)
            wb = np.asarray(self.w_b, dtype=float)
            lam_w = float(np.sum(lam * wb) / (np.sum(wb) + 1e-12))
            alpha = float(os.getenv("QD_ECO_MUT_LAMBDA_DRIFT", "0.1"))
            for pk in g2.absorption_peaks:
                pk.center_nm = float(np.clip(pk.center_nm + alpha * (lam_w - pk.center_nm), 380.0, 780.0))
        except Exception:
            pass
        return g2

    # M3b: provide banded surface albedo A_b^surface (NB×lat×lon) for future shortwave band coupling
    # Returns (A_bands, w_b) or (None, None) if not available.
    def get_surface_albedo_bands(self) -> tuple[np.ndarray | None, np.ndarray | None]:
        try:
            nb = self.bands.nbands
            soil_ref = float(os.getenv("QD_ECO_SOIL_REFLECT", "0.20"))
            # Prefer population-provided cohort/species mixing if available
            if self.pop is not None:
                A = self.pop.get_surface_albedo_bands(nb, soil_ref=soil_ref)  # [NB, lat, lon]
                self._last_A_bands = A
                self._last_w_b = self.w_b.copy()
                return A, self._last_w_b

            # Fallback: single-template (no species/cohort info), land-only
            h, w = self.grid.lat_mesh.shape
            A = np.full((nb, h, w), np.nan, dtype=float)
            f_canopy = np.full((h, w), np.nan, dtype=float)
            f_canopy[self.land_mask] = 1.0
            for b in range(nb):
                leaf_rb = float(self.R_leaf[b])  # [0,1]
                A_b2d = leaf_rb * f_canopy + (1.0 - f_canopy) * soil_ref
                Ab = np.full((h, w), np.nan, dtype=float)
                Ab[self.land_mask] = np.clip(A_b2d[self.land_mask], 0.0, 1.0)
                A[b, :, :] = Ab
            self._last_A_bands = A
            self._last_w_b = self.w_b.copy()
            return A, self._last_w_b
        except Exception:
            return None, None

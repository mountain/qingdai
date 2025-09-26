from __future__ import annotations

import os
import json
import numpy as np
from dataclasses import dataclass
import time
import glob

from .spectral import make_bands, band_weights_from_mode, default_leaf_reflectance
from .population import PopulationManager
from .genes import Genes, Peak, reflectance_from_genes, absorbance_from_genes


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
        # Genotype absorbance cache (M1.5)
        self._absorb_cache: dict[str, np.ndarray] = {}
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
                A_i = self._absorb_from_genes_cached(g)  # [NB] in [0,1]
                R_i = np.clip(1.0 - A_i, 0.0, 1.0)
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

        # Map species identities by spread modes: 'seed' → 'tree', 'diffusion' → 'grass'
        # Respect explicit per-species identity overrides via QD_ECO_SPECIES_{i}_IDENTITY if provided.
        try:
            modes = getattr(self.pop, "species_modes", []) if self.pop is not None else []
            for i, g in enumerate(self.genes_list):
                # Skip if user explicitly provided IDENTITY for this species
                if os.getenv(f"QD_ECO_SPECIES_{i}_IDENTITY"):
                    continue
                mode_i = (modes[i] if (i < len(modes) and modes[i] in ("seed", "diffusion"))
                          else ("seed" if i == 1 else "diffusion"))
                g.identity = ("tree" if mode_i == "seed" else "grass")
            if self._diag:
                n_tree = sum(1 for i, g in enumerate(self.genes_list)
                             if (os.getenv(f"QD_ECO_SPECIES_{i}_IDENTITY") is None) and g.identity == "tree")
                n_grass = sum(1 for i, g in enumerate(self.genes_list)
                              if (os.getenv(f"QD_ECO_SPECIES_{i}_IDENTITY") is None) and g.identity == "grass")
                print(f"[Ecology] Species identities mapped by modes: grass={n_grass}, tree={n_tree}")
        except Exception as _ei:
            if self._diag:
                print(f"[Ecology] Species identity mapping skipped: {_ei}")

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
        Export current species gene table (with weights) to JSON for audit/visualization and sharing.
        File: {out_dir}/genes_day_{day:.1f}.json

        Schema v3 fields:
          - schema_version: int
          - source: string
          - day: float (planetary days)
          - bands: { nbands, lambda_centers_nm[], delta_lambda_nm[], lambda_edges_nm[] }
          - band_weights: w_b[] (normalized)
          - genes: [
              {
                index, identity, provenance,
                alloc_root, alloc_stem, alloc_leaf,
                leaf_area_per_energy,
                drought_tolerance, gdd_germinate, lifespan_days,
                peaks: [ { center_nm, sigma_nm, variance_nm2, height } ],
                weight (optional)
              }, ...
            ]
        """
        try:
            os.makedirs(out_dir, exist_ok=True)
            path = os.path.join(out_dir, f"genes_day_{day_value:05.1f}.json")
            table = []
            # Species weights (if available)
            try:
                weights = [float(x) for x in np.asarray(getattr(self.pop, "species_weights", []), dtype=float).tolist()]
            except Exception:
                weights = None

            for i, g in enumerate(self.genes_list):
                peaks = getattr(g, "absorption_peaks", []) or []
                peaks_out = []
                for pk in peaks:
                    try:
                        sigma = float(pk.width_nm)
                        peaks_out.append({
                            "center_nm": float(pk.center_nm),
                            "sigma_nm": sigma,
                            "variance_nm2": float(sigma * sigma),
                            "height": float(pk.height),
                        })
                    except Exception:
                        # Minimal fallback
                        peaks_out.append({
                            "center_nm": float(getattr(pk, "center_nm", 0.0)),
                            "sigma_nm": float(getattr(pk, "width_nm", 0.0)),
                            "variance_nm2": float(getattr(pk, "width_nm", 0.0)) ** 2,
                            "height": float(getattr(pk, "height", 0.0)),
                        })

                entry = {
                    "index": i,
                    "identity": getattr(g, "identity", f"sp{i}"),
                    "provenance": getattr(g, "provenance", None),
                    "alloc_root": float(getattr(g, "alloc_root", 0.0)),
                    "alloc_stem": float(getattr(g, "alloc_stem", 0.0)),
                    "alloc_leaf": float(getattr(g, "alloc_leaf", 0.0)),
                    "leaf_area_per_energy": float(getattr(g, "leaf_area_per_energy", 0.0)),
                    "drought_tolerance": float(getattr(g, "drought_tolerance", 0.0)),
                    "gdd_germinate": float(getattr(g, "gdd_germinate", 0.0)),
                    "lifespan_days": int(getattr(g, "lifespan_days", 0)),
                    "peaks_model": "gaussian",
                    "peaks": peaks_out,
                    # Embed band definitions per-gene to make the entry self-contained
                    "lambda_centers_nm": [float(x) for x in np.asarray(self.bands.lambda_centers, dtype=float).tolist()],
                    "delta_lambda_nm": [float(x) for x in np.asarray(self.bands.delta_lambda, dtype=float).tolist()],
                    "lambda_edges_nm": [float(x) for x in np.asarray(self.bands.lambda_edges, dtype=float).tolist()],
                }
                if weights is not None and i < len(weights):
                    entry["weight"] = weights[i]
                table.append(entry)

            doc = {
                "schema_version": 3,
                "source": "PyGCM.EcologyAdapter.export_genes",
                "day": float(day_value),
                "bands": {
                    "nbands": int(self.bands.nbands),
                    "band_weights": [float(x) for x in np.asarray(self.w_b, dtype=float).tolist()]
                },
                "genes": table,
            }

            with open(path, "w", encoding="utf-8") as f:
                json.dump(doc, f, ensure_ascii=False, indent=2)

            if self._diag:
                print(f"[Ecology] Genes exported: {path} (Ns={len(table)}, schema=3)")
        except Exception as e:
            if self._diag:
                print(f"[Ecology] Genes export failed: {e}")

    # --- Genes autosave helpers (JSON) ---
    def save_genes_json(self, path: str, day_value: float | None = None) -> bool:
        """
        Save current core gene information (per-species Genes + band weights) to a stable JSON:
        - Intended for autosave/resume under data/genes_autosave.json
        - Compact but self-contained: includes bands nb/weights; each gene embeds peak list
        """
        try:
            os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        except Exception:
            pass
        try:
            # Species weights (optional snapshot)
            try:
                weights = [float(x) for x in np.asarray(getattr(self.pop, "species_weights", []), dtype=float).tolist()]
            except Exception:
                weights = None
            genes_tbl = []
            for i, g in enumerate(self.genes_list):
                peaks = getattr(g, "absorption_peaks", []) or []
                peaks_out = []
                for pk in peaks:
                    try:
                        sigma = float(pk.width_nm)
                        peaks_out.append({
                            "center_nm": float(pk.center_nm),
                            "sigma_nm": sigma,
                            "variance_nm2": float(sigma * sigma),
                            "height": float(pk.height),
                        })
                    except Exception:
                        peaks_out.append({
                            "center_nm": float(getattr(pk, "center_nm", 0.0)),
                            "sigma_nm": float(getattr(pk, "width_nm", 0.0)),
                            "variance_nm2": float(getattr(pk, "width_nm", 0.0)) ** 2,
                            "height": float(getattr(pk, "height", 0.0)),
                        })
                entry = {
                    "index": i,
                    "identity": getattr(g, "identity", f"sp{i}"),
                    "provenance": getattr(g, "provenance", None),
                    "alloc_root": float(getattr(g, "alloc_root", 0.0)),
                    "alloc_stem": float(getattr(g, "alloc_stem", 0.0)),
                    "alloc_leaf": float(getattr(g, "alloc_leaf", 0.0)),
                    "leaf_area_per_energy": float(getattr(g, "leaf_area_per_energy", 0.0)),
                    "drought_tolerance": float(getattr(g, "drought_tolerance", 0.0)),
                    "gdd_germinate": float(getattr(g, "gdd_germinate", 0.0)),
                    "lifespan_days": int(getattr(g, "lifespan_days", 0)),
                    "peaks_model": "gaussian",
                    "peaks": peaks_out,
                }
                genes_tbl.append(entry)
            doc = {
                "schema_version": 3,
                "source": "PyGCM.EcologyAdapter.save_genes_json",
                "day": float(day_value) if day_value is not None else None,
                "bands": {
                    "nbands": int(self.bands.nbands),
                    "band_weights": [float(x) for x in np.asarray(self.w_b, dtype=float).tolist()]
                },
                "genes": genes_tbl,
            }
            if weights is not None and len(weights) > 0:
                doc["species_weights"] = weights
            with open(path, "w", encoding="utf-8") as f:
                json.dump(doc, f, ensure_ascii=False, indent=2)
            if self._diag:
                print(f"[Ecology] Genes autosave written: {path} (Ns={len(genes_tbl)})")
            return True
        except Exception as e:
            if self._diag:
                print(f"[Ecology] Genes autosave save failed: {e}")
            return False

    def load_genes_json(self, path: str, *, on_mismatch: str = "keep") -> bool:
        """
        Load genes autosave JSON (data/genes_autosave.json) and rebuild per-species reflectance.
        - If bands mismatch, we still rebuild reflectance using current bands from stored peaks.
        - If population present, update species reflectance table to new Ns; weights apply later via NPZ.
        """
        try:
            with open(path, "r", encoding="utf-8") as f:
                doc = json.load(f)
        except Exception as e:
            if self._diag:
                print(f"[Ecology] Genes autosave load failed: {e}")
            return False
        try:
            genes_in = []
            for rec in doc.get("genes", []):
                # Reconstruct Peaks and Genes
                peaks = []
                for pk in rec.get("peaks", []) or []:
                    try:
                        # prefer sigma_nm if provided; variance_nm2 kept for back-compat
                        sigma = float(pk.get("sigma_nm", 0.0))
                        if sigma <= 0 and "variance_nm2" in pk:
                            var = float(pk.get("variance_nm2", 0.0))
                            sigma = float(np.sqrt(max(0.0, var)))
                        peaks.append(Peak(float(pk.get("center_nm", 0.0)), sigma, float(pk.get("height", 0.0))))
                    except Exception:
                        continue
                g = Genes(
                    identity=str(rec.get("identity", "sp")),
                    alloc_root=float(rec.get("alloc_root", 0.3)),
                    alloc_stem=float(rec.get("alloc_stem", 0.2)),
                    alloc_leaf=float(rec.get("alloc_leaf", 0.5)),
                    leaf_area_per_energy=float(rec.get("leaf_area_per_energy", 2.0e-3)),
                    absorption_peaks=peaks,
                    drought_tolerance=float(rec.get("drought_tolerance", 0.3)),
                    gdd_germinate=float(rec.get("gdd_germinate", 80.0)),
                    lifespan_days=int(rec.get("lifespan_days", 365)),
                    provenance="autosave:genes_json",
                )
                # normalize allocations
                s = g.alloc_root + g.alloc_stem + g.alloc_leaf
                if s > 0:
                    g.alloc_root /= s; g.alloc_stem /= s; g.alloc_leaf /= s
                genes_in.append(g)

            if len(genes_in) == 0:
                if self._diag:
                    print("[Ecology] Genes autosave contains no genes; keeping current.")
                return False

            # Apply: replace adapter genes list
            self.genes_list = genes_in

            # Rebuild per-species reflectance with current bands
            try:
                R_species_nb = np.stack([reflectance_from_genes(self.bands, g) for g in self.genes_list], axis=0)
                if getattr(self, "pop", None) is not None:
                    self.pop.set_species_reflectance_bands(R_species_nb)
            except Exception as e:
                if self._diag:
                    print(f"[Ecology] Genes reflectance rebuild failed: {e}")

            if self._diag:
                print(f"[Ecology] Genes autosave loaded: Ns={len(self.genes_list)} (band nb={self.bands.nbands})")
            return True
        except Exception as e:
            if self._diag:
                print(f"[Ecology] Genes autosave parse failed: {e}")
            return False

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

    # --- M1.5: Genotype absorbance cache helpers ---
    def _absorb_from_genes_cached(self, g: Genes) -> np.ndarray:
        """
        Return band absorbance A_b[NB] for genes 'g' using cache keyed by
        (identity, peaks, band centers). Falls back to compute if missing.
        """
        try:
            peaks = getattr(g, "absorption_peaks", []) or []
            pk_key = tuple((float(getattr(p, "center_nm", 0.0)),
                            float(getattr(p, "width_nm", 0.0)),
                            float(getattr(p, "height", 0.0))) for p in peaks)
            lam_key = tuple(float(x) for x in np.asarray(self.bands.lambda_centers, dtype=float).ravel().tolist())
            key = f"{getattr(g, 'identity', 'gene')}|{pk_key}|{lam_key}"
        except Exception:
            key = f"{id(g)}|{self.bands.nbands}"
        A = self._absorb_cache.get(key)
        if A is None:
            try:
                A = absorbance_from_genes(self.bands, g)
                A = np.clip(np.asarray(A, dtype=float).ravel(), 0.0, 1.0)
            except Exception:
                # Fallback flat absorbance ~ 0.5
                A = np.full((self.bands.nbands,), 0.5, dtype=float)
            self._absorb_cache[key] = A
        return A

    # --- M1.5: Autosave / Load of extended ecology state ---
    def save_autosave(self, path_npz: str, day_value: float | None = None) -> bool:
        """
        Save extended ecology state to a single file with atomic replace and optional rolling backups.
        Supports:
          - NPZ (legacy/debug)
          - NetCDF (preferred): ecology.nc (standardized)
        Env (optional):
          - QD_ECO_AUTOSAVE_KEEP (default 4): number of timestamped backups to keep (applies to same extension)
        """
        try:
            out_dir = os.path.dirname(path_npz) or "."
            base = os.path.basename(path_npz)
            name, ext = os.path.splitext(base)
            os.makedirs(out_dir, exist_ok=True)
            ts = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
            backup_path = os.path.join(out_dir, f"{name}_{ts}{ext}")
            tmp_path = os.path.join(out_dir, f".{name}.tmp{ext}")

            # Collect data
            LAI_src = None
            species_w = None
            R_species_nb = None
            if getattr(self, "pop", None) is not None:
                LAI_src = getattr(self.pop, "LAI", None)
                if LAI_src is None and hasattr(self.pop, "total_LAI"):
                    LAI_src = self.pop.total_LAI()
                species_w = np.asarray(getattr(self.pop, "species_weights", np.asarray([1.0])), dtype=float)
                if getattr(self.pop, "_species_R_leaf", None) is not None:
                    R_species_nb = np.asarray(self.pop._species_R_leaf, dtype=float)

            # Extension switch
            if ext.lower() == ".nc":
                # Write NetCDF
                from netCDF4 import Dataset
                with Dataset(tmp_path, "w") as ds:
                    # Dimensions
                    nlat, nlon = self.grid.n_lat, self.grid.n_lon
                    ds.createDimension("lat", nlat)
                    ds.createDimension("lon", nlon)
                    if species_w is not None:
                        ds.createDimension("species", int(species_w.size))
                    ds.createDimension("band", int(self.bands.nbands))

                    # Coordinates
                    vlat = ds.createVariable("lat", "f4", ("lat",))
                    vlon = ds.createVariable("lon", "f4", ("lon",))
                    vlat[:] = self.grid.lat
                    vlon[:] = self.grid.lon

                    # Core fields
                    if LAI_src is not None:
                        vLAI = ds.createVariable("LAI", "f4", ("lat", "lon"))
                        vLAI[:] = np.asarray(LAI_src, dtype=np.float32)
                    if species_w is not None:
                        vsw = ds.createVariable("species_weights", "f4", ("species",))
                        vsw[:] = species_w.astype(np.float32)
                    # Bands
                    vcent = ds.createVariable("bands_lambda_centers", "f4", ("band",))
                    vdl = ds.createVariable("bands_delta_lambda", "f4", ("band",))
                    vw = ds.createVariable("w_b", "f4", ("band",))
                    vcent[:] = np.asarray(self.bands.lambda_centers, dtype=np.float32)
                    vdl[:] = np.asarray(self.bands.delta_lambda, dtype=np.float32)
                    vw[:] = np.asarray(self.w_b, dtype=np.float32)

                    # Species reflectance
                    if R_species_nb is not None and species_w is not None:
                        vR = ds.createVariable("R_species_nb", "f4", ("species", "band"))
                        vR[:] = R_species_nb.astype(np.float32)

                    # Optional day
                    if day_value is not None:
                        vday = ds.createVariable("day_value", "f4")
                        vday[...] = float(day_value)

                    # Metadata
                    ds.setncattr("title", "Qingdai Ecology State")
                    ds.setncattr("schema_version", 1)
                    ds.setncattr("source", "EcologyAdapter.save_autosave")
                os.replace(tmp_path, path_npz)
            else:
                # NPZ fallback
                data = {}
                data["schema_version"] = np.int32(1)
                if LAI_src is not None:
                    data["LAI"] = np.asarray(LAI_src, dtype=float)
                if species_w is not None:
                    data["species_weights"] = species_w
                data["bands_lambda_centers"] = np.asarray(self.bands.lambda_centers, dtype=float)
                data["bands_delta_lambda"] = np.asarray(self.bands.delta_lambda, dtype=float)
                data["w_b"] = np.asarray(self.w_b, dtype=float)
                if R_species_nb is not None:
                    data["R_species_nb"] = R_species_nb
                if day_value is not None:
                    data["day_value"] = float(day_value)
                with open(tmp_path, "wb") as f:
                    np.savez(f, **data)
                    f.flush()
                    os.fsync(f.fileno())
                os.replace(tmp_path, path_npz)

            # Timestamped backup (best-effort)
            try:
                import shutil
                shutil.copy2(path_npz, backup_path)
            except Exception:
                backup_path = None

            # Rolling retention
            try:
                keep = int(os.getenv("QD_ECO_AUTOSAVE_KEEP", "4"))
                pattern = os.path.join(out_dir, f"{name}_*{ext}")
                files = sorted(glob.glob(pattern), key=lambda p: os.path.getmtime(p), reverse=True)
                for old in files[keep:]:
                    try:
                        os.remove(old)
                    except Exception:
                        pass
            except Exception:
                pass

            # Also write genes.json next to ecology state (best-effort, standardized name)
            try:
                genes_json = os.path.join(out_dir, "genes.json")
                _ = self.save_genes_json(genes_json, day_value=day_value if day_value is not None else None)
            except Exception:
                pass

            if self._diag:
                msg = f"[Ecology] Autosave written: {path_npz}"
                if backup_path:
                    msg += f"; backup={backup_path}"
                print(msg)
            return True
        except Exception as e:
            if self._diag:
                print(f"[Ecology] Autosave+ load failed: {e}")
            return False

    def load_autosave(self, path_npz: str, *, on_mismatch: str = "fallback") -> bool:
        """
        Load extended ecology state. If bands length mismatches, we only restore LAI/species weights.
        on_mismatch: 'fallback' | 'ignore' – behavior when advanced fields mismatch.
        """
        if getattr(self, "pop", None) is None:
            return False
        try:
            import os as _os
            ext = _os.path.splitext(path_npz)[1].lower()
            centers = None
            R_species_nb = None
            if ext == ".nc":
                from netCDF4 import Dataset as _DS
                with _DS(path_npz, "r") as ds:
                    LAI = np.asarray(ds.variables["LAI"]) if "LAI" in ds.variables else None
                    w = np.asarray(ds.variables["species_weights"]) if "species_weights" in ds.variables else None
                    centers = np.asarray(ds.variables["bands_lambda_centers"]) if "bands_lambda_centers" in ds.variables else None
                    if "R_species_nb" in ds.variables:
                        R_species_nb = np.asarray(ds.variables["R_species_nb"])
            else:
                data = np.load(path_npz)
                LAI = np.asarray(data.get("LAI"))
                w = np.asarray(data.get("species_weights"))
                centers = data.get("bands_lambda_centers")
                R_species_nb = data.get("R_species_nb")
            if LAI is None or LAI.ndim != 2 or w is None or w.ndim != 1:
                if self._diag:
                    print("[Ecology] Autosave+ malformed: require LAI (2D) and species_weights (1D).")
                return False
            # Restore LAI/species weights (like legacy helper)
            pop = self.pop
            pop.LAI = np.clip(LAI, 0.0, pop.params.lai_max)
            S = int(w.size)
            K = int(getattr(pop, "K", 1))
            H, W = pop.shape
            w = np.clip(w, 0.0, None)
            ssum = float(np.sum(w))
            pop.species_weights = (w / ssum) if ssum > 0 else np.full((S,), 1.0 / float(max(S, 1)), dtype=float)
            pop.Ns = int(pop.species_weights.shape[0])
            pop.LAI_layers_SK = np.zeros((pop.Ns, max(1, K), H, W), dtype=float)
            for s_idx in range(pop.Ns):
                frac_s = float(pop.species_weights[s_idx])
                for k in range(max(1, K)):
                    pop.LAI_layers_SK[s_idx, k, :, :] = frac_s * (pop.LAI / float(max(1, K)))
            pop.LAI_layers = np.sum(pop.LAI_layers_SK, axis=0)

            # Try to restore bands & reflectance if shapes compatible
            ok_adv = (centers is not None and len(centers) == self.bands.nbands and
                      R_species_nb is not None and R_species_nb.shape[1] == self.bands.nbands)
            if ok_adv:
                try:
                    pop.set_species_reflectance_bands(np.asarray(R_species_nb, dtype=float))
                except Exception:
                    ok_adv = False
            elif self._diag and on_mismatch != "ignore":
                print("[Ecology] Autosave+: advanced fields (bands/R_species) mismatched – restored base LAI/weights only.")
            if self._diag:
                sdiag = pop.summary()
                print(f"[Ecology] Autosave+ loaded: LAI(min/mean/max)={sdiag['LAI_min']:.2f}/{sdiag['LAI_mean']:.2f}/{sdiag['LAI_max']:.2f}; "
                      f"S={pop.Ns}, K={K}, advanced={'OK' if ok_adv else 'NO'}")
            return True
        except Exception as e:
            if self._diag:
                print(f"[Ecology] Autosave+ load failed: {e}")
            return False

    # --- M1.5: Genotype absorbance cache helpers ---
    def _absorb_from_genes_cached(self, g: Genes) -> np.ndarray:
        """
        Return band absorbance A_b[NB] for genes 'g' using cache keyed by
        (identity, peaks, band centers). Falls back to compute if missing.
        """
        try:
            peaks = getattr(g, "absorption_peaks", []) or []
            pk_key = tuple((float(getattr(p, "center_nm", 0.0)),
                            float(getattr(p, "width_nm", 0.0)),
                            float(getattr(p, "height", 0.0))) for p in peaks)
            lam_key = tuple(float(x) for x in np.asarray(self.bands.lambda_centers, dtype=float).ravel().tolist())
            key = f"{getattr(g, 'identity', 'gene')}|{pk_key}|{lam_key}"
        except Exception:
            key = f"{id(g)}|{self.bands.nbands}"
        A = self._absorb_cache.get(key)
        if A is None:
            try:
                A = absorbance_from_genes(self.bands, g)
                A = np.clip(np.asarray(A, dtype=float).ravel(), 0.0, 1.0)
            except Exception:
                # Fallback flat absorbance ~ 0.5
                A = np.full((self.bands.nbands,), 0.5, dtype=float)
            self._absorb_cache[key] = A
        return A

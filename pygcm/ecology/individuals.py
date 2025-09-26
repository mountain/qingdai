from __future__ import annotations

import os
from dataclasses import dataclass
import numpy as np

try:
    # Dual-star ISR to bands utility and bands structure
    from .spectral import dual_star_insolation_to_bands
except Exception:
    dual_star_insolation_to_bands = None


@dataclass
class IndividualPoolConfig:
    sample_frac: float = 0.02            # fraction of land cells to sample
    per_cell: int = 150                  # individuals per sampled cell (kept consistent with runtime default)
    substeps_per_day: int = 10           # K: subdaily steps per planetary day
    nb_max: int = 16                     # maximum supported spectral bands (for sanity checks)
    diag: bool = True


class IndividualPool:
    """
    Vectorized pool of 'individuals' for a sampled subset of land cells.
    Goal:
      - Capture subdaily spectral adaptation (K substeps/day) and species-specific growth signals
      - Keep cost bounded by (sample_frac * n_cells * per_cell)
      - Aggregate per-cell species weights daily and feed back into PopulationManager (LAI_layers_SK split)

    Notes:
      - Individuals carry (species_id, energy_buffer, water_stress_days, per-band reflectance).
      - Subdaily: accumulate E_day_indiv via dot(Ab_i, I_bands(cell)).
      - Daily: compute per-cell species energy sums → weights; apply to PopulationManager species split in sampled cells.
    """

    def __init__(self,
                 grid,
                 land_mask: np.ndarray,
                 eco_adapter,  # EcologyAdapter (needs .bands, .pop with species_weights, .pop.set_species_reflectance_bands already called)
                 *,
                 sample_frac: float = 0.02,
                 per_cell: int = 100,
                 substeps_per_day: int = 10,
                 diag: bool = True) -> None:
        self.grid = grid
        self.land_mask = (land_mask == 1)
        self.h, self.w = self.grid.lat_mesh.shape
        self.cfg = IndividualPoolConfig(
            sample_frac=float(os.getenv("QD_ECO_INDIV_SAMPLE_FRAC", str(sample_frac))),
            per_cell=int(os.getenv("QD_ECO_INDIV_PER_CELL", str(per_cell))),
            substeps_per_day=max(1, int(os.getenv("QD_ECO_INDIV_SUBSTEPS_PER_DAY", str(substeps_per_day)))),
            nb_max=int(os.getenv("QD_ECO_SPECTRAL_BANDS", "16")) if os.getenv("QD_ECO_SPECTRAL_BANDS") else 16,
            diag=(int(os.getenv("QD_ECO_DIAG", "1")) == 1) and diag
        )
        self.bands = getattr(eco_adapter, "bands", None)
        if self.bands is None:
            raise RuntimeError("IndividualPool requires EcologyAdapter with bands.")
        self.nb = int(getattr(self.bands, "nbands", 0) or getattr(self.bands, "nb", 0) or 0)
        if self.nb <= 0 or self.nb > self.cfg.nb_max:
            # fallback to reasonable band count
            self.nb = min(max(self.nb, 8), self.cfg.nb_max)

        # Species weight distribution from population manager
        pop = getattr(eco_adapter, "pop", None)
        if pop is None:
            raise RuntimeError("IndividualPool requires EcologyAdapter.pop (PopulationManager).")
        sp_weights = np.asarray(getattr(pop, "species_weights", np.asarray([1.0])), dtype=float)
        if sp_weights.ndim != 1 or sp_weights.size <= 0:
            sp_weights = np.asarray([1.0], dtype=float)
        self.ns = int(sp_weights.size)
        ssum = float(np.sum(sp_weights))
        self.sp_weights = (sp_weights / ssum) if ssum > 0 else np.full((self.ns,), 1.0 / float(self.ns), dtype=float)

        # Build sampled cell list (lat, lon indices)
        land_idx = np.flatnonzero(self.land_mask.ravel())
        n_land = int(land_idx.size)
        n_cells_sampled = max(1, int(self.cfg.sample_frac * n_land))
        rng = np.random.default_rng(seed=42)
        if n_cells_sampled >= n_land:
            sampled = land_idx
        else:
            sampled = rng.choice(land_idx, size=n_cells_sampled, replace=False)
        jj = sampled // self.w
        ii = sampled % self.w
        self.sample_j = np.asarray(jj, dtype=np.int32)
        self.sample_i = np.asarray(ii, dtype=np.int32)
        self.n_cells = int(self.sample_j.size)

        # Per-cell individuals
        self.per_cell = int(self.cfg.per_cell)
        self.n_indiv = int(self.n_cells * self.per_cell)

        # Map each individual to its cell index [0..n_cells)
        self.indiv_cell_index = np.repeat(np.arange(self.n_cells, dtype=np.int32), self.per_cell)

        # Draw species_id for each individual based on global species weights
        self.indiv_species_id = rng.choice(np.arange(self.ns, dtype=np.int32), size=self.n_indiv, p=self.sp_weights)

        # Per-individual per-band reflectance Ab_i[b]
        # Prefer species leaf reflectance bands from PopulationManager if available
        species_R = getattr(pop, "_species_R_leaf", None)
        if species_R is None or species_R.shape[0] != self.ns:
            # fallback to equal-neutral reflectance 0.5
            species_R = np.full((self.ns, self.nb), 0.5, dtype=float)
        if species_R.shape[1] != self.nb:
            # band count mismatch: pad/crop
            if species_R.shape[1] > self.nb:
                species_R = species_R[:, :self.nb]
            else:
                species_R = np.pad(species_R, ((0, 0), (0, self.nb - species_R.shape[1])), mode="edge")
        # Small per-individual jitter to encourage micro-diversity
        jitter = 0.02
        Ab = species_R[self.indiv_species_id, :] + rng.normal(0.0, jitter, size=(self.n_indiv, self.nb))
        self.indiv_Ab = np.clip(Ab, 0.0, 1.0)

        # Drought tolerance per species if present on eco_adapter.genes_list; fallback to 0.5
        genes_list = getattr(eco_adapter, "genes_list", None)
        tol = np.full((self.ns,), 0.5, dtype=float)
        if genes_list and len(genes_list) == self.ns:
            for s in range(self.ns):
                try:
                    tol[s] = float(getattr(genes_list[s], "drought_tolerance", 0.5))
                except Exception:
                    tol[s] = 0.5
        self.species_drought_tol = np.clip(tol, 0.0, 1.0)
        self.indiv_tol = self.species_drought_tol[self.indiv_species_id]

        # State buffers
        self.indiv_E_day = np.zeros((self.n_indiv,), dtype=float)
        self.indiv_water_stress_days = np.zeros((self.n_indiv,), dtype=float)

        # Substep scheduler
        self._substep_period = None  # set in try_substep once we know day length
        self._substep_accum = 0.0

        if self.cfg.diag:
            print(f"[EcoIndiv] initialized: sample_frac={self.cfg.sample_frac:.3f}, "
                  f"cells={self.n_cells}, per_cell={self.per_cell}, N={self.n_indiv}, "
                  f"NB={self.nb}, K={self.cfg.substeps_per_day}")

    def try_substep(self,
                    isr_A: np.ndarray,
                    isr_B: np.ndarray,
                    eco_adapter,
                    soil_W_land: np.ndarray | float,
                    dt_seconds: float,
                    day_length_seconds: float) -> None:
        """
        Potentially perform one subdaily step if accumulated dt exceeds period.
        Accumulate per-individual energy from banded irradiance and update stress buffer.
        """
        if dual_star_insolation_to_bands is None:
            return
        if self._substep_period is None:
            self._substep_period = float(day_length_seconds) / float(self.cfg.substeps_per_day)
            self._substep_accum = 0.0
        self._substep_accum += float(dt_seconds)
        if self._substep_accum < self._substep_period:
            return
        # Consume one substep
        self._substep_accum -= self._substep_period

        # Compute banded irradiance [NB, H, W]
        try:
            I_b = dual_star_insolation_to_bands(isr_A, isr_B, eco_adapter.bands)
        except Exception:
            return
        # Gather per sampled cell: [n_cells, NB]
        I_b_cells = I_b[:, self.sample_j, self.sample_i].T  # shape (n_cells, NB)
        # Repeat per individual
        I_b_indiv = I_b_cells[self.indiv_cell_index, :]  # [N, NB]
        # Energy increment: dot(Ab_i, I_b_cell) * substep_dt
        dE = np.einsum("ij,ij->i", self.indiv_Ab, I_b_indiv) * float(self._substep_period)
        self.indiv_E_day += np.maximum(0.0, dE)

        # Water stress: accumulate fractional days if soil < tol
        if soil_W_land is None:
            soil_idx = np.zeros((self.h, self.w), dtype=float)
        elif np.isscalar(soil_W_land):
            soil_idx = np.full((self.h, self.w), float(soil_W_land))
        else:
            soil_idx = np.asarray(soil_W_land, dtype=float)
            if soil_idx.shape != (self.h, self.w):
                soil_idx = np.full((self.h, self.w), float(np.nanmean(soil_idx)))
        # Compute soil index for sampled cells and normalize by a cap (handled upstream)
        soil_cells = soil_idx[self.sample_j, self.sample_i]  # [n_cells]
        soil_indiv = soil_cells[self.indiv_cell_index]       # [N]
        mask_stress = (soil_indiv < self.indiv_tol)
        # accumulate in 'days' units
        self.indiv_water_stress_days[mask_stress] += float(self._substep_period) / float(day_length_seconds)

    def step_daily(self,
                   eco_adapter,
                   soil_W_land: np.ndarray | float,
                   Ts_map: np.ndarray | None = None,
                   day_length_hours: float = 24.0) -> None:
        """
        End-of-day update:
          - Convert per-individual E_day to per-cell species weights (energy proxy)
          - Reset E_day and decay/clear stress days if soil >= tol
          - Apply species weights to PopulationManager by adjusting LAI_layers_SK splits
        """
        pop = getattr(eco_adapter, "pop", None)
        if pop is None or getattr(pop, "LAI_layers_SK", None) is None:
            # Nothing to do
            self.indiv_E_day[:] = 0.0
            self.indiv_water_stress_days[:] = 0.0
            return

        # Build per-cell species energy sum
        N = self.n_indiv
        S = self.ns
        C = self.n_cells
        E = self.indiv_E_day
        sp = self.indiv_species_id
        cell = self.indiv_cell_index

        E_s_c = np.zeros((S, C), dtype=float)
        # accumulate by species and cell using numpy.add.at
        np.add.at(E_s_c, (sp, cell), E)

        # Normalize to weights per cell
        denom = np.sum(E_s_c, axis=0) + 1e-12
        W_s_c = E_s_c / denom[None, :]

        # Optional: include stress penalty (downweight species whose individuals suffered long stress)
        # Here we compute mean stress per (s,c) and apply mild penalty factor
        stress_penalty = float(os.getenv("QD_ECO_INDIV_STRESS_PENALTY", "0.2"))
        if stress_penalty > 0.0:
            stress_s_c = np.zeros((S, C), dtype=float)
            # accumulate stress days per species/cell
            np.add.at(stress_s_c, (sp, cell), self.indiv_water_stress_days)
            cnt_s_c = np.zeros((S, C), dtype=float)
            np.add.at(cnt_s_c, (sp, cell), 1.0)
            # Safe division to avoid RuntimeWarning: compute only where counts > 0
            mean_stress = np.divide(stress_s_c, cnt_s_c, out=np.zeros_like(stress_s_c), where=(cnt_s_c > 0))
            # penalty = 1 / (1 + a * mean_stress), clip
            pen = 1.0 / (1.0 + stress_penalty * mean_stress)
            W_s_c = W_s_c * pen
            # renormalize
            denom2 = np.sum(W_s_c, axis=0) + 1e-12
            W_s_c = W_s_c / denom2[None, :]

        # Apply to PopulationManager: redistribute LAI_layers_SK species splits for sampled cells
        LAI_SK = np.maximum(pop.LAI_layers_SK, 0.0)  # [S,K,H,W]
        K = int(getattr(pop, "K", 1))
        H, W = self.h, self.w

        # --- Simple LAI growth/decay + very local dispersal (daily) ---
        # Energy scaling reference (median of sampled-cell total energy)
        medE = float(np.median(denom[denom > 0])) if np.any(denom > 0) else 1.0
        # Tunable daily rates (dimensionless per day)
        lai_grow = float(os.getenv("QD_ECO_LAI_GROWTH_RATE", "0.002"))   # growth when energy above median
        lai_decay = float(os.getenv("QD_ECO_LAI_DECAY_RATE", "0.001"))   # decay with stress
        recruit_frac = float(os.getenv("QD_ECO_LAI_RECRUIT_FRAC", "0.2"))  # fraction of positive dLAI that spills to neighbors

        # For each sampled cell, compute total per-layer LAI and split across species following W_s_c
        for ci in range(C):
            j = int(self.sample_j[ci])
            i = int(self.sample_i[ci])
            # total LAI per layer at this cell
            total_k = np.sum(LAI_SK[:, :, j, i], axis=0)  # [K]
            wk = W_s_c[:, ci]  # [S]

            # --- magnitude update (growth/decay) ---
            total_old = float(np.sum(total_k))
            # energy-based growth driver (relative to median)
            E_cell = float(denom[ci])
            e_scaled = E_cell / (medE + 1e-12)
            # optional cell mean stress (if computed above)
            try:
                mean_stress_cell = float(np.sum(mean_stress[:, ci] * wk)) if 'mean_stress' in locals() else 0.0
            except Exception:
                mean_stress_cell = 0.0
            dLAI = lai_grow * (e_scaled - 1.0) - lai_decay * mean_stress_cell
            # apply as absolute change (per day), scaled by current canopy magnitude for sensitivity
            dLAI *= max(total_old, 1.0)
            # clamp to physical limits
            lai_max_cell = float(getattr(pop.params, "lai_max", 5.0))
            new_total = float(np.clip(total_old + dLAI, 0.0, lai_max_cell))
            scale = (new_total / (total_old + 1e-12)) if total_old > 0.0 else (new_total / max(lai_max_cell, 1.0))

            # --- rebuild SK with updated magnitude and species split wk ---
            for k in range(K):
                new_k = total_k[k] * scale
                if new_k <= 0.0:
                    LAI_SK[:, k, j, i] = 0.0
                else:
                    LAI_SK[:, k, j, i] = wk * new_k

            # --- very local dispersal: spill a portion of positive growth to 4-neighbors (land only) ---
            recruit = max(0.0, new_total - total_old) * recruit_frac
            if recruit > 0.0:
                # 4-neighbors with periodic longitude
                jn = [max(0, j - 1), min(H - 1, j + 1), j, j]
                in_ = [(i - 1) % W, (i + 1) % W, i, i]
                share = recruit / 4.0
                for (jj, ii) in zip(jn, in_):
                    # skip non-land cells by checking a minimal canopy allowance
                    # (PopulationManager 不直接持有 land_mask，这里用 LAI 总量为 0 近似海洋/不适生)
                    tot_nb = float(np.sum(LAI_SK[:, :, jj, ii]))
                    # Distribute equally over layers with current species weights wk
                    add_each_layer = share / float(max(K, 1))
                    for k in range(K):
                        LAI_SK[:, k, jj, ii] += wk * add_each_layer

        pop.LAI_layers_SK = np.clip(LAI_SK, 0.0, pop.params.lai_max)
        pop.LAI_layers = np.sum(pop.LAI_layers_SK, axis=0)
        pop.LAI = np.sum(pop.LAI_layers, axis=0)
        # Update species weights (global) from area-summed LAI (just for adapter bookkeeping)
        pop.recompute_species_weights_from_LAI()

        # M3+: Couple per-individual reproduction → PopulationManager.seed_bank
        try:
            if int(os.getenv("QD_ECO_INDIV_SEED_COUPLE", "1")) == 1 and hasattr(pop, "seed_bank"):
                # per-cell total energy over the day
                E_tot_cells = denom  # [C]
                repro_frac = float(getattr(pop, "repro_fraction", float(os.getenv("QD_ECO_REPRO_FRACTION", "0.2"))))
                seed_energy = float(getattr(pop, "seed_energy", float(os.getenv("QD_ECO_SEED_ENERGY", "1.0"))))
                retain = float(os.getenv("QD_ECO_SEED_BANK_RETAIN", "0.2"))
                bank_max = float(os.getenv("QD_ECO_SEED_BANK_MAX", "1000.0"))
                # soil gating (reuse soil_idx if available)
                soil_cells = None
                if 'soil_idx' in locals():
                    soil_cells = soil_idx[self.sample_j, self.sample_i]
                # Seeds produced per cell (retained portion)
                seeds_cells = np.maximum(0.0, repro_frac) * np.maximum(0.0, E_tot_cells) / max(seed_energy, 1e-12)
                if soil_cells is not None:
                    seeds_cells = seeds_cells * np.clip(soil_cells, 0.0, 1.0)
                seeds_cells = retain * seeds_cells
                # Scatter-add into pop.seed_bank at sampled locations
                for ci in range(C):
                    j = int(self.sample_j[ci]); i = int(self.sample_i[ci])
                    pop.seed_bank[j, i] = np.clip(pop.seed_bank[j, i] + seeds_cells[ci], 0.0, bank_max)
        except Exception:
            pass

        # Reset daily buffers (energy and, partially, stress)
        self.indiv_E_day[:] = 0.0
        # Stress relief: if soil >= tol for an individual, reduce counter
        if soil_W_land is None:
            soil_idx = np.zeros((H, W), dtype=float)
        elif np.isscalar(soil_W_land):
            soil_idx = np.full((H, W), float(soil_W_land))
        else:
            soil_idx = np.asarray(soil_W_land, dtype=float)
            if soil_idx.shape != (H, W):
                soil_idx = np.full((H, W), float(np.nanmean(soil_idx)))
        soil_cells = soil_idx[self.sample_j, self.sample_i]
        soil_indiv = soil_cells[self.indiv_cell_index]
        ok = (soil_indiv >= self.indiv_tol)
        # decay factor
        decay = float(os.getenv("QD_ECO_INDIV_STRESS_DECAY", "0.5"))
        self.indiv_water_stress_days[ok] *= decay
        self.indiv_water_stress_days[~ok] = np.minimum(self.indiv_water_stress_days[~ok] + 1.0, 365.0)

        if self.cfg.diag:
            beta_hint = float(np.mean(np.max(W_s_c, axis=0)))
            print(f"[EcoIndiv] daily applied to {self.n_cells} cells × {self.per_cell} indiv; "
                  f"mean max species share per cell ~ {beta_hint:.2f} (lower→more even).")

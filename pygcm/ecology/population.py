from __future__ import annotations

import os
import numpy as np
from dataclasses import dataclass
from .spectral import absorbance_from_genes


@dataclass
class LAIParams:
    """Simple LAI prognostic parameters (M2 minimal)."""
    lai_max: float = 5.0                 # maximum LAI
    k_canopy: float = 0.5                # Beer-Lambert like coefficient
    growth_per_j: float = 2.0e-5         # LAI growth per unit "J-equivalent" daily energy
    senesce_per_day: float = 0.01        # daily senescence under stress
    stress_thresh: float = 0.3           # soil water index threshold (0..1)
    stress_strength: float = 1.0         # scaling for senescence when below threshold

    @staticmethod
    def from_env() -> "LAIParams":
        def f(name: str, default: float) -> float:
            try:
                return float(os.getenv(name, str(default)))
            except Exception:
                return default
        return LAIParams(
            lai_max=f("QD_ECO_LAI_MAX", 5.0),
            k_canopy=f("QD_ECO_LAI_K", 0.5),
            growth_per_j=f("QD_ECO_LAI_GROWTH", 2.0e-5),
            senesce_per_day=f("QD_ECO_LAI_SENESCENCE", 0.01),
            stress_thresh=f("QD_ECO_SOIL_STRESS_THRESH", 0.3),
            stress_strength=f("QD_ECO_SOIL_STRESS_GAIN", 1.0),
        )


class PopulationManager:
    """
    Minimal M2 population manager per-grid:
    - Prognostic LAI[lat,lon] on land; initialized small but >0
    - Subdaily: accumulate daily energy buffer E_day (proxy from ISR)
    - Daily: update LAI using growth from E_day and senescence under soil water stress
    - Provide canopy reflectance factor f(LAI) used by adapter to synthesize alpha_map

    Notes:
    - Units are normalized proxies (J-equivalents); growth_per_j should be tuned empirically.
    - soil_water_index is expected in [0,1] (0=dry, 1=wet).
    """
    def __init__(self, land_mask: np.ndarray, *, diag: bool = True):
        self.land = (land_mask == 1)
        self.shape = land_mask.shape
        self.params = LAIParams.from_env()
        # Initialize LAI: small positive over land
        self.LAI = np.zeros(self.shape, dtype=float)
        self.LAI[self.land] = float(os.getenv("QD_ECO_LAI_INIT", "0.2"))
        # Daily energy buffer (proxy J-equivalents)
        self.E_day = np.zeros(self.shape, dtype=float)
        self._diag = diag
        # --- M2: canopy cache and recompute policy ---
        self._hours_accum: float = 0.0
        try:
            self._light_update_every_hours = float(os.getenv("QD_ECO_LIGHT_UPDATE_EVERY_HOURS", "6"))
        except Exception:
            self._light_update_every_hours = 6.0
        try:
            self._lai_recompute_delta = float(os.getenv("QD_ECO_LIGHT_RECOMPUTE_LAI_DELTA", "0.05"))
        except Exception:
            self._lai_recompute_delta = 0.05
        self._canopy_f_cached: np.ndarray | None = None
        self._lai_snapshot: np.ndarray = self.total_LAI().copy()
        self._next_recompute_hours: float = self._light_update_every_hours
        # --- M1: genotype absorbance cache (genes → A_b) ---
        self._genotype_absorb_cache: dict[str, np.ndarray] = {}

        # Cohort layers (K) scaffold: optional vertical layering; extend to species S×K
        try:
            self.K = max(1, int(os.getenv("QD_ECO_COHORT_K", "1")))
        except Exception:
            self.K = 1
        # M4: species (genes) mixture support (weights sum to 1, default 1 species)
        # Parsed by adapter to compute banded leaf reflectance; here only store weights.
        species_weights_env = os.getenv("QD_ECO_SPECIES_WEIGHTS", "").strip()
        # Track whether weights come from env (affects default spread-mode assignment policy)
        self._weights_from_env = bool(species_weights_env)
        if species_weights_env:
            try:
                w = [float(x) for x in species_weights_env.split(",") if x.strip() != ""]
            except Exception:
                w = [1.0]
        else:
            try:
                ns_default = max(1, int(os.getenv("QD_ECO_NS", "20")))
            except Exception:
                ns_default = 20
            # default to Ns=ns_default equal weights
            w = [1.0 / float(ns_default)] * ns_default
        s = sum(w) if len(w) > 0 else 1.0
        self.species_weights = np.asarray([max(0.0, wi) for wi in w], dtype=float)
        if s <= 0:
            try:
                ns_default = max(1, int(os.getenv("QD_ECO_NS", "20")))
            except Exception:
                ns_default = 20
            self.species_weights = np.full((ns_default,), 1.0 / float(ns_default), dtype=float)
        else:
            self.species_weights /= s

        # Species count (after species_weights is defined)
        try:
            self.Ns = int(self.species_weights.shape[0])
        except Exception:
            self.Ns = 1

        # LAI_layers_SK shape [S, K, lat, lon]: initialize by species weight×equal split over K
        self.LAI_layers_SK = np.zeros((self.Ns, self.K, *self.shape), dtype=float)
        for s_idx in range(self.Ns):
            frac_s = float(self.species_weights[s_idx]) if self.Ns > 0 else 1.0
            for k in range(self.K):
                self.LAI_layers_SK[s_idx, k, :, :] = frac_s * (self.LAI / float(self.K))
        # Backward-compat aggregated [K,lat,lon]
        self.LAI_layers = np.sum(self.LAI_layers_SK, axis=0)

        # Species leaf reflectance cache per band (filled by caller via set_species_reflectance_bands)
        self._species_R_leaf = None  # shape [Ns, NB]

        # --- Spatial spread (colonization) controls (M2.5 minimal) ---
        # Enable a conservative neighbor-exchange ("diffusion") of LAI over land to emulate
        # vegetative expansion/propagule dispersal at grid scale.
        # Env:
        #   QD_ECO_SPREAD_ENABLE  (0/1, default 0)
        #   QD_ECO_SPREAD_RATE    (per-day fraction in 0..0.5, default 0.0)
        #   QD_ECO_SPREAD_NEIGHBORS ('vonNeumann' or 'moore', default 'vonNeumann')
        try:
            self.spread_enable = int(os.getenv("QD_ECO_SPREAD_ENABLE", "0")) == 1
        except Exception:
            self.spread_enable = False
        try:
            self.spread_rate = float(os.getenv("QD_ECO_SPREAD_RATE", "0.0"))
        except Exception:
            self.spread_rate = 0.0
        try:
            self.spread_neighbors = os.getenv("QD_ECO_SPREAD_NEIGHBORS", "vonNeumann").strip().lower()
        except Exception:
            self.spread_neighbors = "vonneumann"

        # Seed-based dispersal parameters (optional, ties speed to seed investment)
        try:
            self.spread_mode = os.getenv("QD_ECO_SPREAD_MODE", "diffusion").strip().lower()  # 'diffusion' | 'seed'
        except Exception:
            self.spread_mode = "diffusion"
        try:
            self.repro_fraction = float(os.getenv("QD_ECO_REPRO_FRACTION", "0.2"))  # fraction of daily energy to reproduction
        except Exception:
            self.repro_fraction = 0.2
        try:
            self.seed_energy = float(os.getenv("QD_ECO_SEED_ENERGY", "1.0"))  # energy per seed
        except Exception:
            self.seed_energy = 1.0
        try:
            self.seed_scale = float(os.getenv("QD_ECO_SEED_SCALE", "1.0"))   # scaling in r_eff = r0*(1-exp(-seeds/seed_scale))
        except Exception:
            self.seed_scale = 1.0
        try:
            self.seedling_lai = float(os.getenv("QD_ECO_SEEDLING_LAI", "0.02"))  # LAI boost per established seed unit
        except Exception:
            self.seedling_lai = 0.02
        # Age map (days since establishment); starts at 0, increments daily where LAI>0
        self.age_days = np.zeros(self.shape, dtype=float)
        # Seed bank (M3): per-cell seed reservoir (arb. units), retained/germinated daily
        self.seed_bank = np.zeros(self.shape, dtype=float)
        # Spread gating map (e.g., from soil moisture); initialize to land mask
        self._spread_gate = self.land.astype(float)

        # Per-species spread modes (default: species 0 'diffusion' (grass), species 1 'seed' (tree) if present)
        self.species_modes: list[str] = []
        self._init_species_modes()

    def _init_species_modes(self) -> None:
        """
        Initialize per-species spread modes with the following policy:
        - If QD_ECO_SPECIES_{i}_MODE is provided for a species i, honor it.
        - Else if weights are provided (QD_ECO_SPECIES_WEIGHTS present), draw ONE species index
          from the given weight distribution to be 'seed' (tree), mark the rest 'diffusion' (grass).
        - Else (no weights provided), assign each unspecified species independently at random with
          equal probability to 'seed' or 'diffusion'.
        A fixed RNG seed can be provided via QD_ECO_RAND_SEED for reproducibility.
        """
        try:
            S = int(getattr(self, "Ns", 1))
        except Exception:
            S = 1
        # Start with any explicit per-species modes from env
        explicit = []
        for i in range(S):
            m = os.getenv(f"QD_ECO_SPECIES_{i}_MODE", "").strip().lower() if os.getenv else ""
            explicit.append(m if m in ("seed", "diffusion") else "")
        modes = [""] * S
        # Fill explicit first
        for i in range(S):
            if explicit[i]:
                modes[i] = explicit[i]
        # RNG
        try:
            seed_val = os.getenv("QD_ECO_RAND_SEED")
            rng = np.random.default_rng(int(seed_val)) if seed_val not in (None, "",) else np.random.default_rng()
        except Exception:
            rng = np.random.default_rng()
        # Unspecified indices
        unspec_idx = [i for i in range(S) if modes[i] == ""]
        if len(unspec_idx) == 0:
            self.species_modes = modes
            return
        # If weights come from env: choose exactly one 'seed' species by weights, rest 'diffusion'
        weights_from_env = bool(getattr(self, "_weights_from_env", False))
        if weights_from_env:
            try:
                w = np.asarray(getattr(self, "species_weights", np.ones((S,), dtype=float)), dtype=float)
                w = np.clip(w, 0.0, None)
                w = w / (np.sum(w) + 1e-12)
                # Draw one index globally across all species by weights
                chosen = int(rng.choice(np.arange(S), p=w))
            except Exception:
                chosen = 1 if S > 1 else 0
            for i in unspec_idx:
                modes[i] = "seed" if i == chosen else "diffusion"
        else:
            # No weights from env: assign per-species uniformly at random
            for i in unspec_idx:
                modes[i] = "seed" if rng.random() < 0.5 else "diffusion"
        self.species_modes = modes

    def set_species_modes(self, modes: list[str]) -> None:
        """
        Set per-species spread modes; list length should be Ns.
        Valid entries: 'seed' | 'diffusion'. Missing entries fallback to defaults.
        """
        try:
            S = int(getattr(self, "Ns", 1))
        except Exception:
            S = 1
        out = []
        for i in range(S):
            if i < len(modes) and str(modes[i]).lower() in ("seed", "diffusion"):
                out.append(str(modes[i]).lower())
            else:
                # fallback same as init
                if i == 1:
                    out.append("seed")
                else:
                    out.append("diffusion")
        self.species_modes = out

    def step_subdaily(self, isr_total: np.ndarray, dt_seconds: float, *, return_bands: bool = False, soil_ref: float = 0.20) -> np.ndarray | None:
        """
        Accumulate daily energy buffer from incoming shortwave proxy.
        Use simple proportionality: dE = isr_total * dt.

        When return_bands=True, also返回 A_b^surface（NB×lat×lon），NB 取自物种反射缓存的列数；
        若未设置物种反射缓存，则返回 None（由上层按需调用 get_surface_albedo_bands）。
        """
        if isr_total is None:
            return None
        # Ensure shapes match
        if isr_total.shape != self.shape:
            isr = np.full(self.shape, float(np.nan))
            isr[:] = float(np.nanmean(isr_total))
        else:
            isr = isr_total
        dE = np.nan_to_num(isr) * float(dt_seconds)
        self.E_day += dE

        # Update canopy cache clock and recompute policy
        self._hours_accum += float(dt_seconds) / 3600.0
        if self._should_recompute_canopy():
            self._recompute_canopy_cache()
            self._lai_snapshot = self.total_LAI().copy()
            self._next_recompute_hours = self._hours_accum + self._light_update_every_hours

        if return_bands:
            if getattr(self, "_species_R_leaf", None) is None:
                return None
            nb = int(self._species_R_leaf.shape[1])
            try:
                return self.get_surface_albedo_bands(nb, soil_ref=soil_ref)
            except Exception:
                return None
        return None

    def total_LAI(self) -> np.ndarray:
        """Return total LAI (sum over species and layers)."""
        if getattr(self, "LAI_layers_SK", None) is not None:
            return np.sum(self.LAI_layers_SK, axis=(0, 1))
        if getattr(self, "LAI_layers", None) is not None:
            return np.sum(self.LAI_layers, axis=0)
        return self.LAI

    def canopy_height_map(self) -> np.ndarray:
        """
        Simple canopy height proxy (m) from layered LAI:
          H = H_scale * Σ_k (h_k * LAI_k) / Σ_k LAI_k,  h_k = (k+1)/K
        """
        try:
            H_scale = float(os.getenv("QD_ECO_HEIGHT_SCALE_M", "10.0"))
        except Exception:
            H_scale = 10.0
        if getattr(self, "LAI_layers_SK", None) is None:
            # Single layer: height ~ H_scale * f(LAI)
            f = 1.0 - np.exp(-self.params.k_canopy * np.maximum(self.LAI, 0.0))
            H = H_scale * f
        else:
            K = self.K
            idx = np.arange(1, K + 1, dtype=float)[None, :, None, None] / float(K)  # [1,K,1,1]
            LAI_layers_pos = np.maximum(self.LAI_layers_SK, 0.0)  # [S,K,lat,lon]
            LAI_by_k = np.sum(LAI_layers_pos, axis=0)  # [K,lat,lon]
            num = np.sum(idx[0] * LAI_by_k, axis=0)
            den = np.sum(LAI_by_k, axis=0) + 1e-12
            H = H_scale * (num / den)
        # Land-only map, NaN over ocean
        out = np.full(self.shape, np.nan, dtype=float)
        out[self.land] = H[self.land]
        return out

    def species_density_maps(self) -> list[np.ndarray]:
        """
        Return per-species density maps as Σ_k LAI_s,k (land-only).
        """
        maps = []
        if getattr(self, "LAI_layers_SK", None) is None:
            # Fallback proxy
            Ltot = self.total_LAI()
            for wi in np.atleast_1d(self.species_weights):
                m = np.full(self.shape, np.nan, dtype=float)
                m[self.land] = (wi * Ltot[self.land])
                maps.append(m)
            return maps
        S = self.Ns
        L_s = np.sum(np.maximum(self.LAI_layers_SK, 0.0), axis=1)  # [S,lat,lon]
        for s in range(S):
            m = np.full(self.shape, np.nan, dtype=float)
            m[self.land] = L_s[s, :, :][self.land]
            maps.append(m)
        return maps

    def recompute_species_weights_from_LAI(self) -> None:
        """
        Reset species_weights by normalizing area-summed LAI per species.
        (Unweighted sum over land; sufficient for reflectance mixing.)
        """
        if getattr(self, "LAI_layers_SK", None) is None:
            return
        S = self.Ns
        L_s = np.sum(np.maximum(self.LAI_layers_SK, 0.0), axis=1)  # [S,lat,lon]
        totals = np.zeros((S,), dtype=float)
        for s in range(S):
            totals[s] = float(np.nansum(L_s[s, :, :][self.land]))
        ssum = float(np.sum(totals))
        if ssum <= 0:
            self.species_weights = np.full((S,), 1.0 / float(S), dtype=float)
        else:
            self.species_weights = np.clip(totals / ssum, 0.0, 1.0)

    def add_species_from_parent(self, parent_idx: int, frac: float = 0.02) -> int:
        """
        Split a fraction of parent species LAI into a new species across all layers (conservative).
        Returns the index of the new species.
        """
        if getattr(self, "LAI_layers_SK", None) is None:
            return 0
        S_old, K, H, W = self.LAI_layers_SK.shape
        p = int(np.clip(parent_idx, 0, S_old - 1))
        f = float(np.clip(frac, 0.0, 0.5))
        if f <= 0.0:
            return p
        # Allocate new array
        new = np.zeros((S_old + 1, K, H, W), dtype=float)
        new[:S_old, :, :, :] = self.LAI_layers_SK
        # Transfer fraction from parent
        transfer = f * self.LAI_layers_SK[p, :, :, :]
        new[p, :, :, :] = self.LAI_layers_SK[p, :, :, :] - transfer
        new[S_old, :, :, :] = transfer
        self.LAI_layers_SK = np.clip(new, 0.0, self.params.lai_max)
        self.Ns = S_old + 1
        # Refresh aggregated views
        self.LAI_layers = np.sum(self.LAI_layers_SK, axis=0)
        self.LAI = np.sum(self.LAI_layers, axis=0)
        # Recompute species weights from LAI
        self.recompute_species_weights_from_LAI()
        return self.Ns - 1

    def step_daily(self, soil_water_index: np.ndarray | float | None) -> None:
        """
        Daily LAI update:
        LAI_next = LAI + growth(E_day) - senescence(stress)
        - growth = growth_per_j * E_day (land only), saturating at lai_max
        - stress = max(0, thresh - soil) * senesce_per_day * gain
        Reset E_day at the end.
        """
        P = self.params
        land = self.land

        # Growth term from daily energy (split into growth vs reproduction globally by repro_fraction)
        E_map = self.E_day
        repro_frac = float(np.clip(getattr(self, "repro_fraction", 0.0), 0.0, 0.95))
        growth_energy = (1.0 - repro_frac) * np.nan_to_num(E_map)
        growth = P.growth_per_j * growth_energy
        growth = np.where(land, growth, 0.0)

        # Soil water stress term
        if soil_water_index is None:
            soil = np.zeros(self.shape, dtype=float)
        elif np.isscalar(soil_water_index):
            soil = np.full(self.shape, float(soil_water_index))
        else:
            soil = np.asarray(soil_water_index, dtype=float)
            if soil.shape != self.shape:
                # fallback to mean
                soil = np.full(self.shape, float(np.nanmean(soil)))

        stress = np.maximum(0.0, P.stress_thresh - np.clip(soil, 0.0, 1.0))
        sen = P.senesce_per_day * P.stress_strength * stress
        sen = np.where(land, sen, 0.0)

        # Build spread gating map from soil (optional) for today's spread step
        try:
            if int(os.getenv("QD_ECO_SPREAD_GATE_SOIL", "1")) == 1:
                exp = float(os.getenv("QD_ECO_SPREAD_SOIL_EXP", "1.0"))
                gate = np.clip(soil, 0.0, 1.0) ** exp
                self._spread_gate = np.where(land, gate, 0.0)
            else:
                self._spread_gate = self.land.astype(float)
        except Exception:
            self._spread_gate = self.land.astype(float)

        # Layered Beer-Lambert allocation (top-down) using daily energy proxy
        K = int(getattr(self, "K", 1))
        if K > 1 and getattr(self, "LAI_layers_SK", None) is not None:
            # 1) Layer capture for total canopy (sum over species)
            I_in = np.nan_to_num(self.E_day)  # scalar "light"
            cap_k = np.zeros((K, *self.shape), dtype=float)
            LAI_k_total = np.sum(np.maximum(self.LAI_layers_SK, 0.0), axis=0)  # [K,lat,lon]
            for k in range(K):
                T_k = np.exp(-P.k_canopy * LAI_k_total[k, :, :])
                cap_k[k, :, :] = I_in * (1.0 - T_k)
                I_in = I_in * T_k
            cap_sum = np.sum(cap_k, axis=0)  # [lat,lon]

            # 2) Distribute growth to species×layers by LAI share within each layer
            growth_total = growth  # [lat,lon]
            growth_layers_SK = np.zeros_like(self.LAI_layers_SK)
            # weights within layer
            LAI_prev_SK = np.maximum(self.LAI_layers_SK, 0.0)  # [S,K,lat,lon]
            LAI_prev_by_k = np.sum(LAI_prev_SK, axis=0)  # [K,lat,lon]
            with np.errstate(invalid="ignore", divide="ignore"):
                w_s_k = np.where(LAI_prev_by_k[None, :, :, :] > 0.0,
                                 LAI_prev_SK / (LAI_prev_by_k[None, :, :, :] + 1e-12),
                                 1.0 / float(self.Ns))
            # growth split across layers first by cap_k, then within layer by species share
            with np.errstate(invalid="ignore", divide="ignore"):
                wcap_k = cap_k / (cap_sum[None, :, :] + 1e-12)  # [K,lat,lon]
            no_cap = (cap_sum <= 0.0)
            has_cap = ~no_cap
            # no capture → equal split across K and species
            if np.any(no_cap):
                eq = (growth_total[no_cap] / float(K) / float(self.Ns))
                for s in range(self.Ns):
                    for k in range(K):
                        growth_layers_SK[s, k, no_cap] = eq
            if np.any(has_cap):
                for s in range(self.Ns):
                    for k in range(K):
                        growth_layers_SK[s, k, has_cap] = w_s_k[s, k, has_cap] * wcap_k[k, has_cap] * growth_total[has_cap]

            # 3) Senescence per species proportional to current LAI share
            LAI_tot_prev = np.sum(LAI_prev_SK, axis=(0, 1))
            with np.errstate(invalid="ignore", divide="ignore"):
                wsen_s_k = np.where(LAI_tot_prev[None, None, :, :] > 0.0,
                                    LAI_prev_SK / (LAI_tot_prev[None, None, :, :] + 1e-12),
                                    1.0 / float(self.Ns * K))
            sen_layers_SK = wsen_s_k * sen[None, None, :, :]

            # 4) Update SK layers and clamp
            self.LAI_layers_SK = np.clip(LAI_prev_SK + growth_layers_SK - sen_layers_SK, 0.0, P.lai_max)

            # 5) Upward transfer species-wise
            try:
                upfrac = float(os.getenv("QD_ECO_LAYER_UPFRAC", "0.1"))
            except Exception:
                upfrac = 0.1
            if upfrac > 0.0:
                for s in range(self.Ns):
                    for k in range(K - 1, 0, -1):
                        excess = np.maximum(0.0, self.LAI_layers_SK[s, k, :, :] - self.LAI_layers_SK[s, k - 1, :, :])
                        delta = upfrac * excess
                        self.LAI_layers_SK[s, k, :, :] -= delta
                        self.LAI_layers_SK[s, k - 1, :, :] += delta

            # Refresh aggregates after SK update
            self.LAI_layers = np.sum(self.LAI_layers_SK, axis=0)
            self.LAI = np.sum(self.LAI_layers, axis=0)
        else:
            # Fallback single-layer update
            self.LAI = np.clip(self.LAI + growth - sen, 0.0, P.lai_max)

        # Optional spatial spread after growth/senescence
        seeded_mask = None
        if getattr(self, "spread_enable", False) and float(getattr(self, "spread_rate", 0.0)) > 0.0:
            try:
                S = int(getattr(self, "Ns", 1))
                modes = getattr(self, "species_modes", [])
                any_seed = False
                for s_idx in range(S):
                    mode_s = (modes[s_idx] if s_idx < len(modes) else ("seed" if s_idx == 1 else "diffusion"))
                    if mode_s == "seed":
                        m = self._seed_based_spread_species(s_idx)
                        if m is not None:
                            any_seed = True
                    else:
                        self._apply_neighbor_spread_species(s_idx, rate=float(getattr(self, "spread_rate", 0.0)))
                if self._diag and np.random.rand() < 0.05:
                    sdiag = self.summary()
                    print(f"[Ecology] per-species spread modes: {modes} | r0={self.spread_rate:.3f}/day | "
                          f"LAI(min/mean/max)={sdiag['LAI_min']:.2f}/{sdiag['LAI_mean']:.2f}/{sdiag['LAI_max']:.2f}")
                if any_seed:
                    # Mark seeded today at cell-level for age reset below
                    # (approximation: where any species seeded)
                    seeded_mask = np.zeros(self.shape, dtype=bool)
                    S = int(getattr(self, "Ns", 1))
                    # Reconstruct from today's additions: we don't track per-call adds, so fallback to heuristic:
                    # Use small threshold on today's LAI increase could be implemented; here simply keep non-None flag.
                    # Keeping seeded_mask as None would skip the age reset special casing; instead, we reset age where total LAI increased.
                    # As a simple proxy, leave seeded_mask=None and rely on age reset policy below to use inc mask when available.
                    seeded_mask = None
            except Exception:
                pass

        # Age update: increment by 1 day where LAI>0, but keep newly seeded pixels at 0 today
        try:
            land_mask = self.land
            has_lai = (np.maximum(self.total_LAI(), 0.0) > 0.0) & land_mask
            if seeded_mask is None:
                self.age_days[has_lai] += 1.0
            else:
                inc_mask = has_lai & (~seeded_mask)
                self.age_days[inc_mask] += 1.0
        except Exception:
            pass

        # Germination & decay of seed bank (M3)
        try:
            germ_frac = float(os.getenv("QD_ECO_SEED_GERMINATE_FRAC", "0.10"))
        except Exception:
            germ_frac = 0.10
        try:
            decay = float(os.getenv("QD_ECO_SEED_BANK_DECAY", "0.02"))
        except Exception:
            decay = 0.02
        # Soil-gated germination
        gate = getattr(self, "_spread_gate", None)
        if gate is None:
            gate = self.land.astype(float)
        else:
            gate = np.where(self.land, np.clip(gate, 0.0, 1.0), 0.0)
        seeds_to_germ = np.maximum(0.0, germ_frac) * self.seed_bank * gate
        # Convert to seedling LAI at lowest layer, per-species by weights
        try:
            s_lai = float(os.getenv("QD_ECO_SEEDLING_LAI", "0.02"))
        except Exception:
            s_lai = 0.02
        if getattr(self, "LAI_layers_SK", None) is not None:
            try:
                S = int(getattr(self, "Ns", 1))
                w = np.asarray(getattr(self, "species_weights", np.ones((S,), dtype=float)), dtype=float)
                w = w / (np.sum(w) + 1e-12)
            except Exception:
                S = 1
                w = np.asarray([1.0], dtype=float)
            add_total = s_lai * seeds_to_germ  # [lat,lon]
            for s in range(S):
                add_s = w[s] * add_total
                self.LAI_layers_SK[s, 0, self.land] = np.clip(
                    self.LAI_layers_SK[s, 0, self.land] + add_s[self.land],
                    0.0, self.params.lai_max
                )
            # refresh aggregates
            self.LAI_layers = np.sum(self.LAI_layers_SK, axis=0)
            self.LAI = np.sum(self.LAI_layers, axis=0)
        else:
            self.LAI[self.land] = np.clip(self.LAI[self.land] + s_lai * seeds_to_germ[self.land], 0.0, self.params.lai_max)
        # Seed bank decay after germination
        try:
            self.seed_bank = np.maximum(0.0, self.seed_bank - seeds_to_germ)
            self.seed_bank *= max(0.0, (1.0 - decay))
        except Exception:
            pass

        # Reset daily energy buffer
        self.E_day[:] = 0.0

    def _apply_neighbor_spread(self) -> None:
        """
        Backward-compatible total-LAI diffusion (kept for global mode).
        """
        self._apply_neighbor_spread_species(None, rate=float(getattr(self, "spread_rate", 0.0)))

    def _apply_neighbor_spread_species(self, s_idx: int | None, rate: float) -> None:
        """
        Per-species conservative neighbor exchange.
        - If s_idx is None: apply to total LAI (legacy).
        - Else: apply to species s only, scaling its layers by per-cell factor.
        """
        rate = float(max(0.0, min(0.5, rate)))
        if rate <= 0.0:
            return
        land = self.land
        neigh_mode = str(getattr(self, "spread_neighbors", "vonneumann")).lower()
        if neigh_mode in ("moore", "8", "8n"):
            offsets = [(-1,-1), (-1,0), (-1,1),
                       ( 0,-1),          ( 0,1),
                       ( 1,-1), ( 1,0),  ( 1,1)]
        else:
            offsets = [(-1,0), (0,-1), (0,1), (1,0)]

        if s_idx is None or getattr(self, "LAI_layers_SK", None) is None:
            # total-LAI path
            LAI_prev = self.total_LAI()
            H, W = LAI_prev.shape
            num_valid = np.zeros((H, W), dtype=float)
            for dy, dx in offsets:
                num_valid += np.roll(land, shift=(-dy, -dx), axis=(0, 1)).astype(float)
            gate = getattr(self, "_spread_gate", None)
            if gate is None:
                gate = land.astype(float)
            else:
                gate = np.where(land, np.clip(gate, 0.0, 1.0), 0.0)
            outflow = rate * np.maximum(LAI_prev, 0.0) * gate
            with np.errstate(invalid="ignore", divide="ignore"):
                share = np.where(num_valid > 0.0, outflow / (num_valid + 1e-12), 0.0)
            inflow = np.zeros((H, W), dtype=float)
            for dy, dx in offsets:
                inflow += np.roll(share, shift=(dy, dx), axis=(0, 1))
            LAI_new = np.full_like(LAI_prev, 0.0)
            LAI_raw = LAI_prev - outflow + inflow
            # Cap daily positive increment to avoid runaway front
            try:
                dmax = float(os.getenv("QD_ECO_SPREAD_DLAI_MAX", "0.02"))
            except Exception:
                dmax = 0.02
            inc = LAI_raw - LAI_prev
            inc_pos = np.minimum(np.maximum(inc, 0.0), dmax)
            dec = np.minimum(inc, 0.0)
            LAI_capped = LAI_prev + inc_pos + dec
            LAI_new[land] = np.clip(LAI_capped[land], 0.0, self.params.lai_max)
            if getattr(self, "LAI_layers_SK", None) is not None:
                LAI_tot_prev = np.maximum(np.sum(self.LAI_layers_SK, axis=(0, 1)), 0.0)
                with np.errstate(invalid="ignore", divide="ignore"):
                    factor = np.where(LAI_tot_prev > 0.0, LAI_new / (LAI_tot_prev + 1e-12), 0.0)
                self.LAI_layers_SK = np.clip(self.LAI_layers_SK * factor[None, None, :, :], 0.0, self.params.lai_max)
                self.LAI_layers = np.sum(self.LAI_layers_SK, axis=0)
                self.LAI = np.sum(self.LAI_layers, axis=0)
            else:
                self.LAI = np.clip(LAI_new, 0.0, self.params.lai_max)
            return

        # species-only path
        S, K, H, W = self.LAI_layers_SK.shape
        s = int(np.clip(s_idx, 0, S - 1))
        LAI_s_prev = np.maximum(np.sum(self.LAI_layers_SK[s, :, :, :], axis=0), 0.0)  # [H,W]
        num_valid = np.zeros((H, W), dtype=float)
        for dy, dx in offsets:
            num_valid += np.roll(land, shift=(-dy, -dx), axis=(0, 1)).astype(float)
        gate = getattr(self, "_spread_gate", None)
        if gate is None:
            gate = land.astype(float)
        else:
            gate = np.where(land, np.clip(gate, 0.0, 1.0), 0.0)
        outflow = rate * LAI_s_prev * gate
        with np.errstate(invalid="ignore", divide="ignore"):
            share = np.where(num_valid > 0.0, outflow / (num_valid + 1e-12), 0.0)
        inflow = np.zeros((H, W), dtype=float)
        for dy, dx in offsets:
            inflow += np.roll(share, shift=(dy, dx), axis=(0, 1))
        LAI_s_new = np.full_like(LAI_s_prev, 0.0)
        LAI_s_raw = LAI_s_prev - outflow + inflow
        # Cap daily positive increment for species s
        try:
            dmax = float(os.getenv("QD_ECO_SPREAD_DLAI_MAX", "0.02"))
        except Exception:
            dmax = 0.02
        inc = LAI_s_raw - LAI_s_prev
        inc_pos = np.minimum(np.maximum(inc, 0.0), dmax)
        dec = np.minimum(inc, 0.0)
        LAI_s_capped = LAI_s_prev + inc_pos + dec
        LAI_s_new[land] = np.clip(LAI_s_capped[land], 0.0, self.params.lai_max)

        with np.errstate(invalid="ignore", divide="ignore"):
            factor = np.where(LAI_s_prev > 0.0, LAI_s_new / (LAI_s_prev + 1e-12), 0.0)
        # scale only species s layers
        self.LAI_layers_SK[s, :, :, :] = np.clip(self.LAI_layers_SK[s, :, :, :] * factor[None, :, :], 0.0, self.params.lai_max)
        # refresh aggregates
        self.LAI_layers = np.sum(self.LAI_layers_SK, axis=0)
        self.LAI = np.sum(self.LAI_layers, axis=0)

    def _seed_based_spread(self, soil: np.ndarray | float | None = None) -> np.ndarray | None:
        """
        Legacy total-LAI seed-based spread (kept for backward compatibility).
        """
        return self._seed_based_spread_species(None)

    def _seed_based_spread_species(self, s_idx: int | None) -> np.ndarray | None:
        """
        Per-species seed-based colonization.
        - If s_idx is None: legacy total-LAI version.
        - Else: use species' LAI share to apportion reproduction energy and seedling establishment.
        """
        r0 = float(max(0.0, min(0.5, getattr(self, "spread_rate", 0.0))))
        if r0 <= 0.0:
            return None

        land = self.land
        H, W = self.shape
        E_map = np.nan_to_num(self.E_day)
        repro_frac = float(np.clip(getattr(self, "repro_fraction", 0.0), 0.0, 0.95))
        seed_energy = float(max(1e-12, getattr(self, "seed_energy", 1.0)))
        seed_scale = float(max(1e-12, getattr(self, "seed_scale", 1.0)))
        s_lai = float(max(0.0, getattr(self, "seedling_lai", 0.02)))

        # Neighbor set
        neigh_mode = str(getattr(self, "spread_neighbors", "vonneumann")).lower()
        if neigh_mode in ("moore", "8", "8n"):
            offsets = [(-1,-1), (-1,0), (-1,1),
                       ( 0,-1),          ( 0,1),
                       ( 1,-1), ( 1,0),  ( 1,1)]
        else:
            offsets = [(-1,0), (0,-1), (0,1), (1,0)]

        # If per-species
        if s_idx is not None and getattr(self, "LAI_layers_SK", None) is not None:
            S, K, _, _ = self.LAI_layers_SK.shape
            s = int(np.clip(s_idx, 0, S - 1))
            LAI_s = np.maximum(np.sum(self.LAI_layers_SK[s, :, :, :], axis=0), 0.0)
            LAI_tot = np.maximum(np.sum(self.LAI_layers_SK, axis=(0, 1)), 0.0)
            with np.errstate(invalid="ignore", divide="ignore"):
                share = np.where(LAI_tot > 0.0, LAI_s / (LAI_tot + 1e-12), 0.0)
            E_repro_s = repro_frac * E_map * share
            Seeds = np.maximum(E_repro_s / seed_energy, 0.0) * land.astype(float)
        else:
            # Legacy total-LAI
            E_repro = repro_frac * E_map
            Seeds = np.maximum(E_repro / seed_energy, 0.0) * land.astype(float)

        # Effective dispersal coefficient
        r_eff = r0 * (1.0 - np.exp(-Seeds / seed_scale))
        # Local seed-bank retention (M3): retain a fraction of produced seeds in current cell
        try:
            retain = float(os.getenv("QD_ECO_SEED_BANK_RETAIN", "0.2"))
        except Exception:
            retain = 0.2
        try:
            bank_max = float(os.getenv("QD_ECO_SEED_BANK_MAX", "1000.0"))
        except Exception:
            bank_max = 1000.0
        try:
            self.seed_bank = np.clip(self.seed_bank + retain * Seeds, 0.0, bank_max)
        except Exception:
            pass
        # Optional soil gating of seed dispersal rate
        gate = getattr(self, "_spread_gate", None)
        if gate is None:
            gate = land.astype(float)
        else:
            gate = np.where(land, np.clip(gate, 0.0, 1.0), 0.0)
        r_eff = r_eff * gate

        # Count valid land neighbors for each source
        num_valid = np.zeros((H, W), dtype=float)
        for dy, dx in offsets:
            num_valid += np.roll(land, shift=(-dy, -dx), axis=(0, 1)).astype(float)

        # Seeds exported per neighbor
        with np.errstate(invalid="ignore", divide="ignore"):
            seeds_share = np.where(num_valid > 0.0, r_eff * Seeds / (num_valid + 1e-12), 0.0)

        # Convert seeds to seedling LAI increments at destination
        add = np.zeros((H, W), dtype=float)
        for dy, dx in offsets:
            add += s_lai * np.roll(seeds_share, shift=(dy, dx), axis=(0, 1))
        # Cap seedling LAI daily addition to avoid explosive spread
        try:
            dmax_seed = float(os.getenv("QD_ECO_SEED_DLAI_MAX", "0.01"))
        except Exception:
            dmax_seed = 0.01
        add = np.minimum(add, dmax_seed)

        # Apply only on land
        seeded_mask = (add > 0.0) & land
        if np.any(seeded_mask):
            if s_idx is not None and getattr(self, "LAI_layers_SK", None) is not None:
                # Inject seedlings to species s at lowest layer
                self.LAI_layers_SK[s_idx, 0, seeded_mask] = np.clip(
                    self.LAI_layers_SK[s_idx, 0, seeded_mask] + add[seeded_mask],
                    0.0, self.params.lai_max
                )
                self.LAI_layers = np.sum(self.LAI_layers_SK, axis=0)
                self.LAI = np.sum(self.LAI_layers, axis=0)
            else:
                # Legacy: distribute to species by weights in k=0
                if getattr(self, "LAI_layers_SK", None) is not None:
                    try:
                        S = int(self.Ns)
                        w = np.asarray(self.species_weights, dtype=float)
                        w = w / (np.sum(w) + 1e-12)
                    except Exception:
                        S = 1
                        w = np.asarray([1.0], dtype=float)
                    for s in range(S):
                        self.LAI_layers_SK[s, 0, seeded_mask] = np.clip(
                            self.LAI_layers_SK[s, 0, seeded_mask] + w[s] * add[seeded_mask],
                            0.0, self.params.lai_max
                        )
                    self.LAI_layers = np.sum(self.LAI_layers_SK, axis=0)
                    self.LAI = np.sum(self.LAI_layers, axis=0)
                else:
                    self.LAI[seeded_mask] = np.clip(self.LAI[seeded_mask] + add[seeded_mask], 0.0, self.params.lai_max)
            # Reset age to 0 at newly seeded cells
            try:
                self.age_days[seeded_mask] = 0.0
            except Exception:
                pass

        return seeded_mask if np.any(seeded_mask) else None

    def canopy_reflectance_factor(self) -> np.ndarray:
        """
        Return f(LAI) in [0,1] used to scale leaf reflectance to canopy-scale reflectance.
        Using a saturating law: f(LAI) = 1 - exp(-k * LAI)
        Employ cached value if可用，按策略周期或LAI变化量触发重算。
        """
        if self._canopy_f_cached is None:
            self._recompute_canopy_cache()
        out = np.full(self.shape, np.nan, dtype=float)
        out[self.land] = self._canopy_f_cached[self.land]
        return out

    # M4: supply species leaf reflectance per band from adapter/genes (shape [Ns, NB])
    def set_species_reflectance_bands(self, R_leaf_species_nb: np.ndarray) -> None:
        """
        R_leaf_species_nb: array of shape [Ns, NB] with values in [0,1]
        """
        try:
            arr = np.asarray(R_leaf_species_nb, dtype=float)
            if arr.ndim != 2:
                return
            self._species_R_leaf = np.clip(arr, 0.0, 1.0)
        except Exception:
            self._species_R_leaf = None

    def effective_leaf_reflectance_bands(self, nb: int) -> np.ndarray:
        """
        Compute effective leaf reflectance per band by mixing species according to species_weights.
        Returns array [NB] in [0,1].
        """
        if self._species_R_leaf is None:
            # Fallback: flat neutral reflectance (0.5) if nothing provided
            return np.full((nb,), 0.5, dtype=float)
        Ns, NB = self._species_R_leaf.shape
        if NB != nb:
            # shape mismatch -> simple fallback
            return np.full((nb,), float(np.nanmean(self._species_R_leaf)), dtype=float)
        w = self.species_weights
        if w.size != Ns:
            # weights mismatch -> equal weights
            w = np.full((Ns,), 1.0 / max(1, Ns), dtype=float)
        # R_eff[b] = Σ_i w_i R_i[b]
        return np.clip(np.tensordot(w, self._species_R_leaf, axes=(0, 0)), 0.0, 1.0)

    def get_surface_albedo_bands(self, nb: int, soil_ref: float = 0.20) -> np.ndarray:
        """
        Build A_b^surface (NB×lat×lon) using canopy factor f(LAI), mixed leaf reflectance per band,
        and soil_ref as background:
            A_b(x,y) = R_eff[b] * f(LAI(x,y)) + (1 - f(LAI(x,y))) * soil_ref
        Land-only; ocean is NaN.
        """
        f_canopy = self.canopy_reflectance_factor()  # [lat,lon] in [0,1]
        R_eff = self.effective_leaf_reflectance_bands(nb)  # [NB]
        h, w = self.shape
        A = np.full((nb, h, w), np.nan, dtype=float)
        for b in range(nb):
            Ab = R_eff[b] * f_canopy + (1.0 - f_canopy) * soil_ref
            # land-only
            Z = np.full((h, w), np.nan, dtype=float)
            Z[self.land] = np.clip(Ab[self.land], 0.0, 1.0)
            A[b, :, :] = Z
        return A

    # --- M2 helpers: canopy cache policy ---
    def _should_recompute_canopy(self) -> bool:
        try:
            if self._canopy_f_cached is None:
                return True
            # Time-based
            if self._hours_accum >= self._next_recompute_hours:
                return True
            # LAI change threshold
            lai_now = self.total_LAI()
            delta = np.nanmean(np.abs(lai_now - self._lai_snapshot))
            base = np.nanmean(np.maximum(self._lai_snapshot, 1e-6))
            ratio = (delta / base) if base > 0 else delta
            return bool(ratio >= self._lai_recompute_delta)
        except Exception:
            return True

    def _recompute_canopy_cache(self) -> None:
        k = self.params.k_canopy
        LAI_tot = self.total_LAI()
        f = 1.0 - np.exp(-k * np.maximum(LAI_tot, 0.0))
        self._canopy_f_cached = f

    def lai_delta_ratio(self) -> float:
        lai_now = self.total_LAI()
        delta = np.nanmean(np.abs(lai_now - self._lai_snapshot))
        base = np.nanmean(np.maximum(self._lai_snapshot, 1e-6))
        return float((delta / base) if base > 0 else delta)

    # --- M1: genotype absorbance cache API ---
    def get_absorbance_for_genes(self, bands, genes) -> np.ndarray:
        """
        返回给定 Genes 在指定 bands 下的 A_b（缓存）。
        """
        try:
            peaks = getattr(genes, "absorption_peaks", []) or []
            pk_key = tuple((float(getattr(p, "center_nm", 0.0)),
                            float(getattr(p, "width_nm", 0.0)),
                            float(getattr(p, "height", 0.0))) for p in peaks)
            lam_key = tuple(float(x) for x in np.asarray(bands.lambda_centers, dtype=float).ravel().tolist())
            key = f"{getattr(genes, 'identity', 'gene')}|{pk_key}|{lam_key}"
        except Exception:
            key = f"{id(genes)}|{getattr(bands, 'nbands', 0)}"
        A = self._genotype_absorb_cache.get(key)
        if A is None:
            try:
                A = absorbance_from_genes(bands, genes)
                A = np.clip(np.asarray(A, dtype=float).ravel(), 0.0, 1.0)
            except Exception:
                A = np.full((getattr(bands, "nbands", 1),), 0.5, dtype=float)
            self._genotype_absorb_cache[key] = A
        return A

    def summary(self) -> dict:
        """Return simple diagnostics for logging."""
        land = self.land
        L = self.total_LAI()[land]
        if L.size == 0:
            return {"LAI_min": 0.0, "LAI_mean": 0.0, "LAI_max": 0.0}
        return {
            "LAI_min": float(np.min(L)),
            "LAI_mean": float(np.mean(L)),
            "LAI_max": float(np.max(L)),
        }

from __future__ import annotations

import os
import numpy as np
from dataclasses import dataclass


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

        # Cohort layers (K) scaffold: optional vertical layering; extend to species S×K
        try:
            self.K = max(1, int(os.getenv("QD_ECO_COHORT_K", "1")))
        except Exception:
            self.K = 1
        # Species count
        try:
            self.Ns = int(self.species_weights.shape[0])
        except Exception:
            self.Ns = 1
        # LAI_layers_SK shape [S, K, lat, lon]: initialize by species weight×equal split over K
        self.LAI_layers_SK = np.zeros((self.Ns, self.K, *self.shape), dtype=float)
        for s in range(self.Ns):
            frac_s = float(self.species_weights[s]) if self.Ns > 0 else 1.0
            for k in range(self.K):
                self.LAI_layers_SK[s, k, :, :] = frac_s * (self.LAI / float(self.K))
        # Backward-compat aggregated [K,lat,lon]
        self.LAI_layers = np.sum(self.LAI_layers_SK, axis=0)

        # M4: species (genes) mixture support (weights sum to 1, default 1 species)
        # Parsed by adapter to compute banded leaf reflectance; here only store weights.
        species_weights_env = os.getenv("QD_ECO_SPECIES_WEIGHTS", "").strip()
        if species_weights_env:
            try:
                w = [float(x) for x in species_weights_env.split(",") if x.strip() != ""]
            except Exception:
                w = [1.0]
        else:
            w = [1.0]
        s = sum(w) if len(w) > 0 else 1.0
        self.species_weights = np.asarray([max(0.0, wi) for wi in w], dtype=float)
        if s <= 0:
            self.species_weights = np.asarray([1.0], dtype=float)
        else:
            self.species_weights /= s
        # Species leaf reflectance cache per band (filled by caller via set_species_reflectance_bands)
        self._species_R_leaf = None  # shape [Ns, NB]

    def step_subdaily(self, isr_total: np.ndarray, dt_seconds: float) -> None:
        """
        Accumulate daily energy buffer from incoming shortwave proxy.
        Use simple proportionality: dE = isr_total * dt.
        """
        if isr_total is None:
            return
        # Ensure shapes match
        if isr_total.shape != self.shape:
            isr = np.full(self.shape, float(np.nan))
            # broadcast center if 1D? Keep simple: clip to mean
            isr[:] = float(np.nanmean(isr_total))
        else:
            isr = isr_total
        dE = np.nan_to_num(isr) * float(dt_seconds)
        self.E_day += dE

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

        # Growth term from daily energy
        growth = P.growth_per_j * self.E_day
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

        # Reset daily energy buffer
        self.E_day[:] = 0.0

    def canopy_reflectance_factor(self) -> np.ndarray:
        """
        Return f(LAI) in [0,1] used to scale leaf reflectance to canopy-scale reflectance.
        Using a saturating law: f(LAI) = 1 - exp(-k * LAI)
        """
        k = self.params.k_canopy
        LAI_tot = self.total_LAI()
        f = 1.0 - np.exp(-k * np.maximum(LAI_tot, 0.0))
        # Keep only land; ocean remains NaN to let caller blend selectively
        out = np.full(self.shape, np.nan, dtype=float)
        out[self.land] = f[self.land]
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

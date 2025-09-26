#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Minimal checks for P018 M1/M2 completion:
- M1:
  * spectral absorbance cache path (PopulationManager.get_absorbance_for_genes)
  * Plant.update_one_day band-integrated path and reflectance_bands output
- M2:
  * canopy cache policy (time and LAI-delta triggers)
  * banded surface albedo aggregation (PopulationManager.get_surface_albedo_bands)
Run:
  python3 -m scripts.test_ecology_m1_m2
"""
import os
import sys
import types
import numpy as np

# Ensure package import
sys.path.append(os.path.dirname(os.path.dirname(__file__)))

from pygcm.ecology.population import PopulationManager
from pygcm.ecology.plant import Plant, PlantState
from pygcm.ecology.genes import Genes

class DummyBands:
    def __init__(self, nb):
        self.nbands = nb
        self.lambda_centers = np.linspace(400.0, 700.0, nb)
        self.delta_lambda = np.full((nb,), (700.0-400.0)/nb)
        self.band_weights = np.full((nb,), 1.0/nb)

def make_dummy_genes():
    # Simple two-peak absorption
    g = Genes.from_env()  # start with defaults
    # build peaks if empty
    try:
        if (not getattr(g, "absorption_peaks", None)) or len(g.absorption_peaks) == 0:
            peak = types.SimpleNamespace(center_nm=450.0, width_nm=40.0, height=0.8)
            peak2 = types.SimpleNamespace(center_nm=680.0, width_nm=30.0, height=0.9)
            g.absorption_peaks = [peak, peak2]
    except Exception:
        pass
    # ensure allocations normalized
    s = g.alloc_root + g.alloc_stem + g.alloc_leaf
    if s <= 0:
        g.alloc_root, g.alloc_stem, g.alloc_leaf = 0.3, 0.2, 0.5
    else:
        g.alloc_root /= s; g.alloc_stem /= s; g.alloc_leaf /= s
    # ensure growth params
    if not hasattr(g, "leaf_area_per_energy"):
        g.leaf_area_per_energy = 0.1
    if not hasattr(g, "gdd_germinate"):
        g.gdd_germinate = 10.0
    if not hasattr(g, "lifespan_days"):
        g.lifespan_days = 200
    if not hasattr(g, "drought_tolerance"):
        g.drought_tolerance = 0.3
    return g

def check_close(name, a, b, tol=1e-6):
    if not (abs(a-b) <= tol):
        raise AssertionError(f"{name} not close: {a} vs {b}")

def main():
    # Set environment for deterministic cache policy
    os.environ.setdefault("QD_ECO_LIGHT_UPDATE_EVERY_HOURS", "6")
    os.environ.setdefault("QD_ECO_LIGHT_RECOMPUTE_LAI_DELTA", "0.05")
    os.environ.setdefault("QD_ECO_LAI_K_EXT", "0.4")

    # Construct land mask (small grid)
    H, W = 8, 12
    land = np.zeros((H, W), dtype=int)
    land[2:7, 3:10] = 1

    pm = PopulationManager(land_mask=land, diag=False)

    # Provide species reflectance bands (Ns, NB)
    NB = 8
    R_leaf = np.clip(np.linspace(0.1, 0.5, NB), 0, 1)  # monotonic increasing reflectance
    R_species = np.stack([R_leaf for _ in range(pm.Ns)], axis=0)
    pm.set_species_reflectance_bands(R_species)

    # M1: genotype absorbance cache
    bands = DummyBands(NB)
    genes = make_dummy_genes()
    A1 = pm.get_absorbance_for_genes(bands, genes)
    A2 = pm.get_absorbance_for_genes(bands, genes)  # cached second call
    if A1.shape[0] != NB:
        raise AssertionError("Absorbance length mismatch")
    diff_cache = float(np.max(np.abs(A1 - A2)))
    if diff_cache != 0.0:
        raise AssertionError("Absorbance cache mismatch between repeated calls")

    # M2: canopy cache recompute by time and delta
    # initial compute is lazy, ensure factor exists
    f0 = pm.canopy_reflectance_factor()
    if not np.all(np.isnan(f0) | ((f0 >= 0.0) & (f0 <= 1.0))):
        raise AssertionError("Initial canopy reflectance factor out of [0,1]")

    # advance subdaily without large LAI change, less than 6h
    isr = np.zeros((H, W)); isr[pm.land] = 200.0
    pm.step_subdaily(isr, dt_seconds=3600)  # +1h
    f1 = pm.canopy_reflectance_factor()
    # not guaranteed to change; check still within range
    if not np.all(np.isnan(f1) | ((f1 >= 0.0) & (f1 <= 1.0))):
        raise AssertionError("Canopy reflectance factor after 1h out of [0,1]")

    # force LAI change beyond threshold to trigger recompute
    pm.LAI[pm.land] += 0.1  # large LAI bump
    pm.step_subdaily(isr, dt_seconds=300)  # tick clock
    f2 = pm.canopy_reflectance_factor()
    if not np.all(np.isnan(f2) | ((f2 >= 0.0) & (f2 <= 1.0))):
        raise AssertionError("Canopy reflectance factor after LAI bump out of [0,1]")

    # M2: banded surface albedo aggregation
    Ab = pm.get_surface_albedo_bands(NB, soil_ref=0.2)
    if Ab.shape != (NB, H, W):
        raise AssertionError(f"A_b^surface shape mismatch {Ab.shape}")
    # check range on land
    land_vals = Ab[:, pm.land]
    if not np.all((land_vals >= -1e-9) & (land_vals <= 1.0 + 1e-9)):
        raise AssertionError("A_b^surface values out of [0,1]")

    # M1: Plant.update_one_day band-integrated path and reflectance_bands
    p = Plant(genes=genes)
    # ensure plant alive and growing
    p.state = PlantState.GROWING
    I_b = np.full((NB,), 100.0)  # W m^-2 nm^-1 band-mean
    dl = bands.delta_lambda
    rep = p.update_one_day(
        Ts_day=290.0,
        day_length_hours=12.0,
        soil_water_index=0.5,
        I_bands_weighted_scalar=None,
        I_bands=I_b,
        A_b_genotype=A1,
        delta_lambda=dl,
        light_availability=0.8,
    )
    # energy gain should be positive
    if rep.energy_gain <= 0:
        raise AssertionError("Plant.update_one_day banded energy should be > 0")
    # reflectance bands provided and in [0,1]
    if rep.reflectance_bands is None or rep.reflectance_bands.shape[0] != NB:
        raise AssertionError("Plant.update_one_day should return reflectance_bands of length NB")
    if not np.all((rep.reflectance_bands >= 0.0) & (rep.reflectance_bands <= 1.0)):
        raise AssertionError("Plant.reflectance_bands out of [0,1]")

    print("[OK] P018 M1/M2 minimal checks passed.")

if __name__ == "__main__":
    main()

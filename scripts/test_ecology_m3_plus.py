#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Smoke test for P018 M3+:
- IndividualPool subdaily energy accumulation
- Daily coupling of per-cell reproduction energy into PopulationManager.seed_bank
- Subsequent germination increases LAI

Pass criteria (heuristic):
- After IndividualPool.step_daily, seed_bank sum > 0 (coupling worked)
- After PopulationManager.step_daily with germination on, LAI mean is non-decreasing vs pre-germination
"""

import os
import numpy as np
import sys

# Ensure package import
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from pygcm.grid import SphericalGrid
from pygcm.ecology.adapter import EcologyAdapter
from pygcm.ecology.individuals import IndividualPool

def main():
    # Deterministic environment and moderate parameters
    os.environ["QD_ECO_DIAG"] = "0"
    os.environ["QD_ECO_LAI_INIT"] = "0.05"
    os.environ["QD_ECO_LAI_MAX"] = "5.0"
    os.environ["QD_ECO_LAI_GROWTH"] = "1.0e-8"        # keep growth small to avoid instant saturation
    os.environ["QD_ECO_SOIL_STRESS_THRESH"] = "0.1"
    os.environ["QD_ECO_SOIL_STRESS_GAIN"] = "1.0"

    # Individual pool sampling
    os.environ["QD_ECO_INDIV_SAMPLE_FRAC"] = "0.5"
    os.environ["QD_ECO_INDIV_PER_CELL"] = "40"
    os.environ["QD_ECO_INDIV_SUBSTEPS_PER_DAY"] = "2"
    os.environ["QD_ECO_INDIV_STRESS_PENALTY"] = "0.0"   # disable penalty to reduce confounders

    # M3+ coupling ON
    os.environ["QD_ECO_INDIV_SEED_COUPLE"] = "1"

    # Reproduction and seed-bank parameters (coupling target)
    os.environ["QD_ECO_REPRO_FRACTION"] = "0.3"
    os.environ["QD_ECO_SEED_ENERGY"] = "1.0"
    os.environ["QD_ECO_SEED_BANK_RETAIN"] = "0.5"
    os.environ["QD_ECO_SEED_BANK_MAX"] = "1e9"

    # Germination/decay for next day
    os.environ["QD_ECO_SEED_GERMINATE_FRAC"] = "0.2"
    os.environ["QD_ECO_SEED_BANK_DECAY"] = "0.02"
    os.environ["QD_ECO_SEEDLING_LAI"] = "0.02"

    # Grid and masks
    H, W = 24, 48
    grid = SphericalGrid(n_lat=H, n_lon=W)
    land = np.ones((H, W), dtype=int)     # all land for simplicity

    # Adapter and pop
    eco = EcologyAdapter(grid, land)
    pop = eco.pop
    assert pop is not None, "PopulationManager not initialized"

    # Record initial LAI and seed bank
    lai0_mean = float(np.nanmean(pop.LAI[land == 1]))
    seed0 = float(np.sum(pop.seed_bank[land == 1]))

    # Individual pool
    pool = IndividualPool(grid, land, eco, sample_frac=0.5, per_cell=40, substeps_per_day=2, diag=False)

    # One "day": perform subdaily steps via try_substep by passing dt >= period once
    isr_A = np.where(land == 1, 300.0, 0.0)
    isr_B = np.where(land == 1, 120.0, 0.0)
    day_len = 24 * 3600.0
    # Two substeps configured: call with dt=day_len/2 twice to ensure 2 substeps consumed
    pool.try_substep(isr_A, isr_B, eco, soil_W_land=np.where(land == 1, 0.7, 0.0), dt_seconds=day_len/2, day_length_seconds=day_len)
    pool.try_substep(isr_A, isr_B, eco, soil_W_land=np.where(land == 1, 0.7, 0.0), dt_seconds=day_len/2, day_length_seconds=day_len)

    # Daily boundary: couple per-cell reproduction energy to seed_bank
    pool.step_daily(eco, soil_W_land=np.where(land == 1, 0.7, 0.0), Ts_map=None, day_length_hours=24.0)
    seed1 = float(np.sum(pop.seed_bank[land == 1]))
    if not (seed1 > seed0):
        raise AssertionError(f"Seed bank did not increase after indiv->seed coupling: {seed1} vs {seed0}")

    # Now perform PopulationManager daily step to germinate part of the seed bank
    lai_pre_germ = float(np.nanmean(pop.LAI[land == 1]))
    pop.step_daily(np.where(land == 1, 0.8, 0.0))
    lai_post_germ = float(np.nanmean(pop.LAI[land == 1]))
    if not (lai_post_germ >= lai_pre_germ):
        raise AssertionError(f"LAI did not stay non-decreasing after germination: {lai_post_germ} vs {lai_pre_germ}")

    print("[OK] P018 M3+ indiv→seed_bank coupling smoke test passed.")
    print(f"LAI mean: init={lai0_mean:.4f}, pre_germ={lai_pre_germ:.4f}, post_germ={lai_post_germ:.4f}")
    print(f"Seed bank sum: before={seed0:.2f}, after indiv-couple={seed1:.2f}")

if __name__ == "__main__":
    main()

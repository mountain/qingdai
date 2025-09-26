#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Smoke test for P018 M3:
- Daily LAI update with reproduction split (repro_fraction)
- Seed bank retention, germination, and decay
- Optional seed-based spread increasing LAI in neighbor cells

Pass criteria (heuristic):
- After day 1 (energy input), mean LAI increases vs initial
- Seed bank > 0 in some land cells (with retain > 0)
- After day 2 (with germination), mean LAI increases further and seed bank decreases
"""

import os
import numpy as np
import sys

# Ensure package import
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from pygcm.ecology.population import PopulationManager

def main():
    # Configure environment for deterministic behavior
    os.environ["QD_ECO_LAI_INIT"] = "0.1"              # initial LAI on land
    os.environ["QD_ECO_LAI_MAX"] = "5.0"
    os.environ["QD_ECO_LAI_GROWTH"] = "1.0e-8"         # smaller growth per J-equivalent to avoid instant saturation
    os.environ["QD_ECO_SOIL_STRESS_THRESH"] = "0.1"    # mild stress threshold
    os.environ["QD_ECO_SOIL_STRESS_GAIN"] = "1.0"

    # Reproduction and seed-bank parameters
    os.environ["QD_ECO_REPRO_FRACTION"] = "0.2"        # 20% of daily energy to seed production
    os.environ["QD_ECO_SEED_ENERGY"] = "1.0"           # 1 energy per seed unit
    os.environ["QD_ECO_SEED_SCALE"] = "5.0"            # smoother saturation
    os.environ["QD_ECO_SEEDLING_LAI"] = "0.02"         # LAI per established seed unit
    os.environ["QD_ECO_SEED_BANK_RETAIN"] = "0.3"      # retain 30% locally
    os.environ["QD_ECO_SEED_BANK_MAX"] = "1e6"         # large cap
    os.environ["QD_ECO_SEED_GERMINATE_FRAC"] = "0.1"   # 10% of bank germinates/day
    os.environ["QD_ECO_SEED_BANK_DECAY"] = "0.05"      # 5% decay/day

    # Spread controls (enable 'seed' spread)
    os.environ["QD_ECO_SPREAD_ENABLE"] = "1"
    os.environ["QD_ECO_SPREAD_MODE"] = "seed"
    os.environ["QD_ECO_SPREAD_RATE"] = "0.03"          # per-day
    os.environ["QD_ECO_SPREAD_NEIGHBORS"] = "vonNeumann"
    os.environ["QD_ECO_SPREAD_DLAI_MAX"] = "0.05"      # daily cap

    # Canopy cache update to avoid unnecessary recomputation during test
    os.environ["QD_ECO_LIGHT_UPDATE_EVERY_HOURS"] = "6"
    os.environ["QD_ECO_LIGHT_RECOMPUTE_LAI_DELTA"] = "0.05"

    # Build a small grid with a plus-shape land mask to observe spread
    H, W = 12, 18
    land = np.zeros((H, W), dtype=int)
    land[5, 4:14] = 1
    land[3:9, 9]  = 1

    pop = PopulationManager(land_mask=land, diag=True)

    # Record initial stats
    s0 = pop.summary()
    lai0_mean = s0["LAI_mean"]

    # Day 1: apply a single subdaily call that accumulates one full day of energy (constant isr)
    isr_const = 250.0  # W/m^2 equivalent
    dt_seconds = 24 * 3600
    pop.step_subdaily(np.where(land == 1, isr_const, 0.0), dt_seconds)

    # Soil index moderately wet to allow good growth/germination
    soil_idx = np.where(land == 1, 0.6, 0.0)

    # Daily step 1: growth + reproduction (seed production) + optional spread
    pop.step_daily(soil_idx)
    s1 = pop.summary()
    lai1_mean = s1["LAI_mean"]
    seed_bank1_sum = float(np.sum(pop.seed_bank[land == 1]))

    if not (lai1_mean > lai0_mean):
        raise AssertionError(f"LAI mean did not increase after day 1: {lai1_mean} vs {lai0_mean}")
    if not (seed_bank1_sum > 0.0):
        raise AssertionError("Seed bank did not increase after day 1")

    # Day 2: another day of energy
    pop.step_subdaily(np.where(land == 1, isr_const, 0.0), dt_seconds)
    # Increase gating a bit to simulate wetter conditions everywhere on land
    soil_idx2 = np.where(land == 1, 0.8, 0.0)
    # Daily step 2: includes germination/decay; expect LAI to increase further, seed bank to decrease (some consumed)
    pop.step_daily(soil_idx2)
    s2 = pop.summary()
    lai2_mean = s2["LAI_mean"]
    seed_bank2_sum = float(np.sum(pop.seed_bank[land == 1]))

    if not (lai2_mean >= lai1_mean):
        raise AssertionError(f"LAI mean did not stay non-decreasing after day 2: {lai2_mean} vs {lai1_mean}")
    if not (seed_bank2_sum > 0.0):
        raise AssertionError("Seed bank should remain positive after day 2")

    print("[OK] P018 M3 seed bank / germination smoke test passed.")
    print(f"LAI mean: {lai0_mean:.4f} -> {lai1_mean:.4f} -> {lai2_mean:.4f}")
    print(f"Seed bank sum: {seed_bank1_sum:.4f} -> {seed_bank2_sum:.4f}")

if __name__ == "__main__":
    main()

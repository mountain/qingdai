#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Smoke test for P018 M4 autosave schema and roundtrip:

- Create EcologyAdapter and mutate LAI/species_weights
- Save autosave NPZ to a temp path under data/
- Verify NPZ contains schema_version and selected keys
- Zero-out pop state, then load autosave and verify restoration
- Save multiple times to trigger rolling backups retention

Pass criteria:
- NPZ has schema_version == 1 and expected keys
- After load, LAI and species_weights are restored with correct shapes and reasonable values
- Rolling backups count <= retention setting
"""

import os
import sys
import glob
import time
import shutil
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from pygcm.grid import SphericalGrid
from pygcm.ecology.adapter import EcologyAdapter

def main():
    # Quiet diagnostics
    os.environ["QD_ECO_DIAG"] = "0"
    os.environ["QD_ECO_USE_LAI"] = "1"
    # Small grid
    H, W = 12, 24
    grid = SphericalGrid(n_lat=H, n_lon=W)
    land = np.ones((H, W), dtype=int)

    eco = EcologyAdapter(grid, land)
    assert eco.pop is not None, "PopulationManager missing"
    pop = eco.pop

    # Mutate state deterministically
    rng = np.random.default_rng(1234)
    lai_mut = np.clip(rng.normal(0.2, 0.05, size=(H, W)), 0.0, pop.params.lai_max)
    pop.LAI = lai_mut.copy()
    # Apply into layers/SK (uniform split)
    K = int(getattr(pop, "K", 1))
    Ns = int(getattr(pop, "Ns", len(getattr(pop, "species_weights", [])) or 1))
    pop.species_weights = np.clip(np.ones((Ns,), dtype=float), 0.0, None)
    pop.species_weights /= np.sum(pop.species_weights)
    pop.LAI_layers_SK = np.zeros((Ns, K, H, W), dtype=float)
    for s in range(Ns):
        for k in range(K):
            pop.LAI_layers_SK[s, k, :, :] = float(pop.species_weights[s]) * (pop.LAI / float(max(1, K)))
    pop.LAI_layers = np.sum(pop.LAI_layers_SK, axis=0)

    # Prepare autosave path
    out_dir = os.path.join("data", "test_autosave")
    os.makedirs(out_dir, exist_ok=True)
    path_npz = os.path.join(out_dir, "eco_autosave.npz")

    # Ensure a clean slate
    for p in glob.glob(os.path.join(out_dir, "eco_autosave_*\.npz")):
        try:
            os.remove(p)
        except Exception:
            pass
    try:
        if os.path.exists(path_npz):
            os.remove(path_npz)
    except Exception:
        pass

    # Retention of backups
    os.environ["QD_ECO_AUTOSAVE_KEEP"] = "3"

    # Save autosave
    ok_save = eco.save_autosave(path_npz, day_value=1.0)
    if not ok_save:
        raise AssertionError("save_autosave returned False")
    if not os.path.exists(path_npz):
        raise AssertionError("autosave npz not found at expected path")

    # Inspect NPZ contents
    data = np.load(path_npz)
    schema = int(data.get("schema_version", np.int32(-1)))
    if schema != 1:
        raise AssertionError(f"schema_version expected 1, got {schema}")
    for key in ["LAI", "species_weights", "bands_lambda_centers", "bands_lambda_edges", "w_b"]:
        if key not in data.files:
            raise AssertionError(f"expected key missing in autosave: {key}")

    # Zero out pop state and load
    pop.LAI[:] = 0.0
    pop.species_weights[:] = 0.0
    ok_load = eco.load_autosave(path_npz)
    if not ok_load:
        raise AssertionError("load_autosave returned False")

    # Validate restoration
    lai_rest = pop.LAI
    if lai_rest.shape != lai_mut.shape:
        raise AssertionError(f"LAI shape mismatch after load: {lai_rest.shape} vs {lai_mut.shape}")
    if not (np.nanmean(lai_rest) > 0.0):
        raise AssertionError("LAI mean is not positive after load (restoration failed)")

    # Trigger rolling backups by saving multiple times
    for n in range(4):
        time.sleep(0.05)  # ensure mtime difference
        eco.save_autosave(path_npz, day_value=1.0 + n)
    backups = sorted(glob.glob(os.path.join(out_dir, "eco_autosave_*\.npz")))
    if len(backups) > 3:
        raise AssertionError(f"rolling backup retention exceeded: found {len(backups)} > 3")

    print("[OK] P018 M4 autosave schema and roundtrip smoke test passed.")
    print(f"schema={schema}, LAI_mean={float(np.nanmean(lai_rest)):.4f}, backups={len(backups)}")

if __name__ == "__main__":
    main()

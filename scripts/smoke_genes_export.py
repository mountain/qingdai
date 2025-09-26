#!/usr/bin/env python3
"""
Minimal smoke test for EcologyAdapter.export_genes.

- Instantiates a small spherical grid and a uniform land mask
- Builds EcologyAdapter (without running the full simulation loop)
- Calls export_genes("output", day_value=3.1)
- Prints a short validation summary:
    * schema_version
    * bands keys
    * whether each gene entry contains per-gene band arrays
    * sample peak keys
"""
import os
import json
import numpy as np

from pygcm.grid import SphericalGrid
from pygcm.ecology.adapter import EcologyAdapter


def find_latest_genes_json(out_dir: str) -> str | None:
    try:
        files = [f for f in os.listdir(out_dir) if f.startswith("genes_day_") and f.endswith(".json")]
        if not files:
            return None
        files.sort(reverse=True)
        return os.path.join(out_dir, files[0])
    except Exception:
        return None


def main():
    os.makedirs("output", exist_ok=True)

    # Small grid; uniform land for deterministic output
    grid = SphericalGrid(n_lat=121, n_lon=240)
    land_mask = np.ones_like(grid.lat_mesh, dtype=int)

    # Build adapter (disable LAI to minimize dependencies for this smoke)
    os.environ.setdefault("QD_ECO_USE_LAI", "0")
    adapter = EcologyAdapter(grid, land_mask)

    # Export and locate the newest file
    adapter.export_genes("output", day_value=3.1)
    path = find_latest_genes_json("output")
    if not path:
        print("ERROR: No output/genes_day_*.json found.")
        raise SystemExit(2)

    with open(path, "r", encoding="utf-8") as f:
        doc = json.load(f)

    print(f"Genes JSON: {path}")
    print("schema_version:", doc.get("schema_version"))
    print("bands keys:", list((doc.get("bands") or {}).keys()))

    genes = doc.get("genes") or []
    if not genes:
        print("WARNING: No genes entries found in JSON.")
        raise SystemExit(3)

    g0 = genes[0]
    per_gene_keys = ["lambda_centers_nm", "delta_lambda_nm", "lambda_edges_nm"]
    has_per_gene_bands = all(k in g0 for k in per_gene_keys)
    print("gene0 has per-gene bands:", has_per_gene_bands)
    print("gene0 keys (subset):", [k for k in g0.keys() if k in {"index","identity","provenance","peaks_model","peaks"}])

    peaks = g0.get("peaks") or []
    if peaks:
        print("gene0 peak keys sample:", list(peaks[0].keys()))
    else:
        print("WARNING: gene0 contains no peaks data.")

    # Success criteria: schema v3 and per-gene bands present
    ok = (doc.get("schema_version") == 3) and has_per_gene_bands
    print("OK:", ok)
    raise SystemExit(0 if ok else 4)


if __name__ == "__main__":
    main()

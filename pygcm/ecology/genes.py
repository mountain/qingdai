from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Tuple, Optional
import os
import math
import numpy as np


@dataclass
class Peak:
    """Gaussian absorption peak parameters in nm-domain."""
    center_nm: float
    width_nm: float
    height: float  # 0..1


@dataclass
class Genes:
    """
    Minimal gene definition for M4 scaffolding.

    Fields:
    - identity: string label
    - alloc_root/stem/leaf: energy allocation ratios (should sum ≈1, normalized internally)
    - leaf_area_per_energy: scalar to convert daily net energy to leaf area increment (m^2 per J-equiv)
    - absorption_peaks: list of gaussian peaks defining leaf absorbance spectrum
    - drought_tolerance: soil-water index threshold (0..1) below which stress accumulates faster
    - gdd_germinate: degree-day threshold for seed->growing transition (K·day)
    - lifespan_days: hard cap on age (days)
    """
    identity: str = "grass"
    alloc_root: float = 0.3
    alloc_stem: float = 0.2
    alloc_leaf: float = 0.5
    leaf_area_per_energy: float = 2.0e-3
    absorption_peaks: List[Peak] = field(default_factory=list)
    drought_tolerance: float = 0.3
    gdd_germinate: float = 80.0
    lifespan_days: int = 365
    provenance: Optional[str] = None

    @staticmethod
    def from_env(prefix: str = "QD_ECO_GENE_") -> "Genes":
        """
        Build a default 'grass-like' gene from environment variables.

        Peak string format (comma-separated peaks, each as center:width:height):
            QD_ECO_GENE_PEAKS="450:40:0.6, 680:30:0.8"
        """
        ident = os.getenv(prefix + "IDENTITY", "grass").strip()
        def f(name: str, default: float) -> float:
            try:
                return float(os.getenv(prefix + name, str(default)))
            except Exception:
                return default

        peaks_env = os.getenv(prefix + "PEAKS", "").strip()
        peaks: List[Peak] = []
        if peaks_env:
            parts = peaks_env.split(",")
            for p in parts:
                try:
                    c, w, h = p.strip().split(":")
                    peaks.append(Peak(float(c), float(w), float(h)))
                except Exception:
                    continue
        # Provide a reasonable default if none specified (two-band absorber)
        if not peaks:
            peaks = [Peak(450.0, 40.0, 0.6), Peak(680.0, 30.0, 0.8)]

        g = Genes(
            identity=ident,
            alloc_root=f("ALLOC_ROOT", 0.3),
            alloc_stem=f("ALLOC_STEM", 0.2),
            alloc_leaf=f("ALLOC_LEAF", 0.5),
            leaf_area_per_energy=f("LEAF_AREA_PER_EN", 2.0e-3),
            absorption_peaks=peaks,
            drought_tolerance=f("DROUGHT_TOL", 0.3),
            gdd_germinate=f("GDD_GERMINATE", 80.0),
            lifespan_days=int(f("LIFESPAN_DAYS", 365)),
        )
        # Normalize allocations
        s = g.alloc_root + g.alloc_stem + g.alloc_leaf
        if s <= 0:
            g.alloc_root, g.alloc_stem, g.alloc_leaf = 0.3, 0.2, 0.5
        else:
            g.alloc_root /= s
            g.alloc_stem /= s
            g.alloc_leaf /= s
        g.provenance = f"env:{prefix}"
        return g


def absorbance_from_genes(bands, genes: Genes) -> np.ndarray:
    """
    Compute band absorbance A_b in [0,1] for a given Genes instance and spectral bands.

    bands: SpectralBands (from ecology.spectral)
      - lambda_centers [NB]
      - delta_lambda [NB]
    """
    lam_b = np.asarray(bands.lambda_centers, dtype=float)  # [NB]
    A = np.zeros_like(lam_b, dtype=float)
    for pk in genes.absorption_peaks:
        if pk.width_nm <= 0 or pk.height <= 0:
            continue
        # Gaussian peak
        A += pk.height * np.exp(-((lam_b - pk.center_nm) ** 2) / (2 * (pk.width_nm ** 2)))
    # Clip to [0,1]
    return np.clip(A, 0.0, 1.0)


def reflectance_from_genes(bands, genes: Genes) -> np.ndarray:
    """
    Leaf reflectance spectrum (band-averaged) from genes:
      R_leaf[b] = 1 - A_leaf[b]
    """
    A = absorbance_from_genes(bands, genes)
    return np.clip(1.0 - A, 0.0, 1.0)

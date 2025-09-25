from __future__ import annotations

# Re-export public API for pygcm.ecology

from .adapter import EcologyAdapter
from .population import PopulationManager
from .spectral import (
    SpectralBands,
    make_bands,
    band_weights_from_mode,
    default_leaf_reflectance,
)
from .types import WeatherInstant, WeatherDaily
from .phyto import PhytoManager

__all__ = [
    "EcologyAdapter",
    "PopulationManager",
    "SpectralBands",
    "make_bands",
    "band_weights_from_mode",
    "default_leaf_reflectance",
    "WeatherInstant",
    "WeatherDaily",
    "PhytoManager",
]

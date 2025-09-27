"""
P020 Phase 0: world/ skeleton (façade) with DI-ready stubs.

This module introduces minimal classes to support the QD_USE_OO switch and to
prepare for Phase 1–5 without changing runtime behavior yet. In Phase 0, the
legacy engine in scripts/run_simulation.py remains the execution path.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import TYPE_CHECKING

# Typing-only import to keep mypy happy without requiring runtime import
if TYPE_CHECKING:  # pragma: no cover
    from pygcm.grid import SphericalGrid as _SphericalGrid
else:  # pragma: no cover
    _SphericalGrid = object  # sentinel for type hints

# Runtime import (optional; keeps Phase 0 lightweight)
try:
    from pygcm.grid import SphericalGrid as _RuntimeSphericalGrid
except Exception:  # pragma: no cover
    _RuntimeSphericalGrid = None


# ---------------------------
# Configuration & Parameters
# ---------------------------


@dataclass(frozen=True)
class SimConfig:
    """Minimal simulation configuration for Phase 0 (env-driven)."""

    n_lat: int = 121
    n_lon: int = 240
    dt_seconds: float = 300.0
    use_ocean: bool = True
    use_ecology: bool = True
    use_routing: bool = True

    @classmethod
    def from_env(cls) -> SimConfig:
        def _ibool(name: str, default: str = "1") -> bool:
            try:
                return int(os.getenv(name, default)) == 1
            except Exception:
                return default == "1"

        def _int(name: str, default: str) -> int:
            try:
                return int(os.getenv(name, default))
            except Exception:
                return int(default)

        def _float(name: str, default: str) -> float:
            try:
                return float(os.getenv(name, default))
            except Exception:
                return float(default)

        return cls(
            n_lat=_int("QD_N_LAT", "121"),
            n_lon=_int("QD_N_LON", "240"),
            dt_seconds=_float("QD_DT_SECONDS", "300"),
            use_ocean=_ibool("QD_USE_OCEAN", "1"),
            use_ecology=_ibool("QD_ECO_ENABLE", "1"),
            use_routing=_ibool("QD_HYDRO_ENABLE", "1"),
        )


@dataclass(frozen=True)
class PhysicsParams:
    """Placeholder for physics parameters (Phase 0)."""

    # Add fields in Phase 1–2 when promoting to Pydantic models.
    pass


@dataclass(frozen=True)
class SpectralBands:
    """Placeholder for spectral band definitions (Phase 0)."""

    # Add fields in Phase 2 when banded optics are factored out.
    pass


@dataclass(frozen=True)
class EcologyParams:
    """Placeholder for ecology parameters (Phase 0)."""

    # Add fields in later phases.
    pass


@dataclass(frozen=True)
class ParamsRegistry:
    """Aggregate of parameter models (Phase 0 placeholder)."""

    physics: PhysicsParams = PhysicsParams()
    bands: SpectralBands = SpectralBands()
    ecology: EcologyParams = EcologyParams()


# -------------
# World State
# -------------


@dataclass
class WorldState:
    """Minimal world state for Phase 0. Extended in later phases."""

    t_seconds: float = 0.0


# ----------------
# Qingdai World
# ----------------


class QingdaiWorld:
    """
    Phase 0 façade:
    - Supports DI (dependency injection) via constructor keyword args.
    - Provides create_default() to assemble from env + grid.
    - Implements nominal double-buffer layout without mutating legacy engine.
    """

    def __init__(
        self,
        config: SimConfig,
        params: ParamsRegistry,
        grid: _SphericalGrid | None,
        *,
        state: WorldState | None = None,
        atmos=None,
        ocean=None,
        surface=None,
        hydrology=None,
        routing=None,
        ecology=None,
        forcing=None,
    ) -> None:
        self.config = config
        self.params = params
        self.grid = grid
        # Double buffer placeholders (Phase 0: no external use)
        self.current_state = state or WorldState(t_seconds=0.0)
        self.next_state = WorldState(t_seconds=0.0)

        # DI slots (kept for future phases; unused in Phase 0)
        self.atmos = atmos
        self.ocean = ocean
        self.surface = surface
        self.hydrology = hydrology
        self.routing = routing
        self.ecology = ecology
        self.forcing = forcing

    @classmethod
    def create_default(cls) -> QingdaiWorld:
        cfg = SimConfig.from_env()
        grd = (
            _RuntimeSphericalGrid(cfg.n_lat, cfg.n_lon)
            if _RuntimeSphericalGrid is not None
            else None
        )
        pr = ParamsRegistry()
        return cls(cfg, pr, grd)

    # Phase 0: value-semantics step() without touching legacy runtime
    def step(self) -> WorldState:
        """Advance internal clock by one dt and swap buffers (placeholder)."""
        self.next_state.t_seconds = self.current_state.t_seconds + float(self.config.dt_seconds)
        # Buffer swap (no copies)
        self.current_state, self.next_state = self.next_state, self.current_state
        return self.current_state

    def run(self, n_steps: int | None = None, duration_days: float | None = None) -> None:
        """
        Phase 0: façade notice. The legacy engine in scripts/run_simulation.py
        remains the execution path in this phase. This method is a stub to allow
        early wiring and testing of imports/switches without behavior changes.
        """
        msg = (
            "[P020 Phase 0] QingdaiWorld skeleton is active (façade installed). "
            "Legacy engine in scripts/run_simulation.py continues to run; "
            "this stub will be replaced in later phases."
        )
        print(msg)


__all__ = [
    "SimConfig",
    "PhysicsParams",
    "SpectralBands",
    "EcologyParams",
    "ParamsRegistry",
    "WorldState",
    "QingdaiWorld",
]

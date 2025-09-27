"""
pytest configuration for Phase 0

Goals:
- keep tests fast and deterministic
- avoid plotting and heavy autosave I/O during imports or quick runs
- shrink default grid unless a test overrides explicitly
"""

import os
import sys

import pytest

# Ensure project root on sys.path for 'pygcm' and 'scripts' imports
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)


@pytest.fixture(autouse=True)
def _phase0_env(monkeypatch):
    # Small grid by default (tests can override via monkeypatch in the test)
    monkeypatch.setenv("QD_N_LAT", os.getenv("QD_N_LAT", "10"))
    monkeypatch.setenv("QD_N_LON", os.getenv("QD_N_LON", "20"))
    # Short-circuit plots
    monkeypatch.setenv("QD_PLOT_EVERY_DAYS", os.getenv("QD_PLOT_EVERY_DAYS", "1000000"))
    # Force non-interactive backend for matplotlib (avoid display requirements)
    monkeypatch.setenv("MPLBACKEND", os.getenv("MPLBACKEND", "Agg"))
    # Disable autosave and loading by default for tests
    monkeypatch.setenv("QD_AUTOSAVE_ENABLE", os.getenv("QD_AUTOSAVE_ENABLE", "0"))
    monkeypatch.setenv("QD_AUTOSAVE_LOAD", os.getenv("QD_AUTOSAVE_LOAD", "0"))
    # Disable optional heavy subsystems unless a test enables them explicitly
    monkeypatch.setenv("QD_USE_OCEAN", os.getenv("QD_USE_OCEAN", "0"))
    monkeypatch.setenv("QD_ECO_ENABLE", os.getenv("QD_ECO_ENABLE", "0"))
    monkeypatch.setenv("QD_PHYTO_ENABLE", os.getenv("QD_PHYTO_ENABLE", "0"))
    monkeypatch.setenv("QD_HYDRO_ENABLE", os.getenv("QD_HYDRO_ENABLE", "0"))
    # Ensure no routing generation in tests
    monkeypatch.setenv("QD_HYDRO_NETCDF", os.getenv("QD_HYDRO_NETCDF", ""))
    # Keep JAX path deterministic; tests validate boolean API only
    monkeypatch.setenv("QD_USE_JAX", os.getenv("QD_USE_JAX", "0"))
    yield

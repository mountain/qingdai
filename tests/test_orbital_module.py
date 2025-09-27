"""
Pytest migration of scripts/test_orbital_module.py

Covers invariants and basic properties without plotting:
- Periods/omegas consistency
- Barycenter geometry (r_A + r_B == A_BINARY; M_A*r_A == M_B*r_B)
- Vectorized API shapes for positions and flux
- Flux series positivity and variability over one planetary year
"""

import numpy as np


def test_periods_and_omegas_consistency():
    from pygcm.orbital import OrbitalSystem

    orb = OrbitalSystem()

    # Periods must be positive/finite
    assert orb.T_binary > 0.0 and np.isfinite(orb.T_binary)
    assert orb.T_planet > 0.0 and np.isfinite(orb.T_planet)

    # Angular velocities consistent: omega = 2π / T
    assert np.isclose(orb.omega_binary, 2.0 * np.pi / orb.T_binary, rtol=1e-12, atol=0.0)
    assert np.isclose(orb.omega_planet, 2.0 * np.pi / orb.T_planet, rtol=1e-12, atol=0.0)


def test_barycenter_geometry():
    from pygcm import constants as const
    from pygcm.orbital import OrbitalSystem

    orb = OrbitalSystem()

    # Geometry about barycenter
    # Sum of orbital radii equals binary semi-major axis
    assert np.isclose(orb.r_A + orb.r_B, const.A_BINARY, rtol=0.0, atol=0.0)
    # Torque balance around barycenter: M_A * r_A == M_B * r_B
    assert np.isclose(const.M_A * orb.r_A, const.M_B * orb.r_B, rtol=1e-12, atol=0.0)


def test_vectorized_positions_and_flux_shapes():
    from pygcm.orbital import OrbitalSystem

    orb = OrbitalSystem()

    # Vectorized time
    t = np.linspace(0.0, orb.T_planet, 128, dtype=float)
    x_A, y_A, x_B, y_B = orb.calculate_stellar_positions(t)

    # All arrays must align by shape
    assert isinstance(x_A, np.ndarray) and x_A.shape == t.shape
    assert isinstance(y_A, np.ndarray) and y_A.shape == t.shape
    assert isinstance(x_B, np.ndarray) and x_B.shape == t.shape
    assert isinstance(y_B, np.ndarray) and y_B.shape == t.shape

    flux = orb.calculate_total_flux(t)
    assert isinstance(flux, np.ndarray) and flux.shape == t.shape


def test_flux_series_properties_over_one_year():
    from pygcm.orbital import OrbitalSystem

    orb = OrbitalSystem()

    # Sample one full planetary year
    n = 1000
    t = np.linspace(0.0, orb.T_planet, n, dtype=float)
    flux = orb.calculate_total_flux(t)

    # Positivity and finiteness
    assert np.all(np.isfinite(flux))
    assert np.all(flux > 0.0)

    # Should show some variability over the year (non-constant series)
    fmin, fmean, fmax = float(np.min(flux)), float(np.mean(flux)), float(np.max(flux))
    assert fmin < fmax
    assert fmin < fmean < fmax

    # Enforce minimal relative variability (very loose threshold)
    # Avoid over-constraining: only check > 0.1% relative spread
    spread = (fmax - fmin) / fmean if fmean > 0 else 0.0
    assert spread > 1e-3

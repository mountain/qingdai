"""
Phase 0: initial pytest test suite

Covers:
- world façade basics (creation, env-driven config, step() clock advance)
- orbital module basic invariants and a t=0 flux check against constants
- jax_compat availability and boolean API
- scripts.run_simulation OO strict short-circuit (QD_USE_OO=1, QD_USE_OO_STRICT=1)
"""

import importlib
import time

import numpy as np


def test_world_create_and_step(monkeypatch):
    # Make grid small for quick instantiation
    monkeypatch.setenv("QD_N_LAT", "10")
    monkeypatch.setenv("QD_N_LON", "20")
    monkeypatch.setenv("QD_DT_SECONDS", "123")

    from pygcm.world import QingdaiWorld, SimConfig

    cfg = SimConfig.from_env()
    assert cfg.n_lat == 10
    assert cfg.n_lon == 20
    assert cfg.dt_seconds == 123.0

    world = QingdaiWorld.create_default()
    # Grid may be None if import error, but in this repo the grid implementation exists
    assert world is not None
    assert hasattr(world, "config")
    assert hasattr(world, "grid")
    # step() should advance time by dt
    t0 = world.current_state.t_seconds
    st = world.step()
    assert st.t_seconds == t0 + cfg.dt_seconds
    # Double buffer swap implied: next step advances again
    st2 = world.step()
    assert st2.t_seconds == t0 + 2 * cfg.dt_seconds


def test_orbital_basic_and_flux_t0():
    # Import orbital system and constants
    from pygcm import constants as const
    from pygcm.orbital import OrbitalSystem

    orb = OrbitalSystem()
    # Periods must be positive and finite
    assert orb.T_binary > 0.0
    assert orb.T_planet > 0.0
    assert np.isfinite(orb.T_binary)
    assert np.isfinite(orb.T_planet)

    # At t=0, positions simplify to: x_A=r_A, x_B=-r_B, x_p=a_p (all y=0)
    t0 = 0.0
    x_A, y_A, x_B, y_B = orb.calculate_stellar_positions(t0)
    assert np.isclose(y_A, 0.0)
    assert np.isclose(y_B, 0.0)
    assert np.isclose(x_A, orb.r_A)
    assert np.isclose(x_B, -orb.r_B)

    # Distances
    x_p = const.A_PLANET * np.cos(orb.omega_planet * t0)
    y_p = const.A_PLANET * np.sin(orb.omega_planet * t0)
    assert np.isclose(y_p, 0.0)
    assert np.isclose(x_p, const.A_PLANET)

    d_A = abs(x_p - x_A)  # = A_PLANET - r_A
    d_B = abs(x_p - x_B)  # = A_PLANET + r_B

    # Expected flux
    S_A_exp = const.L_A / (4.0 * np.pi * d_A**2)
    S_B_exp = const.L_B / (4.0 * np.pi * d_B**2)
    S_tot_exp = S_A_exp + S_B_exp

    S_tot = orb.calculate_total_flux(t0)
    # Within a tight relative tolerance
    assert np.isclose(S_tot, S_tot_exp, rtol=1e-12, atol=0.0)


def test_jax_compat_is_enabled_boolean():
    # is_enabled should be importable and return a boolean, regardless of jax presence
    from pygcm.jax_compat import is_enabled

    val = is_enabled()
    assert isinstance(val, bool | np.bool_)


def test_run_simulation_oo_strict_short_circuit(capsys, monkeypatch):
    """
    Ensure main() returns quickly when QD_USE_OO=1 and QD_USE_OO_STRICT=1,
    exercising the Phase 0 switch without executing the legacy engine.
    """
    # Set switches
    monkeypatch.setenv("QD_USE_OO", "1")
    monkeypatch.setenv("QD_USE_OO_STRICT", "1")
    # Keep plots disabled just in case
    monkeypatch.setenv("QD_PLOT_EVERY_DAYS", "1000000")

    # Import the module and run main
    mod = importlib.import_module("scripts.run_simulation")

    t0 = time.time()
    mod.main()
    elapsed = time.time() - t0

    out = capsys.readouterr().out
    # Should print the façade notice
    assert "QingdaiWorld" in out or "façade" in out or "fa\xc3\xa7ade" in out
    # Should return quickly (well below 1 second on CI)
    assert elapsed < 1.0

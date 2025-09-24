"""
jax_compat.py — Optional JAX acceleration layer for Qingdai PyGCM (P016 M1–M3)

Provides:
- JAX enable switch via env QD_USE_JAX (0/1)
- Compatibility helpers:
    * jax_map_coordinates: JAX-native map_coordinates fallback to SciPy if disabled
    * laplacian_sphere_jax: jitted spherical Laplacian
    * hyperdiffuse_jax: jitted ∇⁴ with optional substeps
    * advect_semilag_jax: jitted semi-Lagrangian advection (bilinear)
- to_numpy: safe conversion from device arrays to numpy
- is_enabled: query flag
"""
from __future__ import annotations

import os
import numpy as _np
from typing import Tuple

# Global flag — do not raise if JAX unavailable; just fall back to NumPy/SciPy
_JAX_ENABLED = False
_JAX = None
_JNP = None
_JAX_SCIPY_NDIMAGE = None

try:
    _JAX_ENABLED = int(os.getenv("QD_USE_JAX", "0")) == 1
except Exception:
    _JAX_ENABLED = False

if _JAX_ENABLED:
    try:
        import jax as _JAX
        import jax.numpy as _JNP
        import jax.scipy.ndimage as _JAX_SCIPY_NDIMAGE
        # Optional: select platform (cpu|gpu|tpu)
        plat = os.getenv("QD_JAX_PLATFORM")
        if plat:
            # Note: Must be set before JAX backend initialization in real systems;
            # here we best-effort via environment variable.
            os.environ.setdefault("JAX_PLATFORM_NAME", plat)
        _JAX_ENABLED = True
    except Exception:
        # If import fails, silently disable JAX
        _JAX = None
        _JNP = None
        _JAX_SCIPY_NDIMAGE = None
        _JAX_ENABLED = False


def is_enabled() -> bool:
    return _JAX_ENABLED


def to_numpy(x):
    """Convert JAX array (if enabled) to NumPy array; else return input."""
    if _JAX_ENABLED:
        try:
            return _np.asarray(x)
        except Exception:
            return x
    return x


def jax_map_coordinates(arr, coords, order: int = 1, mode: str = "wrap", prefilter: bool = False):
    """
    JAX-compatible map_coordinates.
    - Uses jax.scipy.ndimage.map_coordinates when JAX is enabled
    - Falls back to scipy.ndimage.map_coordinates otherwise
    Note: prefilter is ignored in JAX path (no-op).
    """
    if _JAX_ENABLED and (_JAX_SCIPY_NDIMAGE is not None):
        # JAX wants coords as a sequence of arrays
        return _JAX_SCIPY_NDIMAGE.map_coordinates(arr, coords, order=order, mode=mode)
    else:
        from scipy.ndimage import map_coordinates as _sc_map
        return _sc_map(arr, coords, order=order, mode=mode, prefilter=prefilter)


# ---------------- JAX-jitted kernels (with NumPy fallbacks) ---------------- #

def laplacian_sphere(F, dlat: float, dlon: float, coslat, a: float):
    """
    Spherical Laplacian of scalar F using divergence form with cosφ metric.
    If JAX enabled: jitted; else NumPy implementation.
    """
    if _JAX_ENABLED:
        @ _JAX.jit
        def _lap(F_, coslat_):
            F_ = _JNP.nan_to_num(F_)
            dF_dphi = _JNP.gradient(F_, dlat, axis=0)
            term_phi = (1.0 / coslat_) * _JNP.gradient(coslat_ * dF_dphi, dlat, axis=0)
            d2F_dlam2 = (_JNP.roll(F_, -1, axis=1) - 2.0 * F_ + _JNP.roll(F_, 1, axis=1)) / (dlon ** 2)
            term_lam = d2F_dlam2 / (coslat_ ** 2)
            return (term_phi + term_lam) / (a ** 2)
        return _lap(F, coslat)
    else:
        F = _np.nan_to_num(F)
        dF_dphi = _np.gradient(F, dlat, axis=0)
        term_phi = (1.0 / coslat) * _np.gradient(coslat * dF_dphi, dlat, axis=0)
        d2F_dlam2 = (_np.roll(F, -1, axis=1) - 2.0 * F + _np.roll(F, 1, axis=1)) / (dlon ** 2)
        term_lam = d2F_dlam2 / (coslat ** 2)
        return (term_phi + term_lam) / (a ** 2)


def hyperdiffuse(F, k4, dt: float, n_substeps: int, dlat: float, dlon: float, coslat, a: float):
    """
    Apply explicit 4th-order hyperdiffusion dF/dt = -k4 ∇⁴ F
    - k4 can be scalar or array broadcastable to F
    - If JAX enabled: jitted with lax.fori_loop; else NumPy fallback
    """
    if dt <= 0.0:
        return F
    if _JAX_ENABLED:
        k4_is_scalar = False
        try:
            k4_is_scalar = _np.isscalar(k4)
        except Exception:
            k4_is_scalar = False

        @ _JAX.jit
        def _step_once(F_, k4_):
            L = laplacian_sphere(F_, dlat, dlon, coslat, a)
            L2 = laplacian_sphere(L, dlat, dlon, coslat, a)
            return F_ - k4_ * L2 * (dt / _JNP.maximum(1, n_substeps))

        @ _JAX.jit
        def _loop(F_):
            sub_dt = dt / _JNP.maximum(1, n_substeps)
            # Allow scalar or array k4
            k4_ = _JNP.array(k4) if not k4_is_scalar else _JNP.array(float(k4))
            def body(i, val):
                L = laplacian_sphere(val, dlat, dlon, coslat, a)
                L2 = laplacian_sphere(L, dlat, dlon, coslat, a)
                return val - k4_ * L2 * sub_dt
            return _JAX.lax.fori_loop(0, _JNP.maximum(1, n_substeps), body, _JNP.nan_to_num(F_))
        return _loop(F)
    else:
        # NumPy fallback
        try:
            if _np.isscalar(k4):
                k4_arr = float(k4)
                if k4_arr <= 0.0:
                    return F
            else:
                k4_arr = _np.nan_to_num(k4, copy=False)
                if _np.all(k4_arr <= 0.0):
                    return F
        except Exception:
            return F
        n = max(1, int(n_substeps))
        sub_dt = dt / n
        out = _np.nan_to_num(F, copy=True)
        for _ in range(n):
            L = laplacian_sphere(out, dlat, dlon, coslat, a)
            L2 = laplacian_sphere(L, dlat, dlon, coslat, a)
            out = out - k4_arr * L2 * sub_dt
        return _np.nan_to_num(out, copy=False)


def advect_semilag(field, u, v, dt: float, a: float, dlat: float, dlon: float, coslat):
    """
    Semi-Lagrangian advection with bilinear interpolation.
    coords in index space: (row, col) with longitude wrap.
    JAX path uses jax.scipy.ndimage.map_coordinates; fallback uses SciPy.
    """
    # Convert velocities to index displacements
    dlam = u * dt / (a * _np.maximum(1e-6, coslat))
    dphi = v * dt / a
    dx = dlam / dlon
    dy = dphi / dlat

    # Grid index meshes
    lats = _np.arange(field.shape[0])
    lons = _np.arange(field.shape[1])
    JJ, II = _np.meshgrid(lats, lons, indexing="ij")
    dep_J = JJ - dy
    dep_I = II - dx

    if _JAX_ENABLED:
        f = _JNP.asarray(field)
        dep_J_j = _JNP.asarray(dep_J)
        dep_I_j = _JNP.asarray(dep_I)
        return _JAX_SCIPY_NDIMAGE.map_coordinates(f, [dep_J_j, dep_I_j], order=1, mode="wrap")
    else:
        from scipy.ndimage import map_coordinates as _sc_map
        return _sc_map(field, [dep_J, dep_I], order=1, mode="wrap", prefilter=False)

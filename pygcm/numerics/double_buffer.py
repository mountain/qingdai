"""
DoubleBufferingArray (DBA): double-buffered array for atomic step updates.

Goals (P021 M1):
- Provide read/write separation with O(1) swap.
- Behave "array-like" enough for NumPy API via magic methods (non-JIT path):
  * __getitem__ reads from .read
  * __setitem__ writes to .write
  * __array__ exposes .read for NumPy interop
  * __array_ufunc__ routes np.ufunc calls; out=DBA writes .write
- Keep JAX kernels pure: jitted functions must consume .read and write back to .write outside.

Notes:
- Avoid replacing underlying buffers; always use slice writes to .write.
- In hotspot paths, prefer explicit `arr = dba.read` to avoid implicit copies via __array__.

Tests to cover (see tests/test_double_buffering.py):
- read/write isolation and swap() semantics
- __getitem__/__setitem__ behave as specified
- __array__: np.asarray(dba), np.sin(dba) read from .read
- __array_ufunc__: np.add(dba, 1) returns ndarray; np.add(dba, 1, out=dba) writes .write
- self-aliasing protection: dba[...] = dba raises
- JAX smoke: jnp.add(dba.read, 1) executes when JAX enabled
"""

from __future__ import annotations

from typing import Any

import numpy as _np

try:
    # Optional JAX interop (read path only)
    _HAS_JAX = True
except Exception:
    _HAS_JAX = False

try:
    # Local helper: robust conversion to NumPy for __array__
    from pygcm.jax_compat import to_numpy as _to_numpy
except Exception:

    def _to_numpy(x):
        return _np.array(x, copy=False)


class DoubleBufferingArray:
    """
    Double-buffered 2D/ND array with explicit read/write and O(1) swap.

    Construction:
      dba = DoubleBufferingArray((n_lat, n_lon), dtype=float, initial_value=0.0)

    Use:
      # read current state
      r = dba.read
      # write next state
      dba.write[:] = r + 1.0
      # atomically make writes visible to readers
      dba.swap()

    Magic methods (non-JIT convenience):
      - __getitem__(key) -> self.read[key]
      - __setitem__(key, value) -> self.write[key] = value
      - __array__(dtype) -> NumPy array view/copy of .read
      - __array_ufunc__: route np.ufunc; out=DBA writes to .write
    """

    __slots__ = ("_a", "_b", "_read_idx", "_write_synced", "__weakref__")
    __array_priority__ = 1000  # prefer our __array_ufunc__ over ndarray

    def __init__(self, shape: tuple[int, ...], dtype: Any = _np.float64, initial_value: Any = 0.0):
        # Allocate two independent buffers
        self._a = _np.full(shape, initial_value, dtype=dtype)
        self._b = _np.full(shape, initial_value, dtype=dtype)
        self._read_idx = 0  # 0 => _a is read, _b is write; 1 => reversed
        self._write_synced = False  # lazily mirror read->write on first write

    # ---- core properties ----
    @property
    def read(self):
        """Array representing the current readable buffer."""
        return self._a if self._read_idx == 0 else self._b

    @property
    def write(self):
        """Array representing the current writable buffer (next state)."""
        return self._b if self._read_idx == 0 else self._a

    def swap(self) -> None:
        """O(1) pointer flip between read and write buffers."""
        # Flip LSB
        self._read_idx ^= 1
        # After flipping, mark write as unsynced (copy-on-write semantics)
        self._write_synced = False

    # ---- convenience ----
    @property
    def shape(self) -> tuple[int, ...]:
        return self.read.shape

    @property
    def dtype(self) -> _np.dtype:
        return self.read.dtype

    def zero_write(self) -> None:
        """Zero the write buffer."""
        self.write[...] = 0
        self._write_synced = True

    # ---- magic: indexing ----
    def __getitem__(self, key):
        return self.read[key]

    def __setitem__(self, key, value):
        # prevent obvious self-aliasing (dba[...] = dba)
        if value is self:
            raise ValueError(
                "DoubleBufferingArray: self-aliasing write is not allowed (dba[...] = dba)."
            )
        # lazy sync write buffer with current read on first write after swap
        if not self._write_synced:
            self.write[...] = self.read
            self._write_synced = True
        # allow ndarray / scalar / sequences
        self.write[key] = value

    # ---- magic: numpy array coercion ----
    def __array__(self, dtype=None):
        # Return NumPy array of .read (copy if necessary to ensure dtype)
        arr = self.read
        try:
            return _to_numpy(arr if dtype is None else _np.asarray(arr, dtype=dtype))
        except Exception:
            return _np.array(arr, copy=dtype is None, dtype=dtype)

    # ---- magic: ufunc routing ----
    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """
        Minimal rule set:
        - If out is not provided or not DBA: replace DBA inputs with .read, call ufunc, return ndarray(s).
        - If out contains DBA(s): replace those out operands with their .write arrays, call ufunc, and return out.
        """
        if method != "__call__":
            return NotImplemented

        # Extract/replace inputs
        proc_inputs = []
        for x in inputs:
            if isinstance(x, DoubleBufferingArray):
                proc_inputs.append(x.read)
            else:
                proc_inputs.append(x)

        # Handle out=...
        out = kwargs.get("out", None)
        if out is None:
            # Pure functional: return ndarray result(s)
            return ufunc(*proc_inputs, **kwargs)

        # out can be a single array or a tuple of arrays
        if not isinstance(out, tuple):
            out = (out,)
        # Ensure any DBA out targets are synced (copy-on-write) before writing
        for y in out:
            if isinstance(y, DoubleBufferingArray) and not y._write_synced:
                y.write[...] = y.read
                y._write_synced = True

        out_arrays = []
        for y in out:
            if isinstance(y, DoubleBufferingArray):
                out_arrays.append(y.write)
            else:
                out_arrays.append(y)

        kwargs["out"] = tuple(out_arrays)
        res = ufunc(*proc_inputs, **kwargs)
        return res

    # ---- repr/help ----
    def __repr__(self) -> str:
        # Show small preview of read/write pointers
        return f"DoubleBufferingArray(shape={self.shape}, dtype={self.dtype}, read=buf{self._read_idx}, write=buf{1 ^ self._read_idx})"

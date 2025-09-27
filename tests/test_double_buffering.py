import numpy as np
import pytest
from pygcm.numerics.double_buffer import DoubleBufferingArray as DBA


def test_read_write_isolation_and_swap():
    dba = DBA((2, 3), dtype=float, initial_value=0.0)
    # Initially read is zeros
    np.testing.assert_allclose(dba.read, 0.0)
    # __setitem__ writes to write buffer, not read
    dba[...] = 1.0
    np.testing.assert_allclose(dba.read, 0.0)
    # After swap, read reflects previous writes
    dba.swap()
    np.testing.assert_allclose(dba.read, 1.0)
    # Write different pattern
    dba[0, :] = 5.0
    np.testing.assert_allclose(dba.read[0, :], 1.0)  # unchanged until swap
    dba.swap()
    np.testing.assert_allclose(dba.read[0, :], 5.0)
    np.testing.assert_allclose(dba.read[1, :], 1.0)


def test_getitem_and_setitem_semantics():
    dba = DBA((2, 2), dtype=float, initial_value=2.0)
    # read-only via __getitem__
    r00 = dba[0, 0]
    assert r00 == 2.0
    # __setitem__ writes to write
    dba[1, 1] = 9.0
    assert dba.read[1, 1] == 2.0
    dba.swap()
    assert dba.read[1, 1] == 9.0


def test_array_and_numpy_interop():
    dba = DBA((2, 2), dtype=float, initial_value=0.5)
    arr = np.asarray(dba)
    # arr corresponds to .read
    np.testing.assert_allclose(arr, dba.read)
    # ufunc without out returns ndarray; does not mutate write/read
    y = np.sin(dba)
    assert isinstance(y, np.ndarray)
    np.testing.assert_allclose(dba.read, 0.5)  # unchanged
    # out= writes to write
    np.add(dba, 1.0, out=dba)
    # Not visible yet
    np.testing.assert_allclose(dba.read, 0.5)
    dba.swap()
    np.testing.assert_allclose(dba.read, 1.5)


def test_array_ufunc_out_multiple():
    # Test that out as tuple works for multi-output ufuncs (where applicable)
    # Use np.divmod which returns tuple; emulate by writing only to first component via out
    # Note: np.divmod supports out with two arrays; here we check routing logic robustness.
    a = DBA((2, 2), dtype=int, initial_value=9)
    q = DBA((2, 2), dtype=int, initial_value=0)
    r = DBA((2, 2), dtype=int, initial_value=0)
    # Write dividend to write of 'a' and swap to read
    a[...] = 9
    a.swap()
    # out expects numpy arrays; our __array_ufunc__ will pass .write
    np.divmod(a, 4, out=(q, r))
    # Not visible until swap
    np.testing.assert_allclose(q.read, 0)
    np.testing.assert_allclose(r.read, 0)
    q.swap()
    r.swap()
    np.testing.assert_allclose(q.read, 2)
    np.testing.assert_allclose(r.read, 1)


def test_self_aliasing_write_raises():
    dba = DBA((2, 2), dtype=float, initial_value=0.0)
    with pytest.raises(ValueError):
        dba[...] = dba  # self-alias protection


@pytest.mark.parametrize("dtype", [np.float32, np.float64, int])
def test_zero_write_and_repr(dtype):
    dba = DBA((1, 3), dtype=dtype, initial_value=7)
    # Move 7's to read by swapping once after default write -> read pattern
    dba.swap()
    # zero write and then swap
    dba.zero_write()
    dba.swap()
    np.testing.assert_allclose(dba.read, 0)
    # repr does not crash
    assert "DoubleBufferingArray" in repr(dba)


def test_jax_smoke_if_available():
    try:
        import jax.numpy as jnp
    except Exception:
        pytest.skip("JAX not available")
    dba = DBA((2, 2), dtype=float, initial_value=3.0)
    dba.swap()  # make 3.0 visible in read
    z = jnp.add(dba.read, 1.0)
    # Convert to numpy for assertion
    z_np = np.asarray(z)
    np.testing.assert_allclose(z_np, 4.0)

"""The array namespace seam: host stays numpy, device arrays stay on the device, mixing raises.

``array_namespace`` reads the namespace off a kernel's inputs. Host inputs must give
the ``numpy`` module itself (so the numpy path is the old code, bit for bit), CuPy
inputs give CuPy, and a mix raises rather than moving data. The first three tests
run everywhere on a fake ``cupy`` module; the device-parity test of each geometry
kernel needs CuPy and a CUDA device and skips cleanly otherwise, as ``test_gpu.py``
does.
"""

from __future__ import annotations

import sys
from types import ModuleType
from typing import Callable

import numpy as np
import pytest

from tanglepack.numerics.geometry import chord_deviation_batch, curvature_area_batch
from tanglepack.numerics.namespace import array_namespace


def _layer() -> tuple[np.ndarray, ...]:
    """A smooth 400-point refine layer ``(p0, p1, left, right)``, its first left neighbour NaN."""
    t = np.linspace(0.0, 3.0, 400)
    xy = np.column_stack((t, np.sin(3 * t) + 0.1 * t**2))
    left = np.vstack(([np.nan, np.nan], xy[:-3]))
    return xy[:-2], xy[1:-1], left, xy[2:]


@pytest.fixture
def fake_cupy(monkeypatch: pytest.MonkeyPatch) -> ModuleType:
    """A stand-in ``cupy`` in ``sys.modules`` whose ``ndarray`` is a dummy class."""
    fake = ModuleType("cupy")
    fake.ndarray = type("ndarray", (), {})
    monkeypatch.setitem(sys.modules, "cupy", fake)
    return fake


def test_host_arrays_compute_in_numpy_itself(monkeypatch: pytest.MonkeyPatch) -> None:
    host = _layer()
    monkeypatch.delitem(sys.modules, "cupy", raising=False)
    assert array_namespace(*host) is np
    fake = ModuleType("cupy")
    fake.ndarray = type("ndarray", (), {})
    monkeypatch.setitem(sys.modules, "cupy", fake)
    assert array_namespace(*host) is np


def test_device_arrays_compute_in_their_own_namespace(fake_cupy: ModuleType) -> None:
    assert array_namespace(fake_cupy.ndarray(), fake_cupy.ndarray()) is fake_cupy


def test_mixed_host_and_device_arrays_raise(fake_cupy: ModuleType) -> None:
    with pytest.raises(TypeError):
        array_namespace(np.zeros(2), fake_cupy.ndarray())


@pytest.mark.parametrize(
    "kernel", [curvature_area_batch, chord_deviation_batch], ids=lambda k: k.__name__
)
def test_kernel_matches_numpy_on_the_device(kernel: Callable[..., np.ndarray]) -> None:
    cp = pytest.importorskip("cupy")
    try:
        devices = cp.cuda.runtime.getDeviceCount()
    except Exception:  # noqa: BLE001 - any CUDA runtime failure means no device
        devices = 0
    if devices < 1:
        pytest.skip("CuPy is installed but no CUDA device is available")

    host = _layer()
    expected = kernel(*host)
    assert type(expected) is np.ndarray  # CuPy imported: host stays host
    out = kernel(*(cp.asarray(a) for a in host))
    assert isinstance(out, cp.ndarray)
    # CuPy's pow and hypot differ from numpy's in the last bit: parity, never identity.
    np.testing.assert_allclose(cp.asnumpy(out), expected, rtol=1e-6, atol=0)
    with pytest.raises(TypeError):
        kernel(cp.asarray(host[0]), *host[1:])

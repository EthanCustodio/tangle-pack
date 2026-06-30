"""GPU backend equivalence (skipped cleanly when CuPy is unavailable).

When CuPy and a device are present, routing the batched map through the GPU must
grow the same manifold as the CPU path. The test also checks the always-present
API surface: ``enable_gpu`` raises a clear error if CuPy cannot be imported.
"""

from __future__ import annotations

import numpy as np
import pytest

import tanglepack
from tanglepack import DynamicalSystem


def _map(point):
    x, y = point
    return np.stack([y - 10 + x**2, -x], axis=0)


def _imap(point):
    x, y = point
    return np.stack([-y, x + 10 - y**2], axis=0)


def _jac(point):
    x, y = point
    return np.array([[2 * x, 1], [-1, 0]])


def _grow():
    wb = tanglepack.TangleWorkbench(_map, _imap, _jac)
    fp = wb.construct_fixed_point([4, -4])
    wb.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    wb.initialize_both_manifolds(fp)
    wb.grow_n_times(fp, "unstable", num_iterations=5)
    return wb, fp


def test_enable_gpu_without_cupy_raises_clearly():
    try:
        import cupy  # noqa: F401
    except Exception:
        system = DynamicalSystem(_map, _imap)
        with pytest.raises(ImportError, match="CuPy"):
            tanglepack.enable_gpu(system)
    else:
        pytest.skip("CuPy is installed; the no-CuPy error path is not exercised.")


def test_gpu_growth_matches_cpu():
    pytest.importorskip("cupy")

    wb_cpu, fp_cpu = _grow()
    cd_cpu = np.asarray(
        wb_cpu.manifolds[(fp_cpu, "unstable", 0, 0)].get_cdist_array()
    ).ravel()

    wb_gpu = tanglepack.TangleWorkbench(_map, _imap, _jac)
    tanglepack.enable_gpu(wb_gpu, min_batch_points=1)
    fp = wb_gpu.construct_fixed_point([4, -4])
    wb_gpu.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    wb_gpu.initialize_both_manifolds(fp)
    wb_gpu.grow_n_times(fp, "unstable", num_iterations=5)
    cd_gpu = np.asarray(
        wb_gpu.manifolds[(fp, "unstable", 0, 0)].get_cdist_array()
    ).ravel()

    assert len(cd_cpu) == len(cd_gpu)
    assert np.allclose(cd_cpu, cd_gpu, rtol=1e-6, atol=1e-9)

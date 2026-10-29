"""The optional GPU backend: parity with the CPU path, and a clear error without CuPy.

When CuPy and a device are present, routing the batched map through the GPU
(``tanglepack.enable_gpu``) must grow the same manifold as the CPU path,
``disable_gpu`` must restore the CPU callables, and the same-stability sweep
on the device must list exactly the CPU's hits; these tests skip cleanly
otherwise. The always-present API surface runs everywhere: ``enable_gpu``
raises ``ImportError`` when CuPy cannot be imported (forced by hiding the
module, so the path runs even where CuPy is installed) and ``TypeError`` on a
target that holds no dynamical system.
"""

from __future__ import annotations

import sys

import numpy as np
import pytest

import tanglepack
from tanglepack import DynamicalSystem
from tanglepack.examples import HENON_K10, henon_jacobian, henon_map, henon_map_inverse
from tanglepack.numerics.diagnostics import count_same_stability_crossings, same_stability_crossings

import cases
from cases import K10_SADDLE_SEED, PERIOD1_ORIENTATION

#: Unstable growth steps compared between the two backends.
GROWTH_STEPS = 5


def _grown_cdists(session: tanglepack.TangleSession) -> np.ndarray:
    """
    Grow the k=10 saddle's unstable manifold and return its cdists.

    Args:
        session: A fresh k=10 session (GPU enabled or not).

    Returns:
        The canonical distances along the first unstable branch.
    """
    fp = session.construct_fixed_point(list(K10_SADDLE_SEED))
    session.orient_eigenvectors(
        fp, {key: np.array(value) for key, value in PERIOD1_ORIENTATION.items()}
    )
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=GROWTH_STEPS)
    return np.asarray(session.manifolds[(fp, "unstable", 0, 0)].get_cdist_array()).ravel()


def _k10_session() -> tanglepack.TangleSession:
    """A fresh k=10 session with the analytic Jacobian."""
    return tanglepack.TangleSession(
        henon_map(*HENON_K10),
        henon_map_inverse(*HENON_K10),
        henon_jacobian(*HENON_K10),
        area_cutoff=1e-4,
    )


def _require_device() -> None:
    """Skip unless CuPy imports and sees a CUDA device."""
    cp = pytest.importorskip("cupy")
    try:
        devices = cp.cuda.runtime.getDeviceCount()
    except Exception:  # noqa: BLE001 - any CUDA runtime failure means no device
        devices = 0
    if devices < 1:
        pytest.skip("CuPy is installed but no CUDA device is available")


def test_gpu_growth_matches_cpu() -> None:
    _require_device()
    cpu = _grown_cdists(_k10_session())

    session = _k10_session()
    system = session.workbench.dynamical_system
    cpu_map = system.map
    assert tanglepack.enable_gpu(session, min_batch_points=1) is system
    assert system.map is not cpu_map
    gpu = _grown_cdists(session)

    assert len(cpu) == len(gpu)
    assert np.allclose(cpu, gpu, rtol=1e-6, atol=1e-9)
    tanglepack.disable_gpu(session)
    assert system.map is cpu_map


def test_enable_gpu_without_cupy_raises_clearly(monkeypatch: pytest.MonkeyPatch) -> None:
    """With CuPy unimportable the error is an ImportError and the map stays the CPU one."""
    monkeypatch.setitem(sys.modules, "cupy", None)  # ``import cupy`` now raises
    system = DynamicalSystem(henon_map(*HENON_K10), henon_map_inverse(*HENON_K10))
    cpu_map = system.map
    with pytest.raises(ImportError):
        tanglepack.enable_gpu(system)
    assert system.map is cpu_map


def test_enable_gpu_rejects_a_target_without_a_system() -> None:
    with pytest.raises(TypeError):
        tanglepack.enable_gpu(object())


def test_gpu_same_stability_readings_are_the_cpus() -> None:
    """The listings (in order) and the counts are equal on the CPU and on the
    device, on the k28 series blasted to its depth, its last bridges counted
    alone too."""
    _require_device()
    series = cases.SERIES["k28_series"]
    case = series.build()
    zone = max(case.zones, key=lambda zone: zone.area)
    for _ in range(series.blasts):
        cases._blast(case.session, zone, [zone.fixed_point], diagnose=False,
                     min_separation=series.min_separation)
    workbench = case.workbench
    readings = []
    for gpu in (False, True):
        if gpu:
            tanglepack.enable_gpu(case.session)
        readings.append((
            same_stability_crossings(workbench, region="all"),
            same_stability_crossings(workbench),
            count_same_stability_crossings(workbench, workbench.bridges),
            count_same_stability_crossings(workbench, workbench.bridges[-50:]),
        ))
    tanglepack.disable_gpu(case.session)
    cpu, device = readings
    assert cpu[0] and cpu[2].in_region, "the case must hold hits for the test to bite"
    assert device == cpu

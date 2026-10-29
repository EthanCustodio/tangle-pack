"""The ``workers`` option runs on several cores and changes nothing computed.

``workers`` splits a large map batch across a thread or forked process pool,
tests the same-stability sweep's pair blocks on threads and passes the cores
to the blast's separation queries. For an elementwise map every result is the
serial one, number for number: each test compares workers 1 with workers 4.
These are execution properties, not physical laws, so they live in the unit
tier. The builds are small, so where a test must see the map split it lowers
``MIN_PARALLEL_POINTS`` (the size below which a batch stays in the caller).
"""

from __future__ import annotations

import functools
import importlib
import threading

import numpy as np
import pytest

from tanglepack import DynamicalSystem, TangleSession, TangleWorkbench
from tanglepack.examples import HENON_K10, henon_map, henon_map_inverse
from tanglepack.examples.duffing import duffing_map, duffing_map_inverse
from tanglepack.numerics.diagnostics import (
    count_same_stability_crossings,
    same_stability_crossings,
)

import cases

_system_module = importlib.import_module("tanglepack.numerics.DynamicalSystem")

#: The worker count every test compares with the serial run.
WORKERS = 4


def _scalar_only(fn):
    """Wrap a map so it rejects anything that is not a single ``(2,)`` point."""

    def wrapped(point):
        point = np.asarray(point, dtype=float)
        if point.shape != (2,):
            raise ValueError("scalar-only map")
        return fn(point)

    return wrapped


def _thread_names(fn, names: set):
    """Wrap a map so every call records the name of the thread it ran on."""

    def wrapped(point):
        names.add(threading.current_thread().name)
        return fn(point)

    return wrapped


MAPS = {
    "henon": (henon_map(*HENON_K10), henon_map_inverse(*HENON_K10)),
    "duffing": (duffing_map(), duffing_map_inverse()),
}


@pytest.fixture
def split_small(monkeypatch):
    """Split every batch of the builds below, however small."""
    monkeypatch.setattr(_system_module, "MIN_PARALLEL_POINTS", 16)
    monkeypatch.setattr(_system_module, "MIN_PARALLEL_POINTS_EXPENSIVE", 16)


@pytest.mark.parametrize("executor", ["thread", "process"])
@pytest.mark.parametrize("scalar", [False, True], ids=["batched", "scalar_only"])
@pytest.mark.parametrize("name", sorted(MAPS))
def test_split_batch_equals_serial(name, scalar, executor):
    """``map_batch`` / ``map_inv_batch`` are equal at workers 1 and 4."""
    forward, backward = MAPS[name]
    if scalar:
        forward, backward = _scalar_only(forward), _scalar_only(backward)
    size = (
        _system_module.MIN_PARALLEL_POINTS_EXPENSIVE if scalar else _system_module.MIN_PARALLEL_POINTS
    ) + 37
    coords = np.random.default_rng(0).uniform(-1.5, 1.5, (size, 2))
    serial = DynamicalSystem(forward, backward)
    split = DynamicalSystem(forward, backward, workers=WORKERS, map_executor=executor)
    assert np.array_equal(split.map_batch(coords), serial.map_batch(coords))
    assert np.array_equal(split.map_inv_batch(coords), serial.map_inv_batch(coords))


def test_small_batches_stay_in_the_caller():
    """A batch below the threshold runs on the calling thread; a large one does not."""
    names: set = set()
    forward, backward = MAPS["henon"]
    system = DynamicalSystem(_thread_names(forward, names), backward, workers=WORKERS)
    caller = threading.current_thread().name
    system.map_batch(np.zeros((_system_module.MIN_PARALLEL_POINTS - 1, 2)))
    assert names == {caller}
    names.clear()
    system.map_batch(np.zeros((_system_module.MIN_PARALLEL_POINTS, 2)))
    assert caller not in names and len(names) > 1


def test_setters_rebuild_the_pool():
    """Changing ``workers`` or ``map_executor`` maps the next batch in a new pool."""
    names: set = set()
    forward, backward = MAPS["henon"]
    system = DynamicalSystem(_thread_names(forward, names), backward, workers=WORKERS)
    coords = np.zeros((_system_module.MIN_PARALLEL_POINTS, 2))
    system.map_batch(coords)
    first_pool = {name.rsplit("_", 1)[0] for name in names}
    names.clear()
    system.workers = 2
    system.map_batch(coords)
    assert len(names) <= 2
    assert {name.rsplit("_", 1)[0] for name in names}.isdisjoint(first_pool)

    system.map_executor = "process"
    assert np.array_equal(system.map_batch(coords), DynamicalSystem(forward, backward).map_batch(coords))
    with pytest.raises(ValueError):
        system.map_executor = "gpu"
    with pytest.raises(ValueError):
        system.workers = 0


def test_session_and_workbench_forward_the_setting():
    """The session's and the workbench's settings are the system's."""
    session = TangleSession(*MAPS["henon"], workers=WORKERS, map_executor="process")
    system = session.workbench.dynamical_system
    assert (session.workers, session.map_executor) == (system.workers, system.map_executor)
    assert (system.workers, system.map_executor) == (WORKERS, "process")
    session.workers, session.map_executor = 2, "thread"
    assert (system.workers, system.map_executor) == (2, "thread")


def test_load_gives_one_worker(tmp_path):
    """``workers`` is a runtime resource: a loaded workbench has ``workers=1``."""
    case = cases.build_k10(through="bridges")
    case.session.workers = WORKERS
    case.workbench.save(tmp_path / "k10")
    loaded = TangleWorkbench.load(tmp_path / "k10", *MAPS["henon"])
    assert loaded.workers == 1


def _numbers(workbench: TangleWorkbench) -> tuple:
    """Every point, crossing and iterate link of a workbench, as comparable values."""
    store = workbench.point_store
    registry = workbench.intersection_registry
    crossings = [
        (iid, ix.unstable_cdist, ix.stable_cdist, tuple(ix.coords if ix.coords is not None else ()),
         ix.crossing_sign, ix.resolution_ratio)
        for iid, ix in sorted(registry, key=lambda item: item[0])
    ]
    bridges = sorted(
        (bridge.id or (-1, -1), np.asarray(bridge.get_point_array()).tobytes())
        for bridge in workbench.bridges
    )
    return (
        store.xy[: len(store)].tobytes(),
        store.cdist[: len(store)].tobytes(),
        crossings,
        sorted(registry.iterate_table.items()),
        sorted((s, sorted(c)) for s, c in registry.unresolved_links.items()),
        bridges,
    )


def test_k10_build_is_equal(monkeypatch, split_small):
    """A whole k10 build at workers 4, every batch split, equals the serial build."""
    serial = cases.build_k10()
    monkeypatch.setattr(cases, "TangleSession", functools.partial(TangleSession, workers=WORKERS))
    split = cases.build_k10()
    assert split.session.workers == WORKERS
    assert _numbers(split.workbench) == _numbers(serial.workbench)
    assert split.session.topology().dynamics.describe() == serial.session.topology().dynamics.describe()


def test_same_stability_readings_are_equal():
    """The listing (in order) and the count are equal at workers 1 and 4, on
    the k28 series blasted to its depth (the two-blast case holds no hit)."""
    series = cases.SERIES["k28_series"]
    case = series.build()
    zone = max(case.zones, key=lambda zone: zone.area)
    for _ in range(series.blasts):
        cases._blast(case.session, zone, [zone.fixed_point], diagnose=False,
                     min_separation=series.min_separation)
    workbench = case.workbench
    readings = []
    for workers in (1, WORKERS):
        workbench.workers = workers
        readings.append((
            same_stability_crossings(workbench, region="all"),
            same_stability_crossings(workbench),
            count_same_stability_crossings(workbench, workbench.bridges),
        ))
    serial, split = readings
    assert serial[0] and serial[2].in_region, "the case must hold hits for the test to bite"
    assert split == serial

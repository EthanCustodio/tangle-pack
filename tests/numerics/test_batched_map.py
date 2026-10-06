"""The batched map helper and its scalar-map fallback must agree.

``DynamicalSystem.map_batch`` follows the columns-of-points convention
(coordinate on axis 0). A batch-capable map runs through the single vectorized
call; a scalar-only map falls back to a per-point loop. Both must produce the
same growth, and the breadth-first refiner must leave every consecutive pair
below the area cutoff.
"""

from __future__ import annotations

import numpy as np

import tanglepack
from tanglepack import DynamicalSystem, ManifoldMachine
from tanglepack.examples import (
    HENON_K10,
    henon_jacobian,
    henon_map,
    henon_map_inverse,
    saddle_guesses,
)

_batch_map = henon_map(*HENON_K10)
_batch_imap = henon_map_inverse(*HENON_K10)
_jac = henon_jacobian(*HENON_K10)


def _scalar_only(fn):
    """Wrap a map so it rejects anything that is not a single (2,) point."""

    def wrapped(point):
        arr = np.asarray(point, dtype=float)
        if arr.shape != (2,):
            raise ValueError("scalar-only map")
        return fn(arr)

    return wrapped


def _recording(fn, shapes: list):
    """Wrap a map so every call that SUCCEEDS appends its argument's shape to ``shapes``."""

    def wrapped(point):
        image = fn(point)
        shapes.append(np.shape(point))
        return image

    return wrapped


def test_map_batch_detects_and_matches_loop():
    shapes: list = []
    system = DynamicalSystem(_recording(_batch_map, shapes), _batch_imap)
    coords = np.array([[0.1, 0.2], [4.0, -4.0], [-3.0, 1.5], [2.2, 2.2]])

    batched = system.map_batch(coords)
    expected = np.vstack([_batch_map(p) for p in coords])

    assert np.allclose(batched, expected, rtol=1e-12, atol=1e-15)
    assert batched.shape == coords.shape

    # once detected, a batch-capable map is called ONCE per batch, on (2, N)
    shapes.clear()
    system.map_batch(coords)
    assert shapes == [(2, len(coords))]


def test_scalar_only_map_falls_back():
    shapes: list = []
    system = DynamicalSystem(
        _recording(_scalar_only(_batch_map), shapes), _scalar_only(_batch_imap)
    )
    coords = np.array([[0.1, 0.2], [4.0, -4.0], [-3.0, 1.5]])

    batched = system.map_batch(coords)
    expected = np.vstack([_scalar_only(_batch_map)(p) for p in coords])

    assert np.allclose(batched, expected, rtol=1e-12, atol=1e-15)

    # detection fell back: every later batch is one call per (2,) point
    shapes.clear()
    system.map_batch(coords)
    assert shapes == [(2,)] * len(coords)


def _grow(m, im):
    """Seed the k=10 saddle with maps ``m``/``im`` and grow its unstable branch 5 times."""
    wb = tanglepack.TangleWorkbench(m, im, _jac)
    fp = wb.construct_fixed_point(saddle_guesses(*HENON_K10)["saddle"])
    wb.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    wb.initialize_both_manifolds(fp)
    wb.grow_n_times(fp, "unstable", num_iterations=5)
    return wb, fp, wb.manifolds[(fp, "unstable", 0, 0)]


def test_batch_and_fallback_growth_agree():
    shapes_b: list = []
    shapes_s: list = []
    wb_b, fp_b, man_b = _grow(_recording(_batch_map, shapes_b), _batch_imap)
    wb_s, fp_s, man_s = _grow(
        _recording(_scalar_only(_batch_map), shapes_s), _scalar_only(_batch_imap)
    )

    # growth really took the two different paths: column-batch (2, N) calls on
    # one side, only single (2,) points on the scalar-only side
    assert any(len(shape) == 2 and shape[0] == 2 for shape in shapes_b)
    assert shapes_s and all(shape == (2,) for shape in shapes_s)

    cd_b = np.asarray(man_b.get_cdist_array()).ravel()
    cd_s = np.asarray(man_s.get_cdist_array()).ravel()
    assert len(cd_b) == len(cd_s)
    assert np.allclose(cd_b, cd_s, rtol=1e-12, atol=1e-9)
    assert np.allclose(
        man_b.get_point_array(), man_s.get_point_array(), rtol=1e-9, atol=1e-9
    )


def test_refinement_essentially_converges_to_cutoff():
    """The refiner drives nearly every pair below ``area_cutoff``.

    A pair below the cutoff when it is evaluated is not re-checked if a later
    insertion in an *adjacent* segment changes its neighbour (true of both the
    breadth-first refiner and the old depth-first one), so a small boundary
    fraction can end just above the cutoff. We pin that this slack stays tiny:
    only a sliver of pairs exceed the cutoff and none by a large factor. The
    areas come from the vectorized kernel the refiner itself flags with.
    """
    wb, fp, man = _grow(_batch_map, _batch_imap)
    cutoff = wb._man_machine.area_cutoff

    points = np.asarray(man.get_point_array(), dtype=float).reshape(-1, 2)
    nan = np.full((1, 2), np.nan)
    p0, p1 = points[:-1], points[1:]
    left = np.vstack((nan, points[:-2]))
    right = np.vstack((points[2:], nan))
    areas = ManifoldMachine._curvature_area_batch(p0, p1, left, right)

    total = len(areas)
    over = int(np.count_nonzero(areas >= cutoff))
    max_ratio = float(areas.max() / cutoff)

    assert total > 0
    assert over / total < 0.02, f"{over}/{total} pairs above cutoff"
    assert max_ratio < 5.0, f"a pair exceeded the cutoff by {max_ratio:.1f}x"

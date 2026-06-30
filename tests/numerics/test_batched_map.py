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
from tanglepack import DynamicalSystem, ManifoldView


def _batch_map(point):
    x, y = point
    return np.stack([y - 10 + x**2, -x], axis=0)


def _batch_imap(point):
    x, y = point
    return np.stack([-y, x + 10 - y**2], axis=0)


def _jac(point):
    x, y = point
    return np.array([[2 * x, 1], [-1, 0]])


def _scalar_only(fn):
    """Wrap a map so it rejects anything that is not a single (2,) point."""

    def wrapped(point):
        arr = np.asarray(point, dtype=float)
        if arr.shape != (2,):
            raise ValueError("scalar-only map")
        return fn(arr)

    return wrapped


def test_map_batch_detects_and_matches_loop():
    system = DynamicalSystem(_batch_map, _batch_imap)
    coords = np.array([[0.1, 0.2], [4.0, -4.0], [-3.0, 1.5], [2.2, 2.2]])

    batched = system.map_batch(coords)
    expected = np.vstack([_batch_map(p) for p in coords])

    assert system._map_batchable is True
    assert np.allclose(batched, expected, rtol=1e-12, atol=1e-15)
    assert batched.shape == coords.shape


def test_scalar_only_map_falls_back():
    system = DynamicalSystem(_scalar_only(_batch_map), _scalar_only(_batch_imap))
    coords = np.array([[0.1, 0.2], [4.0, -4.0], [-3.0, 1.5]])

    batched = system.map_batch(coords)
    expected = np.vstack([_scalar_only(_batch_map)(p) for p in coords])

    assert system._map_batchable is False  # detection fell back
    assert np.allclose(batched, expected, rtol=1e-12, atol=1e-15)


def _grow(m, im):
    wb = tanglepack.TangleWorkbench(m, im, _jac)
    fp = wb.construct_fixed_point([4, -4])
    wb.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    wb.initialize_both_manifolds(fp)
    wb.grow_n_times(fp, "unstable", num_iterations=5)
    return wb, fp, wb.manifolds[(fp, "unstable", 0, 0)]


def test_batch_and_fallback_growth_agree():
    wb_b, fp_b, man_b = _grow(_batch_map, _batch_imap)
    wb_s, fp_s, man_s = _grow(_scalar_only(_batch_map), _scalar_only(_batch_imap))

    assert wb_b._man_machine.system._map_batchable is True
    assert wb_s._man_machine.system._map_batchable is False

    cd_b = np.asarray(man_b.get_cdist_array()).ravel()
    cd_s = np.asarray(man_s.get_cdist_array()).ravel()
    assert len(cd_b) == len(cd_s)
    assert np.allclose(cd_b, cd_s, rtol=1e-12, atol=1e-9)
    assert np.allclose(
        man_b.get_point_array(), man_s.get_point_array(), rtol=1e-9, atol=1e-9
    )


def test_refinement_essentially_converges_to_cutoff():
    """The refiner drives nearly every pair below the cutoff.

    A pair below the cutoff when it is evaluated is not re-checked if a later
    insertion in an *adjacent* segment changes its neighbour (true of both the
    breadth-first refiner and the old depth-first one), so a small boundary
    fraction can end just above the cutoff. We pin that this slack stays tiny:
    only a sliver of pairs exceed the cutoff and none by a large factor.
    """
    wb, fp, man = _grow(_batch_map, _batch_imap)
    machine = wb._man_machine
    viewer = ManifoldView(man, machine.system)
    cutoff = machine.area_cutoff

    total = over = 0
    max_ratio = 0.0
    prev = man.root
    cur = man.walk_fwd(None, prev)
    while cur is not None:
        area = machine._curvature_area((prev, cur), viewer)
        total += 1
        if area >= cutoff:
            over += 1
            max_ratio = max(max_ratio, area / cutoff)
        if cur is man.tail:
            break
        nxt = man.walk_fwd(prev, cur)
        prev, cur = cur, nxt

    assert over / total < 0.02, f"{over}/{total} pairs above cutoff"
    assert max_ratio < 5.0, f"a pair exceeded the cutoff by {max_ratio:.1f}x"

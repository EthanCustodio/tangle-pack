"""Phase 8A: the single manifold walker, the single branch walker and the
shared Hénon example maps.

The duplications were collapsed here and each one is pinned:

* :meth:`BaseManifold._collect` is now the ONE read behind all five array
  getters. The getters are thin wrappers, so the test compares each of them
  against the curve's rows read straight off its chain -- if a wrapper
  ever passes the wrong ``kind`` / ``value`` / ``stop_at_final`` the comparison
  catches it.
* A curve rooted at a BranchPoint reads on from position 0 of the chain its
  slot holds (``BaseManifold.ids``; the slot is read in one place,
  ``_slot_chain``).
* ``tanglepack.examples.henon`` is the one Hénon definition.
"""

from __future__ import annotations


import numpy as np
import pytest

import tanglepack
from tanglepack.examples import (
    HENON_K10,
    HENON_P3,
    henon_jacobian,
    henon_map,
    henon_map_inverse,
    saddle_guesses,
)
from tanglepack.numerics import BaseManifold, BranchPoint, FixedPointSolver


# --------------------------------------------------------------------------- #
# one walker
# --------------------------------------------------------------------------- #
def _hand_walk(manifold, stop_at_final=True):
    """The curve read straight off its chain: the periodic point (or None) and
    the rows from its first row's position to its tail's (to the chain's end
    without the stop)."""
    root = manifold._root
    head = root if isinstance(root, BranchPoint) else None
    if head is None:
        store = manifold.store
        chain, lo = store.chain_of(root), store.position_of(root)
    else:
        chain, lo = head.chain(manifold.stability, manifold.branch_index), 0
        store = chain.store
    ids = chain.ids
    hi = int(np.flatnonzero(ids == manifold._tail)[0]) + 1 if stop_at_final else len(ids)
    return head, store, ids[lo:hi]


def _xy(head, store, rows):
    return np.vstack(([] if head is None else [head.get_point()]) + [store.xy[r] for r in rows.tolist()])


def _cd(head, store, rows, stability):
    return np.vstack(([] if head is None else [[head.get_cdist(stability)]]) + [[store.cdist[r]] for r in rows.tolist()])


def _iterate_of(manifold, node, n=1):
    return node.iterate(n if manifold.stability == "unstable" else -n)


@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_every_getter_matches_an_independent_walk(grown_both, stability):
    """All five getters are one walk with different filters."""
    workbench, fp = grown_both
    manifold = workbench.manifolds[(fp, stability, 0, 0)]

    head, store, rows = _hand_walk(manifold)
    _all_head, _store, all_rows = _hand_walk(manifold, stop_at_final=False)
    assert len(rows) > 10, "fixture too small to be a real walk"

    assert np.array_equal(manifold.get_point_array(), _xy(head, store, rows))
    np.testing.assert_array_equal(manifold.ids, rows)

    # the cdist read deliberately ignores the tail (documented on the getter)
    assert np.array_equal(manifold.get_cdist_array(), _cd(head, store, all_rows, stability))

    link = store.image if manifold.stability == "unstable" else store.preimage
    images = link[rows]
    # A periodic point always has its image, the next orbit point, so the
    # non-iterated read never holds the head.
    non_iterated = rows[images < 0]
    assert np.array_equal(manifold.get_non_iterated_point_array(), _xy(None, store, non_iterated))
    expected = _cd(None, store, non_iterated, stability)
    assert np.array_equal(
        manifold.get_non_iterated_cdist_array(), expected if len(expected) else np.array([])
    )

    image_head = None if head is None else _iterate_of(manifold, head)
    assert np.array_equal(manifold.get_iterated_point_array(), _xy(image_head, store, images[images >= 0]))


def test_final_node_truncates_the_point_walk(grown_unstable):
    """``final_node`` stops the walk; it is the same walker, bounded earlier."""
    workbench, fp, manifold = grown_unstable
    head, store, rows = _hand_walk(manifold)
    stop = int(rows[len(rows) // 2])

    assert np.array_equal(
        manifold.get_point_array(final_node=stop), _xy(head, store, rows[: len(rows) // 2 + 1])
    )


def test_collect_is_the_only_traversal(grown_unstable):
    """The getters delegate: break ``_collect`` and every one of them breaks."""
    workbench, fp, manifold = grown_unstable

    class _CollectReached(RuntimeError):
        """Raised only by the sabotaged ``_collect``: the getter went through it."""

    def _boom(*args, **kwargs):
        raise _CollectReached

    manifold._collect = _boom
    for call in (
        manifold.get_point_array,
        manifold.get_cdist_array,
        manifold.get_non_iterated_point_array,
        manifold.get_non_iterated_cdist_array,
        manifold.get_iterated_point_array,
    ):
        with pytest.raises(_CollectReached):
            call()


# --------------------------------------------------------------------------- #
# one branch walker
# --------------------------------------------------------------------------- #
def test_first_node_steps_past_the_root_branch_point(initialized):
    """A manifold anchored on a BranchPoint starts at position 0 of the chain
    its slot holds."""
    workbench, fp = initialized

    root = fp.branch_points[0]
    assert isinstance(root, BranchPoint)

    view = BaseManifold(
        root,
        "unstable",
        stretch_param=1,
        fixed_point=fp,
        branch_index=0,
        manifold_key=(fp, "unstable", 0, 0),
    )
    assert len(view.ids) and view.ids[0] == root.forward_branches[0].ids[0]
    for branch, chain in enumerate(root.forward_branches):
        view.branch_index = branch
        assert view.ids.tolist() == ([] if chain is None else chain.ids.tolist())


def test_first_node_on_an_ordinary_root_is_the_root(initialized):
    """No BranchPoint, no step: the helper is a no-op on an ordinary root."""
    workbench, fp = initialized
    manifold = workbench.manifolds[(fp, "unstable", 0, 0)]

    ordinary = int(manifold.ids[0])
    view = BaseManifold(
        ordinary,
        "unstable",
        stretch_param=1,
        fixed_point=fp,
        branch_index=0,
        manifold_key=(fp, "unstable", 0, 0),
        store=manifold.store,
    )
    assert view.ids[0] == ordinary


def test_first_node_needs_a_branch_index(initialized):
    """Leaving a BranchPoint is ambiguous without one; say so instead of guessing."""
    workbench, fp = initialized
    view = BaseManifold(
        fp.branch_points[0],
        "unstable",
        stretch_param=1,
        fixed_point=fp,
        branch_index=0,
        manifold_key=(fp, "unstable", 0, 0),
    )
    view.branch_index = None
    with pytest.raises(ValueError):
        view.ids


# --------------------------------------------------------------------------- #
# one Hénon definition
# --------------------------------------------------------------------------- #
def test_henon_factories_match_the_closed_form():
    k, b = 2.8, 1
    fwd, inv, jac = henon_map(k, b), henon_map_inverse(k, b), henon_jacobian(k, b)
    point = np.array([0.7, -1.3])

    x, y = point
    assert np.allclose(fwd(point), [y - k + x**2, -b * x])
    assert np.allclose(inv(fwd(point)), point)
    assert np.allclose(jac(point), [[2 * x, 1], [-b, 0]])
    # area preservation: |det J| == |b| everywhere
    assert np.isclose(np.linalg.det(jac(point)), b)


def test_henon_factories_are_batch_capable_on_axis_zero():
    """The numerics layer maps a whole refinement layer in one call."""
    fwd = henon_map(*HENON_K10)
    batch = np.array([[0.1, 4.0, -3.0], [0.2, -4.0, 1.5]])

    mapped = fwd(batch)
    assert mapped.shape == batch.shape
    for i in range(batch.shape[1]):
        assert np.allclose(mapped[:, i], fwd(batch[:, i]))


@pytest.mark.parametrize(
    "params, name",
    [(HENON_K10, "saddle"), (HENON_K10, "inversion"),
     (HENON_P3, "period_3"), (HENON_P3, "saddle")],
)
def test_saddle_guesses_converge_to_real_fixed_points(params, name):
    """Every recorded guess converges to a real orbit of its own minimal period."""
    guess = saddle_guesses(*params)[name]
    workbench = tanglepack.TangleWorkbench(
        henon_map(*params), henon_map_inverse(*params), henon_jacobian(*params)
    )
    fp = workbench.construct_fixed_point(guess)

    assert fp.period == len(np.atleast_2d(guess))
    assert FixedPointSolver.minimal_period(fp.coordinates) == fp.period
    orbit = np.asarray(fp.coordinates, dtype=float)
    images = np.array([np.ravel(henon_map(*params)(point))[:2] for point in orbit])
    assert np.allclose(images, np.roll(orbit, -1, axis=0), atol=1e-9)


def test_saddle_guesses_returns_a_fresh_copy():
    """A caller that mutates the guess must not poison the next caller."""
    first = saddle_guesses(*HENON_K10)
    first["saddle"][0] = 99
    assert saddle_guesses(*HENON_K10)["saddle"] == [4, -4]


def test_saddle_guesses_refuses_unknown_parameters():
    with pytest.raises(KeyError):
        saddle_guesses(3.7, 1)



"""Phase 8A: the single manifold walker, the single branch walker and the
shared Hénon example maps.

The duplications were collapsed here and each one is pinned:

* :meth:`BaseManifold._collect` is now the ONE traversal behind all five array
  getters. The getters are thin wrappers, so the test compares each of them
  against an independent hand-rolled walk of the linked list -- if a wrapper
  ever passes the wrong ``kind`` / ``value`` / ``stop_at_final`` the comparison
  catches it.
* :meth:`BaseManifold.first_node` is the ONE "step past the root BranchPoint"
  helper that ``ManifoldMachine._branch_view`` and
  ``ManifoldInitializer.construct_kevin_way`` had each open-coded.
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
from tanglepack.numerics import BaseManifold, BranchPoint


# --------------------------------------------------------------------------- #
# one walker
# --------------------------------------------------------------------------- #
def _hand_walk(manifold, stop_at_final=True):
    """Walk the linked list without touching BaseManifold's getters."""
    nodes = []
    prev = None
    current = manifold.root
    while current is not None:
        nodes.append(current)
        if stop_at_final and current is manifold.tail:
            break
        prev, current = current, manifold.walk_fwd(
            prev, current, branch_index=manifold.branch_index
        )
    return nodes


def _has_iterate(manifold, node, n=1):
    return (
        node.exists_next_iterate(n)
        if manifold.stability == "unstable"
        else node.exists_prev_iterate(n)
    )


def _iterate_of(manifold, node, n=1):
    return (
        node.get_next_iterate(n)
        if manifold.stability == "unstable"
        else node.get_prev_iterate(n)
    )


@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_every_getter_matches_an_independent_walk(grown_both, stability):
    """All five getters are one walk with different filters."""
    workbench, fp = grown_both
    manifold = workbench.manifolds[(fp, stability, 0, 0)]

    nodes = _hand_walk(manifold)
    all_nodes = _hand_walk(manifold, stop_at_final=False)
    assert len(nodes) > 10, "fixture too small to be a real walk"

    assert np.array_equal(
        manifold.get_point_array(), np.vstack([n.get_point() for n in nodes])
    )
    assert manifold.get_point_array(return_nodes=True) == nodes

    # the cdist walk deliberately ignores the tail (documented on the getter)
    assert np.array_equal(
        manifold.get_cdist_array(), np.vstack([n.cdist for n in all_nodes])
    )
    assert manifold.get_cdist_array(return_nodes=True) == all_nodes

    expected_non_iterated = [n for n in nodes if not _has_iterate(manifold, n)]
    assert manifold.get_non_iterated_point_array(return_nodes=True) == (
        expected_non_iterated
    )
    assert np.array_equal(
        manifold.get_non_iterated_cdist_array(),
        (
            np.vstack([n.cdist for n in expected_non_iterated])
            if expected_non_iterated
            else np.array([])
        ),
    )

    expected_iterates = [
        _iterate_of(manifold, n) for n in nodes if _has_iterate(manifold, n)
    ]
    assert manifold.get_iterated_point_array(return_nodes=True) == expected_iterates


def test_final_node_truncates_the_point_walk(grown_unstable):
    """``final_node`` stops the walk; it is the same walker, bounded earlier."""
    workbench, fp, manifold = grown_unstable
    nodes = _hand_walk(manifold)
    stop = nodes[len(nodes) // 2]

    assert manifold.get_point_array(final_node=stop, return_nodes=True) == (
        nodes[: nodes.index(stop) + 1]
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


def test_collect_rejects_a_broken_iterate_link(grown_unstable):
    """A node that claims an iterate it does not have is an error, not a None row."""
    workbench, fp, manifold = grown_unstable
    node = next(n for n in _hand_walk(manifold) if _has_iterate(manifold, n))

    class _Liar:
        def exists_next_iterate(self, num_iterates=1):
            return True

        def get_next_iterate(self, num_iterates=1):
            return None

    node.__class__ = type("LyingPoint", (_Liar, type(node)), {})
    try:
        with pytest.raises(ValueError):
            manifold.get_iterated_point_array()
    finally:
        node.__class__ = type(node).__mro__[2]


# --------------------------------------------------------------------------- #
# one branch walker
# --------------------------------------------------------------------------- #
def test_first_node_steps_past_the_root_branch_point(initialized):
    """A manifold anchored on a BranchPoint starts one walk_fwd step past it."""
    workbench, fp = initialized
    manifold = workbench.manifolds[(fp, "unstable", 0, 0)]

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
    assert view.first_node() is view.walk_fwd(None, root, 0)
    assert view.first_node() is root.forward_branches[0]
    assert view.first_node(1) is root.forward_branches[1]


@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_walk_back_inverts_walk_fwd_through_the_root(initialized, stability):
    """``walk_back(walk_fwd(prev, node), node) is prev`` everywhere, the root too.

    At an ordinary node this is the linked list read both ways; at the root
    BranchPoint the walk enters on one slot and leaves on the matching slot, and
    walking back must undo exactly that toggle.
    """
    workbench, fp = initialized
    manifold = workbench.manifolds[(fp, stability, 0, 0)]
    root = fp.branch_points[0]
    slots = root.backward_branches if stability == "unstable" else root.forward_branches
    entries = [node for node in slots if node is not None]  # unseeded slots are None
    assert entries, "the seeded root must have an entry slot"

    checked = 0
    for prev in entries:
        out = manifold.walk_fwd(prev, root)
        assert manifold.walk_back(out, root) is prev
        checked += 1

    prev, cur = root, manifold.first_node()
    while cur is not None and cur is not manifold.tail:
        nxt = manifold.walk_fwd(prev, cur)
        assert manifold.walk_back(nxt, cur) is prev
        prev, cur = cur, nxt
        checked += 1
    assert checked > len(entries), "the walk must cover ordinary nodes too"


def test_first_node_on_an_ordinary_root_is_the_root(initialized):
    """No BranchPoint, no step: the helper is a no-op on an ordinary root."""
    workbench, fp = initialized
    manifold = workbench.manifolds[(fp, "unstable", 0, 0)]

    ordinary = manifold.first_node()
    view = BaseManifold(
        ordinary,
        "unstable",
        stretch_param=1,
        fixed_point=fp,
        branch_index=0,
        manifold_key=(fp, "unstable", 0, 0),
    )
    assert view.first_node() is ordinary


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
        view.first_node()


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
    [(HENON_K10, "saddle"), (HENON_K10, "inversion"), (HENON_P3, "saddle")],
)
def test_saddle_guesses_converge_to_real_fixed_points(params, name):
    guess = saddle_guesses(*params)[name]
    workbench = tanglepack.TangleWorkbench(
        henon_map(*params), henon_map_inverse(*params), henon_jacobian(*params)
    )
    fp = workbench.construct_fixed_point(guess)

    coord = np.asarray(fp.coordinates[0], dtype=float).ravel()[:2]
    image = np.asarray(henon_map(*params)(coord), dtype=float).ravel()[:2]
    assert np.allclose(image, coord, atol=1e-9)


def test_saddle_guesses_returns_a_fresh_copy():
    """A caller that mutates the guess must not poison the next caller."""
    first = saddle_guesses(*HENON_K10)
    first["saddle"][0] = 99
    assert saddle_guesses(*HENON_K10)["saddle"] == [4, -4]


def test_saddle_guesses_refuses_unknown_parameters():
    with pytest.raises(KeyError):
        saddle_guesses(3.7, 1)



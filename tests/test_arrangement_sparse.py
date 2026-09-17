"""
The sparse mode of :class:`~tanglepack.topology.Arrangement.Arrangement`.

A minimal trellis carries a SUBSET of the bridges. Building the dense
arrangement over such a snapshot manufactures an unstable arc between any two
kept crossings that are consecutive on a branch, whether or not a bridge was
kept there, and stubs every absent unstable ray as a slit. Sparse mode takes
the bridge list literally and lets the sectors either side of an absent ray
merge. These tests pin that on a hand-built line crossed by three lobes.
"""

from __future__ import annotations

import numpy as np

from tanglepack.numerics.Bridge import Bridge
from tanglepack.numerics.Intersection import Intersection
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.numerics.Point import Point
from tanglepack.topology.Arrangement import Arrangement
from tanglepack.topology.Trellis import Trellis
from tanglepack.topology.TrellisBranch import TrellisBranch


class _FakeFixedPoint:
    """The bare minimum a TrellisBranch reads off a fixed point."""

    period = 1
    k_value = 1
    coordinates = np.array([[0.0, 0.0]])
    unstable_eigenvalues = [2.0]


def _lobe(fp, u_key, first, second, y0, bulge):
    """A three-point bridge leaving the line at y0 and returning at y0 + 1."""
    root = Point(0.0, y0 - 0.1, 5.0 * (y0 - 0.1))
    mid = Point(bulge, y0 + 0.5, 5.0 * (y0 + 0.5))
    tail = Point(0.0, y0 + 1.1, 5.0 * (y0 + 1.1))
    root.forward, mid.backward = mid, root
    mid.forward, tail.backward = tail, mid
    return Bridge(
        root=root,
        stability="unstable",
        stretch_param=1.0,
        fixed_point=fp,
        tail=tail,
        branch_index=0,
        manifold_key=u_key,
        first_intersection=first,
        second_intersection=second,
    )


def _three_lobes(keep):
    """
    One stable line ``x = 0`` crossed at y = 0, 1, 2, 3 by one unstable curve.

    The unstable branch visits ``a, b, c, d`` in that order, bulging right,
    left, right (signs alternate ``+ - + -``). ``keep`` selects which of the
    three lobes ``(a,b), (b,c), (c,d)`` are handed to the trellis as bridges;
    every crossing stays on both branches regardless.
    """
    fp = _FakeFixedPoint()
    u_key = (fp, "unstable", 0, 0)
    s_key = (fp, "stable", 0, 0)
    registry = IntersectionRegistry()
    ids = []
    for index, sign in enumerate((1, -1, 1, -1)):
        ids.append(
            registry.add(
                Intersection.synthetic(
                    coords=(0.0, float(index)),
                    unstable_cdist=5.0 * index,
                    stable_cdist=float(index),
                    manifold_a_key=u_key,
                    manifold_b_key=s_key,
                    crossing_sign=sign,
                )
            )
        )
    a, b, c, d = ids
    lobes = {
        (a, b): _lobe(fp, u_key, a, b, 0.0, 1.0),
        (b, c): _lobe(fp, u_key, b, c, 1.0, -1.0),
        (c, d): _lobe(fp, u_key, c, d, 2.0, 1.0),
    }
    names = {"ab": (a, b), "bc": (b, c), "cd": (c, d)}
    bridges = [lobes[names[name]] for name in keep]
    branches = {
        u_key: TrellisBranch(
            key=u_key,
            fixed_point=fp,
            stability="unstable",
            orbit_index=0,
            branch_index=0,
            intersection_ids=list(ids),
        ),
        s_key: TrellisBranch(
            key=s_key,
            fixed_point=fp,
            stability="stable",
            orbit_index=0,
            branch_index=0,
            intersection_ids=list(ids),
        ),
    }
    trellis = Trellis(
        fixed_points=[fp], registry=registry, branches=branches, bridges=bridges
    )
    return trellis, (a, b, c, d)


def _unstable_pairs(arrangement):
    return {
        (arc.lo_id, arc.hi_id)
        for face in arrangement.faces
        for arc in face.arcs
        if arc.kind == "unstable"
    }


def test_dense_build_over_all_three_lobes_has_three_regions():
    trellis, _ = _three_lobes(("ab", "bc", "cd"))
    arrangement = Arrangement.from_trellis(trellis)
    assert len(arrangement.regions) == 3
    assert arrangement.euler_characteristic == 2


def test_sparse_unstable_arcs_are_exactly_the_kept_bridges():
    trellis, (a, b, c, d) = _three_lobes(("ab", "cd"))
    arrangement = Arrangement.from_trellis(trellis, sparse=True)
    assert _unstable_pairs(arrangement) == {(a, b), (c, d)}
    assert arrangement.sparse
    assert "sparse" in arrangement.summary()


def test_dense_build_over_the_subset_manufactures_the_missing_arc():
    """Why sparse mode exists: the dense build fills the gap silently."""
    trellis, (a, b, c, d) = _three_lobes(("ab", "cd"))
    arrangement = Arrangement.from_trellis(trellis)
    assert (b, c) in _unstable_pairs(arrangement)


def test_sparse_kept_endpoints_are_degree_three_without_unstable_stubs():
    trellis, (a, b, c, d) = _three_lobes(("ab", "cd"))
    arrangement = Arrangement.from_trellis(trellis, sparse=True)
    assert set(arrangement._nodes[b].slots) == {"s+", "s-", "u-"}
    assert set(arrangement._nodes[c].slots) == {"s+", "s-", "u+"}
    # The stable ends of the line still dangle.
    assert set(arrangement._nodes[a].slots) == {"s+", "s-", "u+"}
    assert arrangement._nodes[arrangement._half_edges[arrangement._nodes[a].slots["s-"]].head].virtual
    assert set(arrangement._nodes[d].slots) == {"s+", "s-", "u-"}
    assert arrangement._nodes[arrangement._half_edges[arrangement._nodes[d].slots["s+"]].head].virtual
    # Nothing virtual hangs off an unstable slot anywhere.
    for node in arrangement._nodes.values():
        if node.virtual:
            continue
        for slot in ("u+", "u-"):
            if slot in node.slots:
                head = arrangement._half_edges[node.slots[slot]].head
                assert not arrangement._nodes[head].virtual


def test_sparse_faces_close_and_euler_holds():
    trellis, (a, b, c, d) = _three_lobes(("ab", "cd"))
    arrangement = Arrangement.from_trellis(trellis, sparse=True)
    # Two closed lobes and one open outer face; the gap between the lobes is
    # part of the outer face, not a slit.
    assert len(arrangement.regions) == 2
    assert len(arrangement.faces) == 3
    assert arrangement.component_count == 1
    assert arrangement.euler_characteristic == 2
    corners = {region.corners for region in arrangement.regions}
    assert corners == {tuple(sorted((a, b))), tuple(sorted((c, d)))} or all(
        set(region.corners) in ({a, b}, {c, d}) for region in arrangement.regions
    )


def test_sparse_and_dense_agree_when_every_bridge_is_kept():
    trellis, _ = _three_lobes(("ab", "bc", "cd"))
    dense = Arrangement.from_trellis(trellis)
    sparse = Arrangement.from_trellis(trellis, sparse=True)
    assert _unstable_pairs(dense) == _unstable_pairs(sparse)
    assert len(dense.regions) == len(sparse.regions) == 3
    # The dense build stubs the two loose unstable ends; sparse does not. A stub
    # adds one vertex and one edge, so the face count is the same either way.
    assert len(sparse.faces) == len(dense.faces)
    assert "+2 dangling ends" in sparse.summary()
    assert "+4 dangling ends" in dense.summary()

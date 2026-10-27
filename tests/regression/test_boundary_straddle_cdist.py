"""Regression: bridge endpoints are real manifold nodes, not offset straddle points.

The bug: ``create_bridges`` used to bound each bridge with two freshly built
"straddle" points, placed at a 10% offset from the crossing toward a segment
endpoint but *assigned the cdist interpolated at the crossing* — i.e. the cdist of
a location the point does not occupy. Because a bridge's cdists are scaled by
``stretch_param`` every time it is iterated, that offset error was baked in and
compounded; the fresh points also laid a second polyline over the existing curve,
which produced overlapping duplicate bridges and spurious unstable x unstable
detections.

The fix: a bridge is now bounded by two points that already exist on its parent
unstable manifold — the real node just below its first crossing and the real node
just above its second — so an endpoint's cdist is, by construction, the cdist of
its own location.

This pins that fix through the public ``create_bridges`` path:

1. root and tail are (by identity) nodes of the parent unstable manifold;
2. they bracket the bridge's two crossings in unstable cdist;
3. they are the *tightest* such nodes — root is the last manifold node at or below
   the first crossing, tail the first at or above the second — which is exactly
   what fails if an endpoint carries the crossing's cdist instead of its own;
4. each endpoint's cdist agrees with the cdist interpolated across its two
   manifold neighbours at its own position, to within the curvature of the arc.

Note:
    Check 4 uses ``rtol=1e-3``: the two-neighbour chord skips the curve between
    them, so the chord estimate differs from the node's true arclength cdist by the
    local curvature (worst case ~4e-4 on this fixture). The old bug displaced the
    endpoint's cdist by ~10% of the segment's cdist span, an order of magnitude
    above that.
"""

from __future__ import annotations

import numpy as np

from helpers.invariants import node_cdist
from tanglepack import BranchPoint, Tangle


def _key(end):
    """An end as a dict key: its row, or the periodic point by identity."""
    return ("periodic", id(end)) if isinstance(end, BranchPoint) else ("row", int(end))


def _node_index_map(nodes) -> dict:
    return {_key(node): i for i, node in enumerate(nodes)}


def test_bridge_endpoints_are_bracketing_manifold_nodes(henon_tangle_with_bridges):
    workbench, fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry
    bridges = workbench.bridges
    assert bridges, "fixture produced no bridges; the test would be vacuous"

    for bridge in bridges:
        manifold = workbench.manifolds[bridge.manifold_key]
        head, store, rows = manifold._rows()
        nodes = ([] if head is None else [head]) + rows.tolist()
        cdists = [node_cdist(n, "unstable") if isinstance(n, BranchPoint) else float(store.cdist[n])
                  for n in nodes]
        index_of = _node_index_map(nodes)

        # (1) the endpoints are nodes of the parent manifold (rows, or its periodic point)
        rows = bridge.ids
        root_key = _key(bridge.root if isinstance(bridge.root, BranchPoint) else rows[0])
        tail_key = _key(rows[-1])
        assert root_key in index_of, (
            "bridge root is not a node of its parent unstable manifold "
            "(a fresh straddle point was created instead)"
        )
        assert tail_key in index_of, (
            "bridge tail is not a node of its parent unstable manifold "
            "(a fresh straddle point was created instead)"
        )
        root_index = index_of[root_key]
        tail_index = index_of[tail_key]
        assert root_index < tail_index, "bridge endpoints are out of order"

        first_cdist = registry[bridge.first_intersection].unstable_cdist
        second_cdist = registry[bridge.second_intersection].unstable_cdist
        root_cdist = node_cdist(bridge.root, "unstable")
        tail_cdist = node_cdist(bridge.tail, "unstable")

        # (2) the endpoints bracket both crossings
        assert root_cdist <= first_cdist, (
            f"root cdist {root_cdist} is above its first crossing {first_cdist}"
        )
        assert first_cdist <= second_cdist, (
            "bridge crossings are not in unstable dynamical order: "
            f"{first_cdist} then {second_cdist}"
        )
        assert tail_cdist >= second_cdist, (
            f"tail cdist {tail_cdist} is below its second crossing {second_cdist}"
        )

        # (3) they are the tightest bracketing nodes
        last_below = max(i for i, c in enumerate(cdists) if c <= first_cdist)
        first_above = min(i for i, c in enumerate(cdists) if c >= second_cdist)
        assert root_index == last_below, (
            f"root is node {root_index} but the last node at or below the first "
            f"crossing is {last_below}"
        )
        assert tail_index == first_above, (
            f"tail is node {tail_index} but the first node at or above the second "
            f"crossing is {first_above}"
        )

        # (4) each endpoint's cdist describes its own location
        for name, node, node_index in (
            ("root", bridge.root, root_index),
            ("tail", bridge.tail, tail_index),
        ):
            if node_index == 0 or node_index == len(nodes) - 1:
                continue  # no two-sided neighbourhood to interpolate across
            chord_cdist = Tangle._cdist_between(
                store,
                nodes[node_index - 1],
                nodes[node_index + 1],
                "unstable",
                store.xy[nodes[node_index]],
            )
            actual = node_cdist(node, "unstable")
            assert np.isclose(actual, chord_cdist, rtol=1e-3), (
                f"{name} cdist {actual} disagrees with the cdist at its own "
                f"location {chord_cdist}"
            )

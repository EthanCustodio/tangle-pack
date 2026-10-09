"""One stable trim: ``TangleWorkbench.trim_stable_at`` (SIMP-10, WP11)."""

from __future__ import annotations

import logging

from cases import build_k10
from helpers.logs import assert_logged


def _first_node_at_or_past(manifold, cdist):
    return next(
        node
        for node in manifold._walk_nodes(stop_at_final=False)
        if node is not manifold.root and node.get_cdist("stable") >= cdist
    )


def _stable_branches(workbench, fp):
    return [key for key in workbench.manifolds if key[0] is fp and key[1] == "stable"]


def test_the_trim_keeps_the_first_node_at_or_past_the_cdist(small_tangle):
    workbench, fp = small_tangle
    key = _stable_branches(workbench, fp)[0]
    manifold = workbench.manifolds[key]
    nodes = list(manifold._walk_nodes(stop_at_final=False))
    lo, hi = nodes[len(nodes) // 3].get_cdist("stable"), nodes[len(nodes) // 3 + 1].get_cdist("stable")
    cdist = 0.5 * (lo + hi)  # strictly between two nodes
    expected = _first_node_at_or_past(manifold, cdist)

    assert workbench.trim_stable_at(key, cdist) is manifold
    assert manifold.tail is expected
    assert expected.get_cdist("stable") >= cdist > lo


def test_a_cdist_on_a_node_keeps_that_node(small_tangle):
    """"At or beyond": a cdist exactly on a node makes that node the tail, not the next."""
    workbench, fp = small_tangle
    key = _stable_branches(workbench, fp)[0]
    manifold = workbench.manifolds[key]
    nodes = list(manifold._walk_nodes(stop_at_final=False))
    node = nodes[len(nodes) // 3]

    workbench.trim_stable_at(key, node.get_cdist("stable"))
    assert manifold.tail is node


def test_a_cdist_past_the_end_leaves_the_branch_untrimmed(small_tangle, caplog):
    workbench, fp = small_tangle
    key = _stable_branches(workbench, fp)[0]
    manifold = workbench.manifolds[key]
    tail = manifold.tail
    end = max(node.get_cdist("stable") for node in manifold._walk_nodes(stop_at_final=False))

    with caplog.at_level(logging.WARNING, logger="tanglepack.numerics.TangleWorkbench"):
        workbench.trim_stable_at(key, 2.0 * end)
    assert manifold.tail is tail
    assert_logged(caplog, logging.WARNING, "tanglepack.numerics.TangleWorkbench", count=1)


def test_trimming_the_manifolds_is_the_trim_at_the_outermost_crossing():
    """Two identical builds: one trimmed by ``trim_stable_manifolds``, the other by
    ``trim_stable_at`` at each branch's outermost registered crossing."""
    trimmed, by_hand = (build_k10(unstable_steps=7, through="intersected") for _ in range(2))
    workbench, fp = by_hand.workbench, by_hand.fixed_point
    outermost = {}
    for _, ix in workbench.intersection_registry:
        outermost[ix.manifold_b_key] = max(outermost.get(ix.manifold_b_key, 0.0), ix.stable_cdist)
    keys = [key for key in _stable_branches(workbench, fp) if outermost.get(key, 0.0) > 0.0]
    assert keys
    for key in keys:
        workbench.trim_stable_at(key, outermost[key])

    trimmed.workbench.trim_stable_manifolds(trimmed.fixed_point)
    for key in keys:
        other = (trimmed.fixed_point, *key[1:])
        expected = workbench.manifolds[key].tail.get_cdist("stable")
        assert trimmed.workbench.manifolds[other].tail.get_cdist("stable") == expected

"""One stable trim: ``TangleWorkbench.trim_stable_at`` (SIMP-10, WP11)."""

from __future__ import annotations

import logging

from cases import build_k10
from helpers.fakes import bare_fixed_point
from helpers.logs import assert_logged
from helpers.polylines import polyline_workbench
from tanglepack import BranchPoint


def _chain_rows(manifold) -> list[int]:
    """The rows from the root to the end of its chain (a periodic-point root has none)."""
    return manifold._rows(stop_at_final=False)[2].tolist()


def _first_row_at_or_past(manifold, cdist):
    """The first row past the root whose cdist is at or beyond ``cdist``."""
    store, rows = manifold.store, _chain_rows(manifold)
    if not isinstance(manifold.root, BranchPoint):
        rows = rows[1:]
    return next(row for row in rows if store.cdist[row] >= cdist)


def _stable_branches(workbench, fp):
    return [key for key in workbench.manifolds if key[0] is fp and key[1] == "stable"]


def test_the_trim_keeps_the_first_node_at_or_past_the_cdist(small_tangle):
    workbench, fp = small_tangle
    key = _stable_branches(workbench, fp)[0]
    manifold = workbench.manifolds[key]
    nodes = [manifold.root, *_chain_rows(manifold)]
    cdist_of = lambda node: node.get_cdist("stable") if node is manifold.root else manifold.store.cdist[node]  # noqa: E731
    lo, hi = cdist_of(nodes[len(nodes) // 3]), cdist_of(nodes[len(nodes) // 3 + 1])
    cdist = 0.5 * (lo + hi)  # strictly between two nodes
    expected = _first_row_at_or_past(manifold, cdist)

    assert workbench.trim_stable_at(key, cdist) is manifold
    assert manifold.ids[-1] == expected
    assert manifold.store.cdist[expected] >= cdist > lo


def test_a_cdist_on_a_node_keeps_that_node(small_tangle):
    """"At or beyond": a cdist exactly on a node makes that node the tail, not the next."""
    workbench, fp = small_tangle
    key = _stable_branches(workbench, fp)[0]
    manifold = workbench.manifolds[key]
    nodes = [manifold.root, *_chain_rows(manifold)]
    node = nodes[len(nodes) // 3]

    workbench.trim_stable_at(key, manifold.store.cdist[node])
    assert manifold.ids[-1] == node


def test_the_root_is_never_the_new_tail():
    """A branch rooted at a row, trimmed at or below the root's cdist, keeps the
    next row: the trim never shrinks a curve to its root alone."""
    fp = bare_fixed_point()
    key = (fp, "stable", 0, 0)
    workbench = polyline_workbench(fp, {key: [(0.0, 0.0), (0.0, 1.0), (0.0, 2.0)]})
    manifold = workbench.manifolds[key]
    _head, _store, rows = manifold._rows()

    for cdist in (0.0, -1.0):  # the root's cdist is 0.0
        workbench.trim_stable_at(key, cdist)
        assert len(manifold.ids) > 1
        assert manifold._rows()[2].tolist() == rows[:2].tolist()


def test_a_cdist_past_the_end_leaves_the_branch_untrimmed(small_tangle, caplog):
    workbench, fp = small_tangle
    key = _stable_branches(workbench, fp)[0]
    manifold = workbench.manifolds[key]
    tail = manifold.ids[-1]
    end = max(manifold.store.cdist[_chain_rows(manifold)])

    with caplog.at_level(logging.WARNING, logger="tanglepack.numerics.TangleWorkbench"):
        workbench.trim_stable_at(key, 2.0 * end)
    assert manifold.ids[-1] == tail
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
        expected = workbench.manifolds[key].tail.cdist
        assert trimmed.workbench.manifolds[other].tail.cdist == expected

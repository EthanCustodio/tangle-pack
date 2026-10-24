"""Tangle segment-index bookkeeping and point-insertion side effects.

Pins:
1. A query-only curve (a blast image curve, a bridge cut from it, a partial
   piece) holds no segment record of its own (N1, 2026-10-24): it owns no
   segment, and the only query-only ``_Segment`` the index keeps is one with a
   crossing pair still to read.
2. The edge dedup still decides which edges of a new image curve are queried:
   none an earlier image curve queried, none of a separator spliced into a
   bridge cut from one, and a curve that lies wholly on a branch is dropped on
   a fresh recompute. These held before N1 too; the identity probe's builds
   never reach them.
3. An indexed curve added twice records no ``None`` segment id: since N1
   ``_insert_segment`` returns None for an edge it only claims, and that
   return must never reach ``_manifold_segs``.
4. Inserting a point into a linked list must never overwrite the host
   point's ``stretch_param``. The inserted node inherits a missing value
   instead (a fresh crossing separator carries ``stretch_param=None``).
"""

from __future__ import annotations

import numpy as np

from cases import build_k10
from helpers.fakes import bare_fixed_point
from helpers.polylines import polyline_manifold, polyline_workbench
from tanglepack import BaseManifold, Bridge, Point, define_resonance_zone
from tanglepack.numerics.Tangle import Tangle


def _nodes(curve) -> list:
    return list(curve.get_point_array(return_nodes=True))


def _query_only_records_are_paired(workbench) -> None:
    """No image curve, nor a bridge cut from one, owns a segment; every
    query-only record left is half of a crossing pair."""
    tangle = workbench.Tangle
    curves = workbench._image_curves
    image_nodes = {id(node) for curve in curves for node in _nodes(curve)}
    owners = list(tangle._manifold_segs)
    assert not [curve for curve in curves if curve in owners]
    assert not [owner for owner in owners if id(owner.root) in image_nodes]
    assert any(id(bridge.root) in image_nodes for bridge in workbench.bridges if bridge.id is not None)

    paired = set().union(*tangle._intersecting_segments)
    query_only = [sid for sid, seg in tangle._seg_lookup.items() if any(seg.manifold is c for c in curves)]
    assert query_only and set(query_only) <= paired


def test_query_only_curves_hold_only_paired_records():
    case = build_k10(through="pips")
    session, fixed_point = case.session, case.fixed_points[0]
    zone = define_resonance_zone(session.workbench, session.strong_pip(fixed_point))
    session.blast_zone(zone, num_iterations=2, fixed_point=[fixed_point], diagnose=False)
    workbench = session.workbench
    _query_only_records_are_paired(workbench)

    workbench.compute_intersections(workbench.fixed_points, preserve_ids=True, diagnose=False)
    workbench.rebuild_bridges()
    _query_only_records_are_paired(workbench)


class _RecordingIndex:
    """An rtree stand-in that records the bounds it is queried with."""

    def __init__(self, rtree) -> None:
        self.rtree, self.queried = rtree, []

    def intersection(self, bounds):
        self.queried.append(tuple(bounds))
        return self.rtree.intersection(bounds)


def _bounds(p0, p1) -> tuple:
    (x0, y0), (x1, y1) = p0.get_point(), p1.get_point()
    return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))


def _stable_index(fp, stable: list) -> Tangle:
    """A Tangle whose stable rtree holds one hand-placed stable curve, recording queries."""
    tangle = Tangle()
    tangle.add_manifolds([polyline_manifold(stable, (fp, "stable", 0, 0))])
    tangle._rtrees["stable"] = _RecordingIndex(tangle._rtrees["stable"])
    return tangle


def test_an_image_over_an_earlier_image_queries_none_of_the_shared_edges():
    fp = bare_fixed_point(label="A")
    key = (fp, "unstable", 1, 0)
    tangle = _stable_index(fp, [(0.5, -1.0), (0.5, 1.0)])
    earlier = polyline_manifold([(x, 0.0) for x in range(4)], key)
    tangle.add_manifold(earlier, index_segments=False)

    # A later image starts on the earlier curve's nodes and runs on past its tail,
    # as merge_manifolds leaves two images that share straddle nodes.
    shared = _nodes(earlier)[1:]
    tail = shared[-1]
    for x in (4.0, 5.0):
        node = Point(x, 0.0, x)
        tail.insert_point_forward(node)
        tail = node
    later = BaseManifold(shared[0], "unstable", 1.0, fixed_point=fp, tail=tail,
                         branch_index=0, manifold_key=key)
    tangle._rtrees["stable"].queried.clear()
    tangle.add_manifold(later, index_segments=False)

    shared_bounds = {_bounds(a, b) for a, b in zip(shared, shared[1:])}
    queried = tangle._rtrees["stable"].queried
    assert len(shared_bounds) == 2 and len(queried) == 2
    assert not shared_bounds & set(queried)


def test_a_separator_in_a_query_only_bridge_is_not_queried_again():
    fp = bare_fixed_point(label="A")
    key = (fp, "unstable", 1, 0)
    # A stable V crossing ONE unstable segment twice: the cut splices a separator.
    tangle = _stable_index(fp, [(-0.5, 1.0), (-0.25, -1.0), (0.0, 1.0)])
    image = polyline_manifold([(-2.0, 0.0), (-1.0, 0.0), (1.0, 0.0), (2.0, 0.0)], key)
    tangle.add_manifold(image, index_segments=False)
    crossings = tangle.resolve_crossings()
    for i, crossing in enumerate(crossings):
        crossing.id = i
    (bridge,) = tangle.create_bridges(crossings, for_manifold=image)
    assert len(_nodes(bridge)) == 3, "the two crossings share a segment: one separator"

    over = Bridge(root=bridge.root, stability="unstable", stretch_param=1.0, fixed_point=fp,
                  tail=bridge.tail, branch_index=0, manifold_key=key)
    tangle._rtrees["stable"].queried.clear()
    assert not tangle.add_manifold(over, index_segments=False)
    assert not tangle._rtrees["stable"].queried


def test_a_curve_lying_on_a_branch_is_dropped_on_a_fresh_recompute():
    fp = bare_fixed_point(label="A")
    branch_key = (fp, "unstable", 0, 0)
    workbench = polyline_workbench(fp, {
        branch_key: [(x, 0.0) for x in np.arange(-3.0, 4.0)],
        (fp, "stable", 0, 0): [(0.5, -2.0), (0.5, 2.0)],
    })
    nodes = _nodes(workbench.manifolds[branch_key])

    def image(root, tail) -> Bridge:
        return Bridge(root=root, stability="unstable", stretch_param=1.0, fixed_point=fp,
                      tail=tail, branch_index=0, manifold_key=branch_key)

    # Growth merged this image into the branch: every edge of it is a branch edge.
    merged = image(nodes[2], nodes[5])
    own = polyline_manifold([(0.0, 1.0), (1.0, 1.5), (2.0, 1.0)], branch_key)
    standing = image(own.root, own.tail)
    workbench._image_curves[:] = [merged, standing]

    workbench.compute_intersections([fp], infer_iterates=False, diagnose=False)
    assert workbench._image_curves == [standing]


def test_add_manifold_twice_records_no_none_ids(henon_tangle_with_bridges):
    workbench, _ = henon_tangle_with_bridges
    tangle = workbench.Tangle

    # Re-adding purges and re-indexes; the second add of a branch claims the
    # edges its bridges still hold, where _insert_segment returns None.
    for manifold in workbench.manifolds.values():
        tangle.add_manifold(manifold)
        tangle.add_manifold(manifold)

    for manifold, seg_ids in tangle._manifold_segs.items():
        assert None not in seg_ids, f"None segment id recorded for {manifold}"
        assert seg_ids <= tangle._seg_lookup.keys()


def test_insert_point_keeps_host_stretch_param():
    host = Point(0.0, 0.0, cdist=0.0, stretch_param=2.5)
    separator = Point(1.0, 0.0, cdist=1.0)  # no stretch_param, like a crossing separator

    host.insert_point_forward(separator)

    assert host.stretch_param == 2.5, "host stretch_param was overwritten"
    assert separator.stretch_param == 2.5, "inserted node did not inherit"


def test_insert_point_does_not_clobber_existing_node_param():
    host = Point(0.0, 0.0, cdist=0.0, stretch_param=2.5)
    node = Point(1.0, 0.0, cdist=1.0, stretch_param=3.5)

    host.insert_point_backward(node)

    assert host.stretch_param == 2.5
    assert node.stretch_param == 3.5


def test_insert_iterate_touches_no_stretch_param():
    """Iterate-list inserts leave stretch_param alone on both sides.

    Callers construct iterate points with their stretch already set (or fill
    it explicitly afterwards), and the host may even be a BranchPoint, which
    has no scalar stretch_param at all.
    """
    host = Point(0.0, 0.0, cdist=1.0, stretch_param=2.5)
    image = Point(0.5, 0.5, cdist=2.0)

    host.insert_next_iterate(image)

    assert host.stretch_param == 2.5
    assert image.stretch_param is None

    preimage = Point(-0.5, -0.5, cdist=0.5)
    host.insert_prev_iterate(preimage)

    assert host.stretch_param == 2.5
    assert preimage.stretch_param is None

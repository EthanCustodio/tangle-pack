"""Tangle segment-index bookkeeping and point-insertion side effects.

Pins:
1. A query-only curve (a blast image curve, a bridge cut from it, a partial
   piece) holds no segment record of its own (N1, 2026-10-24): it owns no
   segment, and the only query-only ``_Segment`` the index keeps is one with a
   crossing pair still to read.
2. The edge dedup still decides which edges of a new image curve are queried:
   none an earlier image curve queried, none of a separator spliced into a
   bridge cut from one, and a curve that lies wholly on a branch is dropped on
   a fresh recompute; and a branch added over a query-only curve (the
   ``reset=False`` recompute) neither indexes nor queries the shared edges.
   These held before N1 too; the identity probe's builds never reach them.
3. An indexed curve added twice records no ``None`` segment id: since N1
   ``_insert_segment`` returns None for an edge it only claims, and that
   return must never reach ``_manifold_segs``.
4. The sid values decide float bits (N5b, 2026-10-26): a pair's two segments
   are resolved in the order ``tuple(frozenset(pair))`` gives, so a query-only
   add builds every pair with the NEW sid first, and the crossing's cdists are
   Python floats. A query-only curve takes one sid per new edge, in curve
   order, whether the edge pairs or not: a curve that crosses nothing still
   returns its sids, so a fresh recompute keeps it.
5. The edge leaving a periodic point is queried from the periodic point's
   coordinates (N6, 2026-10-26): the point has no row, and its stand-in index
   must never read a real row's coordinates. N8 keeps this spelling: the
   periodic point is a view of its fixed point, still on no row. The decoy
   row 0 gives the pin its teeth.

Curves that share points are views of one curve between different ends, as
the images a blast merged and the branch growth merged them into are.
"""

from __future__ import annotations

from itertools import count

import numpy as np

from cases import build_k10
from helpers.fakes import bare_fixed_point
from helpers.polylines import polyline_manifold, polyline_workbench
from tanglepack import BaseManifold, Bridge, BranchPoint, define_resonance_zone
from tanglepack.numerics import ManifoldMachine
from tanglepack.numerics.PointStore import PointStore
from tanglepack.numerics.Tangle import Tangle


def _rows(curve) -> list[int]:
    return curve.ids.tolist()


def _root_in(curve, rows: set, heads: set) -> bool:
    """Whether the curve's root is one of ``rows`` or one of the periodic points ``heads``."""
    root = curve.root
    return id(root) in heads if isinstance(root, BranchPoint) else int(curve.ids[0]) in rows


def _query_only_records_are_paired(workbench) -> None:
    """No image curve, nor a bridge cut from one, owns a segment; every
    query-only record left is half of a crossing pair."""
    tangle = workbench.Tangle
    curves = workbench._image_curves
    image_rows = {row for curve in curves for row in _rows(curve)}
    image_heads = {id(curve.root) for curve in curves if isinstance(curve.root, BranchPoint)}
    owners = list(tangle._manifold_segs)
    assert not [curve for curve in curves if curve in owners]
    assert not [owner for owner in owners if _root_in(owner, image_rows, image_heads)]
    assert any(_root_in(bridge, image_rows, image_heads) for bridge in workbench.bridges if bridge.id is not None)

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

    def intersection_v(self, mins, maxs):
        self.queried.extend(tuple(box) for box in np.hstack((mins, maxs)).tolist())
        return self.rtree.intersection_v(mins, maxs)


def _bounds(p0, p1) -> tuple:
    (x0, y0), (x1, y1) = p0, p1
    return (min(x0, x1), min(y0, y1), max(x0, x1), max(y0, y1))


def _stable_index(fp, stable: list) -> Tangle:
    """A Tangle whose stable rtree holds one hand-placed stable curve, recording queries."""
    tangle = Tangle()
    tangle.add_manifolds([polyline_manifold(stable, (fp, "stable", 0, 0))])
    tangle._rtrees["stable"] = _RecordingIndex(tangle._rtrees["stable"])
    return tangle


def _view(curve: BaseManifold, root: int, tail: int, key: tuple) -> BaseManifold:
    """The curve over ``curve``'s rows from index ``root`` to index ``tail``."""
    rows = curve.ids
    return BaseManifold(int(rows[root]), "unstable", 1.0, fixed_point=key[0], tail=int(rows[tail]),
                        branch_index=0, manifold_key=key, store=curve.store)


def test_an_image_over_an_earlier_image_queries_none_of_the_shared_edges():
    fp = bare_fixed_point(label="A")
    key = (fp, "unstable", 1, 0)
    tangle = _stable_index(fp, [(0.5, -1.0), (0.5, 1.0)])
    whole = polyline_manifold([(x, 0.0) for x in range(6)], key)
    xy = whole.store.xy[whole.ids]
    earlier = _view(whole, 0, 3, key)
    tangle.add_manifold(earlier, index_segments=False)

    # A later image starts on the earlier curve's nodes and runs on past its tail,
    # as merge_manifolds leaves two images that share straddle nodes.
    later = _view(whole, 1, 5, key)
    tangle._rtrees["stable"].queried.clear()
    tangle.add_manifold(later, index_segments=False)

    shared = xy[1:4]
    shared_bounds = {_bounds(a, b) for a, b in zip(shared, shared[1:])}
    queried = tangle._rtrees["stable"].queried
    assert len(shared_bounds) == 2 and len(queried) == 2
    assert not shared_bounds & set(queried)


def test_a_branch_over_a_query_only_curve_neither_indexes_nor_queries_its_edges():
    # The reset=False recompute adds a branch to an index that already holds
    # image curves; as before N1, an edge one of them queried is left alone.
    fp = bare_fixed_point(label="A")
    key = (fp, "unstable", 0, 0)
    tangle = _stable_index(fp, [(0.5, -1.0), (0.5, 1.0)])
    whole = polyline_manifold([(x, 0.0) for x in range(5)], key)
    xy = whole.store.xy[whole.ids]
    tangle.add_manifold(_view(whole, 0, 3, key), index_segments=False)

    branch = _view(whole, 0, 4, key)
    tangle._rtrees["stable"].queried.clear()
    tangle.add_manifold(branch)

    assert len(tangle._rtrees["unstable"]) == 1
    assert tangle._rtrees["stable"].queried == [_bounds(xy[3], xy[4])]


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
    assert len(_rows(bridge)) == 3, "the two crossings share a segment: one separator"

    over = Bridge(root=bridge.root, stability="unstable", stretch_param=1.0, fixed_point=fp,
                  tail=bridge.tail, branch_index=0, manifold_key=key, store=bridge.store)
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
    rows = _rows(workbench.manifolds[branch_key])
    store = workbench.manifolds[branch_key].store

    def image(root, tail) -> Bridge:
        return Bridge(root=root, stability="unstable", stretch_param=1.0, fixed_point=fp,
                      tail=tail, branch_index=0, manifold_key=branch_key, store=store)

    # Growth merged this image into the branch: every edge of it is a branch edge.
    merged = image(rows[2], rows[5])
    own = polyline_manifold([(0.0, 1.0), (1.0, 1.5), (2.0, 1.0)], branch_key)
    standing = image(own.root, own.tail)
    workbench._image_curves[:] = [merged, standing]

    workbench.compute_intersections([fp], infer_iterates=False, diagnose=False)
    assert workbench._image_curves == [standing]


def test_a_crossing_pair_keeps_the_new_sid_first(monkeypatch):
    fp = bare_fixed_point(label="A")
    start = 10**6
    monkeypatch.setattr(Tangle, "_ids", count(start))
    tangle = Tangle()
    tangle.add_manifolds([polyline_manifold([(0.5, -1.0), (0.5, 1.0)], (fp, "stable", 0, 0))])
    (cand_sid,) = tangle._seg_lookup
    # The new sid lands in the candidate's hash slot, where the frozenset's
    # element order depends on which element went in first.
    new_sid = cand_sid + 8
    assert tuple(frozenset((new_sid, cand_sid))) != tuple(frozenset((cand_sid, new_sid)))
    monkeypatch.setattr(Tangle, "_ids", count(new_sid))

    image = polyline_manifold([(0.0, 0.0), (1.0, 0.0)], (fp, "unstable", 1, 0))
    assert tangle.add_manifold(image, index_segments=False) == {new_sid}
    (pair,) = tangle._intersecting_segments
    assert tuple(pair) == (new_sid, cand_sid)

    (crossing,) = tangle.resolve_crossings()
    assert type(crossing.unstable_cdist) is float and type(crossing.stable_cdist) is float


def test_a_curve_crossing_nothing_takes_one_sid_per_new_edge(monkeypatch):
    fp = bare_fixed_point(label="A")
    branch_key = (fp, "unstable", 0, 0)
    workbench = polyline_workbench(fp, {
        branch_key: [(x, 0.0) for x in np.arange(-3.0, 4.0)],
        (fp, "stable", 0, 0): [(0.5, -2.0), (0.5, 2.0)],
    })
    workbench.compute_intersections([fp], infer_iterates=False, diagnose=False)
    image = polyline_manifold([(x, 5.0) for x in range(5)], branch_key)

    start = next(Tangle._ids)
    monkeypatch.setattr(Tangle, "_ids", count(start))
    sids = workbench.Tangle.add_manifold(image, index_segments=False)
    assert sids == set(range(start, start + 4))
    assert next(Tangle._ids) == start + 4
    assert not any(pair & sids for pair in workbench.Tangle._intersecting_segments)

    standing = Bridge(root=image.root, stability="unstable", stretch_param=1.0, fixed_point=fp,
                      tail=image.tail, branch_index=0, manifold_key=branch_key, store=image.store)
    workbench._image_curves[:] = [standing]
    workbench.compute_intersections([fp], infer_iterates=False, diagnose=False)
    assert workbench._image_curves == [standing]


def test_a_query_only_curve_queries_the_edge_leaving_its_periodic_point():
    fp = bare_fixed_point(label="A")
    root = fp.branch_points[0]
    store = PointStore()
    store.append_one(5.0, 5.0, None)  # row 0, a decoy far from every curve

    def curve(points, key, root=None) -> BaseManifold:
        rows = store.append(points, [float(i + 1) for i in range(len(points))])
        ManifoldMachine._manifold_from_rows(store, rows, key[1], 1.0, fp)
        if root is not None:
            root.attach(store.chain_of(rows[0]), key[1], key[3])
        return BaseManifold(int(rows[0]) if root is None else root, key[1], 1.0, fixed_point=fp,
                            tail=int(rows[-1]), branch_index=key[3], manifold_key=key, store=store)

    tangle = Tangle()
    tangle.add_manifold(curve([(0.5, -1.0), (0.5, 1.0)], (fp, "stable", 0, 0)))
    image = curve([(1.0, 0.0), (2.0, 0.0)], (fp, "unstable", 0, 0), root=root)
    (crossing,) = tangle.populate_intersections(tangle.add_manifold(image, index_segments=False))
    assert np.allclose(crossing.coords, (0.5, 0.0))
    assert crossing.unstable_segment[0] is root


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

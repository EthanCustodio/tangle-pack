"""The same-stability polyline diagnostic (WP1, 2026-10-07).

Two unstable curves never cross, nor do two stable ones, and no curve crosses
itself (CLAUDE.md fundamental invariant), so a crossing of two same-stability
POLYLINES is a polygon artifact of under-resolution.
:func:`tanglepack.numerics.diagnostics.same_stability_crossings` finds them --
self-crossings included, which the Tangle index never tests -- and the
workbench reports them: WARNING inside the tangle region, INFO outside. The
test-tier law ``helpers.laws.check_no_same_stability_crossing`` stays the
independent oracle; here the two are only compared on synthetic polylines.
"""

from __future__ import annotations

import itertools
import logging
import time

import numpy as np
import pytest

from cases import build_k10_inversion
from helpers import laws
from helpers.fakes import bare_fixed_point
from helpers.logs import assert_logged
from helpers.polylines import TEST_STORE, periodic_point, store_rows, polyline_manifold, polyline_workbench
from tanglepack import BaseManifold, Bridge
from tanglepack.numerics import ManifoldMachine
from tanglepack.numerics.PointStore import PointStore
from tanglepack.numerics.diagnostics import (
    count_same_stability_crossings,
    same_stability_crossings,
)

WORKBENCH = "tanglepack.numerics.TangleWorkbench"

#: Segment 0 and segment 2 cross; segment 1 is the fold between them.
Z_KINK = np.array([[0.0, 0.0], [2.0, 0.0], [0.5, 0.5], [1.5, -0.5], [3.0, -0.5]])


@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_a_z_kink_is_reported_once(stability: str):
    fp = bare_fixed_point(label="A")
    workbench = polyline_workbench(fp, {(fp, stability, 0, 0): Z_KINK})
    (hit,) = same_stability_crossings(workbench, region="all")
    assert hit.stability == stability
    assert hit.curve_a == hit.curve_b
    assert {hit.seg_a, hit.seg_b} == {0, 2}
    assert 0.0 < hit.suggested_area_cutoff < hit.local_area
    assert 0.0 < hit.sin <= 1.0


def test_neighbouring_segments_and_a_shared_root_are_not_crossings():
    fp = bare_fixed_point(label="A")
    # A hairpin: every segment meets the next at a sharp turn, none crosses.
    hairpin = np.array([[0.0, 0.0], [2.0, 0.0], [0.0, 0.01], [2.0, 0.02], [0.0, 0.03]])
    workbench = polyline_workbench(fp, {(fp, "stable", 0, 0): hairpin})
    root = periodic_point(5.0, 5.0)
    for branch, direction in ((0, 1.0), (1, -1.0)):
        key = (fp, "unstable", 0, branch)
        points = [(5.0 + direction * step, 5.0 + 0.1 * step) for step in (1, 2, 3)]
        workbench.register_manifold(key, polyline_manifold(points, key, root=root))
    assert same_stability_crossings(workbench, region="all") == []


def test_a_bridge_cut_from_a_manifold_is_that_manifold():
    fp = bare_fixed_point(label="A")
    unstable = np.column_stack([np.linspace(-3.0, 3.0, 13), np.zeros(13)])
    stable = np.array([[-2.0, -1.0], [-2.0, 1.0], [0.25, 1.0], [0.25, -1.0], [2.0, -1.0], [2.0, 1.0]])
    workbench = polyline_workbench(
        fp, {(fp, "unstable", 0, 0): unstable, (fp, "stable", 0, 0): stable}
    )
    workbench.compute_intersections([fp], infer_iterates=False)
    assert workbench.create_bridges(fp)
    assert same_stability_crossings(workbench, region="all") == []


def test_a_bridge_off_every_manifold_leaving_a_manifold_s_periodic_point_is_that_manifold():
    """A bridge leaving the periodic point a manifold leaves from, on the other
    slot, is cut from that manifold even when its rows lie on no manifold: their
    crossing is neither listed nor counted (N6's mC / mH, pinned in N8)."""
    fp = bare_fixed_point(label="A")
    workbench = polyline_workbench(fp, {})
    periodic = periodic_point()
    key = (fp, "unstable", 0, 0)
    workbench.register_manifold(key, polyline_manifold([(1.0, 0.0), (2.0, 0.0)], key, root=periodic))
    rows = store_rows([(1.5, 1.0), (1.5, -1.0)], [1.0, 2.0])
    periodic.attach(TEST_STORE.new_chain(rows, 1.0), "unstable", 1)
    bridge = Bridge(root=periodic, stability="unstable", stretch_param=1.0, fixed_point=fp,
                    tail=int(rows[-1]), branch_index=1, manifold_key=(fp, "unstable", 0, 1),
                    store=TEST_STORE)
    workbench._register_bridge(bridge)

    assert same_stability_crossings(workbench, region="all") == []
    count = count_same_stability_crossings(workbench, [bridge])
    assert (count.in_region, count.outside) == (0, 0)


def test_two_curves_of_one_stability_cross_for_both_stabilities_never_across():
    fp = bare_fixed_point(label="A")
    rising = np.array([[-1.0, -1.0], [0.1, 0.1], [1.0, 1.0]])
    falling = np.array([[-1.0, 1.0], [0.1, -0.1], [1.0, -1.0]])
    workbench = polyline_workbench(
        fp,
        {
            (fp, "unstable", 0, 0): rising,
            (fp, "unstable", 0, 1): falling,
            (fp, "stable", 0, 0): rising + 5.0,
            (fp, "stable", 0, 1): falling + 5.0,
        },
    )
    # The unstable pair also crosses nothing stable, and vice versa.
    hits = same_stability_crossings(workbench, region="all")
    assert sorted(hit.stability for hit in hits) == ["stable", "unstable"]
    assert all(hit.curve_a != hit.curve_b for hit in hits)


def test_an_out_of_region_artifact_is_info_and_an_in_region_one_a_warning(caplog):
    fp = bare_fixed_point(label="A")
    line = np.column_stack([np.linspace(-1.0, 1.0, 5), np.zeros(5)])
    cross = np.column_stack([np.full(5, 0.1), np.linspace(-1.0, 1.0, 5)])
    far_kink = Z_KINK + 100.0
    workbench = polyline_workbench(
        fp,
        {
            (fp, "unstable", 0, 0): line,
            (fp, "stable", 0, 0): cross,
            (fp, "stable", 0, 1): far_kink,
        },
    )
    caplog.set_level(logging.DEBUG, logger=WORKBENCH)
    workbench.compute_intersections([fp], infer_iterates=False)
    (hit,) = workbench.check_resolution().same_stability_crossings
    assert not hit.in_tangle_region
    assert_logged(caplog, logging.INFO, WORKBENCH)
    assert_logged(caplog, logging.WARNING, WORKBENCH, count=0)

    near_kink = 0.2 * Z_KINK + np.array([-0.7, 0.05])  # crosses itself at (-0.5, 0.05)
    workbench.register_manifold(
        (fp, "stable", 0, 1), polyline_manifold(near_kink, (fp, "stable", 0, 1))
    )
    caplog.clear()
    workbench.compute_intersections([fp], infer_iterates=False)
    (hit,) = workbench.check_resolution().same_stability_crossings
    assert hit.in_tangle_region
    assert not workbench.check_resolution().clean
    assert_logged(caplog, logging.WARNING, WORKBENCH, count=1)

    caplog.clear()
    workbench.compute_intersections([fp], infer_iterates=False)
    assert_logged(caplog, logging.WARNING, WORKBENCH, count=0)


def test_agrees_with_the_test_tier_oracle_on_random_polylines():
    rng = np.random.default_rng(20261007)
    fp = bare_fixed_point(label="A")
    polylines = [np.cumsum(rng.normal(size=(40, 2)), axis=0) for _ in range(4)]
    keys = [(fp, "unstable", 0, index) for index in range(len(polylines))]
    workbench = polyline_workbench(fp, dict(zip(keys, polylines)))
    names = {f"A unstable 0.{index}": index for index in range(len(polylines))}

    found = {
        frozenset({(names[hit.curve_a], hit.seg_a), (names[hit.curve_b], hit.seg_b)})
        for hit in same_stability_crossings(workbench, region="all")
    }
    everything = np.vstack(polylines)
    lo, hi = everything.min(axis=0) - 1.0, everything.max(axis=0) + 1.0
    expected = {
        frozenset({(int(a), int(i)), (int(b), int(j))})
        for a, i, b, j in laws._transversal_hits(polylines, lo, hi)
    }
    assert expected and found == expected


def _image_chain(fp, points, cuts, store=None):
    """Partial image bridges along one curve through ``points``, cut at the
    node windows ``cuts``; overlapping windows share edges, as a blast's do."""
    store = TEST_STORE if store is None else store
    rows = store_rows(points, store=store)
    ManifoldMachine._manifold_from_rows(store, rows, "unstable", 1.0, fp)
    return [
        Bridge(
            root=int(rows[lo]), stability="unstable", stretch_param=1.0, fixed_point=fp,
            tail=int(rows[hi]), branch_index=0, manifold_key=(fp, "unstable", 0, 0), store=store,
        )
        for lo, hi in cuts
    ]


def test_two_periodic_points_are_two_nodes():
    """The first segments of curves leaving two different periodic points share
    no node, so their crossing is listed and counted (N6, 2026-10-26: a periodic
    point has no row, and each one is keyed apart)."""
    fp = bare_fixed_point(label="A")
    workbench = polyline_workbench(fp, {})
    here, there = periodic_point(0.0, 0.0), periodic_point(0.0, 2.0)
    key, image_key = (fp, "unstable", 0, 0), (fp, "unstable", 1, 0)
    workbench.register_manifold(key, polyline_manifold([(2.0, 2.0), (3.0, 2.0)], key, root=here))
    image = polyline_manifold([(2.0, 0.0), (3.0, 0.0)], image_key, root=there)
    bridge = Bridge(root=there, stability="unstable", stretch_param=1.0, fixed_point=fp,
                    tail=image.tail, branch_index=0, manifold_key=image_key)
    workbench._register_bridge(bridge)

    (hit,) = same_stability_crossings(workbench, region="all")
    assert (hit.seg_a, hit.seg_b) == (0, 0) and np.allclose(hit.point, (1.0, 1.0))
    count = count_same_stability_crossings(workbench, [bridge])
    assert (count.in_region, count.outside) == (0, 1)


def test_equal_rows_of_two_stores_are_two_nodes():
    """Two curves held in two stores, crossing on segments with EQUAL rows in
    their stores, share no node: the crossing is listed and counted (N6,
    2026-10-26: a row is keyed with its store)."""
    fp = bare_fixed_point(label="A")
    key = (fp, "unstable", 0, 0)
    store = PointStore()
    rows = store_rows([(0.0, 0.0), (2.0, 2.0), (3.0, 2.0)], store=store)
    ManifoldMachine._manifold_from_rows(store, rows, "unstable", 1.0, fp)
    workbench = polyline_workbench(fp, {})
    workbench.register_manifold(key, BaseManifold(int(rows[0]), "unstable", 1.0, fixed_point=fp, tail=int(rows[-1]),
                                                  branch_index=0, manifold_key=key, store=store))
    (bridge,) = _image_chain(fp, [(0.0, 2.0), (2.0, 0.0), (3.0, 0.0)], [(0, 2)], store=PointStore())
    workbench._register_bridge(bridge)
    assert (rows[0], rows[1]) == (bridge.ids[0], bridge.ids[0] + 1) == (0, 1)

    (hit,) = same_stability_crossings(workbench, region="all")
    assert (hit.seg_a, hit.seg_b) == (0, 0) and np.allclose(hit.point, (1.0, 1.0))
    count = count_same_stability_crossings(workbench, [bridge])
    assert (count.in_region, count.outside) == (0, 1)


@pytest.mark.parametrize("seed", [1, 2, 3, 4, 5])
def test_the_count_is_the_full_listing_of_the_new_curves_on_random_curves(seed: int):
    """C5 on random curves: the count builds pairs for the new segments only (N2,
    2026-10-24) and still finds exactly the full sweep's hits naming a new curve,
    in the same orientation -- self-crossings of a new curve, new against new,
    edges a new bridge shares with an old one and tied box minima (the odd seeds
    round every vertex) included."""
    rng = np.random.default_rng(seed)
    grid = np.round if seed % 2 else (lambda points: points)
    walk = lambda count: grid(np.cumsum(rng.normal(size=(count, 2)), axis=0))  # noqa: E731
    fp = bare_fixed_point(label="A")
    workbench = polyline_workbench(
        fp,
        {
            (fp, "unstable", 0, 0): walk(30),
            (fp, "unstable", 0, 1): walk(30),
            (fp, "stable", 0, 0): [(0.5, -40.0), (0.5, 40.0)],
        },
    )
    workbench.compute_intersections([fp], infer_iterates=False, diagnose=False)
    cuts = [(0, 12), (11, 24), (23, 40), (39, 59)]
    first = _image_chain(fp, walk(60), cuts)
    kinked = np.vstack([walk(20), 3.0 * Z_KINK])  # its last bridge crosses itself
    second = _image_chain(fp, kinked + 0.25, [(0, 10), (9, len(kinked) - 1)])
    bridges = first + second
    for bridge in bridges:
        workbench._register_bridge(bridge)
    new = [bridge for bridge in bridges if rng.random() < 0.5 or bridge is second[-1]]

    names = {f"partial image bridge {id(bridge):#x}" for bridge in new}
    involved = [
        hit
        for hit in same_stability_crossings(workbench, [fp])
        if hit.stability == "unstable" and {hit.curve_a, hit.curve_b} & names
    ]
    inside = [hit for hit in involved if hit.in_tangle_region]
    count = count_same_stability_crossings(workbench, new)
    assert any(hit.curve_a == hit.curve_b for hit in involved)
    assert count.in_region == len(inside)
    assert count.outside == len(involved) - len(inside)
    assert count.curve_pairs == tuple(sorted({f"{hit.curve_a} x {hit.curve_b}" for hit in inside}))
    assert count.suggested_area_cutoff == (
        min(hit.suggested_area_cutoff for hit in inside) if inside else None
    )


def _sweep_oracle(curves, box):
    """The sweep's hits by plain loops, in its order. ``curves`` are ``(name,
    nodes, points)`` in listing order; returns ``(name_a, seg_a, name_b, seg_b,
    in_box)`` with the earlier segment along the sweep axis first."""
    segments, seen = [], set()
    for name, nodes, points in curves:
        for seg, (ends, a, b) in enumerate(zip(zip(nodes, nodes[1:]), points, points[1:])):
            if frozenset(ends) not in seen:  # one edge reached from two curves is one segment
                seen.add(frozenset(ends))
                segments.append((name, seg, set(ends), np.asarray(a, float), np.asarray(b, float)))
    lo = [np.minimum(a, b) for *_, a, b in segments]
    hi = [np.maximum(a, b) for *_, a, b in segments]
    overlap = lambda s, t, k: lo[s][k] <= hi[t][k] and lo[t][k] <= hi[s][k]  # noqa: E731
    pairs = list(itertools.combinations(range(len(segments)), 2))
    axis = min((0, 1), key=lambda k: sum(overlap(s, t, k) for s, t in pairs))  # fewer overlaps, 0 on a tie
    sweep = sorted(range(len(segments)), key=lambda s: lo[s][axis])
    hits = []
    for s, t in itertools.combinations(sweep, 2):
        (name_s, seg_s, ends_s, a, b), (name_t, seg_t, ends_t, c, d) = segments[s], segments[t]
        if not (overlap(s, t, 0) and overlap(s, t, 1)) or ends_s & ends_t:
            continue
        d1, d2, r = b - a, d - c, c - a
        den = d1[0] * d2[1] - d1[1] * d2[0]
        if abs(den) <= 1e-12 * np.hypot(*d1) * np.hypot(*d2):  # the module's parallel threshold
            continue
        along_s, along_t = (r[0] * d2[1] - r[1] * d2[0]) / den, (r[0] * d1[1] - r[1] * d1[0]) / den
        if 1e-9 < along_s < 1 - 1e-9 and 1e-9 < along_t < 1 - 1e-9:  # interior, as the module's margin
            point = a + along_s * d1
            hits.append((name_s, seg_s, name_t, seg_t, bool(np.all((box[0] <= point) & (point <= box[1])))))
    return hits


@pytest.mark.parametrize("seed", [1, 2, 3, 4])
def test_the_listing_and_the_count_are_the_sweep_s_in_its_order_and_orientation(seed: int):
    """The full listing is the sweep's hits in sweep order, each pair oriented
    earlier-along-the-sweep first (not curve order), and the count names its pairs
    in that orientation (D1, 2026-10-28: pinned before the strip sweep replaced
    the global one). Shared periodic root, overlapping image bridges, a stable
    field and tied box minima (odd seeds round every vertex) included."""
    rng = np.random.default_rng(seed)
    grid = np.round if seed % 2 else (lambda points: points)
    walk = lambda count: grid(np.cumsum(rng.normal(size=(count, 2)), axis=0))  # noqa: E731
    fp = bare_fixed_point(label="A", coordinates=[(0.0, 0.0)])
    root = fp.branch_points[0]
    workbench = polyline_workbench(fp, {})
    curves = {"unstable": [], "stable": []}
    for key, points, start in [
        ((fp, "unstable", 0, 0), walk(30), root),
        ((fp, "unstable", 0, 1), walk(30), root),
        ((fp, "unstable", 1, 0), walk(30), None),
        ((fp, "stable", 0, 0), walk(30), None),
        ((fp, "stable", 0, 1), [(0.5, -40.0), (0.5, 40.0)], None),
    ]:
        manifold = polyline_manifold(points, key, root=start)
        workbench.register_manifold(key, manifold)
        head = [root] if start else []
        nodes, xy = head + list(manifold.ids), [p.get_point() for p in head] + list(points)
        curves[key[1]].append((f"A {key[1]} {key[2]}.{key[3]}", nodes, xy))
    workbench.compute_intersections([fp], infer_iterates=False, diagnose=False)
    chain = walk(50)
    bridges = _image_chain(fp, chain, [(0, 14), (13, 30), (29, 49)])
    for bridge in bridges:
        workbench._register_bridge(bridge)
        rows = list(bridge.ids)
        start = rows[0] - bridges[0].ids[0]
        curves["unstable"].append((f"partial image bridge {id(bridge):#x}", rows, chain[start : start + len(rows)]))
    everywhere = (np.full(2, -np.inf), np.full(2, np.inf))

    listed = [(h.curve_a, h.seg_a, h.curve_b, h.seg_b) for h in same_stability_crossings(workbench, region="all")]
    expected = [hit[:4] for stability in ("unstable", "stable") for hit in _sweep_oracle(curves[stability], everywhere)]
    curve_order = {name: k for k, (name, *_) in enumerate(curves["unstable"] + curves["stable"])}
    assert any(curve_order[a] > curve_order[b] for a, _, b, _ in expected) and listed == expected

    crossings = np.array([ix.coords for _id, ix in workbench.intersection_registry])
    lo, hi = crossings.min(axis=0), crossings.max(axis=0)
    box = (lo - 0.1 * (hi - lo).max(), hi + 0.1 * (hi - lo).max())
    new = [bridges[0], bridges[1 + seed % 2]]
    names = {f"partial image bridge {id(bridge):#x}" for bridge in new}
    involved = [hit for hit in _sweep_oracle(curves["unstable"], box) if {hit[0], hit[2]} & names]
    count = count_same_stability_crossings(workbench, new)
    assert (count.in_region, count.outside) == (sum(h[4] for h in involved), sum(not h[4] for h in involved))
    assert count.curve_pairs == tuple(sorted({f"{a} x {b}" for a, _, b, _, inside in involved if inside}))


@pytest.mark.perf
def test_the_diagnostic_costs_at_most_half_a_compute_intersections():
    """On the largest law build (the inversion case's ~5e4 points, its escaping
    arm crossing itself thousands of times outside the tangle box) the pass
    stays vectorised: at most half the time of the crossing computation."""
    case = build_k10_inversion(through="grown")
    workbench, fixed_points = case.workbench, case.fixed_points

    def best_of_three(call) -> float:
        times = []
        for _ in range(3):
            start = time.perf_counter()
            call()
            times.append(time.perf_counter() - start)
        return min(times)

    compute = best_of_three(
        lambda: workbench.compute_intersections(
            fixed_points, infer_iterates=False, diagnose=False
        )
    )
    diagnostic = best_of_three(lambda: same_stability_crossings(workbench))
    assert diagnostic <= 0.5 * compute, (diagnostic, compute)

"""The crossing resolution ratio and its report (WP1, 2026-10-07).

Every detected crossing carries ``resolution_ratio = min(l_u, l_s) * |sin| /
max(h_u, h_s)``: how far the two segments separate across the shorter one,
over how far the curves they stand for may stray from them. The ratio only
flags -- nothing is discarded -- and the workbench reports new flags in one
WARNING on its own logger. Rules only: no case's ratio or count is pinned.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from cases import build_k28
from helpers.fakes import bare_fixed_point
from helpers.logs import assert_logged
from helpers.polylines import polyline_manifold, polyline_workbench
from tanglepack import Intersection, IntersectionRegistry
from tanglepack.numerics.Tangle import Tangle

WORKBENCH = "tanglepack.numerics.TangleWorkbench"
THRESHOLD = 3.0


def _ratios(unstable: np.ndarray, stable: np.ndarray) -> list[float]:
    """The resolution ratios of every crossing of two hand-placed curves."""
    fp = bare_fixed_point(label="A")
    tangle = Tangle()
    tangle.add_manifolds(
        [
            polyline_manifold(unstable, (fp, "unstable", 0, 0)),
            polyline_manifold(stable, (fp, "stable", 0, 0)),
        ]
    )
    return sorted(ix.resolution_ratio for ix in tangle.resolve_crossings())


def _line(spacing: float = 1.0, count: int = 6) -> np.ndarray:
    """A straight horizontal polyline whose middle segment straddles the origin."""
    xs = (np.arange(count) - count / 2 + 0.5) * spacing
    return np.column_stack([xs, np.zeros(count)])


def _arc(angle: float) -> np.ndarray:
    """A gently curved polyline through the origin, rotated by ``angle``."""
    s = np.array([-1.5, -0.5, 0.5, 1.5])
    local = np.column_stack([s, 0.05 * s**2])
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])
    return local @ rotation.T


def _shallow() -> np.ndarray:
    """A coarse parabola grazing the x-axis: one chord crosses it almost parallel."""
    xs = np.arange(-3, 4) + 0.45
    return np.column_stack([xs, 0.05 * xs**2 - 0.0125])


# --------------------------------------------------------------------------- #
# The ratio
# --------------------------------------------------------------------------- #
def test_ratio_is_similarity_invariant():
    xs = np.arange(-2.05, 2.0, 0.5)
    unstable = np.column_stack([xs, 0.1 * xs**2])
    stable = _arc(1.0) + np.array([0.1, 0.2])
    reference = _ratios(unstable, stable)
    assert reference

    angle, scale, shift = 0.7, 1e-4, np.array([3.0, -2.0])
    rotation = np.array([[np.cos(angle), -np.sin(angle)], [np.sin(angle), np.cos(angle)]])

    def moved(points: np.ndarray) -> np.ndarray:
        return scale * points @ rotation.T + shift

    assert _ratios(moved(unstable), moved(stable)) == pytest.approx(reference, rel=1e-6)


def test_ratio_grows_with_the_crossing_angle():
    ratios = [_ratios(_line(), _arc(angle)) for angle in (0.3, 0.6, 0.9, 1.2, 1.5)]
    assert all(len(r) == 1 for r in ratios)
    flat = [r[0] for r in ratios]
    assert flat == sorted(flat) and len(set(flat)) == len(flat)


def test_ratio_grows_with_segment_length_at_fixed_deviation():
    # The unstable line is straight (no deviation of its own) and shorter than the
    # stable segment, so lengthening it raises min(l) and nothing else.
    ratios = [_ratios(_line(spacing), _arc(np.pi / 2)) for spacing in (0.1, 0.3, 0.9)]
    assert all(len(r) == 1 for r in ratios)
    flat = [r[0] for r in ratios]
    assert flat == sorted(flat) and len(set(flat)) == len(flat)


def test_shallow_grazing_pair_scores_below_a_transversal_one_and_the_threshold():
    shallow = _ratios(_line(1.0, 8), _shallow())
    transversal = _ratios(_line(1.0, 8), _arc(np.pi / 2))
    assert shallow and transversal
    assert min(shallow) < THRESHOLD < min(transversal)


def test_straight_segments_give_a_finite_ratio():
    (ratio,) = _ratios(_line(), np.array([[0.1, -2.0], [0.1, -0.7], [0.1, 0.6], [0.1, 1.9]]))
    assert np.isfinite(ratio) and ratio > THRESHOLD


# --------------------------------------------------------------------------- #
# The registry and the workbench
# --------------------------------------------------------------------------- #
def test_unresolved_crossings_lists_the_ids_below_the_threshold_in_unstable_order():
    registry = IntersectionRegistry()
    ratios = {0.5: 4.0, 0.1: 0.2, 0.9: 2.9, 0.3: 3.0, 0.7: None}
    ids = {}
    for seg, (cdist, ratio) in enumerate(ratios.items()):
        ids[cdist] = registry.add(
            Intersection.from_segments(
                coords=(cdist, 0.0),
                unstable_cdist=cdist,
                stable_cdist=cdist,
                seg1_id=2 * seg,
                seg2_id=2 * seg + 1,
                resolution_ratio=ratio,
            )
        )
    assert registry.unresolved_crossings(THRESHOLD) == [ids[0.1], ids[0.9]]


def _crossing_workbench(stable: np.ndarray, threshold=THRESHOLD):
    fp = bare_fixed_point(label="A")
    workbench = polyline_workbench(
        fp,
        {(fp, "unstable", 0, 0): _line(1.0, 8), (fp, "stable", 0, 0): stable},
        resolution_ratio_threshold=threshold,
    )
    return workbench, fp


def test_unresolved_crossings_raise_one_warning_on_the_workbench_logger(caplog):
    workbench, fp = _crossing_workbench(_shallow())
    caplog.set_level(logging.DEBUG, logger=WORKBENCH)
    workbench.compute_intersections([fp], infer_iterates=False)
    assert_logged(caplog, logging.WARNING, WORKBENCH, count=1)
    assert workbench.pending_resolution_flags() == []
    assert workbench.check_resolution().unresolved_crossings

    caplog.clear()
    workbench.report_resolution()
    assert_logged(caplog, logging.WARNING, WORKBENCH, count=0)


def test_a_clean_crossing_raises_no_warning(caplog):
    workbench, fp = _crossing_workbench(_arc(np.pi / 2))
    caplog.set_level(logging.DEBUG, logger=WORKBENCH)
    workbench.compute_intersections([fp], infer_iterates=False)
    assert_logged(caplog, logging.WARNING, WORKBENCH, count=0)
    assert workbench.check_resolution().clean


def test_no_threshold_no_warning(caplog):
    workbench, fp = _crossing_workbench(_shallow(), threshold=None)
    caplog.set_level(logging.DEBUG, logger=WORKBENCH)
    workbench.compute_intersections([fp], infer_iterates=False)
    assert_logged(caplog, logging.WARNING, WORKBENCH, count=0)
    assert workbench.pending_resolution_flags() == []
    assert workbench.check_resolution().unresolved_crossings == ()


def test_diagnose_false_defers_the_warning(caplog):
    workbench, fp = _crossing_workbench(_shallow())
    caplog.set_level(logging.DEBUG, logger=WORKBENCH)
    workbench.compute_intersections([fp], infer_iterates=False, diagnose=False)
    assert_logged(caplog, logging.WARNING, WORKBENCH, count=0)
    assert workbench.pending_resolution_flags()
    workbench.report_resolution()
    assert_logged(caplog, logging.WARNING, WORKBENCH, count=1)


def test_anchors_and_synthetic_crossings_carry_no_ratio():
    workbench, fp = _crossing_workbench(_arc(np.pi / 2))
    workbench.compute_intersections([fp], infer_iterates=False)
    registry = workbench.intersection_registry
    assert any(ix.is_synthetic for _iid, ix in registry)
    for _iid, ix in registry:
        assert (ix.resolution_ratio is None) == ix.is_synthetic
    assert Intersection.synthetic((0.0, 0.0), 0.0, 0.0).resolution_ratio is None


def test_a_multi_bridge_blast_warns_at_most_once_per_call(caplog):
    case = build_k28(blasts=0, through="zones")
    session = case.session
    (fp,) = case.fixed_points
    (zone,) = case.zones
    # Flag every crossing, so every bridge the blast maps adds unresolved ones.
    session.workbench.resolution_ratio_threshold = float("inf")
    caplog.set_level(logging.WARNING, logger=WORKBENCH)
    result = session.blast_zone(zone, num_iterations=3, fixed_point=[fp], min_separation=1e-5)
    assert sum(len(step.interior_parents) for step in result.steps) > 1
    assert_logged(caplog, logging.WARNING, WORKBENCH, count=1)

    caplog.clear()
    session.blast_zone(zone, num_iterations=1, fixed_point=[fp], diagnose=False)
    assert_logged(caplog, logging.WARNING, WORKBENCH, count=0)



def _image_bridge(fp, points):
    """A partial image bridge through ``points`` on ``fp``'s unstable branch 0.0."""
    from tanglepack import Bridge, Point

    nodes = [Point(x, y, float(k)) for k, (x, y) in enumerate(points)]
    for node, following in zip(nodes, nodes[1:]):
        node.insert_point_forward(following)
    return Bridge(
        root=nodes[0], stability="unstable", stretch_param=1.0, fixed_point=fp,
        tail=nodes[-1], branch_index=0, manifold_key=(fp, "unstable", 0, 0),
    )


def _involved_hits(workbench, fixed_points, bridges):
    """The full listing's unstable hits naming one of ``bridges``."""
    from tanglepack.numerics.diagnostics import same_stability_crossings

    names = {f"image bridge {bridge.id}" for bridge in bridges if bridge.id is not None}
    names |= {f"partial image bridge {id(bridge):#x}" for bridge in bridges if bridge.id is None}
    return [
        hit
        for hit in same_stability_crossings(workbench, fixed_points)
        if hit.stability == "unstable" and {hit.curve_a, hit.curve_b} & names
    ]


def test_the_count_is_the_listing_of_the_new_curves():
    """C5: ``count_same_stability_crossings`` counts exactly the full listing's hits
    that name one of the new curves, in and out of the tangle region; the hits among
    old curves only are not counted."""
    from tanglepack.numerics.diagnostics import count_same_stability_crossings

    fp = bare_fixed_point(label="A")
    workbench = polyline_workbench(fp, {
        (fp, "unstable", 0, 0): [(0, 0), (30, 0)],
        (fp, "stable", 0, 0): [(1, -1), (1, 1)],   # the two crossings span the tangle box
        (fp, "stable", 0, 1): [(9, -1), (9, 1)],
    })
    workbench.compute_intersections([fp], infer_iterates=False)
    old = _image_bridge(fp, [(2, -0.5), (2.5, 0.5)])
    new = [_image_bridge(fp, [(3, -0.5), (3.5, 0.5)]), _image_bridge(fp, [(20, -1), (20, 1)])]
    for bridge in [old, *new]:
        workbench._register_bridge(bridge)

    count = count_same_stability_crossings(workbench, new)
    involved = _involved_hits(workbench, [fp], new)
    assert count.in_region == sum(hit.in_tangle_region for hit in involved) == 1
    assert count.outside == sum(not hit.in_tangle_region for hit in involved) == 1
    assert count.suggested_area_cutoff is not None
    assert count_same_stability_crossings(workbench, [old]).in_region == 1


def test_a_blast_counts_and_remembers_nothing(monkeypatch):
    """A blast's report carries the count over the bridges it added and stores no
    same-stability hit (every hit it can find involves a curve new to it)."""
    from tanglepack.numerics.diagnostics import count_same_stability_crossings

    case = build_k28(blasts=0, through="zones")
    workbench = case.workbench
    before = {id(bridge) for bridge in workbench.bridges}
    remembered = set(workbench._reported_artifacts)
    reports = []
    report_resolution = workbench.report_resolution
    monkeypatch.setattr(
        workbench,
        "report_resolution",
        lambda *args, **kwargs: reports.append(report_resolution(*args, **kwargs)) or reports[-1],
    )

    case.session.blast_zone(case.zones[0], 1, fixed_point=case.fixed_points, diagnose=True)

    added = [bridge for bridge in workbench.bridges if id(bridge) not in before]
    (report,) = reports
    assert added and report.same_stability_crossings == ()
    assert report.same_stability_count == count_same_stability_crossings(workbench, added)
    assert workbench._reported_artifacts == remembered


@pytest.mark.perf
def test_blast_count_is_the_new_bridges_hits():
    """C5 on a real deep blast (k=2.8, the ninth blast, where the hits appear): the
    count is the full listing's hits naming the bridges the blast added (~10 s)."""
    from tanglepack.numerics.diagnostics import count_same_stability_crossings

    case = build_k28(blasts=8, through="zones")
    workbench, fixed_points = case.workbench, case.fixed_points
    before = {id(bridge) for bridge in workbench.bridges}
    case.session.blast_zone(case.zones[0], 1, fixed_point=fixed_points, diagnose=False)

    added = [bridge for bridge in workbench.bridges if id(bridge) not in before]
    count = count_same_stability_crossings(workbench, added)
    involved = _involved_hits(workbench, fixed_points, added)
    assert count.in_region == sum(hit.in_tangle_region for hit in involved)
    assert count.outside == sum(not hit.in_tangle_region for hit in involved)
    assert count.in_region or count.outside, "nine blasts of k=2.8 must leave a hit to count"

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

import logging
import time

import numpy as np
import pytest

from cases import build_k10_inversion
from helpers import laws
from helpers.fakes import bare_fixed_point
from helpers.logs import assert_logged
from helpers.polylines import polyline_manifold, polyline_workbench
from tanglepack import BranchPoint
from tanglepack.numerics.diagnostics import same_stability_crossings

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
    root = BranchPoint(2, (0.0, 0.0), 5.0, 5.0)
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

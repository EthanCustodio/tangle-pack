"""NUM-12 (D1): grow, blast and recompute in any order without losing a computed crossing.

A blast's image curves are computed manifold beyond the grown branches. Before
WP10b every recompute indexed the branches only, so a recompute after a blast,
with or without growth in between, dropped every crossing born on an image
curve, the bridges cut out of it and the blast frontier, and the class count
fell (k=2.8 two blasts 4 -> 3, nested 13 -> 12).
"""

from __future__ import annotations

import re

import numpy as np
import pytest

import cases


def _recompute(case: cases.Case) -> None:
    workbench = case.session.workbench
    workbench.compute_intersections(case.fixed_points, preserve_ids=True)
    workbench.rebuild_bridges()


def _crossings(case: cases.Case) -> dict[int, tuple[np.ndarray, float]]:
    """``{id: (coords, unstable segment length)}`` of every detected crossing."""
    return {
        iid: (
            np.ravel(ix.coords).astype(float),
            float(np.linalg.norm(np.subtract(*(node.get_point() for node in ix.unstable_segment)))),
        )
        for iid, ix in case.session.workbench.intersection_registry
        if not ix.is_synthetic
    }


@pytest.mark.parametrize("stability", ["unstable", "stable"])
@pytest.mark.parametrize("name", ["k28_two_blasts", "nested"])
def test_growth_keeps_every_crossing_and_class(name: str, stability: str) -> None:
    """One growth step after the blasts: every crossing still lies on a computed curve,
    so each is registered again (by position, within its own unstable segment), and
    neither the class table nor the closed alphabet shrinks."""
    case = cases.BUILDERS[name]()
    session = case.session
    before = _crossings(case)
    classes, closed = len(session.bridge_classes()), len(session.symbolic_dynamics().classes)

    for fp in case.fixed_points:
        session.workbench.grow_n_times(fp, stability, num_iterations=1)
    _recompute(case)

    after = np.array([coords for coords, _ in _crossings(case).values()])
    lost = [
        iid for iid, (coords, length) in before.items()
        if np.min(np.linalg.norm(after - coords, axis=1)) > length
    ]
    assert not lost, f"crossings {lost} vanished although their curves are still computed"
    assert len(session.bridge_classes()) >= classes
    assert len(session.symbolic_dynamics().classes) >= closed


def test_blasting_does_not_depend_on_recomputes() -> None:
    """Two blasts, a recompute, one more blast: the same tangle as three blasts straight."""
    def signature(case: cases.Case) -> tuple:
        session = case.session
        dynamics = session.symbolic_dynamics()
        return (
            len(session.workbench.intersection_registry),
            len(session.workbench.bridges),
            len(session.bridge_classes()),
            re.sub(r"\d+", "#", dynamics.describe()),
        )

    straight = cases.build_k28(blasts=3)
    mixed = cases.build_k28(blasts=2)
    _recompute(mixed)
    cases._blast(mixed.session, mixed.zones[0], mixed.fixed_points)
    assert signature(mixed) == signature(straight)

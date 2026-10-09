"""
Regression (PK-5, D10): choosing a strong pip cuts the trellis at its orbit.

The partition and the dual graph read the stable manifold cut at the chosen
pip's orbit (D10), so the topology at a pip is the topology of a numerical trim
at that pip: same holes (by the cdists of their bounding crossings, iterate and
side), same class count, same words in element names, same reliability, or the
same exception both ways. Where the cut removes a crossing between two kept
crossings of one unstable curve, the trellis joins the two bridges across it
(an arc of ``T(q0)``, not a workbench bridge). The inversion x period-2 orbit at
10 steps, whose stable manifold runs past its pip orbit, walks every class.
No count, id or word of a build is pinned.
"""

from __future__ import annotations

import operator

import numpy as np
import pytest

from cases import BUILDERS, build_inversion_period2, build_k10
from helpers import laws
from helpers.names import homotopy_pair, word_in_names


def _signature(session) -> tuple:
    """Holes by cdists, class count, words in element names, reliability."""
    topology = session.topology()
    registry = session.workbench.intersection_registry

    def at(crossing_id: int) -> tuple:
        ix = registry[crossing_id]
        return round(float(ix.unstable_cdist), 6), round(float(ix.stable_cdist), 6)

    holes = tuple(sorted(
        (tuple(sorted(map(at, hole.bounding_ids or ()))), hole.iterate, hole.bridge_side)
        for hole in topology.holes
    ))
    dyn = topology.dynamics
    words = tuple(sorted(
        (homotopy_pair(dyn, cd), tuple(word_in_names(dyn, cd))) for cd in dyn.classes.values()
    ))
    return holes, len(topology.classes), words, dyn.is_reliable


def _outcome(session) -> object:
    """The signature, or the type of what the topology raises."""
    try:
        return _signature(session)
    except Exception as exc:  # noqa: BLE001 - compared by type, both ways
        return type(exc)


def _cut_vs_trim(build, pip: int) -> tuple:
    """The topology at ``pip`` and the topology after a numerical trim at the same pip
    on a fresh build (registry ids are reproducible within a process)."""
    view, trim = build(), build()
    view.session.classify_strong_pips()
    view.session.set_strong_pip(view.fixed_points[0], pip)
    assert trim.registry[pip].stable_cdist == view.registry[pip].stable_cdist
    trim.session.resonance_zone(pip)
    return _outcome(view.session), _outcome(trim.session)


def test_k10_cut_at_every_candidate_equals_a_trim() -> None:
    case = build_k10(through="pips")
    candidates = case.session.strong_pip_candidates(case.fixed_point)
    differ = [
        pip
        for pip in candidates
        if operator.ne(*_cut_vs_trim(lambda: build_k10(through="pips"), pip))
    ]
    assert not differ, (differ, candidates)


def _p2inv10(through: str = "bridges"):
    """Inversion x period 2 at 10 / 10 steps: its stable manifold runs past the pip orbit."""
    with np.errstate(over="ignore", invalid="ignore"):
        return build_inversion_period2(unstable_steps=10, stable_steps=10, through=through)


def test_p2inv10_cut_equals_a_trim() -> None:
    probe = _p2inv10()
    probe.session.classify_strong_pips()
    with np.errstate(over="ignore", invalid="ignore"):
        view, trim = _cut_vs_trim(_p2inv10, probe.session.strong_pip(probe.fixed_point))
    assert view == trim


def test_p2inv10_every_class_walks() -> None:
    """D3: with S cut at the pip orbit every active class walks and the words are reliable."""
    case = _p2inv10("partitioned")
    with np.errstate(over="ignore", invalid="ignore"):
        assert laws.check_every_active_class_walks(case) > 0
        assert case.session.symbolic_dynamics().is_reliable


def _joined_runs(session, fp) -> list[tuple[int, int]]:
    """``(a, b)`` for every run ``a, x_1, .., b`` of crossings consecutive along ONE
    computed unstable curve whose inner crossings the trellis dropped and whose
    ends it kept."""
    trellis = session.trellis(fp)
    kept = set(trellis.own_intersection_ids)
    curves: dict = {}
    for crossing_id, ix in session.workbench.intersection_registry:
        if ix.unstable_manifold is not None:
            curves.setdefault(id(ix.unstable_manifold), []).append(crossing_id)
    runs = []
    for ids in curves.values():
        ids.sort(key=lambda i: float(trellis.intersection(i).unstable_cdist))
        ends = [i for i, crossing_id in enumerate(ids) if crossing_id in kept]
        runs += [(ids[i], ids[j]) for i, j in zip(ends, ends[1:]) if j > i + 1]
    return runs


def test_a_crossing_cut_away_mid_arc_joins_its_bridges() -> None:
    """The first k10 candidate whose cut drops a crossing between two kept crossings of
    one curve: the trellis spans those two by one bridge the workbench does not hold."""
    case = build_k10(through="pips")
    session, fp = case.session, case.fixed_point
    for candidate in session.strong_pip_candidates(fp):
        session.set_strong_pip(fp, candidate)
        runs = _joined_runs(session, fp)
        if runs:
            break
    assert runs, "no candidate cuts a crossing away mid-arc"
    trellis = session.trellis(fp)
    computed = {bridge.id for bridge in session.workbench.bridges}
    for a, b in runs:
        bridge = trellis.bridge_between(a, b)
        assert bridge is not None and bridge.id == (a, b) and bridge.id not in computed

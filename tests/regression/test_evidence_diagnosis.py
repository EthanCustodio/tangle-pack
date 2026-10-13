"""Regression (INT-1, WP7): a walk the registered image disagrees with is diagnosed.

D3: the shortest walk IS the word and the registered image chain is the
cross-check. When the two disagree the evidence carries an
:class:`~tanglepack.topology.SymbolicDynamics.EvidenceDiagnosis`: ``"stable
end"`` when two consecutive crossings of the chain on ONE stable branch share a
crossing sign (the image passed around the end of the computed stable branch,
so its crossings with the trimmed-away part are not registered); else
``"order flip"`` when a member end and a bound of its own element have
registered images in the reverse stable order (the map preserves the order
along a branch, so the cdists there are noise, C21, WP12); else ``"walk"``. A
``"stable end"`` or ``"order flip"`` image gives no verdict (C11, C21:
``word_agrees`` None, logged at INFO, never a reliability failure by itself); a
``"walk"`` disagreement is a WARNING and makes the dynamics unreliable. The
k=2.8 zone blasted six times has a stable-end member (~3 s), blasted nine times
an order flip (~8 s, INT-11). Rule-based: the cause is checked against an
independent scan of the chain and the table, nothing of the build is pinned.

The order flip needs the trellis of the pip at its full depth. At nine blasts
the deepest point of the pip's orbit has resolution ratio 2.99, below the
threshold, so under N it shortens ``T(q0)`` by one iterate and the flipped copy
lies past the cut (Q13 a, B1); refining the blasts resolves it but also moves
the copy off the hole-bridge end. The order-flip build therefore blasts with
clause 1 of N off; the crossings of the flip itself are resolved (ratios above
5000).
"""

from __future__ import annotations

import logging

import pytest

from cases import ZONE_MIN_SEPARATION, build_k28
from tanglepack.topology.BridgeClass import _image_chain


@pytest.fixture(scope="module")
def k28_six_blasts():
    return build_k28(blasts=6).session


def _passes_a_stable_end(trellis, chain) -> bool:
    """Two consecutive chain crossings on one stable branch with one sign."""
    crossings = [chain[0][0], *(pair[1] for pair in chain)]
    for x, y in zip(crossings, crossings[1:]):
        a, b = trellis.intersection(x), trellis.intersection(y)
        if a.manifold_b_key == b.manifold_b_key and a.crossing_sign == b.crossing_sign:
            return True
    return False


def _order_flip(trellis, homotopy, bridge_class, evidence):
    """A member end and a bound of its element whose registered images reverse
    their stable order, or None."""
    ends = evidence.bridge_id if evidence.direction >= 0 else evidence.bridge_id[::-1]
    for end, ref in zip(ends, (bridge_class.source, bridge_class.target)):
        element = homotopy.element(ref)
        for bound in {element.lo_id, element.hi_id} - {None, end}:
            images = trellis.iterate(end, 1), trellis.iterate(bound, 1)
            if None in images:
                continue
            x, y, fx, fy = (trellis.intersection(i) for i in (end, bound, *images))
            if x.manifold_b_key != y.manifold_b_key or fx.manifold_b_key != fy.manifold_b_key:
                continue
            if (x.stable_cdist < y.stable_cdist) != (fx.stable_cdist < fy.stable_cdist):
                return end, bound
    return None


def _expected_cause(trellis, homotopy, cd, evidence) -> str:
    """Stable end first, then order flip, else walk."""
    if _passes_a_stable_end(trellis, _image_chain(trellis, evidence.bridge_id)):
        return "stable end"
    if _order_flip(trellis, homotopy, cd.entry.bridge_class, evidence):
        return "order flip"
    return "walk"


def test_disagreeing_evidence_carries_its_cause(k28_six_blasts, caplog):
    session = k28_six_blasts
    trellis = session.minimal_trellis().trellis
    homotopy = session.homotopy_partition()
    with caplog.at_level(logging.INFO, logger="tanglepack.topology.SymbolicDynamics"):
        dyn = session.symbolic_dynamics(rebuild=True)
    evidences = [evidence for cd in dyn.classes.values() for evidence in cd.evidence]
    diagnosed = [(cd, ev) for cd in dyn.classes.values() for ev in cd.evidence if ev.diagnosis]
    assert diagnosed, "no registered image disagrees with its walk on this build"
    for cd, evidence in diagnosed:
        expected = _expected_cause(trellis, homotopy, cd, evidence)
        assert evidence.diagnosis.cause == expected, evidence.bridge_id
        assert bool(evidence.diagnosis.crossings) == (expected != "walk")
        assert evidence.word_agrees is (None if expected != "walk" else False)
    assert all(evidence.diagnosis is None for evidence in evidences if evidence.word_agrees)

    walk_blamed = [evidence for _cd, evidence in diagnosed if evidence.word_agrees is False]
    warnings = [r for r in caplog.records if r.levelno == logging.WARNING and "registered image" in r.getMessage()]
    assert len(warnings) == len(walk_blamed)
    for cd in dyn.classes.values():
        assert (cd.verified is False) == any(evidence.word_agrees is False for evidence in cd.evidence)
    if walk_blamed:
        assert not dyn.is_reliable


def test_an_order_flip_gives_no_verdict(caplog):
    """INT-11 kind 1 (C21): at nine blasts a deep lambda-lemma copy hugs a
    hole-bridge end so closely that the registered images of a member end and its
    element's bound come out in the reverse stable order. That evidence is noise:
    no verdict, and the dynamics stays reliable."""
    case = build_k28(blasts=0, through="zones")
    session, (fp,), (zone,) = case.session, case.fixed_points, case.zones
    session.workbench.resolution_ratio_threshold = None
    for _ in range(9):
        session.blast_zone(
            zone, num_iterations=1, fixed_point=[fp], min_separation=ZONE_MIN_SEPARATION
        )
    trellis = session.minimal_trellis().trellis
    homotopy = session.homotopy_partition()
    with caplog.at_level(logging.INFO, logger="tanglepack.topology.SymbolicDynamics"):
        dyn = session.symbolic_dynamics()
    flipped = [
        (cd, ev)
        for cd in dyn.classes.values()
        for ev in cd.evidence
        if ev.diagnosis is not None
        and _expected_cause(trellis, homotopy, cd, ev) == "order flip"
    ]
    assert flipped, "precondition: a registered image reverses a member end's stable order"
    for cd, evidence in flipped:
        assert evidence.diagnosis.cause == "order flip"
        assert evidence.word_agrees is None
        end, bound = evidence.diagnosis.crossings
        x, y = (trellis.intersection(i) for i in (end, bound))
        fx, fy = (trellis.intersection(trellis.iterate(i, 1)) for i in (end, bound))
        assert (x.stable_cdist < y.stable_cdist) != (fx.stable_cdist < fy.stable_cdist)
    assert dyn.is_reliable

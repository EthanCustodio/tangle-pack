"""Regression (INT-1, WP7): a walk the registered image disagrees with is diagnosed.

D3: the shortest walk IS the word and the registered image chain is the
cross-check. When the two disagree the evidence carries an
:class:`~tanglepack.topology.SymbolicDynamics.EvidenceDiagnosis`: ``"stable
end"`` when two consecutive crossings of the chain on ONE stable branch share a
crossing sign (the image passed around the end of the computed stable branch,
so its crossings with the trimmed-away part are not registered), else
``"walk"``. Either way the dynamics is not reliable. The k=2.8 zone blasted six
times has such a member (~3 s). Rule-based: the cause is checked against an
independent scan of the chain, nothing of the build is pinned.
"""

from __future__ import annotations

import pytest

from cases import build_k28
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


def test_disagreeing_evidence_carries_its_cause(k28_six_blasts):
    session = k28_six_blasts
    trellis = session.minimal_trellis().trellis
    dyn = session.symbolic_dynamics()
    disagreeing = [
        evidence
        for cd in dyn.classes.values()
        for evidence in cd.evidence
        if evidence.word_agrees is False
    ]
    assert disagreeing, "no registered image disagrees with its walk on this build"
    for evidence in disagreeing:
        diagnosis = evidence.diagnosis
        assert diagnosis is not None, evidence.bridge_id
        chain = _image_chain(trellis, evidence.bridge_id)
        expected = "stable end" if _passes_a_stable_end(trellis, chain) else "walk"
        assert diagnosis.cause == expected, evidence.bridge_id
        assert bool(diagnosis.crossings) == (expected == "stable end")
    agreeing = [
        evidence
        for cd in dyn.classes.values()
        for evidence in cd.evidence
        if evidence.word_agrees
    ]
    assert all(evidence.diagnosis is None for evidence in agreeing)
    assert not dyn.is_reliable

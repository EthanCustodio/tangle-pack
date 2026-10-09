"""Regression (INT-1, WP7): a walk the registered image disagrees with is diagnosed.

D3: the shortest walk IS the word and the registered image chain is the
cross-check. When the two disagree the evidence carries an
:class:`~tanglepack.topology.SymbolicDynamics.EvidenceDiagnosis`: ``"stable
end"`` when two consecutive crossings of the chain on ONE stable branch share a
crossing sign (the image passed around the end of the computed stable branch,
so its crossings with the trimmed-away part are not registered), else
``"walk"``. A ``"stable end"`` image is INCOMPLETE evidence and gives no verdict
(C11, WP10: ``word_agrees`` None, logged at INFO, never a reliability failure by
itself); a ``"walk"`` disagreement is a WARNING and makes the dynamics
unreliable. The k=2.8 zone blasted six times has such a member (~3 s).
Rule-based: the cause is checked against an independent scan of the chain,
nothing of the build is pinned.
"""

from __future__ import annotations

import logging

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


def test_disagreeing_evidence_carries_its_cause(k28_six_blasts, caplog):
    session = k28_six_blasts
    trellis = session.minimal_trellis().trellis
    with caplog.at_level(logging.INFO, logger="tanglepack.topology.SymbolicDynamics"):
        dyn = session.symbolic_dynamics(rebuild=True)
    evidences = [evidence for cd in dyn.classes.values() for evidence in cd.evidence]
    diagnosed = [evidence for evidence in evidences if evidence.diagnosis is not None]
    assert diagnosed, "no registered image disagrees with its walk on this build"
    for evidence in diagnosed:
        chain = _image_chain(trellis, evidence.bridge_id)
        expected = "stable end" if _passes_a_stable_end(trellis, chain) else "walk"
        assert evidence.diagnosis.cause == expected, evidence.bridge_id
        assert bool(evidence.diagnosis.crossings) == (expected == "stable end")
        assert evidence.word_agrees is (None if expected == "stable end" else False)
    assert all(evidence.diagnosis is None for evidence in evidences if evidence.word_agrees)

    walk_blamed = [evidence for evidence in diagnosed if evidence.word_agrees is False]
    warnings = [r for r in caplog.records if r.levelno == logging.WARNING and "registered image" in r.getMessage()]
    assert len(warnings) == len(walk_blamed)
    for cd in dyn.classes.values():
        assert (cd.verified is False) == any(evidence.word_agrees is False for evidence in cd.evidence)
    if walk_blamed:
        assert not dyn.is_reliable

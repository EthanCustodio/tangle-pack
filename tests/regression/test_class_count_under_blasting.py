"""The D1 loop: blast, recompute the topology from scratch, compare.

The author's decisions of 2026-10-07: topology is a pure function of the
numerical state, recomputed from scratch after every blast (D1), and growth of
the class count while nothing else changes is a bug. The alphabet is the
minimal trellis's bridges only (D2). Since D12 was reversed (2026-10-31: use
every computed manifold) the class count may fall; the entropy may not.

Each series (:data:`cases.SERIES`) blasts one zone once at a time and records
a :class:`cases.SeriesStep` per recompute. The laws are rule-based, never a
pinned count:

* the topological entropy ``ln rho`` of the refined transition matrix never
  decreases (2026-10-31). The class count may fall: no unstable branch is cut
  at the pip orbit, so a deeper crossing can refute a reference pair and merge
  classes (k32_series blast 10), which is discovery, not a defect;
* while the hole set stays the same between two recomputes, so does the class
  count (at least one such pair must exist, or the law checked nothing), and so
  does the alphabet closed under the word map: the closure is a function of the
  dual graph alone, never of which bridges a blast happened to compute (D11);
* every class has a member the minimal trellis keeps;
* a reference pair is only ever lost to a witness the definition names,
  read off the table alone (BCG-4: a switch is a discovery, not a defect);
* a class keeps its letter while its two elements keep their bounding
  crossings, and no two such identities share a letter (C2, BCG-7);
* ``session.trellis()`` read before the recompute already shows the
  recomputed topology's trellis: the selection trellis is cut at the pips of
  its fixed points when it is read, never by whichever product runs first (D1);
* no class is ambiguous: the walk leaves from the edges its landing span meets
  and the non-crossing rule picks the side of every class with a member lobe
  (D13; the ties of k=2.8 and k=3.2 at 8+ blasts were the span's);
* a plain recompute keeps the crossing set N and the bridges a blast left: N
  is one definition on every path (B1, Q13 (a); the blast path's vote on
  near-coincident detections made k=2.8 drift from blast 6);
* no bridge ends at a detection outside N, as the blast left it or after the
  recompute: an unresolved detection lies inside a bridge (B1);
* the depth ``m`` of ``T(q0)`` never decreases (B3, Q2): a blast only adds
  links to the pip's orbit. It need not deepen either: it stalls where
  ``f^{m+1}(q0)`` is detected unresolved or not detected at all, which the
  measured series record but no law pins;
* Theorem T (A.8d, B7): after every blast, the members of a class whose
  images are registered share one ``P1``-level itinerary (stable-end chains
  skipped); the per-step count of classes tested is printed;
* walking a lobe image's base branch from its anchor enters and leaves the
  image in turn, and the iterated cut's empty stretches there are its inside
  gaps (B7b, ``helpers.laws.base_chord_misreadings``): the lobe whose chain and
  stable orders differ, ``k32_series``'s hole ``(6, 7)``, exists only at a
  blast step; the per-step count of lobes read is printed.

The per-kind counts are printed in every failure message as diagnostics.

Every law runs twice per series: on the full blast and on the representative
blast (``representatives="fewest_points"``, one bridge per refined class, the
rest held). A representative series that holds nothing is the full series
again, so its laws are a visible skip, and some series must hold a bridge (else
the policy never reached the blast). Theorem T is the premise the representative rule rests on, so it
is checked on the representative records too. Equivalence with the full series
(``pip_depth`` included) is not asserted: the representative blast holds
bridges, so its X is thinner and its topology may arrive later;
``scripts/blast_representatives.py`` measures it.

``_series`` caches the immutable records of each series and policy once per
run; no session survives the call, so no build state is shared between tests.
"""

from __future__ import annotations

from functools import cache
from typing import Optional

import pytest

from cases import SERIES, SeriesStep, blast_series, issue_marks


#: The blasts' representative policies: None is the full blast.
POLICIES = (None, "fewest_points")


@cache
def _series(name: str, representatives: Optional[str]) -> tuple[SeriesStep, ...]:
    """The records of one series under one policy, built once."""
    return blast_series(name, representatives)


def _law_series(name: str, representatives: Optional[str]) -> tuple[SeriesStep, ...]:
    """The records a law reads: a representative series that held nothing is skipped."""
    steps = _series(name, representatives)
    if representatives is not None and not any(step.held for step in steps):
        pytest.skip(f"{name}: the representative series holds nothing, so it is the full series")
    return steps


def _table(steps: tuple[SeriesStep, ...]) -> str:
    """``blasts: holes, classes (active/inert), closed, ambiguous, outside minimal, h`` per step."""
    return "; ".join(
        f"{s.blasts}: {len(s.holes)} holes, {s.classes} ({s.active}/{s.inert}), "
        f"{s.closed} closed, {s.ambiguous} ambiguous, {s.outside_minimal} outside, "
        f"h = {s.entropy:.6f}"
        for s in steps
    )


def _params(law: str) -> list:
    """One param per series, carrying the law's known issue on it."""
    return [pytest.param(name, marks=issue_marks(law, name), id=name) for name in SERIES]


@pytest.fixture(params=POLICIES, ids=lambda policy: policy or "full")
def representatives(request: pytest.FixtureRequest) -> Optional[str]:
    """Every law runs on the full blast and on the representative blast."""
    return request.param


def test_some_representative_series_holds_a_bridge() -> None:
    """The policy reaches the blast: a representative series holds a bridge at
    some step, and a full series never does."""
    assert not any(step.held for name in SERIES for step in _series(name, None))
    held = {name: [step.held for step in _series(name, "fewest_points")] for name in SERIES}
    assert any(any(counts) for counts in held.values()), held


@pytest.mark.parametrize("name", _params("series_entropy_never_decreases"))
def test_entropy_never_decreases(name: str, representatives: Optional[str]) -> None:
    """One more blast never lowers the topological entropy (the class count may fall)."""
    steps = _law_series(name, representatives)
    for before, after in zip(steps, steps[1:]):
        assert after.entropy >= before.entropy - 1e-9, _table(steps)


@pytest.mark.parametrize("name", _params("series_class_count_constant_while_holes_constant"))
def test_class_count_constant_while_holes_constant(name: str, representatives: Optional[str]) -> None:
    """D1 + D2: the same hole set gives the same class count."""
    steps = _law_series(name, representatives)
    same = [(a, b) for a, b in zip(steps, steps[1:]) if a.holes == b.holes]
    assert same, f"no two consecutive recomputes share a hole set: {_table(steps)}"
    for before, after in same:
        assert before.classes == after.classes, _table(steps)


@pytest.mark.parametrize("name", _params("series_closed_alphabet_constant_while_holes_constant"))
def test_closed_alphabet_constant_while_holes_constant(name: str, representatives: Optional[str]) -> None:
    """D11: the same hole set gives the same alphabet closed under the word map; a
    class a word names counts whether or not a blast computed a bridge for it."""
    steps = _law_series(name, representatives)
    same = [(a, b) for a, b in zip(steps, steps[1:]) if a.holes == b.holes]
    assert same, f"no two consecutive recomputes share a hole set: {_table(steps)}"
    for before, after in same:
        assert before.closed == after.closed, _table(steps)


@pytest.mark.parametrize("name", _params("series_alphabet_inside_minimal_trellis"))
def test_alphabet_inside_minimal_trellis_under_blasting(name: str, representatives: Optional[str]) -> None:
    """D2 after every blast: every class has a member the minimal trellis keeps."""
    steps = _law_series(name, representatives)
    assert all(s.outside_minimal == 0 for s in steps), _table(steps)


@pytest.mark.parametrize("name", _params("series_lost_references_are_witnessed"))
def test_lost_references_are_witnessed(name: str, representatives: Optional[str]) -> None:
    """BCG-4: every reference a blast removes is killed by a point of X the table names."""
    steps = _law_series(name, representatives)
    for step in steps:
        assert not step.unwitnessed, (step.blasts, sorted(step.unwitnessed), _table(steps))


@pytest.mark.parametrize("name", _params("series_letters_follow_element_identity"))
def test_letters_follow_element_identity(name: str, representatives: Optional[str]) -> None:
    """C2: a class identified by its elements' bounding crossings keeps its letter
    across recomputes, and a letter names one identity."""
    steps = _law_series(name, representatives)
    seen: dict = {}
    for step in steps:
        for identity, letter in step.letters:
            assert seen.setdefault(identity, letter) == letter, (step.blasts, letter, seen[identity])
    assert len(set(seen.values())) == len(seen), "two identities share a letter"


@pytest.mark.parametrize("name", _params("series_session_trellis_is_history_free"))
def test_session_trellis_is_history_free(name: str, representatives: Optional[str]) -> None:
    """D1, D10: what ``session.trellis()`` shows after a blast does not depend on
    whether a product was computed first."""
    steps = _law_series(name, representatives)
    assert all(s.history_free for s in steps), [s.blasts for s in steps if not s.history_free]


@pytest.mark.parametrize("name", _params("series_no_walk_is_ambiguous"))
def test_no_walk_is_ambiguous(name: str, representatives: Optional[str]) -> None:
    """D3 + D13: after every blast each class has one word (no class is ambiguous)."""
    steps = _law_series(name, representatives)
    assert all(s.ambiguous == 0 for s in steps), _table(steps)


@pytest.mark.parametrize("name", _params("series_recompute_keeps_the_crossing_set"))
def test_recompute_keeps_the_crossing_set(name: str, representatives: Optional[str]) -> None:
    """B1: ``compute_intersections(preserve_ids=True)`` and ``rebuild_bridges()``
    after a blast find the same N and cut the same bridges."""
    steps = _law_series(name, representatives)
    drift = [
        (s.blasts, *(f"{len(blast)} -> {len(again)}" for blast, again in zip(s.crossings, s.recomputed)))
        for s in steps
        if s.crossings != s.recomputed
    ]
    assert not drift, f"(blasts, |N|, bridges) changed by a recompute: {drift}"


@pytest.mark.parametrize("name", _params("series_bridges_end_at_crossings_of_n"))
def test_every_bridge_ends_at_crossings_of_n(name: str, representatives: Optional[str]) -> None:
    """B1: the blast's cut and the recompute's cut are both at crossings of N only."""
    steps = _law_series(name, representatives)
    stray = [
        (s.blasts, which, bridge_id)
        for s in steps
        for which, (crossings, bridges) in (("blast", s.crossings), ("recompute", s.recomputed))
        for bridge_id in bridges
        if not set(bridge_id) <= crossings
    ]
    assert not stray, f"(blasts, cut by, bridge) ending outside N: {stray}"


@pytest.mark.parametrize("name", _params("series_pip_orbit_never_shallower"))
def test_pip_orbit_never_shallower(name: str, representatives: Optional[str]) -> None:
    """B3: one more blast never shortens the registered orbit of the chosen pip."""
    depths = [s.pip_depth for s in _law_series(name, representatives)]
    assert depths == sorted(depths), depths


@pytest.mark.parametrize("name", _params("series_trusted_depth_within_pip_depth"))
def test_trusted_depth_never_exceeds_pip_depth(name: str, representatives: Optional[str]) -> None:
    """2026-10-29: every blast trusts its registered pip orbit no deeper than it
    registers it, and its ``pip_depth`` is the depth composed here."""
    steps = _law_series(name, representatives)
    depths = [(s.blasts, s.pip_depth, s.blast_depths) for s in steps if s.blast_depths]
    assert depths
    assert all(trusted <= pip == composed for _b, composed, (pip, trusted) in depths), depths


@pytest.mark.parametrize("name", _params("series_members_share_the_P1_itinerary"))
def test_members_share_the_P1_itinerary_under_blasting(name: str, representatives: Optional[str]) -> None:
    """Theorem T after every blast: no class's members read different ``P1``
    itineraries (``helpers.laws.p1_itinerary_disagreements``)."""
    steps = _law_series(name, representatives)
    print(f"{name} p1_tested per step:", [s.p1_tested for s in steps])
    disagreements = [(s.blasts, text) for s in steps for text in s.p1_disagreements]
    assert not disagreements, "\n".join(f"step {blasts}: {text}" for blasts, text in disagreements)


@pytest.mark.parametrize("name", _params("series_base_chords_inside_the_lobe_image"))
def test_base_chords_are_inside_the_lobe_image_under_blasting(
    name: str, representatives: Optional[str]
) -> None:
    """B7b after every blast: the base-branch gaps of every lobe image alternate
    out/in from the anchor, and the empty ones are exactly its inside gaps
    (``helpers.laws.base_chord_misreadings``)."""
    steps = _law_series(name, representatives)
    print(f"{name} chords_checked per step:", [s.chords_checked for s in steps])
    misreadings = [(s.blasts, text) for s in steps for text in s.chord_misreadings]
    assert not misreadings, "\n".join(f"step {blasts}: {text}" for blasts, text in misreadings)

"""
Dual-graph walks, element landings and the trellis itinerary.

Synthetic graphs (real ``StableNode`` / ``FaceNode`` objects, arcs with no
geometry) pin the walk: trivial cases, entry-then-exit recording, walls never
crossed, merged-face self-adjacency, ties, dedupe of parallel nodes,
unreachability, missing nodes, several start faces, the ``max_walks`` cap, the
landing span and the non-crossing rule's search (its marker verdict stubbed,
or only its enclosure, for the ``sigma'`` rule; the enclosure itself is in
``test_enclosure.py``). Hand-made
partition families and a duck-typed trellis pin the landing and the trellis
itinerary. Real walks are checked by the law tier (even itineraries, same-side
pairs, contained landings), never pinned.
"""

from __future__ import annotations

import logging
from typing import Callable, Optional

import pytest

from helpers.logs import assert_logged
from tanglepack.topology import DualWalk
from tanglepack.topology.DualWalk import (
    MAX_DETOUR,
    ElementLanding,
    Walk,
    WalkStep,
    land_element,
    shortest_walks,
    trellis_itinerary,
)
from tanglepack.topology.PartitionFamily import (
    HomotopyPartition,
    IteratedHomotopyPartition,
    PartitionFamily,
)
from helpers.fakes import (
    FakeDual,
    FakeTrellis,
    IntervalSpec,
    bare_fixed_point,
    make_result,
    ref,
    stable_key,
)

LOGGER = "tanglepack.topology.DualWalk"


# --------------------------------------------------------------------------- #
# Builders
# --------------------------------------------------------------------------- #
@pytest.fixture
def key():
    return stable_key(bare_fixed_point())


def L(key, index):
    return ref(key, "left", index)


def R(key, index):
    return ref(key, "right", index)


def _pairs(itinerary):
    return [(itinerary[i], itinerary[i + 1]) for i in range(0, len(itinerary), 2)]


# --------------------------------------------------------------------------- #
# shortest_walks: synthetic
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "goal_index, faces",
    [(0, "start"), (1, "shared")],
    ids=["start_is_goal", "start_face_is_goal_face"],
)
def test_trivial_walks_stay_on_the_start_face(key, goal_index, faces):
    """A goal equal to the start, or one on the start's face, needs no crossing."""
    dual = FakeDual()
    f0, f1, f2 = dual.face(), dual.face(), dual.face()
    dual.unified(f0, f1, L(key, 0), R(key, 0))
    if faces == "shared":
        dual.unified(f0, f2, L(key, 1), R(key, 1))
    search = shortest_walks(dual, L(key, 0), L(key, goal_index))
    assert search.status == "trivial"
    assert search.walk is not None
    assert search.walk.itinerary == (L(key, 0), L(key, goal_index))
    assert search.itinerary == search.walk.itinerary
    assert search.walk.steps == [] and search.walk.length == 0
    assert search.walk.faces == [f0]
    assert search.start_faces == [f0]


@pytest.mark.parametrize("entry_side", ["left", "right"])
def test_single_crossing_records_entry_then_exit(key, entry_side):
    """One unified node crossed from either side: entry element, then exit element.

    The ``right`` case is the mirror image of the ``left`` one (every face
    pair swapped), so the itinerary reads the right row, then the left.
    """
    near, far = (L, R) if entry_side == "left" else (R, L)
    dual = FakeDual()
    f0, f1, fx, fy = dual.face(), dual.face(), dual.face(), dual.face()

    def edge(kind, start_face, other_face, index):
        faces = (start_face, other_face) if entry_side == "left" else (other_face, start_face)
        return getattr(dual, kind)(*faces, L(key, index), R(key, index))

    edge("wall", f0, fx, 0)
    node = edge("unified", f0, f1, 1)
    edge("wall", f1, fy, 2)
    search = shortest_walks(dual, near(key, 0), near(key, 2))
    assert search.status == "unique" and not search.truncated
    walk = search.walk
    assert isinstance(walk, Walk)
    assert walk.itinerary == (near(key, 0), near(key, 1), far(key, 1), near(key, 2))
    assert len(walk.itinerary) % 2 == 0
    assert walk.faces == [f0, f1]
    (step,) = walk.steps
    assert isinstance(step, WalkStep)
    assert step.node is node
    assert step.entry_side == entry_side and step.exit_side != entry_side
    assert step.from_face is f0 and step.to_face is f1
    assert step.entry_element == near(key, 1) and step.exit_element == far(key, 1)
    assert walk.pairs == [(near(key, 0), near(key, 1)), (far(key, 1), near(key, 2))]
    assert walk.multiplicity == 1


def test_walls_are_never_crossed(key, caplog):
    dual = FakeDual()
    f0, f1 = dual.face(), dual.face()
    dual.wall(f0, f1, L(key, 0), R(key, 0))
    dual.wall(f1, f0, L(key, 1), R(key, 1))
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        search = shortest_walks(dual, L(key, 0), L(key, 1))
    assert search.status == "unreachable"
    assert search.walks == [] and search.walk is None
    assert search.start_faces == [f0]
    assert_logged(caplog, logging.WARNING, LOGGER)


def test_unreachable_isolated_goal(key):
    dual = FakeDual()
    f0, f1, f2, f3 = dual.face(), dual.face(), dual.face(), dual.face()
    dual.unified(f0, f1, L(key, 0), R(key, 0))
    dual.unified(f2, f3, L(key, 1), R(key, 1))
    assert shortest_walks(dual, L(key, 0), L(key, 1)).status == "unreachable"


def test_merged_face_self_adjacency_is_skipped(key):
    dual = FakeDual()
    f0, f1, fy = dual.face(), dual.face(), dual.face()
    # One merged face on BOTH sides of a unified edge: a step that goes nowhere.
    dual.unified(f0, f0, L(key, 0), R(key, 0))
    dual.unified(f0, f1, L(key, 1), R(key, 1))
    dual.wall(f1, fy, L(key, 2), R(key, 2))
    search = shortest_walks(dual, L(key, 0), L(key, 2))
    assert search.status == "unique"
    assert search.walk.itinerary == (L(key, 0), L(key, 1), R(key, 1), L(key, 2))
    assert search.walk.length == 1


def test_only_the_shortest_layer_is_kept(key):
    dual = FakeDual()
    f0, f1, f2, fy = dual.face(), dual.face(), dual.face(), dual.face()
    dual.unified(f0, f1, L(key, 0), R(key, 0))  # direct
    dual.unified(f0, f2, L(key, 1), R(key, 1))  # detour via f2
    dual.unified(f2, f1, L(key, 2), R(key, 2))
    dual.wall(f1, fy, L(key, 3), R(key, 3))
    search = shortest_walks(dual, L(key, 0), L(key, 3))
    assert search.status == "unique"
    assert search.walk.length == 1
    assert search.walk.faces == [f0, f1]


def test_tie_with_different_itineraries_is_ambiguous(key, caplog):
    dual = FakeDual()
    f0, f1, fx, fy = dual.face(), dual.face(), dual.face(), dual.face()
    dual.wall(f0, fx, L(key, 0), R(key, 0))
    dual.unified(f0, f1, L(key, 1), R(key, 1))
    dual.unified(f0, f1, L(key, 2), R(key, 2))
    dual.wall(f1, fy, L(key, 3), R(key, 3))
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        search = shortest_walks(dual, L(key, 0), L(key, 3))
    assert search.status == "ambiguous" and not search.truncated
    assert [walk.itinerary for walk in search.walks] == [
        (L(key, 0), L(key, 1), R(key, 1), L(key, 3)),
        (L(key, 0), L(key, 2), R(key, 2), L(key, 3)),
    ]
    assert all(walk.multiplicity == 1 for walk in search.walks)
    assert_logged(caplog, logging.WARNING, LOGGER)
    assert search.walk.itinerary[1] == L(key, 1)  # deterministic first candidate


def test_parallel_nodes_with_equal_elements_dedupe_to_unique(key, caplog):
    dual = FakeDual()
    f0, f1, fx, fy = dual.face(), dual.face(), dual.face(), dual.face()
    dual.wall(f0, fx, L(key, 0), R(key, 0))
    dual.unified(f0, f1, L(key, 1), R(key, 1))
    dual.unified(f0, f1, L(key, 1), R(key, 1))
    dual.wall(f1, fy, L(key, 3), R(key, 3))
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        search = shortest_walks(dual, L(key, 0), L(key, 3))
    assert search.status == "unique"
    assert len(search.walks) == 1
    assert search.walk.multiplicity == 2
    assert search.walk.itinerary == (L(key, 0), L(key, 1), R(key, 1), L(key, 3))
    assert_logged(caplog, logging.WARNING, LOGGER, count=0)


def test_no_start_node_and_no_goal_node(key, caplog):
    dual = FakeDual()
    f0, f1 = dual.face(), dual.face()
    dual.unified(f0, f1, L(key, 0), R(key, 0))
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        missing_start = shortest_walks(dual, L(key, 7), L(key, 0))
        missing_goal = shortest_walks(dual, L(key, 0), R(key, 9))
    assert missing_start.status == "no_start_node"
    assert missing_goal.status == "no_goal_node"
    assert missing_start.walks == [] and missing_goal.walks == []
    assert missing_start.start_faces == [] and missing_goal.start_faces == []
    assert_logged(caplog, logging.WARNING, LOGGER, count=2)


def test_multiple_start_faces_are_all_searched(key, caplog):
    dual = FakeDual()
    f0, f3, f1, fx, fy = dual.face(), dual.face(), dual.face(), dual.face(), dual.face()
    # The start element owns two edges, one on face f0 and one on face f3.
    dual.wall(f0, fx, L(key, 0), R(key, 0))
    dual.wall(f3, fx, L(key, 0), R(key, 5))
    dual.unified(f3, f1, L(key, 1), R(key, 1))  # only f3 leads to the goal
    dual.wall(f1, fy, L(key, 2), R(key, 2))
    with caplog.at_level(logging.INFO, logger=LOGGER):
        search = shortest_walks(dual, L(key, 0), L(key, 2))
    assert search.status == "unique"
    assert [face.index for face in search.start_faces] == [f0.index, f3.index]
    assert search.walk.faces == [f3, f1]
    assert_logged(caplog, logging.INFO, LOGGER)


def test_max_walks_cap_sets_truncated(key, caplog):
    dual = FakeDual()
    f0, f1, f2, fx, fy = dual.face(), dual.face(), dual.face(), dual.face(), dual.face()
    dual.wall(f0, fx, L(key, 0), R(key, 0))
    for index in (1, 2, 3):
        dual.unified(f0, f1, L(key, index), R(key, index))
    for index in (4, 5, 6):
        dual.unified(f1, f2, L(key, index), R(key, index))
    dual.wall(f2, fy, L(key, 9), R(key, 9))
    full = shortest_walks(dual, L(key, 0), L(key, 9), max_walks=9)
    assert full.status == "ambiguous" and not full.truncated
    assert len(full.walks) == 9
    assert all(len(walk.itinerary) == 6 for walk in full.walks)
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        capped = shortest_walks(dual, L(key, 0), L(key, 9), max_walks=4)
    assert capped.truncated
    assert capped.status == "ambiguous"
    assert len(capped.walks) == 4
    assert_logged(caplog, logging.WARNING, LOGGER)


def test_walks_are_sorted_deterministically(key):
    dual = FakeDual()
    f0, f1, fx, fy = dual.face(), dual.face(), dual.face(), dual.face()
    dual.wall(f0, fx, L(key, 0), R(key, 0))
    # Registered in descending element order; the result must still ascend.
    for index in (3, 1, 2):
        dual.unified(f0, f1, L(key, index), R(key, index))
    dual.wall(f1, fy, L(key, 4), R(key, 4))
    search = shortest_walks(dual, L(key, 0), L(key, 4))
    assert [walk.itinerary[1].element_id for walk in search.walks] == [1, 2, 3]
    explicit = shortest_walks(dual, L(key, 0), L(key, 4), fixed_points=[key[0]])
    assert [walk.itinerary for walk in explicit.walks] == [
        walk.itinerary for walk in search.walks
    ]


# --------------------------------------------------------------------------- #
# shortest_walks: the landing span and the non-crossing rule (D13)
# --------------------------------------------------------------------------- #
def _landing(target, span, from_table=(True, True)):
    """A contained landing on ``target`` with the given image span."""
    return ElementLanding(
        source=target,
        image_key=target.branch_key,
        image_side=target.side,
        span=span,
        from_table=from_table,
        target=target,
        contained=True,
        singleton=False,
        covering=(target,),
        reason=None,
    )


@pytest.fixture
def two_starts(key):
    """The start element owns the edges [0, 1] (face f0) and [1, 2] (face f3); each face
    reaches the goal's face through its own unified node, so two shortest walks tie."""
    dual = FakeDual()
    f0, f3, f1, fx, fy = (dual.face() for _ in range(5))
    dual.wall(f0, fx, L(key, 0), R(key, 0), lo_cdist=0.0, hi_cdist=1.0)
    dual.wall(f3, fx, L(key, 0), R(key, 5), lo_cdist=1.0, hi_cdist=2.0)
    dual.unified(f0, f1, L(key, 1), R(key, 1))
    dual.unified(f3, f1, L(key, 2), R(key, 2))
    dual.wall(f1, fy, L(key, 3), R(key, 3))
    return dual, f0, f3


@pytest.mark.parametrize(
    "span, from_table",
    [((0.25, 0.75), (True, True)), ((0.25, 1.0), (True, True)), ((0.25, 1.005), (True, False))],
    ids=["inside", "touching_the_neighbour", "scaled_end_past_the_neighbour"],
)
def test_the_walk_leaves_from_the_edges_the_landing_span_meets(key, two_starts, span, from_table):
    """The image of an element is its landing span, so the walk leaves from an edge the span
    overlaps: a neighbour merely touched, or overlapped by less than a scaled end's error,
    is no start."""
    dual, f0, _f3 = two_starts
    assert shortest_walks(dual, L(key, 0), L(key, 3)).status == "ambiguous"
    landings = (_landing(L(key, 0), span, from_table), _landing(L(key, 3), (0.0, 1.0)))
    search = shortest_walks(dual, L(key, 0), L(key, 3), landings=landings)
    assert search.status == "unique" and search.start_faces == [f0]
    assert search.walk.faces[0] is f0
    assert search.walk.start_node.face_on("left") is f0
    assert search.walk.goal_node in dual.stable_nodes_of(L(key, 3))


def test_a_span_meeting_no_edge_falls_back_to_every_edge(key, caplog):
    dual = FakeDual()
    f0, f3, f1, fx, fy = (dual.face() for _ in range(5))
    dual.wall(f0, fx, L(key, 0), R(key, 0))
    dual.wall(f3, fx, L(key, 0), R(key, 5), lo_cdist=1.0, hi_cdist=2.0)
    dual.unified(f3, f1, L(key, 1), R(key, 1))  # only f3 leads to the goal
    dual.wall(f1, fy, L(key, 2), R(key, 2))
    landings = (_landing(L(key, 0), (5.0, 6.0)), _landing(L(key, 2), (0.0, 1.0)))
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        search = shortest_walks(dual, L(key, 0), L(key, 2), landings=landings)
    assert search.status == "unique" and [face.index for face in search.start_faces] == [f0.index, f3.index]
    assert_logged(caplog, logging.WARNING, LOGGER, count=1)


@pytest.fixture
def tie(key):
    """Two unified nodes from f0 to f1: two shortest walks with different itineraries."""
    dual = FakeDual()
    f0, f1, fx, fy = dual.face(), dual.face(), dual.face(), dual.face()
    dual.wall(f0, fx, L(key, 0), R(key, 0))
    dual.unified(f0, f1, L(key, 1), R(key, 1))
    dual.unified(f0, f1, L(key, 2), R(key, 2))
    dual.wall(f1, fy, L(key, 3), R(key, 3))
    return dual


def test_the_markers_break_a_tie(key, tie, monkeypatch):
    monkeypatch.setattr(DualWalk, "holds_markers", lambda dual, walk, wanted: walk.itinerary[1] == L(key, 2))
    search = shortest_walks(tie, L(key, 0), L(key, 3), wanted={})
    assert search.status == "unique" and search.constrained and search.rejected == 1
    assert search.walk.itinerary == (L(key, 0), L(key, 2), R(key, 2), L(key, 3))
    assert not shortest_walks(tie, L(key, 0), L(key, 3)).constrained


def test_a_walk_that_cannot_be_judged_is_never_picked_by_order(key, tie, monkeypatch, caplog):
    """An admissible walk and one whose curve the model cannot judge (None) stay a tie."""
    monkeypatch.setattr(
        DualWalk, "holds_markers", lambda dual, walk, wanted: True if walk.itinerary[1] == L(key, 2) else None
    )
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        search = shortest_walks(tie, L(key, 0), L(key, 3), wanted={})
    assert search.status == "ambiguous" and search.rejected == 0
    assert len(search.walks) == 2
    assert_logged(caplog, logging.WARNING, LOGGER)


def test_the_markers_may_ask_for_a_detour_through_a_face_twice(key, monkeypatch):
    """The direct walk breaks the markers; the admissible one loops out of f0 through f2
    and back (faces may repeat, a node is crossed once), two crossings longer."""
    dual = FakeDual()
    f0, f1, f2, fx, fy = (dual.face() for _ in range(5))
    dual.wall(f0, fx, L(key, 0), R(key, 0))
    dual.unified(f0, f1, L(key, 1), R(key, 1))
    out = dual.unified(f0, f2, L(key, 4), R(key, 4))
    dual.unified(f2, f0, L(key, 5), R(key, 5))
    dual.wall(f1, fy, L(key, 3), R(key, 3))
    monkeypatch.setattr(
        DualWalk, "holds_markers", lambda dual, walk, wanted: bool(walk.steps) and walk.steps[0].node is out
    )
    search = shortest_walks(dual, L(key, 0), L(key, 3), wanted={})
    assert search.status == "unique" and search.rejected >= 1
    assert search.walk.length == 3
    assert search.walk.faces == [f0, f2, f0, f1]


class _HoldsNothing(FakeDual):
    """A dual graph whose every closed curve holds no marker: only ``sigma'`` can reject."""

    def enclosed(self, along, chords=(), along_halves=()):
        return [{}]


def test_a_walk_crossing_sigma_prime_is_rejected(key):
    """Two walks from the edge ``[0, 1]`` to the edge ``[2, 3]``: the one crossing the
    edge ``[1, 2]`` between them crosses ``f(sigma_b)``, which ``f(b)`` meets only at its
    ends (D13); the one going round beyond the goal does not."""
    dual = _HoldsNothing()
    fa, fb, fb2, fc, fx, fy = (dual.face() for _ in range(6))
    dual.wall(fa, fx, L(key, 0), R(key, 0), lo_cdist=0.0, hi_cdist=1.0)
    inside = dual.unified(fa, fb, L(key, 1), R(key, 1), lo_cdist=1.0, hi_cdist=2.0)
    dual.wall(fc, fy, L(key, 3), R(key, 3), lo_cdist=2.0, hi_cdist=3.0)
    beyond = dual.unified(fa, fb2, L(key, 4), R(key, 4), lo_cdist=3.0, hi_cdist=4.0)
    dual.unified(fb, fc, L(key, 5), R(key, 5), lo_cdist=5.0, hi_cdist=6.0)
    dual.unified(fb2, fc, L(key, 6), R(key, 6), lo_cdist=6.0, hi_cdist=7.0)
    assert shortest_walks(dual, L(key, 0), L(key, 3)).status == "ambiguous"
    search = shortest_walks(dual, L(key, 0), L(key, 3), wanted={})
    across, around = sorted(
        shortest_walks(dual, L(key, 0), L(key, 3)).walks, key=lambda walk: walk.steps[0].node is beyond
    )
    assert across.steps[0].node is inside and DualWalk.holds_markers(dual, across, {}) is False
    assert DualWalk.holds_markers(dual, around, {}) is True
    assert search.status == "unique" and search.rejected == 1
    assert search.walk.steps[0].node is beyond


def test_a_budget_spent_before_any_walk_holds_the_markers_is_truncated(key, tie, monkeypatch, caplog):
    """Out of budget with nothing admitted is no verdict on the walks never judged."""
    monkeypatch.setattr(DualWalk, "holds_markers", lambda dual, walk, wanted: False)
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        search = shortest_walks(tie, L(key, 0), L(key, 3), wanted={}, max_walks=1)
    assert search.status == "truncated" and search.truncated and search.walk is None
    assert_logged(caplog, logging.WARNING, LOGGER)


def test_walks_of_an_admitted_itinerary_are_counted_not_judged(key, monkeypatch):
    """Parallel nodes carrying the same elements: one judgement admits the itinerary,
    the other walk adds to its multiplicity and spends no budget."""
    dual = FakeDual()
    f0, f1, fx, fy = dual.face(), dual.face(), dual.face(), dual.face()
    dual.wall(f0, fx, L(key, 0), R(key, 0))
    dual.unified(f0, f1, L(key, 1), R(key, 1))
    dual.unified(f0, f1, L(key, 1), R(key, 1))
    dual.wall(f1, fy, L(key, 3), R(key, 3))
    judged = []
    monkeypatch.setattr(DualWalk, "holds_markers", lambda dual, walk, wanted: judged.append(walk) or True)
    search = shortest_walks(dual, L(key, 0), L(key, 3), wanted={}, max_walks=1)
    assert search.status == "unique" and not search.truncated
    assert len(judged) == 1 and search.walk.multiplicity == 2


def test_no_admissible_walk_within_the_detour_cap(key, tie, monkeypatch, caplog):
    monkeypatch.setattr(DualWalk, "holds_markers", lambda dual, walk, wanted: False)
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        search = shortest_walks(tie, L(key, 0), L(key, 3), wanted={})
    assert search.status == "no_admissible_walk" and search.constrained
    assert search.walks == [] and search.walk is None and search.rejected == 2
    assert MAX_DETOUR >= 1
    assert_logged(caplog, logging.WARNING, LOGGER)


def _crossing_into(node) -> Callable:
    """A ``chords_cross`` stub: a chord crosses any wall iff it runs into ``node``'s edge
    (a wall's own chord is never read)."""
    return lambda face, chord, wall: chord[2] == node.edge_key


def test_a_wall_breaks_a_tie(key, tie, caplog):
    """D13 / D14: images of disjoint arcs are disjoint, so a walk crossing the image walk
    of an arc on another unstable branch is rejected; the clear walk is the word."""
    (through_1,) = tie.stable_nodes_of(L(key, 1))
    tie.chords_cross = _crossing_into(through_1)
    walls = [(tie.faces[0], ("wall chord",))]
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        search = shortest_walks(tie, L(key, 0), L(key, 3), walls=walls)
    assert search.status == "unique" and search.rejected == 1 and not search.constrained
    assert search.walk.itinerary == (L(key, 0), L(key, 2), R(key, 2), L(key, 3))
    assert not [record for record in caplog.records if record.levelno >= logging.WARNING]


def test_a_wall_crossing_every_walk_keeps_the_tie(key, tie, caplog):
    """No clear walk: nothing is rejected and the tie stays flagged."""
    tie.chords_cross = lambda face, chord, wall: True
    walls = [(tie.faces[0], ("wall chord",))]
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        search = shortest_walks(tie, L(key, 0), L(key, 3), walls=walls)
    assert search.status == "ambiguous" and search.rejected == 0 and len(search.walks) == 2
    assert_logged(caplog, logging.WARNING, LOGGER)


def test_a_walk_standing_for_several_face_paths_is_never_rejected_by_a_wall(key, tie, caplog):
    """A second node carrying ``L_1`` gives that itinerary two face paths: which one the
    image runs through is unknown, so a wall crossing both rejects neither."""
    tie.unified(tie.faces[0], tie.faces[1], L(key, 1), R(key, 1))
    through_1 = {node.edge_key for node in tie.stable_nodes_of(L(key, 1))}
    tie.chords_cross = lambda face, chord, wall: chord[2] in through_1
    walls = [(tie.faces[0], ("wall chord",))]
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        search = shortest_walks(tie, L(key, 0), L(key, 3), walls=walls)
    assert search.status == "ambiguous" and search.rejected == 0
    assert sorted(walk.multiplicity for walk in search.walks) == [1, 2]
    assert_logged(caplog, logging.WARNING, LOGGER)


def test_no_walls_and_walls_elsewhere_change_nothing(key, tie):
    """Without walls, or with walls only in a face no walk visits, the tie is the old one."""
    (through_1,) = tie.stable_nodes_of(L(key, 1))
    tie.chords_cross = _crossing_into(through_1)
    plain = shortest_walks(tie, L(key, 0), L(key, 3))
    elsewhere = shortest_walks(tie, L(key, 0), L(key, 3), walls=[(tie.faces[2], ("wall chord",))])
    for search in (plain, elsewhere):
        assert search.status == "ambiguous" and search.rejected == 0
        assert [walk.itinerary for walk in search.walks] == [walk.itinerary for walk in plain.walks]


def test_walk_chords_run_from_edge_to_edge_in_each_face(key, tie):
    """One chord per face visited: in by the start edge (or the edge crossed), out by the
    edge crossed next (or the goal edge), each end on the side the face lies on."""
    (walk, _other) = shortest_walks(tie, L(key, 0), L(key, 3)).walks
    (start,) = tie.stable_nodes_of(L(key, 0))
    (goal,) = tie.stable_nodes_of(L(key, 3))
    (step,) = walk.steps
    assert DualWalk.walk_chords(walk) == [
        (walk.faces[0], (start.edge_key, "left", step.node.edge_key, step.entry_side)),
        (walk.faces[1], (step.node.edge_key, step.exit_side, goal.edge_key, "left")),
    ]

# --------------------------------------------------------------------------- #
# land_element: synthetic
# --------------------------------------------------------------------------- #
#: Crossings on the one stable branch: id -> stable cdist.
CDISTS = {11: 0.25, 12: 0.5, 1: 1.0, 2: 2.0, 3: 4.0, 14: 0.375}
#: The +1 iterates (per_step_beta("stable") = 0.25).
ITERATES = {1: 11, 2: 12, 3: 1}

#: One family's intervals (anchor outward) and its crossing owners, spelled out.
Family = tuple[list[IntervalSpec], dict[int, int]]

HOMOTOPY: Family = (
    [
        (None, 1, 0.0, 1.0, True, False),
        (1, 2, 1.0, 2.0, True, True),
        (2, 3, 2.0, 4.0, False, True),
    ],
    {1: 1, 2: 1, 3: 2},
)
ITERATED: Family = (
    [
        (None, 11, 0.0, 0.25, True, False),
        (11, 12, 0.25, 0.5, True, True),
        (12, 1, 0.5, 1.0, False, False),
        (1, 2, 1.0, 2.0, True, True),
        (2, 3, 2.0, 4.0, False, True),
    ],
    {11: 1, 12: 1, 1: 3, 2: 3, 3: 4},
)
ITERATED_PARENTS = [0, 0, 0, 1, 2]
#: A homotopy family with a singleton pinched at crossing 2.
PINCHED_HOMOTOPY: Family = (
    [
        (None, 1, 0.0, 1.0, True, False),
        (1, 2, 1.0, 2.0, True, False),
        (2, 2, 2.0, 2.0, True, True),
        (2, 3, 2.0, 4.0, False, True),
    ],
    {1: 1, 2: 2, 3: 3},
)


def _families(
    trellis: FakeTrellis,
    key: tuple,
    homotopy: Family = HOMOTOPY,
    iterated: Family = ITERATED,
    parents: Optional[list[int]] = None,
) -> tuple[HomotopyPartition, IteratedHomotopyPartition]:
    """Both families on both sides of one branch (mirrored)."""
    homotopy_family = HomotopyPartition.from_results(
        [make_result(key, side, *homotopy) for side in ("left", "right")],
        trellis=trellis,
    )
    iterated_family = IteratedHomotopyPartition(
        [make_result(key, side, *iterated, parents=parents) for side in ("left", "right")],
        trellis=trellis,
        homotopy=homotopy_family,
    )
    return homotopy_family, iterated_family


def _trellis(
    cdists: dict[int, float] = CDISTS,
    iterates: dict[int, int] = ITERATES,
    keys: "tuple | Callable[[tuple], dict[int, tuple]]" = (),
    **kwargs,
) -> tuple[FakeTrellis, tuple]:
    """A fake trellis on one stable branch of a real ``beta = 0.25`` fixed point.

    Args:
        cdists: Stable cdist per crossing id.
        iterates: ``+1`` iterate per crossing id.
        keys: The branch key of every crossing (the fixed point's own branch
            when empty), or a function of that branch key returning the
            per-crossing keys.
        **kwargs: Passed to :class:`FakeTrellis`.

    Returns:
        ``(trellis, the fixed point's stable branch key)``.
    """
    fixed_point = bare_fixed_point(beta=0.25)
    key = stable_key(fixed_point)
    if callable(keys):
        keys = keys(key)
    return FakeTrellis(cdists, iterates, keys or key, fixed_points=[fixed_point], **kwargs), key


@pytest.fixture
def landing_setup():
    trellis, key = _trellis()
    homotopy, iterated = _families(trellis, key, parents=ITERATED_PARENTS)
    return trellis, key, homotopy, iterated


@pytest.mark.parametrize(
    "element, span",
    [(0, (0.0, 0.25)), (1, (0.25, 0.5)), (2, (0.5, 1.0))],
    ids=["unbounded_low_end_maps_to_zero", "interior", "outermost"],
)
def test_landing_from_the_table_is_contained(landing_setup, element, span):
    """Every homotopy element lands, by the table alone, inside one iterated element."""
    trellis, key, homotopy, iterated = landing_setup
    landing = land_element(trellis, homotopy, iterated, L(key, element))
    assert isinstance(landing, ElementLanding)
    assert landing.resolved and landing.contained and not landing.singleton
    assert landing.source == L(key, element)
    assert landing.image_key == key and landing.image_side == "left"
    assert landing.span == span
    assert landing.from_table == (True, True)
    assert landing.target == L(key, element)
    assert landing.covering == (L(key, element),)
    assert landing.reason is None
    assert repr(landing)


def test_unbounded_high_end_is_scaled():
    trellis, key = _trellis()
    open_ended: Family = (HOMOTOPY[0][:2] + [(2, None, 2.0, 4.0, False, True)], {1: 1, 2: 1})
    homotopy, iterated = _families(trellis, key, homotopy=open_ended)
    landing = land_element(trellis, homotopy, iterated, L(key, 2))
    assert landing.span == (0.5, 1.0)
    assert landing.from_table == (True, False)  # 2 -> 12 registered; hi = 4.0 * 0.25 scaled
    assert landing.target == L(key, 2) and landing.contained


def test_scaled_end_tolerance_is_relative():
    # The unregistered crossing 4 sits at 4.004: its scaled image 1.001 overshoots
    # the iterated element (0.5, 1.0) by 1e-3 -- inside SCALING_RTOL, so contained.
    trellis, key = _trellis(cdists={**CDISTS, 4: 4.004})
    open_ended: Family = (
        HOMOTOPY[0][:2] + [(2, 4, 2.0, 4.004, False, True)],
        {1: 1, 2: 1, 4: 2},
    )
    homotopy, iterated = _families(trellis, key, homotopy=open_ended)
    landing = land_element(trellis, homotopy, iterated, L(key, 2))
    assert landing.from_table == (True, False)
    assert landing.target == L(key, 2) and landing.contained


def test_orientation_reversing_map_flips_the_side():
    trellis, key = _trellis(orientation_preserving=False)
    homotopy, iterated = _families(trellis, key)
    landing = land_element(trellis, homotopy, iterated, L(key, 1))
    assert landing.image_side == "right"
    assert landing.target == R(key, 1) and landing.contained
    preserving, key = _trellis()
    assert land_element(preserving, *_families(preserving, key), L(key, 1)).target == L(key, 1)


def test_zero_width_span_probes_its_low_end():
    trellis, key = _trellis(cdists={**CDISTS, 9: 1.0, 19: 0.25}, iterates={**ITERATES, 9: 19})
    # A hand-made non-singleton element whose two ends coincide in cdist.
    degenerate: Family = (
        [(None, 1, 0.0, 1.0, True, False), (1, 9, 1.0, 1.0, True, True)],
        {1: 1, 9: 1},
    )
    homotopy, iterated = _families(trellis, key, homotopy=degenerate)
    landing = land_element(trellis, homotopy, iterated, L(key, 1))
    assert landing.span == (0.25, 0.25)
    assert not landing.singleton
    assert landing.target == L(key, 1) and landing.contained


def test_singleton_lands_on_the_registered_iterate():
    trellis, key = _trellis()
    pinched_iterated: Family = (
        [
            (None, 11, 0.0, 0.25, True, False),
            (11, 12, 0.25, 0.5, True, False),
            (12, 12, 0.5, 0.5, True, True),
            (12, 1, 0.5, 1.0, False, False),
            (1, 2, 1.0, 2.0, True, False),
            (2, 2, 2.0, 2.0, True, True),
            (2, 3, 2.0, 4.0, False, True),
        ],
        {11: 1, 12: 2, 1: 4, 2: 5, 3: 6},
    )
    homotopy, iterated = _families(trellis, key, PINCHED_HOMOTOPY, pinched_iterated)
    landing = land_element(trellis, homotopy, iterated, L(key, 2))
    assert landing.singleton and landing.resolved and landing.contained
    assert landing.target == L(key, 2)
    assert landing.span == (0.5, 0.5) and landing.from_table == (True, True)
    assert landing.covering == (L(key, 2),)
    assert repr(landing)


def test_singleton_with_unregistered_image_is_unresolved(caplog):
    trellis, key = _trellis(iterates={1: 11, 3: 1})  # no 2 -> 12
    homotopy, iterated = _families(trellis, key, PINCHED_HOMOTOPY)
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        landing = land_element(trellis, homotopy, iterated, L(key, 2))
    assert landing.singleton and not landing.resolved
    assert landing.target is None and landing.covering == ()
    assert landing.reason
    assert landing.span == (0.5, 0.5) and landing.from_table == (False, False)
    assert_logged(caplog, logging.WARNING, LOGGER)
    assert repr(landing)


def test_singleton_owned_by_a_non_singleton_is_info_not_warning(caplog):
    """A singleton whose image is not pinched is the backward-only-holes rule at
    work (the forward hole is never punched), so it is reported at INFO and the
    landing is resolved normally."""
    trellis, key = _trellis()
    # The iterated partition does not pinch 12: element #1 owns it.
    homotopy, iterated = _families(trellis, key, PINCHED_HOMOTOPY, ITERATED)
    with caplog.at_level(logging.INFO, logger=LOGGER):
        landing = land_element(trellis, homotopy, iterated, L(key, 2))
    assert landing.singleton and landing.resolved and landing.contained
    assert landing.target == L(key, 1)
    assert_logged(caplog, logging.WARNING, LOGGER, count=0)
    assert_logged(caplog, logging.INFO, LOGGER)


def test_span_straddling_two_elements_fills_covering(caplog):
    trellis, key = _trellis()
    finer: Family = (
        [
            (None, 11, 0.0, 0.25, True, False),
            (11, 14, 0.25, 0.375, True, True),
            (14, 12, 0.375, 0.5, False, True),
            (12, 1, 0.5, 1.0, False, False),
            (1, 2, 1.0, 2.0, True, True),
            (2, 3, 2.0, 4.0, False, True),
        ],
        {11: 1, 14: 1, 12: 2, 1: 4, 2: 4, 3: 5},
    )
    homotopy, iterated = _families(trellis, key, iterated=finer)
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        landing = land_element(trellis, homotopy, iterated, L(key, 1))
    assert landing.resolved and not landing.contained
    assert landing.span == (0.25, 0.5)
    assert landing.target == L(key, 1)  # the midpoint 0.375 is owned by the closed element
    assert landing.covering == (L(key, 1), L(key, 2))
    assert landing.reason is None
    assert_logged(caplog, logging.WARNING, LOGGER)
    assert repr(landing)


def test_image_on_an_unexpected_branch_asserts():
    """The iterate table links only onto the advanced branches and the scaled
    fallback is advance_key itself, so an image elsewhere is a broken table:
    asserted, never a resolvable state (WP11)."""
    other = stable_key(bare_fixed_point(2))

    def keys(key: tuple) -> dict[int, tuple]:
        # The image of crossing 2 sits on a foreign branch.
        return {**{iid: key for iid in CDISTS}, 12: other}

    trellis, key = _trellis(keys=keys)
    homotopy, iterated = _families(trellis, key)
    with pytest.raises(AssertionError, match="imaged onto"):
        land_element(trellis, homotopy, iterated, L(key, 1))


def test_probe_outside_the_iterated_partition_is_unresolved(caplog):
    trellis, key = _trellis()
    short: Family = (ITERATED[0][:2], {11: 1, 12: 1})  # the iterated partition stops at 0.5
    homotopy, iterated = _families(trellis, key, iterated=short)
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        landing = land_element(trellis, homotopy, iterated, L(key, 2))
    assert not landing.resolved and landing.target is None
    assert landing.reason
    assert landing.span == (0.5, 1.0)


def test_land_element_refuses_an_unknown_ref(landing_setup):
    trellis, key, homotopy, iterated = landing_setup
    with pytest.raises(IndexError):
        land_element(trellis, homotopy, iterated, L(key, 42))


# --------------------------------------------------------------------------- #
# trellis_itinerary: synthetic
# --------------------------------------------------------------------------- #
class _OwnerFamily:
    """Owner of crossing ``c`` on side ``s`` is element ``c`` on that side."""

    #: The production rule: the owner on the end's row, on the branch it meets.
    owner_at_end = PartitionFamily.owner_at_end

    def __init__(self, key, missing=()):
        self.key = key
        self.missing = set(missing)

    def owner_of_intersection(self, crossing_id, side, branch_key):
        assert branch_key == self.key
        if crossing_id in self.missing:
            raise ValueError(f"crossing {crossing_id} is on no {side!r}-partitioned branch")
        return ref(self.key, side, crossing_id)


@pytest.fixture
def chain_setup(key):
    # Signs alternate along the image arc: rows are then consistent per bridge.
    signs = {0: +1, 1: -1, 2: +1, 3: -1}
    trellis = FakeTrellis({i: float(i) for i in signs}, {}, key, signs=signs)
    chain = [(0, 1), (1, 2), (2, 3)]
    return trellis, chain


def test_trellis_itinerary_is_even_and_alternates_rows(key, chain_setup):
    trellis, chain = chain_setup
    itinerary = trellis_itinerary(trellis, _OwnerFamily(key), chain, +1)
    assert len(itinerary) == 2 * len(chain) == 6
    assert itinerary == (
        L(key, 0),  # bridge (0,1) leaves 0 along u+: sign(0) > 0 -> left
        L(key, 1),  # enters 1 along u-: sign(1) < 0 -> left
        R(key, 1),  # bridge (1,2) leaves 1: sign(1) > 0 fails -> right
        R(key, 2),  # enters 2: sign(2) < 0 fails -> right
        L(key, 2),  # bridge (2,3) leaves 2: sign(2) > 0 -> left
        L(key, 3),  # enters 3: sign(3) < 0 -> left
    )
    # Interior crossings: entry row then exit row, on opposite sides.
    for index in (1, 3):
        assert itinerary[index].element_id == itinerary[index + 1].element_id
        assert itinerary[index].side != itinerary[index + 1].side
    # Disjoint pairs are single bridges: same side at both ends.
    assert all(a.side == b.side for a, b in _pairs(itinerary))


def test_trellis_itinerary_direction_minus_one_reverses(key, chain_setup):
    trellis, chain = chain_setup
    forward = trellis_itinerary(trellis, _OwnerFamily(key), chain, +1)
    backward = trellis_itinerary(trellis, _OwnerFamily(key), chain, -1)
    assert backward == tuple(reversed(forward))
    # Read backward the walk starts at c3 on the row bridge (2,3) arrived on,
    # and every pair is the swapped forward pair in reverse order.
    assert backward[0] == L(key, 3) and backward[-1] == L(key, 0)
    assert _pairs(backward) == [(b, a) for a, b in reversed(_pairs(forward))]


def test_trellis_itinerary_single_pair(key):
    signs = {5: +1, 6: -1}
    trellis = FakeTrellis({5: 1.0, 6: 2.0}, {}, key, signs=signs)
    assert trellis_itinerary(trellis, _OwnerFamily(key), [(5, 6)], +1) == (L(key, 5), L(key, 6))


def test_trellis_itinerary_errors(key, chain_setup):
    trellis, chain = chain_setup
    with pytest.raises(ValueError):
        trellis_itinerary(trellis, _OwnerFamily(key), [], +1)
    with pytest.raises(ValueError):
        trellis_itinerary(trellis, _OwnerFamily(key), [(0, 1), (2, 3)], +1)
    with pytest.raises(ValueError):
        trellis_itinerary(trellis, _OwnerFamily(key, missing={2}), chain, +1)
    unsigned = FakeTrellis({0: 0.0, 1: 1.0}, {}, key, signs={0: 0, 1: -1})
    with pytest.raises(ValueError):
        trellis_itinerary(unsigned, _OwnerFamily(key), [(0, 1)], +1)

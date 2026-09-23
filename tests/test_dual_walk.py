"""
Dual-graph walks, element landings and the trellis itinerary.

Synthetic graphs (real ``StableNode`` / ``FaceNode`` objects, ``arc=None``)
pin the walk: trivial cases, entry-then-exit recording, walls never crossed,
merged-face self-adjacency, ties, dedupe of parallel nodes, unreachability,
missing nodes, several start faces and the ``max_walks`` cap. Hand-made
partition families and a duck-typed trellis pin the landing and the trellis
itinerary. One fixture test on k=10 pins the real thing without registry ids.
"""

from __future__ import annotations

import logging

import pytest

from minimal_helpers import build_pieces
from tanglepack.topology.BridgeClass import _image_chain
from tanglepack.topology.DualGraph import DualGraph
from tanglepack.topology.DualWalk import (
    ElementLanding,
    Walk,
    WalkSearch,
    WalkStep,
    land_element,
    shortest_walks,
    trellis_itinerary,
)
from tanglepack.topology.PartitionFamily import (
    HomotopyPartition,
    IteratedHomotopyPartition,
)
from walk_helpers import (
    FakeDual,
    FakeFixedPoint,
    FakeTrellis,
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
    return stable_key(FakeFixedPoint())


def L(key, index):
    return ref(key, "left", index)


def R(key, index):
    return ref(key, "right", index)


def _pairs(itinerary):
    return [(itinerary[i], itinerary[i + 1]) for i in range(0, len(itinerary), 2)]


def _warnings(caplog):
    return [record for record in caplog.records if record.levelno == logging.WARNING]


# --------------------------------------------------------------------------- #
# shortest_walks: synthetic
# --------------------------------------------------------------------------- #
def test_start_equals_goal_is_trivial(key):
    dual = FakeDual()
    f0, f1 = dual.face(), dual.face()
    dual.unified(f0, f1, L(key, 0), R(key, 0))
    search = shortest_walks(dual, L(key, 0), L(key, 0))
    assert search.status == "trivial"
    assert search.walk is not None
    assert search.walk.itinerary == (L(key, 0), L(key, 0))
    assert search.walk.steps == [] and search.walk.length == 0
    assert search.walk.faces == [f0]
    assert search.start_faces == [f0]


def test_start_face_equal_to_goal_face_is_trivial(key):
    dual = FakeDual()
    f0, f1, f2 = dual.face(), dual.face(), dual.face()
    dual.unified(f0, f1, L(key, 0), R(key, 0))
    dual.unified(f0, f2, L(key, 1), R(key, 1))
    search = shortest_walks(dual, L(key, 0), L(key, 1))
    assert search.status == "trivial"
    assert search.walk.itinerary == (L(key, 0), L(key, 1))
    assert search.walk.faces == [f0]
    assert search.itinerary == search.walk.itinerary


def test_single_crossing_records_entry_then_exit(key):
    dual = FakeDual()
    f0, f1, fx, fy = dual.face(), dual.face(), dual.face(), dual.face()
    dual.wall(f0, fx, L(key, 0), R(key, 0))
    node = dual.unified(f0, f1, L(key, 1), R(key, 1))
    dual.wall(f1, fy, L(key, 2), R(key, 2))
    search = shortest_walks(dual, L(key, 0), L(key, 2))
    assert search.status == "unique" and not search.truncated
    walk = search.walk
    assert isinstance(walk, Walk)
    assert walk.itinerary == (L(key, 0), L(key, 1), R(key, 1), L(key, 2))
    assert len(walk.itinerary) % 2 == 0
    assert walk.faces == [f0, f1]
    (step,) = walk.steps
    assert isinstance(step, WalkStep)
    assert step.node is node
    assert step.entry_side == "left" and step.exit_side == "right"
    assert step.from_face is f0 and step.to_face is f1
    assert step.entry_element == L(key, 1) and step.exit_element == R(key, 1)
    assert walk.pairs == [(L(key, 0), L(key, 1)), (R(key, 1), L(key, 2))]
    assert walk.multiplicity == 1


def test_crossing_from_the_right_side_records_right_then_left(key):
    dual = FakeDual()
    f0, f1, fx, fy = dual.face(), dual.face(), dual.face(), dual.face()
    dual.wall(fx, f0, L(key, 0), R(key, 0))
    dual.unified(f1, f0, L(key, 1), R(key, 1))
    dual.wall(fy, f1, L(key, 2), R(key, 2))
    search = shortest_walks(dual, R(key, 0), R(key, 2))
    assert search.status == "unique"
    assert search.walk.itinerary == (R(key, 0), R(key, 1), L(key, 1), R(key, 2))
    assert search.walk.steps[0].entry_side == "right"


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
    assert any("unreachable" in record.getMessage() for record in _warnings(caplog))


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
    assert any("ambiguous" in record.getMessage() for record in _warnings(caplog))
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
    assert not _warnings(caplog)


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
    messages = [record.getMessage() for record in _warnings(caplog)]
    assert any("start element" in message for message in messages)
    assert any("goal element" in message for message in messages)


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
    assert any(
        record.levelno == logging.INFO and "faces 2 faces" in record.getMessage()
        for record in caplog.records
    )


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
    assert any("truncated" in record.getMessage() for record in _warnings(caplog))


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


def test_walk_search_repr_and_walk_repr(key):
    dual = FakeDual()
    f0, f1 = dual.face(), dual.face()
    dual.unified(f0, f1, L(key, 0), R(key, 0))
    search = shortest_walks(dual, L(key, 0), R(key, 0))
    assert isinstance(search, WalkSearch)
    assert "unique" in repr(search)
    assert "1 step(s)" in repr(search.walk)


# --------------------------------------------------------------------------- #
# land_element: synthetic
# --------------------------------------------------------------------------- #
#: Crossings on the one stable branch: id -> stable cdist.
CDISTS = {11: 0.25, 12: 0.5, 1: 1.0, 2: 2.0, 3: 4.0, 14: 0.375}
#: The +1 iterates (per_step_beta("stable") = 0.25).
ITERATES = {1: 11, 2: 12, 3: 1}

HOMOTOPY = [
    (None, 1, 0.0, 1.0, True, False),
    (1, 2, 1.0, 2.0, True, True),
    (2, 3, 2.0, 4.0, False, True),
]
ITERATED = [
    (None, 11, 0.0, 0.25, True, False),
    (11, 12, 0.25, 0.5, True, True),
    (12, 1, 0.5, 1.0, False, False),
    (1, 2, 1.0, 2.0, True, True),
    (2, 3, 2.0, 4.0, False, True),
]
ITERATED_PARENTS = [0, 0, 0, 1, 2]


def _families(trellis, key, homotopy_specs=HOMOTOPY, iterated_specs=ITERATED, parents=None):
    """Both families on both sides of one branch (mirrored)."""
    homotopy = HomotopyPartition.from_results(
        [make_result(key, side, homotopy_specs) for side in ("left", "right")],
        trellis=trellis,
    )
    iterated = IteratedHomotopyPartition(
        [make_result(key, side, iterated_specs, parents=parents) for side in ("left", "right")],
        trellis=trellis,
        homotopy=homotopy,
    )
    return homotopy, iterated


@pytest.fixture
def landing_setup():
    fixed_point = FakeFixedPoint(beta=0.25)
    key = stable_key(fixed_point)
    trellis = FakeTrellis(CDISTS, ITERATES, key, fixed_points=[fixed_point])
    homotopy, iterated = _families(trellis, key, parents=ITERATED_PARENTS)
    return trellis, key, homotopy, iterated


def test_contained_landing_from_the_table(landing_setup):
    trellis, key, homotopy, iterated = landing_setup
    landing = land_element(trellis, homotopy, iterated, L(key, 1))
    assert isinstance(landing, ElementLanding)
    assert landing.resolved and landing.contained and not landing.singleton
    assert landing.source == L(key, 1)
    assert landing.image_key == key and landing.image_side == "left"
    assert landing.span == (0.25, 0.5)
    assert landing.from_table == (True, True)
    assert landing.target == L(key, 1)
    assert landing.covering == (L(key, 1),)
    assert landing.reason is None
    assert "L#1 -> " in repr(landing)


def test_unbounded_low_end_maps_to_zero(landing_setup):
    trellis, key, homotopy, iterated = landing_setup
    landing = land_element(trellis, homotopy, iterated, L(key, 0))
    assert landing.span == (0.0, 0.25)
    assert landing.from_table == (True, True)
    assert landing.target == L(key, 0) and landing.contained


def test_outermost_element_lands_by_the_table(landing_setup):
    trellis, key, homotopy, iterated = landing_setup
    landing = land_element(trellis, homotopy, iterated, L(key, 2))
    assert landing.span == (0.5, 1.0)
    assert landing.target == L(key, 2) and landing.contained


def test_unbounded_high_end_is_scaled():
    fixed_point = FakeFixedPoint(beta=0.25)
    key = stable_key(fixed_point)
    trellis = FakeTrellis(CDISTS, ITERATES, key, fixed_points=[fixed_point])
    open_ended = HOMOTOPY[:2] + [(2, None, 2.0, 4.0, False, True)]
    homotopy, iterated = _families(trellis, key, homotopy_specs=open_ended)
    landing = land_element(trellis, homotopy, iterated, L(key, 2))
    assert landing.span == (0.5, 1.0)
    assert landing.from_table == (True, False)  # 2 -> 12 registered; hi = 4.0 * 0.25 scaled
    assert landing.target == L(key, 2) and landing.contained


def test_scaled_end_tolerance_is_relative():
    fixed_point = FakeFixedPoint(beta=0.25)
    key = stable_key(fixed_point)
    # The unregistered crossing 4 sits at 4.004: its scaled image 1.001 overshoots
    # the iterated element (0.5, 1.0) by 1e-3 -- inside SCALING_RTOL, so contained.
    cdists = {**CDISTS, 4: 4.004}
    trellis = FakeTrellis(cdists, ITERATES, key, fixed_points=[fixed_point])
    open_ended = HOMOTOPY[:2] + [(2, 4, 2.0, 4.004, False, True)]
    homotopy, iterated = _families(trellis, key, homotopy_specs=open_ended)
    landing = land_element(trellis, homotopy, iterated, L(key, 2))
    assert landing.from_table == (True, False)
    assert landing.target == L(key, 2) and landing.contained


def test_orientation_reversing_map_flips_the_side():
    fixed_point = FakeFixedPoint(beta=0.25)
    key = stable_key(fixed_point)
    trellis = FakeTrellis(
        CDISTS, ITERATES, key, orientation_preserving=False, fixed_points=[fixed_point]
    )
    homotopy, iterated = _families(trellis, key)
    landing = land_element(trellis, homotopy, iterated, L(key, 1))
    assert landing.image_side == "right"
    assert landing.target == R(key, 1) and landing.contained
    preserving = FakeTrellis(CDISTS, ITERATES, key, fixed_points=[fixed_point])
    assert land_element(preserving, *_families(preserving, key), L(key, 1)).target == L(key, 1)


def test_zero_width_span_probes_its_low_end():
    fixed_point = FakeFixedPoint(beta=0.25)
    key = stable_key(fixed_point)
    cdists = {**CDISTS, 9: 1.0, 19: 0.25}
    trellis = FakeTrellis(cdists, {**ITERATES, 9: 19}, key, fixed_points=[fixed_point])
    # A hand-made non-singleton element whose two ends coincide in cdist.
    degenerate = [(None, 1, 0.0, 1.0, True, False), (1, 9, 1.0, 1.0, True, True)]
    homotopy, iterated = _families(trellis, key, homotopy_specs=degenerate)
    landing = land_element(trellis, homotopy, iterated, L(key, 1))
    assert landing.span == (0.25, 0.25)
    assert not landing.singleton
    assert landing.target == L(key, 1) and landing.contained


def test_singleton_lands_on_the_registered_iterate():
    fixed_point = FakeFixedPoint(beta=0.25)
    key = stable_key(fixed_point)
    trellis = FakeTrellis(CDISTS, ITERATES, key, fixed_points=[fixed_point])
    homotopy_specs = [
        (None, 1, 0.0, 1.0, True, False),
        (1, 2, 1.0, 2.0, True, False),
        (2, 2, 2.0, 2.0, True, True),
        (2, 3, 2.0, 4.0, False, True),
    ]
    iterated_specs = [
        (None, 11, 0.0, 0.25, True, False),
        (11, 12, 0.25, 0.5, True, False),
        (12, 12, 0.5, 0.5, True, True),
        (12, 1, 0.5, 1.0, False, False),
        (1, 2, 1.0, 2.0, True, False),
        (2, 2, 2.0, 2.0, True, True),
        (2, 3, 2.0, 4.0, False, True),
    ]
    homotopy, iterated = _families(trellis, key, homotopy_specs, iterated_specs)
    landing = land_element(trellis, homotopy, iterated, L(key, 2))
    assert landing.singleton and landing.resolved and landing.contained
    assert landing.target == L(key, 2)
    assert landing.span == (0.5, 0.5) and landing.from_table == (True, True)
    assert landing.covering == (L(key, 2),)
    assert "singleton" in repr(landing)


def test_singleton_with_unregistered_image_is_unresolved(caplog):
    fixed_point = FakeFixedPoint(beta=0.25)
    key = stable_key(fixed_point)
    trellis = FakeTrellis(CDISTS, {1: 11, 3: 1}, key, fixed_points=[fixed_point])  # no 2 -> 12
    homotopy_specs = [
        (None, 1, 0.0, 1.0, True, False),
        (1, 2, 1.0, 2.0, True, False),
        (2, 2, 2.0, 2.0, True, True),
        (2, 3, 2.0, 4.0, False, True),
    ]
    homotopy, iterated = _families(trellis, key, homotopy_specs)
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        landing = land_element(trellis, homotopy, iterated, L(key, 2))
    assert landing.singleton and not landing.resolved
    assert landing.target is None and landing.covering == ()
    assert "not registered" in landing.reason
    assert landing.span == (0.5, 0.5) and landing.from_table == (False, False)
    assert any("not registered" in record.getMessage() for record in _warnings(caplog))
    assert "unresolved" in repr(landing)


def test_singleton_owned_by_a_non_singleton_is_info_not_warning(caplog):
    """A singleton whose image is not pinched is the backward-only-holes rule at
    work (the forward hole is never punched), so it is reported at INFO and the
    landing is resolved normally."""
    fixed_point = FakeFixedPoint(beta=0.25)
    key = stable_key(fixed_point)
    trellis = FakeTrellis(CDISTS, ITERATES, key, fixed_points=[fixed_point])
    homotopy_specs = [
        (None, 1, 0.0, 1.0, True, False),
        (1, 2, 1.0, 2.0, True, False),
        (2, 2, 2.0, 2.0, True, True),
        (2, 3, 2.0, 4.0, False, True),
    ]
    # The iterated partition does not pinch 12: element #1 owns it.
    homotopy, iterated = _families(trellis, key, homotopy_specs, ITERATED)
    with caplog.at_level(logging.INFO, logger=LOGGER):
        landing = land_element(trellis, homotopy, iterated, L(key, 2))
    assert landing.singleton and landing.resolved and landing.contained
    assert landing.target == L(key, 1)
    assert not _warnings(caplog)
    infos = [r for r in caplog.records if r.levelno == logging.INFO]
    assert any(
        "NON-singleton" in r.getMessage() and "backward only" in r.getMessage() for r in infos
    )


def test_span_straddling_two_elements_fills_covering(caplog):
    fixed_point = FakeFixedPoint(beta=0.25)
    key = stable_key(fixed_point)
    trellis = FakeTrellis(CDISTS, ITERATES, key, fixed_points=[fixed_point])
    finer = [
        (None, 11, 0.0, 0.25, True, False),
        (11, 14, 0.25, 0.375, True, True),
        (14, 12, 0.375, 0.5, False, True),
        (12, 1, 0.5, 1.0, False, False),
        (1, 2, 1.0, 2.0, True, True),
        (2, 3, 2.0, 4.0, False, True),
    ]
    homotopy, iterated = _families(trellis, key, iterated_specs=finer)
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        landing = land_element(trellis, homotopy, iterated, L(key, 1))
    assert landing.resolved and not landing.contained
    assert landing.span == (0.25, 0.5)
    assert landing.target == L(key, 1)  # the midpoint 0.375 is owned by the closed element
    assert landing.covering == (L(key, 1), L(key, 2))
    assert landing.reason is None
    assert any("straddles" in record.getMessage() for record in _warnings(caplog))
    assert "straddles" in repr(landing)


def test_image_on_an_unexpected_branch_is_unresolved(caplog):
    fixed_point = FakeFixedPoint(beta=0.25)
    key = stable_key(fixed_point)
    other = stable_key(FakeFixedPoint(period=2))
    keys = {iid: key for iid in CDISTS}
    keys[12] = other  # the image of crossing 2 sits on a foreign branch
    trellis = FakeTrellis(CDISTS, ITERATES, keys, fixed_points=[fixed_point])
    homotopy, iterated = _families(trellis, key)
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        landing = land_element(trellis, homotopy, iterated, L(key, 1))
    assert not landing.resolved
    assert "instead of" in landing.reason
    assert any("instead of" in record.getMessage() for record in _warnings(caplog))


def test_probe_outside_the_iterated_partition_is_unresolved(caplog):
    fixed_point = FakeFixedPoint(beta=0.25)
    key = stable_key(fixed_point)
    trellis = FakeTrellis(CDISTS, ITERATES, key, fixed_points=[fixed_point])
    short = ITERATED[:2]  # the iterated partition stops at 0.5
    homotopy, iterated = _families(trellis, key, iterated_specs=short)
    with caplog.at_level(logging.WARNING, logger=LOGGER):
        landing = land_element(trellis, homotopy, iterated, L(key, 2))
    assert not landing.resolved and landing.target is None
    assert "no single iterated element" in landing.reason
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

    def __init__(self, key, missing=()):
        self.key = key
        self.missing = set(missing)

    def owner_of_intersection(self, crossing_id, side):
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
    with pytest.raises(ValueError, match="at least one"):
        trellis_itinerary(trellis, _OwnerFamily(key), [], +1)
    with pytest.raises(ValueError, match="not consecutive"):
        trellis_itinerary(trellis, _OwnerFamily(key), [(0, 1), (2, 3)], +1)
    with pytest.raises(ValueError, match="grow or blast"):
        trellis_itinerary(trellis, _OwnerFamily(key, missing={2}), chain, +1)
    unsigned = FakeTrellis({0: 0.0, 1: 1.0}, {}, key, signs={0: 0, 1: -1})
    with pytest.raises(ValueError, match="no crossing sign"):
        trellis_itinerary(unsigned, _OwnerFamily(key), [(0, 1)], +1)


# --------------------------------------------------------------------------- #
# k=10 fixture
# --------------------------------------------------------------------------- #
def test_k10_landings_and_the_active_class_walk(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    dual = DualGraph(pieces.minimal, pieces.iterated, strong_pips=pieces.strong_pips)
    trellis = pieces.full

    landings = {}
    for result in pieces.homotopy:
        for interval in result.intervals:
            element = result.ref(interval.element_id)
            landing = land_element(trellis, pieces.homotopy, pieces.iterated, element)
            assert landing.resolved, landing.reason
            assert landing.contained and not landing.singleton
            assert landing.image_key == element.branch_key
            assert landing.image_side == element.side  # Henon with b = 1 preserves orientation
            landings[element] = landing
    # The anchor-side element of each side lands in the anchor-side iterated element.
    for result in pieces.homotopy:
        anchor = landings[result.ref(0)]
        assert anchor.target.element_id == 0 and anchor.span[0] == 0.0

    active = pieces.table.active
    assert len(active) == 1
    entry = active[0]
    x, y = entry.bridge_class.source, entry.bridge_class.target
    search = shortest_walks(dual, landings[x].target, landings[y].target)
    assert search.status == "unique", search
    assert not search.truncated
    walk = search.walk
    assert len(walk.itinerary) == 6
    assert walk.itinerary[0] == landings[x].target and walk.itinerary[-1] == landings[y].target
    for first, second in _pairs(walk.itinerary):
        assert first.side == second.side
    assert all(step.node.is_unified for step in walk.steps)

    # Every registered member image reads the same route off the regular
    # trellis, up to the homotopy parents of each pair.
    def parents(itinerary):
        return [
            pieces.iterated.element(element).parent_element_id for element in itinerary
        ]

    checked = 0
    for member in entry.members:
        if member.is_loop:
            continue
        chain = _image_chain(trellis, member.bridge_id)
        if chain is None:
            continue
        itinerary = trellis_itinerary(trellis, pieces.iterated, chain, member.direction)
        assert len(itinerary) % 2 == 0
        assert len(itinerary) == len(walk.itinerary)
        assert parents(itinerary) == parents(walk.itinerary)
        checked += 1
    assert checked > 0

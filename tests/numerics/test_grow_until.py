"""Phase 7 -- the predicate-driven growth drivers.

``grow_until`` is the generic loop the two named wrappers are built on: grow one
iteration on each requested stability, recompute the crossings keeping their ids,
re-infer the iterate table, stop when the caller's predicate says the tangle is
big enough.

What is pinned here:

* the predicate is consulted BEFORE the first growth, so an already-satisfied
  request costs nothing (a driver that always grew once would silently double
  the size of every tangle it was pointed at);
* the cap raises ``ValueError`` like every other ``grow_until_*`` sibling, and
  it raises only after having actually grown;
* ``grow`` selects which manifolds move, and the wrappers pick it from the
  direction they are closing (a forward image lives further out on the UNSTABLE
  branch, a backward image further out on the STABLE one);
* ``ids="all"`` freezes the id set at call time, so the crossings born during
  the loop do not move the goalposts;
* crossing ids survive the loop -- the recompute runs with ``preserve_ids=True``,
  which is the only reason a caller may hold an id across a growth round at all.

Note:
    ``grow_until_faces_closed`` is pinned on its two reachable outcomes rather
    than on an open-to-closed transition. The unbounded face of a tangle touches
    the outermost crossings permanently: growing a branch pushes its tip further
    out but the crossings already on the outer envelope stay on it, and every
    experiment on the k=10 fixture (four growth patterns, both stabilities) shows
    the open-face corner set only ever GROWING, with each round's new crossings
    born already interior. The cap test uses the anchor, which never closes for a
    different and more repairable reason: the arrangement registers one anchor per
    (unstable branch, stable branch) pair and has no notion of a branch continuing
    THROUGH a node, so a periodic point's ``u-``/``s-`` slots get virtual stubs
    instead of the other branch's ``u+``/``s+``. That is the deferred item in the
    "Inversion caveat" Dev Note of ``tanglepack.topology.Arrangement``, not a
    consequence of what an anchor is.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack.topology.Trellis import Trellis


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _point_counts(workbench) -> dict:
    """Point count of every registered manifold, keyed by manifold key."""
    return {
        key: len(manifold.get_point_array())
        for key, manifold in workbench.manifolds.items()
    }


def _non_anchor_ids(workbench) -> list[int]:
    """Registry ids of the detected (non-synthetic) crossings."""
    return sorted(
        iid for iid, ix in workbench.intersection_registry if not ix.is_synthetic
    )


def _interior_ids(workbench) -> list[int]:
    """Crossings whose every incident face in the all-fixed-points arrangement is closed."""
    arrangement = Trellis.from_workbench(workbench, None).arrangement
    seen: set[int] = set()
    open_corners: set[int] = set()
    for face in arrangement.faces:
        for corner in face.corners:
            seen.add(corner)
            if not face.is_closed:
                open_corners.add(corner)
    return sorted(seen - open_corners)


# --------------------------------------------------------------------------- #
# grow_until -- the generic loop
# --------------------------------------------------------------------------- #
def test_grow_until_checks_the_predicate_before_growing(small_tangle):
    """A predicate that is already true must not move a single point."""
    workbench, fp = small_tangle
    before = _point_counts(workbench)

    rounds = workbench.grow_until(fp, lambda wb: True, max_iterations=5)

    assert rounds == 0
    assert _point_counts(workbench) == before


def test_grow_until_raises_at_the_cap(small_tangle):
    """An impossible predicate exhausts the cap and raises, like the siblings."""
    workbench, fp = small_tangle
    before = _point_counts(workbench)

    with pytest.raises(ValueError, match="[Mm]ax iterations"):
        workbench.grow_until(fp, lambda wb: False, max_iterations=1)

    after = _point_counts(workbench)
    assert any(after[key] > before[key] for key in before), (
        "the cap must be reached by growing, not by refusing to start"
    )


def test_grow_until_grows_only_the_requested_stabilities(small_tangle):
    """``grow=('unstable',)`` leaves the stable manifold exactly where it was."""
    workbench, fp = small_tangle
    before = _point_counts(workbench)

    with pytest.raises(ValueError):
        workbench.grow_until(
            fp, lambda wb: False, grow=("unstable",), max_iterations=1
        )

    after = _point_counts(workbench)
    for key in before:
        if key[1] == "stable":
            assert after[key] == before[key]
    assert after[(fp, "unstable", 0, 0)] > before[(fp, "unstable", 0, 0)]


def test_grow_until_rejects_an_empty_grow_set(small_tangle):
    workbench, fp = small_tangle
    with pytest.raises(ValueError, match="at least one"):
        workbench.grow_until(fp, lambda wb: False, grow=())


def test_grow_until_requires_an_initialized_manifold(fixed_point):
    workbench, fp = fixed_point  # no manifolds seeded yet
    with pytest.raises(ValueError, match="has not been initialized"):
        workbench.grow_until(fp, lambda wb: True)


def test_grow_until_preserves_crossing_ids_across_rounds(small_tangle):
    """The loop recomputes with ``preserve_ids=True``: held ids stay valid."""
    workbench, fp = small_tangle
    held = _non_anchor_ids(workbench)
    assert held, "the fixture must already have a crossing to hold an id of"
    coords = {iid: workbench.intersection_registry[iid].coords for iid in held}

    with pytest.raises(ValueError):
        workbench.grow_until(
            fp, lambda wb: False, grow=("unstable",), max_iterations=2
        )

    for iid, before in coords.items():
        assert iid in workbench.intersection_registry
        after = workbench.intersection_registry[iid].coords
        assert np.allclose(after, before, atol=1e-6), (
            f"crossing {iid} was renumbered onto a different point"
        )


# --------------------------------------------------------------------------- #
# grow_until_iterates_closed
# --------------------------------------------------------------------------- #
def test_iterates_closed_forward_terminates(small_tangle):
    """The forward image of the fixture's only crossing appears after one round."""
    workbench, fp = small_tangle
    targets = _non_anchor_ids(workbench)
    table = workbench.intersection_registry.iterate_table
    assert not any((iid, 1) in table for iid in targets), (
        "the fixture must start with no forward image, or the driver is untested"
    )

    rounds = workbench.grow_until_iterates_closed(
        fp, ids=targets, direction="forward", max_iterations=4
    )

    assert rounds >= 1
    table = workbench.intersection_registry.iterate_table
    assert all((iid, 1) in table for iid in targets)


def test_iterates_closed_backward_terminates(small_tangle):
    """A backward image lives further out on the STABLE branch."""
    workbench, fp = small_tangle
    targets = _non_anchor_ids(workbench)
    before = _point_counts(workbench)

    rounds = workbench.grow_until_iterates_closed(
        fp, ids=targets, direction="backward", max_iterations=4
    )

    assert rounds >= 1
    table = workbench.intersection_registry.iterate_table
    assert all((iid, -1) in table for iid in targets)
    # backward closes on the stable side: the unstable manifold must not move
    after = _point_counts(workbench)
    assert after[(fp, "unstable", 0, 0)] == before[(fp, "unstable", 0, 0)]
    assert after[(fp, "stable", 0, 0)] > before[(fp, "stable", 0, 0)]


def test_iterates_closed_all_freezes_the_id_set(small_tangle):
    """Crossings born during the loop do not move the goalposts."""
    workbench, fp = small_tangle
    frozen = _non_anchor_ids(workbench)

    rounds = workbench.grow_until_iterates_closed(fp, max_iterations=4)

    assert rounds >= 1
    table = workbench.intersection_registry.iterate_table
    assert all((iid, 1) in table for iid in frozen)
    born = [iid for iid in _non_anchor_ids(workbench) if iid not in frozen]
    assert born, "the round must have produced new crossings for this to mean anything"
    assert any((iid, 1) not in table for iid in born), (
        "the driver stopped only because the id set was frozen at call time"
    )


def test_iterates_closed_is_a_no_op_when_already_closed(small_tangle):
    """Asking again immediately grows nothing."""
    workbench, fp = small_tangle
    workbench.grow_until_iterates_closed(fp, max_iterations=4)
    frozen = [
        iid
        for iid in _non_anchor_ids(workbench)
        if (iid, 1) in workbench.intersection_registry.iterate_table
    ]
    before = _point_counts(workbench)

    assert workbench.grow_until_iterates_closed(fp, ids=frozen) == 0
    assert _point_counts(workbench) == before


def test_iterates_closed_rejects_an_unknown_id(small_tangle):
    workbench, fp = small_tangle
    unknown = max(iid for iid, _ix in workbench.intersection_registry) + 100
    with pytest.raises(ValueError, match="not in the registry"):
        workbench.grow_until_iterates_closed(fp, ids=[unknown])


def test_iterates_closed_rejects_an_unknown_direction(small_tangle):
    workbench, fp = small_tangle
    with pytest.raises(ValueError, match="direction"):
        workbench.grow_until_iterates_closed(fp, direction="sideways")


# --------------------------------------------------------------------------- #
# grow_until_faces_closed
# --------------------------------------------------------------------------- #
def test_faces_closed_returns_immediately_for_interior_crossings(small_tangle):
    """Crossings already surrounded by closed faces cost no growth at all."""
    workbench, fp = small_tangle
    workbench.grow_n_times(fp, "unstable", num_iterations=2)
    workbench.compute_intersections([fp], preserve_ids=True)
    interior = _interior_ids(workbench)
    assert interior, "the fixture must have an interior crossing to ask about"
    before = _point_counts(workbench)

    assert workbench.grow_until_faces_closed(fp, ids=interior) == 0
    assert _point_counts(workbench) == before


def test_faces_closed_raises_at_the_cap_for_the_anchor(small_tangle):
    """An anchor cannot close under the arrangement's current anchor model.

    A periodic point is topologically ONE degree-four node whose ``u-`` and ``s-``
    rays are the other branch's ``u+`` and ``s+``. The arrangement does not model
    that: ``_register_anchors`` registers one anchor per (unstable branch, stable
    branch) pair and nothing joins them, so those two slots are filled with
    virtual stubs and every face through one of them is open, however far the
    manifolds are grown. Until that deferred item is closed (the "Inversion
    caveat" Dev Note in ``tanglepack.topology.Arrangement``) an anchor is exactly
    the request the cap exists to stop.
    """
    workbench, fp = small_tangle
    anchors = [iid for iid, ix in workbench.intersection_registry if ix.is_synthetic]
    assert anchors, "compute_intersections must have registered the anchor"
    assert anchors[0] not in _interior_ids(workbench)

    with pytest.raises(ValueError, match="[Mm]ax iterations"):
        workbench.grow_until_faces_closed(fp, ids=anchors, max_iterations=2)


def test_faces_closed_rejects_an_empty_id_set(small_tangle):
    workbench, fp = small_tangle
    with pytest.raises(ValueError, match="at least one"):
        workbench.grow_until_faces_closed(fp, ids=[])


# --------------------------------------------------------------------------- #
# session exposure
# --------------------------------------------------------------------------- #
def test_session_exposes_the_drivers(henon_map, henon_map_inverse):
    """The drivers are session methods, and the trellis cache follows along."""
    from tanglepack import TangleSession

    session = TangleSession(henon_map, henon_map_inverse)
    fp = session.construct_fixed_point([4, -4])
    session.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=7)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])

    stale = session.trellis(fp)
    rounds = session.grow_until_iterates_closed(fp, max_iterations=4)

    assert rounds >= 1
    assert session.trellis(fp) is not stale, (
        "growth moved the generation; the cached trellis must have been rebuilt"
    )
    assert session.grow_until(fp, lambda wb: True) == 0
    # reachable as session methods (delegated to the workbench, not re-declared)
    for name in ("grow_until", "grow_until_iterates_closed", "grow_until_faces_closed"):
        assert callable(getattr(session, name))


# --------------------------------------------------------------------------- #
# bridges stay consistent with the grown curves
# --------------------------------------------------------------------------- #
def test_grow_until_recuts_existing_bridges(small_tangle):
    """Growth inserts crossings between old bridge endpoints; the cut must follow.

    An arrangement asserts that every bridge is a pair of CONSECUTIVE crossings on
    an unstable branch. A round of growth can drop a new crossing inside an
    already-cut bridge, so a driver that left the bridge set alone would leave the
    workbench in a state where the topological layer refuses to build.
    """
    workbench, fp = small_tangle
    workbench.grow_n_times(fp, "unstable", num_iterations=2)
    workbench.compute_intersections([fp])
    workbench.create_bridges(fp)
    assert workbench.bridges, "the fixture must start with bridges to recut"

    # the STABLE side is the one that does it: growing the stable manifold drops
    # new crossings BETWEEN the endpoints of bridges already cut on the unstable
    # branch, while growing the unstable manifold only extends past the last one.
    with pytest.raises(ValueError, match="[Mm]ax iterations"):
        workbench.grow_until(fp, lambda wb: False, grow=("stable",), max_iterations=1)

    # builds at all == every bridge is still consecutive on its branch
    arrangement = Trellis.from_workbench(workbench, None).arrangement
    registry = workbench.intersection_registry
    for bridge in workbench.bridges:
        assert bridge.id is not None
        assert all(endpoint in registry for endpoint in bridge.id)
    assert arrangement.faces


# --------------------------------------------------------------------------- #
# "all" is scoped to the fixed point being grown
# --------------------------------------------------------------------------- #
def _two_saddle_workbench(henon_map, henon_map_inverse, henon_jacobian):
    """The k=10 map's two saddles, each with a small tangle, in one workbench.

    The second is the inversion saddle at (-2.3166, 2.3166) that
    ``tests/conftest.py`` documents -- the cheapest second tangle available on the
    fixture map. Grown just far enough that each saddle owns a crossing the other
    does not.
    """
    from tanglepack import TangleWorkbench

    workbench = TangleWorkbench(henon_map, henon_map_inverse, henon_jacobian)

    simple = workbench.construct_fixed_point([4, -4])
    workbench.orient_eigenvectors(
        simple, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    workbench.initialize_both_manifolds(simple)
    workbench.grow_n_times(simple, "unstable", num_iterations=7)
    workbench.grow_until_turnaround(simple, "stable")

    other = workbench.construct_fixed_point([-2.3166, 2.3166])
    assert other.check_inversion()
    workbench.initialize_both_manifolds(other)
    workbench.grow_n_times(other, "unstable", num_iterations=5)
    workbench.grow_n_times(other, "stable", num_iterations=5)

    workbench.compute_intersections([simple, other])
    return workbench, simple, other


def _ids_of(workbench, fixed_point) -> set:
    """The non-anchor registry ids of the crossings involving one fixed point."""
    registry = workbench.intersection_registry
    id_of = {id(ix): iid for iid, ix in registry}
    return {
        id_of[id(ix)]
        for ix in registry.from_fixed_point(fixed_point)
        if not ix.is_synthetic
    }


def test_iterates_closed_all_is_scoped_to_the_grown_fixed_point(
    henon_map, henon_map_inverse, henon_jacobian
):
    """``ids="all"`` must not hold another tangle's crossings hostage.

    Only ``fixed_point``'s manifolds grow, so a set resolved against the whole
    REGISTRY would contain crossings this call can never close -- and would run to
    the cap by construction in any session with more than one saddle.
    """
    workbench, simple, other = _two_saddle_workbench(
        henon_map, henon_map_inverse, henon_jacobian
    )
    mine = _ids_of(workbench, simple)
    theirs = _ids_of(workbench, other)
    every = set(_non_anchor_ids(workbench))
    assert mine and theirs - mine, (
        "each saddle must own a crossing the other does not for this to bite"
    )

    frozen = set(workbench._frozen_ids("all", fixed_point=simple, allow_all=True))

    assert frozen == mine
    assert not (frozen & (theirs - mine))
    assert frozen != every, "the whole-registry resolution is the bug being pinned"

    # and the excluded ids really are unreachable from this call: growing only
    # `simple` leaves them without a forward image, so the loop would have capped.
    workbench.grow_n_times(simple, "unstable", num_iterations=1)
    workbench.compute_intersections([simple, other], preserve_ids=True)
    table = workbench.intersection_registry.iterate_table
    assert any((iid, 1) not in table for iid in theirs - mine)


# --------------------------------------------------------------------------- #
# the per-round recut does not corrupt blast state
# --------------------------------------------------------------------------- #
def test_grow_until_clears_the_iterated_flag_it_invalidates(small_tangle, caplog):
    """A recut that drops the images must not leave the sources marked iterated.

    ``rebuild_bridges`` carries ``iterated`` across by id, but the recut throws
    away the image bridges and partial pieces that flag was recording. A bridge
    left marked ``iterated`` with no derivable image is invisible to
    ``uniiterated_bridges``, so a later blast silently skips it forever.
    """
    import logging

    workbench, fp = small_tangle
    workbench.grow_n_times(fp, "unstable", num_iterations=2)
    workbench.compute_intersections([fp])
    workbench.create_bridges(fp)
    workbench.iterate_all_bridges()
    assert any(bridge.iterated for bridge in workbench.bridges), (
        "the fixture must have iterated something for the repair to be tested"
    )

    with caplog.at_level(logging.WARNING, logger="tanglepack.numerics.TangleWorkbench"):
        with pytest.raises(ValueError, match="[Mm]ax iterations"):
            workbench.grow_until(
                fp, lambda wb: False, grow=("stable",), max_iterations=1
            )

    stranded = [
        bid
        for bid, bridge in workbench._bridges.items()
        if bridge.iterated and workbench.image_bridges(bid) is None
    ]
    assert not stranded, (
        f"bridges {stranded} are marked iterated but their image cannot be "
        "derived: a later blast will never re-iterate them"
    )
    assert any("recut the bridge set" in record.message for record in caplog.records), (
        "losing a blast frontier must be logged, not silent"
    )

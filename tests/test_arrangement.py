"""Phase 6.4/6.5/6.6 -- the planar arrangement and its regions.

**The arrangement** (6.4) is checked on a hand-built two-crossing fixture whose
face count can be worked out by hand and by Euler's formula; on the real
fixtures the faces that merely ENCLOSE another component are pinned as
non-regions, open faces refuse geometric questions, and the lookups agree with
the regions.

**The regions** (6.5, 6.6): a region whose mapped corner cycle is carried by a
proper sub-face of its true image is rejected by area preservation, a region's
preimage is checked against a hand-derived backward image, and a blast (arcs
on iterated bridges) still builds.

The crossing-sign laws and the per-case arrangement laws (Euler per component,
regions sound and disjoint, ``image_of`` against the dynamics, preimage of the
image) run on every law case in ``tests/invariants/`` (``test_law_crossings``,
``test_law_anchors``, ``test_law_arrangement``).
"""

from __future__ import annotations

import numpy as np
import pytest

from helpers.fakes import bare_fixed_point
from tanglepack.numerics.Bridge import Bridge
from tanglepack.numerics.Intersection import Intersection
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry
from tanglepack.numerics.Point import Point
from tanglepack.topology.Arrangement import Arrangement
from tanglepack.topology.TopologyResults import Arc, canonical_corners
from tanglepack.topology.Trellis import Trellis
from tanglepack.topology.TrellisBranch import TrellisBranch


# --------------------------------------------------------------------------- #
# 6.4 -- a hand-built arrangement whose faces can be counted by hand
# --------------------------------------------------------------------------- #
def _hand_built_trellis():
    """One stable line crossed twice by one unstable lobe.

    Two crossings ``a`` and ``b``; the stable branch runs a -> b and the unstable
    branch runs a -> b as well (the lobe leaves the stable line at a and returns at
    b). Signs are opposite, as the alternation test above establishes they must be.
    That is 2 real nodes and 2 real edges, plus 4 dangling ends; Euler gives
    ``F = E - V + 2 = 6 - 6 + 2 = 2`` faces, of which exactly one -- the lobe -- is
    bounded entirely by manifold.
    """
    fp = bare_fixed_point(beta=0.5, coordinates=[(0.0, 0.0)])
    u_key = (fp, "unstable", 0, 0)
    s_key = (fp, "stable", 0, 0)

    registry = IntersectionRegistry()
    a = registry.add(
        Intersection.synthetic(
            coords=(0.0, 0.0),
            unstable_cdist=0.0,
            stable_cdist=0.0,
            manifold_a_key=u_key,
            manifold_b_key=s_key,
            crossing_sign=1,
        )
    )
    b = registry.add(
        Intersection.synthetic(
            coords=(0.0, 2.0),
            unstable_cdist=5.0,
            stable_cdist=2.0,
            manifold_a_key=u_key,
            manifold_b_key=s_key,
            crossing_sign=-1,
        )
    )

    # A real Bridge for (a, b): the lobe, bulging out to the right of the line.
    root = Point(0.0, -0.2, -0.5)
    mid = Point(3.0, 1.0, 2.5)
    tail = Point(0.0, 2.2, 5.5)
    root.forward, mid.backward = mid, root
    mid.forward, tail.backward = tail, mid
    bridge = Bridge(
        root=root,
        stability="unstable",
        stretch_param=1.0,
        fixed_point=fp,
        tail=tail,
        branch_index=0,
        manifold_key=u_key,
        first_intersection=a,
        second_intersection=b,
    )

    branches = {
        u_key: TrellisBranch(
            key=u_key,
            fixed_point=fp,
            stability="unstable",
            orbit_index=0,
            branch_index=0,
            intersection_ids=[a, b],
        ),
        s_key: TrellisBranch(
            key=s_key,
            fixed_point=fp,
            stability="stable",
            orbit_index=0,
            branch_index=0,
            intersection_ids=[a, b],
        ),
    }
    trellis = Trellis(
        fixed_points=[fp],
        registry=registry,
        branches=branches,
        bridges=[bridge],
    )
    return trellis, a, b


def test_hand_built_arrangement_has_exactly_one_closed_face():
    trellis, a, b = _hand_built_trellis()
    arrangement = Arrangement.from_trellis(trellis)

    assert len(arrangement.faces) == 2
    assert len(arrangement.regions) == 1
    assert len(arrangement.open_faces) == 1

    lobe = arrangement.regions[0]
    assert lobe.is_closed
    assert lobe.corners == (a, b)
    assert {arc.kind for arc in lobe.arcs} == {"stable", "unstable"}
    assert lobe.bridge_ids == [(a, b)]


def test_hand_built_region_geometry_closes_and_contains_its_lobe():
    trellis, _a, _b = _hand_built_trellis()
    lobe = Arrangement.from_trellis(trellis).regions[0]

    ring = lobe.boundary_points
    assert len(ring) >= 3
    # The bridge bulges to x = 3 and the stable line sits on x = 0.
    assert lobe.contains((1.0, 1.0))
    assert not lobe.contains((-1.0, 1.0))
    assert lobe.verify_representative_point()
    # Bounded faces are traversed clockwise (see the Arrangement Dev Notes).
    assert lobe.area < 0


def test_hand_built_lookups():
    trellis, a, b = _hand_built_trellis()
    arrangement = Arrangement.from_trellis(trellis)
    lobe = arrangement.regions[0]

    assert arrangement.region((b, a)) is lobe  # any rotation/direction
    assert arrangement.regions_at(a) == [lobe]
    assert arrangement.regions_at(b) == [lobe]
    for arc in lobe.arcs:
        assert arrangement.regions_bounded_by(arc) == [lobe]
        # The face on the other side of every arc is the open one.
        assert lobe.neighbor_across(arc) is None
    assert arrangement.region((a,)) is None


def test_every_bridge_must_be_consecutive_on_its_branch():
    """The consistency assertion fires on a bridge that is not an arc."""
    trellis, a, b = _hand_built_trellis()
    stray = trellis.bridges[0]
    stray.second_intersection = b + 99
    with pytest.raises(AssertionError):
        Arrangement.from_trellis(trellis)


def test_an_arc_without_a_bridge_is_still_an_arc():
    """The converse does NOT hold: a blast leaves arcs with no bridge object.

    Dropping the bridge list leaves the branch orderings untouched, so the lobe is
    still an arc and still closes a face -- it simply carries no
    :data:`~tanglepack.numerics.Bridge.BridgeId`. With no bridge AND no parent
    manifold to fall back on, its geometry collapses to the chord between the two
    crossings, which is the documented degradation.
    """
    trellis, a, b = _hand_built_trellis()
    trellis.bridges = []
    arrangement = Arrangement.from_trellis(trellis)

    assert len(arrangement.regions) == 1
    lobe = arrangement.regions[0]
    assert lobe.corners == (a, b)
    assert lobe.bridge_ids == []
    assert np.allclose(lobe.boundary_points, [[0.0, 0.0], [0.0, 2.0]])
    assert lobe.area == 0.0


def test_canonical_corners_preserves_cyclic_order():
    assert canonical_corners([5, 2, 9]) == (2, 9, 5)
    assert canonical_corners([]) == ()


def test_arc_reverse_round_trips():
    arc = Arc(kind="stable", lo_id=1, hi_id=4, branch_key=("k",), bridge_id=None)
    assert arc.tail_id == 1 and arc.head_id == 4
    back = arc.reversed()
    assert back.tail_id == 4 and back.head_id == 1
    assert back.edge_key == arc.edge_key
    assert back.reversed() == arc


# --------------------------------------------------------------------------- #
# 6.4/6.5 -- the real fixtures
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_p3_has_a_closed_face_that_swallows_the_inner_tangle(henon_p3_session):
    """The outer tangle's enclosing face is closed but is NOT a region.

    It is a perfectly good cycle of crossings, and it geometrically contains the
    whole other component, so by the definition of a region it is not one. It is
    kept on ``containing_faces`` rather than discarded, and every region really
    does sit inside it.
    """
    session, _fp3, _fp1, _zone = henon_p3_session
    arrangement = session.arrangement()

    assert arrangement.containing_faces, (
        "the nested fixture must produce at least one enclosing face"
    )
    for face in arrangement.containing_faces:
        assert face.is_closed and not face.is_minimal
        assert face not in arrangement.regions
        swallowed = [
            region
            for region in arrangement.regions
            if region.representative_point is not None
            and face.contains(region.representative_point)
        ]
        assert len(swallowed) > 1, (
            "an enclosing face must contain more than one region, or it would be "
            "a region itself"
        )


def test_open_faces_refuse_to_answer_geometric_questions(k10_session):
    session, _fp = k10_session
    arrangement = session.arrangement()
    assert arrangement.open_faces
    face = arrangement.open_faces[0]
    for attempt in (
        lambda: face.boundary_points,
        lambda: face.area,
        lambda: face.representative_point,
        lambda: face.contains((0.0, 0.0)),
    ):
        with pytest.raises(ValueError):
            attempt()


def test_regions_at_and_bounded_by_agree_with_the_regions(k10_session):
    session, _fp = k10_session
    arrangement = session.arrangement()

    for region in arrangement.regions:
        for corner in region.corners:
            assert region in arrangement.regions_at(corner)
        for arc in region.arcs:
            bounding = arrangement.regions_bounded_by(arc)
            assert region in bounding
            assert len(bounding) <= 2
            neighbor = region.neighbor_across(arc)
            if neighbor is not None:
                assert neighbor is not region
                assert region in [
                    other
                    for other in arrangement.regions_bounded_by(arc)
                    if other is not neighbor
                ]


def test_stable_arcs_report_their_partition_elements(k10_session):
    session, fp = k10_session
    session.classify_strong_pips()
    session.compute_pseudoneighbors(fp)
    session.punch_holes(fp)
    session.partition_stable_manifold(fp)

    trellis = session.trellis(fp)
    arrangement = session.arrangement()
    assert trellis.stable_partitions, "the partition must have produced a result"
    labelled = 0
    for region in arrangement.regions:
        for arc in region.stable_arcs:
            for element in arc.elements(trellis):
                assert element.branch_key == arc.branch_key
                assert element.side in ("left", "right")
                labelled += 1
        for arc in region.arcs:
            if arc.kind == "unstable":
                assert arc.elements(trellis) == []
    assert labelled > 0, "no stable arc reported a partition element"


# --------------------------------------------------------------------------- #
# 6.6 -- the combinatorial image agrees with the dynamics
# --------------------------------------------------------------------------- #
def test_k10_rejects_an_image_that_is_only_a_sub_face(k10_session):
    """A face carrying the mapped corner cycle need not BE the image.

    On k=10 exactly one region's mapped corners are carried by a closed face, and
    that face is a proper sub-face of the true image: an extra stable arc runs
    across the image's interior and splits it, and the mapped corners sit on the
    larger piece. Area preservation catches it — the piece is smaller than the
    source — and the union of the piece and its neighbour across the splitting arc
    recovers the source area exactly. ``image_of`` must return None rather than the
    piece.
    """
    session, _fp = k10_session
    arrangement = session.arrangement()
    trellis = arrangement.trellis

    carriers = 0
    for region in arrangement.regions:
        images = [trellis.iterate(corner, 1) for corner in region.corners]
        if any(image is None for image in images):
            continue
        faces = [
            face
            for face in arrangement.regions
            if all(corner in face.corners for corner in images)
        ]
        if len(faces) != 1:
            continue
        carrier = faces[0]
        carriers += 1

        assert abs(carrier.area) < abs(region.area), (
            "this test is about a carrier that is too SMALL to be the image"
        )
        assert arrangement.image_of(region, 1) is None

        # The missing area is one neighbouring face: the image is the union.
        missing = abs(region.area) - abs(carrier.area)
        assert any(
            abs(other.area) == pytest.approx(missing, rel=2e-2)
            for other in arrangement.regions
            if other is not carrier
        ), "the sub-face plus one neighbour should recover the source area"
    assert carriers >= 1, "the fixture must exercise the sub-face rejection"


@pytest.mark.slow
def test_preimage_of_a_region_is_its_hand_derived_backward_image(henon_p3_session):
    """A worked example: take one region, map its corners BACK by hand, look up.

    Not a restatement of ``preimage_of``'s implementation -- the expected answer is
    built here from the registry's iterate table and the plain corner lookup, and
    then checked to be a genuine preimage of the region under the inverse map by
    carrying its representative point forward one step.
    """
    session, _fp3, _fp1, _zone = henon_p3_session
    arrangement = session.arrangement()
    trellis = arrangement.trellis
    forward = session.workbench.dynamical_system.map

    checked = 0
    for region in arrangement.regions:
        expected_corners = [trellis.iterate(corner, -1) for corner in region.corners]
        if any(corner is None for corner in expected_corners):
            continue
        expected = arrangement.region(expected_corners)
        if expected is None:
            continue
        if abs(abs(expected.area) - abs(region.area)) > 2e-2 * abs(region.area):
            continue  # a sub-face, rejected by image_of -- covered by its own test

        assert arrangement.preimage_of(region, 1) is expected
        # ...and it really is a preimage: map it forward for real.
        mapped = np.asarray(
            forward(np.asarray(expected.representative_point, dtype=float)),
            dtype=float,
        )
        assert region.contains(mapped)
        checked += 1
    assert checked > 0, "no region had a hand-derivable preimage"


def test_arrangement_survives_a_blast(k10_session):
    """A blast adds crossings on ITERATED bridges, whose arcs were never cut.

    Those arcs have no ``BridgeId``; the arrangement must still build (the
    every-bridge-is-an-arc invariant is untouched) and still close faces. This is
    the case that made the two-way bridge/arc equality assertion untenable.
    """
    session, fp = k10_session
    trellis = session.trellis(fp)
    trellis.classify_strong_pips()
    session.add_resonance_zones([trellis.strong_pip])
    (zone,) = session.resonance_zones.values()
    session.blast_zone(zone, 2, min_separation=1e-3)

    arrangement = session.arrangement()
    assert arrangement.regions
    assert any(
        arc.bridge_id is None
        for region in arrangement.regions
        for arc in region.arcs
        if arc.kind == "unstable"
    ), "the blast must leave at least one unstable arc with no bridge object"

    representable = 0
    for region in arrangement.regions:
        if len(region.boundary_points) < 3:
            # A blast puts crossings closer together than the manifold has nodes,
            # so a face can end up with both of its arcs collapsed to the same
            # chord. Such a face is real combinatorially and unrepresentable
            # geometrically at this resolution; it has no representative point,
            # and that is the honest answer rather than a fabricated one.
            assert region.representative_point is None
            assert region.area == 0.0
            continue
        assert region.verify_representative_point(), (
            f"region {region.corners} has a representative point outside itself"
        )
        representable += 1
    assert representable > 0, "the blast closed no representable face"

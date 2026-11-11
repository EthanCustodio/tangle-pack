"""
Geometric self-diagnostics of the hole and partition algorithm: flags, never decisions.

The crossing signs decide every row and every opening, and a hole's side is
decided once at its reference hole and carried along its orbit
(``StablePartition``). This module re-reads the same quantities off the
polylines, where a coarse polyline can disagree with them:

* **Rows.** :func:`bridge_row_diagnosis` reads the side of the stable branch a
  bridge approaches each end from on its polyline. On a 2-4 node bridge the
  "approach" is the vector to the far node, so a disagreement with the
  sign is an under-resolved bridge.
* **Hole sides.** A propagated hole's side MEASURED on the curves: its orbit's
  reference hole carried back by the real inverse map and read on the image
  sub-arc, which escapes the lobe when the reference lobe is under-resolved.

* **Words.** A class's word is its dual-graph walk; a member's registered
  image chain, read with the same translation and reduction, only cross-checks
  it (:func:`verify_words`). A disagreement is diagnosed (:class:`EvidenceDiagnosis`):
  ``"stable end"`` evidence is incomplete and ``"order flip"`` evidence noise,
  neither gives a verdict; ``"walk"`` counts against ``is_reliable``.

:func:`geometry_flags` collects rows and sides and :func:`warn_geometry_flags` logs one
WARNING per kind (``Trellis.punch_holes(diagnose=True)``). The two invariants
the algorithm asserts, I1 (``check_holes_share_bridge_side``) and I2
(``check_bridge_rows_consistent``), stay in ``StablePartition`` and only borrow
their messages from here (:func:`bridge_side_violations`,
:func:`bridge_row_violation`), so nothing here feeds a decision, mirroring
``numerics/diagnostics.py``.
"""

from __future__ import annotations

import logging
from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from itertools import groupby
from typing import Iterable, Literal, Optional, Sequence, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ..numerics.geometry import oriented_bridge_polyline
from .BridgeClass import _image_chain
from .DualWalk import trellis_itinerary
from .StablePartition import bridge_for_pair, row_of_end, rows_of_bridge
from .TopologyResults import ElementRef, Hole, Side, side_after

if TYPE_CHECKING:
    from .Trellis import Trellis
    from ..numerics.Bridge import Bridge, BridgeId
    from ..numerics.Intersection import Intersection, ManifoldKey
    from .BridgeClass import BridgeClass, BridgeClassEntry, BridgeClassTable
    from .ElementNaming import ElementNaming
    from .PartitionFamily import PartitionFamily
    from .SymbolicDynamics import ClassDynamics, Symbol

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: Where a hole's ``bridge_side`` came from (:func:`side_source`).
SideSource = Literal["signs", "carried"]


def _side_of(tangent: NDArray[np.float64], displacement: NDArray[np.float64]) -> Optional[Side]:
    """Side of the displacement relative to the tangent, standard orientation.

    Standing on the manifold facing along ``tangent`` (toward the anchor),
    positive cross = left hand — the author's dynamical-direction convention.
    """
    cross = float(tangent[0] * displacement[1] - tangent[1] * displacement[0])
    if cross > 0.0:
        return "left"
    if cross < 0.0:
        return "right"
    return None


def _anchorward_look(
    trellis: "Trellis", ix: "Intersection"
) -> Optional[NDArray[np.float64]]:
    """Unit stable direction toward the anchor at a crossing, read off the manifold.

    Toward the nearest node below the crossing's stable cdist, else away from
    the nearest one above (node cdists are monotone along a branch, so both are
    a bisect); None without manifold nodes.
    """
    manifold = trellis.manifolds.get(ix.manifold_b_key)
    if manifold is None:
        return None
    points, cdists = manifold.points, manifold.cdists
    here = np.asarray(ix.get_point(), dtype=np.float64)
    below = bisect_left(cdists, ix.stable_cdist) - 1
    above = bisect_right(cdists, ix.stable_cdist)
    if below >= 0:
        look = points[below] - here
    elif above < len(points):
        look = here - points[above]
    else:
        return None
    norm = float(np.linalg.norm(look))
    return look / norm if norm > 0.0 else None


def _approach_displacement(
    poly: NDArray[np.float64], q: NDArray[np.float64]
) -> Optional[NDArray[np.float64]]:
    """How the arc approaches its end nearest ``q``: the vector from ``q`` into the bridge.

    The end node sits just past the crossing, so the vector runs to the second
    node clearly distinct from ``q`` (the first when there is only one). On a
    two-node polyline that is the FAR end, which is why this reading only
    flags.
    """
    eps = 1e-9 * (1.0 + float(np.linalg.norm(poly[-1] - poly[0])))
    if np.linalg.norm(poly[0] - q) > np.linalg.norm(poly[-1] - q):
        poly = poly[::-1]
    distinct = np.flatnonzero(np.linalg.norm(poly[1:] - q, axis=1) > eps)[:2]
    return poly[1 + distinct[-1]] - q if len(distinct) else None


def _oriented_bridge_polyline(
    trellis: "Trellis", bridge: "Bridge"
) -> Optional[NDArray[np.float64]]:
    """The bridge polyline oriented by the unstable dynamical direction.

    A bridge is unstable manifold, so its dynamical direction — the direction
    of forward flow — points away from the fixed point, i.e. along increasing
    unstable canonical distance. Storage order (root→tail) need not match;
    this is the single place that mismatch is normalized: the polyline is
    reversed when the first endpoint's unstable cdist exceeds the second's.
    Endpoint cdists always exist, even for keyless iterated-bridge children,
    so no manifold key is consulted.

    Args:
        trellis: The Trellis resolving the bridge's endpoint intersections.
        bridge: The bridge whose polyline to orient.

    Returns:
        The (N, 2) polyline in dynamical orientation, or None when the bridge
        has fewer than two points, an unresolved endpoint, or endpoints of
        equal unstable cdist (orientation undecidable — logged as a warning).

    Note:
        Memoised per trellis and per bridge version
        (:meth:`Trellis.cached_bridge_polyline`); the result is shared, so it
        must not be written to.
    """
    return trellis.cached_bridge_polyline(
        bridge, lambda: _build_oriented_bridge_polyline(trellis, bridge)
    )


def _build_oriented_bridge_polyline(
    trellis: "Trellis", bridge: "Bridge"
) -> Optional[NDArray[np.float64]]:
    """The uncached body of :func:`_oriented_bridge_polyline`.

    Resolves the bridge's two endpoint canonical distances against the trellis
    and hands the orientation itself to
    :func:`~tanglepack.numerics.geometry.oriented_bridge_polyline`, which is the
    single implementation shared with the region layer.
    """
    if bridge.first_intersection is None or bridge.second_intersection is None:
        return None
    return oriented_bridge_polyline(
        bridge,
        trellis.intersection(bridge.first_intersection).unstable_cdist,
        trellis.intersection(bridge.second_intersection).unstable_cdist,
        tol=trellis.registry.cdist_tol,
    )


def _bridge_side_of(
    trellis: "Trellis",
    bridge: "Bridge",
    point: NDArray[np.float64],
    span: Optional[tuple[float, float]] = None,
) -> Optional[Side]:
    """Which side of ``bridge`` the point sits on, in dynamical orientation.

    Standing on the bridge looking along the unstable dynamical direction
    (away from the fixed point, increasing unstable cdist), positive
    cross(tangent, displacement) = left — the same convention as
    :func:`_side_of` on the stable manifold. The side is read at the nearest
    SEGMENT (the foot of the perpendicular, clipped to the segment), so a
    sub-arc lying inside one segment of a coarse bridge still answers.

    Args:
        trellis: The Trellis resolving the bridge's endpoint intersections.
        bridge: The bridge to classify against.
        point: Phase-space point to classify.
        span: Optional unstable-cdist extent of the sub-arc to read the side
            against (a propagated hole's image sub-arc); the whole bridge
            when omitted.

    Returns:
        ``"left"`` or ``"right"``, or None when the bridge cannot be oriented
        or the point sits on the arc.
    """
    poly = _oriented_bridge_polyline(trellis, bridge)
    if poly is None:
        return None
    if span is not None:
        # The image sub-arc only: a far fold of the same bridge can sit nearer
        # the point than the sub-arc it belongs to.
        cdists = np.sort(bridge.cdists)
        lo = max(bisect_left(cdists, span[0]) - 1, 0)
        poly = poly[lo : bisect_right(cdists, span[1]) + 1]
    point = np.asarray(point, dtype=np.float64)
    start, step = poly[:-1], np.diff(poly, axis=0)
    along = np.einsum("ij,ij->i", point - start, step)
    length2 = np.einsum("ij,ij->i", step, step)
    t = np.clip(along / np.maximum(length2, np.finfo(np.float64).tiny), 0.0, 1.0)
    foot = start + t[:, None] * step
    nearest = int(np.argmin(np.linalg.norm(point - foot, axis=1)))
    return _side_of(step[nearest], point - foot[nearest])


@dataclass(frozen=True)
class RowEnd:
    """
    One end of a bridge, read by its crossing sign and by its polyline.

    Attributes:
        intersection_id: Registry id of the crossing.
        crossing_sign: Its ``crossing_sign``.
        combinatorial_row: :func:`row_of_end` (None at a zero sign), the row
            that decides.
        geometric_row: The side of the stable branch, looking toward the
            anchor, the polyline approaches the crossing from (None at an
            anchor, which has no stable frame of its own, or when undecidable).
        approach_angle: Degrees between the arc's approach displacement and
            the stable tangent, in ``[0, 90]`` (None when undecidable). A
            fraction of a degree is a near-tangent end.
        displacement: Length of that approach displacement (None when
            undecidable).
        unstable_cdist: The crossing's unstable canonical distance.
        stable_cdist: The crossing's stable canonical distance.
    """

    intersection_id: int
    crossing_sign: float
    combinatorial_row: Optional[Side]
    geometric_row: Optional[Side]
    approach_angle: Optional[float]
    displacement: Optional[float]
    unstable_cdist: float
    stable_cdist: float

    def describe(self) -> str:
        """One line naming every reading of this end."""
        angle = "-" if self.approach_angle is None else f"{self.approach_angle:.3g} deg"
        disp = "-" if self.displacement is None else f"{self.displacement:.3g}"
        return (
            f"crossing {self.intersection_id}: sign {self.crossing_sign:+g}, "
            f"combinatorial row {self.combinatorial_row}, geometric row "
            f"{self.geometric_row}, approach {angle} over |d| {disp}, cdists "
            f"(u {self.unstable_cdist!r}, s {self.stable_cdist!r})"
        )


@dataclass(frozen=True)
class RowDiagnosis:
    """
    Every reading of a bridge's two rows, for the row assertion and the flag.

    Attributes:
        bridge_id: The two crossing ids, ordered by unstable cdist (the
            :data:`~tanglepack.numerics.Bridge.BridgeId` order).
        branch_key: The stable branch holding the first crossing (None when
            it lies on no branch of the trellis).
        ends: The ``"first"`` and ``"second"`` ends, in ``bridge_id`` order.
        n_nodes: Nodes in the bridge polyline the geometric rows are read on.
        chord: Distance between the two crossings.
    """

    bridge_id: tuple[int, int]
    branch_key: Optional["ManifoldKey"]
    ends: tuple[RowEnd, RowEnd]
    n_nodes: int
    chord: float

    @property
    def signs_alternate(self) -> bool:
        """True when the two crossing signs differ (the combinatorial rows then agree)."""
        return self.ends[0].crossing_sign * self.ends[1].crossing_sign < 0

    @property
    def geometry_disagrees(self) -> bool:
        """Some end's polyline reads a row other than its crossing sign: under-resolved."""
        return any(
            end.geometric_row is not None and end.geometric_row != end.combinatorial_row
            for end in self.ends
        )

    @property
    def verdict(self) -> str:
        """Which reading to distrust."""
        if self.signs_alternate:
            angles = [e.approach_angle for e in self.ends if e.approach_angle is not None]
            angle = f"{min(angles):.3g} deg" if angles else "an undecidable angle"
            return (
                "signs alternate: the combinatorial rows agree; the geometric "
                f"reading on a {self.n_nodes}-node polyline at {angle} is likely "
                "below resolution"
            )
        return "signs equal: the combinatorial rows disagree as well"

    def describe(self) -> str:
        """The diagnosis as text: the bridge, both ends, and the verdict."""
        key = self.branch_key
        branch = "None" if key is None else f"({key[1]}, {key[2]}, {key[3]})"
        return "\n".join(
            [
                f"bridge {self.bridge_id} on stable branch {branch}: "
                f"{self.n_nodes} node(s), chord {self.chord:.3g}",
                *(f"  {end.describe()}" for end in self.ends),
                f"  verdict: {self.verdict}",
            ]
        )


def bridge_row_diagnosis(
    trellis: "Trellis", bridge: "Bridge"
) -> Optional[RowDiagnosis]:
    """
    Read both ends of a bridge by their crossing signs and by the polylines.

    The geometric row of an end is the side of the stable branch, looking
    toward the anchor (:func:`_anchorward_look`), of the arc's approach
    displacement there (:func:`_approach_displacement`). It is the
    cross-check of :func:`row_of_end`, never the decision.

    Args:
        trellis: The Trellis carrying the intersections, branches and manifolds.
        bridge: The bridge to diagnose.

    Returns:
        The :class:`RowDiagnosis`, or None for a partial piece or a bridge
        that cannot be oriented.
    """
    poly = _oriented_bridge_polyline(trellis, bridge)
    if poly is None or bridge.id is None:
        return None
    ends = []
    for endpoint, iid in zip(("first", "second"), bridge.id):
        ix = trellis.intersection(iid)
        disp = look = None
        if not ix.is_anchor:
            disp = _approach_displacement(poly, np.asarray(ix.get_point(), dtype=np.float64))
            look = None if disp is None else _anchorward_look(trellis, ix)
        geometric_row = angle = displacement = None
        if look is not None:
            geometric_row = _side_of(look, disp)
            displacement = float(np.linalg.norm(disp))
            angle = float(np.degrees(np.arctan2(
                abs(float(look[0] * disp[1] - look[1] * disp[0])),
                abs(float(np.dot(look, disp))),
            )))
        ends.append(RowEnd(
            intersection_id=iid,
            crossing_sign=float(ix.crossing_sign),
            combinatorial_row=(
                row_of_end(trellis, bridge.id, endpoint) if ix.crossing_sign else None
            ),
            geometric_row=geometric_row,
            approach_angle=angle,
            displacement=displacement,
            unstable_cdist=ix.unstable_cdist,
            stable_cdist=ix.stable_cdist,
        ))
    branch = trellis.branch_containing(bridge.id[0], "stable")
    first, second = (
        np.asarray(trellis.intersection(i).get_point(), dtype=np.float64) for i in bridge.id
    )
    return RowDiagnosis(
        bridge_id=bridge.id,
        branch_key=None if branch is None else branch.key,
        ends=tuple(ends),
        n_nodes=len(poly),
        chord=float(np.linalg.norm(second - first)),
    )


def bridge_row_violation(
    trellis: "Trellis", bridge: "Bridge"
) -> Optional[str]:
    """
    Describe how a bridge breaks the row invariant, or None if it holds.

    A bridge whose two defining crossings lie on the SAME stable branch must
    approach both of them from the same side of that branch: the lobe lies
    entirely on one side of the stable manifold, because the endpoints are
    consecutive crossings and the arc between them therefore never crosses it
    again. The rows are the crossing signs' (:func:`rows_of_bridge`), so
    the invariant fails exactly when the two signs are EQUAL. Bridges whose
    crossings sit on different stable branches (heteroclinic, or a period > 1
    orbit) are exempt, and so is a partial piece. A polyline that reads other
    rows while the signs alternate is under-resolved; :func:`geometry_flags`
    reports it.

    Args:
        trellis: The Trellis carrying the intersections, branches and manifolds.
        bridge: The bridge to check.

    Returns:
        A human-readable violation message followed by the bridge's
        :class:`RowDiagnosis` (built only on failure), or None when the
        invariant holds or the check does not apply.
    """
    if bridge.id is None:
        return None
    first, second = bridge.id
    unstable_key = trellis.unstable_key_of(bridge.id)
    branch = trellis.branch(trellis.stable_key_at(first, unstable_key))
    if branch is None or branch is not trellis.branch(trellis.stable_key_at(second, unstable_key)):
        return None  # different stable branches (heteroclinic / period > 1)
    row_first, row_second = rows_of_bridge(trellis, bridge.id)
    if row_first == row_second:
        return None
    diagnosis = bridge_row_diagnosis(trellis, bridge)
    return (
        f"Bridge {bridge.id} on stable branch (orbit {branch.orbit_index}, branch "
        f"{branch.branch_index}) approaches its two defining crossings from opposite "
        f"sides of that branch (equal crossing signs): {first} -> {row_first}, "
        f"{second} -> {row_second}. A bridge's endpoints are consecutive crossings, "
        f"so the arc between them cannot cross the stable manifold again."
        + ("" if diagnosis is None else f"\n{diagnosis.describe()}")
    )


def side_source(hole: Hole) -> SideSource:
    """A direct hole reads its side off its crossing signs; a propagated one carries it."""
    return "signs" if hole.pair is not None else "carried"


def measured_sides(
    trellis: "Trellis", holes: Iterable[Hole]
) -> dict[tuple[tuple[int, int], int], Side]:
    """
    Each propagated hole's side MEASURED on the curves: the flag on its carried side.

    Replays the backward chain of ``propagate_reference_holes`` on the curves:
    the orbit's reference hole is carried back by the real inverse map, one
    step at a time, and read at the nearest segment of the hole bridge's
    image sub-arc (the reference bridge's crossings moved by
    :meth:`Trellis.image_cdist`). The carried point escapes the lobe when the
    reference lobe is under-resolved, so a side other than
    ``bridge_side`` is a flag, never an error.

    Args:
        trellis: The Trellis the holes were punched on.
        holes: The holes to measure, typically ``trellis.holes``.

    Returns:
        ``{(origin, iterate): side}`` for every propagated hole whose orbit has
        a reference hole among ``holes``; ``(origin, iterate)`` is unique per
        propagation. Empty without an inverse map; a hole whose carried point
        sits on the arc has no entry.

    Note:
        It calls the map, so it runs on the flag path only
        (``Trellis.punch_holes(diagnose=True)``), never on an assertion. The
        reference point is the one figures draw (:func:`plotting.hole_point`).
    """
    from .plotting import hole_point

    if trellis.dynamical_system is None:
        return {}
    holes = list(holes)
    references = {h.origin: h for h in holes if h.pair is not None and h.pair.is_reference}
    sides = {}
    for hole in holes:
        reference = references.get(hole.origin)
        if hole.pair is not None or reference is None:
            continue
        steps = -hole.iterate
        carried = hole_point(trellis, reference)
        for _ in range(steps):
            carried = np.asarray(
                trellis.dynamical_system.map_inv(carried), dtype=np.float64
            ).ravel()
        bridge = bridge_for_pair(trellis, reference.pair)
        span = tuple(sorted(
            trellis.image_cdist(end_id, -steps, "unstable", branch=bridge.manifold_key)[1]
            for end_id in bridge.id
        ))
        side = _bridge_side_of(trellis, trellis.bridge_between(*hole.bounding_ids), carried, span)
        if side is not None:
            sides[(hole.origin, hole.iterate)] = side
    return sides


def _orbit_side(hole: Hole, orientation_preserving: bool) -> Optional[Side]:
    """A hole's bridge side carried back to iterate 0 of its orbit (None without
    a side, or without an iterate to carry it by under a reversing map)."""
    if hole.bridge_side is None or (hole.iterate is None and not orientation_preserving):
        return None
    return side_after(hole.bridge_side, -(hole.iterate or 0), orientation_preserving)


@dataclass(frozen=True)
class SideDiagnosis:
    """
    One hole that fails a hole-side check, with every reading of its side.

    Attributes:
        bounding_ids: The hole's bridge (its two defining crossings).
        origin: The hole's orbit (``Hole.origin``).
        iterate: The hole's place on it (``Hole.iterate``).
        side_source: Where its ``bridge_side`` came from (:func:`side_source`).
        bridge_side: The stored side.
        measured_side: The side measured on the curves (:func:`measured_sides`;
            None when not measured).
        expected_side: The orbit's side at this iterate, from the first hole of
            the origin (parity-flipped under an orientation-reversing map);
            None when the orbit check does not apply.
        reference_iterate: The iterate of that first hole.
    """

    bounding_ids: Optional[tuple[int, int]]
    origin: Optional[tuple[int, int]]
    iterate: Optional[int]
    side_source: Optional[SideSource]
    bridge_side: Optional[Side]
    measured_side: Optional[Side]
    expected_side: Optional[Side]
    reference_iterate: Optional[int]

    @property
    def orbit_disagrees(self) -> bool:
        """The stored side breaks the orbit invariant."""
        return self.expected_side is not None and self.bridge_side != self.expected_side

    @property
    def measured_disagrees(self) -> bool:
        """The measured side contradicts the stored (carried) one."""
        return self.measured_side is not None and self.measured_side != self.bridge_side

    def describe(self) -> str:
        """One line naming the hole and each reading that disagrees."""
        parts = [
            f"hole {self.bounding_ids} of origin {self.origin} at iterate "
            f"{self.iterate} ({self.side_source}): bridge_side {self.bridge_side}"
        ]
        if self.orbit_disagrees:
            parts.append(
                f"its orbit expects {self.expected_side} (from iterate "
                f"{self.reference_iterate})"
            )
        if self.measured_disagrees:
            parts.append(f"measured {self.measured_side} on its image sub-arc")
        return ", ".join(parts)


def side_diagnoses(
    holes: Iterable[Hole],
    *,
    orientation_preserving: bool = True,
    measured: Optional[dict[tuple[tuple[int, int], int], Side]] = None,
) -> list[SideDiagnosis]:
    """
    Diagnose every hole that disagrees with its orbit or with its measured side.

    The ORBIT check: a propagated hole is the backward image of its orbit's
    reference hole, and the map carries a bridge's dynamical orientation to
    the image bridge's, so an orientation-PRESERVING map (det J > 0 — both
    Hénon fixtures have b = 1) leaves the left/right classification fixed
    along the whole backward chain. Under an orientation-REVERSING map the
    side alternates with the parity of :attr:`Hole.iterate`, which
    ``orientation_preserving=False`` encodes. Since propagated holes CARRY
    their side (author's rule), this guards the direct holes of one
    orbit (its references, each read off its own crossing signs) and the
    carrying itself.

    The MEASURED reading is the geometric cross-check: a hole whose measured
    side (:func:`measured_sides`: the reference hole carried back by the real
    inverse map, read on the image sub-arc) contradicts its stored side is
    under-resolved. It is a flag (:func:`geometry_flags`), never an error.

    Holes with no ``origin`` or no comparable side skip the orbit check, and
    propagated holes with no entry in ``measured`` skip the measured one; both
    counts are logged at INFO (the second only when ``measured`` is given).

    Args:
        holes: The holes to check, typically ``trellis.holes``.
        orientation_preserving: True when the map has det J > 0.
        measured: The measured sides by ``(origin, iterate)``
            (:func:`measured_sides`); None for the orbit check alone, which
            never calls the map.

    Returns:
        One :class:`SideDiagnosis` per disagreeing hole, in input order.
    """
    holes = list(holes)
    orbit_side: dict[tuple[int, int], tuple[Side, Optional[int]]] = {}
    unplaced = 0
    for hole in holes:
        side = _orbit_side(hole, orientation_preserving)
        if hole.origin is None or side is None:
            unplaced += 1
        else:
            orbit_side.setdefault(hole.origin, (side, hole.iterate))
    unmeasured = 0 if measured is None else sum(
        1
        for hole in holes
        if side_source(hole) == "carried" and (hole.origin, hole.iterate) not in measured
    )
    if unplaced or unmeasured:
        logger.info(
            "Hole-side check skips %d hole(s) with no comparable orbit side and "
            "%d propagated hole(s) with no measured side", unplaced, unmeasured
        )

    diagnoses: list[SideDiagnosis] = []
    for hole in holes:
        expected = reference_iterate = None
        if hole.origin in orbit_side and _orbit_side(hole, orientation_preserving):
            expected, reference_iterate = orbit_side[hole.origin]
            expected = side_after(expected, hole.iterate or 0, orientation_preserving)
        diagnosis = SideDiagnosis(
            bounding_ids=hole.bounding_ids,
            origin=hole.origin,
            iterate=hole.iterate,
            side_source=side_source(hole),
            bridge_side=hole.bridge_side,
            measured_side=None if measured is None else measured.get((hole.origin, hole.iterate)),
            expected_side=expected,
            reference_iterate=reference_iterate,
        )
        if diagnosis.orbit_disagrees or diagnosis.measured_disagrees:
            diagnoses.append(diagnosis)
    return diagnoses


def bridge_side_violations(
    holes: Iterable[Hole], *, orientation_preserving: bool = True
) -> list[str]:
    """
    Describe every hole whose side of its bridge breaks the orbit invariant (I1).

    See :func:`side_diagnoses` for the check and which holes it skips; a
    measured side that disagrees is a flag, not a violation.

    Args:
        holes: The holes to check, typically ``trellis.holes``.
        orientation_preserving: True when the map has det J > 0.

    Returns:
        One line per hole whose orbit disagrees (its :class:`SideDiagnosis`).
        Empty when I1 holds.
    """
    return [
        diagnosis.describe()
        for diagnosis in side_diagnoses(holes, orientation_preserving=orientation_preserving)
        if diagnosis.orbit_disagrees
    ]


def geometry_flags(
    trellis: "Trellis", bridges: Iterable["Bridge"], holes: Iterable[Hole]
) -> tuple[list[RowDiagnosis], list[SideDiagnosis]]:
    """
    The geometric readings that disagree with the sign rule.

    The crossing signs decide rows and openings and the orbit carries hole
    sides; the polylines only cross-check them. A bridge whose polyline reads
    another row than its signs at some end, or a propagated hole whose
    carried point lies on the other side of its image sub-arc, marks an
    under-resolved polyline (or a carried point that left its lobe): a reason
    to refine, never an error.

    Args:
        trellis: The Trellis carrying the bridges and manifolds.
        bridges: The bridges to read, typically the hole bridges.
        holes: The holes whose measured sides (:func:`measured_sides`) to
            compare.

    Returns:
        ``(rows, sides)``: the :class:`RowDiagnosis` of every bridge whose
        ``geometry_disagrees``, and the :class:`SideDiagnosis` of every hole
        whose ``measured_disagrees``.
    """
    rows = [
        diagnosis
        for diagnosis in (bridge_row_diagnosis(trellis, bridge) for bridge in bridges)
        if diagnosis is not None and diagnosis.geometry_disagrees
    ]
    holes = list(holes)
    sides = [
        diagnosis
        for diagnosis in side_diagnoses(
            holes,
            orientation_preserving=trellis.orientation_preserving,
            measured=measured_sides(trellis, holes),
        )
        if diagnosis.measured_disagrees
    ]
    return rows, sides


def warn_geometry_flags(
    trellis: "Trellis", bridges: Iterable["Bridge"], holes: Iterable[Hole]
) -> None:
    """One WARNING per kind of :func:`geometry_flags` (a resolution flag)."""
    rows, sides = geometry_flags(trellis, bridges, holes)
    if rows:
        logger.warning(
            "%d hole bridge(s) read a geometric row other than their crossing "
            "signs: under-resolved polylines, flag only (D6); %s\n%s",
            len(rows),
            [diagnosis.bridge_id for diagnosis in rows[:8]],
            rows[0].describe(),
        )
    if sides:
        logger.warning(
            "%d propagated hole(s) measured on the other side of their bridge than "
            "their orbit carries: under-resolved lobes, flag only (D6); %s",
            len(sides),
            "; ".join(diagnosis.describe() for diagnosis in sides[:5]),
        )


# ── the words against the registered images ─────────────────────────────────


@dataclass(frozen=True)
class EvidenceDiagnosis:
    """
    Why a registered image disagrees with the walked word.

    Attributes:
        position: The first token where the two words differ.
        walked: The walked token's occurrence there (iterated pair, class
            source end first), ``None`` past the end of the walked word.
        registered: The registered token's occurrence there, ``None`` past the
            end of the registered word.
        cause: ``"stable end"`` when two consecutive crossings of the
            registered chain on one stable branch share a crossing sign: the
            image passed around the end of the computed stable branch, so its
            crossings with the trimmed-away part are not registered and the
            registered word is incomplete. ``"order flip"`` when a member
            end and a bound of its own element have registered images in the
            reverse stable order: the map preserves the order along a branch,
            so the registered cdists contradict the dynamics there and the
            evidence is noise. ``"walk"`` otherwise: no numerical
            explanation, the shortest walk is not the image's class.
        crossings: The equal-sign run behind ``"stable end"``, the
            ``(end, bound)`` behind ``"order flip"``, else ``()``.
    """

    position: int
    walked: Optional[tuple[ElementRef, ElementRef]]
    registered: Optional[tuple[ElementRef, ElementRef]]
    cause: Literal["stable end", "order flip", "walk"]
    crossings: tuple[int, ...] = ()


@dataclass
class ImageEvidence:
    """
    The registered forward image of one member bridge, as a word.

    Attributes:
        bridge_id: The member.
        direction: Its direction in its class.
        itinerary: The iterated elements along its image chain, in the same
            even-length, disjoint-pair format as a dual-graph walk, running
            from the image of the class source to the image of the target.
        symbols: The itinerary translated (loops dropped), freely reduced.
        chain: The registered image chain the itinerary was read off.
        word_agrees: True when the non-loop ``(class, direction)`` sequence
            equals the class's primary word, False when it does not; ``None``
            when the class has no itinerary to compare with, or when the
            evidence is incomplete or noise (a ``"stable end"`` or
            ``"order flip"`` diagnosis).
        itinerary_agrees: True when the iterated refs agree exactly.
        diagnosis: Why the word disagrees, when it does.
    """

    bridge_id: "BridgeId"
    direction: int
    itinerary: tuple[ElementRef, ...]
    symbols: list["Symbol"]
    chain: tuple["BridgeId", ...] = ()
    word_agrees: Optional[bool] = None
    itinerary_agrees: Optional[bool] = None
    diagnosis: Optional[EvidenceDiagnosis] = None

    @property
    def word(self) -> str:
        """The registered word as text."""
        return " ".join(symbol.text for symbol in self.symbols)


def registered_evidence(
    trellis: "Trellis",
    iterated: "PartitionFamily",
    naming: "ElementNaming",
    table: "BridgeClassTable",
    letters: dict["BridgeClass", str],
    entry: "BridgeClassEntry",
) -> list[ImageEvidence]:
    """
    The registered forward images of a class's members, as words.

    Args:
        trellis: The full trellis (iterate table, crossings).
        iterated: The iterated partition family owning the image crossings.
        naming: Resolves iterated refs to their homotopy parents.
        table: The bridge classes.
        letters: The letter map (see :func:`translate_itinerary`).
        entry: The class whose members to read.

    Returns:
        One :class:`ImageEvidence` per non-loop member with a registered image
        chain, in member order, its word freely reduced like the class word,
        with ``word_agrees`` / ``itinerary_agrees`` left ``None``
        (:func:`verify_words` fills them). Loop members and members without a
        registered chain contribute nothing.
    """
    from .SymbolicDynamics import _signature, chain_pieces, free_reduce, translate_itinerary

    evidence: list[ImageEvidence] = []
    for member in entry.members:
        if member.is_loop:
            continue
        chain = _image_chain(trellis, member.bridge_id)
        if not chain:
            continue
        try:
            itinerary = tuple(trellis_itinerary(trellis, iterated, chain, member.direction))
        except ValueError as error:
            logger.warning(
                "registered image of member %s of class %s cannot be read as an "
                "itinerary; it contributes no evidence (%s)",
                member.bridge_id,
                entry.letter if entry.letter is not None else entry.bridge_class.label,
                error,
            )
            continue
        symbols, _loops = translate_itinerary(
            itinerary, naming, table, letters, trellis.fixed_points,
            chain_pieces(table, chain, member.direction) if any(c.bridge_class.lobe for c in table) else None,
        )
        word = [symbols[i] for i in free_reduce(_signature(symbols))]
        evidence.append(
            ImageEvidence(member.bridge_id, member.direction, itinerary, word, tuple(chain))
        )
    return evidence


def _diagnose(
    trellis: "Trellis",
    homotopy: "PartitionFamily",
    cls: "BridgeClass",
    evidence: ImageEvidence,
    symbols: Sequence["Symbol"],
    position: int,
) -> EvidenceDiagnosis:
    """Why a registered word parts from the walked one at ``position``."""
    walked = symbols[position].occurrence if position < len(symbols) else None
    registered = (
        evidence.symbols[position].occurrence if position < len(evidence.symbols) else None
    )
    chain = evidence.chain
    crossings = [chain[0][0], *(pair[1] for pair in chain)] if chain else []
    # The first two or more consecutive crossings on one stable branch with one sign.
    def mark(crossing: int) -> tuple:
        ix = trellis.intersection(crossing)
        return ix.manifold_b_key, ix.crossing_sign

    runs = (tuple(run) for _, run in groupby(crossings, key=mark))
    run = next((run for run in runs if len(run) > 1), None)
    if run:
        return EvidenceDiagnosis(position, walked, registered, "stable end", run)
    flip = _order_flip(trellis, homotopy, cls, evidence)
    if flip:
        return EvidenceDiagnosis(position, walked, registered, "order flip", flip)
    return EvidenceDiagnosis(position, walked, registered, "walk")


def _order_flip(
    trellis: "Trellis", homotopy: "PartitionFamily", cls: "BridgeClass", evidence: ImageEvidence
) -> Optional[tuple[int, int]]:
    """A member end and a bound of its own element whose registered images run in
    the reverse stable order: the map preserves the order along a branch, so
    the registered cdists there are noise. Lookup only."""
    ends = evidence.bridge_id if evidence.direction >= 0 else evidence.bridge_id[::-1]
    for end, ref in zip(ends, (cls.source, cls.target)):
        element = homotopy.element(ref)
        for bound in (element.lo_id, element.hi_id):
            if bound is None or bound == end:
                continue
            images = trellis.iterate(end, 1), trellis.iterate(bound, 1)
            if None in images:
                continue
            x, y = (trellis.intersection(i) for i in (end, bound))
            fx, fy = (trellis.intersection(i) for i in images)
            if x.manifold_b_key == y.manifold_b_key and fx.manifold_b_key == fy.manifold_b_key and (
                (x.stable_cdist < y.stable_cdist) != (fx.stable_cdist < fy.stable_cdist)
            ):
                return end, bound
    return None


def verify_words(
    records: Iterable["ClassDynamics"],
    trellis: "Trellis",
    naming: "ElementNaming",
    table: "BridgeClassTable",
    letters: dict["BridgeClass", str],
) -> None:
    """
    Check every class's word against the registered images of its members.

    Runs once every word is final, never changes one, and fills each record's
    ``evidence``, ``verified`` and ``itinerary_verified``.

    Args:
        records: The classes, their ``symbols`` and ``itinerary`` final.
        trellis: The full trellis (iterate table, crossings).
        naming: The element naming (its iterated family owns the image
            crossings, its homotopy family the members' elements).
        table: The closed alphabet.
        letters: ``{class: letter}`` over the table.
    """
    from .SymbolicDynamics import _signature

    iterated = naming.iterated
    for cd in records:
        cd.evidence = registered_evidence(trellis, iterated, naming, table, letters, cd.entry)
        if not cd.evidence or cd.itinerary is None:
            continue
        signature = _signature(cd.symbols)
        for evidence in cd.evidence:
            evidence.itinerary_agrees = tuple(evidence.itinerary) == tuple(cd.itinerary)
            registered = _signature(evidence.symbols)
            if registered == signature:
                evidence.word_agrees = True
                continue
            position = next(
                (i for i, (w, r) in enumerate(zip(signature, registered)) if w != r),
                min(len(signature), len(registered)),
            )
            evidence.diagnosis = diagnosis = _diagnose(
                trellis, naming.homotopy, cd.entry.bridge_class, evidence, cd.symbols, position
            )
            incomplete = diagnosis.cause in ("stable end", "order flip")
            evidence.word_agrees = None if incomplete else False
            if diagnosis.cause == "stable end":
                why = (
                    f"incomplete: the image passes the end of the computed stable branch at "
                    f"crossings {list(diagnosis.crossings)} (same branch, same sign), so it "
                    "gives no verdict"
                )
            elif diagnosis.cause == "order flip":
                end, bound = diagnosis.crossings
                gaps = [
                    trellis.intersection(a).stable_cdist - trellis.intersection(b).stable_cdist
                    for a, b in ((end, bound), (trellis.iterate(end, 1), trellis.iterate(bound, 1)))
                ]
                why = (
                    f"noise: member end {end} and its element's bound {bound} are "
                    f"{gaps[0]:+.2e} apart in stable cdist but their registered images "
                    f"{gaps[1]:+.2e}, an order the map cannot reverse, so it gives no verdict"
                )
            else:
                why = "no numerical explanation: the shortest walk is not the image's class"
            (logger.info if incomplete else logger.warning)(
                "class %s: the registered image of member %s reads %r but the %s itinerary "
                "reads %r; they part at token %d: %s",
                cd.letter, evidence.bridge_id, evidence.word,
                cd.source, cd.word, diagnosis.position, why,
            )
        verdicts = [e.word_agrees for e in cd.evidence if e.word_agrees is not None]
        cd.verified = all(verdicts) if verdicts else None
        cd.itinerary_verified = all(evidence.itinerary_agrees for evidence in cd.evidence)
        if cd.verified and not cd.itinerary_verified:
            logger.info(
                "class %s: word verified but the registered itinerary differs from the "
                "%s one at element level",
                cd.letter,
                cd.source,
            )


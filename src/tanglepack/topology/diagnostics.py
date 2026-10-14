"""
Geometric self-diagnostics of the hole and partition algorithm: flags, never decisions.

The crossing signs decide every row and every opening, and a hole's side is
decided once at its reference hole and carried along its orbit (D6,
``StablePartition``). This module re-reads the same quantities off the
polylines, where a coarse polyline can disagree with them:

* **Rows.** :func:`bridge_row_diagnosis` reads the side of the stable branch a
  bridge approaches each end from on its polyline. On a 2-4 node bridge the
  "approach" is the vector to the far node (NUM-1), so a disagreement with the
  sign is an under-resolved bridge.
* **Hole sides.** A propagated hole's side MEASURED on the curves: its orbit's
  reference hole carried back by the real inverse map and read on the image
  sub-arc, which escapes the lobe when the reference lobe is under-resolved
  (NUM-2).

:func:`geometry_flags` collects both and :func:`warn_geometry_flags` logs one
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
from typing import Iterable, Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from ..numerics.geometry import oriented_bridge_polyline
from .StablePartition import row_of_end, rows_of_bridge
from .TopologyResults import Hole, Side, SideSource, side_after

if TYPE_CHECKING:
    from .Trellis import Trellis
    from ..numerics.Bridge import Bridge
    from ..numerics.Intersection import Intersection, ManifoldKey

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def _side_of(tangent: NDArray[np.float64], displacement: NDArray[np.float64]) -> Optional[Side]:
    """Side of the displacement relative to the tangent, standard orientation.

    Standing on the manifold facing along ``tangent`` (toward the anchor),
    positive cross = left hand — the author's dynamical-direction convention
    (2026-07-09), which supersedes the earlier flipped calibration.
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
    nodes = manifold.get_point_array(return_nodes=True)
    cdists = [node.cdist for node in nodes]
    here = np.asarray(ix.get_point(), dtype=np.float64)
    below = bisect_left(cdists, ix.stable_cdist) - 1
    above = bisect_right(cdists, ix.stable_cdist)
    if below >= 0:
        look = np.asarray(nodes[below].get_point(), dtype=np.float64) - here
    elif above < len(nodes):
        look = here - np.asarray(nodes[above].get_point(), dtype=np.float64)
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
    flags (NUM-1).
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
        # the point than the sub-arc it belongs to (2026-10-02).
        nodes = bridge.get_point_array(return_nodes=True)
        cdists = sorted(node.cdist for node in nodes)
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
            that decides (D6).
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
    cross-check of :func:`row_of_end`, never the decision (D6).

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
    again. The rows are the crossing signs' (:func:`rows_of_bridge`, D6), so
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
        side_source: Where its ``bridge_side`` came from (``Hole.side_source``).
        bridge_side: The stored side.
        measured_side: The side measured on the curves (``Hole.measured_side``).
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
    holes: Iterable[Hole], *, orientation_preserving: bool = True
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
    their side (author, 2026-10-07) this guards the direct holes of one
    orbit (its references, each read off its own crossing signs) and the
    carrying itself.

    The MEASURED reading is the geometric cross-check: a hole whose
    ``measured_side`` (the reference hole carried back by the real inverse
    map, read on the image sub-arc) contradicts its stored side is
    under-resolved. It is a flag (:func:`geometry_flags`), never an error (D6).

    Holes with no ``origin`` or no comparable side skip the orbit check, and
    holes with no ``measured_side`` skip the measured one; both counts are
    logged at INFO.

    Args:
        holes: The holes to check, typically ``trellis.holes``.
        orientation_preserving: True when the map has det J > 0.

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
    unmeasured = sum(
        1 for hole in holes if hole.side_source == "carried" and hole.measured_side is None
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
            side_source=hole.side_source,
            bridge_side=hole.bridge_side,
            measured_side=hole.measured_side,
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
    The geometric readings that disagree with the sign rule (D6).

    The crossing signs decide rows and openings and the orbit carries hole
    sides; the polylines only cross-check them. A bridge whose polyline reads
    another row than its signs at some end, or a propagated hole whose
    carried point lies on the other side of its image sub-arc, marks an
    under-resolved polyline (or a carried point that left its lobe): a reason
    to refine, never an error.

    Args:
        trellis: The Trellis carrying the bridges and manifolds.
        bridges: The bridges to read, typically the hole bridges.
        holes: The holes whose measured sides to compare.

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
    sides = [
        diagnosis
        for diagnosis in side_diagnoses(
            holes, orientation_preserving=trellis.orientation_preserving
        )
        if diagnosis.measured_disagrees
    ]
    return rows, sides


def _warn_geometry_flags(
    trellis: "Trellis", bridges: Iterable["Bridge"], holes: Iterable[Hole]
) -> None:
    """One WARNING per kind of :func:`geometry_flags` (a resolution flag, D6)."""
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

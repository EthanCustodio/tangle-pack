"""
Walks over the dual graph: where a homotopy element lands, and the shortest
route between two iterated elements.

This module is the walking half of the symbolic dynamics. Given a bridge
class ``{X, Y}`` of HOMOTOPY elements, :func:`land_element` says which
ITERATED element each of ``X`` and ``Y`` maps into under one map step
(:class:`ElementLanding`), and :func:`shortest_walks` walks the
:class:`~tanglepack.topology.DualGraph.DualGraph` from one landing to the
other, crossing the stable manifold only at unified (traversable) nodes and
recording, at every crossing, the element on the side it entered from and
then the element on the side it left to (:class:`Walk`, :class:`WalkSearch`).
The itinerary has even length and splits into disjoint consecutive pairs,
each of which is one bridge class (or its inverse) -- the caller
(``SymbolicDynamics``) turns those pairs into symbols. Between the landing
spans, the walk obeys the non-crossing rule (D13): :func:`lobe_image_markers`
reads which markers a member's lobe holds and :func:`holds_markers` whether a
walk's closed curve holds their images, and the decided walks of classes on
other unstable branches are walls a tie may not cross (:func:`walk_chords`).

A singleton element ``[x, x]`` owns no stable edge and has no dual-graph
node, so a class with a singleton landing cannot be walked. For it
:func:`trellis_itinerary` reads the same even-length itinerary off the
REGULAR trellis: the chain of consecutive registered crossings along the
image of a member bridge, each interior crossing contributing its entry-row
element and then its exit-row element.

Everything here is lookup only: the map is never called. Images come from
:meth:`~tanglepack.topology.Trellis.Trellis.iterate` and
:meth:`~tanglepack.topology.Trellis.Trellis.image_cdist` (iterate table
first, scaled canonical distance as the fallback) and rows from
:func:`~tanglepack.topology.StablePartition.row_of_end` (crossing sign
alone).

Dev Notes:

* ``dual`` in :func:`shortest_walks` is duck-typed: it needs only
  ``stable_nodes_of(ref)``. The walk reads faces off the nodes themselves
  (``StableNode.faces``, ``FaceNode.exits``), so tests build synthetic
  :class:`~tanglepack.topology.DualGraph.StableNode` /
  :class:`~tanglepack.topology.DualGraph.FaceNode` objects with ``arc=None``.
* The walk iterates a node's sides EXPLICITLY (``node.faces[side] is F``)
  rather than calling ``StableNode.side_of`` / ``other_face``: when one
  merged face lies on both sides of a unified edge ``side_of`` answers the
  first side and ``other_face`` returns the same face, which is exactly the
  self-adjacency the walk must skip.
* Landing does NOT cross-check the midpoint owner against the owner of a
  member endpoint's ``+1`` iterate: an image endpoint on an existing hole
  boundary (never re-cut) can be owned by the neighbour element, and the
  owner of a boundary point says nothing about where the span's interior
  lands. The comparison is logged at DEBUG only.
* The reverse-direction trellis itinerary is the element-wise REVERSE of
  the forward one. The plan phrases it as "reverse the chain and swap each
  pair", but ``row_of_end`` reads a crossing's row from its position in a
  :data:`~tanglepack.numerics.Bridge.BridgeId` ORDERED BY UNSTABLE CDIST
  (``row(first) = left iff sign > 0``, ``row(second) = left iff sign < 0``),
  so applying it to a swapped pair flips every row. Reading the rows on the
  forward chain and reversing the tuple gives the walk that traverses the
  same bridges the other way: it enters each crossing on the row the forward
  walk left it on, and leaves on the row it entered.
* The non-crossing rule (D13, author 2026-10-09; WP9b). ``f`` is a
  homeomorphism and an unstable arc never crosses another one, so the image
  ``f(b)`` of a member bridge cannot cross the images of the unstable arcs
  around it. Two consequences make it combinatorial. (1) The walk runs
  between the landing SPANS (C14): the hole bridges bounding a homotopy
  element ``X`` map onto the kept image bridges at the ends of ``f(X)``, so
  ``f(b)`` leaves the stable manifold from an edge ``f(X)`` meets. (2) The
  marker law: ``b`` and the stable segment ``sigma_b`` between its ends bound
  a disk, and ``f(b) + f(sigma_b)`` bounds its image, so it holds the image of
  every marker (periodic point, hole) the lobe holds and of no other. The
  word is the shortest walk whose curve, closed by ``sigma'`` along the
  stable manifold, does so. The rule applies to every class with a member
  lobe, not only to ties (C13): a unique shortest walk can break it. A walk
  is looked for up to :data:`MAX_DETOUR` crossings past the shortest one
  (faces may then repeat); the detours found so far are two crossings long.
  A walk crossing ``sigma'``, or whose chords cross in a face, is rejected
  (``f(b)`` is simple and meets ``f(sigma_b)`` only at its ends; pulling it
  tight keeps it so); one the midpoint model cannot draw (an edge crossed
  twice) is kept but never preferred by order, so it can only make a tie.
* Walls (D13, D14; WP9c). Images of disjoint unstable arcs are disjoint, so
  ``f(b)`` cannot cross the image walk of an arc on another unstable branch:
  given the decided walks of such classes (``walls``, :func:`walk_chords`), a
  tie loses every walk with a chord crossing a wall in its face
  (:meth:`~tanglepack.topology.DualGraph.DualGraph.chords_cross`); it counts
  in :attr:`WalkSearch.rejected`. This is the rule for a class with no member
  lobe (member-less, D11, or with its ends on two branches or two tangles),
  and it applies to member-lobe ties too. A walk standing for several face
  paths (``multiplicity > 1``) is never rejected, and its class makes no
  wall: which face the image runs through is unknown. When every walk crosses
  a wall the tie stands (never observed).
* When the ``max_walks`` cap cuts the enumeration short the status still
  says ``"unique"`` or ``"ambiguous"`` for the itineraries that survived and
  :attr:`WalkSearch.truncated` is set (with a WARNING), because a truncated
  search still has a best candidate to hand on. Only a constrained search
  (D13) that runs out before any judged walk holds the markers reports
  ``"truncated"`` as its status, with no walk: it has no candidate, and
  ``"no_admissible_walk"`` would claim every walk up to the detour cap broke
  the rule. The cap counts judgements there, and the other walks of an
  itinerary already holding the markers are counted, not judged.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Iterable, Iterator, Literal, Optional, Sequence, TYPE_CHECKING

from .BridgeClass import element_sort_key
from .StablePartition import owns_cdist, row_of_end
from .TopologyResults import (
    ElementRef,
    Endpoint,
    OPPOSITE_SIDE,
    PartitionInterval,
    Side,
    endpoint_index,
    side_after,
)
from .Trellis import SCALING_RTOL

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ..numerics.Bridge import BridgeId
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.Intersection import ManifoldKey
    from .BridgeClass import BridgeClassEntry
    from .DualGraph import Chord, FaceNode, Marker, StableNode
    from .PartitionFamily import PartitionFamily
    from .Trellis import Trellis

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


# --------------------------------------------------------------------------- #
# Landing: where a homotopy element goes under one map step
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class ElementLanding:
    """
    Where one HOMOTOPY element lands, under one map step, in the ITERATED partition.

    Attributes:
        source: The homotopy element mapped (a ref into the homotopy family).
        image_key: The stable branch its image lies on
            (``FixedPoint.advance_key(branch_key, 1)``).
        image_side: The side its image lies on: the same side on an
            orientation-preserving map, the opposite one otherwise.
        span: ``(lo, hi)`` stable canonical distances of the image on
            ``image_key``. Both equal for a singleton.
        from_table: Per end of ``span``, whether the distance came from the
            iterate table (``True``) or from the scaling law (``False``). The
            anchor end (cdist 0) counts as a table answer: it is exact.
        target: The iterated element the image lands in, or None when
            unresolved (see ``reason``).
        contained: True when the whole span lies inside ``target`` (to the
            per-end tolerance); False when it straddles a cut, in which case
            ``target`` is the owner of the span's midpoint and ``covering``
            lists every element the span meets.
        singleton: True when ``source`` is a singleton ``[x, x]``: its image
            is the registered iterate of ``x`` and ``target`` is that
            crossing's owner (no dual-graph node exists for it).
        covering: Every iterated element meeting the span, anchor outward.
            ``(target,)`` when contained; empty when unresolved.
        reason: Why the landing is unresolved, or None.
    """

    source: ElementRef
    image_key: "ManifoldKey"
    image_side: Side
    span: tuple[float, float]
    from_table: tuple[bool, bool]
    target: Optional[ElementRef]
    contained: bool
    singleton: bool
    covering: tuple[ElementRef, ...]
    reason: Optional[str]

    @property
    def resolved(self) -> bool:
        """True when the landing names a target element."""
        return self.target is not None

    def __repr__(self) -> str:
        target = "unresolved" if self.target is None else self.target.label
        flags = []
        if self.singleton:
            flags.append("singleton")
        if not self.contained and self.target is not None:
            flags.append("straddles")
        extra = f" [{', '.join(flags)}]" if flags else ""
        return (
            f"ElementLanding({self.source.label} -> {target}, "
            f"span=({self.span[0]:.4g}, {self.span[1]:.4g}){extra})"
        )


def _is_singleton(interval: PartitionInterval) -> bool:
    """A pinched element ``[x, x]``: both ends the same registered crossing."""
    return interval.lo_id is not None and interval.lo_id == interval.hi_id


def _end_tolerance(cdist: float, from_table: bool) -> float:
    """The containment tolerance at one end of an image span: none at a registered
    end, the scaling law's relative error at a scaled one."""
    return 0.0 if from_table else SCALING_RTOL * abs(cdist)


def _meets_span(
    interval: PartitionInterval, lo: float, hi: float, tol_lo: float, tol_hi: float
) -> bool:
    """
    Whether an element overlaps ``[lo, hi]`` with positive length, or owns one
    of its ends (a neighbour merely touching an open end does not count).
    """
    overlaps = interval.lo_cdist < hi - tol_hi and interval.hi_cdist > lo + tol_lo
    return overlaps or owns_cdist(interval, lo, tol_lo) or owns_cdist(interval, hi, tol_hi)


def _unresolved(
    ref: ElementRef,
    image_key: "ManifoldKey",
    image_side: Side,
    span: tuple[float, float],
    from_table: tuple[bool, bool],
    reason: str,
    *,
    singleton: bool = False,
) -> ElementLanding:
    """An unresolved landing with its reason (logged at WARNING by the caller)."""
    return ElementLanding(
        source=ref,
        image_key=image_key,
        image_side=image_side,
        span=span,
        from_table=from_table,
        target=None,
        contained=False,
        singleton=singleton,
        covering=(),
        reason=reason,
    )


def _land_singleton(
    trellis: "Trellis",
    iterated: "PartitionFamily",
    ref: ElementRef,
    interval: PartitionInterval,
    image_key: "ManifoldKey",
    image_side: Side,
) -> ElementLanding:
    """The singleton path: the registered iterate of the pinched crossing."""
    crossing_id = interval.lo_id
    assert crossing_id is not None
    image_id = trellis.iterate(crossing_id, 1)
    key, cdist, from_table = trellis.image_cdist(crossing_id, 1, "stable", branch=ref.branch_key)
    span = (float(cdist), float(cdist))
    flags = (bool(from_table), bool(from_table))
    if image_id not in trellis:
        reason = (
            f"image crossing of singleton element {ref.label} (crossing "
            f"{crossing_id}) is no crossing of the trellis; grow or blast"
        )
        logger.warning("%s", reason)
        return _unresolved(ref, image_key, image_side, span, flags, reason, singleton=True)
    # infer_iterates links a crossing only to one on the advanced branches, and
    # the scaled fallback IS advance_key, so the image is always on image_key.
    assert key == image_key, f"singleton {ref.label} imaged onto {key[1:]}, not {image_key[1:]}"
    try:
        target = iterated.owner_of_intersection(image_id, image_side, image_key)
    except ValueError as error:
        reason = (
            f"image crossing {image_id} of singleton element {ref.label} has no "
            f"owner on the {image_side} side of the iterated partition: {error}"
        )
        logger.warning("%s", reason)
        return _unresolved(ref, image_key, image_side, span, flags, reason, singleton=True)
    target_interval = iterated.element(target)
    if not _is_singleton(target_interval):
        # Expected, not a lookup error: holes propagate backward only, so the
        # forward image of the hole pinching this singleton is never punched
        # and its image crossing is a boundary point of a full-width element
        # (or an interior one); a singleton at the trimmed end of the manifold
        # likewise has an image with nothing trimmed beyond it.
        logger.info(
            "singleton element %s maps to crossing %d, owned on the %s side by "
            "NON-singleton element %s [%.6g, %.6g]: the hole that would pinch the "
            "image is a forward iterate and holes propagate backward only (or the "
            "singleton sits at the trimmed end of the manifold)",
            ref.label,
            image_id,
            image_side,
            target.label,
            target_interval.lo_cdist,
            target_interval.hi_cdist,
        )
    logger.debug("landed singleton %s at crossing %d in %s", ref.label, image_id, target.label)
    return ElementLanding(
        source=ref,
        image_key=image_key,
        image_side=image_side,
        span=span,
        from_table=flags,
        target=target,
        contained=True,
        singleton=True,
        covering=(target,),
        reason=None,
    )


def land_element(
    trellis: "Trellis",
    homotopy: "PartitionFamily",
    iterated: "PartitionFamily",
    ref: ElementRef,
) -> ElementLanding:
    """
    Map one homotopy element forward one step into the iterated partition.

    The image branch is ``FixedPoint.advance_key(branch_key, 1)`` and the
    image side is the element's own side, or the opposite one on an
    orientation-reversing map. A singleton element goes through the iterate
    table (``Trellis.iterate(x, 1)``) and lands in the owner of that crossing.
    Any other element has its two ends mapped -- the anchor end to cdist 0,
    a registered end through :meth:`~.Trellis.Trellis.image_cdist` (table
    first, scaled fallback), an unbounded high end by the scaling law -- and
    lands in the iterated element owning the midpoint of the image span.
    Containment of the whole span is exact at a table end and checked to
    ``SCALING_RTOL * |c|`` at a scaled end.

    Args:
        trellis: The trellis resolving iterates, image distances and the
            orientation of the map.
        homotopy: The family ``ref`` belongs to.
        iterated: The family the image is located in (normally the iterated
            homotopy partition; any family covering the image branch and side
            works).
        ref: The homotopy element to map.

    Returns:
        The :class:`ElementLanding`. It is unresolved (``target is None``,
        ``reason`` set, WARNING logged) when a singleton's image crossing is
        not a crossing of the trellis, or when the iterated partition has no
        single owner for the probe point. A span that straddles a cut is
        resolved to its midpoint owner with ``contained=False`` and
        ``covering`` filled (WARNING).

    Raises:
        ValueError: If ``ref`` names no element of ``homotopy``.
        AssertionError: If an image end lands on a branch other than
            ``advance_key``'s (the iterate table links only there).
    """
    interval = homotopy.element(ref)
    fixed_point: "FixedPoint" = ref.branch_key[0]
    image_key = fixed_point.advance_key(ref.branch_key, 1)
    image_side = side_after(ref.side, 1, trellis.orientation_preserving)

    if _is_singleton(interval):
        return _land_singleton(trellis, iterated, ref, interval, image_key, image_side)

    ends: list[tuple[float, bool]] = []
    for end_id, end_cdist, is_low in (
        (interval.lo_id, interval.lo_cdist, True),
        (interval.hi_id, interval.hi_cdist, False),
    ):
        if end_id is None:
            if is_low:
                ends.append((0.0, True))
            else:
                ends.append((float(end_cdist) * fixed_point.per_step_beta("stable"), False))
            continue
        key, cdist, from_table = trellis.image_cdist(end_id, 1, "stable", branch=ref.branch_key)
        # Always the advanced branch: the table links only there, and the
        # scaled fallback is advance_key itself.
        assert key == image_key, f"end {end_id} of {ref.label} imaged onto {key[1:]}, not {image_key[1:]}"
        ends.append((float(cdist), bool(from_table)))

    (lo, lo_table), (hi, hi_table) = ends
    if lo > hi:
        (lo, lo_table), (hi, hi_table) = (hi, hi_table), (lo, lo_table)
    span = (lo, hi)
    flags = (lo_table, hi_table)
    probe = lo if hi <= lo else 0.5 * (lo + hi)

    try:
        owner = iterated.interval_at(image_key, image_side, probe)
    except ValueError as error:
        reason = (
            f"no single iterated element owns the image of {ref.label} at "
            f"cdist {probe:.6g} on the {image_side} side of branch "
            f"{image_key[1:]}: {error}"
        )
        logger.warning("%s", reason)
        return _unresolved(ref, image_key, image_side, span, flags, reason)
    target = iterated.ref(image_key, image_side, owner.element_id)

    tol_lo = _end_tolerance(lo, lo_table)
    tol_hi = _end_tolerance(hi, hi_table)
    contained = owner.lo_cdist - tol_lo <= lo and hi <= owner.hi_cdist + tol_hi
    if contained:
        covering: tuple[ElementRef, ...] = (target,)
    else:
        result = iterated.result(image_key, image_side)
        covering = tuple(
            result.ref(candidate.element_id)
            for candidate in result.intervals
            if _meets_span(candidate, lo, hi, tol_lo, tol_hi)
        )
        logger.warning(
            "image span [%.6g, %.6g] of element %s straddles %d iterated elements "
            "(%s); keeping the midpoint owner %s",
            lo,
            hi,
            ref.label,
            len(covering),
            ", ".join(candidate.label for candidate in covering),
            target.label,
        )

    if logger.isEnabledFor(logging.DEBUG):
        _debug_endpoint_owners(trellis, iterated, ref, interval, image_side, target)
    return ElementLanding(
        source=ref,
        image_key=image_key,
        image_side=image_side,
        span=span,
        from_table=flags,
        target=target,
        contained=contained,
        singleton=False,
        covering=covering,
        reason=None,
    )


def _debug_endpoint_owners(
    trellis: "Trellis",
    iterated: "PartitionFamily",
    ref: ElementRef,
    interval: PartitionInterval,
    image_side: Side,
    target: ElementRef,
) -> None:
    """DEBUG only: report the owners of the registered end images (never used)."""
    for end_id in (interval.lo_id, interval.hi_id):
        if end_id is None:
            continue
        image_id = trellis.iterate(end_id, 1)
        if image_id not in trellis:
            continue
        try:
            owner = iterated.owner_of_intersection(image_id, image_side, target.branch_key)
        except ValueError:
            continue
        if owner != target:
            logger.debug(
                "element %s lands in %s while the image %d of its end %d is owned "
                "by %s (a boundary handed to the neighbour by a closed cut)",
                ref.label,
                target.label,
                image_id,
                end_id,
                owner.label,
            )


# --------------------------------------------------------------------------- #
# Walks
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class WalkStep:
    """
    One crossing of the stable manifold along a walk.

    Attributes:
        node: The unified stable node crossed.
        entry_side: The side of the edge the walk came from.
        exit_side: The side it left to (the opposite side).
        from_face: The face the step started in.
        to_face: The face it ended in.
    """

    node: "StableNode"
    entry_side: Side
    exit_side: Side
    from_face: "FaceNode"
    to_face: "FaceNode"

    @property
    def entry_element(self) -> ElementRef:
        """The element on the side the walk entered from."""
        return self.node.elements[self.entry_side]

    @property
    def exit_element(self) -> ElementRef:
        """The element on the side the walk left to."""
        return self.node.elements[self.exit_side]


@dataclass
class Walk:
    """
    One shortest walk between two iterated elements.

    Attributes:
        steps: The crossings, in order (empty for a trivial walk).
        itinerary: ``(start, e_1, e_2, ..., goal)``: the start element, then
            per step the entry-side element and the exit-side element, then
            the goal. Even length; the disjoint consecutive pairs
            ``(itinerary[2i], itinerary[2i + 1])`` are each one bridge class
            or its inverse.
        faces: The faces visited, ``len(steps) + 1`` of them.
        multiplicity: How many distinct face paths produced this itinerary
            (parallel unified nodes carrying the same elements collapse into
            one walk); under the non-crossing rule, the walks of the
            itinerary not rejected.
        start_node: The stable node of the start element the walk leaves from
            (on its first face).
        goal_node: The stable node of the goal element it arrives at (on its
            last face).
    """

    steps: list[WalkStep]
    itinerary: tuple[ElementRef, ...]
    faces: list["FaceNode"]
    multiplicity: int = 1
    start_node: Optional["StableNode"] = None
    goal_node: Optional["StableNode"] = None

    @property
    def length(self) -> int:
        """The number of stable-manifold crossings."""
        return len(self.steps)

    @property
    def pairs(self) -> list[tuple[ElementRef, ElementRef]]:
        """The disjoint consecutive pairs of the itinerary."""
        return [
            (self.itinerary[index], self.itinerary[index + 1])
            for index in range(0, len(self.itinerary), 2)
        ]

    def __repr__(self) -> str:
        route = " ".join(ref.label for ref in self.itinerary)
        return f"Walk({self.length} step(s), x{self.multiplicity}: {route})"


#: The outcome of a shortest-walk search. :func:`shortest_walks` reports the
#: cap through :attr:`WalkSearch.truncated`; ``"truncated"`` is the status only
#: when a constrained search ran out with no walk (see the module Dev Notes).
WalkStatus = Literal[
    "unique",
    "trivial",
    "ambiguous",
    "unreachable",
    "no_start_node",
    "no_goal_node",
    "no_admissible_walk",
    "truncated",
    "span_meets_no_node",
]

#: How many crossings longer than the shortest walk the non-crossing rule may
#: go looking for an admissible one (D13).
MAX_DETOUR = 4


@dataclass
class WalkSearch:
    """
    The result of :func:`shortest_walks`.

    Attributes:
        start: The iterated element the walk starts in.
        goal: The iterated element it ends in.
        status: ``"trivial"`` (start and goal share a face, or are equal),
            ``"unique"`` (one distinct shortest itinerary), ``"ambiguous"``
            (several distinct shortest itineraries, all kept),
            ``"unreachable"``, ``"no_start_node"``, ``"no_goal_node"``,
            ``"no_admissible_walk"`` (the non-crossing rule rejected every
            walk up to :data:`MAX_DETOUR` crossings past the shortest),
            ``"truncated"`` (the ``max_walks`` cap ran out before any walk
            judged under that rule held the markers) or
            ``"span_meets_no_node"`` (a landing span overlaps none of its
            element's stable nodes).
        walks: The distinct shortest walks in deterministic order (by the
            :func:`~tanglepack.topology.BridgeClass.element_sort_key` of
            their itineraries); empty unless the status is ``"trivial"``,
            ``"unique"`` or ``"ambiguous"``.
        start_faces: The distinct faces on ``start``'s side of the nodes the
            walk may leave from.
        truncated: True when the ``max_walks`` cap stopped the enumeration of
            shortest face paths (of judged walks, when constrained), so
            ``walks`` (and their multiplicities) may be incomplete.
        constrained: True when the markers of a member lobe were imposed
            (the non-crossing rule, D13).
        rejected: How many itineraries, no longer than the one chosen, broke
            the markers or crossed the image of a disjoint arc (a wall).
    """

    start: ElementRef
    goal: ElementRef
    status: WalkStatus
    walks: list[Walk] = field(default_factory=list)
    start_faces: list["FaceNode"] = field(default_factory=list)
    truncated: bool = False
    constrained: bool = False
    rejected: int = 0

    @property
    def walk(self) -> Optional[Walk]:
        """The first walk in deterministic order, or None when there is none."""
        return self.walks[0] if self.walks else None

    @property
    def itinerary(self) -> Optional[tuple[ElementRef, ...]]:
        """The first walk's itinerary, or None."""
        walk = self.walk
        return None if walk is None else walk.itinerary

    def __repr__(self) -> str:
        return (
            f"WalkSearch({self.start.label} -> {self.goal.label}: {self.status}, "
            f"{len(self.walks)} walk(s){', truncated' if self.truncated else ''})"
        )


def _distinct_faces(faces: Iterable["FaceNode"]) -> list["FaceNode"]:
    """Faces deduplicated by index, in index order."""
    by_index = {face.index: face for face in faces}
    return [by_index[index] for index in sorted(by_index)]


def _fixed_points_of(refs: Iterable[ElementRef]) -> list["FixedPoint"]:
    """The refs' fixed points in first-seen order (identity-matched)."""
    found: list["FixedPoint"] = []
    for ref in refs:
        fixed_point = ref.fixed_point
        if not any(candidate is fixed_point for candidate in found):
            found.append(fixed_point)
    return found


_Predecessor = tuple["FaceNode", "StableNode", Side, Side]


def _bfs_layers(
    start_faces: list["FaceNode"], goal_indices: set[int]
) -> tuple[dict[int, list[_Predecessor]], list["FaceNode"], int]:
    """
    Layered BFS over faces through unified nodes.

    Returns:
        ``(predecessors, goals_reached, depth)``: every predecessor
        ``(from_face, node, entry_side, exit_side)`` of each face at its
        minimal depth, the goal faces found in the first layer holding one
        (empty when the search is exhausted) and that layer's depth.
    """
    depth: dict[int, int] = {face.index: 0 for face in start_faces}
    predecessors: dict[int, list[_Predecessor]] = {face.index: [] for face in start_faces}
    frontier = list(start_faces)
    level = 0
    while frontier:
        layer: dict[int, tuple["FaceNode", list[_Predecessor]]] = {}
        for face in frontier:
            for node in face.exits:
                for entry_side in node.sides:
                    if node.faces.get(entry_side) is not face:
                        continue
                    exit_side = OPPOSITE_SIDE[entry_side]
                    other = node.faces.get(exit_side)
                    if other is None or other is face:
                        continue
                    if other.index in depth:
                        continue
                    layer.setdefault(other.index, (other, []))[1].append(
                        (face, node, entry_side, exit_side)
                    )
        if not layer:
            break
        level += 1
        for index, (face, preds) in layer.items():
            depth[index] = level
            predecessors[index] = preds
        reached = [layer[index][0] for index in sorted(layer) if index in goal_indices]
        if reached:
            return predecessors, reached, level
        frontier = [layer[index][0] for index in sorted(layer)]
    return predecessors, [], level


class _Budget:
    """A shared enumeration budget; ``truncated`` flips when a path is refused."""

    def __init__(self, remaining: int) -> None:
        self.remaining = max(int(remaining), 0)
        self.truncated = False

    def take(self) -> bool:
        """Spend one unit; False (and ``truncated``) once the budget is gone."""
        if self.remaining <= 0:
            self.truncated = True
            return False
        self.remaining -= 1
        return True


def _backtrack(
    goal_face: "FaceNode",
    predecessors: dict[int, list[_Predecessor]],
    budget: _Budget,
) -> Iterable[list[_Predecessor]]:
    """
    Every shortest face path into ``goal_face``, as predecessor lists in walk
    order, while the budget lasts.
    """
    stack: list[tuple["FaceNode", list[_Predecessor]]] = [(goal_face, [])]
    while stack:
        face, suffix = stack.pop()
        preds = predecessors.get(face.index, [])
        if not preds:
            if not budget.take():
                return
            yield suffix
            continue
        # Push in reverse so the first-recorded predecessor is explored first.
        for pred in reversed(preds):
            stack.append((pred[0], [pred] + suffix))


def _nodes_meeting(
    nodes: list["StableNode"], landing: ElementLanding, ref: ElementRef
) -> list["StableNode"]:
    """
    The nodes whose edge the landing span overlaps; empty (WARNING) when none.

    The span should overlap an edge by more than its end tolerance (nothing at
    a table end, ``SCALING_RTOL * |c|`` at a scaled one), so a neighbour
    touched by a scaled end's error is no start. A span narrower than its own
    error overlaps no edge by that much, so then the edges it overlaps at all
    are kept: the tolerance drops neighbours, never the one edge a span lies in.
    """
    lo, hi = landing.span
    tol = max(
        (SCALING_RTOL * abs(c) for c, exact in zip(landing.span, landing.from_table) if not exact),
        default=0.0,
    )
    overlaps = [(node, min(hi, node.hi_cdist) - max(lo, node.lo_cdist)) for node in nodes]
    meeting = [n for n, o in overlaps if o > tol] or [n for n, o in overlaps if o > 0]
    if not meeting:
        logger.warning(
            "the landing span [%.6g, %.6g] of %s meets none of its %d stable node(s) "
            "(resolution flag): the class is unresolved",
            lo,
            hi,
            ref.label,
            len(nodes),
        )
    return meeting


def _stable_edges_between(dual, branch_key: "ManifoldKey", lo: float, hi: float) -> set[tuple]:
    """The edges of one stable branch inside ``[lo, hi]`` (registry cdists, compared exactly)."""
    return {
        node.edge_key
        for node in dual.stable_nodes_on(branch_key)
        if lo <= node.lo_cdist and node.hi_cdist <= hi
    }


def _lobe(dual, bridge_id: "BridgeId") -> Optional[set[tuple]]:
    """A kept bridge and the stable segment between its ends, or None across two branches."""
    trellis = dual.trellis
    unstable_key = trellis.unstable_key_of(bridge_id)
    branches = {trellis.stable_key_at(end, unstable_key) for end in bridge_id}
    if len(branches) != 1:
        return None
    lo, hi = sorted(float(trellis.intersection(end).stable_cdist) for end in bridge_id)
    return _stable_edges_between(dual, branches.pop(), lo, hi) | {dual.unstable_edges[bridge_id]}


class MembersDisagree(ValueError):
    """Two member lobes of one class hold different marker images: the class is
    not one homotopy class (a resolution flag)."""


def lobe_image_markers(dual, entry: "BridgeClassEntry") -> Optional[dict["Marker", bool]]:
    """
    The images of the markers a member's lobe holds and leaves out (D13).

    A member bridge ``b`` the minimal trellis keeps, with both ends on one
    stable branch, closes with the stable segment ``sigma_b`` between them into
    a Jordan curve: ``b`` meets no stable arc inside, and ``sigma_b`` no part
    of ``b``. The map is a homeomorphism, so ``f(b) + f(sigma_b)`` holds
    ``f(m)`` exactly when ``b + sigma_b`` holds ``m``, for every marker ``m``.

    Args:
        dual: The dual graph (its :attr:`~.DualGraph.DualGraph.markers`).
        entry: The class whose members to read.

    Returns:
        ``{image marker: held}`` for every marker the lobe decides and whose
        image is placed, read off the class's first such member (every other
        lobe agrees); None when the class has no member lobe (its walk is held
        only by the walls).

    Raises:
        MembersDisagree: If two member lobes disagree on a marker.
    """
    lobes: list[tuple["BridgeId", dict["Marker", bool]]] = []
    for member in entry.members:
        if member.is_loop or member.bridge_id not in dual.unstable_edges:  # not kept
            continue
        curve = _lobe(dual, member.bridge_id)
        if curve is None:
            continue
        (held,) = dual.enclosed(curve)
        images = {
            dual.image_marker(marker): inside
            for marker, inside in held.items()
            if inside is not None and dual.image_marker(marker) is not None
        }
        lobes.append((member.bridge_id, images))
    if not lobes:
        return None
    (first_id, first), others = lobes[0], lobes[1:]
    for bridge_id, images in others:
        clash = [marker for marker, inside in images.items() if first.get(marker, inside) != inside]
        if clash:
            raise MembersDisagree(
                f"members {first_id} and {bridge_id} of class {entry.bridge_class.label} "
                f"hold different markers ({', '.join(map(str, clash))})"
            )
    return first


def walk_chords(walk: Walk) -> list[tuple["FaceNode", "Chord"]]:
    """
    The walk's chord in each face it visits.

    Each chord runs from the edge the walk entered the face by (its start edge
    first) to the edge it leaves by (its goal edge last), each end named
    ``(edge key, side)`` on the side the face lies on.
    """
    ends_in = [(walk.start_node.edge_key, walk.itinerary[0].side)] + [
        (step.node.edge_key, step.exit_side) for step in walk.steps
    ]
    ends_out = [(step.node.edge_key, step.entry_side) for step in walk.steps] + [
        (walk.goal_node.edge_key, walk.itinerary[-1].side)
    ]
    return [(face, end_in + end_out) for face, end_in, end_out in zip(walk.faces, ends_in, ends_out)]


def _half_toward(node: "StableNode", other: "StableNode") -> tuple:
    """The half of a node's edge facing another node of its branch, named by its end vertex."""
    return (node.edge_key, node.arc.hi_id if other.lo_cdist >= node.hi_cdist else node.arc.lo_id)


def holds_markers(dual, walk: Walk, wanted: dict["Marker", bool]) -> Optional[bool]:
    """
    Whether a walk, closed along the stable manifold, holds the wanted markers (D13).

    The walk runs from the middle of its start edge through its faces to the
    middle of its goal edge; ``sigma'``, the stable edges between the two plus
    the half of each end edge facing the other, closes it. The curve stands
    for ``f(b) + f(sigma_b)``, so each marker in ``wanted`` must be held or
    left out as ``wanted`` says; markers outside ``wanted`` are free.

    Args:
        dual: The dual graph.
        walk: A walk with its ``start_node`` and ``goal_node``.
        wanted: ``{marker: held}`` (:func:`lobe_image_markers`).

    Returns:
        True when SOME routing through the faces keeps the curve simple and
        agrees with ``wanted`` on every marker of ``wanted`` the curve decides
        (held where ``wanted`` says held, left out where it says left out).
        False when none does, which includes a walk crossing ``sigma'``
        (``f(b)`` meets ``f(sigma_b)`` only at its ends) and one whose chords
        cross in a face (``f(b)`` is simple, and pulling it tight keeps it
        so). None when the curve is beyond the model: the ends on two
        branches, an edge crossed twice or a walk leaving and re-entering one
        edge (one crossing per edge, at its middle). A trivial walk inside one
        edge holds nothing.
    """
    start, goal = walk.start_node, walk.goal_node
    if start.branch_key != goal.branch_key:
        return None
    if start is goal:
        return None if walk.steps else not any(wanted.values())
    lo, hi = sorted((start.lo_cdist, start.hi_cdist, goal.lo_cdist, goal.hi_cdist))[1:3]
    between = _stable_edges_between(dual, start.branch_key, lo, hi) - {start.edge_key, goal.edge_key}
    crossed = [step.node.edge_key for step in walk.steps]
    if set(crossed) & between:
        return False
    if len(set(crossed)) < len(crossed) or {start.edge_key, goal.edge_key} & set(crossed):
        return None
    halves = {_half_toward(start, goal), _half_toward(goal, start)}
    routings = dual.enclosed(between, walk_chords(walk), halves)
    if routings is None:
        return None
    return any(
        all(held.get(marker) in (None, inside) for marker, inside in wanted.items())
        for held in routings
    )


def _distances_to(goal_faces: list["FaceNode"]) -> dict[int, int]:
    """The fewest crossings from each face to a goal face (face adjacency is symmetric)."""
    distance = {face.index: 0 for face in goal_faces}
    frontier = list(goal_faces)
    while frontier:
        layer = []
        for face in frontier:
            for node in face.exits:
                for other in node.faces.values():
                    if other.index not in distance:
                        distance[other.index] = distance[face.index] + 1
                        layer.append(other)
        frontier = layer
    return distance


def _walks(
    start: ElementRef,
    goal: ElementRef,
    start_nodes: list["StableNode"],
    goal_nodes: list["StableNode"],
    length: int,
) -> Iterator[Walk]:
    """
    Every walk with exactly ``length`` crossings.

    Faces may repeat; a node is crossed at most once and a step never stays
    in its face. A branch is cut once the goal is out of reach in the steps
    left.
    """
    goal_nodes_on: dict[int, list["StableNode"]] = {}
    for node in goal_nodes:
        goal_nodes_on.setdefault(node.face_on(goal.side).index, []).append(node)
    distance = _distances_to([node.face_on(goal.side) for node in goal_nodes])

    def extend(face: "FaceNode", steps: list[WalkStep]) -> Iterator[tuple[list[WalkStep], "StableNode"]]:
        if distance.get(face.index, length + 1) > length - len(steps):
            return
        if len(steps) == length:
            for goal_node in goal_nodes_on.get(face.index, ()):
                yield steps, goal_node
            return
        for node in face.exits:
            if any(step.node is node for step in steps):
                continue
            for entry_side in node.sides:
                other = node.faces.get(OPPOSITE_SIDE[entry_side])
                if node.faces.get(entry_side) is not face or other is None or other is face:
                    continue
                step = WalkStep(node, entry_side, OPPOSITE_SIDE[entry_side], face, other)
                yield from extend(other, steps + [step])

    for start_node in start_nodes:
        first = start_node.face_on(start.side)
        for steps, goal_node in extend(first, []):
            itinerary = [start]
            for step in steps:
                itinerary += [step.entry_element, step.exit_element]
            itinerary.append(goal)
            faces = [first] + [step.to_face for step in steps]
            yield Walk(steps, tuple(itinerary), faces, 1, start_node, goal_node)


def _admissible_walks(
    dual,
    start: ElementRef,
    goal: ElementRef,
    start_nodes: list["StableNode"],
    goal_nodes: list["StableNode"],
    shortest: int,
    wanted: dict["Marker", bool],
    budget: _Budget,
) -> tuple[list[Walk], int, Optional[int]]:
    """
    The admissible walks of the least length in ``[shortest, shortest + MAX_DETOUR]``.

    An itinerary is admissible when one of its walks holds the wanted markers,
    or cannot be judged (kept, so it is never passed over by order alone).
    Each judgement spends one unit of the budget; once an itinerary has a
    walk holding the markers its other walks are counted, not judged.

    Returns:
        ``(walks, rejected itineraries, their length)``; no walks and length
        None when every itinerary up to the detour cap broke the markers, no
        walks and the length reached when the budget ran out first.
    """
    rejected = 0
    for length in range(shortest, shortest + MAX_DETOUR + 1):
        judged: dict[tuple[ElementRef, ...], list[tuple[Optional[bool], Walk]]] = {}
        for walk in _walks(start, goal, start_nodes, goal_nodes, length):
            verdicts = judged.get(walk.itinerary)
            if verdicts and verdicts[0][0]:  # held: one walk is enough
                verdicts[0][1].multiplicity += 1
                continue
            if not budget.take():
                break
            verdict = holds_markers(dual, walk, wanted)
            verdicts = judged.setdefault(walk.itinerary, [])
            verdicts.insert(0 if verdict else len(verdicts), (verdict, walk))
        kept = []
        for verdicts in judged.values():
            admissible = [walk for verdict, walk in verdicts if verdict is not False]
            if admissible:
                admissible[0].multiplicity += len(admissible) - 1
                kept.append(admissible[0])
        rejected += len(judged) - len(kept)
        if kept or budget.truncated:
            return kept, rejected, length
    return [], rejected, None


def shortest_walks(
    dual,
    start: ElementRef,
    goal: ElementRef,
    *,
    max_walks: int = 64,
    fixed_points: Optional[Sequence["FixedPoint"]] = None,
    landings: Optional[tuple[ElementLanding, ElementLanding]] = None,
    wanted: Optional[dict["Marker", bool]] = None,
    walls: Sequence[tuple["FaceNode", "Chord"]] = (),
) -> WalkSearch:
    """
    Every shortest walk of the dual graph from one iterated element to another.

    A walk lives in the faces of the plane and crosses the stable manifold
    only at unified nodes (the strong pip's fundamental segment); walls are
    never crossed. It starts in a face on ``start``'s side of one of its
    stable nodes and ends in a face on ``goal``'s side of one of its nodes.
    A layered breadth-first search over faces records every predecessor at
    minimal depth and stops after the first layer holding a goal face; all
    shortest face paths are then backtracked, turned into itineraries and
    deduplicated.

    With ``landings``, the walk leaves from and arrives at the nodes whose
    edges the landing SPANS meet, not every node of the iterated element: the
    hole bridges bounding a homotopy element ``X`` map onto the kept image
    bridges at the ends of ``f(X)``, which ``f(b)`` cannot cross (D13). With
    ``wanted``, the non-crossing rule: the word is the shortest walk whose
    closed curve holds exactly the wanted markers (:func:`holds_markers`),
    looked for up to :data:`MAX_DETOUR` crossings past the shortest walk
    (faces may then repeat). With ``walls``, a tie loses the walks whose
    chords cross a wall in its face (:meth:`~.DualGraph.DualGraph.chords_cross`):
    the walls stand for the images of arcs disjoint from this class's, and
    images of disjoint arcs are disjoint (D13, D14). When every walk crosses
    one, the tie stands.

    Args:
        dual: Anything with ``stable_nodes_of(ref) -> list[StableNode]``
            (a :class:`~tanglepack.topology.DualGraph.DualGraph`, or a test
            double over synthetic nodes); with ``wanted``, the dual graph's
            markers and enclosure too.
        start: The iterated element to start from.
        goal: The iterated element to reach.
        max_walks: Cap on the number of shortest face paths enumerated (on
            the walks judged against the markers, with ``wanted``: the other
            walks of an itinerary already holding them are counted, not
            judged). When it is hit ``truncated`` is set and a WARNING
            logged; the status still reports the itineraries that survived,
            or is ``"truncated"`` when none did.
        fixed_points: The ordering context for
            :func:`~tanglepack.topology.BridgeClass.element_sort_key`. By
            default it is derived from the refs themselves (start, goal, then
            every itinerary element, in first-seen order), which is
            deterministic for a given graph.
        landings: The landings ``start`` and ``goal`` are the targets of.
        wanted: ``{marker: held}`` from a member lobe
            (:func:`lobe_image_markers`), or None for no constraint.
        walls: ``(face, chord)`` per face of the decided walks of classes
            whose arcs lie on other unstable branches (:func:`walk_chords`);
            with them, ``dual`` needs ``chords_cross`` too. A walk of
            ``multiplicity > 1`` is never rejected by a wall.

    Returns:
        The :class:`WalkSearch`. ``"no_start_node"`` / ``"no_goal_node"``
        when an element owns no stable node (a singleton, or a ref of the
        wrong family); ``"trivial"`` when ``start == goal`` or a start face
        is a goal face (itinerary ``(start, goal)``); ``"unique"`` for one
        distinct shortest itinerary; ``"ambiguous"`` for several (all kept,
        WARNING); ``"unreachable"`` when the search exhausts the faces;
        ``"no_admissible_walk"`` when the non-crossing rule rejects every
        walk up to the detour cap; ``"truncated"`` when the cap ran out
        before any judged walk held the markers; ``"span_meets_no_node"``
        when a landing span overlaps none of its element's stable nodes.
    """
    start_nodes = list(dual.stable_nodes_of(start))
    if not start_nodes:
        logger.warning(
            "walk %s -> %s: the start element owns no stable node (a singleton "
            "element, or a ref outside the dual graph's family)",
            start.label,
            goal.label,
        )
        return WalkSearch(start, goal, "no_start_node")
    goal_nodes = list(dual.stable_nodes_of(goal))
    if not goal_nodes:
        logger.warning(
            "walk %s -> %s: the goal element owns no stable node (a singleton "
            "element, or a ref outside the dual graph's family)",
            start.label,
            goal.label,
        )
        return WalkSearch(start, goal, "no_goal_node")
    if landings is not None:
        start_nodes = _nodes_meeting(start_nodes, landings[0], start)
        goal_nodes = _nodes_meeting(goal_nodes, landings[1], goal)
        if not start_nodes or not goal_nodes:
            return WalkSearch(start, goal, "span_meets_no_node")

    start_faces = _distinct_faces(node.face_on(start.side) for node in start_nodes)
    goal_faces = _distinct_faces(node.face_on(goal.side) for node in goal_nodes)
    if len(start_faces) > 1:
        logger.info(
            "walk %s -> %s: the start element faces %d faces (%s); all are "
            "searched from at depth 0",
            start.label,
            goal.label,
            len(start_faces),
            ", ".join(str(face.index) for face in start_faces),
        )
    goal_indices = {face.index for face in goal_faces}

    def node_on(nodes: list["StableNode"], side: Side, face: "FaceNode") -> "StableNode":
        return next(node for node in nodes if node.face_on(side) is face)

    # start == goal shares every start face (without landings, which may
    # differ between the two ends of a class landing in one element).
    shared = [face for face in start_faces if face.index in goal_indices]
    if shared and wanted is None:
        walk = Walk(
            [],
            (start, goal),
            [shared[0]],
            multiplicity=len(shared),
            start_node=node_on(start_nodes, start.side, shared[0]),
            goal_node=node_on(goal_nodes, goal.side, shared[0]),
        )
        return WalkSearch(start, goal, "trivial", [walk], start_faces)
    shortest = 0
    if not shared:
        predecessors, reached, shortest = _bfs_layers(start_faces, goal_indices)
        if not reached:
            logger.warning(
                "walk %s -> %s: unreachable; every route crosses a wall (the stable "
                "manifold off the fundamental segment)",
                start.label,
                goal.label,
            )
            return WalkSearch(start, goal, "unreachable", [], start_faces)

    budget = _Budget(max_walks)
    rejected = 0
    if wanted is None:
        by_itinerary: dict[tuple[ElementRef, ...], Walk] = {}
        for goal_face in reached:
            for path in _backtrack(goal_face, predecessors, budget):
                steps = [
                    WalkStep(node, entry_side, exit_side, from_face, to_face)
                    for (from_face, node, entry_side, exit_side), to_face in zip(
                        path, [pred[0] for pred in path[1:]] + [goal_face]
                    )
                ]
                itinerary: list[ElementRef] = [start]
                for step in steps:
                    itinerary.append(step.entry_element)
                    itinerary.append(step.exit_element)
                itinerary.append(goal)
                key = tuple(itinerary)
                found = by_itinerary.get(key)
                if found is None:
                    faces = [path[0][0]] + [step.to_face for step in steps]
                    by_itinerary[key] = Walk(
                        steps,
                        key,
                        faces,
                        1,
                        node_on(start_nodes, start.side, faces[0]),
                        node_on(goal_nodes, goal.side, goal_face),
                    )
                else:
                    found.multiplicity += 1
            if budget.truncated:
                break
        found_walks = list(by_itinerary.values())
        length: Optional[int] = shortest
    else:
        found_walks, rejected, length = _admissible_walks(
            dual, start, goal, start_nodes, goal_nodes, shortest, wanted, budget
        )
        if rejected:
            logger.info(
                "walk %s -> %s: the non-crossing rule rejected %d itinerar%s",
                start.label,
                goal.label,
                rejected,
                "y" if rejected == 1 else "ies",
            )
        if not found_walks and budget.truncated:
            logger.warning(
                "walk %s -> %s: the %d walks judged under the non-crossing rule ran "
                "out at %d crossings with no walk holding the markers of the class's "
                "member lobe (%d itineraries rejected); the search is truncated",
                start.label,
                goal.label,
                max_walks,
                length,
                rejected,
            )
            return WalkSearch(start, goal, "truncated", [], start_faces, True, True, rejected)
        if not found_walks:
            logger.warning(
                "walk %s -> %s: no walk up to %d crossings holds the markers of the "
                "class's member lobe (%d itineraries rejected)",
                start.label,
                goal.label,
                shortest + MAX_DETOUR,
                rejected,
            )
            return WalkSearch(start, goal, "no_admissible_walk", [], start_faces, False, True, rejected)
    truncated = budget.truncated

    context = (
        list(fixed_points)
        if fixed_points is not None
        else _fixed_points_of([start, goal] + [ref for walk in found_walks for ref in walk.itinerary])
    )
    walks = sorted(
        found_walks,
        key=lambda walk: tuple(element_sort_key(ref, context) for ref in walk.itinerary),
    )
    if len(walks) > 1 and walls:
        # A walk standing for several face paths is never judged on one of them.
        clear = [
            walk
            for walk in walks
            if walk.multiplicity > 1
            or not any(
                wall_face is face and dual.chords_cross(face, chord, wall)
                for face, chord in walk_chords(walk)
                for wall_face, wall in walls
            )
        ]
        if clear and len(clear) < len(walks):
            crossing = len(walks) - len(clear)
            logger.info(
                "walk %s -> %s: %d %s the image of a disjoint arc",
                start.label,
                goal.label,
                crossing,
                "itinerary crosses" if crossing == 1 else "itineraries cross",
            )
            rejected += crossing
            walks = clear
    if truncated:
        logger.warning(
            "walk %s -> %s: more than %d %s; enumeration truncated after %d "
            "distinct itinerar%s",
            start.label,
            goal.label,
            max_walks,
            "shortest face paths" if wanted is None else "walks judged",
            len(walks),
            "y" if len(walks) == 1 else "ies",
        )
    status: WalkStatus = (
        "trivial" if length == 0 else "unique" if len(walks) == 1 else "ambiguous"
    )
    if status == "ambiguous":
        logger.warning(
            "walk %s -> %s: %d distinct shortest itineraries of %d crossing(s); "
            "the class is ambiguous:\n  %s",
            start.label,
            goal.label,
            len(walks),
            walks[0].length,
            "\n  ".join(" ".join(ref.label for ref in walk.itinerary) for walk in walks),
        )
    return WalkSearch(
        start, goal, status, walks, start_faces, truncated, wanted is not None, rejected
    )


# --------------------------------------------------------------------------- #
# The singleton path: an itinerary read off the regular trellis
# --------------------------------------------------------------------------- #
def trellis_itinerary(
    trellis: "Trellis",
    iterated: "PartitionFamily",
    chain: Sequence["BridgeId"],
    direction: int,
) -> tuple[ElementRef, ...]:
    """
    The itinerary of a member bridge's image, read off the regular trellis.

    ``chain`` is the image chain of a member bridge
    (``BridgeClass._image_chain``): the consecutive registered crossing pairs
    ``[(c0, c1), (c1, c2), ...]`` tiling the image arc in increasing unstable
    canonical distance. Each pair is one bridge; the row of each of its ends
    comes from :func:`~tanglepack.topology.StablePartition.row_of_end` and
    the element on that row from the iterated family. The itinerary is the
    element at ``c0`` on the row the first bridge leaves along, then for each
    interior crossing ``c_i`` the element on the row the bridge ARRIVES on
    followed by the element on the row the next bridge LEAVES on, then the
    element at ``c_n`` on the row the last bridge arrives on -- the same
    even-length, disjoint-pair format as a dual-graph :class:`Walk`.

    Args:
        trellis: The trellis resolving the crossings' signs.
        iterated: The family owning the crossings (the iterated homotopy
            partition).
        chain: The image chain, in increasing unstable canonical distance.
        direction: The member's direction: ``+1`` reads the chain as given
            (source end first); ``-1`` traverses the image the other way, so
            the itinerary is the element-wise REVERSE of the forward one (see
            the module Dev Notes for why the rows are still read on the
            forward chain).

    Returns:
        The itinerary, ``2 * len(chain)`` elements long.

    Raises:
        ValueError: If ``chain`` is empty or not consecutive, if a crossing
            has no crossing sign, or if a crossing has no owner on the needed
            side of the iterated partition (the image is not partitioned
            there; grow or blast).
    """
    pairs = [tuple(pair) for pair in chain]
    if not pairs:
        raise ValueError("an image chain must hold at least one crossing pair")
    for previous, following in zip(pairs, pairs[1:]):
        if previous[1] != following[0]:
            raise ValueError(
                f"image chain is not consecutive: {previous} is followed by {following}"
            )

    def owner(pair: tuple[int, int], endpoint: Endpoint, role: str) -> ElementRef:
        try:
            return iterated.owner_at_end(trellis, pair, endpoint)
        except ValueError as error:
            raise ValueError(
                f"crossing {pair[endpoint_index(endpoint)]} of the image chain has no "
                f"owner on the {row_of_end(trellis, pair, endpoint)} side of the "
                f"iterated partition ({role}); the image is not partitioned there "
                f"-- grow or blast. {error}"
            ) from error

    itinerary: list[ElementRef] = [owner(pairs[0], "first", "start of the chain")]
    for previous, following in zip(pairs, pairs[1:]):
        itinerary.append(owner(previous, "second", "entry row"))
        itinerary.append(owner(following, "first", "exit row"))
    itinerary.append(owner(pairs[-1], "second", "end of the chain"))
    if direction < 0:
        itinerary.reverse()
    assert len(itinerary) == 2 * len(pairs)
    return tuple(itinerary)

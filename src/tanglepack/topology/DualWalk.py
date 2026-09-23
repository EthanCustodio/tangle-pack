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
(``SymbolicDynamics``) turns those pairs into symbols.

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
* ``StableNode.images`` are images of ITERATED elements; the landing needs
  images of HOMOTOPY elements, so it does not read them.
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
* ``"truncated"`` is a member of :data:`WalkStatus` for callers, but
  :func:`shortest_walks` never reports it as the status: when the
  ``max_walks`` cap cuts the enumeration short the status still says
  ``"unique"`` or ``"ambiguous"`` for the itineraries that survived and
  :attr:`WalkSearch.truncated` is set (with a WARNING), because a truncated
  search still has a best candidate to hand on.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Iterable, Literal, Optional, Sequence, TYPE_CHECKING

from .BridgeClass import element_sort_key
from .StablePartition import owns_cdist, row_of_end
from .TopologyResults import ElementRef, OPPOSITE_SIDE, PartitionInterval, Side
from .Trellis import SCALING_RTOL

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ..numerics.Bridge import BridgeId
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.Intersection import ManifoldKey
    from .DualGraph import FaceNode, StableNode
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


def _end_tolerance(cdist: float, from_table: bool, cdist_tol: float) -> float:
    """The containment tolerance at one end of an image span."""
    if from_table:
        return cdist_tol
    return max(cdist_tol, SCALING_RTOL * abs(cdist))


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
    key, cdist, from_table = trellis.image_cdist(crossing_id, 1, "stable")
    span = (float(cdist), float(cdist))
    flags = (bool(from_table), bool(from_table))
    if image_id is None:
        reason = (
            f"image crossing of singleton element {ref.label} (crossing "
            f"{crossing_id}) is not registered; grow or blast"
        )
        logger.warning("%s", reason)
        return _unresolved(ref, image_key, image_side, span, flags, reason, singleton=True)
    if key != image_key:
        reason = (
            f"image of singleton element {ref.label} landed on branch {key[1:]} "
            f"instead of {image_key[1:]}"
        )
        logger.warning("%s", reason)
        return _unresolved(ref, image_key, image_side, span, flags, reason, singleton=True)
    try:
        target = iterated.owner_of_intersection(image_id, image_side)
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
    Containment of the whole span is checked with the registry's
    ``cdist_tol`` at a table end and ``SCALING_RTOL * |c|`` at a scaled end.

    Args:
        trellis: The trellis resolving iterates, image distances, the
            orientation of the map and the canonical-distance tolerance.
        homotopy: The family ``ref`` belongs to.
        iterated: The family the image is located in (normally the iterated
            homotopy partition; any family covering the image branch and side
            works).
        ref: The homotopy element to map.

    Returns:
        The :class:`ElementLanding`. It is unresolved (``target is None``,
        ``reason`` set, WARNING logged) when a singleton's image crossing is
        not registered, when an image end lands on a branch other than the
        expected one, or when the iterated partition has no single owner for
        the probe point. A span that straddles a cut is resolved to its
        midpoint owner with ``contained=False`` and ``covering`` filled
        (WARNING).

    Raises:
        ValueError: If ``ref`` names no element of ``homotopy``.
    """
    interval = homotopy.element(ref)
    fixed_point: "FixedPoint" = ref.branch_key[0]
    image_key = fixed_point.advance_key(ref.branch_key, 1)
    image_side: Side = ref.side if trellis.orientation_preserving else OPPOSITE_SIDE[ref.side]
    cdist_tol = float(trellis.registry.cdist_tol)

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
        key, cdist, from_table = trellis.image_cdist(end_id, 1, "stable")
        if key != image_key:
            partial = (float(cdist), float(cdist))
            reason = (
                f"image of crossing {end_id} bounding element {ref.label} landed "
                f"on branch {key[1:]} instead of {image_key[1:]}"
            )
            logger.warning("%s", reason)
            return _unresolved(ref, image_key, image_side, partial, (from_table, from_table), reason)
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

    tol_lo = _end_tolerance(lo, lo_table, cdist_tol)
    tol_hi = _end_tolerance(hi, hi_table, cdist_tol)
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
        if image_id is None:
            continue
        try:
            owner = iterated.owner_of_intersection(image_id, image_side)
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
            one walk).
    """

    steps: list[WalkStep]
    itinerary: tuple[ElementRef, ...]
    faces: list["FaceNode"]
    multiplicity: int = 1

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


#: The outcome of a shortest-walk search. ``"truncated"`` is reserved for
#: callers; :func:`shortest_walks` reports the cap through
#: :attr:`WalkSearch.truncated` instead (see the module Dev Notes).
WalkStatus = Literal[
    "unique",
    "trivial",
    "ambiguous",
    "unreachable",
    "no_start_node",
    "no_goal_node",
    "truncated",
]


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
            ``"unreachable"``, ``"no_start_node"`` or ``"no_goal_node"``.
        walks: The distinct shortest walks in deterministic order (by the
            :func:`~tanglepack.topology.BridgeClass.element_sort_key` of
            their itineraries); empty unless the status is ``"trivial"``,
            ``"unique"`` or ``"ambiguous"``.
        start_faces: The distinct faces on ``start``'s side of its nodes.
        truncated: True when the ``max_walks`` cap stopped the enumeration of
            shortest face paths, so ``walks`` (and their multiplicities) may
            be incomplete.
    """

    start: ElementRef
    goal: ElementRef
    status: WalkStatus
    walks: list[Walk] = field(default_factory=list)
    start_faces: list["FaceNode"] = field(default_factory=list)
    truncated: bool = False

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
) -> tuple[dict[int, list[_Predecessor]], list["FaceNode"]]:
    """
    Layered BFS over faces through unified nodes.

    Returns:
        ``(predecessors, goals_reached)``: every predecessor
        ``(from_face, node, entry_side, exit_side)`` of each face at its
        minimal depth, and the goal faces found in the first layer holding
        one (empty when the search is exhausted).
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
            return predecessors, reached
        frontier = [layer[index][0] for index in sorted(layer)]
    return predecessors, []


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


def shortest_walks(
    dual,
    start: ElementRef,
    goal: ElementRef,
    *,
    max_walks: int = 64,
    fixed_points: Optional[Sequence["FixedPoint"]] = None,
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

    Args:
        dual: Anything with ``stable_nodes_of(ref) -> list[StableNode]``
            (a :class:`~tanglepack.topology.DualGraph.DualGraph`, or a test
            double over synthetic nodes).
        start: The iterated element to start from.
        goal: The iterated element to reach.
        max_walks: Cap on the number of shortest face paths enumerated. When
            it is hit ``truncated`` is set and a WARNING logged; the status
            still reports the itineraries that survived.
        fixed_points: The ordering context for
            :func:`~tanglepack.topology.BridgeClass.element_sort_key`. By
            default it is derived from the refs themselves (start, goal, then
            every itinerary element, in first-seen order), which is
            deterministic for a given graph.

    Returns:
        The :class:`WalkSearch`. ``"no_start_node"`` / ``"no_goal_node"``
        when an element owns no stable node (a singleton, or a ref of the
        wrong family); ``"trivial"`` when ``start == goal`` or a start face
        is a goal face (itinerary ``(start, goal)``); ``"unique"`` for one
        distinct shortest itinerary; ``"ambiguous"`` for several (all kept,
        WARNING); ``"unreachable"`` when the search exhausts the faces.
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

    if start == goal:
        walk = Walk([], (start, goal), [start_faces[0]], multiplicity=len(start_faces))
        return WalkSearch(start, goal, "trivial", [walk], start_faces)
    shared = [face for face in start_faces if face.index in goal_indices]
    if shared:
        walk = Walk([], (start, goal), [shared[0]], multiplicity=len(shared))
        return WalkSearch(start, goal, "trivial", [walk], start_faces)

    predecessors, reached = _bfs_layers(start_faces, goal_indices)
    if not reached:
        logger.warning(
            "walk %s -> %s: unreachable; every route crosses a wall (the stable "
            "manifold off the fundamental segment)",
            start.label,
            goal.label,
        )
        return WalkSearch(start, goal, "unreachable", [], start_faces)

    budget = _Budget(max_walks)
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
                by_itinerary[key] = Walk(steps, key, faces, multiplicity=1)
            else:
                found.multiplicity += 1
        if budget.truncated:
            break
    truncated = budget.truncated

    context = (
        list(fixed_points)
        if fixed_points is not None
        else _fixed_points_of([start, goal] + [ref for key in by_itinerary for ref in key])
    )
    walks = sorted(
        by_itinerary.values(),
        key=lambda walk: tuple(element_sort_key(ref, context) for ref in walk.itinerary),
    )
    if truncated:
        logger.warning(
            "walk %s -> %s: more than %d shortest face paths; enumeration "
            "truncated after %d distinct itinerar%s",
            start.label,
            goal.label,
            max_walks,
            len(walks),
            "y" if len(walks) == 1 else "ies",
        )
    status: WalkStatus = "unique" if len(walks) == 1 else "ambiguous"
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
    return WalkSearch(start, goal, status, walks, start_faces, truncated)


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

    def owner(crossing_id: int, row: Side, role: str) -> ElementRef:
        try:
            return iterated.owner_of_intersection(crossing_id, row)
        except ValueError as error:
            raise ValueError(
                f"crossing {crossing_id} of the image chain has no owner on the "
                f"{row} side of the iterated partition ({role}); the image is "
                f"not partitioned there -- grow or blast. {error}"
            ) from error

    itinerary: list[ElementRef] = [
        owner(pairs[0][0], row_of_end(trellis, pairs[0], "first"), "start of the chain")
    ]
    for previous, following in zip(pairs, pairs[1:]):
        crossing_id = previous[1]
        itinerary.append(
            owner(crossing_id, row_of_end(trellis, previous, "second"), "entry row")
        )
        itinerary.append(
            owner(crossing_id, row_of_end(trellis, following, "first"), "exit row")
        )
    itinerary.append(
        owner(pairs[-1][1], row_of_end(trellis, pairs[-1], "second"), "end of the chain")
    )
    if direction < 0:
        itinerary.reverse()
    assert len(itinerary) == 2 * len(pairs)
    return tuple(itinerary)

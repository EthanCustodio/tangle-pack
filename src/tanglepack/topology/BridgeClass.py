"""
Bridge classes: the homotopy class of arcs a bridge belongs to.

A bridge runs from one crossing to the next along an unstable branch, and each
of its two ends sits against a stable branch on a definite side (its *row*, see
:func:`~.StablePartition.row_of_end`). The stable partition cuts that branch
into elements, so each end names one element. The UNORDERED pair of those two
:class:`~.TopologyResults.ElementRef` s is the bridge's CLASS: every bridge of a
class is an arc between the same two elements, and up to homotopy relative to
the stable manifold they are all the same arc. That is the coarsest identity
the dual graph can tell apart, and the symbol itineraries are written in.

A class has a canonical DIRECTION, ``source -> target``, oriented anchor outward:
``source`` is the element nearer its branch's anchor. A bridge traversed in the
unstable dynamical direction (its :data:`~tanglepack.numerics.Bridge.BridgeId`
order) either runs ``source -> target`` (symbol ``a``, direction ``+1``) or
``target -> source`` (symbol ``a^-1``, direction ``-1``). Both directions belong
to the one class.

Dev Notes:

* The class is read off the row, not off both sides: the bridge lies on ONE
  side of the stable branch at each end, and that is the side whose partition
  the region it bounds belongs to. Taking the other side would name the element
  across the stable manifold, which no bridge of that class touches.
* The row here is the combinatorial one (``crossing_sign`` alone), so an anchor
  bridge — periodic point out to the first crossing, whose first endpoint is a
  synthetic crossing with no manifold nodes around it — is classed like any
  other. Partial bridges have no :data:`~tanglepack.numerics.Bridge.BridgeId`
  and are skipped.
* LOOPS FOLD (user rule, 2026-09-15). A bridge whose two ends land in the SAME
  element is not a class of its own: it is the forward image of a bridge that
  did connect two elements and has been pushed into one of them, and it will
  keep mapping into that element. Such a loop joins the class of the bridge it
  iterated from — its preimage chain, walked with registry lookups only
  (:meth:`~.Trellis.Trellis.iterate` and
  :meth:`~.Trellis.Trellis.bridge_between`; the map is never called) until a
  non-loop ancestor is found. A class holding a loop is INERT: its arcs escape
  into one element and never take part in an itinerary. A loop whose chain is
  broken (no registered preimage, or a partial one) stands alone as a
  ``BridgeClass(x, x)`` entry, flagged :attr:`BridgeClass.is_loop` and logged.
* Inertness is evidence-based. A class is reported active until the trellis has
  been grown far enough for one of its images to be registered as a loop; on
  k=10 the exterior class is active at 9 unstable steps and inert at 10. A loop
  produced by a too-coarse partition (an interior bridge whose ends the
  partition failed to separate) folds by the same rule — the folding is the
  diagnostic, not an error.
* Partitions are passed in rather than read off the trellis because they do not
  live where the bridges do: a trellis holds a bridge under its UNSTABLE fixed
  point, while the partition owning that bridge's endpoints lives on the trellis
  of the STABLE branch's fixed point. Running this on the all-fixed-points
  trellis with every per-fixed-point trellis's partitions is what covers the
  nested and heteroclinic cases in one pass.
* Ordering is deliberate and total (:func:`class_sort_key`): the table is built
  in sorted order and each class's members are sorted by
  :data:`~tanglepack.numerics.Bridge.BridgeId`, so two runs of the same trellis
  enumerate the classes identically.
* Letters (``a``, ``b``, ...) and resonance-zone membership are session-level
  annotations (the loom owns zones and the persistent alphabet); the slots for
  them live on :class:`BridgeClassEntry` so one object carries the whole record.
* Refined symbols — subdividing one class into several symbols that all still
  belong to it — are future work. They will hang off :class:`BridgeClassEntry`
  (a finer grouping of its members); :class:`BridgeClass` itself stays the
  coarse key.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Iterable, Iterator, Optional, Sequence, TYPE_CHECKING

from .StablePartition import row_of_end
from .TopologyResults import (
    ElementRef,
    Endpoint,
    endpoint_index,
    Side,
    StablePartitionResult,
)

if TYPE_CHECKING:
    from ..numerics.Bridge import BridgeId
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.Intersection import ManifoldKey
    from .Trellis import Trellis

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


@dataclass(frozen=True)
class BridgeClass:
    """
    The unordered pair of partition elements a bridge connects, canonically oriented.

    ``source`` is the element nearer its branch anchor (smaller element id on one
    branch and side; across branches the :func:`element_sort_key` order), so the
    pair ``{X, Y}`` has exactly one representation and the class hashes as one
    symbol whichever way its bridges run. Frozen, so it is usable as a dict key
    and as the key of a session alphabet.

    Attributes:
        source: The anchor-nearer element; the start of the class's direction.
        target: The other element. Equal to ``source`` only for an unresolved loop
            (see the module Dev Notes).
    """

    source: ElementRef
    target: ElementRef

    @property
    def is_loop(self) -> bool:
        """True when both elements coincide — an unresolved loop standing alone."""
        return self.source == self.target

    @property
    def label(self) -> str:
        """A short, deterministic name: the two element labels, ``<->``-joined."""
        if self.is_loop:
            return f"{self.source.label} loop"
        return f"{self.source.label} <-> {self.target.label}"

    def __str__(self) -> str:
        return self.label


@dataclass(frozen=True)
class BridgeMember:
    """
    One bridge's membership in a class: which way it runs, and how it got there.

    Attributes:
        bridge_id: The bridge.
        direction: ``+1`` when the bridge runs ``source -> target`` in the
            unstable dynamical direction, ``-1`` for ``target -> source``, ``0``
            for a loop (both ends in one element; ``a`` and ``a^-1`` coincide).
        folded_from: For a folded loop, the non-loop ancestor bridge whose class
            it joined. ``None`` for a direct member and for an unresolved loop.
        loop_element: For a loop, the element both of its ends sit in.
    """

    bridge_id: "BridgeId"
    direction: int
    folded_from: Optional["BridgeId"] = None
    loop_element: Optional[ElementRef] = None

    @property
    def is_loop(self) -> bool:
        """True when the bridge's two ends sit in the same element."""
        return self.direction == 0


@dataclass
class BridgeClassEntry:
    """
    One bridge class with its members and the session-level annotations.

    Attributes:
        bridge_class: The class.
        members: The member bridges, sorted by
            :data:`~tanglepack.numerics.Bridge.BridgeId`.
        letter: The class's symbol letter (``"a"``, ``"b"``, ...), set by the
            session alphabet for active classes; ``None`` when unlettered.
        zone_key: The ``(fixed_point, branch_index)`` key of the resonance zone
            the class lies in, set by the session; ``None`` when it lies outside
            every zone or no zones are defined.
    """

    bridge_class: BridgeClass
    members: list[BridgeMember] = field(default_factory=list)
    letter: Optional[str] = None
    zone_key: Optional[tuple] = None

    @property
    def inert(self) -> bool:
        """True when the class holds a loop: its arcs map into one element."""
        return any(member.is_loop for member in self.members)

    @property
    def active(self) -> bool:
        """The opposite of :attr:`inert`."""
        return not self.inert

    @property
    def bridge_ids(self) -> list["BridgeId"]:
        """The member bridge ids, in member order."""
        return [member.bridge_id for member in self.members]

    @property
    def loops(self) -> list[BridgeMember]:
        """The loop members (folded or unresolved)."""
        return [member for member in self.members if member.is_loop]

    @property
    def name(self) -> str:
        """The letter when there is one, else the class label."""
        return self.letter if self.letter is not None else self.bridge_class.label

    def member(self, bridge_id: "BridgeId") -> BridgeMember:
        """
        The membership record of one bridge.

        Args:
            bridge_id: The bridge.

        Returns:
            Its :class:`BridgeMember`.

        Raises:
            KeyError: If the bridge is not a member of this class.
        """
        for member in self.members:
            if member.bridge_id == bridge_id:
                return member
        raise KeyError(f"bridge {bridge_id} is not a member of class {self.name}")

    def symbol(self, bridge_id: "BridgeId") -> str:
        """
        The symbol a bridge contributes to an itinerary.

        Args:
            bridge_id: A member bridge.

        Returns:
            ``letter`` for a ``source -> target`` traversal, ``letter^-1`` for
            ``target -> source``, and the bare ``letter`` for a loop.

        Raises:
            ValueError: If the class has no letter (inert classes are unlettered).
            KeyError: If the bridge is not a member.
        """
        member = self.member(bridge_id)
        if self.letter is None:
            raise ValueError(
                f"class {self.bridge_class.label} has no letter "
                f"({'inert' if self.inert else 'unlettered'}), so bridge "
                f"{bridge_id} has no symbol"
            )
        return f"{self.letter}^-1" if member.direction < 0 else self.letter

    def describe(self) -> str:
        """One line: name, elements, activity, zone and members with symbols."""
        zone = "no zone" if self.zone_key is None else zone_label(self.zone_key)
        parts = []
        for member in self.members:
            if member.is_loop:
                origin = (
                    f"folded from {member.folded_from}"
                    if member.folded_from is not None
                    else "unresolved"
                )
                parts.append(f"{member.bridge_id} loop in {member.loop_element} ({origin})")
            elif self.letter is not None:
                parts.append(f"{member.bridge_id} {self.symbol(member.bridge_id)}")
            else:
                parts.append(f"{member.bridge_id} {'+' if member.direction > 0 else '-'}")
        state = "inert" if self.inert else "active"
        head = (
            f"{self.letter}: {self.bridge_class.label}"
            if self.letter is not None
            else self.bridge_class.label
        )
        return f"{head} [{state}, {zone}] {len(self.members)} bridge(s): " + ", ".join(parts)


@dataclass
class BridgeClassTable:
    """
    Every bridge class of a trellis, in :func:`class_sort_key` order.

    Attributes:
        entries: The classes, one :class:`BridgeClassEntry` each.
    """

    entries: list[BridgeClassEntry] = field(default_factory=list)

    def __iter__(self) -> Iterator[BridgeClassEntry]:
        return iter(self.entries)

    def __len__(self) -> int:
        return len(self.entries)

    def __getitem__(self, bridge_class: BridgeClass) -> BridgeClassEntry:
        for entry in self.entries:
            if entry.bridge_class == bridge_class:
                return entry
        raise KeyError(f"no entry for class {bridge_class.label}")

    def __contains__(self, bridge_class: object) -> bool:
        return any(entry.bridge_class == bridge_class for entry in self.entries)

    @property
    def classes(self) -> list[BridgeClass]:
        """The classes in table order."""
        return [entry.bridge_class for entry in self.entries]

    @property
    def active(self) -> list[BridgeClassEntry]:
        """The entries that hold no loop."""
        return [entry for entry in self.entries if entry.active]

    @property
    def inert(self) -> list[BridgeClassEntry]:
        """The entries that hold a loop."""
        return [entry for entry in self.entries if entry.inert]

    @property
    def bridge_ids(self) -> list["BridgeId"]:
        """Every classed bridge id, sorted."""
        return sorted(bid for entry in self.entries for bid in entry.bridge_ids)

    def entry_of(self, bridge_id: "BridgeId") -> BridgeClassEntry:
        """
        The entry a bridge belongs to.

        Args:
            bridge_id: The bridge.

        Returns:
            Its :class:`BridgeClassEntry`.

        Raises:
            KeyError: If no class holds the bridge.
        """
        for entry in self.entries:
            for member in entry.members:
                if member.bridge_id == bridge_id:
                    return entry
        raise KeyError(f"bridge {bridge_id} is in no class")

    def member_of(self, bridge_id: "BridgeId") -> BridgeMember:
        """
        The membership record of a bridge.

        Args:
            bridge_id: The bridge.

        Returns:
            Its :class:`BridgeMember`.

        Raises:
            KeyError: If no class holds the bridge.
        """
        return self.entry_of(bridge_id).member(bridge_id)

    def symbol(self, bridge_id: "BridgeId") -> str:
        """The itinerary symbol of a bridge; see :meth:`BridgeClassEntry.symbol`."""
        return self.entry_of(bridge_id).symbol(bridge_id)

    def as_dict(self) -> dict[BridgeClass, list["BridgeId"]]:
        """The table as a plain ``{class: sorted member ids}`` mapping."""
        return {entry.bridge_class: list(entry.bridge_ids) for entry in self.entries}

    def describe(self) -> str:
        """A multi-line report, one :meth:`BridgeClassEntry.describe` line per class."""
        head = (
            f"{len(self.entries)} bridge class(es): {len(self.active)} active, "
            f"{len(self.inert)} inert"
        )
        return "\n".join([head, *(f"  {entry.describe()}" for entry in self.entries)])

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, BridgeClassTable):
            return NotImplemented
        return self.entries == other.entries


def zone_label(zone_key: tuple) -> str:
    """
    A readable name for a ``(fixed_point, branch_index)`` resonance-zone key.

    Args:
        zone_key: The key, as stored on :attr:`BridgeClassEntry.zone_key`.

    Returns:
        ``"zone p{period} branch {index}"``.
    """
    fixed_point, branch_index = zone_key
    period = getattr(fixed_point, "period", "?")
    return f"zone p{period} branch {branch_index}"


def element_sort_key(
    ref: ElementRef, fixed_points: Sequence["FixedPoint"]
) -> tuple[int, int, int, Side, int]:
    """
    A total, run-stable ordering key for one element reference.

    Args:
        ref: The element reference to order.
        fixed_points: The fixed points in the order that defines the outermost
            component of the key (normally ``trellis.fixed_points``). A ref
            whose fixed point is not among them sorts after all of them.

    Returns:
        ``(fixed-point position, orbit index, branch index, side, element id)``.

    Note:
        Fixed points are matched by IDENTITY, not equality — a FixedPoint is a
        live orbit object, and two runs of the same trellis reuse the same one.
        On one branch and side the key reduces to the element id, which counts
        from the anchor outward: that is what makes ``source`` the anchor-nearer
        element of a class.
    """
    position = len(fixed_points)
    for index, fixed_point in enumerate(fixed_points):
        if fixed_point is ref.fixed_point:
            position = index
            break
    return (position, ref.orbit_index, ref.branch_index, ref.side, ref.element_id)


def class_sort_key(
    bridge_class: BridgeClass, fixed_points: Sequence["FixedPoint"]
) -> tuple:
    """
    A total, run-stable ordering key for a bridge class.

    Args:
        bridge_class: The class to order.
        fixed_points: The fixed points defining the outermost ordering (see
            :func:`element_sort_key`).

    Returns:
        The concatenation of :func:`element_sort_key` for ``source`` and then for
        ``target`` — a 10-tuple.
    """
    return element_sort_key(bridge_class.source, fixed_points) + element_sort_key(
        bridge_class.target, fixed_points
    )


def oriented_class(
    x: ElementRef, y: ElementRef, fixed_points: Sequence["FixedPoint"]
) -> tuple[BridgeClass, int]:
    """
    The canonical class of a bridge running ``x -> y`` and the bridge's direction.

    Args:
        x: The element at the bridge's ``"first"`` end (smaller unstable cdist).
        y: The element at its ``"second"`` end.
        fixed_points: The ordering context (see :func:`element_sort_key`).

    Returns:
        ``(BridgeClass(source, target), direction)`` with ``direction = +1`` when
        ``x`` is the source, ``-1`` when ``y`` is, and ``0`` when ``x == y``.
    """
    if x == y:
        return BridgeClass(x, x), 0
    if element_sort_key(x, fixed_points) <= element_sort_key(y, fixed_points):
        return BridgeClass(x, y), +1
    return BridgeClass(y, x), -1


def bridge_classes(
    trellis: "Trellis",
    partitions: Iterable[StablePartitionResult],
    *,
    bridge_ids: Optional[Iterable["BridgeId"]] = None,
) -> BridgeClassTable:
    """
    Class a trellis's bridges by the pair of partition elements they connect.

    For each bridge id, the row at each end (:func:`~.StablePartition.row_of_end`)
    picks the side, the endpoint's ``manifold_b_key`` picks the stable branch,
    and the partition result for that ``(branch, side)`` names the element that
    owns the endpoint. Two different elements give the bridge's
    :class:`BridgeClass` (oriented anchor outward) and its direction; two equal
    elements make it a loop, which is folded into the class of its first non-loop
    ancestor along the registry's preimage chain (see the module Dev Notes).

    Args:
        trellis: The trellis whose crossings and bridges to read. Normally the
            all-fixed-points trellis, so that nested and heteroclinic bridges
            are covered in one pass.
        partitions: The stable partitions to resolve elements against — every
            per-fixed-point trellis's :attr:`~.Trellis.Trellis.stable_partitions`
            concatenated. Indexed by ``(branch_key, side)``.
        bridge_ids: The bridges to class. Defaults to every non-partial bridge
            of ``trellis``. A loop's ancestor is looked up on the trellis whether
            or not it is among these.

    Returns:
        A :class:`BridgeClassTable` ordered by :func:`class_sort_key`, each
        entry's members sorted by :data:`~tanglepack.numerics.Bridge.BridgeId`.
        Letters and zone keys are left unset (the session fills them).

    Raises:
        ValueError: If two of ``partitions`` cover the same ``(branch, side)``;
            if no partition covers the ``(branch, side)`` a bridge end needs
            (the message names the branch); if a crossing carries no stable
            manifold key; or if a crossing has no ``crossing_sign``.

    Note:
        A crossing that lies on a partitioned branch but is absent from that
        result's ``element_of_intersection`` is an AssertionError, not a
        ValueError: a partition covers every crossing of its branch by
        construction, so a miss is a broken invariant rather than a caller
        error.
    """
    by_branch_side = _index_partitions(partitions)
    ids = (
        list(bridge_ids)
        if bridge_ids is not None
        else [bridge.id for bridge in trellis.bridges if bridge.id is not None]
    )
    fixed_points = trellis.fixed_points

    def ends(bridge_id: "BridgeId") -> tuple[ElementRef, ElementRef]:
        return (
            _element_at_end(trellis, by_branch_side, bridge_id, "first"),
            _element_at_end(trellis, by_branch_side, bridge_id, "second"),
        )

    grouped: dict[BridgeClass, list[BridgeMember]] = {}
    for bridge_id in ids:
        x, y = ends(bridge_id)
        bridge_class, direction = oriented_class(x, y, fixed_points)
        if direction != 0:
            grouped.setdefault(bridge_class, []).append(BridgeMember(bridge_id, direction))
            continue
        ancestor = _non_loop_ancestor(trellis, bridge_id, ends)
        if ancestor is None:
            logger.warning(
                "loop bridge %s (both ends in element %s) has no registered "
                "non-loop preimage; classed on its own as an inert loop",
                bridge_id,
                x.label,
            )
            grouped.setdefault(bridge_class, []).append(
                BridgeMember(bridge_id, 0, loop_element=x)
            )
            continue
        ancestor_id, ancestor_class = ancestor
        grouped.setdefault(ancestor_class, []).append(
            BridgeMember(bridge_id, 0, folded_from=ancestor_id, loop_element=x)
        )
        logger.debug(
            "loop bridge %s folded into class %s via ancestor %s",
            bridge_id,
            ancestor_class.label,
            ancestor_id,
        )

    entries = [
        BridgeClassEntry(
            bridge_class, sorted(members, key=lambda member: member.bridge_id)
        )
        for bridge_class, members in sorted(
            grouped.items(), key=lambda item: class_sort_key(item[0], fixed_points)
        )
    ]
    logger.debug("classed %d bridges into %d classes", len(ids), len(entries))
    return BridgeClassTable(entries)


def _non_loop_ancestor(
    trellis: "Trellis",
    bridge_id: "BridgeId",
    ends,
) -> Optional[tuple["BridgeId", BridgeClass]]:
    """
    Walk a loop's preimage chain to the first bridge that connects two elements.

    Registry lookups only: each step takes the ``-1`` iterate of both endpoints
    (:meth:`~.Trellis.Trellis.iterate`) and the bridge spanning them
    (:meth:`~.Trellis.Trellis.bridge_between`). A bridge's preimage is at most
    one bridge, because the forward map can only add crossings to an arc.

    Returns:
        ``(ancestor id, its class)`` or ``None`` when the chain breaks (a missing
        iterate, no bridge over the preimage pair, a partial one) or cycles.
    """
    fixed_points = trellis.fixed_points
    current = bridge_id
    visited = {bridge_id}
    while True:
        first = trellis.iterate(current[0], -1)
        second = trellis.iterate(current[1], -1)
        if first is None or second is None:
            return None
        bridge = trellis.bridge_between(first, second)
        if bridge is None or bridge.id is None or bridge.id in visited:
            return None
        visited.add(bridge.id)
        x, y = ends(bridge.id)
        bridge_class, direction = oriented_class(x, y, fixed_points)
        if direction != 0:
            return bridge.id, bridge_class
        current = bridge.id


def _index_partitions(
    partitions: Iterable[StablePartitionResult],
) -> dict[tuple["ManifoldKey", Side], StablePartitionResult]:
    """Index partition results by ``(branch_key, side)``, rejecting duplicates."""
    index: dict[tuple["ManifoldKey", Side], StablePartitionResult] = {}
    for result in partitions:
        key = (result.branch_key, result.side)
        if key in index:
            raise ValueError(
                f"two partitions were given for branch {result.branch_key[1:]} "
                f"side {result.side!r}; a (branch, side) names one partition"
            )
        index[key] = result
    return index


def _element_at_end(
    trellis: "Trellis",
    by_branch_side: dict[tuple["ManifoldKey", Side], StablePartitionResult],
    bridge_id: "BridgeId",
    endpoint: Endpoint,
) -> ElementRef:
    """The element a bridge sits against at one of its ends, on that end's row."""
    row = row_of_end(trellis, bridge_id, endpoint)
    intersection_id = bridge_id[endpoint_index(endpoint)]
    branch_key = trellis.intersection(intersection_id).manifold_b_key
    if branch_key is None:
        raise ValueError(
            f"crossing {intersection_id} of bridge {bridge_id} carries no "
            "stable manifold key, so no partition can name its element"
        )
    result = by_branch_side.get((branch_key, row))
    if result is None:
        raise ValueError(
            f"no {row!r} partition covers stable branch {branch_key[1:]}, which "
            f"carries crossing {intersection_id} of bridge {bridge_id}; "
            "partition that branch (or pass its trellis's stable_partitions)"
        )
    element_id = result.element_of_intersection.get(intersection_id)
    assert element_id is not None, (
        f"crossing {intersection_id} lies on partitioned branch "
        f"{branch_key[1:]} side {row!r} but no element of that partition owns "
        "it; a partition covers every crossing of its branch"
    )
    return result.ref(element_id)

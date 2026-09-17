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
  non-loop ancestor is found. A loop whose chain is broken (no registered
  preimage, or a partial one) stands alone as a ``BridgeClass(x, x)`` entry,
  flagged :attr:`BridgeClass.is_loop` and logged.
* INERT = MAPS TO NO ACTIVE CLASS (user definition, 2026-09-16). Each non-loop
  member's forward image is read from the iterate table alone
  (:func:`_image_chain`: the ``+1`` iterates of its two endpoints and every
  registered crossing between them on the image's unstable branch, as a chain
  of consecutive pairs). Each image pair is classed by the same row/element
  rule as a bridge — it needs only crossing ids, so a pair that NO ``Bridge``
  object spans (the parent manifold ends before the image crossing, as on the
  k=2.8 blasted tangle) still has a class or is still a loop. A class is INERT
  when it has at least one resolvable image and every resolvable image is a
  loop or lies in an inert class; it is ACTIVE when any resolvable image lies
  in an active class, its own class included (the anchor bridge always maps
  over itself). Computed as a fixed point over the table
  (:func:`_resolve_inertness`); a class with no resolvable image at all stays
  active for lack of evidence, and loop members contribute no images of their
  own (a loop's image is a loop) but do count as loop evidence for the class
  they folded into, so "holds a loop" still implies inert. On k=10 the
  exterior class {(1,12), (9,11), loop (15,6)} is inert although (1,12) has no
  registered image; on k=2.8 the class of (1,2) is inert through the virtual
  loop (5,6) and the class {(2,5), loop (6,7)} through its folded loop, while
  the anchor class stays active by mapping over itself. Inertness is still
  evidence-based: grow the trellis (or blast) until the images are registered.
  A loop produced by a too-coarse partition (an interior bridge whose ends the
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
        inert: True when the class maps to no active class (see the module Dev
            Notes); set by :func:`bridge_classes` from the image evidence below.
        image_classes: The classes the resolvable forward images of the non-loop
            members land in (this class included when a member maps over
            itself), in first-seen order without repeats. A class listed here
            need not have an entry in the table: an image pair no bridge spans
            still names one.
        image_loops: The image pairs that are loops — both image ends in one
            element — as ``(image pair, element)``, folded loop members included.
            Any such pair is inert evidence.
        unresolved: The non-loop members whose forward image the iterate table
            cannot resolve (an endpoint without a registered ``+1`` iterate).
            They contribute no evidence either way.
        letter: The class's symbol letter (``"a"``, ``"b"``, ...), set by the
            session alphabet for active classes; ``None`` when unlettered.
        zone_key: The ``(fixed_point, branch_index)`` key of the resonance zone
            the class lies in, set by the session; ``None`` when it lies outside
            every zone or no zones are defined.
    """

    bridge_class: BridgeClass
    members: list[BridgeMember] = field(default_factory=list)
    inert: bool = False
    image_classes: list[BridgeClass] = field(default_factory=list)
    image_loops: list[tuple["BridgeId", ElementRef]] = field(default_factory=list)
    unresolved: list["BridgeId"] = field(default_factory=list)
    letter: Optional[str] = None
    zone_key: Optional[tuple] = None

    @property
    def active(self) -> bool:
        """The opposite of :attr:`inert`."""
        return not self.inert

    @property
    def has_image_evidence(self) -> bool:
        """True when at least one member's forward image was resolved."""
        return bool(self.image_classes or self.image_loops)

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

    def describe(self, names: Optional[dict[BridgeClass, str]] = None) -> str:
        """
        One line: name, elements, activity, zone, members with symbols, and images.

        Args:
            names: How to print the classes in :attr:`image_classes` — normally
                the table's ``{class: entry.name}`` so lettered classes show
                their letter. A class missing from it prints its label.

        Returns:
            The description line.
        """
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
        line = f"{head} [{state}, {zone}] {len(self.members)} bridge(s): " + ", ".join(parts)

        names = names or {}
        images = [
            names.get(image_class, image_class.label) for image_class in self.image_classes
        ] + [f"{pair} loop in {element}" for pair, element in self.image_loops]
        if images:
            line += "; images: " + ", ".join(images)
        else:
            line += "; images: none registered"
        if self.unresolved:
            line += f" ({len(self.unresolved)} member(s) with no registered image)"
        return line


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
        """The entries that map to some active class (or have no image evidence)."""
        return [entry for entry in self.entries if entry.active]

    @property
    def inert(self) -> list[BridgeClassEntry]:
        """The entries that map to no active class."""
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
        names = {entry.bridge_class: entry.name for entry in self.entries}
        return "\n".join(
            [head, *(f"  {entry.describe(names)}" for entry in self.entries)]
        )

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
    Each class is then marked inert or active from its members' forward images
    (:func:`_image_chain`, :func:`_resolve_inertness`): inert when it maps only
    to loops and inert classes, active when it maps to an active class (itself
    included) or when no image is registered at all.

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
        entry's members sorted by :data:`~tanglepack.numerics.Bridge.BridgeId`
        and its ``inert`` flag and image evidence filled in. Letters and zone
        keys are left unset (the session fills them).

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
    for entry in entries:
        _collect_image_evidence(trellis, entry, ends, fixed_points)
    _resolve_inertness(entries)
    logger.debug(
        "classed %d bridges into %d classes (%d inert)",
        len(ids),
        len(entries),
        sum(entry.inert for entry in entries),
    )
    return BridgeClassTable(entries)


def _image_chain(trellis: "Trellis", bridge_id: "BridgeId") -> Optional[list["BridgeId"]]:
    """
    The consecutive crossing pairs tiling a bridge's forward image, by lookup only.

    The image of the arc between crossings ``a`` and ``b`` is the arc between
    ``f(a)`` and ``f(b)``, read from the iterate table
    (:meth:`~.Trellis.Trellis.iterate`); the crossings on it are every
    registered crossing of the image's unstable branch whose unstable
    canonical distance lies between the two. No ``Bridge`` object is required:
    an image pair the trellis never cut (the grown manifold stops short of it)
    is returned all the same, because its class is readable from its ids.

    Args:
        trellis: The Trellis whose registry and branches to read.
        bridge_id: The bridge whose image is wanted.

    Returns:
        The image's consecutive pairs in increasing unstable canonical
        distance (each ordered like a
        :data:`~tanglepack.numerics.Bridge.BridgeId`), or ``None`` when an
        endpoint has no registered ``+1`` iterate, the two images carry no
        common unstable branch key, or that branch is not in the trellis.
    """
    images = [trellis.iterate(endpoint, 1) for endpoint in bridge_id]
    if any(image is None for image in images):
        return None
    ia, ib = (trellis.intersection(image) for image in images)
    key = ia.manifold_a_key
    if key is None or ib.manifold_a_key != key:
        logger.warning(
            "image of bridge %s has endpoints %s on unstable branches %s and %s; "
            "cannot place its image arc",
            bridge_id,
            tuple(images),
            None if key is None else key[1:],
            None if ib.manifold_a_key is None else ib.manifold_a_key[1:],
        )
        return None
    branch = trellis.branch(key)
    if branch is None:
        return None

    lo, hi = sorted((ia.unstable_cdist, ib.unstable_cdist))
    on_arc = [
        iid
        for iid in branch.intersection_ids
        if lo <= trellis.intersection(iid).unstable_cdist <= hi
    ]
    for image in images:
        if image not in on_arc:
            on_arc.append(image)
    on_arc.sort(key=lambda iid: trellis.intersection(iid).unstable_cdist)
    return list(zip(on_arc, on_arc[1:]))


def _collect_image_evidence(
    trellis: "Trellis",
    entry: BridgeClassEntry,
    ends,
    fixed_points: Sequence["FixedPoint"],
) -> None:
    """
    Fill an entry's ``image_classes``, ``image_loops`` and ``unresolved``.

    Every non-loop member's image chain (:func:`_image_chain`) is classed pair
    by pair with the same row/element rule as the members themselves; a pair
    with two different elements adds its class, a pair with one element adds a
    loop. Loop members add themselves as loop evidence and no images of their
    own (a loop's image is a loop). An image pair whose end lies on a branch no
    partition covers is logged and skipped.
    """
    seen_classes: set[BridgeClass] = set()
    seen_loops: set["BridgeId"] = set()

    def add_loop(pair: "BridgeId", element: ElementRef) -> None:
        if pair not in seen_loops:
            seen_loops.add(pair)
            entry.image_loops.append((pair, element))

    for member in entry.members:
        if member.is_loop:
            if member.loop_element is not None:
                add_loop(member.bridge_id, member.loop_element)
            continue
        chain = _image_chain(trellis, member.bridge_id)
        if chain is None:
            entry.unresolved.append(member.bridge_id)
            continue
        for pair in chain:
            try:
                x, y = ends(pair)
            except ValueError as error:
                logger.warning(
                    "image pair %s of bridge %s cannot be classed (%s); ignoring it",
                    pair,
                    member.bridge_id,
                    error,
                )
                continue
            image_class, direction = oriented_class(x, y, fixed_points)
            if direction == 0:
                add_loop(pair, x)
            elif image_class not in seen_classes:
                seen_classes.add(image_class)
                entry.image_classes.append(image_class)


def _resolve_inertness(entries: Iterable[BridgeClassEntry]) -> None:
    """
    Mark the inert entries: those that map to no active class.

    A fixed-point iteration from "every class active": an entry becomes inert
    once it has image evidence and every class in its ``image_classes`` is
    already inert (a loop image is never active; an image class with no entry
    of its own — a virtual, uncut pair between two elements — is unknown and
    therefore counts as active). A class among its own images keeps itself
    active. An entry whose class is itself a loop (an unresolved loop standing
    alone) is inert outright.
    """
    entries = list(entries)
    inert: set[BridgeClass] = {
        entry.bridge_class for entry in entries if entry.bridge_class.is_loop
    }
    changed = True
    while changed:
        changed = False
        for entry in entries:
            if entry.bridge_class in inert or not entry.has_image_evidence:
                continue
            if all(image_class in inert for image_class in entry.image_classes):
                inert.add(entry.bridge_class)
                changed = True
    for entry in entries:
        entry.inert = entry.bridge_class in inert


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

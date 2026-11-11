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
* LOOPS FOLD FORWARD (the author's rule). A
  bridge whose two ends land in the SAME element is not a class of its own: it
  is the forward image of a bridge that did connect two elements and has been
  pushed into one of them, and it will keep mapping into that element. Such a
  loop joins the class of its IMAGE-CHAIN PARENT — the one bridge whose
  :func:`_image_chain` holds it (the map is a bijection, so the images of
  disjoint arcs are disjoint and a pair has at most one parent; asserted) —
  climbing through parents that are loops themselves until one connects two
  elements. Lookups only, never the map. Walking PREIMAGES instead would break
  at every fold-born crossing, which has no ``-1`` link. A loop no image
  chain reaches is an ORPHAN: logged once, recorded on
  :attr:`BridgeClassTable.orphans` and left unclassed, never lettered.
* EVIDENCE vs ALPHABET. :func:`bridge_classes` classes every
  bridge of the trellis: that is the EVIDENCE table (the minimal trellis
  reads it). The ALPHABET the session letters, the words and the
  transition matrix use is its restriction to the classes with a member the
  minimal trellis keeps (``MinimalTrellis.alphabet``,
  :meth:`BridgeClassTable.restricted_to`); a blast child is image evidence
  only.
* INERT = TRANSIENT IN THE WORD GRAPH (design note §A.8). A class
  is inert iff no cycle is reachable from it over the class-level words: its
  word is empty or names only inert classes, and an unresolved class (no
  word) is active, as is every class whose word names one. The verdict is the
  dynamics' (``SymbolicDynamics.inert_classes``), so :attr:`BridgeClassEntry.inert`
  is ``None`` here until the dynamics sets it. The minimal trellis decides
  before the words with the author's sufficient condition: a
  bridge carrying a ZEROTH hole is inert, since its lobe holds one marker
  whose image is none, so its word is empty.
* Partitions are passed in rather than read off the trellis because they do not
  live where the bridges do: a trellis holds a bridge under its UNSTABLE fixed
  point, while the partition owning that bridge's endpoints lives on the trellis
  of the STABLE branch's fixed point. Running this on the all-fixed-points
  trellis with every per-fixed-point trellis's partitions is what covers the
  nested and heteroclinic cases in one pass.
* Ordering is deliberate and total: the table is grouped by TANGLE first --
  the classes of ``trellis.fixed_points[0]``, then of the next fixed point,
  then the classes connecting two tangles (heteroclinic) -- and within a group
  built in ascending order of each class's ``min_unstable_cdist`` (the
  smallest unstable canonical distance of a member's first crossing, so the
  class holding a tangle's anchor bridge comes first in its group; author's
  rules; unstable cdists of
  different saddles are not comparable, so they are never interleaved),
  ties broken by :func:`class_sort_key`, and each class's members are sorted
  by :data:`~tanglepack.numerics.Bridge.BridgeId`, so two runs of the same
  trellis enumerate the classes identically. Letters, legends and
  transition-graph slots all follow this one order.
* Letters (``a``, ``b``, ...) and resonance-zone membership are session-level
  annotations (the loom owns zones and the persistent alphabet); the slots for
  them live on :class:`BridgeClassEntry` so one object carries the whole record.
* Refined symbols — subdividing one class into several symbols (``a_1, a_2``)
  that all still belong to it — live in ``SymbolicDynamics.py``
  (``RefinedClass``, ``SymbolicDynamics.refined`` / ``member_refinement``),
  derived from the walked itineraries of all classes, never from this table.
  :class:`BridgeClass` stays the coarse key: a refinement groups an entry's
  members more finely and changes no class identity.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field, replace
from typing import Iterable, Iterator, Mapping, Optional, Sequence, TYPE_CHECKING

from .TopologyResults import (
    ElementRef,
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

    ``source`` is the element nearer its branch anchor (smaller element id, on
    one branch or across branches; :func:`element_sort_key` breaks a tie of
    ids, see :func:`anchor_outward_key`), so the
    pair ``{X, Y}`` has exactly one representation and the class hashes as one
    symbol whichever way its bridges run. Frozen, so it is usable as a dict key
    and as the key of a session alphabet.

    Attributes:
        source: The anchor-nearer element; the start of the class's direction.
        target: The other element. Equal to ``source`` only for the pseudo-class
            :func:`oriented_class` returns for a loop pair; never a table entry.
        lobe: Empty unless the pair holds several classes. A pair with both
            elements on one stable branch names one class only when every arc
            between them closes with the stable segment into a lobe holding the
            same markers; where the lobes differ (a cone point lets an arc pass
            either side of an island; :mod:`.LobeSplit`), each class is the pair
            plus the markers its lobes hold.
    """

    source: ElementRef
    target: ElementRef
    lobe: frozenset = frozenset()

    @property
    def is_loop(self) -> bool:
        """True when both elements coincide: the pseudo-class of a loop pair."""
        return self.source == self.target

    @property
    def label(self) -> str:
        """A short, deterministic name: the two element labels, ``<->``-joined."""
        if self.is_loop:
            return f"{self.source.label} loop"
        label = f"{self.source.label} <-> {self.target.label}"
        if self.lobe:
            label += " {" + ", ".join(sorted(marker_text(m) for m in self.lobe)) + "}"
        return label

    def with_lobe(self, lobe: Optional[frozenset]) -> "BridgeClass":
        """The same pair, told apart by ``lobe`` (unchanged for None)."""
        return self if lobe is None else replace(self, lobe=frozenset(lobe))


def marker_text(marker: tuple) -> str:
    """A short name for a dual-graph marker: ``hole a.b^i`` or ``orbit p.i``."""
    if marker and marker[0] == "hole":
        origin, iterate = marker[1], marker[2]
        return f"hole {'.'.join(map(str, origin))}^{iterate}"
    if marker and marker[0] == "orbit":
        return f"orbit {marker[1]}.{marker[2]}"
    return str(marker)

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
        folded_from: For a loop, the first non-loop bridge up its image-chain
            parents, whose class it joined. ``None`` for a direct member.
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
        inert: True when the class is transient in the word graph (see the
            module Dev Notes); ``None`` until the symbolic dynamics sets it.
        letter: The class's symbol letter (``"a"``, ``"b"``, ...), set by the
            session alphabet for active classes; ``None`` when unlettered.
        zone_key: The key of the resonance zone the class lies in (its fixed
            point, one zone per fixed point), set by the session; ``None`` when it lies outside
            every zone or no zones are defined.
        min_unstable_cdist: The smallest unstable canonical distance of a
            member's first crossing — the table's sort key within a tangle
            (the anchor class has 0); ``inf`` for an entry assembled by hand.
        tangle: The position in ``trellis.fixed_points`` of the one fixed
            point the class belongs to (every member's unstable branch and
            both elements' stable branches on it), or None for a class that
            connects two tangles (heteroclinic). The table groups by it.
    """

    bridge_class: BridgeClass
    members: list[BridgeMember] = field(default_factory=list)
    inert: Optional[bool] = None
    letter: Optional[str] = None
    zone_key: Optional["FixedPoint"] = None
    min_unstable_cdist: float = float("inf")
    tangle: Optional[int] = None

    @property
    def active(self) -> bool:
        """The opposite of :attr:`inert`."""
        return not self.inert

    @property
    def spanned(self) -> bool:
        """False for a class a word names that no computed bridge spans yet:
        lettered and walked like any other, informational only."""
        return bool(self.members)

    @property
    def bridge_ids(self) -> list["BridgeId"]:
        """The member bridge ids, in member order."""
        return [member.bridge_id for member in self.members]

    @property
    def loops(self) -> list[BridgeMember]:
        """The loop members (each folded in from its image-chain parent)."""
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
        """
        One line: name, elements, activity, zone and members with symbols.

        Returns:
            The description line.
        """
        zone = "no zone" if self.zone_key is None else zone_label(self.zone_key)
        parts = []
        for member in self.members:
            if member.is_loop:
                parts.append(
                    f"{member.bridge_id} loop in {member.loop_element} "
                    f"(folded from {member.folded_from})"
                )
            elif self.letter is not None:
                parts.append(f"{member.bridge_id} {self.symbol(member.bridge_id)}")
            else:
                parts.append(f"{member.bridge_id} {'+' if member.direction > 0 else '-'}")
        state = {None: "undecided", True: "inert", False: "active"}[self.inert]
        head = (
            f"{self.letter}: {self.bridge_class.label}"
            if self.letter is not None
            else self.bridge_class.label
        )
        return f"{head} [{state}, {zone}] {len(self.members)} bridge(s): " + ", ".join(parts)


@dataclass
class BridgeClassTable:
    """
    Bridge classes, grouped by tangle (``entry.tangle``, heteroclinic last) and
    in ascending ``min_unstable_cdist`` order within a group (ties by
    :func:`class_sort_key`).

    Built by :func:`bridge_classes` over every bridge (the evidence table), or
    by :meth:`restricted_to` (the alphabet).

    Attributes:
        entries: The classes, one :class:`BridgeClassEntry` each.
        orphans: The loop bridges no image chain reaches, left unclassed (see
            the module Dev Notes); empty on a restricted table.
    """

    entries: list[BridgeClassEntry] = field(default_factory=list)
    orphans: list["BridgeId"] = field(default_factory=list)

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
        """The entries not known to be inert (every entry before the dynamics runs)."""
        return [entry for entry in self.entries if entry.active]

    @property
    def inert(self) -> list[BridgeClassEntry]:
        """The entries the dynamics found transient."""
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

    def symbol(self, bridge_id: "BridgeId") -> str:
        """The itinerary symbol of a bridge; see :meth:`BridgeClassEntry.symbol`."""
        return self.entry_of(bridge_id).symbol(bridge_id)

    def restricted_to(self, bridge_ids: Iterable["BridgeId"]) -> "BridgeClassTable":
        """
        The classes with a member among some bridges: the alphabet of a bridge set.

        Args:
            bridge_ids: The bridges, in either endpoint order (normally
                ``MinimalTrellis.kept_bridge_ids``).

        Returns:
            A table of the SAME entry objects, in table order, members kept
            whole (a class's other members are its image evidence).
        """
        wanted = {frozenset(bid) for bid in bridge_ids}
        return BridgeClassTable(
            [
                entry
                for entry in self.entries
                if any(frozenset(bid) in wanted for bid in entry.bridge_ids)
            ]
        )

    def as_dict(self) -> dict[BridgeClass, list["BridgeId"]]:
        """The table as a plain ``{class: sorted member ids}`` mapping."""
        return {entry.bridge_class: list(entry.bridge_ids) for entry in self.entries}

    def describe(self) -> str:
        """A multi-line report, one :meth:`BridgeClassEntry.describe` line per class."""
        head = (
            f"{len(self.entries)} bridge class(es): {len(self.active)} active, "
            f"{len(self.inert)} inert"
        )
        lines = [head, *(f"  {entry.describe()}" for entry in self.entries)]
        if self.orphans:
            lines.append(f"  orphan loop(s), unclassed: {self.orphans}")
        return "\n".join(lines)


def zone_label(zone_key: "FixedPoint") -> str:
    """
    A readable name for a resonance-zone key (its fixed point).

    Args:
        zone_key: The key, as stored on :attr:`BridgeClassEntry.zone_key`.

    Returns:
        ``"zone p{period}"``.
    """
    return f"zone p{getattr(zone_key, 'period', '?')}"


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


def anchor_outward_key(
    ref: ElementRef, fixed_points: Sequence["FixedPoint"]
) -> tuple:
    """
    The key that orients a class anchor outward, across branches too.

    Args:
        ref: One end of a class.
        fixed_points: The ordering context (see :func:`element_sort_key`).

    Returns:
        ``(element_id,) + element_sort_key(ref)``: anchor-nearness first, so a
        class between two branches of one orbit starts at its anchor-nearer
        element whichever branch that is (the classes of a period-k orbit are
        then oriented alike on every branch, and a word does not pick up an
        inverse from the orbit index wrapping), then the run-stable order.

    Note:
        On one branch and side this is the element-id order, exactly what
        :func:`element_sort_key` gives, so same-branch classes are unchanged.
    """
    return (ref.element_id,) + element_sort_key(ref, fixed_points)


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
        The source is the anchor-nearer end (:func:`anchor_outward_key`).
    """
    if x == y:
        return BridgeClass(x, x), 0
    if anchor_outward_key(x, fixed_points) <= anchor_outward_key(y, fixed_points):
        return BridgeClass(x, y), +1
    return BridgeClass(y, x), -1


def bridge_classes(
    trellis: "Trellis",
    partitions: Iterable[StablePartitionResult],
    lobes: Optional[Mapping["BridgeId", frozenset]] = None,
) -> BridgeClassTable:
    """
    Class a trellis's bridges by the pair of partition elements they connect.

    The result is the EVIDENCE table: every bridge's class. The alphabet is its
    restriction to the minimal trellis (``MinimalTrellis.alphabet``).

    Each end of a bridge id is read by the homotopy family's
    :meth:`~.PartitionFamily.PartitionFamily.owner_at_end`: the row at that end
    (:func:`~.StablePartition.row_of_end`) picks the side,
    :meth:`~.Trellis.Trellis.stable_key_at` picks the stable branch, and
    the partition result for that ``(branch, side)`` names the element that
    owns the endpoint. Two different elements give the bridge's
    :class:`BridgeClass` (oriented anchor outward) and its direction; two equal
    elements make it a loop, which is folded FORWARD into the class of the first
    non-loop bridge up its image-chain parents; a loop with none is an orphan
    (see the module Dev Notes).
    Inertness is left ``None``: it is a verdict on the words (the module Dev
    Notes), which the symbolic dynamics sets.

    Args:
        trellis: The trellis whose crossings and bridges to read. Normally the
            all-fixed-points trellis, so that nested and heteroclinic bridges
            are covered in one pass.
        partitions: The stable partitions to resolve elements against — every
            per-fixed-point trellis's :attr:`~.Trellis.Trellis.stable_partitions`
            concatenated. Indexed by ``(branch_key, side)``.
        lobes: ``{bridge id: lobe}`` for the members of a pair that holds
            several classes (:func:`~.LobeSplit.split_by_lobes`); each such
            bridge's class is its pair with that lobe. A folded loop takes its
            parent's.

    Returns:
        A :class:`BridgeClassTable` ordered by each class's smallest member
        first-crossing unstable canonical distance (``min_unstable_cdist``,
        ties by :func:`class_sort_key`), each entry's members sorted by
        :data:`~tanglepack.numerics.Bridge.BridgeId`, and the orphan loops
        listed on :attr:`BridgeClassTable.orphans`. Letters, verdicts and zone
        keys are left unset (the dynamics and the session fill them).

    Raises:
        ValueError: If two of ``partitions`` cover the same ``(branch, side)``;
            if no partition covers the ``(branch, side)`` a bridge end needs
            (the message names the branch); if a covering result does not own
            a crossing; or if a crossing has no ``crossing_sign``.

    Note:
        A crossing a covering result does not own breaks an invariant the
        partition asserts when it is built, but here it is
        :meth:`~.PartitionFamily.PartitionFamily.owner_of_intersection`'s
        ValueError, as for every other caller. An image pair held by two
        bridges' image chains is an AssertionError (the map is a bijection).
    """
    # PartitionFamily imports this module (``_image_chain``), so not at module level.
    from .PartitionFamily import HomotopyPartition

    family = HomotopyPartition.from_results(partitions, trellis=trellis)
    ids = [bridge.id for bridge in trellis.bridges if bridge.id is not None]
    fixed_points = trellis.fixed_points
    chains = {bid: _image_chain(trellis, bid) for bid in ids}
    image_parent: dict[frozenset, "BridgeId"] = {}
    for bid, chain in chains.items():
        for pair in chain or ():
            assert image_parent.setdefault(frozenset(pair), bid) == bid, (
                f"image pair {pair} lies in the image chains of both {bid} and "
                f"{image_parent[frozenset(pair)]}; the images of disjoint arcs are disjoint"
            )

    def ends(bridge_id: "BridgeId") -> tuple[ElementRef, ElementRef]:
        return (
            family.owner_at_end(trellis, bridge_id, "first"),
            family.owner_at_end(trellis, bridge_id, "second"),
        )

    lobes = lobes or {}
    grouped: dict[BridgeClass, list[BridgeMember]] = {}
    orphans: list["BridgeId"] = []
    for bridge_id in ids:
        x, y = ends(bridge_id)
        bridge_class, direction = oriented_class(x, y, fixed_points)
        if direction != 0:
            bridge_class = bridge_class.with_lobe(lobes.get(tuple(bridge_id)))
            grouped.setdefault(bridge_class, []).append(BridgeMember(bridge_id, direction))
            continue
        folded = _non_loop_parent(bridge_id, image_parent, ends, fixed_points)
        if folded is None:
            orphans.append(bridge_id)
            continue
        parent_id, parent_class = folded
        parent_class = parent_class.with_lobe(lobes.get(tuple(parent_id)))
        grouped.setdefault(parent_class, []).append(
            BridgeMember(bridge_id, 0, folded_from=parent_id, loop_element=x)
        )
    if orphans:
        logger.warning(
            "loop bridge(s) %s lie in no non-loop bridge's image chain; left "
            "unclassed and never lettered (D2)",
            orphans,
        )

    def min_cdist(members: list[BridgeMember]) -> float:
        return min(
            float(trellis.intersection(member.bridge_id[0]).unstable_cdist)
            for member in members
        )

    def tangle_of(bridge_class: BridgeClass, members: list[BridgeMember]) -> Optional[int]:
        owners = {id(bridge_class.source.fixed_point), id(bridge_class.target.fixed_point)}
        for member in members:
            key = trellis.intersection(member.bridge_id[0]).manifold_a_key
            if key is not None:
                owners.add(id(key[0]))
        if len(owners) != 1:
            return None
        (owner,) = owners
        return next((i for i, fp in enumerate(fixed_points) if id(fp) == owner), None)

    def table_key(item) -> tuple:
        bridge_class, members = item
        tangle = tangle_of(bridge_class, members)
        return (
            len(fixed_points) if tangle is None else tangle,
            min_cdist(members),
            class_sort_key(bridge_class, fixed_points),
        )

    entries = [
        BridgeClassEntry(
            bridge_class,
            sorted(members, key=lambda member: member.bridge_id),
            min_unstable_cdist=min_cdist(members),
            tangle=tangle_of(bridge_class, members),
        )
        for bridge_class, members in sorted(grouped.items(), key=table_key)
    ]
    logger.debug(
        "classed %d bridges into %d classes (%d orphan loops)", len(ids), len(entries), len(orphans)
    )
    return BridgeClassTable(entries, orphans)


def _image_chain(trellis: "Trellis", bridge_id: "BridgeId") -> Optional[list["BridgeId"]]:
    """
    The consecutive crossing pairs tiling a bridge's forward image, by lookup only.

    The image of the arc between crossings ``a`` and ``b`` is the arc between
    ``f(a)`` and ``f(b)``, read from the iterate table
    (:meth:`~.Trellis.Trellis.iterate`), on the unstable branch one map step
    after the arc's own (an anchor end is on every branch of its point); the
    crossings on it are every registered crossing of that branch whose
    unstable canonical distance lies between the two. No ``Bridge`` object is required:
    an image pair the trellis never cut (the grown manifold stops short of it)
    is returned all the same, because its class is readable from its ids.

    Args:
        trellis: The Trellis whose registry and branches to read.
        bridge_id: The bridge whose image is wanted.

    Returns:
        The image's consecutive pairs in increasing unstable canonical
        distance (each ordered like a
        :data:`~tanglepack.numerics.Bridge.BridgeId`), or ``None`` when an
        endpoint's ``+1`` iterate is no crossing of the trellis (unregistered,
        or past its pip-orbit cut: the image leaves T(q0), no evidence either
        way), or the advanced branch is not in the trellis.

    Raises:
        AssertionError: If an image end lies off the advanced branch (a table
            link always passes the map gate onto ``advance_key``).
    """
    images = [trellis.iterate(endpoint, 1) for endpoint in bridge_id]
    if any(image not in trellis for image in images):
        return None
    ia, ib = (trellis.intersection(image) for image in images)
    own_key = trellis.unstable_key_of(bridge_id)
    key = own_key[0].advance_key(own_key, 1)
    # A table link passes the map gate onto advance_key, so an image end
    # lies on the advanced branch (an anchor lies on every branch of its point).
    for image in (ia, ib):
        assert image.is_anchor or image.manifold_a_key == key, (
            f"image of bridge {bridge_id} ends off its image branch {key[1:]}: {images}"
        )
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


def _non_loop_parent(
    bridge_id: "BridgeId",
    image_parent: dict[frozenset, "BridgeId"],
    ends,
    fixed_points: Sequence["FixedPoint"],
) -> Optional[tuple["BridgeId", BridgeClass]]:
    """
    Climb a loop's image-chain parents to the first bridge that connects two elements.

    Returns:
        ``(parent id, its class)``, or ``None`` when no image chain holds the
        current bridge or the climb cycles (an orphan).
    """
    visited = {frozenset(bridge_id)}
    current = image_parent.get(frozenset(bridge_id))
    while current is not None and frozenset(current) not in visited:
        visited.add(frozenset(current))
        bridge_class, direction = oriented_class(*ends(current), fixed_points)
        if direction:
            return current, bridge_class
        current = image_parent.get(frozenset(current))
    return None


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

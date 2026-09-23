"""
Partition families: the homotopy partition and its iterated refinement.

A *family* is one partition of the stable manifold, held as one
:class:`~tanglepack.topology.TopologyResults.StablePartitionResult` per
``(branch_key, side)``. :class:`PartitionFamily` is the base every kind
inherits: indexed access, midpoint ownership, a structural signature for cache
keys, and a report. Two kinds live here:

* :class:`HomotopyPartition` — the partition the holes define (the algorithm
  in :mod:`~tanglepack.topology.StablePartition`, unchanged; this class only
  wraps its results).
* :class:`IteratedHomotopyPartition` — the homotopy partition further cut by
  the image bridges of a :class:`~tanglepack.topology.MinimalTrellis.MinimalTrellis`.

A third kind is expected later and will inherit :class:`PartitionFamily` too.

The iterated cut, in bridges
----------------------------
Each NEW bridge of the minimal trellis (an image pair that is not itself a hole
bridge) lies on one side of the stable manifold — its row, from
:func:`~tanglepack.topology.StablePartition.row_of_end` — and cuts the
partition on THAT side only, by the FOOTPRINT rule (author, 2026-09-21). The
image ``f(H)`` of a hole ``H`` meets the stable manifold along the stretch
between the two END crossings of its image chain (the image of the hole's own
stretch) and along one CHORD per pair of consecutive interior crossings (where
the chain crosses over and back). Those footprint stretches are OPEN at both
ends, like a hole's stretch, and the elements outside them are CLOSED there:
a homotopy element ``[ ]`` cut by one lobe reads ``[ ] ( ) [ ]``. Per
crossing: pair the chain's crossings ``x_0 .. x_m`` (unstable order) as
``(x_0, x_m), (x_1, x_2), (x_3, x_4), ...``; at ``x_i`` the flank toward its
PARTNER is opened on the row of the piece being cut. A boundary the homotopy
partition already has (a hole boundary) is never re-cut; its closedness
stands. The other side is not cut at that crossing. When an image is
subdivided (its pieces lie on alternating sides), each piece cuts its own
side, so the crossing shared by two pieces is cut on both sides, once each,
toward the same partner.

Dev Notes:

* Marks. The homotopy result is turned back into the marks vocabulary of
  ``StablePartition._build_intervals`` — boundary cdists plus the sets of
  boundaries whose OUTWARD / ANCHORWARD interval is open — the new bridges
  add their marks, and the intervals are rebuilt by a local copy of that
  function's tail (:func:`_intervals_from_marks`). The copy exists because
  the original is private to a module under concurrent edit and re-derives
  the anchor and outermost boundaries from the trellis, which a refinement
  must not do; ``tests/test_partition_family.py`` pins that rebuilding every
  homotopy result from its own marks reproduces it exactly. In this
  vocabulary a cut is ONE mark per crossing (the flank toward the partner:
  ``open_anchorward`` when the partner is anchorward, ``open_outward`` when
  it is outward), so a single piece can never produce the ``][`` collision
  the singleton pass would otherwise have to arbitrate; two pieces asking for
  opposite flanks at one crossing (possible only with disagreeing chains) are
  caught before the marks change (``"cut conflict"``).
* Partners need no geometry. A crossing's partner lies on ITS stable branch
  whether or not the piece cutting it is cross-branch (a period-k orbit): the
  pair bounds a connected stable arc of the disc ``f(H)`` — its boundary
  stretch or an interior chord — so one rule cuts every piece, and a partner
  on another branch is a broken chain (WARNING, ``"partner on another
  branch"``, no cut). This replaced the geometric reading through
  ``StablePartition._hole_openings`` for cross-branch bridges (2026-09-21);
  no fixture exercises a cross-branch image bridge beyond the ownership
  invariants.
* The crossings come from the FULL chain (:func:`~.BridgeClass._image_chain`,
  lookup only), not from the minimal trellis's kept pieces, which omit pairs
  without a ``Bridge`` object and would break consecutiveness. A chain with
  an even number of pieces ends on different rows, so its end stretch is not
  a footprint on one row; it is paired as stated anyway (WARNING) and its
  leftover crossing gets no partner (``"no partner"``). A piece shared by two
  chains takes its partners from each; a disagreement is a WARNING and no cut
  at that end (``"partner disagreement"``).
* :class:`~tanglepack.topology.TopologyResults.ElementRef` stays
  ``(branch_key, side, element_id)`` with no family tag. The two families of
  one trellis therefore produce refs that collide by value; never mix refs
  from different families in one structure. Iterated elements point back at
  their homotopy parent through ``PartitionInterval.parent_element_id``.
"""

from __future__ import annotations

import logging
from typing import ClassVar, Iterable, Iterator, Optional, TYPE_CHECKING

from .BridgeClass import _image_chain, _index_partitions
from .StablePartition import (
    _element_of_intersection,
    _elements_at_bridge,
    owns_cdist,
    partition_stable_manifold,
    row_of_end,
)
from .TopologyResults import (
    ElementRef,
    PartitionInterval,
    Side,
    StablePartitionResult,
)

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ..numerics.Bridge import BridgeId
    from ..numerics.Intersection import ManifoldKey
    from .MinimalTrellis import MinimalTrellis
    from .Trellis import Trellis

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: Partner marker for a crossing two image chains pair differently (crossing
#: ids are non-negative, so this never collides with one).
_DISAGREEMENT = -1

#: One boundary set of a partition on one (branch, side): the cdist of every
#: boundary id (None = the anchor sentinel), and the ids whose outward /
#: anchorward interval is open there.
Marks = tuple[dict[Optional[int], float], set[Optional[int]], set[Optional[int]]]


class PartitionFamily:
    """
    One partition of the stable manifold, one result per ``(branch, side)``.

    Attributes:
        results: The results indexed by ``(branch_key, side)``, in insertion
            order.
        trellis: The trellis the results were computed over (supplies the
            canonical-distance tolerance and crossing lookups). Optional for
            a family assembled by hand.
    """

    #: The kind of partition, one word, distinct per subclass.
    kind: ClassVar[str] = "base"

    def __init__(
        self,
        results: Iterable[StablePartitionResult],
        *,
        trellis: Optional["Trellis"] = None,
    ) -> None:
        """
        Args:
            results: The per-``(branch, side)`` results of this family.
            trellis: The trellis they were computed over.

        Raises:
            ValueError: If two results name the same ``(branch_key, side)``.
        """
        self.results: dict[tuple["ManifoldKey", Side], StablePartitionResult] = (
            _index_partitions(results)
        )
        self.trellis = trellis

    # ── access ──────────────────────────────────────────────────────────────

    @property
    def tol(self) -> float:
        """The canonical-distance tolerance (0 without a trellis)."""
        return 0.0 if self.trellis is None else float(self.trellis.registry.cdist_tol)

    def result(self, branch_key: "ManifoldKey", side: Side) -> StablePartitionResult:
        """
        The result for one branch and side.

        Args:
            branch_key: Manifold key of a stable branch.
            side: ``"left"`` or ``"right"``.

        Returns:
            The :class:`~tanglepack.topology.TopologyResults.StablePartitionResult`.

        Raises:
            ValueError: If this family has no result for that pair.
        """
        result = self.results.get((branch_key, side))
        if result is None:
            raise ValueError(
                f"no {side!r} partition covers stable branch {branch_key[1:]} in "
                f"this {self.kind} partition family"
            )
        return result

    def as_list(self) -> list[StablePartitionResult]:
        """The results as a list, in insertion order (what the trellis and
        :func:`~tanglepack.topology.BridgeClass.bridge_classes` take)."""
        return list(self.results.values())

    @property
    def branch_keys(self) -> list["ManifoldKey"]:
        """The distinct stable branch keys this family covers, in order."""
        seen: list["ManifoldKey"] = []
        for key, _side in self.results:
            if key not in seen:
                seen.append(key)
        return seen

    def sides(self, branch_key: "ManifoldKey") -> list[Side]:
        """The sides this family covers on one branch."""
        return [side for key, side in self.results if key == branch_key]

    def ref(self, branch_key: "ManifoldKey", side: Side, element_id: int) -> ElementRef:
        """
        The reference naming one element.

        Args:
            branch_key: Manifold key of the stable branch.
            side: The side.
            element_id: The element's index within that result.

        Returns:
            The :class:`~tanglepack.topology.TopologyResults.ElementRef`.

        Raises:
            ValueError: If no result covers ``(branch_key, side)``.
            IndexError: If the id is not an element of that result.
        """
        return self.result(branch_key, side).ref(element_id)

    def element(self, ref: ElementRef) -> PartitionInterval:
        """
        The element one reference names.

        Args:
            ref: The element reference.

        Returns:
            Its :class:`~tanglepack.topology.TopologyResults.PartitionInterval`.

        Raises:
            ValueError: If no result covers the ref's ``(branch, side)``.
            IndexError: If the id is not an element of that result.
        """
        return self.result(ref.branch_key, ref.side).element(ref.element_id)

    def interval_at(
        self, branch_key: "ManifoldKey", side: Side, cdist: float
    ) -> PartitionInterval:
        """
        The element owning one stable canonical distance on one side of a branch.

        Args:
            branch_key: Manifold key of the stable branch.
            side: Which side's partition to read.
            cdist: The stable canonical distance to locate.

        Returns:
            The owning :class:`~tanglepack.topology.TopologyResults.PartitionInterval`.

        Raises:
            ValueError: If no result covers ``(branch_key, side)``, if the
                distance lies outside the partitioned stretch, or (a broken
                partition) if more than one element owns it.
        """
        result = self.result(branch_key, side)
        tol = self.tol
        owners = [
            interval
            for interval in result.intervals
            if owns_cdist(interval, cdist, tol)
        ]
        if len(owners) != 1:
            span = (
                f"{result.intervals[0].lo_cdist:.6g}.."
                f"{result.intervals[-1].hi_cdist:.6g}"
                if result.intervals
                else "empty"
            )
            raise ValueError(
                f"{len(owners)} elements of the {side!r} {self.kind} partition of "
                f"branch {branch_key[1:]} own stable cdist {cdist:.6g} (expected "
                f"exactly one); that partition spans {span}"
            )
        return owners[0]

    def element_at(
        self, branch_key: "ManifoldKey", side: Side, cdist: float
    ) -> ElementRef:
        """
        The reference of the element owning one canonical distance.

        Args:
            branch_key: Manifold key of the stable branch.
            side: Which side's partition to read.
            cdist: The stable canonical distance to locate.

        Returns:
            The owning element's :class:`~tanglepack.topology.TopologyResults.ElementRef`.

        Raises:
            ValueError: As :meth:`interval_at`.
        """
        return self.result(branch_key, side).ref(
            self.interval_at(branch_key, side, cdist).element_id
        )

    def owner_of_intersection(self, intersection_id: int, side: Side) -> ElementRef:
        """
        The element owning one registered crossing on one side.

        Args:
            intersection_id: Registry id of a crossing on a covered branch.
            side: Which side's partition to read.

        Returns:
            The owning element's reference, read from the result's
            ``element_of_intersection`` table.

        Raises:
            ValueError: If no covered result on that side knows the crossing.
        """
        for (key, result_side), result in self.results.items():
            if result_side != side:
                continue
            element_id = result.element_of_intersection.get(intersection_id)
            if element_id is not None:
                return result.ref(element_id)
        raise ValueError(
            f"crossing {intersection_id} is on no {side!r}-partitioned branch of "
            f"this {self.kind} partition family"
        )

    # ── identity and reporting ──────────────────────────────────────────────

    def signature(self) -> tuple:
        """
        A structural signature of this family, for cache keying.

        Returns:
            One ``(branch_key, side, ((lo_id, hi_id, closed_lo, closed_hi),
            ...))`` entry per result, in order — the same tuple
            ``TangleSession._partition_signature`` builds for a gathered list.
        """
        return tuple(
            (
                result.branch_key,
                result.side,
                tuple(
                    (iv.lo_id, iv.hi_id, iv.closed_lo, iv.closed_hi)
                    for iv in result.intervals
                ),
            )
            for result in self.results.values()
        )

    def describe(self) -> str:
        """A human-readable report, one line per element in interval notation."""
        lines = [f"{self.kind} partition family ({len(self.results)} result(s))"]
        for result in self.results.values():
            lines.append(
                f"{result.side} partition (p{result.branch_key[0].period}, "
                f"orbit {result.branch_key[2]}.{result.branch_key[3]}):"
            )
            for iv in result.intervals:
                lo = "[" if iv.closed_lo else "("
                hi = "]" if iv.closed_hi else ")"
                lo_name = "anchor" if iv.lo_id is None else f"id {iv.lo_id}"
                hi_name = "end" if iv.hi_id is None else f"id {iv.hi_id}"
                extra = "  (singleton)" if iv.lo_cdist == iv.hi_cdist else ""
                if iv.parent_element_id is not None:
                    extra += f"  (parent #{iv.parent_element_id}"
                    if iv.cut_by is not None:
                        extra += f", cut by {iv.cut_by}"
                    extra += ")"
                lines.append(
                    f"  #{iv.element_id} {lo}{iv.lo_cdist:.4g}, {iv.hi_cdist:.4g}{hi}  "
                    f"({lo_name} -> {hi_name}){extra}"
                )
        return "\n".join(lines)

    def __len__(self) -> int:
        return len(self.results)

    def __iter__(self) -> Iterator[StablePartitionResult]:
        return iter(self.results.values())

    def __contains__(self, key: object) -> bool:
        return key in self.results

    def __repr__(self) -> str:
        elements = sum(len(result.intervals) for result in self.results.values())
        return (
            f"<{type(self).__name__} kind={self.kind!r} results={len(self.results)} "
            f"elements={elements}>"
        )


class HomotopyPartition(PartitionFamily):
    """
    The partition the punched holes define (the Stable Manifold Partition
    Algorithm, unchanged).

    Build one from a trellis that has holes, or wrap results already stored on
    trellises with :meth:`from_results`.
    """

    kind: ClassVar[str] = "homotopy"

    @classmethod
    def build(
        cls,
        trellis: "Trellis",
        branch_keys: Optional[Iterable["ManifoldKey"]] = None,
        sides: tuple[Side, ...] = ("left", "right"),
    ) -> "HomotopyPartition":
        """
        Partition every stable branch of a trellis by its holes.

        Args:
            trellis: A trellis whose holes are punched.
            branch_keys: The stable branches to partition; every stable branch
                by default.
            sides: The sides to partition; both by default.

        Returns:
            The family, one result per ``(branch, side)``, computed by
            :func:`~tanglepack.topology.StablePartition.partition_stable_manifold`.
        """
        keys = (
            [branch.key for branch in trellis.stable_branches]
            if branch_keys is None
            else list(branch_keys)
        )
        results = [
            partition_stable_manifold(trellis, key, side) for key in keys for side in sides
        ]
        return cls(results, trellis=trellis)

    @classmethod
    def from_results(
        cls,
        results: Iterable[StablePartitionResult],
        *,
        trellis: Optional["Trellis"] = None,
    ) -> "HomotopyPartition":
        """
        Wrap already-computed results (a trellis's ``stable_partitions``, or a
        session's gathered ones) without recomputing anything.

        Args:
            results: The results.
            trellis: The trellis to read tolerances and crossings from.

        Returns:
            The family.
        """
        return cls(results, trellis=trellis)


# ── marks: the boundary vocabulary of StablePartition._build_intervals ─────


def _marks_of(result: StablePartitionResult) -> Marks:
    """
    Turn a result back into its marks.

    A non-singleton interval open at an end marks that end's boundary as open
    on that flank; singletons contribute nothing (they are regenerated by the
    rebuild from the two open flanks around them).
    """
    boundaries: dict[Optional[int], float] = {}
    open_outward: set[Optional[int]] = set()
    open_anchorward: set[Optional[int]] = set()
    for iv in result.intervals:
        boundaries.setdefault(iv.lo_id, float(iv.lo_cdist))
        boundaries.setdefault(iv.hi_id, float(iv.hi_cdist))
        if iv.lo_cdist == iv.hi_cdist and iv.lo_id == iv.hi_id:
            continue
        if not iv.closed_lo:
            open_outward.add(iv.lo_id)
        if not iv.closed_hi:
            open_anchorward.add(iv.hi_id)
    return boundaries, open_outward, open_anchorward


def _intervals_from_marks(marks: Marks) -> list[PartitionInterval]:
    """
    Assemble ordered intervals from marks (the tail of
    ``StablePartition._build_intervals``, verbatim in spirit).

    A boundary excluded by both neighbouring intervals owns itself as a closed
    singleton. Element ids are stamped by the caller.
    """
    boundaries, open_outward, open_anchorward = marks
    unique: dict = {}
    for bid, cdist in boundaries.items():
        unique[bid if bid is not None else ("sentinel", cdist)] = (bid, cdist)
    ordered = sorted(unique.values(), key=lambda item: item[1])

    intervals: list[PartitionInterval] = []
    for (lo_id, lo_c), (hi_id, hi_c) in zip(ordered, ordered[1:]):
        intervals.append(
            PartitionInterval(
                lo_id,
                hi_id,
                lo_c,
                hi_c,
                closed_lo=lo_id not in open_outward,
                closed_hi=hi_id not in open_anchorward,
            )
        )

    with_singletons: list[PartitionInterval] = []
    for position, (bid, cdist) in enumerate(ordered):
        before = intervals[position - 1] if position > 0 else None
        after = intervals[position] if position < len(intervals) else None
        owned = (before is not None and before.closed_hi) or (
            after is not None and after.closed_lo
        )
        if not owned and (before is not None or after is not None):
            with_singletons.append(
                PartitionInterval(bid, bid, cdist, cdist, True, True)
            )
        if after is not None:
            with_singletons.append(after)
    return with_singletons


class Cut:
    """
    One audit record of the iterated cut: what one bridge end did to one row.

    Attributes:
        bridge_id: The new bridge.
        intersection_id: The crossing this record is about.
        branch_key: The stable branch of that crossing.
        side: The row the bridge cut (its own).
        opened: ``"anchorward"`` / ``"outward"`` for the neighbouring interval
            that was opened at the crossing (the flank toward the partner),
            or None when nothing changed.
        reason: Why nothing changed (``"existing boundary"``, ``"row
            invariant"``, ``"no partition"``, ``"no partner"``, ``"partner
            disagreement"``, ``"partner on another branch"``, ``"cut
            conflict"``), or None.
        partner: The crossing paired with this one in its image chain (the
            other end of the footprint stretch), or None when unknown.
    """

    def __init__(
        self,
        bridge_id: "BridgeId",
        intersection_id: int,
        branch_key: Optional["ManifoldKey"],
        side: Optional[Side],
        opened: Optional[str],
        reason: Optional[str] = None,
        partner: Optional[int] = None,
    ) -> None:
        """Record one bridge end's effect (see the class attributes)."""
        self.bridge_id = bridge_id
        self.intersection_id = intersection_id
        self.branch_key = branch_key
        self.side = side
        self.opened = opened
        self.reason = reason
        self.partner = partner

    def __repr__(self) -> str:
        effect = f"opened {self.opened}" if self.opened else f"skipped ({self.reason})"
        if self.partner is not None:
            effect += f", toward {self.partner}"
        branch = None if self.branch_key is None else self.branch_key[1:]
        return (
            f"Cut({self.bridge_id} at {self.intersection_id} on {branch}/"
            f"{self.side}: {effect})"
        )


class IteratedHomotopyPartition(PartitionFamily):
    """
    The homotopy partition refined by the image bridges of a minimal trellis.

    Attributes:
        homotopy: The family this one refines.
        minimal: The minimal trellis whose image bridges cut it.
        cuts: One :class:`Cut` per new-bridge end, applied or skipped.
    """

    kind: ClassVar[str] = "iterated_homotopy"

    def __init__(
        self,
        results: Iterable[StablePartitionResult],
        *,
        trellis: Optional["Trellis"] = None,
        homotopy: Optional[HomotopyPartition] = None,
        minimal: Optional["MinimalTrellis"] = None,
        cuts: Optional[list[Cut]] = None,
    ) -> None:
        """
        Args:
            results: The refined results (normally from :meth:`from_minimal`).
            trellis: The full trellis.
            homotopy: The family refined.
            minimal: The minimal trellis that cut it.
            cuts: The audit records.
        """
        super().__init__(results, trellis=trellis)
        self.homotopy = homotopy
        self.minimal = minimal
        self.cuts: list[Cut] = list(cuts or [])

    @classmethod
    def from_minimal(
        cls, minimal: "MinimalTrellis", homotopy: HomotopyPartition
    ) -> "IteratedHomotopyPartition":
        """
        Cut a homotopy partition by every image bridge of a minimal trellis.

        Args:
            minimal: The minimal trellis; its ``image_bridge_ids`` are the
                cutting bridges and its full ``trellis`` resolves them.
            homotopy: The family to refine (every stable branch the image
                bridges touch must be covered on both sides).

        Returns:
            The refined family. With no image bridges it reproduces the
            homotopy family element for element.

        Raises:
            AssertionError: If a crossing of a cut branch ends up owned by
                zero or several elements, or an iterated element does not lie
                inside its parent.
        """
        trellis = minimal.trellis
        marks: dict[tuple["ManifoldKey", Side], Marks] = {
            key: _marks_of(result) for key, result in homotopy.results.items()
        }
        existing: dict[tuple["ManifoldKey", Side], set[Optional[int]]] = {
            key: set(boundaries) for key, (boundaries, _o, _a) in marks.items()
        }
        spans: dict[tuple["ManifoldKey", Side], list[tuple[float, float, "BridgeId"]]] = {}
        cuts: list[Cut] = []
        decided: dict[tuple[tuple["ManifoldKey", Side], int], str] = {}

        def open_at(
            key: tuple["ManifoldKey", Side],
            iid: int,
            which: str,
            bid: "BridgeId",
            partner: int,
        ) -> None:
            """Open the ``which`` neighbour at ``iid`` on ``key``, unless it is
            an existing boundary or another piece already opened the other flank."""
            if key not in marks:
                cuts.append(Cut(bid, iid, key[0], key[1], None, "no partition", partner))
                logger.warning(
                    "image bridge %s ends on branch %s side %s, which the "
                    "homotopy partition does not cover; no cut there",
                    bid,
                    key[0][1:],
                    key[1],
                )
                return
            if iid in existing[key]:
                cuts.append(Cut(bid, iid, key[0], key[1], None, "existing boundary", partner))
                logger.debug(
                    "image bridge %s meets existing boundary %d on %s/%s; "
                    "its closedness stands",
                    bid,
                    iid,
                    key[0][1:],
                    key[1],
                )
                return
            previous = decided.get((key, iid))
            if previous is not None and previous != which:
                cuts.append(Cut(bid, iid, key[0], key[1], None, "cut conflict", partner))
                logger.warning(
                    "image bridge %s would open the %s flank of crossing %d on "
                    "%s/%s, already opened %s by another piece; no cut there",
                    bid,
                    which,
                    iid,
                    key[0][1:],
                    key[1],
                    previous,
                )
                return
            decided[(key, iid)] = which
            boundaries, open_outward, open_anchorward = marks[key]
            boundaries[iid] = float(trellis.intersection(iid).stable_cdist)
            (open_anchorward if which == "anchorward" else open_outward).add(iid)
            cuts.append(Cut(bid, iid, key[0], key[1], which, None, partner))

        partners = cls._chain_partners(minimal)
        for bid in minimal.image_bridge_ids:
            first, second = bid
            branch = {iid: trellis.intersection(iid).manifold_b_key for iid in bid}
            row = {
                first: row_of_end(trellis, bid, "first"),
                second: row_of_end(trellis, bid, "second"),
            }
            same_branch = branch[first] == branch[second]
            if same_branch and row[first] != row[second]:
                logger.warning(
                    "image bridge %s approaches its two crossings from "
                    "different sides of one stable branch (row invariant); "
                    "it cuts nothing",
                    bid,
                )
                cuts.append(Cut(bid, first, branch[first], None, None, "row invariant"))
                cuts.append(Cut(bid, second, branch[second], None, None, "row invariant"))
                continue
            slot = partners.get(frozenset(bid), {})
            for iid in bid:
                key = (branch[iid], row[iid])
                partner = slot.get(iid)
                if partner is None or partner == _DISAGREEMENT:
                    reason = "no partner" if partner is None else "partner disagreement"
                    cuts.append(Cut(bid, iid, key[0], key[1], None, reason))
                    continue
                if trellis.intersection(partner).manifold_b_key != branch[iid]:
                    logger.warning(
                        "image bridge %s: crossing %d is paired with %d on another "
                        "stable branch (a footprint stretch lies on one branch); "
                        "no cut there",
                        bid,
                        iid,
                        partner,
                    )
                    cuts.append(
                        Cut(bid, iid, key[0], key[1], None, "partner on another branch", partner)
                    )
                    continue
                own_c = float(trellis.intersection(iid).stable_cdist)
                partner_c = float(trellis.intersection(partner).stable_cdist)
                # The footprint stretch between the two lies open; open its flank.
                open_at(key, iid, "anchorward" if partner_c < own_c else "outward", bid, partner)
            if same_branch:
                near_c, far_c = sorted(
                    float(trellis.intersection(iid).stable_cdist) for iid in bid
                )
                spans.setdefault((branch[first], row[first]), []).append((near_c, far_c, bid))

        results: list[StablePartitionResult] = []
        for (branch_key, side), result_marks in marks.items():
            branch = trellis.branch(branch_key)
            assert branch is not None, f"{branch_key[1:]} is not a branch of the trellis"
            intervals = _intervals_from_marks(result_marks)
            for element_id, interval in enumerate(intervals):
                interval.element_id = element_id
                interval.branch_key = branch_key
                interval.side = side
                mid = (
                    interval.lo_cdist
                    if interval.lo_cdist == interval.hi_cdist
                    else 0.5 * (interval.lo_cdist + interval.hi_cdist)
                )
                parent = homotopy.interval_at(branch_key, side, mid)
                interval.parent_element_id = parent.element_id
                assert (
                    interval.lo_cdist >= parent.lo_cdist - homotopy.tol
                    and interval.hi_cdist <= parent.hi_cdist + homotopy.tol
                ), (
                    f"iterated element {branch_key[1:]}/{side}#{element_id} "
                    f"[{interval.lo_cdist:.6g}, {interval.hi_cdist:.6g}] is not "
                    f"inside its parent #{parent.element_id} "
                    f"[{parent.lo_cdist:.6g}, {parent.hi_cdist:.6g}]"
                )
                covering = [
                    (far_c - near_c, bid)
                    for near_c, far_c, bid in spans.get((branch_key, side), [])
                    if near_c - homotopy.tol <= interval.lo_cdist
                    and interval.hi_cdist <= far_c + homotopy.tol
                ]
                if covering:
                    interval.cut_by = min(covering, key=lambda item: item[0])[1]
            owners = _element_of_intersection(trellis, branch, intervals)
            results.append(
                StablePartitionResult(
                    branch_key=branch_key,
                    side=side,
                    intervals=intervals,
                    element_of_intersection=owners,
                    elements_at_bridge=_elements_at_bridge(
                        minimal.sparse, branch_key, owners
                    ),
                )
            )

        family = cls(
            results, trellis=trellis, homotopy=homotopy, minimal=minimal, cuts=cuts
        )
        applied = sum(1 for cut in cuts if cut.opened is not None)
        logger.info(
            "iterated homotopy partition: %d image bridge(s) applied %d cut(s) "
            "(%d skipped); %d -> %d elements",
            len(minimal.image_bridge_ids),
            applied,
            len(cuts) - applied,
            sum(len(r.intervals) for r in homotopy.results.values()),
            sum(len(r.intervals) for r in family.results.values()),
        )
        return family

    @staticmethod
    def _chain_partners(
        minimal: "MinimalTrellis",
    ) -> dict[frozenset, dict[int, Optional[int]]]:
        """
        The footprint partner of every crossing of every image chain.

        Args:
            minimal: The minimal trellis whose ``image_chains`` name the hole
                bridges that were mapped forward; each chain is re-read in
                full from the trellis (lookup only).

        Returns:
            ``{frozenset(piece): {crossing id: partner id}}`` over every
            piece of every chain. A crossing maps to None when it is the
            leftover of an even chain and to :data:`_DISAGREEMENT` when two
            chains pair it differently.

        Note:
            Crossings ``x_0 .. x_m`` are paired ``(x_0, x_m), (x_1, x_2),
            (x_3, x_4), ...``: the end stretch and the chords of the image
            lobe. A chain that is not consecutive cuts nothing (WARNING).
        """
        trellis = minimal.trellis
        partners: dict[frozenset, dict[int, Optional[int]]] = {}
        for parent, kept in minimal.image_chains.items():
            chain = _image_chain(trellis, parent) or list(kept)
            if not chain:
                continue
            if any(chain[i][1] != chain[i + 1][0] for i in range(len(chain) - 1)):
                logger.warning(
                    "image chain of %s is not consecutive (%s); its pieces cut nothing",
                    parent,
                    chain,
                )
                for piece in chain:
                    partners.setdefault(frozenset(piece), {})
                continue
            crossings = [chain[0][0]] + [pair[1] for pair in chain]
            pieces = len(chain)
            if pieces % 2 == 0:
                logger.warning(
                    "image chain of %s has an even number of pieces (%d), so its "
                    "ends lie on different rows; pairing (x_0, x_m), (x_1, x_2), "
                    "... regardless and leaving x_%d unpaired",
                    parent,
                    pieces,
                    pieces - 1,
                )
            pairs = [(crossings[0], crossings[-1])] + [
                (crossings[i], crossings[i + 1]) for i in range(1, pieces - 1, 2)
            ]
            partner_of: dict[int, int] = {}
            for x, y in pairs:
                partner_of[x] = y
                partner_of[y] = x
            for piece in chain:
                slot = partners.setdefault(frozenset(piece), {})
                for iid in piece:
                    proposed = partner_of.get(iid)
                    if iid in slot and slot[iid] != proposed:
                        logger.warning(
                            "piece %s: chains disagree on the partner of crossing %d "
                            "(%s vs %s); no cut there",
                            piece,
                            iid,
                            slot[iid],
                            proposed,
                        )
                        slot[iid] = _DISAGREEMENT
                    else:
                        slot.setdefault(iid, proposed)
        return partners

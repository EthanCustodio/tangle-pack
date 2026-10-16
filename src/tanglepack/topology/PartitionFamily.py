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

The iterated partition ``P1 = P ∨ f(∂P)``
-----------------------------------------
``P`` is the homotopy partition. ``P1`` cuts each stable branch, in addition
to ``P``'s cuts, at every crossing of the minimal trellis's image bridges -- the
first images of the hole bridges -- with the stable manifold: the images of
``∂P`` at a chain's ends and the folds of the hole images at its interior
crossings (design note §A.8d). Each NEW bridge of the minimal trellis (an
image pair that is not itself a hole bridge) lies on one side of the stable
manifold — its row, from
:func:`~tanglepack.topology.StablePartition.row_of_end` — and cuts the
partition on THAT side only, by the EMPTY-STRETCH rule (the author's).
The stable manifold has EMPTY stretches and the rest is the inclusive shape:

* a hole's stretch is empty (the open interval of the homotopy partition);
* the image of a hole is empty where it meets the stable manifold: a hole
  ``H`` bounded by one bridge iterates forward into a lobe ``f(H)`` that is
  stretched and folded; it sits on the manifold along its BASE (between the
  images of the hole's two ends, the chain's ``(x_0, x_m)``) and folds through
  the manifold along one CHORD per pair of consecutive interior crossings
  (``(x_1, x_2), (x_3, x_4), ...``, the doorways to its folds on the other
  side). Base and chords are read off the hole bridge's full image chain
  (:func:`~.BridgeClass._image_chain`, iterate table only): on every branch,
  the base branch included, the interior crossings pair in STABLE order.
  The anchor's hole, which
  opens only beyond the manifold end, has no stretch and its image contributes
  nothing. An odd count of interior crossings on one branch means the chain
  left the computed stable manifold there: its outermost crossing stays
  unpaired, no stretch beyond the pairs formed is declared empty, and the
  footprint is recorded in
  :attr:`IteratedHomotopyPartition.incomplete_footprints` (data, INFO).

Empty stretches are empty on BOTH sides of the manifold. An image bridge's
crossing ``x`` belongs, on the row the bridge cuts, to the element AWAY from
the empty stretch it abuts: the flank toward the empty stretch is opened
(``open_anchorward`` when the empty stretch is anchorward of ``x``,
``open_outward`` when it is outward). k=2.8 two blasts: the lobe of the hole
``(6, 7)`` has base ``(8, 11)`` and chord ``(9, 10)`` (its fold is the
reference hole on the left), so the right row reads ``R_1^1 = [0, 11]``,
``R_1^2 = (11, 8)``, ``R_1^3 = [8, 7]`` and ``R_5^1 = [4, 9]``, ``R_5^2 =
(9, 10)``, ``R_5^3 = [10, 3]``; the anchor lobe's fold ``(3, 4)`` on the left
abuts the hole ``(4, 5)`` at ``4``, so ``L_1^1 = [0, 4)``, ``L_1^2 = [4, 9]``.
Every crossing then lies in the element of the stable edge the dual-graph
walk crosses at it (a walk never enters a hole face), so a refined class's
occurrence names exactly the elements its member bridges land in. A boundary
the homotopy partition already has is never re-cut; the other side is not cut
at that crossing; a subdivided image cuts each row once per piece.

Dev Notes:

* Marks. The homotopy result is turned back into the marks vocabulary of
  ``StablePartition`` — boundary cdists plus the sets of boundaries whose
  OUTWARD / ANCHORWARD interval is open — the new bridges add their marks,
  and the intervals are rebuilt by the homotopy partition's own
  ``StablePartition.intervals_from_marks`` (the branch's ends are already
  among the marks, so nothing is re-derived from the trellis);
  ``test_partition_family.py`` pins that rebuilding every homotopy result
  from its own marks reproduces it exactly. A cut is one
  mark per crossing and flank; a crossing that abuts empty stretches on BOTH
  flanks (a hole against a base) is pinched into a closed singleton, as the
  homotopy partition does between two holes (INFO). Two pieces sharing a
  crossing on one row ask for the same flanks (they depend on the crossing
  alone); the second is recorded as ``"already cut"``.
* Parent by exact containment. An iterated element's homotopy parent is
  the homotopy span containing it, compared exactly
  (:func:`_parent_interval`; a singleton by its crossing's id), never by a
  midpoint owner under a slack: a real sliver of width 1.1e-6 beyond the
  zone's pip (k=2.8, 9 blasts) lies within any such slack of the closed pip
  singleton. The ``cut_by`` covering test is exact
  for the same reason; :meth:`PartitionFamily.interval_at` is exact.
* The rule is the author's: "all of the stretched and folded shape is
  inclusive, the empty parts between the bends are open", read only off
  holes with a stretch.
* A crossing that abuts no empty stretch gets no cut. When it lies on the
  image chain of a hole bridge whose hole has no stretch of its own (the
  anchor's one-sided hole, or a lobe whose stable span holds other holes'
  crossings), that is the author's rule, by design: only holes with
  a stretch mark a base and chords (INFO, ``"unmarked lobe"``; every occurrence on
  the law cases and the deep builders is this kind:
  inversion_two_blasts, 13 on k28_deep_blast_sliver, 13 / 47 per recompute on
  k28 / k32 at depth). Anything else is a WARNING, ``"no empty stretch"``.
* Chords. An image chain can leave its base branch
  for an excursion across another stable branch (another orbit branch or fixed
  point: p3_6_blasts, p3_15, p3_16, heteroclinic), or spiral back across its
  own (``k32_series`` from 2 blasts). A lobe ``H`` (its bridge and the stable
  segment between the bridge's ends) holds no point of any stable branch of
  the trellis, so ``f(H)`` holds no anchor; walking ANY branch from its anchor
  therefore enters and leaves ``f(H)`` at successive crossings with the image
  chain, the base counting as inside, and the stretches inside are its
  crossings paired in STABLE order (checked by point-in-polygon on both
  branch kinds: ``tests/regression/test_foreign_chords.py``, and per blast
  step in ``tests/regression/test_class_count_under_blasting.py``). Every
  empty stretch lies on one branch, and the base is read at a non-anchor end
  (an anchor lies on every branch); both are asserted.
* Cross-branch pieces (a period-k orbit only; a bridge whose two crossings
  lie on different stable branches) are cut like any other: the rule reads
  the crossing's own branch.
* Only pieces with a ``Bridge`` object cut (``minimal.image_bridge_ids``); a
  chain pair the minimal trellis skipped (no bridge object) cuts nothing,
  as ``MinimalTrellis`` already warns, but it still counts as a base or chord.
* :class:`~tanglepack.topology.TopologyResults.ElementRef` stays
  ``(branch_key, side, element_id)`` with no family tag. The two families of
  one trellis therefore produce refs that collide by value; never mix refs
  from different families in one structure. Iterated elements point back at
  their homotopy parent through ``PartitionInterval.parent_element_id``.
"""

from __future__ import annotations

import logging
from collections import Counter
from typing import ClassVar, Iterable, Iterator, Optional, TYPE_CHECKING

from .BridgeClass import _image_chain, _index_partitions
from .StablePartition import (
    Marks,
    _element_of_intersection,
    intervals_from_marks,
    owns_cdist,
    partition_stable_manifold,
    row_of_end,
)
from .TopologyResults import (
    ElementRef,
    Endpoint,
    PartitionInterval,
    Side,
    StablePartitionResult,
    endpoint_index,
)

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ..numerics.Bridge import BridgeId
    from ..numerics.Intersection import ManifoldKey
    from .MinimalTrellis import MinimalTrellis
    from .Trellis import Trellis

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class PartitionFamily:
    """
    One partition of the stable manifold, one result per ``(branch, side)``.

    Attributes:
        results: The results indexed by ``(branch_key, side)``, in insertion
            order.
        trellis: The trellis the results were computed over (supplies
            crossing lookups). Optional for
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
        self,
        branch_key: "ManifoldKey",
        side: Side,
        cdist: float,
    ) -> PartitionInterval:
        """
        The element owning one stable canonical distance on one side of a branch, exactly.

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

        Note:
            Partition boundaries ARE registry floats, so a registry cdist
            resolves exactly, however close its neighbours sit. An ESTIMATED
            cdist (a midpoint of ``land_element``'s image span) is located
            exactly too: no absolute slack is read.
        """
        result = self.result(branch_key, side)
        owners = [
            interval for interval in result.intervals if owns_cdist(interval, cdist, 0.0)
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

    def owner_at_end(
        self, trellis: "Trellis", bridge_id: "BridgeId", endpoint: Endpoint
    ) -> ElementRef:
        """
        The element one end of an arc sits against: on its row, on its stable branch.

        Args:
            trellis: The trellis the arc's crossings are registered in.
            bridge_id: The arc's two crossings, ordered like a
                :data:`~tanglepack.numerics.Bridge.BridgeId`.
            endpoint: ``"first"`` or ``"second"``.

        Returns:
            The owner on :func:`row_of_end`'s side of the branch
            :meth:`~.Trellis.Trellis.stable_key_at` names (at an anchor, the
            branch of the arc's own orbit point and branch index).

        Raises:
            ValueError: As :meth:`owner_of_intersection`.
        """
        intersection_id = bridge_id[endpoint_index(endpoint)]
        branch_key = trellis.stable_key_at(intersection_id, trellis.unstable_key_of(bridge_id))
        return self.owner_of_intersection(
            intersection_id, row_of_end(trellis, bridge_id, endpoint), branch_key
        )

    def owner_of_intersection(
        self,
        intersection_id: int,
        side: Side,
        branch_key: "ManifoldKey",
    ) -> ElementRef:
        """
        The element owning one registered crossing on one side of one branch.

        Args:
            intersection_id: Registry id of a crossing on ``branch_key``.
            side: Which side's partition to read.
            branch_key: The stable branch to read it on (always named: an
                anchor lies on every branch of its point).

        Returns:
            The owning element's reference, read from the result's
            ``element_of_intersection`` table.

        Raises:
            ValueError: If that branch and side is not covered, or its result
                does not know the crossing.
        """
        result = self.results.get((branch_key, side))
        element_id = None if result is None else result.element_of_intersection.get(intersection_id)
        if element_id is None:
            raise ValueError(
                f"crossing {intersection_id} is on no {side!r}-partitioned branch "
                f"{branch_key[1:]} of this {self.kind} partition family"
            )
        return result.ref(element_id)

    # ── identity and reporting ──────────────────────────────────────────────

    def signature(self) -> tuple:
        """
        A structural signature of this family, for comparing two builds.

        Returns:
            One ``(branch_key, side, ((lo_id, hi_id, closed_lo, closed_hi),
            ...))`` entry per result, in order.
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


# ── marks: the boundary vocabulary of StablePartition.intervals_from_marks ──


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


def _parent_interval(
    child: PartitionInterval, parent_result: StablePartitionResult
) -> PartitionInterval:
    """
    The homotopy element an iterated element refines, by exact containment.

    A homotopy boundary is never re-cut, so the iterated boundaries are a
    superset of the homotopy ones and a span lies inside exactly one homotopy
    span, compared exactly (all boundaries are registry floats). A singleton
    ``[x, x]`` — on a homotopy boundary, or pinched at a new one — belongs to
    the element owning the crossing ``x`` (by id); the synthetic anchor's to
    the element owning cdist 0.0.
    """
    if child.lo_id == child.hi_id:
        element_id = parent_result.element_of_intersection.get(child.lo_id)
        if element_id is not None:
            return parent_result.element(element_id)
        candidates = [
            iv for iv in parent_result.intervals if owns_cdist(iv, child.lo_cdist, 0.0)
        ]
    else:
        candidates = [
            iv
            for iv in parent_result.intervals
            if iv.lo_id != iv.hi_id
            and iv.lo_cdist <= child.lo_cdist
            and child.hi_cdist <= iv.hi_cdist
        ]
    assert len(candidates) == 1, (
        f"iterated element {child.branch_key[1:]}/{child.side}#{child.element_id} "
        f"(ids {child.lo_id}, {child.hi_id}) [{child.lo_cdist!r}, {child.hi_cdist!r}] "
        f"lies inside {len(candidates)} homotopy elements "
        f"{[(iv.element_id, iv.lo_id, iv.hi_id) for iv in candidates]}, expected one"
    )
    return candidates[0]


class Cut:
    """
    One audit record of the iterated cut: what one bridge end did to one row.

    Attributes:
        bridge_id: The new bridge.
        intersection_id: The crossing this record is about.
        branch_key: The stable branch of that crossing.
        side: The row the bridge cut (its own).
        opened: ``"anchorward"`` / ``"outward"`` for the neighbouring interval
            that was opened at the crossing (the flank toward the empty
            stretch, so the crossing stays with the element away from it),
            or None when nothing changed. A crossing pinched between two
            empty stretches has one record per flank.
        reason: Why nothing changed (``"existing boundary"``, ``"row
            invariant"``, ``"no partition"``, ``"unmarked lobe"`` (the crossing
            lies on the image lobe of a hole with no stretch, which marks
            nothing by design), ``"no empty stretch"``, ``"already cut"``), or
            None.
        partner: The far end of the empty stretch this crossing abuts (the
            other end of the hole, base or chord), or None when unknown.
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
        incomplete_footprints: ``(hole bridge, stable branch, count)`` for every
            branch where a hole's image chain crosses an ODD number of times:
            the chain left the computed stable manifold there, so its last
            crossing is unpaired (data, never a reliability failure).
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
        incomplete_footprints: Iterable[tuple["BridgeId", "ManifoldKey", int]] = (),
    ) -> None:
        """
        Args:
            results: The refined results (normally from :meth:`from_minimal`).
            trellis: The full trellis.
            homotopy: The family refined.
            minimal: The minimal trellis that cut it.
            cuts: The audit records.
            incomplete_footprints: The hole images with an odd crossing count
                on a branch, ``(hole bridge, branch, count)``.
        """
        super().__init__(results, trellis=trellis)
        self.homotopy = homotopy
        self.minimal = minimal
        self.cuts: list[Cut] = list(cuts or [])
        self.incomplete_footprints = list(incomplete_footprints)

    def describe(self) -> str:
        """The family's report, then one line per incomplete footprint."""
        lines = [super().describe()]
        lines.extend(
            f"incomplete footprint: hole bridge {bridge} crosses stable branch "
            f"{branch[1:]} {count} time(s), the last crossing unpaired"
            for bridge, branch, count in self.incomplete_footprints
        )
        return "\n".join(lines)

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
        decided: set[tuple[tuple["ManifoldKey", Side], int, str]] = set()

        def open_at(
            key: tuple["ManifoldKey", Side],
            iid: int,
            which: str,
            bid: "BridgeId",
            partner: int,
        ) -> None:
            """Open the ``which`` neighbour at ``iid`` on ``key``, unless it is
            an existing boundary or another piece already opened that flank."""
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
            if (key, iid, which) in decided:
                cuts.append(Cut(bid, iid, key[0], key[1], None, "already cut", partner))
                logger.debug(
                    "image bridge %s: the %s flank of crossing %d on %s/%s was already "
                    "opened by another piece",
                    bid,
                    which,
                    iid,
                    key[0][1:],
                    key[1],
                )
                return
            decided.add((key, iid, which))
            boundaries, open_outward, open_anchorward = marks[key]
            boundaries[iid] = float(trellis.intersection(iid).stable_cdist)
            (open_anchorward if which == "anchorward" else open_outward).add(iid)
            cuts.append(Cut(bid, iid, key[0], key[1], which, None, partner))

        empty = cls._empty_stretches(minimal, homotopy)
        unmarked: Optional[dict[int, "BridgeId"]] = None
        for bid in minimal.image_bridge_ids:
            first, second = bid
            unstable_key = trellis.unstable_key_of(bid)
            branch = {iid: trellis.stable_key_at(iid, unstable_key) for iid in bid}
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
            for iid in bid:
                key = (branch[iid], row[iid])
                if key in marks and iid in existing[key]:
                    cuts.append(Cut(bid, iid, key[0], key[1], None, "existing boundary"))
                    logger.debug(
                        "image bridge %s meets existing boundary %d on %s/%s; "
                        "its closedness stands",
                        bid,
                        iid,
                        key[0][1:],
                        key[1],
                    )
                    continue
                abutting = empty.get(iid, [])
                if not abutting:
                    if unmarked is None:
                        unmarked = cls._unmarked_lobes(minimal, homotopy)
                    if iid in unmarked:
                        cuts.append(Cut(bid, iid, key[0], key[1], None, "unmarked lobe"))
                        logger.info(
                            "image bridge %s: crossing %d lies on the image lobe of hole "
                            "bridge %s, whose hole has no stretch of its own (the "
                            "anchor's one-sided hole, or a lobe holding other holes): by "
                            "the empty-stretch rule (author, 2026-09-23) it marks "
                            "nothing, so its element is not cut there",
                            bid,
                            iid,
                            unmarked[iid],
                        )
                    else:
                        cuts.append(Cut(bid, iid, key[0], key[1], None, "no empty stretch"))
                        logger.warning(
                            "image bridge %s: crossing %d abuts no hole and no image of "
                            "a hole, so its element is not cut there",
                            bid,
                            iid,
                        )
                    continue
                own_c = float(trellis.intersection(iid).stable_cdist)
                flanks: dict[str, int] = {}
                for other, kind in abutting:
                    if other not in trellis.branch(key[0]).intersection_ids:
                        # Every empty stretch lies on one branch; only an anchor,
                        # which lies on all of them, abuts one along another branch.
                        assert trellis.intersection(iid).is_anchor, (bid, iid, other, kind)
                        continue
                    other_c = float(trellis.intersection(other).stable_cdist)
                    flanks.setdefault("anchorward" if other_c < own_c else "outward", other)
                if len(flanks) == 2:
                    logger.info(
                        "image bridge %s: crossing %d is pinched between empty "
                        "stretches on both flanks (%s); it becomes a singleton",
                        bid,
                        iid,
                        abutting,
                    )
                # The crossing stays with the element AWAY from the empty stretch.
                for which, other in flanks.items():
                    open_at(key, iid, which, bid, other)
            if same_branch:
                near_c, far_c = sorted(
                    float(trellis.intersection(iid).stable_cdist) for iid in bid
                )
                spans.setdefault((branch[first], row[first]), []).append((near_c, far_c, bid))

        results: list[StablePartitionResult] = []
        for (branch_key, side), result_marks in marks.items():
            branch = trellis.branch(branch_key)
            assert branch is not None, f"{branch_key[1:]} is not a branch of the trellis"
            intervals = intervals_from_marks(result_marks)
            for element_id, interval in enumerate(intervals):
                interval.element_id = element_id
                interval.branch_key = branch_key
                interval.side = side
                parent = _parent_interval(interval, homotopy.result(branch_key, side))
                interval.parent_element_id = parent.element_id
                assert (
                    parent.lo_cdist <= interval.lo_cdist
                    and interval.hi_cdist <= parent.hi_cdist
                ), (
                    f"iterated element {branch_key[1:]}/{side}#{element_id} "
                    f"(ids {interval.lo_id}, {interval.hi_id}) "
                    f"[{interval.lo_cdist!r}, {interval.hi_cdist!r}] of width "
                    f"{interval.hi_cdist - interval.lo_cdist!r} is not inside its "
                    f"parent #{parent.element_id} (ids {parent.lo_id}, {parent.hi_id}) "
                    f"[{parent.lo_cdist!r}, {parent.hi_cdist!r}]"
                )
                covering = [
                    (far_c - near_c, bid)
                    for near_c, far_c, bid in spans.get((branch_key, side), [])
                    if near_c <= interval.lo_cdist and interval.hi_cdist <= far_c
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
                )
            )

        family = cls(
            results,
            trellis=trellis,
            homotopy=homotopy,
            minimal=minimal,
            cuts=cuts,
            incomplete_footprints=cls._incomplete_footprints(minimal, homotopy),
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
    def _empty_stretches(
        minimal: "MinimalTrellis", homotopy: HomotopyPartition
    ) -> dict[int, list[tuple[int, str]]]:
        """
        The empty stretches of the stable manifold, indexed by their ends.

        Args:
            minimal: The minimal trellis (its hole bridges and full trellis).
            homotopy: The homotopy partition, whose open intervals are the
                holes' stretches.

        Returns:
            ``{crossing id: [(other end, kind), ...]}`` for every end of an
            empty stretch, ``kind`` being ``"hole"`` (a hole's own stretch),
            ``"base"`` (where a hole's image lobe sits on the manifold) or
            ``"chord"`` (where it folds through to the other side). Stretches
            are recorded without a side: they are empty on both.

        Note:
            A hole's lobe is its bounding bridge's full image chain (lookup
            only), read for every hole bridge whose stretch is an open
            interval of the homotopy partition; the anchor's one-sided hole
            has none. Its ``(x_0, x_m)`` is the base, on one stable branch.
            Interior crossings pair as chords per stable branch, in stable
            order on every branch, the base branch included (a lobe
            image holds no anchor, so walking a branch from its anchor enters
            and leaves it in turn). An odd count on a branch leaves its
            outermost crossing unpaired: the footprint is incomplete
            (:meth:`_incomplete_footprints`).
        """
        trellis = minimal.trellis
        stretches = dict.fromkeys(_hole_stretches(homotopy), "hole")
        for hole_bridge in minimal.hole_bridge_ids:
            if frozenset(hole_bridge) not in stretches:
                continue
            chain = _image_chain(trellis, hole_bridge)
            if not chain:
                continue
            if any(chain[i][1] != chain[i + 1][0] for i in range(len(chain) - 1)):
                logger.warning(
                    "image chain of hole bridge %s is not consecutive (%s); its lobe "
                    "marks no empty stretch",
                    hole_bridge,
                    chain,
                )
                continue
            crossings = [chain[0][0]] + [pair[1] for pair in chain]
            # The base f(sigma) lies on ONE stable branch: its ends are the +1
            # images of one branch's crossings. An anchor lies on every branch of
            # its point, so the branch is read at a non-anchor end.
            ends = {
                trellis.intersection(c).manifold_b_key
                for c in (crossings[0], crossings[-1])
                if not trellis.intersection(c).is_anchor
            }
            assert len(ends) == 1, (hole_bridge, chain)
            stretches.setdefault(frozenset((crossings[0], crossings[-1])), "base")
            # Chords: f(H) holds no anchor, so walking any branch from its anchor
            # enters and leaves it in turn, the base branch no exception (Jordan):
            # crossings pair in STABLE order.
            folds: dict = {}
            for c in crossings[1:-1]:
                folds.setdefault(trellis.intersection(c).manifold_b_key, []).append(c)
            for on_branch in folds.values():
                on_branch.sort(key=lambda c: trellis.intersection(c).stable_cdist)
                for x, y in zip(on_branch[::2], on_branch[1::2]):
                    stretches.setdefault(frozenset((x, y)), "chord")
        empty: dict[int, list[tuple[int, str]]] = {}
        for stretch, kind in stretches.items():
            if len(stretch) != 2:
                continue
            a, b = sorted(stretch)
            empty.setdefault(a, []).append((b, kind))
            empty.setdefault(b, []).append((a, kind))
        return empty

    @staticmethod
    def _incomplete_footprints(
        minimal: "MinimalTrellis", homotopy: HomotopyPartition
    ) -> list[tuple["BridgeId", "ManifoldKey", int]]:
        """
        The lobe images whose footprint on some stable branch is incomplete.

        Args:
            minimal: The minimal trellis (its hole bridges and full trellis).
            homotopy: The homotopy partition, whose open intervals are the
                holes' stretches.

        Returns:
            ``(hole bridge, stable branch, count)`` for every branch the image
            chain of a hole with a stretch crosses an ODD number of times at
            its interior crossings: the chain left the computed stable
            manifold there (a "stable end"), so :meth:`_empty_stretches`
            leaves the outermost crossing unpaired. Logged once at INFO when there
            is any: data, never a reliability failure.
        """
        trellis = minimal.trellis
        stretched = _hole_stretches(homotopy)
        incomplete = []
        for hole_bridge in minimal.hole_bridge_ids:
            if frozenset(hole_bridge) not in stretched:
                continue
            chain = _image_chain(trellis, hole_bridge) or []
            if any(pair[1] != after[0] for pair, after in zip(chain, chain[1:])):
                continue  # marks no empty stretch (``_empty_stretches`` warns)
            folds = Counter(
                trellis.intersection(second).manifold_b_key for _first, second in chain[:-1]
            )
            incomplete.extend(
                (hole_bridge, branch_key, count)
                for branch_key, count in folds.items()
                if count % 2
            )
        if incomplete:
            logger.info(
                "an odd number of fold crossings, last one unpaired (incomplete "
                "footprint): %s",
                [(bridge, branch[1:], count) for bridge, branch, count in incomplete],
            )
        return incomplete

    @staticmethod
    def _unmarked_lobes(
        minimal: "MinimalTrellis", homotopy: HomotopyPartition
    ) -> dict[int, "BridgeId"]:
        """Every crossing on the image chain of a hole bridge whose hole has no
        stretch (the anchor's one-sided hole, or a lobe whose stable span holds
        other holes' crossings), mapped to that hole bridge: by the author's
        rule such a lobe marks no empty stretch."""
        trellis = minimal.trellis
        stretched = _hole_stretches(homotopy)
        return {
            c: hole_bridge
            for hole_bridge in minimal.hole_bridge_ids
            if frozenset(hole_bridge) not in stretched
            for pair in _image_chain(trellis, hole_bridge) or ()
            for c in pair
        }


def _hole_stretches(homotopy: HomotopyPartition) -> set[frozenset]:
    """The holes' own stretches: the open intervals of the homotopy partition,
    each as the frozenset of its two bounding crossing ids."""
    return {
        frozenset((iv.lo_id, iv.hi_id))
        for result in homotopy.results.values()
        for iv in result.intervals
        if iv.lo_id is not None
        and iv.hi_id is not None
        and iv.lo_id != iv.hi_id
        and not iv.closed_lo
        and not iv.closed_hi
    }

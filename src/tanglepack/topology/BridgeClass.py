"""
Bridge classes: the pair of partition elements a bridge connects.

A bridge runs from one crossing to the next along an unstable branch, and each
of its two ends sits against a stable branch on a definite side (its *row*, see
:func:`~.StablePartition.row_of_end`). The stable partition cuts that branch
into elements, so each end names one element; the ordered pair of those two
:class:`~.TopologyResults.ElementRef` s is the bridge's CLASS. Two bridges of
the same class connect the same two elements of the stable partition, which
is the coarsest identity the dual graph can tell apart.

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
* Partitions are passed in rather than read off the trellis because they do not
  live where the bridges do: a trellis holds a bridge under its UNSTABLE fixed
  point, while the partition owning that bridge's endpoints lives on the trellis
  of the STABLE branch's fixed point. Running this on the all-fixed-points
  trellis with every per-fixed-point trellis's partitions is what covers the
  nested and heteroclinic cases in one pass.
* Ordering is deliberate and total (:func:`class_sort_key`): the returned dict
  is built in sorted order and each class's members are sorted by
  :data:`~tanglepack.numerics.Bridge.BridgeId`, so two runs of the same trellis
  enumerate the classes identically.

Open question: a class is currently the finest invariant we have — two bridges
sharing a class are assumed to behave alike under the map. Nothing checks
that yet; a consumer that needs it must verify it.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Iterable, Optional, Sequence, TYPE_CHECKING

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
    The ordered pair of partition elements a bridge connects.

    Ordered like a :data:`~tanglepack.numerics.Bridge.BridgeId`, i.e. by the
    unstable dynamical direction: :attr:`first` is the element at the end with
    the smaller unstable canonical distance. Frozen, so it hashes and is usable
    as a symbol.

    Attributes:
        first: Element at the bridge's ``"first"`` end, on that end's row.
        second: Element at the bridge's ``"second"`` end, on that end's row.
    """

    first: ElementRef
    second: ElementRef

    @property
    def label(self) -> str:
        """A short, deterministic name: the two elements' labels, arrow-joined."""
        return f"{self.first.label} -> {self.second.label}"

    def __str__(self) -> str:
        return self.label


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
        The concatenation of :func:`element_sort_key` for ``first`` and then for
        ``second`` — a 10-tuple.
    """
    return element_sort_key(bridge_class.first, fixed_points) + element_sort_key(
        bridge_class.second, fixed_points
    )


def bridge_classes(
    trellis: "Trellis",
    partitions: Iterable[StablePartitionResult],
    *,
    bridge_ids: Optional[Iterable["BridgeId"]] = None,
) -> dict[BridgeClass, list["BridgeId"]]:
    """
    Group a trellis's bridges by the pair of partition elements they connect.

    For each bridge id, the row at each end (:func:`~.StablePartition.row_of_end`)
    picks the side, the endpoint's ``manifold_b_key`` picks the stable branch,
    and the partition result for that ``(branch, side)`` names the element that
    owns the endpoint. The two elements are the bridge's
    :class:`BridgeClass`.

    Args:
        trellis: The trellis whose crossings and bridges to read. Normally the
            all-fixed-points trellis, so that nested and heteroclinic bridges
            are covered in one pass.
        partitions: The stable partitions to resolve elements against — every
            per-fixed-point trellis's :attr:`~.Trellis.Trellis.stable_partitions`
            concatenated. Indexed by ``(branch_key, side)``.
        bridge_ids: The bridges to class. Defaults to every non-partial bridge
            of ``trellis``.

    Returns:
        A dict from class to the sorted list of that class's
        :data:`~tanglepack.numerics.Bridge.BridgeId` s, itself ordered by
        :func:`class_sort_key`.

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

    grouped: dict[BridgeClass, list["BridgeId"]] = {}
    for bridge_id in ids:
        first = _element_at_end(trellis, by_branch_side, bridge_id, "first")
        second = _element_at_end(trellis, by_branch_side, bridge_id, "second")
        grouped.setdefault(BridgeClass(first, second), []).append(bridge_id)

    fixed_points = trellis.fixed_points
    logger.debug("classed %d bridges into %d classes", len(ids), len(grouped))
    return {
        bridge_class: sorted(members)
        for bridge_class, members in sorted(
            grouped.items(), key=lambda item: class_sort_key(item[0], fixed_points)
        )
    }


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

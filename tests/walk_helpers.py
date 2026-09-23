"""
Synthetic builders for dual-walk tests: fake fixed points, stable branches,
partition results, a duck-typed trellis, and real ``StableNode`` /
``FaceNode`` objects wired into a tiny ``FakeDual``.

Nothing here touches the map or the registry; every builder produces the
same dataclasses the real :class:`~tanglepack.topology.DualGraph.DualGraph`
produces, with ``arc=None`` where the geometry is irrelevant.
"""

from __future__ import annotations

from types import SimpleNamespace
from typing import Iterable, Optional

from tanglepack.topology.DualGraph import FaceNode, StableNode
from tanglepack.topology.TopologyResults import (
    ElementRef,
    PartitionInterval,
    Side,
    StablePartitionResult,
)


# --------------------------------------------------------------------------- #
# Fixed points and branch keys
# --------------------------------------------------------------------------- #
class FakeFixedPoint:
    """
    A stand-in for :class:`~tanglepack.numerics.FixedPoint.FixedPoint`.

    Attributes:
        period: The orbit period (used by ``ElementRef.label``).
        k_value: Steps of a branch return.
        beta: The per-step STABLE canonical-distance factor (< 1); the
            unstable factor is its reciprocal.
    """

    def __init__(self, period: int = 1, k_value: int = 1, beta: float = 0.25) -> None:
        self.period = period
        self.k_value = k_value
        self.beta = beta
        self.name = f"fp(p{period})"

    @property
    def num_branches(self) -> int:
        return 2 if self.k_value == 2 * self.period else 1

    def per_step_beta(self, stability: str) -> float:
        return self.beta if stability == "stable" else 1.0 / self.beta

    def advance_key(self, key: tuple, n: int = 1) -> tuple:
        """One orbit index step per map step; branch flips on a wrap with inversion."""
        fixed_point, stability, orbit, branch = key
        total = orbit + n
        wraps = total // self.period
        new_orbit = total % self.period
        if self.num_branches == 2 and wraps % 2:
            branch = 1 - branch
        return (fixed_point, stability, new_orbit, branch)

    def __repr__(self) -> str:
        return self.name


def stable_key(fixed_point: FakeFixedPoint, orbit: int = 0, branch: int = 0) -> tuple:
    """A stable manifold key on a fake fixed point."""
    return (fixed_point, "stable", orbit, branch)


def ref(branch_key: tuple, side: Side, element_id: int) -> ElementRef:
    """Shorthand for an :class:`ElementRef`."""
    return ElementRef(branch_key, side, element_id)


# --------------------------------------------------------------------------- #
# Partition results by hand
# --------------------------------------------------------------------------- #
#: One interval spec: ``(lo_id, hi_id, lo_cdist, hi_cdist, closed_lo, closed_hi)``.
IntervalSpec = tuple[Optional[int], Optional[int], float, float, bool, bool]


def make_result(
    branch_key: tuple,
    side: Side,
    intervals: Iterable[IntervalSpec],
    owners: Optional[dict[int, int]] = None,
    parents: Optional[Iterable[Optional[int]]] = None,
) -> StablePartitionResult:
    """
    A stamped :class:`StablePartitionResult` from interval specs.

    Args:
        branch_key: The stable branch.
        side: The side.
        intervals: Anchor-outward interval specs.
        owners: ``element_of_intersection`` (crossing id -> element id). When
            omitted it is derived: every registered end of every interval is
            owned by the interval closed there (a singleton owns itself).
        parents: Optional ``parent_element_id`` per interval.

    Returns:
        The result with element ids, branch key and side stamped.
    """
    built: list[PartitionInterval] = []
    parent_list = list(parents) if parents is not None else None
    for index, (lo_id, hi_id, lo_c, hi_c, closed_lo, closed_hi) in enumerate(intervals):
        interval = PartitionInterval(lo_id, hi_id, float(lo_c), float(hi_c), closed_lo, closed_hi)
        interval.element_id = index
        interval.branch_key = branch_key
        interval.side = side
        if parent_list is not None:
            interval.parent_element_id = parent_list[index]
        built.append(interval)
    if owners is None:
        owners = {}
        for interval in built:
            if interval.lo_id is not None and interval.closed_lo:
                owners.setdefault(interval.lo_id, interval.element_id)
            if interval.hi_id is not None and interval.closed_hi:
                owners.setdefault(interval.hi_id, interval.element_id)
    return StablePartitionResult(
        branch_key=branch_key,
        side=side,
        intervals=built,
        element_of_intersection=dict(owners),
    )


# --------------------------------------------------------------------------- #
# A duck-typed trellis
# --------------------------------------------------------------------------- #
class FakeTrellis:
    """
    The slice of :class:`~tanglepack.topology.Trellis.Trellis` the walk layer reads.

    Attributes:
        cdists: Stable canonical distance per crossing id.
        iterates: ``+1`` iterate per crossing id (the inverse is derived).
        keys: Stable branch key per crossing id (one key for all when a single
            key is given).
        signs: ``crossing_sign`` per crossing id (default +1).
        orientation_preserving: The map's orientation.
        fixed_points: The ordering context.
        registry: Carries ``cdist_tol``.
    """

    def __init__(
        self,
        cdists: dict[int, float],
        iterates: Optional[dict[int, int]] = None,
        keys: "tuple | dict[int, tuple]" = (),
        *,
        signs: Optional[dict[int, int]] = None,
        orientation_preserving: bool = True,
        cdist_tol: float = 1e-6,
        fixed_points: Iterable = (),
    ) -> None:
        self.cdists = dict(cdists)
        self.iterates = dict(iterates or {})
        self.preimages = {image: source for source, image in self.iterates.items()}
        self.keys = keys
        self.signs = dict(signs or {})
        self.orientation_preserving = orientation_preserving
        self.registry = SimpleNamespace(cdist_tol=cdist_tol)
        self.fixed_points = list(fixed_points)

    def key_of(self, crossing_id: int) -> tuple:
        if isinstance(self.keys, dict):
            return self.keys[crossing_id]
        return self.keys

    def intersection(self, crossing_id: int) -> SimpleNamespace:
        return SimpleNamespace(
            id=crossing_id,
            stable_cdist=self.cdists[crossing_id],
            crossing_sign=self.signs.get(crossing_id, 1),
            manifold_b_key=self.key_of(crossing_id),
        )

    def iterate(self, crossing_id: int, n: int) -> Optional[int]:
        current: Optional[int] = crossing_id
        table = self.iterates if n > 0 else self.preimages
        for _ in range(abs(n)):
            if current is None:
                return None
            current = table.get(current)
        return current

    def image_cdist(self, crossing_id: int, n: int = 1, stability: str = "stable"):
        assert stability == "stable", "the fake only follows the stable side"
        key = self.key_of(crossing_id)
        fixed_point = key[0]
        image = self.iterate(crossing_id, n)
        if image is not None:
            return self.key_of(image), self.cdists[image], True
        return (
            fixed_point.advance_key(key, n),
            self.cdists[crossing_id] * fixed_point.per_step_beta("stable") ** n,
            False,
        )


# --------------------------------------------------------------------------- #
# Synthetic dual graphs
# --------------------------------------------------------------------------- #
class FakeDual:
    """
    The one method :func:`~tanglepack.topology.DualWalk.shortest_walks` needs.

    Attributes:
        faces: The face nodes, ``faces[i].index == i``.
        nodes: Every stable node, in creation order.
    """

    def __init__(self) -> None:
        self.faces: list[FaceNode] = []
        self.nodes: list[StableNode] = []
        self._next_edge = 0

    def face(self, kind: str = "region", *, unbounded: bool = False) -> FaceNode:
        """Add a face node."""
        node = FaceNode(index=len(self.faces), faces=[], kind=kind, is_unbounded=unbounded)
        self.faces.append(node)
        return node

    def _edge_key(self, branch_key: tuple) -> tuple:
        lo = 2 * self._next_edge
        self._next_edge += 1
        return ("stable", lo, lo + 1, branch_key)

    def unified(
        self,
        left_face: FaceNode,
        right_face: FaceNode,
        left_element: ElementRef,
        right_element: ElementRef,
        *,
        lo_cdist: float = 0.0,
        hi_cdist: float = 1.0,
    ) -> StableNode:
        """
        Add one solid (traversable) node between two faces.

        The same face may be passed on both sides (merged-face self-adjacency).
        """
        branch_key = left_element.branch_key
        edge_key = self._edge_key(branch_key)
        node = StableNode(
            arc=None,
            edge_key=edge_key,
            key=(edge_key, "both"),
            branch_key=branch_key,
            lo_cdist=lo_cdist,
            hi_cdist=hi_cdist,
            sides=("left", "right"),
            traversable=True,
            fundamental_span=(lo_cdist, hi_cdist),
        )
        node.elements = {"left": left_element, "right": right_element}
        node.faces = {"left": left_face, "right": right_face}
        left_face.stable_nodes.append((node, "left"))
        right_face.stable_nodes.append((node, "right"))
        self.nodes.append(node)
        return node

    def wall(
        self,
        left_face: FaceNode,
        right_face: FaceNode,
        left_element: ElementRef,
        right_element: ElementRef,
        *,
        lo_cdist: float = 0.0,
        hi_cdist: float = 1.0,
    ) -> tuple[StableNode, StableNode]:
        """Add the two open side nodes of one wall edge."""
        branch_key = left_element.branch_key
        edge_key = self._edge_key(branch_key)
        made: list[StableNode] = []
        for side, face, element in (
            ("left", left_face, left_element),
            ("right", right_face, right_element),
        ):
            node = StableNode(
                arc=None,
                edge_key=edge_key,
                key=(edge_key, side),
                branch_key=branch_key,
                lo_cdist=lo_cdist,
                hi_cdist=hi_cdist,
                sides=(side,),
                traversable=False,
                fundamental_span=None,
            )
            node.elements = {side: element}
            node.faces = {side: face}
            face.stable_nodes.append((node, side))
            self.nodes.append(node)
            made.append(node)
        return made[0], made[1]

    def stable_nodes_of(self, element: ElementRef) -> list[StableNode]:
        """Every node whose element on ``element.side`` is ``element``."""
        return [
            node
            for node in self.nodes
            if element.side in node.sides and node.elements[element.side] == element
        ]

"""
The minimal trellis: the bridges needed to express a trellis's topology, no more.

A grown trellis carries every bridge the numerics happened to cut. Most of them
are redundant for the combinatorics: the stable-manifold partition is defined by
the bridges that carry holes, and the only other bridges the dual graph needs are
the first forward images of those hole bridges whose bridge class is still
active. So the minimal trellis is

* every bridge with a hole attached (``Hole.bounding_ids`` names it), and
* for each such bridge whose class is ACTIVE, the bridges tiling its first
  forward image, read from the iterate table alone.

An inert hole bridge is kept but not mapped forward. The stable manifold is kept
whole; the crossings that survive as nodes are the kept bridges' endpoints, the
anchors and the outermost crossing of each stable branch (the ends of the
partitioned stretch).

The result is a synthetic :class:`~tanglepack.topology.Trellis.Trellis` sharing
the parent's registry, manifolds and generation but carrying only the kept
crossings on its branches and only the kept bridges. Its arrangement is built in
sparse mode (:meth:`~tanglepack.topology.Arrangement.Arrangement.from_trellis`
with ``sparse=True``) so that an absent bridge is an absent wall, not a slit.

Dev Notes:

* :func:`image_chain` is a verbatim copy of ``BridgeClass._image_chain`` from
  the concurrent ``fix-backward-hole-propagation`` work (2026-09-16), which is
  not on this branch. Delete this copy and import that function once the two
  branches merge; ``tests/test_minimal_trellis.py`` compares the two whenever
  both exist.
* A subdivided image (the parent bridge crossed the trimmed-away part of the
  stable manifold, so its image crosses the computed part in the middle)
  contributes EVERY consecutive pair of the chain. A pair no ``Bridge`` object
  spans (the grown manifold stops short of it, or a blast child was never cut
  there) is skipped with a WARNING: the minimal trellis is evidence-based, like
  inertness, and grows as the trellis does.
* "Active" is read from the supplied
  :class:`~tanglepack.topology.BridgeClass.BridgeClassTable` and nothing else.
  Its definition is changing in the concurrent work (from "holds a loop" to
  "maps to no active class"); only ``entry.active`` is consumed here so both
  definitions plug in.
* Anchors are kept as nodes even when no kept bridge ends there: a periodic
  point is the crossing of its own two manifolds and the stable branch's
  partition starts at it. The outermost crossing is kept because the
  partition ends at it; a dropped outer crossing would lose the last stretch
  of stable manifold from every face.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Iterable, Optional, TYPE_CHECKING

from .Arrangement import Arrangement
from .Trellis import Trellis
from .TrellisBranch import TrellisBranch

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ..numerics.Bridge import Bridge, BridgeId
    from ..numerics.Intersection import ManifoldKey
    from .BridgeClass import BridgeClassTable
    from .TopologyResults import Hole

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


def image_chain(trellis: Trellis, bridge_id: "BridgeId") -> Optional[list["BridgeId"]]:
    """
    The consecutive crossing pairs tiling a bridge's forward image, by lookup only.

    The image of the arc between crossings ``a`` and ``b`` is the arc between
    ``f(a)`` and ``f(b)``, read from the iterate table
    (:meth:`~.Trellis.Trellis.iterate`); the crossings on it are every
    registered crossing of the image's unstable branch whose unstable
    canonical distance lies between the two. No ``Bridge`` object is required:
    an image pair the trellis never cut (the grown manifold stops short of it)
    is returned all the same, because its identity is readable from its ids.

    Args:
        trellis: The Trellis whose registry and branches to read.
        bridge_id: The bridge whose image is wanted.

    Returns:
        The image's consecutive pairs in increasing unstable canonical
        distance (each ordered like a
        :data:`~tanglepack.numerics.Bridge.BridgeId`), or ``None`` when an
        endpoint has no registered ``+1`` iterate, the two images carry no
        common unstable branch key, or that branch is not in the trellis.

    Note:
        Verbatim copy of ``BridgeClass._image_chain`` (see the module Dev
        Notes); to be replaced by an import after the branches merge.
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


@dataclass
class MinimalTrellis:
    """
    A trellis reduced to its hole bridges and the images of the active ones.

    Built by :func:`minimal_trellis`. The full snapshot stays reachable as
    :attr:`trellis` (registry lookups, iterates, strong pips all read from it);
    :attr:`sparse` is the synthetic snapshot the arrangement is built over.

    Attributes:
        trellis: The full trellis this was reduced from.
        hole_bridge_ids: Every bridge with a hole attached, in branch-then-cdist
            order. These define the homotopy partition.
        inert_hole_bridge_ids: The hole bridges whose class is inert: kept, not
            mapped forward.
        unclassed_hole_bridge_ids: Hole bridges the class table does not know
            (logged as a WARNING at build time): kept, not mapped forward.
        image_chains: For each hole bridge that WAS mapped forward, the pairs of
            its image chain that were kept (hole bridges included when the image
            lands on one).
        unmapped: The active hole bridges whose image could not be read from
            the iterate table (no ``+1`` entry for an endpoint).
        skipped_pairs: ``(hole bridge, image pair)`` for every image pair with
            no ``Bridge`` object behind it (logged as a WARNING).
        image_bridge_ids: The kept image pairs that are NOT hole bridges: the
            bridges that refine the homotopy partition into the iterated one.
        kept: The kept ``Bridge`` objects (hole bridges and image bridges).
        dropped_bridge_ids: The ids of every non-partial bridge of the full
            trellis that was not kept.
        sparse: The synthetic filtered trellis (kept crossings and kept
            bridges only, sharing the parent's registry and manifolds).
    """

    trellis: Trellis
    hole_bridge_ids: list["BridgeId"]
    inert_hole_bridge_ids: list["BridgeId"]
    unclassed_hole_bridge_ids: list["BridgeId"]
    image_chains: dict["BridgeId", list["BridgeId"]]
    unmapped: list["BridgeId"]
    skipped_pairs: list[tuple["BridgeId", "BridgeId"]]
    image_bridge_ids: list["BridgeId"]
    kept: list["Bridge"]
    dropped_bridge_ids: list["BridgeId"]
    sparse: Trellis
    _arrangement: Optional[Arrangement] = field(default=None, repr=False)

    # ── membership ──────────────────────────────────────────────────────────

    @property
    def kept_bridge_ids(self) -> list["BridgeId"]:
        """The ids of every kept bridge, hole bridges first then image bridges."""
        return list(self.hole_bridge_ids) + list(self.image_bridge_ids)

    def is_hole_bridge(self, bridge_id: "BridgeId") -> bool:
        """Whether a bridge (in either endpoint order) carries a hole."""
        return frozenset(bridge_id) in self._hole_set

    def is_image_bridge(self, bridge_id: "BridgeId") -> bool:
        """Whether a bridge (in either endpoint order) is a kept image pair."""
        return frozenset(bridge_id) in self._image_set

    def is_kept(self, bridge_id: "BridgeId") -> bool:
        """Whether a bridge (in either endpoint order) is in the minimal trellis."""
        return self.is_hole_bridge(bridge_id) or self.is_image_bridge(bridge_id)

    def parents_of(self, image_bridge_id: "BridgeId") -> list["BridgeId"]:
        """
        The hole bridges whose image chain contains one image bridge.

        Args:
            image_bridge_id: A kept image pair, in either endpoint order.

        Returns:
            The parent hole bridge ids (normally one; two when two hole
            bridges have images tiling the same pair).
        """
        wanted = frozenset(image_bridge_id)
        return [
            parent
            for parent, chain in self.image_chains.items()
            if any(frozenset(pair) == wanted for pair in chain)
        ]

    def kept_ids(self, branch_key: "ManifoldKey") -> list[int]:
        """
        The crossings kept as nodes on one branch, anchor outward.

        Args:
            branch_key: A stable or unstable branch key of the full trellis.

        Returns:
            The kept registry ids in canonical-distance order (empty for a
            branch the full trellis does not carry).
        """
        branch = self.sparse.branch(branch_key)
        return [] if branch is None else list(branch.intersection_ids)

    @property
    def _hole_set(self) -> set[frozenset]:
        return {frozenset(bid) for bid in self.hole_bridge_ids}

    @property
    def _image_set(self) -> set[frozenset]:
        return {frozenset(bid) for bid in self.image_bridge_ids}

    # ── arrangement ─────────────────────────────────────────────────────────

    @property
    def arrangement(self) -> Arrangement:
        """
        The sparse arrangement of the kept crossings and bridges (built once).

        Returns:
            :meth:`~tanglepack.topology.Arrangement.Arrangement.from_trellis`
            over :attr:`sparse` with ``sparse=True``.
        """
        if self._arrangement is None:
            self._arrangement = Arrangement.from_trellis(self.sparse, sparse=True)
        return self._arrangement

    # ── reporting ───────────────────────────────────────────────────────────

    def summary(self) -> str:
        """A one-line description of what was kept and what was dropped."""
        return (
            f"MinimalTrellis: {len(self.hole_bridge_ids)} hole bridge(s) "
            f"({len(self.inert_hole_bridge_ids)} inert, "
            f"{len(self.unclassed_hole_bridge_ids)} unclassed), "
            f"{len(self.image_bridge_ids)} image bridge(s) from "
            f"{len(self.image_chains)} chain(s); {len(self.unmapped)} unmapped, "
            f"{len(self.skipped_pairs)} pair(s) skipped; "
            f"{len(self.dropped_bridge_ids)} bridge(s) dropped"
        )

    def describe(self) -> str:
        """A multi-line report: each kept bridge with its role and provenance."""
        lines = [self.summary()]
        for bid in self.hole_bridge_ids:
            role = (
                "inert"
                if frozenset(bid) in {frozenset(b) for b in self.inert_hole_bridge_ids}
                else "unclassed"
                if frozenset(bid) in {frozenset(b) for b in self.unclassed_hole_bridge_ids}
                else "active"
            )
            chain = self.image_chains.get(bid)
            image = "" if chain is None else f" -> image {chain}"
            lines.append(f"  hole bridge {bid} ({role}){image}")
        for bid in self.image_bridge_ids:
            lines.append(f"  image bridge {bid} from {self.parents_of(bid)}")
        for parent, pair in self.skipped_pairs:
            lines.append(f"  skipped image pair {pair} of {parent}: no Bridge object")
        for bid in self.unmapped:
            lines.append(f"  unmapped hole bridge {bid}: no registered +1 iterate")
        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"<{self.summary()}>"


def _bridge_sort_key(trellis: Trellis, bridge_id: "BridgeId") -> tuple:
    """Order bridges by unstable branch (trellis order) then unstable cdist."""
    order = {branch.key: index for index, branch in enumerate(trellis.unstable_branches)}
    lo = trellis.intersection(bridge_id[0])
    return (order.get(lo.manifold_a_key, len(order)), float(lo.unstable_cdist))


def minimal_trellis(
    trellis: Trellis,
    holes: Iterable["Hole"],
    classes: "BridgeClassTable",
) -> MinimalTrellis:
    """
    Reduce a trellis to its hole bridges and the first images of the active ones.

    Args:
        trellis: The full trellis (normally the ALL-fixed-points snapshot, so
            that heteroclinic bridges and every stable branch are present).
        holes: The punched holes whose ``bounding_ids`` name the hole bridges;
            normally every cached per-fixed-point trellis's ``holes``.
        classes: The bridge-class table over the homotopy partition; each hole
            bridge's ``entry.active`` decides whether it is mapped forward.

    Returns:
        The :class:`MinimalTrellis`.

    Raises:
        AssertionError: If a kept bridge's endpoints are not kept crossings,
            or a partial bridge slipped into the kept set.

    Note:
        Pure lookups: the iterate table and the registry are read, the map is
        never called. A hole whose bounding pair no bridge spans, a hole
        bridge the class table does not know, and an image pair no bridge
        spans are each logged as a WARNING and skipped, never invented.
    """
    tol = trellis.registry.cdist_tol

    # 1. Hole bridges.
    hole_ids: dict[frozenset, "BridgeId"] = {}
    for hole in holes:
        if hole.bounding_ids is None:
            continue
        bridge = trellis.bridge_between(*hole.bounding_ids)
        if bridge is None or bridge.id is None:
            logger.warning(
                "hole near crossing %d names bounding pair %s, which no bridge of "
                "this trellis spans; it contributes no bridge to the minimal trellis",
                hole.near_intersection_id,
                hole.bounding_ids,
            )
            continue
        hole_ids.setdefault(frozenset(bridge.id), bridge.id)
    hole_bridge_ids = sorted(hole_ids.values(), key=lambda b: _bridge_sort_key(trellis, b))

    # 2. Images of the active hole bridges.
    inert: list["BridgeId"] = []
    unclassed: list["BridgeId"] = []
    chains: dict["BridgeId", list["BridgeId"]] = {}
    unmapped: list["BridgeId"] = []
    skipped: list[tuple["BridgeId", "BridgeId"]] = []
    image_ids: dict[frozenset, "BridgeId"] = {}
    for bid in hole_bridge_ids:
        try:
            entry = classes.entry_of(bid)
        except KeyError:
            logger.warning(
                "hole bridge %s is in no bridge class; kept, not mapped forward", bid
            )
            unclassed.append(bid)
            continue
        if not entry.active:
            inert.append(bid)
            continue
        chain = image_chain(trellis, bid)
        if chain is None:
            logger.debug(
                "active hole bridge %s has no registered image; not mapped", bid
            )
            unmapped.append(bid)
            continue
        kept_pairs: list["BridgeId"] = []
        for pair in chain:
            if frozenset(pair) in hole_ids:
                kept_pairs.append(hole_ids[frozenset(pair)])
                continue
            bridge = trellis.bridge_between(*pair)
            if bridge is None or bridge.id is None:
                logger.warning(
                    "image pair %s of hole bridge %s has no Bridge object; skipped "
                    "(grow or blast the trellis to register it)",
                    pair,
                    bid,
                )
                skipped.append((bid, pair))
                continue
            image_ids.setdefault(frozenset(bridge.id), bridge.id)
            kept_pairs.append(bridge.id)
        chains[bid] = kept_pairs
    image_bridge_ids = sorted(
        image_ids.values(), key=lambda b: _bridge_sort_key(trellis, b)
    )

    # 3. Kept bridges and kept crossings.
    kept_ids = hole_bridge_ids + image_bridge_ids
    kept: list["Bridge"] = []
    for bid in kept_ids:
        bridge = trellis.bridge_between(*bid)
        assert bridge is not None and bridge.id is not None
        kept.append(bridge)
    kept_set = {frozenset(bid) for bid in kept_ids}
    dropped = sorted(
        (
            bridge.id
            for bridge in trellis.bridges
            if bridge.id is not None and frozenset(bridge.id) not in kept_set
        ),
        key=lambda b: _bridge_sort_key(trellis, b),
    )

    endpoints = {endpoint for bid in kept_ids for endpoint in bid}
    branches: dict["ManifoldKey", TrellisBranch] = {}
    for key, branch in trellis.branches.items():
        if branch.stability == "stable":
            attr, key_attr = "stable_cdist", "manifold_b_key"
        else:
            attr, key_attr = "unstable_cdist", "manifold_a_key"
        ids = {
            iid
            for iid in endpoints
            if getattr(trellis.intersection(iid), key_attr) == key
        }
        for iid in branch.intersection_ids:
            if float(getattr(trellis.intersection(iid), attr)) <= tol:
                ids.add(iid)  # the anchor: a periodic point is its own crossing
        if branch.stability == "stable" and len(branch):
            ids.add(branch.intersection_ids[-1])  # the partition ends here
        assert ids <= set(branch.intersection_ids), (
            f"kept crossings {sorted(ids - set(branch.intersection_ids))} are "
            f"not on branch {key[1:]} of the full trellis"
        )
        branches[key] = TrellisBranch(
            key=key,
            fixed_point=branch.fixed_point,
            stability=branch.stability,
            orbit_index=branch.orbit_index,
            branch_index=branch.branch_index,
            intersection_ids=sorted(
                ids, key=lambda iid: float(getattr(trellis.intersection(iid), attr))
            ),
        )
    node_ids = {iid for branch in branches.values() for iid in branch.intersection_ids}
    assert endpoints <= node_ids, (
        f"kept bridge endpoints {sorted(endpoints - node_ids)} are on no branch"
    )

    sparse = Trellis(
        fixed_points=list(trellis.fixed_points),
        registry=trellis.registry,
        branches=branches,
        bridges=kept,
        dynamical_system=trellis.dynamical_system,
        manifolds=trellis.manifolds,
        generation=trellis._built_generation,
    )
    sparse.stable_partitions = list(trellis.stable_partitions)

    minimal = MinimalTrellis(
        trellis=trellis,
        hole_bridge_ids=hole_bridge_ids,
        inert_hole_bridge_ids=inert,
        unclassed_hole_bridge_ids=unclassed,
        image_chains=chains,
        unmapped=unmapped,
        skipped_pairs=skipped,
        image_bridge_ids=image_bridge_ids,
        kept=kept,
        dropped_bridge_ids=dropped,
        sparse=sparse,
    )
    logger.info("%s", minimal.summary())
    return minimal

"""
Resolution diagnostics: where the polylines are too coarse to trust.

A computed tangle is a set of polylines standing for smooth curves, and two
things go silently wrong where the polylines are too coarse for the geometry:

* **Unresolved crossings.** Near a tangency an unstable and a stable polyline
  can straddle each other without the curves crossing there (or cross without
  the polylines showing it). ``Tangle`` scores every detected crossing with
  its :attr:`~.Intersection.Intersection.resolution_ratio` and chord deviation;
  :meth:`~.IntersectionRegistry.IntersectionRegistry.unresolved_crossings`
  lists the detections that are not crossings of N: below the threshold
  (clause 1), or in a run of detections closer along the unstable polyline
  than their chord deviation (clause 2).
* **Same-stability crossings.** By the fundamental invariant two unstable
  curves never cross, nor do two stable ones, and no curve crosses itself. A
  polyline that does -- a Z-kink of one manifold, two branches, a blast's image
  bridge over another curve -- is a polygon ARTIFACT, never a real crossing.
  ``Tangle`` never pairs two segments of one stability, so it cannot see
  them; :func:`same_stability_crossings` finds them, and
  :func:`count_same_stability_crossings` counts the ones a blast's new image
  bridges take part in (a blast adds no stable curve, and every hit it can
  add involves one of its images).

Nothing here changes the tangle. An unresolved detection is out of N by the
registry's definition, not by this report, which only names it;
a same-stability crossing is a FLAG. :meth:`~.TangleWorkbench.TangleWorkbench.check_resolution`
bundles the two into a :class:`ResolutionReport`.

Dev Notes:

The tangle region is the box of the registered crossings padded by ``pad`` of
its larger side, the rule of the test-tier law
``tests/helpers/laws.check_no_same_stability_crossing`` (which stays the
independent oracle, deliberately not built on this module). Outside it an
escaping arm is under-resolved by design (the inversion saddle's escaping arm
crosses itself thousands of times) and says nothing about the dynamics.

``suggested_area_cutoff`` is half the larger curvature area of the two
crossing segments: refinement only splits a segment whose area reaches
``area_cutoff``, so that is NECESSARY for the refiner to touch the kink. It is
not sufficient -- a kink can survive one level of refinement.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, Literal, Optional

import numpy as np
from numpy.typing import NDArray

from .geometry import curvature_area_batch
from .Intersection import Stability

if TYPE_CHECKING:
    from .BranchPoint import BranchPoint
    from .Bridge import Bridge
    from .FixedPoint import FixedPoint
    from .PointStore import PointStore
    from .TangleWorkbench import TangleWorkbench

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Interior parameter margin: two segments that only touch at an endpoint do not
# cross (the Tangle's own predicate counts touching endpoints).
_MARGIN = 1e-9
# Relative parallel threshold on sin(angle), as in Tangle.
_EPS_REL = 1e-12
# Candidate pairs generated at once by the sweep.
_PAIR_BLOCK = 2_000_000


@dataclass(frozen=True)
class SameStabilityCrossing:
    """
    Two same-stability polyline segments crossing: a polygon artifact.

    Attributes:
        stability: The stability of both curves.
        curve_a: The first curve's name (a manifold key, or an image bridge).
        curve_b: The second curve's name; equal to ``curve_a`` for a self-crossing.
        seg_a: The segment index along ``curve_a`` (from its root).
        seg_b: The segment index along ``curve_b``.
        point: Where the two segments cross.
        sin: ``|sin|`` of the crossing angle.
        local_area: The larger curvature area of the two segments.
        suggested_area_cutoff: ``0.5 * local_area``, the most the refiner's
            ``area_cutoff`` can be for it to split these segments at all
            (necessary, not sufficient; see the module Dev Notes).
        in_tangle_region: Whether ``point`` lies in the padded box of the
            registered crossings.
    """

    stability: Stability
    curve_a: str
    curve_b: str
    seg_a: int
    seg_b: int
    point: tuple[float, float]
    sin: float
    local_area: float
    suggested_area_cutoff: float
    in_tangle_region: bool


@dataclass(frozen=True)
class SameStabilityCount:
    """
    The unstable polyline crossings some new curves take part in, counted.

    Attributes:
        in_region: Hits inside the padded box of the registered crossings.
        outside: Hits outside it.
        suggested_area_cutoff: Half the smallest in-region local area (see
            :class:`SameStabilityCrossing`); None without an in-region hit.
        curve_pairs: The distinct in-region curve pairs, ``"a x b"``, sorted.
    """

    in_region: int
    outside: int
    suggested_area_cutoff: Optional[float]
    curve_pairs: tuple[str, ...]


@dataclass(frozen=True)
class ResolutionReport:
    """
    Everything the polylines fail to resolve, on demand.

    Attributes:
        threshold: The resolution-ratio threshold used (None: clause 1 off).
        unresolved_crossings: Registry ids of the detections that are not
            crossings of N -- ratio below ``threshold`` (clause 1) or in a run
            (clause 2) -- in unstable cdist order.
        same_stability_crossings: Every same-stability polyline crossing (empty
            when only a count was taken).
        same_stability_count: The count of the hits a blast's new image bridges
            take part in, when the report was a blast's; else None.
    """

    threshold: Optional[float]
    unresolved_crossings: tuple[int, ...]
    same_stability_crossings: tuple[SameStabilityCrossing, ...]
    same_stability_count: Optional[SameStabilityCount] = None

    @property
    def clean(self) -> bool:
        """No unresolved crossing and no same-stability crossing in the tangle region."""
        return (
            not self.unresolved_crossings
            and not any(hit.in_tangle_region for hit in self.same_stability_crossings)
            and not (self.same_stability_count and self.same_stability_count.in_region)
        )


def same_stability_crossings(
    workbench: "TangleWorkbench",
    fixed_points: Optional[Iterable["FixedPoint"]] = None,
    *,
    region: Literal["tangle", "all"] = "tangle",
    pad: float = 0.1,
) -> list[SameStabilityCrossing]:
    """
    Every crossing of two unstable, or two stable, polylines.

    The curves are every manifold of the selected fixed points plus the
    unstable image bridges that share no node with a manifold (a blast's
    images; a bridge cut out of a manifold IS that manifold). Two segments that
    share a node -- neighbours on one curve, two branches leaving one root, two
    image bridges overlapping at a crossing -- never count, nor do segments that
    only touch.

    Args:
        workbench: The workbench holding the curves.
        fixed_points: Whose curves to check; None checks every fixed point.
        region: ``"tangle"`` marks a hit ``in_tangle_region`` only inside the
            padded box of the registered crossings (none when nothing is
            registered); ``"all"`` counts the whole plane as the tangle region.
        pad: The box padding, as a fraction of its larger side.

    Returns:
        The hits, unstable ones first, each curve pair in curve order.
    """
    selected = (
        list(workbench.fixed_points) if fixed_points is None else list(fixed_points)
    )
    box = (
        _tangle_box(workbench, selected, pad)
        if region == "tangle"
        else (np.full(2, -np.inf), np.full(2, np.inf))
    )
    return [
        hit
        for stability in ("unstable", "stable")
        for hit in _crossings_among(_curves(workbench, selected, stability), stability, box)
    ]


def count_same_stability_crossings(
    workbench: "TangleWorkbench", bridges: Iterable["Bridge"], *, pad: float = 0.1
) -> SameStabilityCount:
    """
    Count the unstable polyline crossings that involve one of ``bridges``.

    The field is every unstable curve of the bridges' fixed points, as in
    :func:`same_stability_crossings`; only the segment pairs with a segment of
    one of ``bridges`` are built and tested, and nothing is materialised per
    hit. A blast still pays for the whole field once: the segment rows, the
    edge dedup and the sweep's sort run over every unstable segment of those
    fixed points; the pairs are built for the new segments only.

    Args:
        workbench: The workbench holding the curves.
        bridges: The new curves (a blast's image bridges).
        pad: The tangle-box padding, as in :func:`same_stability_crossings`.

    Returns:
        The count; the hits are exactly those of the full listing naming one of
        ``bridges``.
    """
    bridges = list(bridges)
    fixed_points = list({id(bridge.fixed_point): bridge.fixed_point for bridge in bridges}.values())
    curves = _curves(workbench, fixed_points, "unstable")
    new = {id(bridge) for bridge in bridges}
    owner_a, owner_b, *_rest, local_area, inside = _hits(
        curves,
        _tangle_box(workbench, fixed_points, pad),
        new=np.array([id(bridge) in new for _name, _read, bridge in curves], dtype=bool),
    )
    pairs = np.unique(np.stack([owner_a[inside], owner_b[inside]], axis=1), axis=0)
    return SameStabilityCount(
        in_region=int(inside.sum()),
        outside=int((~inside).sum()),
        suggested_area_cutoff=0.5 * float(local_area[inside].min()) if inside.any() else None,
        curve_pairs=tuple(sorted(f"{curves[p][0]} x {curves[q][0]}" for p, q in pairs)),
    )


def _tangle_box(
    workbench: "TangleWorkbench", fixed_points: list["FixedPoint"], pad: float
) -> Optional[tuple[NDArray[np.float64], NDArray[np.float64]]]:
    """The padded box of the registered crossings of these fixed points."""
    points = np.array(
        [
            ix.coords
            for _iid, ix in workbench.intersection_registry
            if any(fp is chosen for fp in ix.fixed_points for chosen in fixed_points)
        ],
        dtype=float,
    )
    if not len(points):
        return None
    lo, hi = points.min(axis=0), points.max(axis=0)
    margin = pad * float((hi - lo).max())
    return lo - margin, hi + margin


def _curves(
    workbench: "TangleWorkbench", fixed_points: list["FixedPoint"], stability: Stability
) -> list[tuple[str, tuple, Optional["Bridge"]]]:
    """The named reads ``(periodic point or None, store, rows)`` of one stability's
    curves, root to tail, each with its image bridge (None for a manifold)."""
    reads = {
        key: manifold._rows()
        for key, manifold in workbench.manifolds.items()
        if key[1] == stability
    }
    curves = [
        (_manifold_name(key), read, None)
        for key, read in reads.items()
        if any(key[0] is fp for fp in fixed_points)
    ]
    if stability == "stable":
        return curves
    # A bridge sharing a node with a manifold is cut from it: a row some manifold
    # of its store holds, or the periodic point some manifold leaves from.
    held: dict[PointStore, NDArray[np.bool_]] = {}
    for _head, store, rows in reads.values():
        if store is not None:
            held.setdefault(store, np.zeros(len(store), dtype=bool))[rows] = True
    heads = {head for head, _store, _rows in reads.values() if head is not None}
    for bridge in workbench.bridges:
        if not any(bridge.fixed_point is fp for fp in fixed_points):
            continue
        head, store, rows = read = bridge._rows()
        if head in heads or (store in held and held[store][rows].any()):
            continue
        name = f"image bridge {bridge.id}" if bridge.id is not None else f"partial image bridge {id(bridge):#x}"
        curves.append((name, read, bridge))
    return curves


def _manifold_name(key: tuple) -> str:
    """A manifold key as text: ``A unstable 0.1``."""
    fixed_point, stability, orbit, branch = key
    return f"{fixed_point.label or '?'} {stability} {orbit}.{branch}"


def _crossings_among(
    curves: list[tuple[str, tuple, Optional["Bridge"]]],
    stability: Stability,
    box: Optional[tuple[NDArray[np.float64], NDArray[np.float64]]],
) -> list[SameStabilityCrossing]:
    """The transversal interior crossings among the segments of some curves."""
    return [
        SameStabilityCrossing(
            stability=stability,
            curve_a=curves[p][0],
            curve_b=curves[q][0],
            seg_a=int(seg_p),
            seg_b=int(seg_q),
            point=(float(x), float(y)),
            sin=float(sin),
            local_area=float(area),
            suggested_area_cutoff=0.5 * float(area),
            in_tangle_region=bool(inside),
        )
        for p, q, seg_p, seg_q, (x, y), sin, area, inside in zip(*_hits(curves, box))
    ]


def _hits(
    curves: list[tuple[str, tuple, Optional["Bridge"]]],
    box: Optional[tuple[NDArray[np.float64], NDArray[np.float64]]],
    new: Optional[NDArray[np.bool_]] = None,
) -> tuple[NDArray, ...]:
    """The transversal interior crossings among the segments of some curves, as arrays.

    Returns, one entry per hit: the two segments' curve indices and positions
    along their curves, the crossing point, ``|sin|`` of the angle, the larger
    curvature area of the two segments, and whether the point lies in ``box``
    (never, for None). With ``new`` (one flag per curve) only the segment pairs
    with a segment of a new curve are tested.
    """
    stores = list(dict.fromkeys(store for _name, (_head, store, rows), _bridge in curves if len(rows)))
    offsets = dict(zip(stores, np.cumsum([0, *(len(store) for store in stores[:-1])]).tolist()))
    periodic: dict[BranchPoint, int] = {}
    rows = [_segment_rows(k, read, offsets, periodic) for k, (_name, read, _bridge) in enumerate(curves)]
    rows = [row for row in rows if row is not None]
    if not rows:
        none = np.empty(0, dtype=np.int64)
        return none, none, none, none, np.empty((0, 2)), np.empty(0), np.empty(0), np.empty(0, dtype=bool)
    a, b, left, right, node_a, node_b, owner, position = (
        np.concatenate(column) for column in zip(*rows)
    )
    # One edge reached from two curves (overlapping image bridges) is one segment.
    edges = np.sort(np.stack([node_a, node_b], axis=1), axis=1)
    _, first = np.unique(edges, axis=0, return_index=True)
    keep = np.sort(first)
    a, b, left, right, node_a, node_b, owner, position = (
        column[keep] for column in (a, b, left, right, node_a, node_b, owner, position)
    )

    i, j = _candidate_pairs(a, b, None if new is None else new[owner])
    shared = (
        (node_a[i] == node_a[j])
        | (node_a[i] == node_b[j])
        | (node_b[i] == node_a[j])
        | (node_b[i] == node_b[j])
    )
    i, j = i[~shared], j[~shared]

    d1, d2 = b[i] - a[i], b[j] - a[j]
    den = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
    scale = np.hypot(*d1.T) * np.hypot(*d2.T)
    r = a[j] - a[i]
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (r[:, 0] * d2[:, 1] - r[:, 1] * d2[:, 0]) / den
        u = (r[:, 0] * d1[:, 1] - r[:, 1] * d1[:, 0]) / den
    hit = (
        (np.abs(den) > _EPS_REL * scale)
        & (t > _MARGIN)
        & (t < 1 - _MARGIN)
        & (u > _MARGIN)
        & (u < 1 - _MARGIN)
    )
    i, j, t = i[hit], j[hit], t[hit]
    sin = np.abs(den[hit]) / scale[hit]
    point = a[i] + t[:, None] * (b[i] - a[i])
    local_area = np.maximum(
        curvature_area_batch(a[i], b[i], left[i], right[i]),
        curvature_area_batch(a[j], b[j], left[j], right[j]),
    )
    in_region = (
        np.zeros(len(i), dtype=bool)
        if box is None
        else np.all((box[0] <= point) & (point <= box[1]), axis=1)
    )
    return owner[i], owner[j], position[i], position[j], point, sin, local_area, in_region


def _segment_rows(
    owner: int,
    read: tuple,
    offsets: dict["PointStore", int],
    periodic: dict["BranchPoint", int],
) -> Optional[tuple[NDArray, ...]]:
    """Per-segment arrays of one curve: ends, outer neighbours, node keys, owner, index.

    A node's key is its row plus its store's offset; a periodic point, which has
    no row, is -1, -2, ... in the order ``periodic`` first meets it (never a row).
    """
    head, store, rows = read
    if len(rows) + (head is not None) < 2:
        return None
    points, keys = store.xy[rows], rows + offsets[store]
    if head is not None:
        points = np.vstack((head.get_point(), points))
        keys = np.concatenate(([-1 - periodic.setdefault(head, len(periodic))], keys))
    missing = np.full((1, 2), np.nan)
    count = len(points) - 1
    return (
        points[:-1],
        points[1:],
        np.vstack([missing, points[:-2]]),
        np.vstack([points[2:], missing]),
        keys[:-1],
        keys[1:],
        np.full(count, owner),
        np.arange(count),
    )


def _candidate_pairs(
    a: NDArray[np.float64], b: NDArray[np.float64], new: Optional[NDArray[np.bool_]] = None
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Every pair of segments whose bounding boxes overlap, by sweep and prune
    (with ``new``, one flag per segment, only the pairs holding a new segment).

    The boxes are sorted along the axis on which fewer of them overlap, each is
    paired with the boxes that start before it ends, and the pairs that also
    overlap on the other axis are kept. With ``new`` an old segment is paired
    only with the new segments inside its forward range, so the pairs are built
    for the new segments alone, in the orientation of the full sweep. A
    tangle's segments are short next to their spacing, so this yields a few
    candidates per segment; it runs in blocks of ``_PAIR_BLOCK`` pairs so a
    dense fold cannot exhaust memory.
    (rtree's bulk ``intersection_v`` is not used: in 1.4.0 it loses its result
    buffer once the results outgrow twice the query count, which every polyline
    does -- each segment meets its two neighbours.)
    """
    mins, maxs = np.minimum(a, b), np.maximum(a, b)
    sweeps = []
    for axis in (0, 1):
        order = np.argsort(mins[:, axis], kind="stable")
        ends = np.searchsorted(mins[order, axis], maxs[order, axis], side="right")
        sweeps.append((axis, order, np.maximum(ends - np.arange(len(order)) - 1, 0)))
    axis, order, counts = min(sweeps, key=lambda sweep: sweep[2].sum())
    other = 1 - axis

    # The partners of sorted position p are pool[start[p] : start[p] + count[p]]:
    # its whole forward range, or for an old p only the new positions in it.
    positions = np.arange(len(order))
    if new is None:
        pool, start, count = positions, positions + 1, counts
    else:
        is_new = new[order]
        new_positions = np.flatnonzero(is_new)
        lo = np.searchsorted(new_positions, positions + 1)
        hi = np.searchsorted(new_positions, positions + 1 + counts)
        pool = np.concatenate([positions, new_positions])
        start = np.where(is_new, positions + 1, len(order) + lo)
        count = np.where(is_new, counts, hi - lo)

    totals = np.cumsum(count)
    cuts = np.searchsorted(totals, np.arange(_PAIR_BLOCK, totals[-1], _PAIR_BLOCK))
    firsts, seconds = [], []
    for block in np.split(positions, cuts):
        block_counts = count[block]
        first = np.repeat(block, block_counts)
        offset = np.arange(block_counts.sum()) - np.repeat(
            np.cumsum(block_counts) - block_counts, block_counts
        )
        i, j = order[first], order[pool[np.repeat(start[block], block_counts) + offset]]
        overlap = (mins[i, other] <= maxs[j, other]) & (mins[j, other] <= maxs[i, other])
        firsts.append(i[overlap])
        seconds.append(j[overlap])
    return np.concatenate(firsts), np.concatenate(seconds)

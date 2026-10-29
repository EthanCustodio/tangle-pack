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
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from typing import TYPE_CHECKING, Callable, Iterable, Literal, Optional

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
# Sweep strips are this many median extents of the new segments wide, at most _MAX_STRIPS.
_STRIP_WIDTH = 4.0
_MAX_STRIPS = 4096
# Node keys shift by this into an edge key: the periodic points are -1, -2, ...
_KEY_OFFSET = 1 << 20


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
        The hits, unstable ones first, in sweep order: each names first the
        segment earlier along the sweep axis (not always the earlier curve's),
        and they run by that segment, then the other.
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
        for hit in _crossings_among(
            _curves(workbench, selected, stability), stability, box, workbench.workers
        )
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
    edge dedup, the sweep axis and sort and the choice of strips run over every
    unstable segment of those fixed points; the pairs are built, and tested,
    for the new segments only, in the strips they span.

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
        workers=workbench.workers,
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
    workers: int = 1,
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
        for p, q, seg_p, seg_q, (x, y), sin, area, inside in zip(*_hits(curves, box, workers=workers))
    ]


def _hits(
    curves: list[tuple[str, tuple, Optional["Bridge"]]],
    box: Optional[tuple[NDArray[np.float64], NDArray[np.float64]]],
    new: Optional[NDArray[np.bool_]] = None,
    workers: int = 1,
) -> tuple[NDArray, ...]:
    """The transversal interior crossings among the segments of some curves, as arrays.

    Returns, one entry per hit: the two segments' curve indices and positions
    along their curves, the crossing point, ``|sin|`` of the angle, the larger
    curvature area of the two segments, and whether the point lies in ``box``
    (never, for None). A hit names first the segment earlier along the sweep
    axis. Without ``new`` the hits come in sweep order (by first, then second
    segment); with ``new`` (one flag per curve) only the segment pairs with a
    segment of a new curve are tested, and the hits come strip by strip.
    ``workers`` threads test the candidate blocks; the hits are the same.
    """
    stores = list(dict.fromkeys(store for _name, (_head, store, rows), _bridge in curves if len(rows)))
    offsets = dict(zip(stores, np.cumsum([0, *(len(store) for store in stores[:-1])]).tolist()))
    periodic: dict[BranchPoint, int] = {}
    rows = [_segment_rows(k, read, offsets, periodic) for k, (_name, read, _bridge) in enumerate(curves)]
    rows = [row for row in rows if row is not None]
    if not rows:
        none = np.empty(0, dtype=np.int64)
        return none, none, none, none, np.empty((0, 2)), np.empty(0), np.empty(0), np.empty(0, dtype=bool)
    a, b, left, right, node_a, node_b, owner, position = (np.concatenate(column) for column in zip(*rows))
    # One edge reached from two curves (overlapping image bridges) is one segment.
    _, first = np.unique(_edge_keys(node_a, node_b), return_index=True)
    segment = np.sort(first)
    mins, maxs = np.minimum(a[segment], b[segment]), np.maximum(a[segment], b[segment])
    axis = _sweep_axis(mins, maxs)
    sweep = np.argsort(mins[:, axis], kind="stable")
    # From here on a segment is its position along the sweep; ``order`` gives its row.
    order = segment[sweep]
    mins, maxs = mins[sweep], maxs[sweep]
    a, b, node_a, node_b = a[order], b[order], node_a[order], node_b[order]
    is_new = np.ones(len(order), dtype=bool) if new is None else new[owner[order]]

    blocks, pairs_in = _candidate_pairs(mins, maxs, is_new, axis, workers)

    def crossings_in(block: NDArray[np.int64]) -> tuple[NDArray, ...]:
        """The transversal interior crossings among one block's candidate pairs."""
        i, j = pairs_in(block)
        shared = (node_a[i] == node_a[j]) | (node_a[i] == node_b[j])
        shared |= (node_b[i] == node_a[j]) | (node_b[i] == node_b[j])
        i, j = i[~shared], j[~shared]
        d1, d2 = b[i] - a[i], b[j] - a[j]
        den = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
        scale = np.hypot(*d1.T) * np.hypot(*d2.T)
        r = a[j] - a[i]
        with np.errstate(divide="ignore", invalid="ignore"):
            t = (r[:, 0] * d2[:, 1] - r[:, 1] * d2[:, 0]) / den
            u = (r[:, 0] * d1[:, 1] - r[:, 1] * d1[:, 0]) / den
        hit = (np.abs(den) > _EPS_REL * scale) & (t > _MARGIN) & (t < 1 - _MARGIN)
        hit &= (u > _MARGIN) & (u < 1 - _MARGIN)
        return i[hit], j[hit], t[hit], np.abs(den[hit]) / scale[hit]

    # The blocks are independent, so they run on ``workers`` threads (numpy
    # releases the GIL in their array work), in block order either way.
    found = [(np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64), np.empty(0), np.empty(0))]
    if workers > 1 and len(blocks) > 1:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            found += pool.map(crossings_in, blocks)
    else:
        found += map(crossings_in, blocks)
    i, j, t, sin = (np.concatenate(column) for column in zip(*found))
    if new is None:
        # The strips find the pairs strip by strip; the listing is in sweep order.
        by_sweep = np.lexsort((j, i))
        i, j, t, sin = i[by_sweep], j[by_sweep], t[by_sweep], sin[by_sweep]

    point = a[i] + t[:, None] * (b[i] - a[i])
    local_area = np.maximum(
        curvature_area_batch(a[i], b[i], left[order[i]], right[order[i]]),
        curvature_area_batch(a[j], b[j], left[order[j]], right[order[j]]),
    )
    in_region = (
        np.zeros(len(i), dtype=bool)
        if box is None
        else np.all((box[0] <= point) & (point <= box[1]), axis=1)
    )
    i, j = order[i], order[j]
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


def _edge_keys(node_a: NDArray[np.int64], node_b: NDArray[np.int64]) -> NDArray[np.int64]:
    """One integer per edge, the same in both directions: its two node keys,
    shifted past the periodic points' -1, -2, ..., packed low key high."""
    low = np.minimum(node_a, node_b) + _KEY_OFFSET
    high = np.maximum(node_a, node_b) + _KEY_OFFSET
    assert low.min() >= 0 and high.max() < 1 << 31, "node keys outside the edge-key range"
    return (low << 32) | high


def _sweep_axis(mins: NDArray[np.float64], maxs: NDArray[np.float64]) -> int:
    """The axis along which fewer boxes overlap, the first on a tie.

    Sorted by their minima, each box overlaps the boxes after it whose minimum
    is at most its maximum, so with ``n`` boxes the overlapping pairs along an
    axis number ``sum(searchsorted(sorted minima, maxima, right)) - n(n+1)/2``,
    the maxima taken in any order.
    """
    n = len(mins)
    overlaps = [
        int(np.searchsorted(np.sort(mins[:, k]), np.sort(maxs[:, k]), side="right").sum()) - n * (n + 1) // 2
        for k in (0, 1)
    ]
    return int(overlaps[1] < overlaps[0])


def _candidate_pairs(
    mins: NDArray[np.float64],
    maxs: NDArray[np.float64],
    is_new: NDArray[np.bool_],
    axis: int,
    workers: int = 1,
) -> tuple[list[NDArray[np.int64]], Callable[[NDArray[np.int64]], tuple[NDArray, NDArray]]]:
    """Every pair of boxes that overlap and hold a new one, by sweep and prune
    inside strips, as blocks of about ``_PAIR_BLOCK`` candidates (a quarter of
    each worker's share with ``workers`` above 1) and the function that builds
    one block's pairs.
    The boxes come sorted by their minima along ``axis`` (the sweep), and a
    pair ``(p, q)`` of sweep positions has ``p < q``. The other axis is cut
    into strips ``_STRIP_WIDTH`` median extents of the new boxes wide (at most
    ``_MAX_STRIPS``), and only the strips some new box spans are swept. In each,
    a box is paired with the boxes after it that start no later than it ends --
    an old box only with the new ones -- and the pairs that also overlap across
    are kept. Two overlapping boxes both span the strip of the larger of their
    minima across, so each pair is kept there and only there. A blast's images
    lie in dense folds: swept over the whole plane, a new box met about a
    thousand candidates for some 75 overlaps (D1, 2026-10-28); a strip holds
    only the folds near it. The strip width changes the cost, never the pairs.
    (rtree's bulk ``intersection_v`` is no help: in 1.4.0 it loses its result
    buffer once the results outgrow twice the query count, which every polyline
    does -- each segment meets its two neighbours --, and the workbench's index
    holds no blast's image curves.)
    """
    if not is_new.any():
        return [], None
    across = 1 - axis
    low, high = float(mins[:, across].min()), float(maxs[:, across].max())
    width = _STRIP_WIDTH * float(np.median(maxs[is_new, across] - mins[is_new, across]))
    strips = int(min(_MAX_STRIPS, np.ceil((high - low) / width))) if width > 0 and high > low else 1
    width = (high - low) / strips if high > low else 1.0
    first_strip, last_strip = (
        np.clip(np.floor((ends[:, across] - low) / width), 0, strips - 1).astype(np.int64)
        for ends in (mins, maxs)
    )
    spanned = np.cumsum(
        np.bincount(first_strip[is_new], minlength=strips + 1)
        - np.bincount(last_strip[is_new] + 1, minlength=strips + 1)
    )[:strips] > 0
    reached = np.concatenate(([0], np.cumsum(spanned)))
    boxes = np.flatnonzero(reached[last_strip + 1] > reached[first_strip])

    # One entry per box and spanned strip it spans, by strip, in sweep order within one.
    spans = last_strip[boxes] - first_strip[boxes] + 1
    box = np.repeat(boxes, spans)
    strip = np.repeat(first_strip[boxes] - np.cumsum(spans) + spans, spans) + np.arange(spans.sum())
    box, strip = box[spanned[strip]], strip[spanned[strip]]
    by_strip = np.argsort(strip, kind="stable")
    box, strip = box[by_strip], strip[by_strip]
    bounds = np.searchsorted(strip, np.arange(strips + 1))
    ends = np.empty(len(box), dtype=np.int64)
    for s in np.flatnonzero(np.diff(bounds)).tolist():
        lo, hi = bounds[s], bounds[s + 1]
        ends[lo:hi] = lo + np.searchsorted(mins[box[lo:hi], axis], maxs[box[lo:hi], axis], side="right")

    # The partners of entry e are pool[start[e] : start[e] + count[e]]: its whole
    # forward range in its strip, or for an old box only the new entries in it.
    entries = np.arange(len(box))
    entry_new = is_new[box]
    new_entries = np.flatnonzero(entry_new)
    after, until = np.searchsorted(new_entries, entries + 1), np.searchsorted(new_entries, ends)
    pool = np.concatenate([entries, new_entries])
    start = np.where(entry_new, entries + 1, len(box) + after)
    count = np.where(entry_new, ends - entries - 1, until - after)

    totals = np.cumsum(count)
    size = _PAIR_BLOCK if workers == 1 else max(1, min(_PAIR_BLOCK, int(totals[-1]) // (4 * workers)))
    blocks = np.split(entries, np.searchsorted(totals, np.arange(size, totals[-1], size)))

    def pairs_in(block: NDArray[np.int64]) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
        """The overlapping pairs of one block's candidates."""
        block_counts = count[block]
        entry = np.repeat(block, block_counts)
        offset = np.arange(block_counts.sum()) - np.repeat(
            np.cumsum(block_counts) - block_counts, block_counts
        )
        p, q = box[entry], box[pool[np.repeat(start[block], block_counts) + offset]]
        keep = (
            (mins[p, across] <= maxs[q, across])
            & (mins[q, across] <= maxs[p, across])
            & (np.maximum(first_strip[p], first_strip[q]) == strip[entry])
        )
        return p[keep], q[keep]

    return blocks, pairs_in

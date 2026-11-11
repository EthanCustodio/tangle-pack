"""
The trellis cartoon: a trellis redrawn in standard form (Mitchell, Physica D
238, 2009, Fig. 3b).

The computed manifolds hide the topology behind the geometry: long thin lobes,
crossings bunched near the stable branch, huge size ratios. The cartoon keeps
only the ORDER of the crossings. Each fixed point is a circle of the zone
cartoon (:func:`~.plotting._zone_circles`): its stable branches are arcs from
the anchor to the pip end, alternating with empty arcs (the unstable sides),
and a period-1 tangle is the upper half circle. On each branch the crossings
the figure shows get consecutive ordinal positions in stable-cdist order
(:meth:`~tanglepack.numerics.Intersection.Intersection.cdist_key`), the anchor
at rank 0; every bridge is a smooth curve between its two crossings, each end
leaving the arc on the side its row gives (:func:`~.StablePartition.row_of_end`,
crossing signs only). Nothing here calls the map or reads a polyline.

Three routes, by the bridge's two ends:

* ``"u"`` -- one branch, one row: a U-arc on that side;
* ``"hairpin"`` -- one branch, the rows differ: the arc goes round the
  branch's pip end through the empty arc (it cannot wrap the anchor: there it
  would cross the unstable root arc, and two unstable curves never cross);
* ``"cross"`` -- the ends on two branches (period k, heteroclinic): a cubic
  leaving each end along its normal, routed round a nested circle like the
  dual-graph cartoon (best effort, as there).

Dev Notes:
    U-arcs and hairpins of one branch are drawn in ONE chart, so they nest
    exactly as the real arcs do. Unfold the two rows into a line,
    ``w = E - r`` on the outside row and ``w = -(E - r)`` on the inside row
    (``r`` the rank, ``E`` one rank past the last: the pip end is ``w = 0``),
    draw every arc as a flat half ellipse over its interval in the upper half
    ``w``-plane (the semicircles squeezed vertically, so they still do not
    cross, and only an arc near ``w = 0`` folds round the pip end), and fold
    the half plane onto the plane slit along the branch by doubling the angle
    (``|w| e^{i arg w} -> |w| e^{2 i arg w}``). The fold
    is a bijection of the open half plane onto the slit plane and every later
    step (angle and radius along the circle, the monotone compressions that
    keep a curve within its circle's room) is monotone, so arcs that do not
    cross give curves that do not cross: a U-arc is the image of a
    one-row interval, a hairpin of an interval round ``w = 0``. Within the
    unfolded line the depth of an arc grows with its span (the ellipse's
    height is :data:`FOLD_FLATTEN` times half its span): ``base + slope *
    span`` with base 0.
    A lobe of a nested zone's parent enclosing the nested orbit
    (:func:`~.plotting.cartoon_enclosures`) is a U-arc routed the long way
    round instead (:func:`~.plotting._polar_route`), outside the chart.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Iterable, Literal, Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

from .BridgeClass import _image_chain
from .StablePartition import row_of_end
from .plotting import (
    CARTOON_ARC_BASE,
    CARTOON_ARC_SLOPE,
    CARTOON_LINE_STYLE,
    CLASS_UNKNOWN_COLOR,
    HOLE_COLORS,
    HOLE_MARKERS,
    STRONG_PIP_STYLE,
    ZONE_EMPTY_STYLE,
    ZONE_INWARD_CAP,
    ZoneLayout,
    _arc_points,
    _normal_bezier,
    _path_point,
    _polar_route,
    _zone_circles,
    class_colors,
    name_mathtext,
)

if TYPE_CHECKING:
    from ..numerics.Bridge import BridgeId
    from .BridgeClass import BridgeClassTable
    from .MinimalTrellis import MinimalTrellis
    from .SymbolicDynamics import SymbolicDynamics
    from .Trellis import Trellis

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: Height of an unfolded arc over its half span: flat, so an arc folds round
#: the pip end only near it.
FOLD_FLATTEN = 0.15

#: Radial offset of one rank unit of the folded chart, in arc lengths of a rank.
FOLD_HEIGHT = 2.5

#: Ranks between a branch's last ranked crossing and its pip end (the fold).
FOLD_MARGIN = 1.0

#: How far a hairpin may run past its pip end, as a fraction of the empty arc.
FOLD_OVERRUN = 0.45

#: Samples per drawn curve.
CURVE_SAMPLES = 96

#: Line widths of a drawn bridge and of its forward iterates.
BRIDGE_LINEWIDTH = 2.0
ITERATE_LINEWIDTH = 0.9


@dataclass
class TrellisCartoonLayout:
    """
    What :func:`plot_trellis_cartoon` drew, for tests and callers.

    Attributes:
        zone: The circles and arcs (:class:`~.plotting.ZoneLayout`); its
            ``ranks`` hold each branch's ranked stable cdists in rank order.
        points: ``{point key: (x, y)}`` for every ranked crossing; the key is
            the crossing id, or ``(anchor id, branch_key)`` for an anchor
            (one point on every branch of its periodic point).
        order: ``{branch_key: [point key, ...]}``, the ranked crossings of
            each stable branch, anchor first.
        bridges: ``{BridgeId: "u" | "hairpin" | "cross"}``, the route drawn.
        paths: ``{BridgeId: (n, 2) array}``, the sampled curve of each drawn
            bridge, from its first crossing to its second.
        iterates: ``{BridgeId: [[pair, ...] per step]}``, the image pairs
            drawn for each drawn bridge (``iterates=n``).
        representatives: ``{class name: BridgeId}`` of the drawn classes.
    """

    zone: ZoneLayout
    points: dict = field(default_factory=dict)
    order: dict = field(default_factory=dict)
    bridges: dict = field(default_factory=dict)
    paths: dict = field(default_factory=dict)
    iterates: dict = field(default_factory=dict)
    representatives: dict = field(default_factory=dict)


# ── which bridges ───────────────────────────────────────────────────────────


def _class_names(
    classes: "BridgeClassTable", dynamics: Optional["SymbolicDynamics"]
) -> dict:
    """``{BridgeClass: name}``: the dynamics' letter (inert ones too), else the entry's name."""
    lettered = {} if dynamics is None else dynamics.classes
    return {
        entry.bridge_class: (
            lettered[entry.bridge_class].letter if entry.bridge_class in lettered else entry.name
        )
        for entry in classes
    }


def class_representatives(
    classes: "BridgeClassTable",
    minimal: "MinimalTrellis",
    dynamics: Optional["SymbolicDynamics"] = None,
) -> dict[str, "BridgeId"]:
    """
    One member per class: a minimal-trellis hole bridge, else image bridge, else
    the member with the smallest unstable cdist.

    Ties within a tier go to the smallest unstable cdist of the first crossing.

    Args:
        classes: The class table (a topology's ``classes``).
        minimal: The minimal trellis of the same recompute.
        dynamics: When given, a class is keyed by its letter in the dynamics
            (inert classes have one there); otherwise by the entry's
            :attr:`~.BridgeClass.BridgeClassEntry.name`.

    Returns:
        ``{class name: BridgeId}`` in table order, for every class with a
        member (a class no computed bridge spans has none).
    """
    names = _class_names(classes, dynamics)
    holes, images = set(minimal.hole_bridge_ids), set(minimal.image_bridge_ids)
    registry = minimal.trellis.registry

    def preference(bridge_id) -> tuple:
        tier = 0 if bridge_id in holes else 1 if bridge_id in images else 2
        return (tier, registry[bridge_id[0]].cdist_key("unstable"), bridge_id)

    return {
        names[entry.bridge_class]: min(entry.bridge_ids, key=preference)
        for entry in classes
        if entry.spanned
    }


def _chosen_bridges(
    trellis: "Trellis",
    classes: "BridgeClassTable",
    names: dict,
    representatives: dict,
    show: str,
    letters: Optional[Iterable[str]],
    bridge_ids: Optional[Iterable["BridgeId"]],
) -> list:
    """The bridge ids to draw (see :func:`plot_trellis_cartoon`)."""
    if bridge_ids is not None:
        chosen = list(dict.fromkeys(bridge_ids))
        outside = [b for b in chosen if any(i not in trellis for i in b)]
        if outside:
            raise ValueError(f"bridges {outside} end at crossings outside the trellis")
        return chosen
    if show not in ("all", "representatives"):
        raise ValueError(f"show must be 'all' or 'representatives', not {show!r}")
    entries = list(classes)
    if letters is not None:
        wanted = set(letters)
        unknown = wanted - set(names.values())
        if unknown:
            raise ValueError(f"no class named {sorted(unknown)}; the classes are {list(names.values())}")
        entries = [entry for entry in entries if names[entry.bridge_class] in wanted]
    if show == "representatives":
        return [representatives[names[e.bridge_class]] for e in entries if e.spanned]
    chosen = [bridge_id for entry in entries for bridge_id in entry.bridge_ids]
    if letters is None:
        chosen += [bridge.id for bridge in trellis.bridges if bridge.id is not None]
    return list(dict.fromkeys(chosen))


def _image_pairs(trellis: "Trellis", bridge_id, steps: int) -> list:
    """The image pairs of ``bridge_id`` after 1..``steps`` steps (``_image_chain`` iterated)."""
    images, current = [], [bridge_id]
    for step in range(1, steps + 1):
        following = []
        for pair in current:
            chain = _image_chain(trellis, pair)
            if chain is None:
                logger.info(
                    "trellis cartoon: the image of %s (step %d from bridge %s) leaves the "
                    "trellis; its iterates stop there", pair, step, bridge_id,
                )
                continue
            following.extend(chain)
        if not following:
            break
        images.append(following)
        current = following
    return images


# ── ranks and points ────────────────────────────────────────────────────────


def _end(trellis: "Trellis", pair, endpoint: Literal["first", "second"]) -> tuple:
    """``(point key, stable branch key, row)`` of one end of an arc."""
    crossing = pair[0 if endpoint == "first" else 1]
    branch_key = trellis.stable_key_at(crossing, trellis.unstable_key_of(pair))
    key = (crossing, branch_key) if trellis.intersection(crossing).is_anchor else crossing
    return key, branch_key, row_of_end(trellis, pair, endpoint)


def pip_orbit(trellis: "Trellis") -> list[tuple[int, int]]:
    """``(j, crossing id)`` of every registered pip-orbit point ``f^j(q0)``, ``j < k``, in the trellis."""
    return [
        (j, image)
        for pip in trellis.pips.values()
        for j in range(trellis.registry[pip].manifold_b_key[0].k_value)
        if (image := trellis.iterate(pip, j)) in trellis
    ]


def _rank_crossings(
    trellis: "Trellis", arcs: Iterable, orbit: set, holes: list, compact: bool
) -> dict:
    """``{branch_key: [point key, ...]}`` in stable-cdist order, anchor first."""
    ranked: dict = {}
    for branch in trellis.stable_branches:
        anchors = [i for i in branch.intersection_ids if trellis.intersection(i).is_anchor]
        shown = {(anchors[0], branch.key)} if anchors else set()
        if not compact:
            shown |= {i for i in branch.intersection_ids if i not in anchors}
        ranked[branch.key] = shown
    for pair in arcs:
        for endpoint in ("first", "second"):
            key, branch_key, _row = _end(trellis, pair, endpoint)
            ranked[branch_key].add(key)
    for crossing in orbit:
        ranked[trellis.intersection(crossing).manifold_b_key].add(crossing)
    for hole in holes:
        pair = _bridge_order(trellis, hole.bounding_ids)
        for endpoint in ("first", "second"):
            key, branch_key, _row = _end(trellis, pair, endpoint)
            ranked[branch_key].add(key)

    def stable_key(point) -> tuple:
        crossing = point[0] if isinstance(point, tuple) else point
        return trellis.intersection(crossing).cdist_key("stable")

    return {branch_key: sorted(keys, key=stable_key) for branch_key, keys in ranked.items()}


def _bridge_order(trellis: "Trellis", ids) -> tuple:
    """Two crossings of one unstable arc, ordered like a ``BridgeId``."""
    return tuple(sorted(ids, key=lambda i: trellis.intersection(i).cdist_key("unstable")))


class _Chart:
    """One stable branch's fold chart: rank and row to the plane (see the Dev Notes)."""

    def __init__(self, zone: ZoneLayout, branch_key, count: int) -> None:
        fp_id, self.alpha, self.sweep = zone.arcs[branch_key]
        self.cx, self.cy, self.radius = zone.circles[fp_id]
        self.end = max(count - 1, 0) + FOLD_MARGIN
        self.unit = self.radius * abs(self.sweep) / self.end
        self.overrun = FOLD_OVERRUN * self.end
        self.outward_room = CARTOON_ARC_BASE + CARTOON_ARC_SLOPE * (count - 1)
        keepout = zone.keepout.get(fp_id)
        self.inward_room = (
            self.radius - keepout if keepout is not None else ZONE_INWARD_CAP * self.radius
        )
        self.sigma = {zone.outside[branch_key]: 1.0}
        self.sigma["right" if zone.outside[branch_key] == "left" else "left"] = -1.0

    def angle(self, rank) -> NDArray[np.float64]:
        """The angle of an ordinal position; past the pip end it slows into the empty arc."""
        rank = np.asarray(rank, dtype=float)
        past = np.maximum(rank - self.end, 0.0)
        rank = np.minimum(rank, self.end) + self.overrun * (1.0 - np.exp(-past / self.overrun))
        return self.alpha + self.sweep * rank / self.end

    def place(self, rank, offset) -> NDArray[np.float64]:
        """Points at ordinal ``rank`` and radial ``offset`` (rank units, + outside)."""
        offset = FOLD_HEIGHT * self.unit * np.asarray(offset, dtype=float)
        out, room_in = self.outward_room, self.inward_room
        offset = np.where(
            offset >= 0.0,
            out * (1.0 - np.exp(-offset / out)),
            -room_in * (1.0 - np.exp(offset / room_in)),
        )
        theta, r = self.angle(rank), self.radius + offset
        return np.column_stack([self.cx + r * np.cos(theta), self.cy + r * np.sin(theta)])

    def unfolded(self, rank: float, row: str) -> float:
        """The end's position on the unfolded line (outside row positive)."""
        return self.sigma[row] * (self.end - rank)

    def arc(self, w_p: float, w_q: float) -> NDArray[np.float64]:
        """The folded half ellipse from ``w_p`` to ``w_q``."""
        centre, half = 0.5 * (w_p + w_q), 0.5 * abs(w_p - w_q)
        phi = np.linspace(0.0, np.pi, CURVE_SAMPLES)
        if w_p < w_q:
            phi = phi[::-1]
        w = centre + half * (np.cos(phi) + 1j * FOLD_FLATTEN * np.sin(phi))
        folded = np.abs(w) * np.exp(2j * np.arctan2(np.maximum(w.imag, 0.0), w.real))
        return self.place(self.end - folded.real, folded.imag)


def _sampled(path) -> NDArray[np.float64]:
    """A matplotlib path (one cubic or a polyline) as ``CURVE_SAMPLES`` points."""
    return np.array([_path_point(path, t) for t in np.linspace(0.0, 1.0, CURVE_SAMPLES)])


def _route(
    trellis: "Trellis", zone: ZoneLayout, charts: dict, ranks: dict, points: dict,
    pair, enclosed: bool,
) -> tuple[str, NDArray[np.float64]]:
    """The route kind and sampled curve of one arc, first crossing to second."""
    (p, key_p, row_p), (q, key_q, row_q) = (_end(trellis, pair, e) for e in ("first", "second"))
    if key_p == key_q:
        chart = charts[key_p]
        kind = "u" if row_p == row_q else "hairpin"
        fp_id = zone.arcs[key_p][0]
        if kind == "u" and enclosed and chart.sigma[row_p] < 0 and fp_id in zone.keepout:
            return kind, _sampled(_polar_route(zone, fp_id, points[p], points[q], around=True))
        return kind, chart.arc(chart.unfolded(ranks[p], row_p), chart.unfolded(ranks[q], row_q))

    def normal(key, branch_key, row):
        cx, cy, _r = zone.circles[zone.arcs[branch_key][0]]
        radial = np.asarray(points[key]) - (cx, cy)
        return charts[branch_key].sigma[row] * radial / np.linalg.norm(radial)

    def cap(branch_key, row):
        return None if charts[branch_key].sigma[row] > 0 else charts[branch_key].inward_room

    fp_p, fp_q = zone.arcs[key_p][0], zone.arcs[key_q][0]
    inward = charts[key_p].sigma[row_p] < 0 and charts[key_q].sigma[row_q] < 0
    if fp_p == fp_q and inward and fp_p in zone.keepout:
        path = _polar_route(zone, fp_p, points[p], points[q], around=enclosed)
    else:
        path = _normal_bezier(
            points[p], normal(p, key_p, row_p), points[q], normal(q, key_q, row_q),
            cap_p=cap(key_p, row_p), cap_q=cap(key_q, row_q),
        )
    return "cross", _sampled(path)


# ── the figure ──────────────────────────────────────────────────────────────


def plot_trellis_cartoon(
    trellis: "Trellis",
    classes: "BridgeClassTable",
    dynamics: Optional["SymbolicDynamics"] = None,
    *,
    minimal: Optional["MinimalTrellis"] = None,
    ax=None,
    show: Literal["all", "representatives"] = "representatives",
    letters: Optional[Iterable[str]] = None,
    bridge_ids: Optional[Iterable["BridgeId"]] = None,
    iterates: int = 0,
    holes: bool = False,
    compact: bool = True,
    interior_side: Optional[dict] = None,
    parents: Optional[dict] = None,
    arc_order: Optional[dict] = None,
    enclosures: Optional[dict] = None,
    label_classes: bool = True,
) -> TrellisCartoonLayout:
    """
    Draw a trellis in standard form: ordinal crossings on circular stable arcs,
    smooth nested bridges.

    Every class is drawn in its colour (:func:`~.plotting.class_colors`,
    unrefined; one grey without ``dynamics``), inert classes dashed, each
    drawn representative with an arrow pointing anchor outward (the class's
    orientation) and its letter. The anchors are filled dots labelled with
    their fixed point, the pip-orbit points open circles.

    Args:
        trellis: The trellis (a topology's: cut at its pips).
        classes: Its class table (a topology's ``classes``).
        dynamics: The symbolic dynamics: letters (inert classes included) and
            colours.
        minimal: The minimal trellis, which picks the representatives
            (:func:`class_representatives`); required for
            ``show="representatives"`` unless ``bridge_ids`` is given.
        ax: Matplotlib axes (the current axes by default).
        show: ``"representatives"`` (one bridge per class) or ``"all"``
            (every member; without ``letters`` every bridge of the trellis,
            those in no class of ``classes`` grey).
        letters: Draw only the classes of these names.
        bridge_ids: Draw exactly these bridges (overrides ``show`` and
            ``letters``).
        iterates: Also draw, thin and dashed, the image pairs of every drawn
            bridge up to this many steps (:func:`~.BridgeClass._image_chain`
            iterated; a chain leaving the trellis stops, with one INFO log).
        holes: Draw the trellis's holes, each just off its crossing inside its
            bridge's lobe, labelled with its iterate.
        compact: Rank only the crossings the figure shows (drawn ends, pip
            orbit, drawn holes); False ranks every crossing of the trellis.
        interior_side, parents, arc_order: The zones, as for the circular
            dual-graph cartoon (:func:`tanglepack.loom.cartoon_zones`).
        enclosures: ``{bridge_id: frozenset(id(nested fixed point))}``
            (:func:`~.plotting.cartoon_enclosures`): the lobes routed the long
            way round a nested circle.
        label_classes: Write each drawn representative's class name.

    Returns:
        The :class:`TrellisCartoonLayout`.

    Raises:
        ValueError: On an unknown ``show`` or class name, a bridge ending
            outside the trellis, or ``show="representatives"`` without
            ``minimal``.
    """
    import matplotlib.pyplot as plt

    target = ax if ax is not None else plt.gca()
    names = _class_names(classes, dynamics)
    entry_of = {member.bridge_id: entry for entry in classes for member in entry.members}
    if minimal is None and bridge_ids is None and show == "representatives":
        raise ValueError("show='representatives' needs the minimal trellis (minimal=...)")
    representatives = {} if minimal is None else class_representatives(classes, minimal, dynamics)
    drawn = _chosen_bridges(trellis, classes, names, representatives, show, letters, bridge_ids)
    images = {bridge_id: _image_pairs(trellis, bridge_id, iterates) for bridge_id in drawn} if iterates else {}
    shown_holes = list(trellis.holes) if holes else []
    arcs = drawn + [pair for steps in images.values() for pairs in steps for pair in pairs]

    orbit = {crossing for _j, crossing in pip_orbit(trellis)}
    zone = ZoneLayout()
    order = _rank_crossings(trellis, arcs, orbit, shown_holes, compact)
    for branch_key, keys in order.items():
        zone.ranks[branch_key] = [
            float(trellis.intersection(k[0] if isinstance(k, tuple) else k).stable_cdist)
            for k in keys
        ]
    _zone_circles(zone, list(order), interior_side, parents, arc_order)
    charts = {key: _Chart(zone, key, len(keys)) for key, keys in order.items()}
    ranks = {point: rank for keys in order.values() for rank, point in enumerate(keys)}
    points = {
        point: tuple(charts[branch_key].place(rank, 0.0)[0])
        for branch_key, keys in order.items()
        for rank, point in enumerate(keys)
    }
    shown = set(drawn)
    layout = TrellisCartoonLayout(
        zone=zone, points=points, order=order, iterates=images,
        representatives={name: b for name, b in representatives.items() if b in shown},
    )
    extent = _draw_frame(target, zone, order, points, orbit, dynamics)

    palette = class_colors(dynamics, refined=False) if dynamics is not None else {}
    arrowed = set(layout.representatives.values())
    enclosed_ids = set(enclosures or {})

    def style_of(bridge_id) -> tuple[str, bool]:
        """The class colour, and whether the class is inert (dashed)."""
        entry = entry_of.get(bridge_id)
        if entry is None:
            return CLASS_UNKNOWN_COLOR, False
        known = dynamics.classes.get(entry.bridge_class) if dynamics is not None else None
        inert = known.kind != "active" if known is not None else bool(entry.inert)
        return palette.get(names[entry.bridge_class], CLASS_UNKNOWN_COLOR), inert

    for bridge_id in drawn:
        color, _inert = style_of(bridge_id)
        for pairs in images.get(bridge_id, []):
            for pair in pairs:
                _kind, curve = _route(trellis, zone, charts, ranks, points, pair, pair in enclosed_ids)
                target.plot(curve[:, 0], curve[:, 1], color=color, linewidth=ITERATE_LINEWIDTH,
                            linestyle="--", zorder=4)
                extent.extend(curve)
    for bridge_id in drawn:
        color, dashed = style_of(bridge_id)
        kind, curve = _route(trellis, zone, charts, ranks, points, bridge_id, bridge_id in enclosed_ids)
        layout.bridges[bridge_id], layout.paths[bridge_id] = kind, curve
        target.plot(curve[:, 0], curve[:, 1], color=color, linewidth=BRIDGE_LINEWIDTH,
                    linestyle="--" if dashed else "-", solid_capstyle="round", zorder=6)
        extent.extend(curve)
        if bridge_id in arrowed:
            entry = entry_of[bridge_id]
            _arrow(target, curve, entry.member(bridge_id).direction, color)
            if label_classes:
                mid = curve[len(curve) // 2]
                target.annotate(
                    name_mathtext(names[entry.bridge_class]), tuple(mid), textcoords="offset points",
                    xytext=(5, 5), fontsize=11, color=color, zorder=9,
                )
    for hole in shown_holes:
        pair = _bridge_order(trellis, hole.bounding_ids)
        extent.append(_draw_hole(target, trellis, zone, charts, ranks, points, hole, pair,
                                 pair in enclosed_ids))

    xy = np.array(extent, dtype=float).reshape(-1, 2)
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    pad = 0.6 + 0.04 * float(np.max(hi - lo))
    target.set_xlim(lo[0] - pad, hi[0] + pad)
    target.set_ylim(lo[1] - pad, hi[1] + pad)
    target.set_aspect("equal", adjustable="box")
    target.set_axis_off()
    return layout


def _draw_frame(target, zone: ZoneLayout, order: dict, points: dict, orbit: set,
                dynamics: Optional["SymbolicDynamics"]) -> list:
    """The circles, stable arcs, branch codes, anchors and pip-orbit points; returns their extent."""
    extent: list = []
    for cx, cy, r in zone.circles.values():
        ring = _arc_points(cx, cy, r, 0.0, 2 * np.pi, n=256)
        target.plot(ring[:, 0], ring[:, 1], **ZONE_EMPTY_STYLE)
        extent.extend([(cx - r, cy - r), (cx + r, cy + r)])
    naming = getattr(dynamics, "naming", None)
    for branch_key, (fp_id, alpha, sweep) in zone.arcs.items():
        cx, cy, r = zone.circles[fp_id]
        stable = _arc_points(cx, cy, r, alpha, alpha + sweep)
        target.plot(stable[:, 0], stable[:, 1], solid_capstyle="round", **CARTOON_LINE_STYLE)
        code = naming.branch_code(branch_key) if naming is not None else f"{branch_key[2]}.{branch_key[3]}"
        at = zone.point(branch_key, alpha - 0.12 * sweep)
        target.annotate(
            code, at, ha="center", va="center", fontsize=8, color="#52514e", zorder=9,
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "edgecolor": "#b4b2a9"},
        )
        anchor = order[branch_key][0]
        if isinstance(anchor, tuple):
            x, y = points[anchor]
            target.plot(x, y, marker="o", markersize=7, color="black", zorder=10)
            fixed_point = branch_key[0]
            target.annotate(
                getattr(fixed_point, "label", None) or "P", (x, y), textcoords="offset points",
                xytext=(-6, -10), ha="right", va="top", fontsize=10, zorder=10,
            )
        for point in order[branch_key]:
            if point in orbit:
                target.scatter(*points[point], s=60, facecolors="none", edgecolors=STRONG_PIP_STYLE["color"],
                               linewidths=1.5, zorder=STRONG_PIP_STYLE["zorder"])
    return extent


def _arrow(target, curve: NDArray[np.float64], direction: int, color: str) -> None:
    """A barbed arrow at the curve's middle, along it (``direction >= 0``) or against it."""
    i = len(curve) // 2
    step = curve[i + 1] - curve[i - 1]
    if direction < 0:
        step = -step
    target.annotate(
        "", xy=tuple(curve[i] + 0.5 * step), xytext=tuple(curve[i] - 0.5 * step),
        arrowprops={"arrowstyle": "-|>", "color": color, "lw": 1.5, "mutation_scale": 16},
        zorder=7,
    )


def _draw_hole(target, trellis: "Trellis", zone: ZoneLayout, charts: dict, ranks: dict,
               points: dict, hole, pair: tuple, enclosed: bool) -> tuple:
    """One hole inside its bridge's lobe (``pair``, in bridge order), near its own crossing,
    labelled with its iterate."""
    _kind, curve = _route(trellis, zone, charts, ranks, points, pair, enclosed)
    near_first = pair[0] == hole.near_intersection_id
    on_bridge = curve[int((0.3 if near_first else 0.7) * (len(curve) - 1))]
    _near, branch_key, _row = _end(trellis, pair, "first" if near_first else "second")
    cx, cy, radius = zone.circles[zone.arcs[branch_key][0]]
    radial = on_bridge - (cx, cy)
    on_stable = np.array([cx, cy]) + radius * radial / np.linalg.norm(radial)
    x, y = 0.5 * (on_bridge + on_stable)
    origins = sorted({h.origin for h in trellis.holes if h.origin is not None})
    index = origins.index(hole.origin) if hole.origin in origins else 0
    color = HOLE_COLORS[index % len(HOLE_COLORS)]
    target.scatter(x, y, marker=HOLE_MARKERS[index % len(HOLE_MARKERS)], color=color, s=40, zorder=14)
    if hole.iterate is not None:
        target.annotate(str(hole.iterate), (x, y), textcoords="offset points", xytext=(4, 4),
                        fontsize=8, color=color, zorder=15)
    return (x, y)

"""
The inner period-3 Hénon tangle as ``k`` increases: one figure per ``k``, headed by its mean lobe area.

Sweeps the nonlinearity ``k`` of the Hénon map ``(x, y) -> (y - k + x**2, -x)``
(``b = 1``, area-preserving), builds the inner period-3 saddle orbit's tangle at
each value -- the orbit of ``henon_cases.build_period3`` at ``k = 2``, continued
in ``k`` from the previous value's orbit so the sweep follows ONE orbit family
-- computes every lobe's area, and writes the frames in ``k`` order to
``figures/henon_period3_k_sweep/`` (``frame_000_k1.300.png``,
``frame_001_k1.400.png``, ...) so they can be flipped through. A ``lobes.csv``
and a ``summary.png`` (mean lobe area against ``k``) close the sweep.

No topology is computed: no strong pip, no holes, no partition, no symbolic
dynamics. Each frame is the manifolds, the crossings and the bridges alone.

A LOBE is a bridge (the unstable arc between two consecutive crossings) closed
by the stable segment between its two ends (the lobe of ``DualWalk._lobe``,
the disk D13 argues over): its area is the shoelace area of that ring, read
off the two polylines. A bridge whose ends lie on two different stable
branches bounds no lobe and is skipped. Images of a lobe have the same area
(area preservation), so the title also gives the mean over the DISTINCT lobe
orbits -- the lobes with no registered preimage lobe, the heads of the
iterate chains. A lobe whose ring reaches beyond the escape radius is an
ESCAPE LOBE: a tongue of the unstable manifold leaving the trapped region
(at ``k = 4`` one reaches ``x ~ 900``), real, but whether the trellis holds it
depends on the growth, so it is counted, drawn grey, and left out of every
statistic and of the frame box. A frame whose median lobe area is within a
decade of ``area_cutoff`` is flagged in its title: the lobes there are at the
resolution floor (the low-``k`` end at the default 1e-7; pass a smaller
``--area-cutoff`` to resolve them).

Growth policy (the one thing that differs from ``build_period3``): the
unstable eigenvalue runs from ~1.1 near the period-tripling bifurcation
(``k ~ 1.25``, where the orbit is born from the elliptic fixed point) to ~56 at
``k = 6``, so no fixed step count serves every ``k`` -- 13 steps blow the point
budget at ``k >= 3`` and leave ``k = 1.5`` with no crossing at all. The step
counts (``--unstable-steps`` / ``--stable-steps``, default the ``k = 2``
recipe's 13 / 9) are therefore read AT ``k = 2`` and scaled at every other
``k`` to the same total stretch: ``steps(k) = round(steps * ln lambda_u(2) /
ln lambda_u(k))`` (35 / 24 at ``k = 1.5``, 7 / 5 at ``k = 3``; ``--no-scale``
uses the counts verbatim). Each manifold is grown one map step at a time
towards that target and stops early when its arm ESCAPES the trapped region
(``max |coordinate|`` beyond ``radius``, by default twice the outer period-1
saddle's coordinate ``1 + sqrt(1 + k)``; the step that crosses it is cheap,
the next one is the explosion) or when the next step is predicted over the
point budget. The title of every frame says which stop it hit. Growing until
escape alone (no target) over-grows the weakly hyperbolic end: at ``k = 1.9``
it runs 20 steps to ~800 lobes, most of them under-resolved at 1e-7.

Run with ``PYTHONPATH=src python3 scripts/henon_period3_k_sweep.py``
(``--k-min``, ``--k-max``, ``--k-step``, ``--area-cutoff``, ``--radius``,
``--unstable-steps``, ``--stable-steps``, ``--no-scale``, ``--no-fill``,
``--axes common``, ``--out``).

Dev Notes:
    * Continuation: the orbit at each ``k`` is solved from the previous
      ``k``'s orbit (``[[0, 1], [-1, 0], [-1, 1]]`` seeds ``k = 2``; the sweep
      runs from the seed outward in both directions so a ``k_min`` below 2
      still continues from the one known orbit), and the eigenvector
      orientation hints are the previous ``k``'s ORIENTED eigenvectors, so
      the branch grown on each orbit point never flips between frames (the
      seed's own hints are ``build_period3``'s). Below
      the period tripling (``k < 1.25``) the continuation lands on the mirror
      image of the orbit through the elliptic point; the default ``k_min``
      stays above it.
    * The axes frame each tangle on its own by default (``--axes common``
      for the union of every frame's box: the tangle grows from ~0.5 wide at
      ``k = 1.3`` to ~5 wide at ``k = 4``, which shrinks the low-``k`` frames
      to a dot). A frame's box holds its orbit points, crossings and inside
      lobes -- never a whole manifold, since a stable branch no crossing trims
      keeps its escaping arm. Both manifolds are
      clipped to a few times the frame's box before they are stored, so a
      frame never holds an explosion step's points.
    * A ``k`` whose solve fails, whose growth exceeds the point budget or
      whose bridges cannot be cut is logged and skipped; the next ``k``
      continues from the last orbit that solved.
    * The low end is a near-tangency, not an under-resolution a finer cutoff
      cures: the splitting of the period-3 manifolds is exponentially small
      near the period tripling, and at ``k <= 1.7`` a cutoff of 1e-10 yields
      MORE and SMALLER lobes than 1e-7 (``k = 1.6``: 141 lobes at 5.6e-8
      become 1551 at 7e-11, nearly every one its own "orbit": polygon noise
      between two curves that nearly coincide), while ``k = 1.8`` converges
      (30 lobes / 14 orbits at 1e-7, 23 / 3 at 1e-10). Trust the frames from
      ``k ~ 1.8`` up at the default cutoff; the ones flagged "at the
      resolution floor" show where the manifolds run, not their lobes.
    * Measured on 2026-10-14 (default sweep, 28 frames, ~15 s): the mean lobe
      area rises monotonically from 1.8e-5 at ``k = 1.9`` through 1.1e-4
      (``k = 2``, 24 lobes of two orbits, all within 1% -- area preservation
      along each chain), 1e-2 at 2.5, 6.6e-2 at 3.0 and 0.25 at 3.7; from
      3.8 the stable branches reach the escaping unstable arm and two outer
      lobes of area ~7 enter the trellis, so the mean jumps to ~1.6 while the
      median stays at 0.36 (the title shows both).
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
from numpy.typing import NDArray

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

import tanglepack  # noqa: E402
from tanglepack.examples.henon import henon_jacobian, henon_map, henon_map_inverse  # noqa: E402
from tanglepack.numerics.geometry import signed_polygon_area  # noqa: E402
from tanglepack.topology.TopologyResults import Arc  # noqa: E402

logger = logging.getLogger("henon_period3_k_sweep")

#: The period-3 orbit of the nested Hénon map at ``k = 2`` (``henon_cases``).
SEED_K = 2.0
SEED_GUESS = [[0, 1], [-1, 0], [-1, 1]]
SEED_ORIENTATION = {"unstable": np.array([0.0, -1.0]), "stable": np.array([-1.0, -1.0])}

OUT_DIR = Path(__file__).resolve().parent.parent / "figures" / "henon_period3_k_sweep"

STABLE_COLOR = "#c0392b"
UNSTABLE_COLOR = "#2c5aa0"
LOBE_COLOR = "#f1c40f"


# ── lobes ──────────────────────────────────────────────────────────────────


@dataclass(frozen=True)
class Lobe:
    """
    One lobe of a trellis: a bridge closed by the stable segment between its ends.

    Attributes:
        bridge_id: The bridge's ``(first_id, second_id)``.
        area: The absolute shoelace area of the ring.
        ring: The ``(N, 2)`` closed polygon (unstable arc, then the stable arc
            back to the start).
        head: True when no registered preimage lobe exists in the trellis: the
            lobe is the head of its iterate chain, one per distinct lobe orbit.
        reach: ``max |coordinate|`` over the ring. A lobe reaching beyond the
            escape radius is a tongue leaving the trapped region (an escape
            lobe): real, but not a lobe of the tangle's interior, so the
            statistics leave it out.
    """

    bridge_id: tuple[int, int]
    area: float
    ring: NDArray[np.float64]
    head: bool
    reach: float


def lobe_areas(trellis, bridges=None) -> list[Lobe]:
    """
    Every lobe of a trellis with its area.

    Args:
        trellis: A :class:`~tanglepack.topology.Trellis.Trellis` of one fixed
            point, bridges cut.
        bridges: The bridges to measure; default every bridge of the trellis.
            The chain-head test reads every bridge of the trellis regardless.

    Returns:
        One :class:`Lobe` per bridge whose two ends lie on ONE stable branch,
        in the given order. A bridge joining two stable branches bounds no
        lobe and is skipped; a partial piece has no id and is skipped.
    """
    bridge_ids = {bridge.id for bridge in trellis.bridges if bridge.id is not None}
    lobes: list[Lobe] = []
    for bridge in (trellis.bridges if bridges is None else bridges):
        if bridge.id is None:
            continue
        first, second = bridge.id
        unstable_key = trellis.unstable_key_of(bridge.id)
        stable_keys = {trellis.stable_key_at(end, unstable_key) for end in bridge.id}
        if len(stable_keys) != 1:
            continue
        stable_key = stable_keys.pop()
        unstable = Arc(
            kind="unstable",
            lo_id=first,
            hi_id=second,
            branch_key=unstable_key,
            bridge_id=bridge.id,
        ).polyline(trellis)
        s_first = trellis.intersection(first).stable_cdist
        s_second = trellis.intersection(second).stable_cdist
        lo_id, hi_id = (first, second) if s_first <= s_second else (second, first)
        # The stable arc is walked back from ``second`` to ``first``.
        stable = Arc(
            kind="stable",
            lo_id=lo_id,
            hi_id=hi_id,
            branch_key=stable_key,
            reverse=(lo_id == first),
        ).polyline(trellis)
        ring = np.vstack([unstable, stable])
        preimage = (trellis.iterate(first, -1), trellis.iterate(second, -1))
        head = not (None not in preimage and tuple(preimage) in bridge_ids)
        lobes.append(
            Lobe(bridge.id, abs(signed_polygon_area(ring)), ring, head, float(np.abs(ring).max()))
        )
    return lobes


# ── one frame ──────────────────────────────────────────────────────────────


@dataclass
class Frame:
    """
    What one ``k`` contributes to the sweep: the curves to draw and the numbers.

    Attributes:
        k: The nonlinearity.
        orbit: The ``(3, 2)`` orbit points.
        lambda_u: The unstable eigenvalue of the 3-cycle.
        unstable_steps: Map steps the unstable manifold was grown.
        stable_steps: Map steps the stable manifold was grown.
        stable: The trimmed stable polylines, one per branch, clipped like
            the unstable ones.
        unstable: The unstable polylines, one per branch, clipped to a few
            times the tangle's box (NaN rows break the line where points were
            dropped).
        crossings: The ``(M, 2)`` registered crossings.
        lobes: Every lobe with its area and ring, escape lobes included.
        radius: The escape radius: a lobe reaching beyond it is an escape lobe.
        area_cutoff: The refinement cutoff the manifolds were grown at.
        escaped: Per stability, whether growth stopped on escape (else it
            reached its target, or the step budget stopped it short).
        target: Per stability, the step target at this ``k``.
    """

    k: float
    orbit: NDArray[np.float64]
    lambda_u: float
    unstable_steps: int
    stable_steps: int
    stable: list[NDArray[np.float64]]
    unstable: list[NDArray[np.float64]]
    crossings: NDArray[np.float64]
    lobes: list[Lobe]
    radius: float
    area_cutoff: float
    escaped: dict[str, bool] = field(default_factory=dict)
    target: dict[str, int] = field(default_factory=dict)

    @property
    def inside(self) -> list[Lobe]:
        """The lobes within the escape radius: the ones every statistic reads."""
        return [lobe for lobe in self.lobes if lobe.reach <= self.radius]

    @property
    def escaping(self) -> list[Lobe]:
        """The escape lobes: tongues reaching beyond the escape radius."""
        return [lobe for lobe in self.lobes if lobe.reach > self.radius]

    @property
    def areas(self) -> NDArray[np.float64]:
        return np.array([lobe.area for lobe in self.inside], dtype=float)

    @property
    def mean_area(self) -> float:
        areas = self.areas
        return float(areas.mean()) if len(areas) else float("nan")

    @property
    def head_areas(self) -> NDArray[np.float64]:
        return np.array([lobe.area for lobe in self.inside if lobe.head], dtype=float)

    @property
    def mean_head_area(self) -> float:
        heads = self.head_areas
        return float(heads.mean()) if len(heads) else float("nan")

    @property
    def mean_escaping_area(self) -> float:
        areas = [lobe.area for lobe in self.escaping]
        return float(np.mean(areas)) if areas else float("nan")

    @property
    def at_resolution_floor(self) -> bool:
        """True when the median lobe area is within a decade of the refinement cutoff."""
        areas = self.areas
        return bool(len(areas)) and float(np.median(areas)) < 10.0 * self.area_cutoff

    def box(self) -> tuple[float, float, float, float]:
        """``(xmin, xmax, ymin, ymax)`` of the tangle: orbit, crossings and inside lobes."""
        return tangle_box(self.orbit, self.crossings, [lobe.ring for lobe in self.inside])


def tangle_box(
    orbit: NDArray[np.float64], crossings: NDArray[np.float64], rings: list[NDArray[np.float64]]
) -> tuple[float, float, float, float]:
    """
    The box of a tangle: its orbit, crossings and lobes -- never a whole manifold.

    A stable branch no crossing trims keeps its escaping arm, and the unstable
    arms always escape, so the manifolds themselves never frame a plot. With no
    crossing the orbit's own box, widened by its diameter, stands in.
    """
    points = np.vstack([orbit, crossings.reshape(-1, 2), *rings])
    xmin, xmax = float(points[:, 0].min()), float(points[:, 0].max())
    ymin, ymax = float(points[:, 1].min()), float(points[:, 1].max())
    if len(crossings) == 0 and not rings:
        margin = max(xmax - xmin, ymax - ymin, 1e-3) * 0.5
        return xmin - margin, xmax + margin, ymin - margin, ymax + margin
    return xmin, xmax, ymin, ymax


def grow_until_escape(
    session,
    fixed_point,
    stability: str,
    *,
    radius: float,
    max_steps: int,
    step_budget: int,
) -> tuple[int, bool]:
    """
    Grow one manifold a step at a time until its arm escapes or the cap is hit.

    Args:
        session: The :class:`~tanglepack.TangleSession`.
        fixed_point: The orbit.
        stability: ``"unstable"`` or ``"stable"``.
        radius: ``max |coordinate|`` that ends the growth: the step crossing
            it is the last.
        max_steps: The cap on the number of steps.
        step_budget: Points one step may be predicted to create (the last
            step's count times its growth factor) before growth stops short.

    Returns:
        ``(steps grown, escaped)``.
    """
    for step in range(1, max_steps + 1):
        stats = session.growth_stats(fixed_point, stability)
        if stats is not None and stats.growth_factor is not None:
            predicted = stats.points_created * stats.growth_factor
            if predicted > step_budget:
                logger.info(
                    "%s manifold: step %d predicted at %.3g points, over the step budget; stopping",
                    stability, step, predicted,
                )
                return step - 1, False
        session.grow_n_times(fixed_point, stability, num_iterations=1)
        stats = session.growth_stats(fixed_point, stability)
        if stats.max_abs_coord > radius:
            return step, True
    return max_steps, False


def clip_polyline(points: NDArray[np.float64], box: tuple[float, float, float, float]) -> NDArray[np.float64]:
    """Drop the points outside ``box``, breaking the line with NaNs where they were."""
    xmin, xmax, ymin, ymax = box
    inside = (points[:, 0] >= xmin) & (points[:, 0] <= xmax) & (points[:, 1] >= ymin) & (points[:, 1] <= ymax)
    if inside.all():
        return points
    out = points.astype(float).copy()
    out[~inside] = np.nan
    # Collapse runs of NaN rows to one so the array stays small.
    keep = inside | np.concatenate([[True], inside[:-1]])
    return out[keep]


def build_frame(
    k: float,
    guess,
    orientation: dict[str, NDArray[np.float64]],
    *,
    area_cutoff: float,
    radius: Optional[float],
    unstable_steps: int,
    stable_steps: int,
    scale: bool,
    lambda_ref: float,
    step_budget: int,
) -> tuple[Frame, list, dict[str, NDArray[np.float64]]]:
    """
    Build the period-3 tangle at one ``k`` and read its lobes.

    Args:
        k: The nonlinearity.
        guess: The ``(3, 2)`` orbit guess (the previous frame's orbit).
        orientation: Eigenvector orientation hints (the previous frame's
            oriented eigenvectors).
        area_cutoff: The manifold machine's refinement cutoff.
        radius: The escape radius; None = ``2 * (1 + sqrt(1 + k))``.
        unstable_steps: Unstable step target at ``lambda_ref``.
        stable_steps: Stable step target at ``lambda_ref``.
        scale: Scale the targets by ``ln(lambda_ref) / ln(lambda_u)`` (equal
            total stretch at every ``k``); False uses them verbatim.
        lambda_ref: The unstable eigenvalue the targets are read at (``k = 2``).
        step_budget: Points a step may be predicted to create.

    Returns:
        ``(frame, orbit, oriented eigenvectors)``, the last two seeding the
        next ``k``.
    """
    session = tanglepack.TangleSession(
        henon_map(k), henon_map_inverse(k), henon_jacobian(k), area_cutoff=area_cutoff
    )
    fp = session.construct_fixed_point(guess)
    session.orient_eigenvectors(fp, orientation)
    session.initialize_both_manifolds(fp)
    lambda_u = float(fp.unstable_eigenvalues[0])
    escape_radius = radius if radius is not None else 2.0 * (1.0 + math.sqrt(1.0 + k))

    factor = math.log(lambda_ref) / math.log(lambda_u) if scale else 1.0
    steps: dict[str, int] = {}
    escaped: dict[str, bool] = {}
    targets: dict[str, int] = {}
    for stability, count in (("unstable", unstable_steps), ("stable", stable_steps)):
        targets[stability] = max(1, int(round(count * factor)))
        steps[stability], escaped[stability] = grow_until_escape(
            session, fp, stability,
            radius=escape_radius, max_steps=targets[stability], step_budget=step_budget,
        )

    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    trellis = session.trellis(fp)
    lobes = lobe_areas(trellis)

    orbit = np.array(fp.coordinates, dtype=float)
    crossings = np.array(
        [trellis.intersection(i).get_point() for i in trellis.intersection_ids if not trellis.intersection(i).is_anchor],
        dtype=float,
    ).reshape(-1, 2)
    # Clip both manifolds' escaping arms to a few times the tangle's own box
    # before storing them: a frame never holds an explosion step's points.
    xmin, xmax, ymin, ymax = tangle_box(
        orbit, crossings, [lobe.ring for lobe in lobes if lobe.reach <= escape_radius]
    )
    cx, cy = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
    half = max(xmax - xmin, ymax - ymin, 1e-3) * 2.0
    clip_box = (cx - half, cx + half, cy - half, cy + half)
    stable = [
        clip_polyline(manifold.get_point_array(), clip_box)
        for (owner, stability, _o, _b), manifold in session.manifolds.items()
        if owner is fp and stability == "stable"
    ]
    unstable = [
        clip_polyline(manifold.get_point_array(), clip_box)
        for (owner, stability, _o, _b), manifold in session.manifolds.items()
        if owner is fp and stability == "unstable"
    ]

    frame = Frame(
        k=k,
        orbit=orbit,
        lambda_u=lambda_u,
        unstable_steps=steps["unstable"],
        stable_steps=steps["stable"],
        stable=stable,
        unstable=unstable,
        crossings=crossings,
        lobes=lobes,
        radius=escape_radius,
        area_cutoff=area_cutoff,
        escaped=escaped,
        target=targets,
    )
    oriented = {
        "unstable": np.array(fp.unstable_eigenvectors[0], dtype=float),
        "stable": np.array(fp.stable_eigenvectors[0], dtype=float),
    }
    return frame, [c.tolist() for c in fp.coordinates], oriented


# ── drawing ────────────────────────────────────────────────────────────────


def draw_frame(frame: Frame, limits: tuple[float, float, float, float], path: Path, *, fill: bool) -> None:
    """Draw one frame to ``path``: manifolds, crossings, lobes, the numbers in the title."""
    fig, ax = plt.subplots(figsize=(9, 9))
    if fill:
        for lobe in frame.escaping:
            ax.fill(lobe.ring[:, 0], lobe.ring[:, 1], color="grey", alpha=0.12, lw=0)
        for lobe in frame.inside:
            ax.fill(lobe.ring[:, 0], lobe.ring[:, 1], color=LOBE_COLOR, alpha=0.25, lw=0)
    for poly in frame.unstable:
        ax.plot(poly[:, 0], poly[:, 1], color=UNSTABLE_COLOR, lw=0.8)
    for poly in frame.stable:
        ax.plot(poly[:, 0], poly[:, 1], color=STABLE_COLOR, lw=0.8)
    if len(frame.crossings):
        ax.scatter(frame.crossings[:, 0], frame.crossings[:, 1], s=6, color="k", zorder=4)
    ax.scatter(frame.orbit[:, 0], frame.orbit[:, 1], s=40, color="k", marker="o", zorder=5)
    for i, (x, y) in enumerate(frame.orbit):
        ax.annotate(f"p{i}", (x, y), textcoords="offset points", xytext=(5, 5), fontsize=9)
    ax.set_xlim(limits[0], limits[1])
    ax.set_ylim(limits[2], limits[3])
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    inside, escaping = frame.inside, frame.escaping
    n_heads = int(sum(lobe.head for lobe in inside))
    lines = [f"Hénon k = {frame.k:.3f}, b = 1 — inner period-3 tangle"]
    if inside:
        areas = frame.areas
        lines.append(
            f"mean lobe area = {frame.mean_area:.4g}   ({len(inside)} lobes; "
            f"median {np.median(areas):.3g}, min {areas.min():.3g}, max {areas.max():.3g})"
        )
        lines.append(f"{n_heads} distinct lobe orbits, mean area {frame.mean_head_area:.4g}")
    else:
        lines.append("mean lobe area = n/a   (no lobe inside the trapped region)")
    if escaping:
        lines.append(
            f"{len(escaping)} escape lobe(s) reaching beyond |coord| = {frame.radius:.3g} "
            f"excluded (mean area {frame.mean_escaping_area:.3g}; grey)"
        )
    if frame.at_resolution_floor:
        lines.append(
            f"⚠ lobe areas within a decade of area_cutoff = {frame.area_cutoff:.0e}: at the resolution floor"
        )
    stop = {
        s: ("escape" if frame.escaped.get(s)
            else "target" if getattr(frame, s + "_steps") == frame.target.get(s)
            else "budget")
        for s in ("unstable", "stable")
    }
    lines.append(
        f"λ_u = {frame.lambda_u:.3g};  grown {frame.unstable_steps} unstable ({stop['unstable']}) / "
        f"{frame.stable_steps} stable ({stop['stable']}) steps;  {len(frame.crossings)} crossings"
    )
    ax.set_title("\n".join(lines), fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def draw_summary(frames: list[Frame], path: Path) -> None:
    """Mean lobe area (over lobes and over distinct lobe orbits) against ``k``."""
    ks = [f.k for f in frames]
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    ax.plot(ks, [f.mean_area for f in frames], "o-", color=UNSTABLE_COLOR, label="mean over lobes")
    ax.plot(ks, [f.mean_head_area for f in frames], "s--", color=STABLE_COLOR, label="mean over distinct lobe orbits")
    ax.set_ylabel("mean lobe area (inside the trapped region)")
    ax.set_yscale("log")
    ax.legend()
    ax.set_title("Inner period-3 Hénon tangle: lobe area against k", loc="left")
    ax2.plot(ks, [len(f.inside) for f in frames], "o-", color="k", label="lobes")
    ax2.plot(ks, [sum(l.head for l in f.inside) for f in frames], "s--", color="grey", label="distinct lobe orbits")
    ax2.plot(ks, [len(f.escaping) for f in frames], "^:", color=LOBE_COLOR, label="escape lobes (excluded)")
    ax2.set_ylabel("count")
    ax2.set_xlabel("k")
    ax2.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def write_csv(frames: list[Frame], path: Path) -> None:
    """One row per ``k``: the lobe statistics and the growth used."""
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["k", "lambda_u", "unstable_steps", "stable_steps", "crossings", "lobes",
             "distinct_lobe_orbits", "mean_lobe_area", "mean_orbit_lobe_area",
             "min_lobe_area", "max_lobe_area", "total_lobe_area",
             "escape_lobes", "mean_escape_lobe_area", "escape_radius", "at_resolution_floor"]
        )
        for f in frames:
            areas = f.areas
            writer.writerow(
                [f"{f.k:.4f}", f"{f.lambda_u:.6g}", f.unstable_steps, f.stable_steps,
                 len(f.crossings), len(f.inside), int(sum(l.head for l in f.inside)),
                 f"{f.mean_area:.6g}", f"{f.mean_head_area:.6g}",
                 f"{areas.min():.6g}" if len(areas) else "", f"{areas.max():.6g}" if len(areas) else "",
                 f"{areas.sum():.6g}" if len(areas) else "",
                 len(f.escaping), f"{f.mean_escaping_area:.6g}" if f.escaping else "",
                 f"{f.radius:.6g}", int(f.at_resolution_floor)]
            )


# ── the sweep ──────────────────────────────────────────────────────────────


def sweep_order(ks: list[float]) -> list[list[float]]:
    """Split ``ks`` into the runs continued from the seed: downward from it, then upward."""
    below = sorted([k for k in ks if k < SEED_K], reverse=True)
    above = sorted([k for k in ks if k >= SEED_K])
    return [run for run in (below, above) if run]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--k-min", type=float, default=1.3, help="first k (default 1.3, just above the period tripling)")
    parser.add_argument("--k-max", type=float, default=4.0, help="last k (default 4.0)")
    parser.add_argument("--k-step", type=float, default=0.1, help="k increment (default 0.1)")
    parser.add_argument("--area-cutoff", type=float, default=1e-7, help="refinement cutoff (default 1e-7)")
    parser.add_argument("--radius", type=float, default=None,
                        help="escape radius, max |coordinate| (default 2 * (1 + sqrt(1 + k)))")
    parser.add_argument("--unstable-steps", type=int, default=13,
                        help="unstable growth steps at k = 2, scaled to equal stretch elsewhere (default 13)")
    parser.add_argument("--stable-steps", type=int, default=9,
                        help="stable growth steps at k = 2, scaled to equal stretch elsewhere (default 9)")
    parser.add_argument("--no-scale", action="store_true", help="use the step counts verbatim at every k")
    parser.add_argument("--step-budget", type=int, default=500_000,
                        help="points one growth step may be predicted to create before growth stops (default 500000)")
    parser.add_argument("--no-fill", action="store_true", help="do not shade the lobes")
    parser.add_argument("--axes", choices=("frame", "common"), default="frame",
                        help="axes limits: each frame's own tangle box (default) or the union of every frame's")
    parser.add_argument("--pad", type=float, default=0.08, help="axes padding as a fraction of the box")
    parser.add_argument("--out", type=Path, default=OUT_DIR, help=f"output folder (default {OUT_DIR})")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)-8s %(name)s: %(message)s",
    )
    logger.setLevel(logging.INFO)

    count = int(round((args.k_max - args.k_min) / args.k_step)) + 1
    ks = [round(args.k_min + i * args.k_step, 10) for i in range(count)]
    args.out.mkdir(parents=True, exist_ok=True)

    # The seed orbit at k = 2: every run continues from it, and the step
    # counts are read at its eigenvalue.
    seed = tanglepack.TangleSession(henon_map(SEED_K), henon_map_inverse(SEED_K), henon_jacobian(SEED_K))
    seed_fp = seed.construct_fixed_point(SEED_GUESS)
    seed.orient_eigenvectors(seed_fp, SEED_ORIENTATION)
    seed_guess = [c.tolist() for c in seed_fp.coordinates]
    seed_orientation = {
        "unstable": np.array(seed_fp.unstable_eigenvectors[0], dtype=float),
        "stable": np.array(seed_fp.stable_eigenvectors[0], dtype=float),
    }
    lambda_ref = float(seed_fp.unstable_eigenvalues[0])

    frames: list[Frame] = []
    for run in sweep_order(ks):
        guess, orientation = seed_guess, seed_orientation
        for k in run:
            started = time.time()
            try:
                frame, guess, orientation = build_frame(
                    k, guess, orientation,
                    area_cutoff=args.area_cutoff,
                    radius=args.radius,
                    unstable_steps=args.unstable_steps,
                    stable_steps=args.stable_steps,
                    scale=not args.no_scale,
                    lambda_ref=lambda_ref,
                    step_budget=args.step_budget,
                )
            except (ValueError, tanglepack.PointBudgetExceeded) as exc:
                logger.warning("k = %.3f skipped: %s: %s", k, type(exc).__name__, str(exc)[:200])
                continue
            frames.append(frame)
            print(
                f"k = {k:.3f}: {len(frame.inside):3d} lobes, mean area {frame.mean_area:.4g}, "
                f"{sum(l.head for l in frame.inside)} distinct orbits (mean {frame.mean_head_area:.4g}), "
                f"{len(frame.escaping)} escape lobes; "
                f"grown {frame.unstable_steps}/{frame.stable_steps} steps, "
                f"{len(frame.crossings)} crossings, {time.time() - started:.1f} s",
                flush=True,
            )

    if not frames:
        print("no frame built", file=sys.stderr)
        sys.exit(1)
    frames.sort(key=lambda f: f.k)

    def padded(box: tuple[float, float, float, float]) -> tuple[float, float, float, float]:
        xmin, xmax, ymin, ymax = box
        pad = args.pad * max(xmax - xmin, ymax - ymin)
        return xmin - pad, xmax + pad, ymin - pad, ymax + pad

    boxes = np.array([f.box() for f in frames])
    common = padded((boxes[:, 0].min(), boxes[:, 1].max(), boxes[:, 2].min(), boxes[:, 3].max()))

    for index, frame in enumerate(frames):
        path = args.out / f"frame_{index:03d}_k{frame.k:.3f}.png"
        limits = common if args.axes == "common" else padded(frame.box())
        draw_frame(frame, limits, path, fill=not args.no_fill)
    write_csv(frames, args.out / "lobes.csv")
    draw_summary(frames, args.out / "summary.png")
    print(f"{len(frames)} frames, lobes.csv and summary.png written to {args.out}")


if __name__ == "__main__":
    main()

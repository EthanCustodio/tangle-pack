"""
The inner period-3 resonance zone as ``k`` increases: one strong pip, tracked through the continuation, its zone blasted.

The companion of ``henon_period3_k_sweep.py`` (whose growth policy, lobe
measure and frame drawing it reuses). The difference: a strong pip is chosen
ONCE, at the seed ``k = 2``, and CARRIED through the parameter continuation
instead of being re-chosen at every ``k``. The seed choice is the CENTRAL pip
(``--seed-pip central``, the default): the pips of the window line the whole
connection from one orbit point to the next (``u(0)`` runs from ``p0`` to
``p1``, ``s(1)`` from ``p1`` back to ``p0``, and they cross all along it,
accumulating at both ends), so "the outermost registered pip" -- the
library's default, smallest unstable cdist -- is whichever one the stable
growth happens to reach, and it moves one return in or out when the growth
changes. The pip nearest the MIDPOINT of the two orbit points is intrinsic:
at ``k = 2`` it is the pip on the map's symmetry line (id 4, ``(-0.577,
0.577)``), where the two manifolds cross farthest from both periodic points.
``--seed-pip default`` takes the library default and ``--seed-pip-id`` a
candidate id (``--list-seed-pips`` prints them). At each
``k`` the resonance zone is cut at the tracked pip, blasted ``--blasts`` times
(single-iteration blasts, as ``build_period3`` does; 15 by default, the
depth of ``henon_blast_period_3.py``), and the frame shows the zone
(shaded), its stable boundary (the stable manifold trimmed at the pip orbit
by the zone), the unstable manifold TRIMMED TO THE ZONE -- every bridge
whose innermost zone is this one (``TangleSession.classify_bridge``, C4):
the boundary arcs from each orbit point to its pip and every interior
bridge, the blast children included; the escaping arm beyond the pips and
the exterior lobes are dropped (``--untrimmed`` adds the whole unstable
manifold faintly) -- the lobes of those bridges and, in its title, the
zone's area, the mean lobe area and their ratio.

The lobes are read on the UNCUT trellis (the computed one, ``strong_pip =
None`` restores it): the trellis of the pip (D10/D12) is cut at the pip
orbit's deepest REGISTERED iterate, and at ``k = 3`` that cut keeps three
lobes of a zone fifteen blasts have filled with dozens of interior bridges
-- the right trellis for the topology, the wrong picture of a blast. The
mean lobe area is therefore the mean over the zone's interior lobes, every
image a blast produced counted (iterates of one lobe have one area, so the
"distinct lobe orbits" mean is the one to compare across ``k``).

The stable manifold is grown to ESCAPE (cap: ``--stable-overgrow`` times
the scaled step target) rather than to the target alone, so the tracked pip
is always registered -- at the scaled 7 stable steps of ``k = 2.2`` the
branch stopped short of it and the tracking fell one return inward; the
zone's trim at the pip orbit cuts the surplus away again, so the trellis,
the blasts and the lobes see no more than they would have.

How the pip is tracked: a pip's SIGNATURE is the orbit index of its unstable
branch, the orbit index of its stable branch and its crossing sign -- all
three persist under continuation, since the orbit is re-solved from the
previous ``k``'s orbit in the same order. Among the strong-pip candidates with
that signature the one NEAREST the previous frame's pip (in the plane, scaled
by the orbit's diameter) is the tracked pip. Same-sign candidates on one
branch pair are one return apart (unstable cdist ratio ``lambda_u``, a
different corner of the zone), so for ``k`` steps of 0.1 the nearest is
unambiguous; a jump of more than a quarter of the orbit diameter is a WARNING
in the log and in the title. A ``k`` with no matching candidate gets a frame
without zone or blast and the tracking carries on from the last pip found.

Outputs, in ``figures/henon_period3_zone_sweep/``: ``frame_k<k>.png`` (named
by ``k``, zero-padded, so they sort in ``k`` order), ``zones.csv`` (zone area,
lobe statistics, blast sizes, the pip per ``k``), ``summary.png`` (zone area,
mean lobe area and their ratio against ``k``) and ``geometry/k<k>.npz`` (every
curve and number of the frame, NaN-separated polylines, for an animation).

Chunked runs: a sweep over a sub-range (``--k-min`` above the seed, or
``--k-max`` below it) first WARMS UP from the seed to the chunk's first ``k``
in steps of ``--warmup-step`` (default 0.1) -- orbit, orientation and tracked
pip carried, nothing blasted or drawn -- so a chunk continues the same pip as
the whole sweep would. With ``--chunk`` the CSV is named after the range and
no summary is drawn; ``--merge`` then concatenates every ``zones_k*.csv`` of
the output folder into ``zones.csv`` and draws the summary. Frames and
geometry files of several chunks share the folder by their ``k`` names. The
blast cost grows ~1.4x per blast at high ``k`` (``k = 4``: ~2 min a frame),
so a dense sweep is best run as parallel chunks of equal cost.

Run with ``PYTHONPATH=src python3 scripts/henon_period3_zone_sweep.py``
(``--k-min``, ``--k-max``, ``--k-step``, ``--blasts``, ``--min-separation``,
``--seed-pip-id``, ``--list-seed-pips``, plus the growth and drawing options of
``henon_period3_k_sweep.py``).

Dev Notes:
    * The zone is defined on the session (``session.resonance_zone``), which
      trims the stable manifold at the pip orbit and remembers the pip, so the
      blasts and the final trellis read the zone's cut (D10: ``trellis(fp)``
      is cut at the pip orbit on every read). The lobes are therefore the
      lobes of the trellis of the pip, blast children included -- the
      forward images of the zone's turnstile lobes inside the zone -- and the
      mean lobe area is reported next to the zone area (``lobe/zone``, the
      share of the zone one lobe carries). Escape lobes are excluded as in the
      companion script.
    * At the seed the candidates with id 3 and 4 are the two pips of the
      window on ``u(0) x s(1)`` (signs ``+1`` / ``-1``); 3 is the library
      default at 9 stable steps (the fixture's pip, zone area 0.538), 4 the
      central one (zone area 0.538 as well, to three digits: zones cut at
      different pips of one window differ by lobes). ``--list-seed-pips``
      prints them all.
    * A blast that fails (point budget, a bridge image that cannot be cut) is
      logged; the frame keeps the zone and the lobes of the blasts that ran.
    * The dense sweep of 2026-10-14 (``figures/henon_period3_zone_sweep_dense/``,
      148 frames at ``k`` step 1/70 from 1.9 to 4.0, 15 blasts, nine parallel
      chunks, ~25 min wall) tracks the central pip within 0.02 orbit diameters
      per step and has no incomplete zone. Its three lowest frames needed a
      finer cutoff (``k = 1.9`` at 1e-9, 1.914 and 1.929 at 1e-8, recorded in
      the ``area_cutoff`` column): at 1e-7 the pip crossing there resolves at
      ratio ~2, below the threshold of 3, so it is no candidate and the
      tracker falls one return inward (zone area continuous, the star jumps).
      ``zone_sweep.mp4`` in that folder is the Manim animation of its
      geometry files (``~/.manim-skill/henon-period3-zone/script.py``: zone,
      stable manifold, orbit and pip morph between frames; bridges and lobes,
      which change combinatorially, cross-fade).
    * Blast cost: the interior frontier grows by ~1.4 per blast at high ``k``
      (``k = 4``: 7 -> 748 interior bridges over 15 blasts, ~2 min; ``k = 3``:
      7 -> 35, 1 s), so the sweep's high end dominates its run time.
    * At ``k = 1.8`` (1e-7, any stable growth 12-24 steps) the tracked pip's
      registered image has the OPPOSITE crossing sign and its second image is
      missing, so the zone gets two cuts of three and its boundary closes with
      a chord (area 0.30 against ~0.43 expected by continuity, no lobe). The
      iterate inference mislinks there because the lobes (~1e-6) sit a decade
      above the cutoff and the crossings crowd: the companion script's
      resolution floor, seen from the iterate table. Such a frame is flagged
      INCOMPLETE in its title, in ``zones.csv`` (``zone_complete``) and on the
      summary; the default range starts at 1.9. Measured on 2026-10-14
      (1.9 to 4.0, 22 frames, 15 blasts, ~7 min, the last five frames ~5 of
      them): the pip moves 0.01-0.02 orbit diameters per step, the zone area
      runs 0.43 -> 3.74, the mean lobe area 1.6e-5 -> 0.49 and their ratio
      3.7e-5 -> 0.13, monotone but for ``k = 3.2`` (the blast discovers
      smaller secondary lobes there: 203 lobes against 38 at 3.1).
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

sys.path.insert(0, str(Path(__file__).resolve().parent))

import tanglepack  # noqa: E402
from tanglepack.examples.henon import henon_jacobian, henon_map, henon_map_inverse  # noqa: E402

from henon_period3_k_sweep import (  # noqa: E402
    LOBE_COLOR,
    SEED_GUESS,
    SEED_K,
    SEED_ORIENTATION,
    STABLE_COLOR,
    UNSTABLE_COLOR,
    Lobe,
    clip_polyline,
    grow_until_escape,
    lobe_areas,
    sweep_order,
    tangle_box,
)

logger = logging.getLogger("henon_period3_zone_sweep")

OUT_DIR = Path(__file__).resolve().parent.parent / "figures" / "henon_period3_zone_sweep"
ZONE_COLOR = "#7fb3d5"
BRIDGE_COLOR = "#1b4f72"
PIP_COLOR = "#8e44ad"

#: A tracked pip may move at most this fraction of the orbit's diameter per ``k`` step.
MAX_PIP_JUMP = 0.25


@dataclass(frozen=True)
class PipSignature:
    """
    What identifies one strong pip across the continuation.

    Attributes:
        unstable_orbit: Orbit index of the pip's unstable branch.
        stable_orbit: Orbit index of the pip's stable branch.
        sign: The crossing sign.
        point: The pip's coordinates at the last ``k`` it was found.
    """

    unstable_orbit: int
    stable_orbit: int
    sign: int
    point: tuple[float, float]

    def label(self) -> str:
        return f"u({self.unstable_orbit}) x s({self.stable_orbit}), sign {self.sign:+d}"


def signature_of(trellis, intersection_id: int) -> PipSignature:
    """The signature of a crossing of a trellis."""
    crossing = trellis.intersection(intersection_id)
    unstable_key = crossing.manifold_a_key if crossing.manifold_a_key[1] == "unstable" else crossing.manifold_b_key
    stable_key = crossing.manifold_b_key if unstable_key is crossing.manifold_a_key else crossing.manifold_a_key
    x, y = crossing.get_point()
    return PipSignature(int(unstable_key[2]), int(stable_key[2]), int(crossing.crossing_sign), (float(x), float(y)))


def central_pip(trellis, orbit: NDArray[np.float64], branch_pair: Optional[tuple[int, int]] = None) -> Optional[int]:
    """
    The strong-pip candidate nearest the midpoint of the orbit points it connects.

    Args:
        trellis: The trellis whose candidates are searched.
        orbit: The ``(period, 2)`` orbit points.
        branch_pair: ``(unstable orbit index, stable orbit index)`` to search
            on; default: the branch pair of the trellis's default pip.

    Returns:
        The candidate id, or None without candidates.
    """
    candidates = list(trellis.strong_pip_candidates)
    if not candidates:
        return None
    if branch_pair is None:
        default = trellis.strong_pip if trellis.strong_pip is not None else min(
            candidates, key=lambda i: trellis.intersection(i).unstable_cdist
        )
        signature = signature_of(trellis, default)
        branch_pair = (signature.unstable_orbit, signature.stable_orbit)
    midpoint = 0.5 * (orbit[branch_pair[0]] + orbit[branch_pair[1]])
    on_pair = [
        (candidate, signature_of(trellis, candidate))
        for candidate in candidates
    ]
    on_pair = [(c, sig) for c, sig in on_pair if (sig.unstable_orbit, sig.stable_orbit) == tuple(branch_pair)]
    if not on_pair:
        return None
    return min(on_pair, key=lambda item: float(np.hypot(*(np.subtract(item[1].point, midpoint)))))[0]


def track_pip(trellis, wanted: PipSignature, orbit: NDArray[np.float64]) -> tuple[Optional[int], float]:
    """
    The strong-pip candidate continuing ``wanted`` on this trellis.

    Args:
        trellis: The trellis whose candidates are searched.
        wanted: The signature from the previous ``k``.
        orbit: The ``(3, 2)`` orbit points, whose diameter scales the jump.

    Returns:
        ``(candidate id or None, jump)`` where ``jump`` is the distance from
        the previous pip as a fraction of the orbit diameter (nan when None).
    """
    diameter = max(float(np.max(np.linalg.norm(orbit[:, None, :] - orbit[None, :, :], axis=-1))), 1e-12)
    best: Optional[int] = None
    best_jump = math.inf
    for candidate in trellis.strong_pip_candidates:
        signature = signature_of(trellis, candidate)
        if (signature.unstable_orbit, signature.stable_orbit, signature.sign) != (
            wanted.unstable_orbit, wanted.stable_orbit, wanted.sign
        ):
            continue
        jump = float(np.hypot(*(np.subtract(signature.point, wanted.point)))) / diameter
        if jump < best_jump:
            best, best_jump = candidate, jump
    return best, (best_jump if best is not None else math.nan)


@dataclass
class ZoneFrame:
    """
    One ``k`` of the zone sweep.

    Attributes:
        k: The nonlinearity.
        orbit: The ``(3, 2)`` orbit points.
        lambda_u: The unstable eigenvalue of the 3-cycle.
        unstable_steps: Map steps the unstable manifold was grown.
        stable_steps: Map steps the stable manifold was grown.
        stable: The stable polylines (trimmed at the pip orbit), clipped.
        unstable: The whole unstable polylines, clipped; empty unless the
            frame was built ``untrimmed``.
        bridges: Every interior bridge polyline (boundary arcs and blast
            children), clipped: the unstable manifold trimmed to the zone.
        crossings: The ``(M, 2)`` crossings of the final trellis.
        lobes: Every lobe of the final trellis.
        radius: The escape radius.
        area_cutoff: The refinement cutoff.
        pip: The tracked pip's signature at this ``k``, None when not found.
        pip_jump: How far the pip moved, as a fraction of the orbit diameter.
        zone_area: The resonance zone's area, nan without a zone.
        zone_boundary: The zone's closed boundary polygon, empty without one.
        zone_cuts: Cut points of the zone, one per branch when complete.
        zone_complete: True when every branch is cut at a registered point
            of the pip's orbit with the pip's own crossing sign (the map
            preserves orientation, so an image with the other sign is a
            mislink at the resolution floor and the zone's area is unreliable).
        blast_sizes: Interior bridges each blast produced.
        blast_error: The message of a failed blast, if any.
        escaped: Per stability, whether growth stopped on escape.
        target: Per stability, the step target.
    """

    k: float
    orbit: NDArray[np.float64]
    lambda_u: float
    unstable_steps: int
    stable_steps: int
    stable: list[NDArray[np.float64]]
    unstable: list[NDArray[np.float64]]
    bridges: list[NDArray[np.float64]]
    crossings: NDArray[np.float64]
    lobes: list[Lobe]
    radius: float
    area_cutoff: float
    pip: Optional[PipSignature]
    pip_jump: float
    zone_area: float
    zone_boundary: NDArray[np.float64]
    zone_cuts: int = 0
    zone_complete: bool = False
    blast_sizes: list[int] = field(default_factory=list)
    blast_error: Optional[str] = None
    escaped: dict[str, bool] = field(default_factory=dict)
    target: dict[str, int] = field(default_factory=dict)

    @property
    def inside(self) -> list[Lobe]:
        return [lobe for lobe in self.lobes if lobe.reach <= self.radius]

    @property
    def escaping(self) -> list[Lobe]:
        return [lobe for lobe in self.lobes if lobe.reach > self.radius]

    @property
    def areas(self) -> NDArray[np.float64]:
        return np.array([lobe.area for lobe in self.inside], dtype=float)

    @property
    def mean_area(self) -> float:
        areas = self.areas
        return float(areas.mean()) if len(areas) else float("nan")

    @property
    def mean_head_area(self) -> float:
        heads = [lobe.area for lobe in self.inside if lobe.head]
        return float(np.mean(heads)) if heads else float("nan")

    @property
    def lobe_over_zone(self) -> float:
        return self.mean_area / self.zone_area if self.zone_area > 0 else float("nan")

    @property
    def at_resolution_floor(self) -> bool:
        areas = self.areas
        return bool(len(areas)) and float(np.median(areas)) < 10.0 * self.area_cutoff

    def box(self) -> tuple[float, float, float, float]:
        rings = [lobe.ring for lobe in self.inside]
        if len(self.zone_boundary):
            rings.append(self.zone_boundary)
        return tangle_box(self.orbit, self.crossings, rings)


def build_zone_frame(
    k: float,
    guess,
    orientation: dict[str, NDArray[np.float64]],
    wanted: Optional[PipSignature],
    *,
    area_cutoff: float,
    radius: Optional[float],
    unstable_steps: int,
    stable_steps: int,
    scale: bool,
    lambda_ref: float,
    step_budget: int,
    blasts: int,
    min_separation: float,
    stable_overgrow: float = 2.0,
    seed_pip: str = "central",
    seed_pip_id: Optional[int] = None,
    untrimmed: bool = False,
) -> tuple[ZoneFrame, list, dict[str, NDArray[np.float64]]]:
    """
    Build the period-3 tangle at one ``k``, cut its zone at the tracked pip and blast it.

    Args:
        k: The nonlinearity.
        guess: The ``(3, 2)`` orbit guess (the previous frame's orbit).
        orientation: Eigenvector orientation hints (the previous frame's).
        wanted: The pip signature to continue; None at the seed, where the
            pip is chosen (``seed_pip_id`` or the trellis's default).
        area_cutoff, radius, unstable_steps, stable_steps, scale, lambda_ref,
            step_budget: As ``henon_period3_k_sweep.build_frame``.
        blasts: Single-iteration blasts of the zone.
        min_separation: The blasts' ``min_separation``.
        stable_overgrow: The stable manifold's step cap as a multiple of its
            scaled target: it is grown to escape, so the tracked pip is
            registered whatever its stable cdist.
        seed_pip: At the seed, ``"central"`` (nearest the midpoint of the
            orbit points it connects) or ``"default"`` (the library's).
        seed_pip_id: At the seed, a candidate id overriding ``seed_pip``.
        untrimmed: Also keep the whole unstable manifold for a faint underlay.

    Returns:
        ``(frame, orbit, oriented eigenvectors)``.
    """
    session = tanglepack.TangleSession(
        henon_map(k), henon_map_inverse(k), henon_jacobian(k), area_cutoff=area_cutoff
    )
    fp = session.construct_fixed_point(guess)
    session.orient_eigenvectors(fp, orientation)
    session.initialize_both_manifolds(fp)
    lambda_u = float(fp.unstable_eigenvalues[0])
    escape_radius = radius if radius is not None else 2.0 * (1.0 + math.sqrt(1.0 + k))
    orbit = np.array(fp.coordinates, dtype=float)

    factor = math.log(lambda_ref) / math.log(lambda_u) if scale else 1.0
    steps: dict[str, int] = {}
    escaped: dict[str, bool] = {}
    targets: dict[str, int] = {}
    for stability, count in (("unstable", unstable_steps), ("stable", stable_steps)):
        targets[stability] = max(1, int(round(count * factor)))
        cap = targets[stability]
        if stability == "stable":
            cap = max(cap, int(math.ceil(cap * stable_overgrow)))
            targets[stability] = cap
        steps[stability], escaped[stability] = grow_until_escape(
            session, fp, stability, radius=escape_radius, max_steps=cap, step_budget=step_budget,
        )
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    trellis = session.trellis(fp)

    # ── the pip ────────────────────────────────────────────────────────
    pip_jump = math.nan
    if wanted is None:
        session.classify_strong_pips(fp)
        if seed_pip_id is not None:
            pip_id = seed_pip_id
        elif seed_pip == "central":
            pip_id = central_pip(trellis, orbit)
        else:
            pip_id = trellis.strong_pip
        if pip_id is not None and pip_id != trellis.strong_pip:
            session.set_strong_pip(fp, pip_id)
    else:
        pip_id, pip_jump = track_pip(trellis, wanted, orbit)
        if pip_id is not None:
            session.set_strong_pip(fp, pip_id)
            if pip_jump > MAX_PIP_JUMP:
                logger.warning(
                    "k = %.3f: the tracked pip jumped %.2f orbit diameters (%s)",
                    k, pip_jump, wanted.label(),
                )
        else:
            logger.warning("k = %.3f: no strong-pip candidate continues %s", k, wanted.label())
    pip = signature_of(trellis, pip_id) if pip_id is not None else None

    # ── the zone and its blasts ───────────────────────────────────────
    zone_area = math.nan
    zone_boundary = np.empty((0, 2))
    zone_cuts = 0
    zone_complete = False
    sizes: list[int] = []
    blast_error: Optional[str] = None
    if pip_id is not None:
        zone = session.resonance_zone(pip_id)
        zone_area = float(zone.area)
        zone_boundary = np.asarray(zone.boundary_polygon(session.workbench), dtype=float).reshape(-1, 2)
        zone_cuts = len(zone.cut_intersections)
        pip_sign = int(trellis.intersection(pip_id).crossing_sign)
        zone_complete = zone_cuts == fp.k_value and all(
            int(cut.crossing_sign) == pip_sign for cut in zone.cut_intersections
        )
        if not zone_complete:
            logger.warning(
                "k = %.3f: the zone is incomplete (%d of %d cuts, signs %s against the pip's %+d): "
                "the pip's orbit is not reliably registered, its area is unreliable",
                k, zone_cuts, fp.k_value, [int(c.crossing_sign) for c in zone.cut_intersections], pip_sign,
            )
        for _ in range(blasts):
            try:
                result = session.blast_zone(
                    zone, num_iterations=1, fixed_point=[fp], min_separation=min_separation
                )
            except (ValueError, AssertionError, tanglepack.PointBudgetExceeded) as exc:
                blast_error = f"{type(exc).__name__}: {str(exc)[:120]}"
                logger.warning("k = %.3f: blast %d failed: %s", k, len(sizes) + 1, blast_error)
                break
            sizes.append(len(result.all_interior_bridges()))
            if result.stop_reason == "frontier_empty":
                break
        trellis = session.trellis(fp)
        chosen = session.chosen_strong_pip(fp)
        if chosen is not None:
            found = trellis.registry.find(chosen)
            if found is not None:
                pip = signature_of(trellis, found)
        # The uncut trellis: the blast's bridges lie past the pip orbit's
        # deepest registered iterate, where the trellis of the pip ends.
        trellis.strong_pip = None
        interior = [b for b in trellis.bridges if b.id is not None and session.classify_bridge(b) is zone]
    else:
        interior = [b for b in trellis.bridges if b.id is not None]

    lobes = lobe_areas(trellis, interior)
    crossings = np.array(
        [trellis.intersection(i).get_point() for i in trellis.intersection_ids if not trellis.intersection(i).is_anchor],
        dtype=float,
    ).reshape(-1, 2)
    rings = [lobe.ring for lobe in lobes if lobe.reach <= escape_radius]
    if len(zone_boundary):
        rings.append(zone_boundary)
    xmin, xmax, ymin, ymax = tangle_box(orbit, crossings, rings)
    cx, cy = 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)
    half = max(xmax - xmin, ymax - ymin, 1e-3) * 2.0
    clip_box = (cx - half, cx + half, cy - half, cy + half)

    def curves(stability: str) -> list[NDArray[np.float64]]:
        return [
            clip_polyline(manifold.get_point_array(), clip_box)
            for (owner, stab, _o, _b), manifold in session.manifolds.items()
            if owner is fp and stab == stability
        ]

    bridges = [
        clip_polyline(np.asarray(bridge.get_point_array(), dtype=float), clip_box)
        for bridge in interior
        if bridge.get_point_array() is not None and len(bridge.get_point_array()) >= 2
    ]
    frame = ZoneFrame(
        k=k,
        orbit=orbit,
        lambda_u=lambda_u,
        unstable_steps=steps["unstable"],
        stable_steps=steps["stable"],
        stable=curves("stable"),
        unstable=curves("unstable") if untrimmed else [],
        bridges=bridges,
        crossings=crossings,
        lobes=lobes,
        radius=escape_radius,
        area_cutoff=area_cutoff,
        pip=pip,
        pip_jump=pip_jump,
        zone_area=zone_area,
        zone_boundary=zone_boundary,
        zone_cuts=zone_cuts,
        zone_complete=zone_complete,
        blast_sizes=sizes,
        blast_error=blast_error,
        escaped=escaped,
        target=targets,
    )
    oriented = {
        "unstable": np.array(fp.unstable_eigenvectors[0], dtype=float),
        "stable": np.array(fp.stable_eigenvectors[0], dtype=float),
    }
    return frame, [c.tolist() for c in fp.coordinates], oriented


def _nan_joined(polylines: list[NDArray[np.float64]]) -> NDArray[np.float64]:
    """Several ``(N, 2)`` polylines as one array, a NaN row between them."""
    if not polylines:
        return np.empty((0, 2))
    parts: list[NDArray[np.float64]] = []
    for poly in polylines:
        parts.append(np.asarray(poly, dtype=float).reshape(-1, 2))
        parts.append(np.full((1, 2), np.nan))
    return np.vstack(parts)[:-1]


def save_geometry(frame: ZoneFrame, path: Path) -> None:
    """Write the frame's curves and numbers to an ``.npz`` (NaN-separated polylines)."""
    inside = frame.inside
    np.savez_compressed(
        path,
        k=frame.k,
        orbit=frame.orbit,
        stable=_nan_joined(frame.stable),
        unstable=_nan_joined(frame.unstable),
        bridges=_nan_joined(frame.bridges),
        crossings=frame.crossings,
        zone_boundary=frame.zone_boundary,
        lobes=_nan_joined([lobe.ring for lobe in inside]),
        lobe_areas=np.array([lobe.area for lobe in inside], dtype=float),
        pip=np.array(frame.pip.point if frame.pip else [np.nan, np.nan], dtype=float),
        zone_area=frame.zone_area,
        zone_complete=int(frame.zone_complete),
        mean_lobe_area=frame.mean_area,
        mean_orbit_lobe_area=frame.mean_head_area,
        lobe_over_zone=frame.lobe_over_zone,
        n_lobes=len(inside),
        n_orbits=int(sum(lobe.head for lobe in inside)),
        blasts=len(frame.blast_sizes),
        interior_bridges=int(sum(frame.blast_sizes)),
        lambda_u=frame.lambda_u,
    )


def list_seed_pips() -> None:
    """Print the seed trellis's strong-pip candidates with their signatures."""
    session = tanglepack.TangleSession(henon_map(SEED_K), henon_map_inverse(SEED_K), henon_jacobian(SEED_K))
    fp = session.construct_fixed_point(SEED_GUESS)
    session.orient_eigenvectors(fp, SEED_ORIENTATION)
    session.initialize_both_manifolds(fp)
    radius = 2.0 * (1.0 + math.sqrt(1.0 + SEED_K))
    session.grow_n_times(fp, "unstable", num_iterations=13)
    grow_until_escape(session, fp, "stable", radius=radius, max_steps=18, step_budget=500_000)
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    session.classify_strong_pips(fp)
    trellis = session.trellis(fp)
    orbit = np.array(fp.coordinates, dtype=float)
    print(f"seed k = {SEED_K:g}: library default pip {trellis.strong_pip}, central pip {central_pip(trellis, orbit)}; candidates:")
    for candidate in trellis.strong_pip_candidates:
        signature = signature_of(trellis, candidate)
        crossing = trellis.intersection(candidate)
        print(
            f"  id {candidate:3d}  {signature.label()}  at ({signature.point[0]:+.4f}, {signature.point[1]:+.4f})"
            f"  unstable cdist {crossing.unstable_cdist:.4g}  stable cdist {crossing.stable_cdist:.4g}"
        )


# ── drawing ────────────────────────────────────────────────────────────────


def draw_zone_frame(frame: ZoneFrame, limits: tuple[float, float, float, float], path: Path, *, fill: bool) -> None:
    """Draw one zone frame to ``path``."""
    fig, ax = plt.subplots(figsize=(9, 9))
    if len(frame.zone_boundary):
        ax.fill(frame.zone_boundary[:, 0], frame.zone_boundary[:, 1], color=ZONE_COLOR, alpha=0.25, lw=0)
    if fill:
        for lobe in frame.escaping:
            ax.fill(lobe.ring[:, 0], lobe.ring[:, 1], color="grey", alpha=0.12, lw=0)
        for lobe in frame.inside:
            ax.fill(lobe.ring[:, 0], lobe.ring[:, 1], color=LOBE_COLOR, alpha=0.3, lw=0)
    for poly in frame.unstable:
        ax.plot(poly[:, 0], poly[:, 1], color=UNSTABLE_COLOR, lw=0.5, alpha=0.3)
    for poly in frame.bridges:
        ax.plot(poly[:, 0], poly[:, 1], color=UNSTABLE_COLOR, lw=0.8)
    for poly in frame.stable:
        ax.plot(poly[:, 0], poly[:, 1], color=STABLE_COLOR, lw=0.9)
    if len(frame.crossings):
        ax.scatter(frame.crossings[:, 0], frame.crossings[:, 1], s=6, color="k", zorder=4)
    ax.scatter(frame.orbit[:, 0], frame.orbit[:, 1], s=40, color="k", zorder=5)
    for i, (x, y) in enumerate(frame.orbit):
        ax.annotate(f"p{i}", (x, y), textcoords="offset points", xytext=(5, 5), fontsize=9)
    if frame.pip is not None:
        ax.scatter(*frame.pip.point, s=140, marker="*", color=PIP_COLOR, edgecolor="k", zorder=6, label="tracked pip")
        ax.legend(loc="lower right", fontsize=9)
    ax.set_xlim(limits[0], limits[1])
    ax.set_ylim(limits[2], limits[3])
    ax.set_aspect("equal")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    inside, escaping = frame.inside, frame.escaping
    lines = [f"Hénon k = {frame.k:.3f}, b = 1 — inner period-3 resonance zone, pip tracked from k = {SEED_K:g}"]
    if frame.pip is None:
        lines.append("no strong-pip candidate continues the tracked pip: no zone, no blast")
    else:
        zone = f"zone area = {frame.zone_area:.4g}"
        if not frame.zone_complete:
            zone += f"  ⚠ INCOMPLETE ({frame.zone_cuts} cut(s)): the pip's orbit is not reliably registered"
        if inside:
            areas = frame.areas
            lines.append(
                f"{zone};  mean lobe area = {frame.mean_area:.4g};  lobe/zone = {frame.lobe_over_zone:.3g}"
            )
            lines.append(
                f"{len(inside)} lobes (median {np.median(areas):.3g}, min {areas.min():.3g}, max {areas.max():.3g});  "
                f"{sum(lobe.head for lobe in inside)} distinct lobe orbits, mean area {frame.mean_head_area:.4g}"
            )
        else:
            lines.append(f"{zone};  no lobe inside the trapped region")
        jump = "" if math.isnan(frame.pip_jump) else f", moved {frame.pip_jump:.2f} orbit diameters"
        flag = "  ⚠ large jump" if frame.pip_jump > MAX_PIP_JUMP else ""
        pip = frame.pip
        lines.append(
            f"pip {pip.label()} at ({pip.point[0]:+.3f}, {pip.point[1]:+.3f}){jump}{flag};  "
            f"{len(frame.blast_sizes)} blast(s), {sum(frame.blast_sizes)} interior bridges"
            + (f" (last {frame.blast_sizes[-1]})" if frame.blast_sizes else "")
            + (f";  blast failed: {frame.blast_error}" if frame.blast_error else "")
        )
    if escaping:
        lines.append(
            f"{len(escaping)} escape lobe(s) beyond |coord| = {frame.radius:.3g} excluded "
            f"(mean area {np.mean([l.area for l in escaping]):.3g}; grey)"
        )
    if frame.at_resolution_floor:
        lines.append(f"⚠ lobe areas within a decade of area_cutoff = {frame.area_cutoff:.0e}: at the resolution floor")
    stop = {
        s: ("escape" if frame.escaped.get(s)
            else "cap" if getattr(frame, s + "_steps") == frame.target.get(s)
            else "budget")
        for s in ("unstable", "stable")
    }
    lines.append(
        f"λ_u = {frame.lambda_u:.3g};  grown {frame.unstable_steps} unstable ({stop['unstable']}) / "
        f"{frame.stable_steps} stable ({stop['stable']}) steps;  {len(frame.crossings)} crossings"
    )
    ax.set_title("\n".join(lines), fontsize=9.5, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


CSV_COLUMNS = [
    "k", "lambda_u", "unstable_steps", "stable_steps", "pip", "pip_x", "pip_y", "pip_jump",
    "zone_area", "zone_cuts", "zone_complete", "crossings", "lobes", "distinct_lobe_orbits", "mean_lobe_area",
    "mean_orbit_lobe_area", "median_lobe_area", "lobe_over_zone", "escape_lobes",
    "blasts", "interior_bridges", "blast_sizes", "blast_error", "at_resolution_floor", "area_cutoff",
]


def csv_row(f: ZoneFrame) -> dict[str, str]:
    """The CSV row of one frame."""
    areas = f.areas
    values = [
        f"{f.k:.4f}", f"{f.lambda_u:.6g}", f.unstable_steps, f.stable_steps,
        f.pip.label() if f.pip else "", f"{f.pip.point[0]:.6g}" if f.pip else "",
        f"{f.pip.point[1]:.6g}" if f.pip else "", f"{f.pip_jump:.4g}",
        f"{f.zone_area:.6g}", f.zone_cuts, int(f.zone_complete), len(f.crossings), len(f.inside),
        int(sum(l.head for l in f.inside)), f"{f.mean_area:.6g}", f"{f.mean_head_area:.6g}",
        f"{np.median(areas):.6g}" if len(areas) else "", f"{f.lobe_over_zone:.6g}",
        len(f.escaping), len(f.blast_sizes), int(sum(f.blast_sizes)), " ".join(map(str, f.blast_sizes)),
        f.blast_error or "", int(f.at_resolution_floor), f"{f.area_cutoff:g}",
    ]
    return {column: str(value) for column, value in zip(CSV_COLUMNS, values)}


def write_csv(rows: list[dict[str, str]], path: Path) -> None:
    """One row per ``k``."""
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=CSV_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def draw_summary(rows: list[dict[str, str]], path: Path) -> None:
    """Zone area, mean lobe area and their ratio against ``k``, from CSV rows."""
    def column(name: str) -> list[float]:
        return [float(row[name]) if row[name] not in ("", "nan") else float("nan") for row in rows]

    ks = column("k")
    fig, (ax, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 9), sharex=True)
    ax.plot(ks, column("zone_area"), "o-", color=ZONE_COLOR, mec="k", label="zone area", ms=4)
    bad = [row for row in rows if row["zone_complete"] == "0"]
    if bad:
        ax.plot([float(r["k"]) for r in bad], [float(r["zone_area"]) for r in bad], "x", color="r", ms=10, label="incomplete zone")
    ax.plot(ks, column("mean_lobe_area"), "s-", color=UNSTABLE_COLOR, label="mean lobe area", ms=4)
    ax.plot(ks, column("mean_orbit_lobe_area"), "^--", color=STABLE_COLOR, label="mean over distinct lobe orbits", ms=4)
    ax.set_yscale("log")
    ax.set_ylabel("area")
    ax.legend()
    ax.set_title("Inner period-3 Hénon resonance zone against k (one pip, tracked)", loc="left")
    ax2.plot(ks, column("lobe_over_zone"), "o-", color="k", ms=4)
    ax2.set_yscale("log")
    ax2.set_ylabel("mean lobe area / zone area")
    ax3.plot(ks, column("lobes"), "o-", color="k", label="lobes", ms=4)
    ax3.plot(ks, column("interior_bridges"), "s--", color=BRIDGE_COLOR, label="interior bridges (all blasts)", ms=4)
    ax3.plot(ks, column("blasts"), "^:", color="grey", label="blasts run (stops when the frontier empties)", ms=4)
    ax3.set_ylabel("count")
    ax3.set_xlabel("k")
    ax3.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def merge_chunks(out: Path) -> None:
    """Concatenate every ``zones_k*.csv`` of ``out`` into ``zones.csv`` and draw the summary."""
    rows: list[dict[str, str]] = []
    for chunk in sorted(out.glob("zones_k*.csv")):
        with chunk.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if "blast_sizes" not in row:  # chunks written before the column existed
                    sizes = row.get("interior_bridges", "")
                    row["blast_sizes"] = sizes
                    row["interior_bridges"] = str(sum(int(x) for x in sizes.split()))
                row.setdefault("area_cutoff", "1e-07")
                rows.append({column: row.get(column, "") for column in CSV_COLUMNS})
    if not rows:
        print(f"no zones_k*.csv in {out}", file=sys.stderr)
        sys.exit(1)
    rows.sort(key=lambda row: float(row["k"]))
    write_csv(rows, out / "zones.csv")
    draw_summary(rows, out / "summary.png")
    print(f"{len(rows)} rows merged into {out / 'zones.csv'}; summary.png drawn")


def frame_name(k: float) -> str:
    """``k02.0143``: zero-padded so names sort in ``k`` order."""
    return f"k{k:07.4f}"


# ── the sweep ──────────────────────────────────────────────────────────────


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--k-min", type=float, default=1.9, help="first k (default 1.9: below it the pip's orbit is not reliably registered)")
    parser.add_argument("--k-max", type=float, default=4.0, help="last k (default 4.0)")
    parser.add_argument("--k-step", type=float, default=0.1, help="k increment (default 0.1)")
    parser.add_argument("--blasts", type=int, default=15, help="single-iteration blasts of the zone per k (default 15)")
    parser.add_argument("--untrimmed", action="store_true",
                        help="also draw the whole unstable manifold faintly (default: only the trellis cut at the pip orbit)")
    parser.add_argument("--min-separation", type=float, default=1e-5, help="the blasts' min_separation (default 1e-5)")
    parser.add_argument("--seed-pip", choices=("central", "default"), default="central",
                        help="the seed's pip: nearest the midpoint of its connection (default) or the library default")
    parser.add_argument("--seed-pip-id", type=int, default=None, help="candidate id to choose at the seed k = 2 instead")
    parser.add_argument("--stable-overgrow", type=float, default=2.0,
                        help="stable step cap as a multiple of the scaled target; grown to escape within it (default 2)")
    parser.add_argument("--list-seed-pips", action="store_true", help="print the seed's candidates and exit")
    parser.add_argument("--warmup-step", type=float, default=0.1,
                        help="k step of the warm-up from the seed to a sub-range's first k (default 0.1)")
    parser.add_argument("--chunk", action="store_true",
                        help="a chunk of a larger sweep: CSV named after the range, no summary (see --merge)")
    parser.add_argument("--merge", action="store_true", help="merge the chunk CSVs in --out into zones.csv and draw the summary")
    parser.add_argument("--no-geometry", action="store_true", help="do not write geometry/<k>.npz per frame")
    parser.add_argument("--area-cutoff", type=float, default=1e-7, help="refinement cutoff (default 1e-7)")
    parser.add_argument("--radius", type=float, default=None, help="escape radius (default 2 * (1 + sqrt(1 + k)))")
    parser.add_argument("--unstable-steps", type=int, default=13, help="unstable steps at k = 2, scaled elsewhere (default 13)")
    parser.add_argument("--stable-steps", type=int, default=9, help="stable steps at k = 2, scaled elsewhere (default 9)")
    parser.add_argument("--no-scale", action="store_true", help="use the step counts verbatim at every k")
    parser.add_argument("--step-budget", type=int, default=500_000, help="predicted points a growth step may create (default 500000)")
    parser.add_argument("--no-fill", action="store_true", help="do not shade the lobes")
    parser.add_argument("--axes", choices=("frame", "common"), default="frame",
                        help="axes limits: each frame's own box (default) or the union of every frame's")
    parser.add_argument("--pad", type=float, default=0.08, help="axes padding as a fraction of the box")
    parser.add_argument("--out", type=Path, default=OUT_DIR, help=f"output folder (default {OUT_DIR})")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO if args.verbose else logging.WARNING,
        format="%(levelname)-8s %(name)s: %(message)s",
    )
    logger.setLevel(logging.INFO)

    if args.list_seed_pips:
        list_seed_pips()
        return
    if args.merge:
        merge_chunks(args.out)
        return

    count = int(round((args.k_max - args.k_min) / args.k_step)) + 1
    ks = [round(args.k_min + i * args.k_step, 10) for i in range(count)]
    args.out.mkdir(parents=True, exist_ok=True)

    growth = dict(
        area_cutoff=args.area_cutoff, radius=args.radius,
        unstable_steps=args.unstable_steps, stable_steps=args.stable_steps,
        scale=not args.no_scale, step_budget=args.step_budget,
        blasts=args.blasts, min_separation=args.min_separation,
        stable_overgrow=args.stable_overgrow, untrimmed=args.untrimmed,
    )

    # The seed: the orbit, its eigenvalue (the step counts are read at it) and
    # the chosen pip's signature, the one thing every run continues.
    seed = tanglepack.TangleSession(henon_map(SEED_K), henon_map_inverse(SEED_K), henon_jacobian(SEED_K))
    seed_fp = seed.construct_fixed_point(SEED_GUESS)
    lambda_ref = float(seed_fp.unstable_eigenvalues[0])
    seed_frame, seed_guess, seed_orientation = build_zone_frame(
        SEED_K, SEED_GUESS, SEED_ORIENTATION, None, lambda_ref=lambda_ref,
        seed_pip=args.seed_pip, seed_pip_id=args.seed_pip_id, **growth
    )
    if seed_frame.pip is None:
        print("the seed has no strong-pip candidate; nothing to track", file=sys.stderr)
        sys.exit(1)
    print(f"seed k = {SEED_K:g}: pip {seed_frame.pip.label()} at ({seed_frame.pip.point[0]:+.4f}, {seed_frame.pip.point[1]:+.4f}), "
          f"zone area {seed_frame.zone_area:.4g}", flush=True)

    frames: list[ZoneFrame] = [seed_frame] if SEED_K in ks else []
    for run in sweep_order([k for k in ks if k != SEED_K]):
        guess, orientation, wanted = seed_guess, seed_orientation, seed_frame.pip
        # Warm up from the seed to the run's first k: continue the orbit and the
        # pip in steps of warmup_step, blasting and drawing nothing.
        direction = 1.0 if run[0] > SEED_K else -1.0
        k_warm = SEED_K + direction * args.warmup_step
        while direction * (run[0] - k_warm) > 1e-9:
            frame, guess, orientation = build_zone_frame(
                round(k_warm, 10), guess, orientation, wanted, lambda_ref=lambda_ref,
                **{**growth, "blasts": 0},
            )
            if frame.pip is not None:
                wanted = frame.pip
            logger.info("warm-up k = %.3f: pip %s, jump %.2f", k_warm, frame.pip.label() if frame.pip else None, frame.pip_jump)
            k_warm += direction * args.warmup_step
        for k in run:
            started = time.time()
            try:
                frame, guess, orientation = build_zone_frame(
                    k, guess, orientation, wanted, lambda_ref=lambda_ref, **growth
                )
            except (ValueError, tanglepack.PointBudgetExceeded) as exc:
                logger.warning("k = %.3f skipped: %s: %s", k, type(exc).__name__, str(exc)[:200])
                continue
            frames.append(frame)
            if frame.pip is not None:
                wanted = frame.pip
            print(
                f"k = {k:.3f}: zone area {frame.zone_area:.4g}, {len(frame.inside):3d} lobes, mean area "
                f"{frame.mean_area:.4g} (lobe/zone {frame.lobe_over_zone:.3g}), pip jump {frame.pip_jump:.2f}, "
                f"blasts {frame.blast_sizes}{' FAILED' if frame.blast_error else ''}; "
                f"grown {frame.unstable_steps}/{frame.stable_steps}, {time.time() - started:.1f} s",
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
    if not args.no_geometry:
        (args.out / "geometry").mkdir(exist_ok=True)
    for frame in frames:
        limits = common if args.axes == "common" else padded(frame.box())
        draw_zone_frame(frame, limits, args.out / f"frame_{frame_name(frame.k)}.png", fill=not args.no_fill)
        if not args.no_geometry:
            save_geometry(frame, args.out / "geometry" / f"{frame_name(frame.k)}.npz")
    rows = [csv_row(f) for f in frames]
    if args.chunk:
        name = f"zones_{frame_name(frames[0].k)}_{frame_name(frames[-1].k)}.csv"
        write_csv(rows, args.out / name)
        print(f"{len(frames)} frames and {name} written to {args.out} (chunk: run --merge for zones.csv and summary.png)")
    else:
        write_csv(rows, args.out / "zones.csv")
        draw_summary(rows, args.out / "summary.png")
        print(f"{len(frames)} frames, zones.csv, summary.png and geometry/ written to {args.out}")


if __name__ == "__main__":
    main()

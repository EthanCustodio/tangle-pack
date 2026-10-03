"""Higher-period and nested Hénon tangles, built up to the stable partition.

Two recipes on the nested Hénon map ``HENON_P3 = (2, 1)``, shared by the
figure scripts and the test fixtures so the two can never drift apart:

* :func:`build_period3` -- the inner period-3 saddle orbit alone;
* :func:`build_nested` -- the outer period-1 saddle and the inner period-3
  orbit together, one resonance zone each.

Both stop after ``partition_stable_manifold``; everything downstream
(minimal trellis, iterated partition, dual graph, symbolic dynamics) is a
session call away.

Dev Notes:

* Parameters follow ``scripts/henon_blast_period_3.py`` and the
  ``henon_p3_session`` fixture: period 3 at ``area_cutoff = 1e-7`` with
  orientation unstable ``[0, -1]`` / stable ``[-1, -1]``, grown 13 unstable /
  9 stable steps; period 1 with unstable ``[-1, 0]`` / stable ``[0, 1]``,
  11 unstable steps and ``grow_until_turnaround`` on the stable side.
* Period 3 is CLOSED at 13 steps: every active hole bridge's image is already
  a hole bridge on the next branch (``a -> b -> c -> a u^-1 w^-1``), so the
  minimal trellis has no image bridges and blasting its zone changes nothing
  topologically. Growing to 15 or 16 unstable steps, or blasting 6 or more
  times, used to trip ``check_holes_share_bridge_side`` (2026-09-30): a
  propagated hole's side was read against the nearest vertex of its whole
  containing bridge, on another fold. Fixed 2026-10-02 (sides are read on the
  image sub-arc); those runs now reach the symbolic dynamics but are not
  reliable (an unreachable class, a virtual ``new1``), so the default stays
  at 13.
* Each nested zone blasts only its own fixed point (2026-10-02): blasting
  both from either zone also iterated the other tangle's bridges, made the
  result depend on blast order and collapsed the period-3 pseudoneighbors at
  6 or more outer blasts.
* The nested case's outer zone is blasted twice (the k=2.8 recipe): with no
  blast the outer class ``e`` is unresolved (a singleton landing with no
  registered chain); after two blasts every class resolves and the outer
  classes reproduce the k=2.8 two-blast words, with real (``unique``) walks.
* Registry ids are not reproducible between builds, so the pips are the
  default strong pips, never pinned by id.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:  # pragma: no cover - typing only
    from ..loom.TangleSession import TangleSession

from .henon import HENON_P3, henon_jacobian, henon_map, henon_map_inverse, saddle_guesses

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

#: Orientation hints of the period-3 orbit's eigenvectors.
PERIOD3_ORIENTATION: dict[str, NDArray[np.float64]] = {
    "unstable": np.array([0, -1]),
    "stable": np.array([-1, -1]),
}
#: Orientation hints of the outer period-1 saddle's eigenvectors.
PERIOD1_ORIENTATION: dict[str, NDArray[np.float64]] = {
    "unstable": np.array([-1, 0]),
    "stable": np.array([0, 1]),
}


@dataclass
class TangleBuild:
    """
    One built tangle, partitioned and ready for the topology layer.

    Attributes:
        session: The :class:`~tanglepack.loom.TangleSession.TangleSession`.
        fixed_points: The fixed points, outermost first.
        pips: The strong pip chosen per fixed point, in the same order.
        stable_points: Every stable-manifold node, stacked ``(N, 2)``, captured
            before the trim at the pips (for framing plots).
        title: A human title, ``"Hénon k=2, period 3"``.
        slug: A file-system name, ``"p3"``.
        blast_sizes: Interior bridges produced by each blast, in order.
    """

    session: "TangleSession"
    fixed_points: list
    pips: list[int]
    stable_points: NDArray[np.float64]
    title: str
    slug: str
    blast_sizes: list[int] = field(default_factory=list)


def _session() -> "TangleSession":
    """A fresh session on the nested Hénon map."""
    from ..loom.TangleSession import TangleSession

    return TangleSession(
        henon_map(*HENON_P3), henon_map_inverse(*HENON_P3), henon_jacobian(*HENON_P3)
    )


def _grow_period3(session, unstable_steps: int, stable_steps: int, area_cutoff: float):
    """Construct, orient, initialise and grow the period-3 orbit."""
    session.workbench._man_machine.area_cutoff = area_cutoff
    fp3 = session.construct_fixed_point(saddle_guesses(*HENON_P3)["period_3"])
    session.orient_eigenvectors(fp3, PERIOD3_ORIENTATION)
    session.initialize_both_manifolds(fp3)
    session.grow_n_times(fp3, "unstable", num_iterations=unstable_steps)
    session.grow_n_times(fp3, "stable", num_iterations=stable_steps)
    return fp3


def _grow_period1(session, unstable_steps: int, area_cutoff: float):
    """Construct, orient, initialise and grow the outer period-1 saddle."""
    session.workbench._man_machine.area_cutoff = area_cutoff
    fp1 = session.construct_fixed_point(saddle_guesses(*HENON_P3)["saddle"])
    session.orient_eigenvectors(fp1, PERIOD1_ORIENTATION)
    session.initialize_both_manifolds(fp1)
    session.grow_n_times(fp1, "unstable", num_iterations=unstable_steps)
    session.grow_until_turnaround(fp1, "stable")
    return fp1


def _stable_points(session, fixed_points) -> NDArray[np.float64]:
    """Every node of the fixed points' stable manifolds, stacked."""
    return np.vstack(
        [
            manifold.get_point_array()
            for (fixed_point, stability, _o, _b), manifold in session.manifolds.items()
            if stability == "stable" and any(fixed_point is fp for fp in fixed_points)
        ]
    )


def _intersect_and_bridge(session, fixed_points) -> None:
    """Crossings over every fixed point, then trim and bridges per fixed point."""
    session.compute_intersections(list(fixed_points))
    for fp in fixed_points:
        session.trim_stable_manifolds(fp)
    for fp in fixed_points:
        session.create_bridges(fp)
    session.infer_iterate_table()


def _repin(session, fixed_points, pips) -> None:
    """Classify again (the trellis was rebuilt) and restore the chosen pips."""
    session.classify_strong_pips()
    for fp, pip in zip(fixed_points, pips):
        session.set_strong_pip(fp, pip)


def _partition(session) -> None:
    """Pseudoneighbors, holes and the stable partition on every trellis."""
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()


def build_period3(
    *,
    unstable_steps: int = 13,
    stable_steps: int = 9,
    area_cutoff: float = 1e-7,
    blasts: int = 0,
    min_separation: float = 1e-5,
) -> TangleBuild:
    """
    The inner period-3 orbit of the nested Hénon map, alone.

    Args:
        unstable_steps: Map steps to grow the unstable manifold.
        stable_steps: Map steps to grow the stable manifold.
        area_cutoff: Refinement area cutoff for the manifold machine.
        blasts: Single-iteration blasts of the resonance zone after the trim.
        min_separation: The blasts' ``min_separation``.

    Returns:
        The :class:`TangleBuild` (one fixed point), partitioned.
    """
    session = _session()
    fp3 = _grow_period3(session, unstable_steps, stable_steps, area_cutoff)
    _intersect_and_bridge(session, [fp3])
    stable_points = _stable_points(session, [fp3])

    session.classify_strong_pips()
    pip = session.trellis(fp3).strong_pip
    zone = session.resonance_zone(pip)
    _repin(session, [fp3], [pip])
    sizes: list[int] = []
    for _ in range(blasts):
        result = session.blast_zone(
            zone, num_iterations=1, fixed_point=[fp3], min_separation=min_separation
        )
        sizes.append(len(result.all_interior_bridges()))
        _repin(session, [fp3], [pip])
    _partition(session)
    blasted = f", {blasts} blast(s)" if blasts else ""
    return TangleBuild(
        session=session,
        fixed_points=[fp3],
        pips=[pip],
        stable_points=stable_points,
        title=f"Hénon k={HENON_P3[0]:g}, period 3{blasted}",
        slug="p3" + (f"_blast{blasts}" if blasts else ""),
        blast_sizes=sizes,
    )


def build_nested(
    *,
    p3_unstable_steps: int = 13,
    p3_stable_steps: int = 9,
    p3_area_cutoff: float = 1e-7,
    p1_unstable_steps: int = 11,
    p1_area_cutoff: float = 1e-7,
    outer_blasts: int = 2,
    inner_blasts: int = 0,
    min_separation: float = 1e-4,
) -> TangleBuild:
    """
    The outer period-1 saddle and the inner period-3 orbit, together.

    Args:
        p3_unstable_steps: Unstable growth of the period-3 orbit.
        p3_stable_steps: Stable growth of the period-3 orbit.
        p3_area_cutoff: Area cutoff while growing the period-3 manifolds.
        p1_unstable_steps: Unstable growth of the period-1 saddle.
        p1_area_cutoff: Area cutoff while growing the period-1 manifolds
            (``henon_blast_period_3.py`` uses 1e-7; the ``henon_p3_session``
            fixture uses 1e-4).
        outer_blasts: Single-iteration blasts of the outer (larger) zone.
        inner_blasts: Single-iteration blasts of the inner (smaller) zone.
            Each zone's blast touches only its own tangle, so the order of
            the inner and outer blasts does not matter.
        min_separation: The blasts' ``min_separation``.

    Returns:
        The :class:`TangleBuild`, fixed points ``[fp1, fp3]`` (outermost
        first), partitioned.
    """
    session = _session()
    fp3 = _grow_period3(session, p3_unstable_steps, p3_stable_steps, p3_area_cutoff)
    fp1 = _grow_period1(session, p1_unstable_steps, p1_area_cutoff)
    fixed_points = [fp1, fp3]
    _intersect_and_bridge(session, [fp3, fp1])
    stable_points = _stable_points(session, fixed_points)

    session.classify_strong_pips()
    pips = [session.trellis(fp).strong_pip for fp in fixed_points]
    session.add_resonance_zones(pips)
    zones = sorted(session.resonance_zones.values(), key=lambda zone: zone.area)
    inner, outer = zones[0], zones[-1]
    _repin(session, fixed_points, pips)
    sizes: list[int] = []
    for zone, count in ((inner, inner_blasts), (outer, outer_blasts)):
        for _ in range(count):
            result = session.blast_zone(
                zone,
                num_iterations=1,
                fixed_point=[zone.fixed_point],
                min_separation=min_separation,
            )
            sizes.append(len(result.all_interior_bridges()))
            _repin(session, fixed_points, pips)
    _partition(session)
    blasted = []
    if inner_blasts:
        blasted.append(f"{inner_blasts} inner")
    if outer_blasts:
        blasted.append(f"{outer_blasts} outer")
    suffix = f", {' + '.join(blasted)} blast(s)" if blasted else ""
    return TangleBuild(
        session=session,
        fixed_points=fixed_points,
        pips=pips,
        stable_points=stable_points,
        title=f"Hénon k={HENON_P3[0]:g}, period 1 + period 3{suffix}",
        slug="nested_p1_p3",
        blast_sizes=sizes,
    )

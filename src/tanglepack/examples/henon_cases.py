"""Higher-period and nested Hénon tangles, built up to their topology.

Three recipes on the nested Hénon map ``HENON_P3 = (2, 1)``, shared by the
figure scripts and the test fixtures so the two can never drift apart:

* :func:`build_period3` -- the inner period-3 saddle orbit alone;
* :func:`build_nested` -- the outer period-1 saddle and the inner period-3
  orbit together, one resonance zone each;
* :func:`build_heteroclinic` -- :func:`build_nested` with period 3 grown far
  enough to cross the outer saddle's stable manifold.

Both end at ``session.topology()`` (holes, partition, classes, minimal
trellis, dual graph and words of every fixed point, from the pips the zones
were cut at).

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
  image sub-arc; since WP4 the signs decide the openings and every class
  walks). Past the reference switch (BCG-4) the word of an image-bridge class
  names the class of its image, a pair no computed bridge spans; since D11
  (WP8b) that pair is a class of its own, flagged "no computed bridge yet",
  and the runs are reliable. The default stays at 13, before the reference
  switch.
* Each nested zone blasts only its own fixed point (2026-10-02): blasting
  both from either zone also iterated the other tangle's bridges, made the
  result depend on blast order and collapsed the period-3 pseudoneighbors at
  6 or more outer blasts.
* The nested case's outer zone is blasted twice (the k=2.8 recipe): with no
  blast the outer class ``e`` is unresolved (a singleton landing with no
  registered chain); after two blasts every class resolves and the outer
  classes reproduce the k=2.8 two-blast words, with real (``unique``) walks.
* Construction order is OUTERMOST first (NEST-6, WP9): the period-1 saddle is
  constructed first, so it is ``A`` and tangle 0 of the session, matching
  ``TangleBuild.fixed_points``; the orbits are still grown period 3 first, each
  at its own ``area_cutoff``, so the session ends at period 1's.
* The heteroclinic recipe (:func:`build_heteroclinic`, WP9; scan of
  ``build_nested`` on 2026-10-08): period 3 at 13-15 unstable steps (any stable
  growth 9-12) has NO heteroclinic crossing; at 16 unstable steps it has four
  ``W^u(p3) x W^s(p1)`` crossings (0.6-0.8 s), two inside both pip-orbit cuts
  (the other two past p3's unstable cut: points of X). 16 / 10 was chosen while
  the law ``dual_face_side_is_geometry`` was a half-plane test that failed at
  16 / 9 and 16 / 12 (WP10 made it a local probe, which passes there); it is
  kept. 17 / 9 has 22 heteroclinic crossings but builds in 34 s; 13 / 15 has
  seven ``W^u(p1) x W^s(p3)`` crossings, two of them on one unstable segment,
  which cut a bridge past its own end until WP10 bracketed the cut (it builds
  now; its iterate links fail the beta and area laws, NUM-20). None of the four
  crossings has an iterate link (the preimage lies past p1's stable trim, the
  image past p3's unstable cut: escape-lobe images, WP10).
* The pips are the default strong pips, never pinned by id. The session
  carries them across every blast (API-1, WP8): a zone remembers the pip it
  was cut at, and each recompute finds it again by its cdists and keys.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

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
    One built tangle and its topology.

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


def _orbit(session, guess: str, orientation: dict[str, NDArray[np.float64]]):
    """Construct and orient one orbit of the nested map (``guess`` names a saddle guess)."""
    fp = session.construct_fixed_point(saddle_guesses(*HENON_P3)[guess])
    session.orient_eigenvectors(fp, orientation)
    return fp


def _grow(
    session, fp, unstable_steps: int, stable_steps: Optional[int], area_cutoff: float
) -> None:
    """Initialise and grow one orbit; ``stable_steps=None`` grows its stable side to turnaround."""
    session.area_cutoff = area_cutoff
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=unstable_steps)
    if stable_steps is None:
        session.grow_until_turnaround(fp, "stable")
    else:
        session.grow_n_times(fp, "stable", num_iterations=stable_steps)


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
        The :class:`TangleBuild` (one fixed point), its topology computed.
    """
    session = _session()
    fp3 = _orbit(session, "period_3", PERIOD3_ORIENTATION)
    _grow(session, fp3, unstable_steps, stable_steps, area_cutoff)
    _intersect_and_bridge(session, [fp3])
    stable_points = _stable_points(session, [fp3])

    session.classify_strong_pips()
    pip = session.trellis(fp3).strong_pip
    zone = session.resonance_zone(pip)
    sizes: list[int] = []
    for _ in range(blasts):
        result = session.blast_zone(
            zone, num_iterations=1, fixed_point=[fp3], min_separation=min_separation
        )
        sizes.append(len(result.all_interior_bridges()))
    session.topology()
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


def _nested_slug(outer_blasts: int, inner_blasts: int) -> str:
    """The figure folder of a nested build: the blast counts other than the
    default (two outer, no inner) are part of it, so runs never overwrite each other."""
    slug = "nested_p1_p3"
    if outer_blasts != 2:
        slug += f"_outer{outer_blasts}"
    if inner_blasts:
        slug += f"_inner{inner_blasts}"
    return slug


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
        first, the session's order too: ``fp1`` is ``A`` and tangle 0), its
        topology computed.
    """
    session = _session()
    fixed_points = [
        _orbit(session, "saddle", PERIOD1_ORIENTATION),
        _orbit(session, "period_3", PERIOD3_ORIENTATION),
    ]
    fp1, fp3 = fixed_points
    _grow(session, fp3, p3_unstable_steps, p3_stable_steps, p3_area_cutoff)
    _grow(session, fp1, p1_unstable_steps, None, p1_area_cutoff)
    _intersect_and_bridge(session, fixed_points)
    stable_points = _stable_points(session, fixed_points)

    session.classify_strong_pips()
    pips = [session.trellis(fp).strong_pip for fp in fixed_points]
    session.add_resonance_zones(pips)
    zones = sorted(session.resonance_zones.values(), key=lambda zone: zone.area)
    inner, outer = zones[0], zones[-1]
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
    session.topology()
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
        slug=_nested_slug(outer_blasts, inner_blasts),
        blast_sizes=sizes,
    )


def build_heteroclinic(**overrides) -> TangleBuild:
    """
    The nested tangles with heteroclinic crossings: period 3 grown to 16 / 10 steps.

    At 16 unstable steps the inner orbit's unstable manifold reaches across the
    outer saddle's stable manifold: four ``W^u(p3) x W^s(p1)`` crossings, two of
    them inside both pip-orbit cuts, a thin lobe across p1's fundamental
    segment (see the Dev Notes for why this recipe).

    Args:
        **overrides: Keywords of :func:`build_nested` replacing this recipe's.

    Returns:
        The :class:`TangleBuild`, fixed points ``[fp1, fp3]``, its topology
        computed, slug ``heteroclinic_p1_p3`` (plus the blast counts, as
        :func:`build_nested`'s).
    """
    build = build_nested(**{"p3_unstable_steps": 16, "p3_stable_steps": 10, **overrides})
    build.title = build.title.replace("period 1 + period 3", "period 1 + period 3, heteroclinic")
    build.slug = build.slug.replace("nested", "heteroclinic", 1)
    return build

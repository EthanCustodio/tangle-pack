"""Higher-period, nested and inversion Hénon tangles, built up to their topology.

Three recipes on the nested Hénon map ``HENON_P3 = (2, 1)``, shared by the
figure scripts and the test fixtures so the two can never drift apart:

* :func:`build_period3` -- the inner period-3 saddle orbit alone;
* :func:`build_nested` -- the outer period-1 saddle and the inner period-3
  orbit together, one resonance zone each;
* :func:`build_heteroclinic` -- :func:`build_nested` with period 3 grown far
  enough to cross the outer saddle's stable manifold.

and one on ``HENON_K4 = (4, 1)``:

* :func:`build_inversion` -- the saddle with both eigenvalues negative, its
  zone at ``f(q0)`` (the recipe of the inversion law cases; why ``f(q0)``:
  Dev Notes).

All end at ``session.topology()`` (holes, partition, classes, minimal
trellis, dual graph and words of every fixed point, from the pips the zones
were cut at).

Dev Notes:

* Parameters follow ``scripts/henon_blast_period_3.py`` and the
  ``henon_p3_session`` fixture: period 3 at ``area_cutoff = 1e-7`` with
  orientation unstable ``[0, -1]`` / stable ``[-1, -1]``, grown 13 unstable /
  9 stable steps; period 1 with unstable ``[-1, 0]`` / stable ``[0, 1]``,
  11 unstable steps and ``grow_until_turnaround`` on the stable side.
* Period 3 at 13 steps is closed up to 5 blasts: every active hole bridge's
  image is already a hole bridge on the next branch (``a -> b -> c -> a u^-1
  v^-1``), so the minimal trellis has no image bridges. At 6 blasts (or at 15
  unstable steps) deep crossings fold onto the fundamental window and the
  reference pair switches (a discovery, not a bug): 8 -> 42 holes, and the
  alphabet plateaus at 23 seed / 43 closed classes from 6 blasts (measured to
  20 blasts), every class walks and the runs are reliable. A pair no computed
  bridge spans is a class of its own, flagged "no computed bridge yet".
  The default stays at 13 steps, unblasted. Every period-3 build has its own
  figure folder (:func:`_period3_slug`).
* Each nested zone blasts only its own fixed point: blasting both from
  either zone would iterate the other tangle's bridges, make the result
  depend on blast order and collapse the period-3 pseudoneighbors at 6 or
  more outer blasts.
* The nested case's outer zone is blasted twice (the k=2.8 recipe), chosen to
  make the case closed: with no blast the outer class ``b`` is unresolved (a
  singleton landing with no registered chain); one blast resolves it through
  the trellis path; two blasts are the fewest at which every outer class
  walks, and the outer classes reproduce the k=2.8 two-blast words with real
  walks (``source == "walk"``).
* Construction order is OUTERMOST first: the period-1 saddle is
  constructed first, so it is ``A`` and tangle 0 of the session, matching
  ``TangleBuild.fixed_points``; the orbits are still grown period 3 first, each
  at its own ``area_cutoff``, so the session ends at period 1's.
* The heteroclinic recipe (:func:`build_heteroclinic`; a scan of
  ``build_nested``): period 3 at 13-15 unstable steps (any stable
  growth 9-12) has NO heteroclinic crossing; at 16 unstable steps it has four
  ``W^u(p3) x W^s(p1)`` crossings (0.6-0.8 s), all four inside p1's stable cut
  (two of them lay past p3's unstable cut until D12 was reversed, 2026-10-31).
  16 / 10 is kept: it passes the law ``dual_face_side_is_geometry`` in every
  form the law has had (a half-plane test failed at 16 / 9 and 16 / 12).
  17 / 9 has 22 heteroclinic crossings but builds in 34 s; 13 / 15 has
  seven ``W^u(p1) x W^s(p3)`` crossings, two of them on one unstable segment
  (it builds; its iterate links fail the beta and area laws). None of the four
  crossings has an iterate link (the preimage lies past p1's stable trim, the
  image is not registered: escape-lobe images).
* The inversion recipe is the law cases' (``tests/cases.build_inversion``):
  11 unstable and 11 stable steps at ``area_cutoff = 1e-5`` (step 12 escapes
  into the point budget), no ``orient_eigenvectors`` (the solver's
  eigenvectors pick branch 0), the zone at ``f(q0)``; the blasts' separation
  is 1e-5. Its class count drops at 6 blasts by a genuine reference switch
  (open).
* Why ``f(q0)`` and not ``q0`` (measured, not designed): on the k=2.8 and
  k=3.2 recipes (``tests/cases.py``) the zone at ``f(q0)`` after n blasts has
  the topology of the zone at ``q0`` after n-1 (holes, classes, words;
  n = 1..3), so ``f(q0)`` costs one blast and keeps the unblasted
  singleton-to-singleton case (the trellis path). On the inversion saddle this
  holds for n = 1 only: its zone at ``q0`` after one blast already has 15
  holes.
* The pips are the default strong pips, never pinned by id. The session
  carries them across every blast: a zone remembers the pip it
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

from .henon import HENON_K4, HENON_P3, henon_jacobian, henon_map, henon_map_inverse, saddle_guesses

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
        slug=_period3_slug(unstable_steps, stable_steps, blasts),
        blast_sizes=sizes,
    )


def _period3_slug(unstable_steps: int, stable_steps: int, blasts: int) -> str:
    """The figure folder of a period-3 build: growth and blast counts other than
    the default (13 / 9 steps, unblasted) are part of it."""
    slug = "p3"
    if unstable_steps != 13:
        slug += f"_u{unstable_steps}"
    if stable_steps != 9:
        slug += f"_s{stable_steps}"
    if blasts:
        slug += f"_blast{blasts}"
    return slug


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
    session.add_resonance_zones(pips, nesting={pips[1]: pips[0]})
    inner, outer = session.resonance_zones[fp3], session.resonance_zones[fp1]
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


def build_inversion(
    *,
    blasts: int = 0,
    unstable_steps: int = 11,
    stable_steps: int = 11,
    area_cutoff: float = 1e-5,
    min_separation: float = 1e-5,
) -> TangleBuild:
    """
    The ``k=4`` saddle with inversion, its zone at ``f(q0)`` blasted.

    Both eigenvalues of the saddle at ``(1 - sqrt5, sqrt5 - 1)`` are negative,
    so ``k_value = 2`` and the two branches of each manifold are one chain.

    Args:
        blasts: Single-iteration blasts of the zone.
        unstable_steps: Map steps to grow the unstable manifold.
        stable_steps: Map steps to grow the stable manifold.
        area_cutoff: Refinement area cutoff for the manifold machine.
        min_separation: The blasts' ``min_separation``.

    Returns:
        The :class:`TangleBuild` (one fixed point, the pip ``f(q0)``), its
        topology computed, slug ``inversion`` (plus ``_blast<n>``).

    Raises:
        ValueError: If the default strong pip has no registered image (too
            little growth: 8 unstable steps).
    """
    from ..loom.TangleSession import TangleSession

    session = TangleSession(
        henon_map(*HENON_K4),
        henon_map_inverse(*HENON_K4),
        henon_jacobian(*HENON_K4),
        area_cutoff=area_cutoff,
    )
    fp = session.construct_fixed_point(saddle_guesses(*HENON_K4)["inversion"])
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=unstable_steps)
    session.grow_n_times(fp, "stable", num_iterations=stable_steps)
    _intersect_and_bridge(session, [fp])
    stable_points = _stable_points(session, [fp])

    session.classify_strong_pips()
    trellis = session.trellis(fp)
    pip = trellis.iterate(trellis.strong_pip, 1)
    if pip is None:
        raise ValueError(
            f"the strong pip has no registered image at {unstable_steps} unstable steps, "
            "so the zone at f(q0) cannot be cut; grow deeper"
        )
    zone = session.resonance_zone(pip)
    sizes: list[int] = []
    for _ in range(blasts):
        result = session.blast_zone(
            zone, num_iterations=1, fixed_point=[fp], min_separation=min_separation
        )
        sizes.append(len(result.all_interior_bridges()))
    session.topology()
    blasted = f", {blasts} blast(s)" if blasts else ""
    return TangleBuild(
        session=session,
        fixed_points=[fp],
        pips=[pip],
        stable_points=stable_points,
        title=f"Hénon k={HENON_K4[0]:g}, inversion{blasted}",
        slug="inversion" + (f"_blast{blasts}" if blasts else ""),
        blast_sizes=sizes,
    )

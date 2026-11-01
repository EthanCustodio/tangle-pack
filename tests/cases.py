"""The tests-only case builders: every Hénon tangle the suite builds, frozen here.

Each builder returns a :class:`Case`: the session, its fixed points (outermost
first), the chosen strong pips, the resonance zones it made and the
:class:`CaseExpect` construction parameters (periods, inversion, tangles,
blasts) the law tier checks the build against. The builder
parameters are FROZEN in this file. They are never read from the defaults of
:mod:`tanglepack.examples.henon_cases` (those serve the figure scripts and may
move), nor from library defaults: every builder passes its ``area_cutoff``
explicitly (``CaseExpect.area_cutoff`` records it), so a change of a figure
recipe or of ``DEFAULT_AREA_CUTOFF`` can never silently change a test case.

Cases:

* ``k10`` -- the ``k=10, b=1`` binary horseshoe, the plain saddle at
  ``(4, -4)``: 9 unstable steps at ``area_cutoff = 1e-4``, the stable manifold
  grown to turnaround, no analytic Jacobian (the conftest recipe since the
  start of the suite).
* ``k28_one_blast`` / ``k28_two_blasts`` -- ``k=2.8``: 10 unstable steps at
  ``area_cutoff = 1e-7``, the stable manifold trimmed at the pip ``f(q0)`` (a
  resonance zone, which the session remembers as the pip), then one or two
  single-iteration blasts of that zone with ``min_separation = 1e-5``
  (:func:`_zone_at_image_pip`).
* ``p3`` -- the inner period-3 orbit of the nested map ``(2, 1)`` alone,
  delegated to :func:`tanglepack.examples.henon_cases.build_period3` with every
  keyword explicit.
* ``nested`` / ``nested_inner_blasts`` -- period 1 + period 3 with two outer
  blasts, and no or two inner blasts, delegated to
  :func:`tanglepack.examples.henon_cases.build_nested` with every keyword
  explicit (``p1_area_cutoff = 1e-7``).
* ``heteroclinic`` -- the nested recipe with period 3 grown 16 unstable / 10
  stable steps (two outer blasts, ``1e-7`` on both): the inner orbit's unstable
  manifold reaches across the outer saddle's stable manifold, four ``A(u) x
  B(s)`` crossings, two inside both pip-orbit cuts (NEST-1, WP9). 16 / 10 was
  chosen while ``dual_face_side_is_geometry`` was a half-plane test that failed
  at 16 / 9 and 16 / 12 (WP10 made it a local probe); the case keeps 16 / 10.
  ``heteroclinic_17u`` grows period 3 one unstable step further, where the
  images of the connecting bridges are computed (a perf regression, WP9c).
* ``nested_unblasted`` -- the nested map with no blast and the outer saddle
  grown at ``p1_area_cutoff = 1e-4`` (the old ``henon_p3_session`` recipe, kept
  as a parameter: Phase 0 probe P6 showed the cutoff changes the outer tangle).
* ``inversion`` / ``inversion_one_blast`` / ``inversion_two_blasts`` -- the
  ``k=4`` map's OTHER saddle at ``(1 - sqrt5, sqrt5 - 1)``: both eigenvalues
  negative (-1.96 / -0.51), ``k_value = 2 * period``, det J = +1 (inversion,
  NOT orientation reversing); 11 unstable / 11 stable steps at ``area_cutoff
  = 1e-5`` (step 12 escapes), then the k=2.8 zone recipe with 0, 1 or 2
  blasts (INV-4).
* ``inversion_period2`` -- inversion on a period-k orbit: the ``k=4.1``
  period-2 orbit ``(-2.0488, -0.0488) <-> (0.0488, 2.0488)``, both ``DM^2``
  eigenvalues negative (-1.863 / -0.537), ``k_value = 4``; 9 unstable / 10
  stable steps at ``area_cutoff = 1e-5`` (Q12), no zone (at 10 / 10 its stable
  manifold runs past the pip orbit; with the trellis cut there (D10, WP8b) the
  only law failing is ``bridge_endpoints_on_own_branch`` on a WORKBENCH bridge,
  WP7-D1, so it is a regression of its own, ``test_pip_orbit_cut``; it escapes
  past ~8 steps at k=4.25).
* ``k10_inversion`` -- the ``k=10`` map's inversion saddle at ``(-2.3166,
  2.3166)``, 6 unstable / 5 stable steps at ``area_cutoff = 1e-4``: the
  numerics fixture behind ``henon_inversion*`` only (it escapes at 7 steps and
  has no holes, so it is no law case).
* ``orientation_reversing`` -- the ``b = -1`` placeholder (``k = 4``, seed
  ``[2, 2]``, ``area_cutoff = 1e-4``): ``construct_fixed_point`` raises
  ``ValueError`` because a single ``k_value`` cannot describe a saddle with
  exactly one negative eigenvalue; det J < 0 is out of scope (C3).
* ``p3_15`` / ``p3_16`` / ``p3_6_blasts`` -- the deep period-3 runs past the
  reference switch (BCG-4): their words name classes outside the minimal
  trellis's alphabet, which the symbolic dynamics closes under the word map
  (PK-1, WP7).
* ``p3_16_12`` / ``inversion_period2_10`` -- period 3 at 16 unstable / 12
  stable steps, and the inversion x period-2 orbit at 10 / 10 steps refined to
  ``area_cutoff = 1e-6`` (at 1e-5 two crossings of its reference window are
  below the resolution threshold, outside N, B1): the deep runs every
  non-mutating law holds on since WP10 (two crossings on one segment used to
  cut a bridge past its own end there; two laws over-claimed on p3).
* ``near_tangency`` -- the period-1 ``b = 1`` saddle at ``k = -0.999``, close to
  the saddle-centre bifurcation at ``k = -1`` where the splitting angle is
  exponentially small, each arm grown to 20 fixed-point separations of
  canonical distance at its own ``area_cutoff = 1e-10``, through
  ``compute_intersections`` only. Its crossings are near-tangent polylines
  straddling each other: the resolution guards' positive case. NOT a law case.
* ``k28_deep_blast_sliver`` -- ``k28`` from the analytic seed, blasted nine
  times: a crossing sits a sliver narrower than ``2 * cdist_tol`` past the pip
  (stress campaign E11, 2026-10-07; the exact-ownership regression). The blasts
  run with clause 3 of N off (the near-anchor guard, 2026-10-29, ends the pip
  orbit before the ninth blast makes the sliver) and with no ``min_separation``
  (the sliver's bridge ends on the pip orbit; since that child lost its
  exemption, 2026-10-31, the 1e-5 guard drops it).

The blast series (:data:`SERIES`, :func:`blast_series`) run the D1 loop --
blast once, recompute the topology from scratch, compare -- from the zone
stage of ``k28`` (8 blasts), ``p3`` (9), ``inversion`` (4), E11's ``k32`` (10;
``build_k32_zone``, the ``k28`` recipe at ``k = 3.2`` from its analytic seed)
and E11's ``nested_outer`` (11 blasts of the outer zone at ``min_separation =
1e-4``).

Dev Notes:

* ``henon_cases`` has no stage parameter (and ``src/`` is not changed by the
  test refactor), so a nested or period-3 build stopped before the partition
  (``through="zones"``) is a local copy of the ``henon_cases`` recipe; only the
  full ``"partitioned"`` build delegates.
* Registry ids are reproducible since WP8 (crossings are registered along
  each unstable branch in cdist order, INT-7; before, three ``build_k10()``
  calls in one process gave three permutations). Tests still never pin an id
  and never pick "the lowest id" as a canonical choice.
* The session carries the strong pips across every blast (API-1, WP8): a
  zone remembers the pip it was cut at, so no builder restores pips by hand,
  and ``session.topology()`` recomputes the holes and the partition from them.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from types import MappingProxyType
from typing import TYPE_CHECKING, Callable, Literal, Mapping, NamedTuple, Optional

import numpy as np
import pytest

from tanglepack import TangleSession
from tanglepack.examples import (
    HENON_K10,
    HENON_P3,
    henon_jacobian,
    henon_map,
    henon_map_inverse,
)
from tanglepack.examples import henon_cases

from helpers.direct import trellis_view
from helpers.laws import base_chord_misreadings, p1_itinerary_disagreements, reference_witness

if TYPE_CHECKING:  # pragma: no cover - typing only
    from tanglepack import FixedPoint, IntersectionRegistry, TangleWorkbench
    from tanglepack.loom.Blast import BlastResult
    from tanglepack.loom.ResonanceZone import ResonanceZone

#: How far a builder runs. Each stage includes every earlier one:
#: ``empty`` (a session, no fixed point), ``fixed_point`` (constructed and
#: oriented), ``seeded`` (fundamental segments), ``grown_unstable`` (unstable
#: growth only), ``grown`` (both manifolds), ``intersected``
#: (``compute_intersections``), ``bridges`` (trim, bridges, iterate table),
#: ``pips`` (strong pips classified), ``zones`` (resonance zones, blasts),
#: ``partitioned`` (``session.topology()``: pseudoneighbors, holes, stable
#: partition and every later stage).
Stage = Literal[
    "empty",
    "fixed_point",
    "seeded",
    "grown_unstable",
    "grown",
    "intersected",
    "bridges",
    "pips",
    "zones",
    "partitioned",
]

_STAGES: tuple[str, ...] = (
    "empty",
    "fixed_point",
    "seeded",
    "grown_unstable",
    "grown",
    "intersected",
    "bridges",
    "pips",
    "zones",
    "partitioned",
)

# --------------------------------------------------------------------------- #
# Frozen parameters
# --------------------------------------------------------------------------- #
HENON_K28: tuple[float, float] = (2.8, 1)
#: The second single saddle of the stress campaign's E11 (``k32_series``).
HENON_K32: tuple[float, float] = (3.2, 1)
#: The map of the inversion law cases (its OTHER saddle, ``K4_INVERSION_SEED``).
HENON_K4: tuple[float, float] = (4.0, 1)
#: The map of the inversion x period-k case (``K41_INVERSION_PERIOD2_SEED``).
HENON_K41: tuple[float, float] = (4.1, 1)
#: The ``b = -1`` orientation-reversing placeholder map (P4).
HENON_ORIENTATION_REVERSING: tuple[float, float] = (4, -1)
#: Close to the saddle-centre bifurcation at k = -1 (the near-tangency case).
HENON_NEAR_TANGENCY: tuple[float, float] = (-0.999, 1)

K10_SADDLE_SEED: tuple[float, float] = (4.0, -4.0)
K10_INVERSION_SEED: tuple[float, float] = (-2.3166, 2.3166)
#: The k=4 OTHER root ``(1 - sqrt5, sqrt5 - 1)``: eigenvalues -1.96 / -0.51 (INV-4).
K4_INVERSION_SEED: tuple[float, float] = (1 - np.sqrt(5), np.sqrt(5) - 1)
#: The ``k=4.1`` period-2 orbit whose ``DM^2`` eigenvalues are both negative
#: (-1.863 / -0.537): inversion on a period-k orbit, ``k_value = 4`` (WP3).
K41_INVERSION_PERIOD2_SEED: tuple[tuple[float, float], ...] = (
    (-2.0488, -0.0488),
    (0.0488, 2.0488),
)
K28_SADDLE_SEED: tuple[float, float] = (4.0, -4.0)
#: The analytic ``k=2.8`` saddle ``x = (2 + sqrt(4 + 4k)) / 2`` (the stress
#: campaign's seed). It converges to the same saddle as ``K28_SADDLE_SEED`` but
#: not to the same last bits, and nine blasts of it leave a real crossing
#: 1.1e-6 past the pip, a sliver whose midpoint lies within ``cdist_tol`` of the
#: pip (the E11 sliver).
K28_ANALYTIC_SADDLE_SEED: tuple[float, float] = (2.9493588689617924, -2.9493588689617924)
#: The analytic ``k=3.2`` saddle ``x = (2 + sqrt(4 + 4k)) / 2`` (E11's seed).
K32_SADDLE_SEED: tuple[float, float] = (3.04939015319192, -3.04939015319192)
ORIENTATION_REVERSING_SEED: tuple[float, float] = (2.0, 2.0)
#: The ``min_separation`` of every single-iteration zone blast (k28, inversion, series).
ZONE_MIN_SEPARATION: float = 1e-5
P3_ORBIT_SEED: tuple[tuple[float, float], ...] = ((0, 1), (-1, 0), (-1, 1))
P3_OUTER_SADDLE_SEED: tuple[float, float] = (4.0, -4.0)

#: Eigenvector orientation hints of every period-1 saddle above.
PERIOD1_ORIENTATION: Mapping[str, tuple[int, int]] = MappingProxyType(
    {"unstable": (-1, 0), "stable": (0, 1)}
)
#: Eigenvector orientation hints of the period-3 orbit.
PERIOD3_ORIENTATION: Mapping[str, tuple[int, int]] = MappingProxyType(
    {"unstable": (0, -1), "stable": (-1, -1)}
)


@dataclass(frozen=True)
class CaseExpect:
    """
    What the law tier may assume about a case.

    Attributes:
        periods: The period of each fixed point, outermost first.
        area_cutoff: The machine's refinement cutoff the build ends with (the
            outer saddle's for nested), always passed explicitly.
        inversion: Whether any fixed point has inversion (``k_value == 2 * period``).
        n_tangles: The number of fixed points.
        blasts: Single-iteration blasts applied (outer blasts for nested).
        heteroclinic: Whether the selection trellis holds a crossing of one
            fixed point's unstable and another's stable branch.
    """

    periods: tuple[int, ...]
    area_cutoff: float
    inversion: bool = False
    n_tangles: int = 1
    blasts: int = 0
    heteroclinic: bool = False


@dataclass
class Case:
    """
    One built Hénon tangle.

    Attributes:
        name: The case name (a key of :data:`BUILDERS`, or a variant label).
        session: The session that owns the workbench.
        fixed_points: The fixed points, outermost first.
        pips: The strong pip chosen per fixed point (None before ``pips``).
        expect: The case expectations.
        zones: The resonance zones made by the build, in creation order.
        stage: The stage the build stopped at.
    """

    name: str
    session: TangleSession
    fixed_points: list["FixedPoint"]
    pips: list[Optional[int]]
    expect: CaseExpect
    zones: list["ResonanceZone"] = field(default_factory=list)
    stage: str = "partitioned"

    @property
    def workbench(self) -> "TangleWorkbench":
        """The session's workbench."""
        return self.session.workbench

    @property
    def registry(self) -> "IntersectionRegistry":
        """The workbench's intersection registry."""
        return self.session.workbench.intersection_registry

    @property
    def fixed_point(self) -> "FixedPoint":
        """The single fixed point of a one-tangle case.

        Raises:
            ValueError: The case has more than one fixed point.
        """
        if len(self.fixed_points) != 1:
            raise ValueError(f"case {self.name!r} has {len(self.fixed_points)} fixed points")
        return self.fixed_points[0]


def _reached(through: str, stage: str) -> bool:
    """True when a build that stops at ``through`` has done ``stage``."""
    if through not in _STAGES:
        raise ValueError(f"unknown stage {through!r}; known stages are {_STAGES}")
    return _STAGES.index(through) >= _STAGES.index(stage)


def _orientation(hints: Mapping[str, tuple[int, int]]) -> dict[str, np.ndarray]:
    """The orientation hints as the arrays ``orient_eigenvectors`` takes."""
    return {stability: np.array(vector) for stability, vector in hints.items()}


def _blast(
    session: TangleSession,
    zone: "ResonanceZone",
    fixed_points: list,
    *,
    diagnose: bool = True,
    min_separation: float = ZONE_MIN_SEPARATION,
    representatives: Optional[str] = None,
) -> "BlastResult":
    """One single-iteration blast of ``zone`` (at ``ZONE_MIN_SEPARATION`` by default)."""
    return session.blast_zone(
        zone,
        num_iterations=1,
        fixed_point=fixed_points,
        min_separation=min_separation,
        diagnose=diagnose,
        representatives=representatives,
    )


def _zone_at_image_pip(
    session: TangleSession, fp: "FixedPoint", blasts: int
) -> tuple[int, "ResonanceZone"]:
    """The k=2.8 zone recipe: the zone trimmed at the pip ``f(q0)``, blasted ``blasts`` times.

    Returns:
        ``(pip, zone)``; the session remembers the pip the zone was cut at.
    """
    session.classify_strong_pips()
    trellis = session.trellis(fp)
    pip = trellis.iterate(trellis.strong_pip, 1)
    assert pip is not None, "the default strong pip must have a registered image"
    zone = session.resonance_zone(pip)
    for _ in range(blasts):
        _blast(session, zone, [fp])
    return pip, zone


# --------------------------------------------------------------------------- #
# k = 10
# --------------------------------------------------------------------------- #
def build_k10(
    *,
    unstable_steps: int = 9,
    through: Stage = "partitioned",
) -> Case:
    """
    The ``k=10`` binary-horseshoe saddle at ``(4, -4)``.

    Args:
        unstable_steps: Unstable growth steps (9; the numerics stage fixtures
            use 7).
        through: The last stage to run (see :data:`Stage`). ``zones`` adds
            nothing for this case.

    Returns:
        The case. No analytic Jacobian is passed (the historical recipe).
    """
    area_cutoff = 1e-4
    session = TangleSession(
        henon_map(*HENON_K10), henon_map_inverse(*HENON_K10), area_cutoff=area_cutoff
    )
    case = Case(
        name="k10",
        session=session,
        fixed_points=[],
        pips=[None],
        expect=CaseExpect(periods=(1,), area_cutoff=area_cutoff),
        stage=through,
    )
    if not _reached(through, "fixed_point"):
        return case
    fp = session.construct_fixed_point(list(K10_SADDLE_SEED))
    session.orient_eigenvectors(fp, _orientation(PERIOD1_ORIENTATION))
    case.fixed_points = [fp]
    if not _reached(through, "seeded"):
        return case
    session.initialize_both_manifolds(fp)
    if not _reached(through, "grown_unstable"):
        return case
    session.grow_n_times(fp, "unstable", num_iterations=unstable_steps)
    if not _reached(through, "grown"):
        return case
    session.grow_until_turnaround(fp, "stable")
    if not _reached(through, "intersected"):
        return case
    session.compute_intersections([fp])
    if not _reached(through, "bridges"):
        return case
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    if not _reached(through, "pips"):
        return case
    session.classify_strong_pips()
    case.pips = [session.strong_pip(fp)]
    if not _reached(through, "partitioned"):
        return case
    session.topology()
    return case


# --------------------------------------------------------------------------- #
# k = 2.8
# --------------------------------------------------------------------------- #
def _saddle_zone(
    henon: tuple[float, float], seed: tuple[float, float], blasts: int
) -> tuple[TangleSession, "FixedPoint", int, "ResonanceZone"]:
    """The period-1 zone recipe of ``k28`` and ``k32``: 10 unstable steps at
    ``area_cutoff = 1e-7``, the stable manifold to turnaround, the zone at
    ``f(q0)`` blasted ``blasts`` times (:func:`_zone_at_image_pip`)."""
    session = TangleSession(
        henon_map(*henon), henon_map_inverse(*henon), henon_jacobian(*henon), area_cutoff=1e-7
    )
    fp = session.construct_fixed_point(list(seed))
    session.orient_eigenvectors(fp, _orientation(PERIOD1_ORIENTATION))
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=10)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    pip, zone = _zone_at_image_pip(session, fp, blasts)
    return session, fp, pip, zone


def build_k28(
    *,
    blasts: int = 1,
    through: Stage = "partitioned",
    seed: tuple[float, float] = K28_SADDLE_SEED,
) -> Case:
    """
    The ``k=2.8`` saddle, its zone trimmed at ``f(q0)`` and blasted.

    Args:
        blasts: Single-iteration blasts of the zone (1 or 2 for the law cases).
        through: ``"zones"`` (trimmed and blasted) or ``"partitioned"``.
        seed: The saddle's solver seed (the deep-blast regression uses
            ``K28_ANALYTIC_SADDLE_SEED``).

    Returns:
        The case, the pip ``f(q0)`` carried by the session through every blast.

    Raises:
        ValueError: ``through`` is neither ``"zones"`` nor ``"partitioned"``.
    """
    if through not in ("zones", "partitioned"):
        raise ValueError("build_k28 stops at 'zones' or 'partitioned'")
    session, fp, pip, zone = _saddle_zone(HENON_K28, seed, blasts)
    if through == "partitioned":
        session.topology()
    return Case(
        name={0: "k28_unblasted", 1: "k28_one_blast", 2: "k28_two_blasts"}.get(
            blasts, f"k28_{blasts}_blasts"
        ),
        session=session,
        fixed_points=[fp],
        pips=[pip],
        expect=CaseExpect(periods=(1,), area_cutoff=1e-7, blasts=blasts),
        zones=[zone],
        stage=through,
    )


def build_k28_deep_blast_sliver() -> Case:
    """
    ``k28`` from the analytic seed, blasted nine times with clause 3 of N off
    and no ``min_separation``.

    The sliver past the pip (stress campaign E11) is made by the ninth blast,
    which the near-anchor guard (2026-10-29) never reaches: it ends the pip
    orbit first. The 1e-5 proximity guard would drop its bridge, which ends
    on the pip orbit, since the pip-orbit exemption went (2026-10-31). The
    build is :func:`build_k28` otherwise.

    Returns:
        The partitioned case.
    """
    case = build_k28(blasts=0, through="zones", seed=K28_ANALYTIC_SADDLE_SEED)
    session, (fp,), (zone,) = case.session, case.fixed_points, case.zones
    session.anchor_resolution_factor = None
    for _ in range(9):
        _blast(session, zone, [fp], min_separation=None)
    session.topology()
    case.name, case.stage = "k28_9_blasts", "partitioned"
    case.expect = CaseExpect(periods=(1,), area_cutoff=1e-7, blasts=9)
    return case


def build_k32_zone() -> Case:
    """
    The ``k=3.2`` saddle by the ``k28`` recipe, from its analytic seed (E11's
    ``k32``), stopped at its zone ``f(q0)``: the start of ``k32_series``.

    Returns:
        The case at the ``zones`` stage, unblasted.
    """
    session, fp, pip, zone = _saddle_zone(HENON_K32, K32_SADDLE_SEED, 0)
    return Case(
        name="k32_unblasted",
        session=session,
        fixed_points=[fp],
        pips=[pip],
        expect=CaseExpect(periods=(1,), area_cutoff=1e-7),
        zones=[zone],
        stage="zones",
    )


# --------------------------------------------------------------------------- #
# The nested map (2, 1): period 3 alone, and period 1 + period 3
# --------------------------------------------------------------------------- #
def build_period3(
    *,
    unstable_steps: int = 13,
    stable_steps: int = 9,
    area_cutoff: float = 1e-7,
    blasts: int = 0,
    min_separation: float = 1e-5,
    through: Literal["zones", "partitioned"] = "partitioned",
) -> Case:
    """
    The inner period-3 orbit alone.

    ``"partitioned"`` delegates to
    :func:`tanglepack.examples.henon_cases.build_period3` with every keyword
    explicit; ``"zones"`` runs the local copy of its recipe and stops at the
    zone, unblasted (the blast series starts there).

    Args:
        unstable_steps: Unstable growth steps (13; 15 and 16 are the known
            open deep runs).
        stable_steps: Stable growth steps.
        area_cutoff: The manifold machine's refinement area cutoff.
        blasts: Single-iteration blasts of the zone.
        min_separation: The blasts' ``min_separation``.
        through: ``"zones"`` or ``"partitioned"``.

    Returns:
        The case.

    Raises:
        ValueError: ``"zones"`` was asked for together with blasts.
    """
    if through == "zones":
        if blasts:
            raise ValueError("a period-3 build stopped at its zone takes no blasts")
        return _period3_zone(unstable_steps, stable_steps, area_cutoff)
    build = henon_cases.build_period3(
        unstable_steps=unstable_steps,
        stable_steps=stable_steps,
        area_cutoff=area_cutoff,
        blasts=blasts,
        min_separation=min_separation,
    )
    name = "p3"
    if unstable_steps != 13:
        name += f"_{unstable_steps}"
    if blasts:
        name += f"_{blasts}_blasts"
    return Case(
        name=name,
        session=build.session,
        fixed_points=list(build.fixed_points),
        pips=list(build.pips),
        expect=CaseExpect(periods=(3,), area_cutoff=area_cutoff, blasts=blasts),
        zones=list(build.session.resonance_zones.values()),
    )


def _period3_zone(unstable_steps: int, stable_steps: int, area_cutoff: float) -> Case:
    """The ``henon_cases.build_period3`` recipe stopped at the zone."""
    session = TangleSession(
        henon_map(*HENON_P3),
        henon_map_inverse(*HENON_P3),
        henon_jacobian(*HENON_P3),
        area_cutoff=area_cutoff,
    )
    fp3 = session.construct_fixed_point([list(p) for p in P3_ORBIT_SEED])
    session.orient_eigenvectors(fp3, _orientation(PERIOD3_ORIENTATION))
    session.initialize_both_manifolds(fp3)
    session.grow_n_times(fp3, "unstable", num_iterations=unstable_steps)
    session.grow_n_times(fp3, "stable", num_iterations=stable_steps)
    session.compute_intersections([fp3])
    session.trim_stable_manifolds(fp3)
    session.create_bridges(fp3)
    session.classify_strong_pips()
    pip = session.trellis(fp3).strong_pip
    zone = session.resonance_zone(pip)
    return Case(
        name="p3_unblasted",
        session=session,
        fixed_points=[fp3],
        pips=[pip],
        expect=CaseExpect(periods=(3,), area_cutoff=area_cutoff),
        zones=[zone],
        stage="zones",
    )


def build_period3_coarse_blasts(*, blasts: int = 6, area_cutoff: float = 1e-5) -> Case:
    """
    The period-3 zone blasted at a coarse refinement cutoff (stress campaign E11).

    The zone is built at 1e-7 (:func:`build_period3` ``through="zones"``), then
    ``blasts`` single-iteration blasts refine at ``area_cutoff`` with
    ``min_separation = 1e-5`` (no same-stability pass), and the topology is
    computed. At 1e-5 and 6 blasts some hole
    bridges hold 2-4 nodes, where a polyline reading of their rows disagrees with
    the crossing signs (NUM-1). Their ends have resolution ratios near 1, so the
    blasts run with clause 1 of N off (``resolution_ratio_threshold = None``):
    under the default threshold those detections are outside N (Q13 a, B1) and
    no hole bridge is that coarse. Clause 3 is off too (2026-10-29): at 1e-5 the
    near-anchor guard stops the coarse blasts before they make such a bridge.

    Args:
        blasts: Single-iteration blasts of the zone.
        area_cutoff: The manifold machine's refinement cutoff during the blasts.

    Returns:
        The partitioned case.
    """
    case = build_period3(through="zones")
    session = case.session
    session.area_cutoff = area_cutoff
    session.workbench.resolution_ratio_threshold = None
    session.anchor_resolution_factor = None
    (zone,) = case.zones
    for _ in range(blasts):
        session.blast_zone(
            zone, num_iterations=1, fixed_point=case.fixed_points,
            min_separation=1e-5, diagnose=False,
        )
    session.topology()
    case.name = "p3_coarse_blasts"
    case.expect = CaseExpect(periods=(3,), area_cutoff=area_cutoff, blasts=blasts)
    case.stage = "partitioned"
    return case


def _nested_through(
    *,
    through: Stage,
    p3_unstable_steps: int,
    p3_stable_steps: int,
    p3_area_cutoff: float,
    p1_unstable_steps: int,
    p1_area_cutoff: float,
) -> Case:
    """The ``henon_cases.build_nested`` recipe stopped before the blasts."""
    session = TangleSession(
        henon_map(*HENON_P3), henon_map_inverse(*HENON_P3), henon_jacobian(*HENON_P3)
    )
    case = Case(
        name="nested_unblasted",
        session=session,
        fixed_points=[],
        pips=[None, None],
        expect=CaseExpect(periods=(1, 3), area_cutoff=p1_area_cutoff, n_tangles=2),
        stage=through,
    )
    # Constructed outermost first (NEST-6: the session's tangle order is the
    # case's), grown period 3 first, each at its own cutoff.
    fp1 = session.construct_fixed_point(list(P3_OUTER_SADDLE_SEED))
    session.orient_eigenvectors(fp1, _orientation(PERIOD1_ORIENTATION))
    fp3 = session.construct_fixed_point([list(p) for p in P3_ORBIT_SEED])
    session.orient_eigenvectors(fp3, _orientation(PERIOD3_ORIENTATION))
    case.fixed_points = [fp1, fp3]

    session.area_cutoff = p3_area_cutoff
    session.initialize_both_manifolds(fp3)
    session.grow_n_times(fp3, "unstable", num_iterations=p3_unstable_steps)
    session.grow_n_times(fp3, "stable", num_iterations=p3_stable_steps)

    session.area_cutoff = p1_area_cutoff
    session.initialize_both_manifolds(fp1)
    session.grow_n_times(fp1, "unstable", num_iterations=p1_unstable_steps)
    session.grow_until_turnaround(fp1, "stable")
    if not _reached(through, "intersected"):
        return case

    session.compute_intersections(case.fixed_points)
    if not _reached(through, "bridges"):
        return case
    for fp in case.fixed_points:
        session.trim_stable_manifolds(fp)
    for fp in case.fixed_points:
        session.create_bridges(fp)
    if not _reached(through, "pips"):
        return case

    # Per trellis, outermost first (the historical henon_p3_session order, which
    # also fixes the order of the session's trellis cache), not the fan-out.
    for fp in case.fixed_points:
        session.trellis(fp).classify_strong_pips()
    case.pips = [session.trellis(fp).strong_pip for fp in case.fixed_points]
    if not _reached(through, "zones"):
        return case
    session.add_resonance_zones(case.pips)
    case.zones = list(session.resonance_zones.values())
    return case


def build_nested(
    *,
    outer_blasts: int = 2,
    inner_blasts: int = 0,
    p3_unstable_steps: int = 13,
    p3_stable_steps: int = 9,
    p3_area_cutoff: float = 1e-7,
    p1_unstable_steps: int = 11,
    p1_area_cutoff: float = 1e-7,
    min_separation: float = 1e-4,
    through: Stage = "partitioned",
) -> Case:
    """
    The outer period-1 saddle and the inner period-3 orbit together.

    Args:
        outer_blasts: Single-iteration blasts of the outer (larger) zone.
        inner_blasts: Single-iteration blasts of the inner (smaller) zone.
        p3_unstable_steps: Unstable growth of the period-3 orbit.
        p3_stable_steps: Stable growth of the period-3 orbit.
        p3_area_cutoff: Area cutoff while growing the period-3 manifolds.
        p1_unstable_steps: Unstable growth of the period-1 saddle.
        p1_area_cutoff: Area cutoff while growing the period-1 manifolds (1e-7;
            the unblasted fixture uses 1e-4, which changes the outer tangle).
        min_separation: The blasts' ``min_separation``.
        through: ``"partitioned"`` delegates to
            :func:`tanglepack.examples.henon_cases.build_nested`; an earlier
            stage (at most ``"zones"``, and only without blasts) runs the
            local copy of its recipe and stops there.

    Returns:
        The case, fixed points ``[fp1, fp3]``.

    Raises:
        ValueError: An early stage was asked for together with blasts.
    """
    if through != "partitioned":
        if outer_blasts or inner_blasts:
            raise ValueError("a nested build stopped before the partition takes no blasts")
        if not _reached(through, "grown"):
            raise ValueError("a nested build starts at the 'grown' stage")
        return _nested_through(
            through=through,
            p3_unstable_steps=p3_unstable_steps,
            p3_stable_steps=p3_stable_steps,
            p3_area_cutoff=p3_area_cutoff,
            p1_unstable_steps=p1_unstable_steps,
            p1_area_cutoff=p1_area_cutoff,
        )
    build = henon_cases.build_nested(
        p3_unstable_steps=p3_unstable_steps,
        p3_stable_steps=p3_stable_steps,
        p3_area_cutoff=p3_area_cutoff,
        p1_unstable_steps=p1_unstable_steps,
        p1_area_cutoff=p1_area_cutoff,
        outer_blasts=outer_blasts,
        inner_blasts=inner_blasts,
        min_separation=min_separation,
    )
    name = {(2, 0): "nested", (2, 2): "nested_inner_blasts"}.get(
        (outer_blasts, inner_blasts), f"nested_{outer_blasts}_{inner_blasts}"
    )
    return Case(
        name=name,
        session=build.session,
        fixed_points=list(build.fixed_points),
        pips=list(build.pips),
        expect=CaseExpect(
            periods=(1, 3),
            area_cutoff=p1_area_cutoff,
            n_tangles=2,
            blasts=outer_blasts,
        ),
        zones=list(build.session.resonance_zones.values()),
    )


def build_heteroclinic(*, p3_unstable_steps: int = 16) -> Case:
    """
    Period 1 + period 3 with heteroclinic crossings: period 3 at 16 / 10 steps.

    Args:
        p3_unstable_steps: Unstable growth of the period-3 orbit (at 17 the
            images of the connecting bridges are computed too).

    Returns:
        The case ``"heteroclinic"``, fixed points ``[fp1, fp3]``.
    """
    case = build_nested(
        outer_blasts=2,
        inner_blasts=0,
        p3_unstable_steps=p3_unstable_steps,
        p3_stable_steps=10,
        p3_area_cutoff=1e-7,
        p1_unstable_steps=11,
        p1_area_cutoff=1e-7,
        min_separation=1e-4,
    )
    case.name = "heteroclinic"
    case.expect = replace(case.expect, heteroclinic=True)
    return case


# --------------------------------------------------------------------------- #
# Inversion: the k = 4 law cases and the k = 10 numerics fixture
# --------------------------------------------------------------------------- #
def build_inversion(*, blasts: int = 0, through: Stage = "partitioned") -> Case:
    """
    The ``k=4`` map's OTHER saddle at ``(1 - sqrt5, sqrt5 - 1)``, its zone at ``f(q0)`` blasted.

    Both eigenvalues are negative (-1.96 and -0.51), so ``k_value = 2 *
    period`` and the two branches of each manifold are one chain; det J = +1,
    so the map is orientation preserving. ``orient_eigenvectors`` is not called
    (the solver's eigenvectors pick branch 0). 11 unstable and 11 stable steps
    at ``area_cutoff = 1e-5`` (step 12 escapes), then the k=2.8 zone recipe
    (:func:`_zone_at_image_pip`).

    Args:
        blasts: Single-iteration blasts of the zone (0, 1 or 2 for the law cases).
        through: ``"zones"`` (the blast series starts there) or ``"partitioned"``.

    Returns:
        The case, the pip ``f(q0)`` carried by the session through every blast.

    Raises:
        ValueError: ``through`` is neither ``"zones"`` nor ``"partitioned"``.
    """
    if through not in ("zones", "partitioned"):
        raise ValueError("build_inversion stops at 'zones' or 'partitioned'")
    area_cutoff = 1e-5
    session = TangleSession(
        henon_map(*HENON_K4),
        henon_map_inverse(*HENON_K4),
        henon_jacobian(*HENON_K4),
        area_cutoff=area_cutoff,
    )
    fp = session.construct_fixed_point(list(K4_INVERSION_SEED))
    assert fp.k_value == 2 * fp.period, "the k=4 other root must have inversion"
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=11)
    session.grow_n_times(fp, "stable", num_iterations=11)
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    pip, zone = _zone_at_image_pip(session, fp, blasts)
    if through == "partitioned":
        session.topology()
    return Case(
        name={0: "inversion", 1: "inversion_one_blast", 2: "inversion_two_blasts"}.get(
            blasts, f"inversion_{blasts}_blasts"
        ),
        session=session,
        fixed_points=[fp],
        pips=[pip],
        expect=CaseExpect(periods=(1,), area_cutoff=area_cutoff, inversion=True, blasts=blasts),
        zones=[zone],
        stage=through,
    )


def build_inversion_period2(
    *,
    unstable_steps: int = 9,
    stable_steps: int = 10,
    area_cutoff: float = 1e-5,
    through: Stage = "partitioned",
) -> Case:
    """
    Inversion on a period-k orbit: the ``k=4.1`` period-2 orbit, no zone.

    Both eigenvalues of ``DM^2`` are negative (-1.863 and -0.537), so
    ``k_value = 2 * period = 4`` and each manifold's four pieces are one chain
    over the two orbit points and both branches; det J = +1. 9 unstable and 10
    stable steps at ``area_cutoff = 1e-5``, then the stable manifold trimmed at
    its outermost crossings and partitioned at the default strong pip. The
    tenth stable step resolves the class a 9-step build leaves unresolved (a
    singleton landing with no registered member image; author Q12).

    Args:
        unstable_steps: Unstable growth steps (the law case: 9).
        stable_steps: Stable growth steps (the law case: 10).
        area_cutoff: The refinement cutoff (the law case: 1e-5).
        through: ``"bridges"`` or ``"partitioned"``.

    Returns:
        The case.

    Raises:
        ValueError: ``through`` is neither ``"bridges"`` nor ``"partitioned"``.
    """
    if through not in ("bridges", "partitioned"):
        raise ValueError("build_inversion_period2 stops at 'bridges' or 'partitioned'")
    session = TangleSession(
        henon_map(*HENON_K41),
        henon_map_inverse(*HENON_K41),
        henon_jacobian(*HENON_K41),
        area_cutoff=area_cutoff,
    )
    fp = session.construct_fixed_point([list(point) for point in K41_INVERSION_PERIOD2_SEED])
    assert (fp.period, fp.k_value) == (2, 4), "the k=4.1 period-2 orbit must have inversion"
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=unstable_steps)
    session.grow_n_times(fp, "stable", num_iterations=stable_steps)
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    case = Case(
        name="inversion_period2",
        session=session,
        fixed_points=[fp],
        pips=[None],
        expect=CaseExpect(periods=(2,), area_cutoff=area_cutoff, inversion=True),
        stage=through,
    )
    if through == "partitioned":
        session.classify_strong_pips()
        case.pips = [session.strong_pip(fp)]
        session.topology()
    return case


def build_k10_inversion(
    *,
    unstable_steps: int = 6,
    stable_steps: int = 5,
    through: Stage = "partitioned",
) -> Case:
    """
    The ``k=10`` map's inversion saddle at ``(-2.3166, 2.3166)``: a numerics fixture only.

    Both eigenvalues are negative (-4.406 and -0.227), so ``k_value = 2 *
    period``; det J = +1. It escapes one step past 6 and has no holes, so it
    is no law case (INV-4); the inversion law cases are :func:`build_inversion`.

    Args:
        unstable_steps: Unstable growth steps (6; 7 does not finish).
        stable_steps: Stable growth steps.
        through: The last stage to run (see :data:`Stage`). ``zones`` adds
            nothing.

    Returns:
        The case.
    """
    area_cutoff = 1e-4
    session = TangleSession(
        henon_map(*HENON_K10),
        henon_map_inverse(*HENON_K10),
        henon_jacobian(*HENON_K10),
        area_cutoff=area_cutoff,
    )
    case = Case(
        name="k10_inversion",
        session=session,
        fixed_points=[],
        pips=[None],
        expect=CaseExpect(periods=(1,), area_cutoff=area_cutoff, inversion=True),
        stage=through,
    )
    if not _reached(through, "fixed_point"):
        return case
    fp = session.construct_fixed_point(list(K10_INVERSION_SEED))
    assert fp.check_inversion(), "the inversion case must have inversion"
    assert fp.k_value == 2 * fp.period
    case.fixed_points = [fp]
    if not _reached(through, "seeded"):
        return case
    session.initialize_both_manifolds(fp)
    if not _reached(through, "grown_unstable"):
        return case
    session.grow_n_times(fp, "unstable", num_iterations=unstable_steps)
    if not _reached(through, "grown"):
        return case
    session.grow_n_times(fp, "stable", num_iterations=stable_steps)
    if not _reached(through, "intersected"):
        return case
    session.compute_intersections([fp])
    if not _reached(through, "bridges"):
        return case
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    if not _reached(through, "pips"):
        return case
    session.classify_strong_pips()
    case.pips = [session.strong_pip(fp)]
    if not _reached(through, "partitioned"):
        return case
    session.topology()
    return case


# --------------------------------------------------------------------------- #
# The b = -1 placeholder
# --------------------------------------------------------------------------- #
def build_orientation_reversing(*, through: Stage = "partitioned") -> Case:
    """
    The ``b = -1`` orientation-reversing placeholder (``k = 4``, seed ``[2, 2]``).

    For ``b = -1`` the Jacobian ``[[2x, 1], [1, 0]]`` has det -1, so exactly
    one eigenvalue is negative at every fixed point and the unstable and
    stable manifolds differ in inversion status. A single ``k_value`` cannot
    express that.

    Args:
        through: The last stage to run, once construction is supported.

    Returns:
        The case (once det J < 0 is supported).

    Raises:
        ValueError: From ``FixedPoint.set_k_value`` via
            ``construct_fixed_point``, today.
    """
    area_cutoff = 1e-4
    session = TangleSession(
        henon_map(*HENON_ORIENTATION_REVERSING),
        henon_map_inverse(*HENON_ORIENTATION_REVERSING),
        henon_jacobian(*HENON_ORIENTATION_REVERSING),
        area_cutoff=area_cutoff,
    )
    fp = session.construct_fixed_point(list(ORIENTATION_REVERSING_SEED))
    case = Case(
        name="orientation_reversing",
        session=session,
        fixed_points=[fp],
        pips=[None],
        expect=CaseExpect(periods=(1,), area_cutoff=area_cutoff),
        stage=through,
    )
    if not _reached(through, "seeded"):
        return case
    session.initialize_both_manifolds(fp)
    if not _reached(through, "grown"):
        return case
    session.grow_n_times(fp, "unstable", num_iterations=6)
    session.grow_until_turnaround(fp, "stable")
    if not _reached(through, "intersected"):
        return case
    session.compute_intersections([fp])
    if not _reached(through, "bridges"):
        return case
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    if not _reached(through, "pips"):
        return case
    session.classify_strong_pips()
    case.pips = [session.strong_pip(fp)]
    if _reached(through, "partitioned"):
        session.topology()
    return case


# --------------------------------------------------------------------------- #
# Near tangency
# --------------------------------------------------------------------------- #
def build_near_tangency() -> Case:
    """
    The period-1 saddle at ``k = -0.999``, grown far into its near-tangent tangle.

    Both arms grow one step at a time until their tail reaches 20 fixed-point
    separations ``sqrt(4 + 4k)`` of canonical distance (lambda is close to 1, so
    a step count would not say how far that is), at ``area_cutoff = 1e-10``
    pinned here; then the crossings are computed. Nothing past
    ``intersected``: the case exists for the resolution guards.

    Returns:
        The case, stage ``intersected``.
    """
    area_cutoff = 1e-10
    session = TangleSession(
        henon_map(*HENON_NEAR_TANGENCY),
        henon_map_inverse(*HENON_NEAR_TANGENCY),
        henon_jacobian(*HENON_NEAR_TANGENCY),
        area_cutoff=area_cutoff,
    )
    k, _b = HENON_NEAR_TANGENCY
    x = 1.0 + 0.5 * np.sqrt(4 + 4 * k)
    fp = session.construct_fixed_point([x, -x])
    session.orient_eigenvectors(fp, _orientation(PERIOD1_ORIENTATION))
    session.initialize_both_manifolds(fp)
    reach = 20 * np.sqrt(4 + 4 * k)
    for stability in ("unstable", "stable"):
        while session.manifolds[(fp, stability, 0, 0)].tail.cdist < reach:
            session.grow_n_times(fp, stability, num_iterations=1)
    session.compute_intersections([fp])
    return Case(
        name="near_tangency",
        session=session,
        fixed_points=[fp],
        pips=[None],
        expect=CaseExpect(periods=(1,), area_cutoff=area_cutoff),
        stage="intersected",
    )


# --------------------------------------------------------------------------- #
# The blast series (D1): blast, recompute the topology from scratch, compare
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SeriesStep:
    """
    One recompute of the D1 loop.

    Attributes:
        blasts: Single-iteration blasts so far.
        holes: ``(bounding_ids, iterate, bridge_side)`` of every hole; ids
            persist across the blasts of one session, so two steps compare.
        classes: ``len(session.bridge_classes())``, the ALPHABET's size (the classes
            with a member the minimal trellis keeps, D2), active and inert together.
        closed: ``len(session.symbolic_dynamics().classes)``, the alphabet closed
            under the word map (C9, D11).
        ambiguous: How many of those classes have several shortest walks with
            different words (D13 picks the side wherever a member lobe exists).
        active: The active classes (a diagnostic, never compared).
        inert: The inert classes (a diagnostic, never compared).
        outside_minimal: Classes with no member in ``minimal_trellis().kept_bridge_ids``.
        references: The reference pseudoneighbor pairs, as id pairs.
        unwitnessed: The references of the previous step that are gone with no
            witness in this one (``helpers.laws.reference_witness``).
        letters: ``(identity, letter)`` of every lettered class, the identity
            being the pair of its elements' ``(branch_key, side, lo_id, hi_id)``
            read off ``session.homotopy_partition()`` here, never through the
            library's own identity.
        pip_depth: ``m``, a report: the deepest registered forward
            iterate of the blasted fixed point's chosen pip, composed here
            through ``trellis.iterate`` (never the library's own count).
        blast_depth: The blast's own ``pip_depth`` of the blasted fixed point
            (``BlastResult``), None for the unblasted step.
        history_free: Whether ``session.trellis()`` read BEFORE the recompute
            shows the crossings and bridges of the recomputed topology's
            trellis (D1: the snapshot never depends on what ran first).
        crossings: ``(N, bridge ids)`` as the blast left them: every registered
            crossing ``unresolved_crossings`` does not list, and the ids of the
            cut bridges.
        recomputed: The same pair after a plain
            ``compute_intersections(preserve_ids=True)`` and
            ``rebuild_bridges()`` (B1: N is one definition on every path).
        p1_tested: How many classes Theorem T was tested on (two or more members
            with readable registered image chains;
            ``helpers.laws.p1_itinerary_disagreements``).
        p1_disagreements: One text per class whose members' ``P1`` itineraries
            differ; recorded, never asserted here.
        chords_checked: How many lobe images cross their base branch at two or
            more interior crossings (``helpers.laws.base_chord_misreadings``).
        chord_misreadings: One text per base-branch gap the iterated cut reads
            otherwise than point-in-polygon; recorded, never asserted here.
        entropy: ``ln`` of the spectral radius of the refined transition matrix
            (:func:`_entropy`; 0 for an empty matrix or ``rho <= 1``).
        held: The bridges the step's blast held (``BlastResult.held``; 0 for the
            unblasted step and for a full blast).
    """

    blasts: int
    holes: frozenset
    classes: int
    closed: int
    ambiguous: int
    active: int
    inert: int
    outside_minimal: int
    references: frozenset
    unwitnessed: frozenset
    letters: frozenset
    pip_depth: int
    blast_depth: Optional[int]
    history_free: bool
    crossings: tuple[frozenset, frozenset]
    recomputed: tuple[frozenset, frozenset]
    p1_tested: int
    p1_disagreements: tuple
    chords_checked: int
    chord_misreadings: tuple
    entropy: float
    held: int = 0


class Series(NamedTuple):
    """One blast series: the zone-stage build, how many blasts, at which separation."""

    build: Callable[[], Case]
    blasts: int
    min_separation: float = ZONE_MIN_SEPARATION


#: ``name -> Series``. Each blasts the build's LARGEST zone (the only one but on
#: ``nested_outer_series``, whose outer zone it is). k28 switches its reference
#: pair at blasts 2 and 8, p3 at blast 6. ``k32_series`` and
#: ``nested_outer_series`` are E11's ``k32`` / ``nested_outer`` (2026-10-07), run
#: to the depth where the base decreased (BCG-5: k32 blast 10, nested_outer blast
#: 11, inversion blast 4). Since D12 was reversed (2026-10-31) the class count
#: may fall there (k32 blast 10: a deeper crossing refutes a reference); the
#: entropy may not.
SERIES: Mapping[str, Series] = MappingProxyType(
    {
        "k28_series": Series(lambda: build_k28(blasts=0, through="zones"), 8),
        "p3_series": Series(lambda: build_period3(through="zones"), 9),
        "inversion_series": Series(lambda: build_inversion(through="zones"), 4),
        "k32_series": Series(build_k32_zone, 10),
        "nested_outer_series": Series(
            lambda: build_nested(outer_blasts=0, through="zones"), 11, min_separation=1e-4
        ),
    }
)


def blast_series(name: str, representatives: Optional[str] = None) -> tuple[SeriesStep, ...]:
    """
    The D1 loop on one session: recompute the topology, blast once, repeat.

    Each blast is one iteration of the build's largest zone, of its own fixed
    point only, with ``diagnose=False``; the session carries the pips. Each
    step ends with a plain recompute of the crossings and bridges, recorded
    beside what the blast left.

    Args:
        name: A key of :data:`SERIES`.
        representatives: The blasts' representative policy (None: every
            interior bridge is iterated, the series as the D1 laws state it).

    Returns:
        One step per recompute, the unblasted one first.
    """
    series = SERIES[name]
    case = series.build()
    session = case.session
    zone = max(case.zones, key=lambda zone: zone.area)
    fp = zone.fixed_point
    steps: list[SeriesStep] = []
    for blast in range(series.blasts + 1):
        result = None
        if blast:
            result = _blast(
                session, zone, [fp], diagnose=False, min_separation=series.min_separation,
                representatives=representatives,
            )
        shown = trellis_view(session.trellis())
        table = session.bridge_classes()
        kept = set(session.minimal_trellis().kept_bridge_ids)
        family = session.homotopy_partition()
        trellis = session.topology().trellis
        dynamics = session.symbolic_dynamics()
        p1_tested, _skipped, p1_disagreements = p1_itinerary_disagreements(session)
        chords_checked, chord_misreadings = base_chord_misreadings(session)
        references = frozenset(p.as_tuple() for p in trellis.reference_pseudoneighbors)
        lost = steps[-1].references - references if steps else frozenset()
        steps.append(
            SeriesStep(
                blasts=blast,
                holes=frozenset(
                    (hole.bounding_ids, hole.iterate, hole.bridge_side) for hole in trellis.holes
                ),
                classes=len(table),
                closed=len(dynamics.classes),
                ambiguous=sum(cd.ambiguous for cd in dynamics.classes.values()),
                active=len(table.active),
                inert=len(table.inert),
                outside_minimal=sum(
                    not any(m.bridge_id in kept for m in entry.members) for entry in table
                ),
                references=references,
                unwitnessed=frozenset(
                    pair for pair in lost if reference_witness(trellis, pair) is None
                ),
                letters=frozenset(
                    (_class_identity(family, entry.bridge_class), entry.letter)
                    for entry in table
                    if entry.letter is not None
                ),
                pip_depth=_pip_depth(session, fp),
                blast_depth=None if result is None else result.pip_depth.get(fp),
                history_free=shown == trellis_view(session.topology().trellis),
                crossings=_crossings_and_bridges(session.workbench),
                recomputed=_recomputed(session.workbench),
                p1_tested=p1_tested,
                p1_disagreements=p1_disagreements,
                chords_checked=chords_checked,
                chord_misreadings=chord_misreadings,
                entropy=_entropy(dynamics.transition_matrix()[1]),
                held=0 if result is None else result.held,
            )
        )
    return tuple(steps)


def _entropy(matrix: np.ndarray) -> float:
    """``ln rho(matrix)``, 0 unless the spectral radius exceeds 1.

    That covers an empty or nilpotent matrix (``rho = 0``) and a defective
    eigenvalue at 1 that eigvals reads just below it.
    """
    rho = np.abs(np.linalg.eigvals(matrix)).max() if matrix.size else 0.0
    return float(np.log(rho)) if rho > 1 else 0.0


def _pip_depth(session: TangleSession, fp: "FixedPoint") -> int:
    """The largest ``m`` with ``f^m(q0)`` registered, ``q0`` the chosen pip of ``fp``."""
    trellis = session.trellis(fp)
    pip = trellis.registry.find(session.chosen_strong_pip(fp))
    m = 0
    while trellis.iterate(pip, m + 1) is not None:
        m += 1
    return m


def _crossings_and_bridges(workbench) -> tuple[frozenset, frozenset]:
    """``(N, bridge ids)``: the registry less what it lists as unresolved, and the cut."""
    registry = workbench.intersection_registry
    unresolved = set(registry.unresolved_crossings(workbench.resolution_ratio_threshold))
    return (
        frozenset(iid for iid in registry.all_ids() if iid not in unresolved),
        frozenset(bridge.id for bridge in workbench.bridges if bridge.id is not None),
    )


def _recomputed(workbench) -> tuple[frozenset, frozenset]:
    """:func:`_crossings_and_bridges` after a plain recompute of the crossings and bridges."""
    workbench.compute_intersections(workbench.fixed_points, preserve_ids=True, diagnose=False)
    workbench.rebuild_bridges()
    return _crossings_and_bridges(workbench)


def _class_identity(family, bridge_class) -> frozenset:
    """The pair of the class's elements, each by its branch, side and bounding crossings."""
    elements = (family.element(ref) for ref in (bridge_class.source, bridge_class.target))
    return frozenset((e.branch_key, e.side, e.lo_id, e.hi_id) for e in elements)


# --------------------------------------------------------------------------- #
# Registries: law cases, known issues, not-applicable pairs
# --------------------------------------------------------------------------- #
#: The cases every physical law runs on (decision 9).
LAW_CASES: tuple[str, ...] = (
    "k10",
    "k28_one_blast",
    "k28_two_blasts",
    "p3",
    "nested",
    "nested_inner_blasts",
    "inversion",
    "inversion_one_blast",
    "inversion_two_blasts",
    "inversion_period2",
    "heteroclinic",
)

def _named(case: Case, name: str) -> Case:
    """``case`` under its BUILDERS name, where its builder's own name would be a law case's."""
    case.name = name
    return case


#: Every named case and its frozen builder.
BUILDERS: Mapping[str, Callable[[], Case]] = MappingProxyType(
    {
        "k10": build_k10,
        "k28_one_blast": lambda: build_k28(blasts=1),
        "k28_two_blasts": lambda: build_k28(blasts=2),
        "p3": build_period3,
        "nested": build_nested,
        "nested_inner_blasts": lambda: build_nested(inner_blasts=2),
        "heteroclinic": build_heteroclinic,
        "heteroclinic_17u": lambda: _named(build_heteroclinic(p3_unstable_steps=17), "heteroclinic_17u"),
        "nested_unblasted": lambda: build_nested(
            outer_blasts=0, p1_area_cutoff=1e-4, through="zones"
        ),
        "inversion": build_inversion,
        "inversion_one_blast": lambda: build_inversion(blasts=1),
        "inversion_two_blasts": lambda: build_inversion(blasts=2),
        "inversion_period2": build_inversion_period2,
        "orientation_reversing": build_orientation_reversing,
        "p3_15": lambda: build_period3(unstable_steps=15),
        "p3_16": lambda: build_period3(unstable_steps=16),
        "p3_6_blasts": lambda: build_period3(blasts=6),
        "p3_16_12": lambda: _named(build_period3(unstable_steps=16, stable_steps=12), "p3_16_12"),
        "inversion_period2_10": lambda: _named(
            build_inversion_period2(unstable_steps=10, stable_steps=10, area_cutoff=1e-6),
            "inversion_period2_10",
        ),
        "p3_coarse_blasts": build_period3_coarse_blasts,
        "near_tangency": build_near_tangency,
        "k28_deep_blast_sliver": build_k28_deep_blast_sliver,
    }
)


@dataclass(frozen=True)
class KnownIssue:
    """
    A known bug recorded as ``xfail(strict=True)``.

    Attributes:
        reason: Why the pair fails today, led by the catalog issue id
            (``docs/open_issues.md``).
        wp: The work package that fixes it (``"WP3"``), the decision that
            keeps it (``"C3"``), or the open author decision it waits on
            (``"DEC-<n>"``, or ``"author: <what>"`` for one not yet numbered).
        raises: The exception type the failure raises, when it is one.
    """

    reason: str
    wp: str
    raises: Optional[type[BaseException]] = None


def _each(laws: tuple[str, ...], names: tuple[str, ...], issue: KnownIssue) -> dict:
    """One ``(law, case)`` entry per law and case, all carrying ``issue``."""
    return {(law, name): issue for law in laws for name in names}


_INV_12_UNBLASTED = KnownIssue(
    "INV-12: with no unstable pip cut the unblasted inversion case keeps one reference "
    "whose hole never propagates (1 hole, 1 class, 0 active): the law checks nothing",
    "author: inversion investigation",
)
_INV_12_SERIES = KnownIssue(
    "INV-12: inversion_series blast 4 has two ambiguous walks and its entropy falls "
    "0.264 -> 0 (reference switch to (56, 57))",
    "author: inversion investigation",
)


#: ``(law_id, case) -> KnownIssue``: every pair fails today and names the WP that
#: fixes it. The law ids are the ``helpers.laws`` check names without
#: ``check_`` (or a standalone test's id); a law-tier pair is left out of its
#: layer's run and gets a set-aside test of its own (``helpers.law_tier``), so
#: a fix flips loudly to XPASS. The registry may only shrink, except for a new
#: law whose entry names the later WP that fixes it, and the INV-12 entries the
#: reversal of D12 exposed (2026-10-31, ``.refactor/decisions.md``).
KNOWN_ISSUES: Mapping[tuple[str, str], KnownIssue] = MappingProxyType(
    {
        ("orientation_reversing_case_builds", "orientation_reversing"): KnownIssue(
            "INV-11: det J < 0 is out of scope (C3); FixedPoint.set_k_value rejects a "
            "saddle with exactly one negative eigenvalue",
            "C3",
            raises=ValueError,
        ),
        **_each(
            (
                "anchor_bridge_class_leads_its_tangle",
                "registered_images_keep_markers",
                "partition_singletons",
                "backward_holes_keep_bridge_side",
                "anchor_opening_faces_its_hole",
                "propagated_holes_land_on_predicted_branch",
                "propagation_terminates",
                "propagated_hole_inside_container",
                "refined_matrix_is_an_in_splitting",
                "matrix_is_token_counts",
                "base_chords_are_inside_the_lobe_image",
            ),
            ("inversion",),
            _INV_12_UNBLASTED,
        ),
        **_each(
            ("series_entropy_never_decreases", "series_no_walk_is_ambiguous"),
            ("inversion_series",),
            _INV_12_SERIES,
        ),
    }
)

_NO_PINCH = "no hole pair pinches a crossing"
_NO_ANCHOR_OPENING = "no hole is bounded by the anchor"
_NO_IMAGE_BRIDGES = "the minimal trellis maps no hole bridge forward (no image bridges)"
_NO_REFINEMENT = "no class splits (no refined children)"
_NO_DISJOINT_WALLS = "no two classes with decided walks have arcs on disjoint unstable branches"
_NO_TWO_CHAINED_MEMBERS = "no class has two members with registered image chains"

#: ``(law_id, case) -> reason``: the law has nothing to check on that case. The
#: layer test leaves it out; its set-aside test runs it, insists it still checks
#: nothing and reports a visible skip.
NOT_APPLICABLE: Mapping[tuple[str, str], str] = MappingProxyType(
    {
        ("arrangement_image_of", "k10"): (
            "the one face carrying a region's mapped corners is a proper sub-face, "
            "which image_of rejects; no region has a closed image"
        ),
        ("arrangement_preimage_inverts_image", "k10"): "no region has a closed image",
        **{
            ("partition_singletons", name): _NO_PINCH
            for name in ("k10", "k28_two_blasts")
        },
        **{
            ("anchor_opening_faces_its_hole", name): _NO_ANCHOR_OPENING
            for name in ("k10", "k28_two_blasts")
        },
        **{
            ("iterated_cut_provenance", name): _NO_IMAGE_BRIDGES
            for name in ("k28_one_blast", "p3", "inversion", "inversion_period2")
        },
        **{
            (law, name): _NO_REFINEMENT
            for law in ("member_matching_consistent",)
            for name in ("k28_one_blast", "p3", "inversion", "inversion_period2")
        },
        **{
            ("every_active_class_walks", name): (
                "every active class lands on a singleton; its word is the trellis "
                "itinerary (D3)"
            )
            for name in ("k28_one_blast", "p3", "inversion", "inversion_period2")
        },
        **{
            ("walks_avoid_disjoint_images", name): _NO_DISJOINT_WALLS
            for name in (
                "k10", "k28_one_blast", "k28_two_blasts", "p3", "inversion", "inversion_period2"
            )
        },
        **{
            ("members_share_the_P1_itinerary", name): _NO_TWO_CHAINED_MEMBERS
            for name in (
                "k28_one_blast", "k28_two_blasts", "p3", "nested", "nested_inner_blasts",
                "inversion_one_blast", "inversion_period2",
            )
        },
        **{
            ("k10_is_the_full_shift", name): "the textbook answer is known for k=10 only"
            for name in LAW_CASES
            if name != "k10"
        },
    }
)


def issue_marks(law_id: str, case: str) -> list:
    """
    The pytest marks of one ``(law, case)`` pair.

    Args:
        law_id: The law's identifier.
        case: The case name.

    Returns:
        ``[xfail(strict=True)]`` (with ``raises`` when given) for a known
        issue, ``[]`` otherwise.
    """
    issue = KNOWN_ISSUES.get((law_id, case))
    if issue is None:
        return []
    kwargs = {"strict": True, "reason": f"{issue.wp}: {issue.reason}"}
    if issue.raises is not None:
        kwargs["raises"] = issue.raises
    return [pytest.mark.xfail(**kwargs)]

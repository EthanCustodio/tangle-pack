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
  resonance zone), then one or two single-iteration blasts of that zone with
  ``min_separation = 1e-5``, the pip restored after each
  (:func:`_zone_at_image_pip`).
* ``p3`` -- the inner period-3 orbit of the nested map ``(2, 1)`` alone,
  delegated to :func:`tanglepack.examples.henon_cases.build_period3` with every
  keyword explicit.
* ``nested`` / ``nested_inner_blasts`` -- period 1 + period 3 with two outer
  blasts, and no or two inner blasts, delegated to
  :func:`tanglepack.examples.henon_cases.build_nested` with every keyword
  explicit (``p1_area_cutoff = 1e-7``).
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
  eigenvalues negative (-1.863 / -0.537), ``k_value = 4``; 9 unstable / 9
  stable steps at ``area_cutoff = 1e-5``, no zone (at 10 steps seven laws fail,
  none of them a period-k hole gap: an element opened by two holes of one
  origin a lap apart, one on a cross-branch container, and landings past the
  next branch's partition -- WP7 / WP8; it escapes past ~8 steps at k=4.25).
* ``k10_inversion`` -- the ``k=10`` map's inversion saddle at ``(-2.3166,
  2.3166)``, 6 unstable / 5 stable steps at ``area_cutoff = 1e-4``: the
  numerics fixture behind ``henon_inversion*`` only (it escapes at 7 steps and
  has no holes, so it is no law case).
* ``orientation_reversing`` -- the ``b = -1`` placeholder (``k = 4``, seed
  ``[2, 2]``, ``area_cutoff = 1e-4``): ``construct_fixed_point`` raises
  ``ValueError`` because a single ``k_value`` cannot describe a saddle with
  exactly one negative eigenvalue; det J < 0 is out of scope (C3).
* ``p3_15`` / ``p3_16`` / ``p3_6_blasts`` -- the deep period-3 runs that are
  known open problems (an unreachable class and a virtual ``new1``, PK-1).
* ``near_tangency`` -- the period-1 ``b = 1`` saddle at ``k = -0.999``, close to
  the saddle-centre bifurcation at ``k = -1`` where the splitting angle is
  exponentially small, each arm grown to 20 fixed-point separations of
  canonical distance at its own ``area_cutoff = 1e-10``, through
  ``compute_intersections`` only. Its crossings are near-tangent polylines
  straddling each other: the resolution guards' positive case. NOT a law case.
* ``k28_deep_blast_sliver`` -- ``k28`` from the analytic seed, blasted nine
  times: a crossing sits a sliver narrower than ``2 * cdist_tol`` past the pip
  (stress campaign E11, 2026-10-07; the exact-ownership regression).

The blast series (:data:`SERIES`, :func:`blast_series`) run the D1 loop --
blast once, recompute the topology from scratch, compare -- from the zone
stage of ``k28`` (8 blasts), ``p3`` (9) and ``inversion`` (2).

Dev Notes:

* ``henon_cases`` has no stage parameter (and ``src/`` is not changed by the
  test refactor), so a nested or period-3 build stopped before the partition
  (``through="zones"``) is a local copy of the ``henon_cases`` recipe; only the
  full ``"partitioned"`` build delegates.
* Registry ids are NOT reproducible between two builds in one process (P6:
  three ``build_k10()`` calls give three permutations of the same crossings),
  so a test never compares ids across builds and never picks "the lowest id"
  as a canonical choice.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING, Callable, Literal, Mapping, Optional

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

if TYPE_CHECKING:  # pragma: no cover - typing only
    from tanglepack import FixedPoint, IntersectionRegistry, TangleWorkbench
    from tanglepack.loom.ResonanceZone import ResonanceZone

#: How far a builder runs. Each stage includes every earlier one:
#: ``empty`` (a session, no fixed point), ``fixed_point`` (constructed and
#: oriented), ``seeded`` (fundamental segments), ``grown_unstable`` (unstable
#: growth only), ``grown`` (both manifolds), ``intersected``
#: (``compute_intersections``), ``bridges`` (trim, bridges, iterate table),
#: ``pips`` (strong pips classified), ``zones`` (resonance zones, blasts),
#: ``partitioned`` (pseudoneighbors, holes, stable partition).
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
    """

    periods: tuple[int, ...]
    area_cutoff: float
    inversion: bool = False
    n_tangles: int = 1
    blasts: int = 0


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


def _partition(session: TangleSession) -> None:
    """Pseudoneighbors, holes and the stable partition through the fan-outs."""
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()


def _repin(session: TangleSession, fixed_points: list, pips: list) -> None:
    """Classify the rebuilt trellis again and restore the chosen pips."""
    session.classify_strong_pips()
    for fp, pip in zip(fixed_points, pips):
        session.set_strong_pip(fp, pip)


def _blast(
    session: TangleSession,
    zone: "ResonanceZone",
    fixed_points: list,
    pips: list,
    *,
    diagnose: bool = True,
) -> None:
    """One single-iteration blast of ``zone`` at ``ZONE_MIN_SEPARATION``, the pips restored."""
    session.blast_zone(
        zone,
        num_iterations=1,
        fixed_point=fixed_points,
        min_separation=ZONE_MIN_SEPARATION,
        diagnose=diagnose,
    )
    _repin(session, fixed_points, pips)


def _zone_at_image_pip(
    session: TangleSession, fp: "FixedPoint", blasts: int
) -> tuple[int, "ResonanceZone"]:
    """The k=2.8 zone recipe: the zone trimmed at the pip ``f(q0)``, blasted ``blasts`` times.

    Returns:
        ``(pip, zone)``; the pip is restored after the trim and after every blast.
    """
    session.classify_strong_pips()
    trellis = session.trellis(fp)
    pip = trellis.iterate(trellis.strong_pip, 1)
    assert pip is not None, "the default strong pip must have a registered image"
    zone = session.resonance_zone(pip)
    _repin(session, [fp], [pip])
    for _ in range(blasts):
        _blast(session, zone, [fp], [pip])
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
    _partition(session)
    return case


# --------------------------------------------------------------------------- #
# k = 2.8
# --------------------------------------------------------------------------- #
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
        The case, the pip ``f(q0)`` restored after every blast.

    Raises:
        ValueError: ``through`` is neither ``"zones"`` nor ``"partitioned"``.
    """
    if through not in ("zones", "partitioned"):
        raise ValueError("build_k28 stops at 'zones' or 'partitioned'")
    area_cutoff = 1e-7
    session = TangleSession(
        henon_map(*HENON_K28),
        henon_map_inverse(*HENON_K28),
        henon_jacobian(*HENON_K28),
        area_cutoff=area_cutoff,
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
    if through == "partitioned":
        _partition(session)
    return Case(
        name={0: "k28_unblasted", 1: "k28_one_blast", 2: "k28_two_blasts"}.get(
            blasts, f"k28_{blasts}_blasts"
        ),
        session=session,
        fixed_points=[fp],
        pips=[pip],
        expect=CaseExpect(periods=(1,), area_cutoff=area_cutoff, blasts=blasts),
        zones=[zone],
        stage=through,
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
    repinned zone, unblasted (the blast series starts there).

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
    """The ``henon_cases.build_period3`` recipe stopped at the repinned zone."""
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
    _repin(session, [fp3], [pip])
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
    ``min_separation = 1e-5`` (no same-stability pass), the pips restored after
    each, and the stable manifold is partitioned. At 1e-5 and 6 blasts some hole
    bridges hold 2-4 nodes, where a polyline reading of their rows disagrees with
    the crossing signs (NUM-1).

    Args:
        blasts: Single-iteration blasts of the zone.
        area_cutoff: The manifold machine's refinement cutoff during the blasts.

    Returns:
        The partitioned case.
    """
    case = build_period3(through="zones")
    session = case.session
    session.area_cutoff = area_cutoff
    (zone,) = case.zones
    for _ in range(blasts):
        session.blast_zone(
            zone, num_iterations=1, fixed_point=case.fixed_points,
            min_separation=1e-5, diagnose=False,
        )
        _repin(session, case.fixed_points, case.pips)
    _partition(session)
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
    """The ``henon_cases.build_nested`` recipe stopped before the repin and blasts."""
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
    session.area_cutoff = p3_area_cutoff
    fp3 = session.construct_fixed_point([list(p) for p in P3_ORBIT_SEED])
    session.orient_eigenvectors(fp3, _orientation(PERIOD3_ORIENTATION))
    session.initialize_both_manifolds(fp3)
    session.grow_n_times(fp3, "unstable", num_iterations=p3_unstable_steps)
    session.grow_n_times(fp3, "stable", num_iterations=p3_stable_steps)

    session.area_cutoff = p1_area_cutoff
    fp1 = session.construct_fixed_point(list(P3_OUTER_SADDLE_SEED))
    session.orient_eigenvectors(fp1, _orientation(PERIOD1_ORIENTATION))
    session.initialize_both_manifolds(fp1)
    session.grow_n_times(fp1, "unstable", num_iterations=p1_unstable_steps)
    session.grow_until_turnaround(fp1, "stable")
    case.fixed_points = [fp1, fp3]
    if not _reached(through, "intersected"):
        return case

    session.compute_intersections([fp3, fp1])
    if not _reached(through, "bridges"):
        return case
    session.trim_stable_manifolds(fp3)
    session.trim_stable_manifolds(fp1)
    session.create_bridges(fp3)
    session.create_bridges(fp1)
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
            local copy of its recipe and stops there, before the repin.

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
        The case, the pip ``f(q0)`` restored after every blast.

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
        _partition(session)
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


def build_inversion_period2(*, through: Stage = "partitioned") -> Case:
    """
    Inversion on a period-k orbit: the ``k=4.1`` period-2 orbit, no zone.

    Both eigenvalues of ``DM^2`` are negative (-1.863 and -0.537), so
    ``k_value = 2 * period = 4`` and each manifold's four pieces are one chain
    over the two orbit points and both branches; det J = +1. 9 unstable and 9
    stable steps at ``area_cutoff = 1e-5``, then the stable manifold trimmed at
    its outermost crossings and partitioned at the default strong pip.

    Args:
        through: ``"bridges"`` or ``"partitioned"``.

    Returns:
        The case.

    Raises:
        ValueError: ``through`` is neither ``"bridges"`` nor ``"partitioned"``.
    """
    if through not in ("bridges", "partitioned"):
        raise ValueError("build_inversion_period2 stops at 'bridges' or 'partitioned'")
    area_cutoff = 1e-5
    session = TangleSession(
        henon_map(*HENON_K41),
        henon_map_inverse(*HENON_K41),
        henon_jacobian(*HENON_K41),
        area_cutoff=area_cutoff,
    )
    fp = session.construct_fixed_point([list(point) for point in K41_INVERSION_PERIOD2_SEED])
    assert (fp.period, fp.k_value) == (2, 4), "the k=4.1 period-2 orbit must have inversion"
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=9)
    session.grow_n_times(fp, "stable", num_iterations=9)
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
        _partition(session)
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
    _partition(session)
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
        _partition(session)
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
        classes: ``len(session.bridge_classes())``, active and inert together.
        active: The active classes (a diagnostic, never compared).
        inert: The inert classes (a diagnostic, never compared).
        outside_minimal: Classes with no member in ``minimal_trellis().kept_bridge_ids``.
    """

    blasts: int
    holes: frozenset
    classes: int
    active: int
    inert: int
    outside_minimal: int


#: ``name -> (zone-stage builder, blasts)``. k28 switches its reference pair at
#: blasts 2 and 8, p3 at blast 6; inversion stops at 2.
SERIES: Mapping[str, tuple[Callable[[], Case], int]] = MappingProxyType(
    {
        "k28_series": (lambda: build_k28(blasts=0, through="zones"), 8),
        "p3_series": (lambda: build_period3(through="zones"), 9),
        "inversion_series": (lambda: build_inversion(through="zones"), 2),
    }
)


def blast_series(name: str) -> tuple[SeriesStep, ...]:
    """
    The D1 loop on one session: recompute the topology, blast once, repeat.

    Each blast is one iteration at ``ZONE_MIN_SEPARATION`` with
    ``diagnose=False``, the pips restored after it.

    Args:
        name: A key of :data:`SERIES`.

    Returns:
        One step per recompute, the unblasted one first.
    """
    build, blasts = SERIES[name]
    case = build()
    session = case.session
    (zone,) = case.zones
    steps = []
    for blast in range(blasts + 1):
        if blast:
            _blast(session, zone, case.fixed_points, case.pips, diagnose=False)
        _partition(session)
        table = session.bridge_classes()
        kept = set(session.minimal_trellis().kept_bridge_ids)
        steps.append(
            SeriesStep(
                blasts=blast,
                holes=frozenset(
                    (hole.bounding_ids, hole.iterate, hole.bridge_side)
                    for fp in case.fixed_points
                    for hole in session.trellis(fp).holes
                ),
                classes=len(table),
                active=len(table.active),
                inert=len(table.inert),
                outside_minimal=sum(
                    not any(m.bridge_id in kept for m in entry.members) for entry in table
                ),
            )
        )
    return tuple(steps)


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
)

#: Every named case and its frozen builder.
BUILDERS: Mapping[str, Callable[[], Case]] = MappingProxyType(
    {
        "k10": build_k10,
        "k28_one_blast": lambda: build_k28(blasts=1),
        "k28_two_blasts": lambda: build_k28(blasts=2),
        "p3": build_period3,
        "nested": build_nested,
        "nested_inner_blasts": lambda: build_nested(inner_blasts=2),
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
        "p3_coarse_blasts": build_period3_coarse_blasts,
        "near_tangency": build_near_tangency,
        "k28_deep_blast_sliver": lambda: build_k28(blasts=9, seed=K28_ANALYTIC_SADDLE_SEED),
    }
)


@dataclass(frozen=True)
class KnownIssue:
    """
    A known bug recorded as ``xfail(strict=True)``.

    Attributes:
        reason: Why the pair fails today, led by the catalog issue id
            (``docs/open_issues.md``).
        wp: The work package that fixes it (``"WP3"``), or the decision that
            keeps it (``"C3"``).
        raises: The exception type the failure raises, when it is one.
    """

    reason: str
    wp: str
    raises: Optional[type[BaseException]] = None


def _each(laws: tuple[str, ...], names: tuple[str, ...], issue: KnownIssue) -> dict:
    """One ``(law, case)`` entry per law and case, all carrying ``issue``."""
    return {(law, name): issue for law in laws for name in names}


_INVERSION = ("inversion", "inversion_one_blast", "inversion_two_blasts")
_NUM16 = KnownIssue(
    "NUM-16: a zero-area region after a blast reports a representative point inside "
    "its neighbour",
    "WP10",
)
_BCG1 = KnownIssue(
    "BCG-1: every trellis bridge is classed, blast children included, so a class can "
    "have no member in the minimal trellis (D2)",
    "WP6",
)

#: ``(law_id, case) -> KnownIssue``: every pair fails today and names the WP that
#: fixes it. The law ids are the ``helpers.laws`` check names without
#: ``check_`` (or a standalone test's id); a law-tier pair is left out of its
#: layer's run and gets a set-aside test of its own (``helpers.law_tier``), so
#: a fix flips loudly to XPASS. The registry may only shrink, except for a new
#: law whose entry names the later WP that fixes it.
KNOWN_ISSUES: Mapping[tuple[str, str], KnownIssue] = MappingProxyType(
    {
        **_each(
            ("arrangement_regions_disjoint",),
            ("k28_one_blast", "k28_two_blasts", "inversion_one_blast", "inversion_two_blasts"),
            _NUM16,
        ),
        **_each(
            ("arrangement_regions_sound", "arrangement_regions_disjoint"),
            ("nested_inner_blasts",),
            KnownIssue(
                "NUM-16: a zero-area region after the inner blasts has no verified "
                "representative point (regions_disjoint then raises TypeError)",
                "WP10",
            ),
        ),
        ("orientation_reversing_case_builds", "orientation_reversing"): KnownIssue(
            "INV-11: det J < 0 is out of scope (C3); FixedPoint.set_k_value rejects a "
            "saddle with exactly one negative eigenvalue",
            "C3",
            raises=ValueError,
        ),
        **_each(
            ("open_p3_deep_runs",),
            ("p3_15", "p3_16", "p3_6_blasts"),
            KnownIssue(
                "PK-1: the deep period-3 runs (15 or 16 steps, 6 blasts) have a "
                "virtual new1 (the unreachable class went with the geometric "
                "openings, WP4)",
                "WP7",
            ),
        ),
        ("arrangement_euler", "inversion_two_blasts"): KnownIssue(
            "INT-9: the dense arrangement draws an arc between consecutive blast-image "
            "crossings no bridge spans; dropping those arcs gives V - E + F = 2",
            "WP7",
        ),
        ("anchor_bridge_class_leads_its_tangle", "inversion"): KnownIssue(
            "PK-3 / DEC-19: the reference pair's +1 image lies past the computed u(0.0) "
            "(the inversion halves are grown one map step apart, a stagger), so s(0.1) "
            "gets no window hole (named in a WARNING since WP5) and the anchor bridge's "
            "image is a loop (an inert anchor class)",
            "WP5",
        ),
        **_each(("alphabet_inside_minimal_trellis",), _INVERSION, _BCG1),
        **_each(
            ("series_class_count_constant_while_holes_constant",),
            ("k28_series", "p3_series"),
            KnownIssue(
                "BCG-1/BCG-2: blast children with no image are classed and orphan loops "
                "stand alone, so the class count grows while the hole set is constant",
                "WP6",
            ),
        ),
        **_each(("series_alphabet_inside_minimal_trellis",), tuple(SERIES), _BCG1),
        ("three_tangles_take_three_families", "hand_built"): KnownIssue(
            "TEST-5(c): TANGLE_COLOR_FAMILIES holds two families, so a third tangle "
            "draws from none",
            "WP9",
        ),
    }
)

_NO_PINCH = "no hole pair pinches a crossing"
_NO_ANCHOR_OPENING = "no hole is bounded by the anchor"
_NO_IMAGE_BRIDGES = "the minimal trellis maps no hole bridge forward (no image bridges)"
_NO_REFINEMENT = "no class splits (no refined children)"

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
            for name in ("k10", "k28_two_blasts", "inversion")
        },
        **{
            ("anchor_opening_faces_its_hole", name): _NO_ANCHOR_OPENING
            for name in ("k10", "k28_two_blasts", "inversion")
        },
        **{
            ("iterated_cut_provenance", name): _NO_IMAGE_BRIDGES
            for name in ("k28_one_blast", "p3", "inversion", "inversion_period2")
        },
        **{
            (law, name): _NO_REFINEMENT
            for law in ("refined_children_inherit_word", "member_matching_consistent")
            for name in ("k28_one_blast", "p3", "inversion", "inversion_period2")
        },
        **{
            ("every_active_class_walks", name): (
                "every active class lands on a singleton; its word is the trellis "
                "itinerary (D3)"
            )
            for name in ("k28_one_blast", "p3")
        },
        **{
            (law, "inversion"): (
                "one reference hole, whose backward image lies past the last crossing of "
                "the other unstable branch"
            )
            for law in (
                "backward_holes_keep_bridge_side",
                "propagated_holes_land_on_predicted_branch",
                "propagation_terminates",
            )
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

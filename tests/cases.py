"""The tests-only case builders: every Hénon tangle the suite builds, frozen here.

Each builder returns a :class:`Case`: the session, its fixed points (outermost
first), the chosen strong pips, the resonance zones it made and the
:class:`CaseExpect` facts the law tier can assume about the case. The builder
parameters are FROZEN in this file. They are never read from the defaults of
:mod:`tanglepack.examples.henon_cases` (those serve the figure scripts and may
move), so a change of a figure recipe can never silently change a test case.

Cases:

* ``k10`` -- the ``k=10, b=1`` binary horseshoe, the plain saddle at
  ``(4, -4)``: 9 unstable steps, the stable manifold grown to turnaround, no
  analytic Jacobian (the conftest recipe since the start of the suite).
* ``k28_one_blast`` / ``k28_two_blasts`` -- ``k=2.8``: 10 unstable steps at
  ``area_cutoff = 1e-7``, the stable manifold trimmed at the pip ``f(q0)`` (a
  resonance zone), then one or two single-iteration blasts of that zone with
  ``min_separation = 1e-5``, the pip restored after each.
* ``p3`` -- the inner period-3 orbit of the nested map ``(2, 1)`` alone,
  delegated to :func:`tanglepack.examples.henon_cases.build_period3` with every
  keyword explicit.
* ``nested`` -- period 1 + period 3 with two outer blasts, delegated to
  :func:`tanglepack.examples.henon_cases.build_nested` with every keyword
  explicit (``p1_area_cutoff = 1e-7``).
* ``nested_unblasted`` -- the nested map with no blast and the outer saddle
  grown at ``p1_area_cutoff = 1e-4`` (the old ``henon_p3_session`` recipe, kept
  as a parameter: Phase 0 probe P6 showed the cutoff changes the outer tangle).
* ``inversion`` -- the ``k=10`` map's OTHER saddle at ``(-2.3166, 2.3166)``:
  both eigenvalues negative, ``k_value = 2 * period``, det J = +1 (inversion,
  NOT orientation reversing); 6 unstable / 5 stable steps (7 unstable steps
  does not finish, P2).
* ``orientation_reversing`` -- the ``b = -1`` placeholder (``k = 4``, seed
  ``[2, 2]``): ``construct_fixed_point`` raises ``ValueError`` today because a
  single ``k_value`` cannot describe a saddle with exactly one negative
  eigenvalue (P4).
* ``p3_15`` / ``p3_16`` / ``p3_6_blasts`` -- the deep period-3 runs that are
  known open problems (an unreachable class and a virtual ``new1``, P3).

Dev Notes:

* ``henon_cases`` has no stage parameter (and ``src/`` is not changed by the
  test refactor), so a nested build stopped before the partition
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
#: The ``b = -1`` orientation-reversing placeholder map (P4).
HENON_ORIENTATION_REVERSING: tuple[float, float] = (4, -1)

K10_SADDLE_SEED: tuple[float, float] = (4.0, -4.0)
K10_INVERSION_SEED: tuple[float, float] = (-2.3166, 2.3166)
K28_SADDLE_SEED: tuple[float, float] = (4.0, -4.0)
ORIENTATION_REVERSING_SEED: tuple[float, float] = (2.0, 2.0)
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
        inversion: Whether any fixed point has inversion (``k_value == 2 * period``).
        n_tangles: The number of fixed points.
        blasts: Single-iteration blasts applied (outer blasts for nested).
        has_image_bridges: Whether the minimal trellis is expected to carry
            image bridges (None when not stated by the author).
    """

    periods: tuple[int, ...]
    inversion: bool = False
    n_tangles: int = 1
    blasts: int = 0
    has_image_bridges: Optional[bool] = None


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
    session = TangleSession(henon_map(*HENON_K10), henon_map_inverse(*HENON_K10))
    case = Case(
        name="k10",
        session=session,
        fixed_points=[],
        pips=[None],
        expect=CaseExpect(periods=(1,), has_image_bridges=True),
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
    session.infer_iterate_table()
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
def build_k28(*, blasts: Literal[0, 1, 2] = 1, through: Stage = "partitioned") -> Case:
    """
    The ``k=2.8`` saddle, its zone trimmed at ``f(q0)`` and blasted.

    Args:
        blasts: Single-iteration blasts of the zone (1 or 2 for the law cases).
        through: The last stage to run; stages before ``bridges`` are not
            offered (nothing uses them).

    Returns:
        The case, the pip ``f(q0)`` restored after every blast.

    Raises:
        ValueError: ``through`` is earlier than ``bridges``.
    """
    if not _reached(through, "bridges"):
        raise ValueError("build_k28 starts at the 'bridges' stage")
    session = TangleSession(
        henon_map(*HENON_K28), henon_map_inverse(*HENON_K28), henon_jacobian(*HENON_K28)
    )
    session.workbench._man_machine.area_cutoff = 1e-7
    fp = session.construct_fixed_point(list(K28_SADDLE_SEED))
    session.orient_eigenvectors(fp, _orientation(PERIOD1_ORIENTATION))
    session.initialize_both_manifolds(fp)
    session.grow_n_times(fp, "unstable", num_iterations=10)
    session.grow_until_turnaround(fp, "stable")
    session.compute_intersections([fp])
    session.trim_stable_manifolds(fp)
    session.create_bridges(fp)
    session.infer_iterate_table()
    name = {0: "k28_unblasted", 1: "k28_one_blast", 2: "k28_two_blasts"}.get(
        blasts, f"k28_{blasts}_blasts"
    )
    case = Case(
        name=name,
        session=session,
        fixed_points=[fp],
        pips=[None],
        expect=CaseExpect(periods=(1,), blasts=blasts),
        stage=through,
    )
    if not _reached(through, "pips"):
        return case

    session.classify_strong_pips()
    trellis = session.trellis(fp)
    pip = trellis.iterate(trellis.strong_pip, 1)
    assert pip is not None, "the default strong pip must have a registered image"
    case.pips = [pip]
    if not _reached(through, "zones"):
        return case

    zone = session.resonance_zone(pip)
    case.zones = [zone]
    session.classify_strong_pips()
    session.set_strong_pip(fp, pip)
    for _ in range(blasts):
        session.blast_zone(zone, num_iterations=1, fixed_point=[fp], min_separation=1e-5)
        session.classify_strong_pips()
        session.set_strong_pip(fp, pip)
    if not _reached(through, "partitioned"):
        return case
    _partition(session)
    return case


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
) -> Case:
    """
    The inner period-3 orbit alone, partitioned.

    Delegates to :func:`tanglepack.examples.henon_cases.build_period3` with
    every keyword explicit.

    Args:
        unstable_steps: Unstable growth steps (13; 15 and 16 are the known
            open deep runs).
        stable_steps: Stable growth steps.
        area_cutoff: The manifold machine's refinement area cutoff.
        blasts: Single-iteration blasts of the zone.
        min_separation: The blasts' ``min_separation``.

    Returns:
        The case.
    """
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
        expect=CaseExpect(periods=(3,), blasts=blasts, has_image_bridges=False),
        zones=list(build.session.resonance_zones.values()),
    )


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
        expect=CaseExpect(periods=(1, 3), n_tangles=2),
        stage=through,
    )
    session.workbench._man_machine.area_cutoff = p3_area_cutoff
    fp3 = session.construct_fixed_point([list(p) for p in P3_ORBIT_SEED])
    session.orient_eigenvectors(fp3, _orientation(PERIOD3_ORIENTATION))
    session.initialize_both_manifolds(fp3)
    session.grow_n_times(fp3, "unstable", num_iterations=p3_unstable_steps)
    session.grow_n_times(fp3, "stable", num_iterations=p3_stable_steps)

    session.workbench._man_machine.area_cutoff = p1_area_cutoff
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
    session.infer_iterate_table()
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
    return Case(
        name="nested" if outer_blasts == 2 and not inner_blasts else "nested_variant",
        session=build.session,
        fixed_points=list(build.fixed_points),
        pips=list(build.pips),
        expect=CaseExpect(
            periods=(1, 3),
            n_tangles=2,
            blasts=outer_blasts,
            has_image_bridges=True if outer_blasts == 2 and not inner_blasts else None,
        ),
        zones=list(build.session.resonance_zones.values()),
    )


# --------------------------------------------------------------------------- #
# The k = 10 inversion saddle and the b = -1 placeholder
# --------------------------------------------------------------------------- #
def build_inversion(
    *,
    unstable_steps: int = 6,
    stable_steps: int = 5,
    through: Stage = "partitioned",
) -> Case:
    """
    The ``k=10`` map's inversion saddle at ``(-2.3166, 2.3166)``.

    Both eigenvalues are negative (-4.406 and -0.227), so ``k_value = 2 *
    period`` and the two branches of each manifold are one chain. det J = +1:
    the map is orientation preserving. ``orient_eigenvectors`` is a no-op here
    (the second branch is the image of the first), so it is not called.

    Args:
        unstable_steps: Unstable growth steps (6; 7 does not finish, P2).
        stable_steps: Stable growth steps.
        through: The last stage to run (see :data:`Stage`). ``zones`` adds
            nothing.

    Returns:
        The case.
    """
    session = TangleSession(
        henon_map(*HENON_K10), henon_map_inverse(*HENON_K10), henon_jacobian(*HENON_K10)
    )
    case = Case(
        name="inversion",
        session=session,
        fixed_points=[],
        pips=[None],
        expect=CaseExpect(periods=(1,), inversion=True),
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
    session.infer_iterate_table()
    if not _reached(through, "pips"):
        return case
    session.classify_strong_pips()
    case.pips = [session.strong_pip(fp)]
    if not _reached(through, "partitioned"):
        return case
    _partition(session)
    return case


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
    session = TangleSession(
        henon_map(*HENON_ORIENTATION_REVERSING),
        henon_map_inverse(*HENON_ORIENTATION_REVERSING),
        henon_jacobian(*HENON_ORIENTATION_REVERSING),
    )
    fp = session.construct_fixed_point(list(ORIENTATION_REVERSING_SEED))
    case = Case(
        name="orientation_reversing",
        session=session,
        fixed_points=[fp],
        pips=[None],
        expect=CaseExpect(periods=(1,)),
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
    session.infer_iterate_table()
    if not _reached(through, "pips"):
        return case
    session.classify_strong_pips()
    case.pips = [session.strong_pip(fp)]
    if _reached(through, "partitioned"):
        _partition(session)
    return case


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
    "inversion",
)

#: Every named case and its frozen builder.
BUILDERS: Mapping[str, Callable[[], Case]] = MappingProxyType(
    {
        "k10": build_k10,
        "k28_one_blast": lambda: build_k28(blasts=1),
        "k28_two_blasts": lambda: build_k28(blasts=2),
        "p3": build_period3,
        "nested": build_nested,
        "nested_unblasted": lambda: build_nested(
            outer_blasts=0, p1_area_cutoff=1e-4, through="zones"
        ),
        "inversion": build_inversion,
        "orientation_reversing": build_orientation_reversing,
        "p3_15": lambda: build_period3(unstable_steps=15),
        "p3_16": lambda: build_period3(unstable_steps=16),
        "p3_6_blasts": lambda: build_period3(blasts=6),
    }
)


@dataclass(frozen=True)
class KnownIssue:
    """
    A known bug recorded as ``xfail(strict=True)``.

    Attributes:
        reason: Why the pair fails today.
        raises: The exception type the failure raises, when it is one.
    """

    reason: str
    raises: Optional[type[BaseException]] = None


#: ``(law_id, case) -> KnownIssue``; fixing these is out of scope (decision 10).
#: The law ids are the ``helpers.laws`` check names without ``check_``. A
#: law-tier pair is left out of its layer's run and gets a test of its own
#: (``helpers.law_tier.known_issue_test``), so a fix flips loudly to XPASS.
KNOWN_ISSUES: Mapping[tuple[str, str], KnownIssue] = MappingProxyType(
    {
        ("one_anchor_per_unstable_branch", "inversion"): KnownIssue(
            "an inversion point registers one (0,0) anchor per (unstable, stable) "
            "branch pair, i.e. two per unstable branch; the author's rule is one (P2)"
        ),
        ("anchor_bridge_class_leads_its_tangle", "inversion"): KnownIssue(
            "with two anchors per unstable branch, the bridge leaving the second "
            "anchor sits in an inert class (Phase 4 finding, tied to the anchor issue)"
        ),
        ("arrangement_euler", "inversion"): KnownIssue(
            "the inversion arrangement has V - E + F = 0 on its one component (four "
            "coincident anchors; Phase 4 finding)"
        ),
        ("itinerary_pairs_same_side", "inversion"): KnownIssue(
            "an inversion class's walked itinerary pairs a left element of branch 0.1 "
            "with a right element of branch 0.0 (a cross-side pair; Phase 4 finding)"
        ),
        ("landings_contained", "inversion"): KnownIssue(
            "an inversion class's one-step landing straddles a cut of the image "
            "branch (contained=False; Phase 4 finding)"
        ),
        ("arrangement_regions_disjoint", "k28_one_blast"): KnownIssue(
            "a zero-area region bounded by an unstable arc no Bridge object spans "
            "(after a blast) reports a representative point inside a neighbouring "
            "region (Phase 4 finding)"
        ),
        ("arrangement_regions_disjoint", "k28_two_blasts"): KnownIssue(
            "a zero-area region bounded by an unstable arc no Bridge object spans "
            "(after a blast) reports a representative point inside a neighbouring "
            "region (Phase 4 finding)"
        ),
        ("orientation_reversing_case_builds", "orientation_reversing"): KnownIssue(
            "det J < 0 is not supported: FixedPoint.set_k_value rejects a saddle "
            "with exactly one negative eigenvalue (P4)",
            raises=ValueError,
        ),
        ("open_p3_deep_runs", "p3_15"): KnownIssue(
            "15 unstable steps: an unreachable class and a virtual new1 (P3)"
        ),
        ("open_p3_deep_runs", "p3_16"): KnownIssue(
            "16 unstable steps: an unreachable class and a virtual new1 (P3)"
        ),
        ("open_p3_deep_runs", "p3_6_blasts"): KnownIssue(
            "6 blasts: an unreachable class and a virtual new1 (P3)"
        ),
    }
)

_NO_HOLES = "no pseudoneighbors, hence no holes, at the only feasible depth (P2)"
_NO_IMAGE_BRIDGES = "the minimal trellis maps no hole bridge forward (no image bridges)"
_NO_REFINEMENT = "no class splits (no refined children)"

#: ``(law_id, case) -> reason``: the law has nothing to check on that case; the
#: layer test skips it (and does not count it) on that case.
NOT_APPLICABLE: Mapping[tuple[str, str], str] = MappingProxyType(
    {
        **{
            (law, "inversion"): _NO_HOLES
            for law in (
                "holes_are_classified",
                "backward_holes_keep_bridge_side",
                "direct_hole_side_is_coordinate_side",
                "direct_hole_opens_inward_pair",
                "openings_on_own_bridge_row",
                "openings_missing_only_at_anchor_or_tail",
                "propagated_holes_land_on_predicted_branch",
                "propagation_terminates",
                "no_direct_hole_beyond_fundamental",
                "reference_pairs_valid",
                "partition_singletons",
                "minimal_trellis_bridges",
                "iterated_cut_provenance",
                "refined_children_inherit_word",
                "member_matching_consistent",
            )
        },
        ("dual_face_side_is_geometry", "inversion"): (
            "the minimal trellis keeps no bridge, so its arrangement closes no region (P2)"
        ),
        ("arrangement_image_of", "inversion"): "image_of closes no region (P2)",
        ("arrangement_preimage_inverts_image", "inversion"): "image_of closes no region (P2)",
        ("arrangement_image_of", "k10"): (
            "the one face carrying a region's mapped corners is a proper sub-face, "
            "which image_of rejects; no region has a closed image"
        ),
        ("arrangement_preimage_inverts_image", "k10"): "no region has a closed image",
        ("partition_singletons", "k10"): "no hole pair pinches a crossing",
        ("partition_singletons", "k28_two_blasts"): "no hole pair pinches a crossing",
        ("iterated_cut_provenance", "k28_one_blast"): _NO_IMAGE_BRIDGES,
        ("iterated_cut_provenance", "p3"): _NO_IMAGE_BRIDGES,
        ("refined_children_inherit_word", "k28_one_blast"): _NO_REFINEMENT,
        ("refined_children_inherit_word", "p3"): _NO_REFINEMENT,
        ("member_matching_consistent", "k28_one_blast"): _NO_REFINEMENT,
        ("member_matching_consistent", "p3"): _NO_REFINEMENT,
    }
)


def issue_marks(law_id: str, case: str) -> list:
    """
    The pytest marks one ``(law, case)`` pair carries.

    Args:
        law_id: The law's identifier.
        case: The case name.

    Returns:
        ``xfail(strict=True)`` (with ``raises`` when given) for a known issue,
        a ``skip`` for a not-applicable pair, nothing otherwise.
    """
    marks = []
    issue = KNOWN_ISSUES.get((law_id, case))
    if issue is not None:
        kwargs = {"strict": True, "reason": issue.reason}
        if issue.raises is not None:
            kwargs["raises"] = issue.raises
        marks.append(pytest.mark.xfail(**kwargs))
    reason = NOT_APPLICABLE.get((law_id, case))
    if reason is not None:
        marks.append(pytest.mark.skip(reason=f"not applicable: {reason}"))
    return marks

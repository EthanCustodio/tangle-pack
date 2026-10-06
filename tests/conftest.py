"""Shared fixtures for the tanglepack test suite.

Every fixture here is a THIN, FUNCTION-SCOPED call to a builder of the
tests-only :mod:`cases` module, where every case's parameters are frozen
(author decision 7, 2026-10-05): a test gets a fresh build it may mutate, and
no fixture shares state between tests. The one exception is the physical-law
tier (``tests/invariants/conftest.py``), which builds one read-only case per
module (``law_case``) behind a fingerprint guard.

Fixtures:

* ``henon_map`` / ``henon_map_inverse`` / ``henon_jacobian`` -- the ``k=10,
  b=1`` binary-horseshoe map (batch capable, coordinate on axis 0).
* k=10 numerics stages, each ``(workbench, fp)`` unless stated:
  ``workbench`` (empty), ``fixed_point``, ``initialized``, ``grown_unstable``
  (``(workbench, fp, unstable manifold)``, 7 unstable steps), ``grown_both``
  (7 steps + stable turnaround), ``small_tangle`` (+ crossings),
  ``henon_tangle_with_bridges`` (9 steps, trimmed, bridges, iterate table;
  Phase 0 probe P6 showed the old 7 + 2 recipe without the table is
  interchangeable).
* ``k10_session`` / ``k10_partitioned`` -- ``(session, fp)`` at the bridges
  and the partitioned stage.
* ``henon_inversion_initialized`` / ``henon_inversion`` -- the inversion saddle
  at ``(-2.3166, 2.3166)``, seeded / grown 6 + 5 steps with crossings.
* ``henon_p3_session`` -- the unblasted nested tangle stopped at its resonance
  zones (outer saddle at ``area_cutoff = 1e-4``), ``(session, fp3, fp1,
  largest zone)``; ``p3_partitioned`` -- the same map partitioned through
  :func:`cases.build_nested` (no blasts, outer saddle at ``area_cutoff =
  1e-4``), ``(session, fp3, fp1)``.
* ``k28_partitioned`` / ``k28_two_blasts_partitioned`` -- ``(session, fp)``,
  the k=2.8 zone blasted once / twice, partitioned.

Registry ids are not reproducible between builds (even of one case in one
process): locate crossings by cdist order or structure, never by id.

Layout (tiers; every test basename is unique, no ``__init__.py``):

* ``invariants/`` -- the physical-law tier: one test per (law, case) over
  :data:`cases.LAW_CASES`, a module-scoped read-only build per case behind a
  fingerprint guard; the checks live in ``helpers/laws.py``.
* ``unit/numerics``, ``unit/topology``, ``unit/loom`` -- rules and kernels of
  each layer, public API first (private access only for the pure kernels the
  2026-10-05 policy allows).
* ``facade/`` -- ``TangleSession``: the one cache-contract table, session =
  direct build, fan-outs, delegates and reports.
* ``golden/`` -- the CLAUDE.md fixture facts, each pinned ONCE
  (``@pytest.mark.golden``, runs by default; author-gated, no re-record).
* ``regression/`` -- fixed bugs that need their own build, and the open deep
  period-3 runs.
* ``plotting/`` -- smoke plus topological properties of the drawings.
* ``helpers/`` (checks, fakes, letter-free name spellers, log assertion),
  ``cases.py`` (frozen builders, ``KNOWN_ISSUES``, ``NOT_APPLICABLE``) and
  ``_tools/coverage_guard.py`` are not collected.

Markers: ``golden`` (runs by default; ``pytest -m golden`` selects the golden
tier) and ``perf`` (opt-in, ``pytest -m perf``; ``addopts`` deselects it).

Dev Notes:
    Provisional rules, tested only in their firm part (author decision 3,
    2026-10-05). Revisit these tests when the author settles a rule:

    * The ``+1 .. +(k_value - 1)`` forward-hole exemption (holes propagate
      backward only; ``_is_forward_beyond_fundamental``). Tested: a recorded
      pair at forward iterate ``>= k_value`` is never punched
      (``invariants/test_law_partition.py``,
      ``unit/topology/test_stable_partition.py``). Untested: that the
      ``+1 .. +(k-1)`` iterates ARE punched.
    * The refined-point cdist. Tested: a refined node's cdist lies strictly
      between its neighbours' (``unit/numerics/test_refinement.py``).
      Untested: that it is the mean of the neighbours.
    * The anchor faces-closed limitation of ``grow_until_faces_closed`` (the
      anchor's faces never close in the arrangement). Untested; noted in
      ``unit/numerics/test_grow_until.py``.

    Dead API: the code stays, its tests were dropped (author decision 4), and
    ``_tools/coverage_guard.py`` whitelists it:
    ``IntersectionRegistry.on_interval`` and ``_get_lambda_u``,
    ``TangleSession.invalidate_trellises`` (deprecated no-op), the scalar
    ``ManifoldMachine._curvature_area`` (superseded by
    ``_curvature_area_batch``), the ``Pseudoneighbor.forward_unstable_branch_cycle``
    / ``StrongPip.forward_stable_branch_cycle`` wrappers (use
    ``FixedPoint.branch_cycle``), and ``BoundaryArc`` (already gone from
    ``src/``).

    KNOWN_ISSUES index (:data:`cases.KNOWN_ISSUES`, every one
    ``xfail(strict=True)``; an XPASS fails the run, so a fix must remove its
    entry). Fixing them is a source change, out of the test refactor's scope:

    * ``one_anchor_per_unstable_branch`` on ``inversion``: the k=10 inversion
      saddle registers one (0,0) anchor per (unstable, stable) branch pair,
      i.e. two per unstable branch; the author's rule is one.
    * Probably downstream of it, all on ``inversion``:
      ``anchor_bridge_class_leads_its_tangle``, ``arrangement_euler``,
      ``itinerary_pairs_same_side``, ``landings_contained``.
    * ``arrangement_regions_disjoint`` on ``k28_one_blast`` and
      ``k28_two_blasts``: a zero-area region bounded by an unstable arc no
      ``Bridge`` spans reports a representative point inside its neighbour.
    * ``orientation_reversing_case_builds`` (``raises=ValueError``): the b = -1
      placeholder; ``det J < 0`` is not supported yet. (The inversion saddle
      has ``det J = +1``; it is NOT orientation reversing.)
    * ``open_p3_deep_runs`` on ``p3_15``, ``p3_16``, ``p3_6_blasts``
      (``regression/test_open_p3_deep_runs.py``): an unreachable class and a
      virtual ``new1``.

    Not-applicable pairs (:data:`cases.NOT_APPLICABLE`) are skips with a
    reason, never xfails. The running record of every deleted test and every
    deviation is ``docs/test_suite_refactor_ledger.md``.
"""

from __future__ import annotations

from typing import Callable

import matplotlib

matplotlib.use("Agg")

import numpy as np  # noqa: E402
import pytest  # noqa: E402

from tanglepack import TangleSession, TangleWorkbench  # noqa: E402
from tanglepack.examples import (  # noqa: E402
    HENON_K10,
    henon_jacobian as _henon_jacobian_factory,
    henon_map as _henon_map_factory,
    henon_map_inverse as _henon_map_inverse_factory,
)

from cases import (  # noqa: E402
    build_inversion,
    build_k10,
    build_k28,
    build_nested,
)

# --------------------------------------------------------------------------- #
# Maps
# --------------------------------------------------------------------------- #
_henon_map = _henon_map_factory(*HENON_K10)
_henon_map_inverse = _henon_map_inverse_factory(*HENON_K10)
_henon_jacobian = _henon_jacobian_factory(*HENON_K10)


@pytest.fixture(name="henon_map")
def _k10_map_fixture() -> Callable[[np.ndarray], np.ndarray]:
    """The ``k=10, b=1`` Hénon map."""
    return _henon_map


@pytest.fixture(name="henon_map_inverse")
def _k10_map_inverse_fixture() -> Callable[[np.ndarray], np.ndarray]:
    """The inverse of the ``k=10, b=1`` Hénon map."""
    return _henon_map_inverse


@pytest.fixture(name="henon_jacobian")
def _k10_jacobian_fixture() -> Callable[[np.ndarray], np.ndarray]:
    """The analytic Jacobian of the ``k=10, b=1`` Hénon map."""
    return _henon_jacobian


def _unstable_manifold(workbench: TangleWorkbench, fp, branch_index: int = 0):
    """The orbit-0 unstable manifold on one branch."""
    return workbench.manifolds[(fp, "unstable", 0, branch_index)]


# --------------------------------------------------------------------------- #
# k = 10: the numerics stages
# --------------------------------------------------------------------------- #
@pytest.fixture
def workbench() -> TangleWorkbench:
    """A fresh k=10 workbench with no fixed points yet."""
    return build_k10(through="empty").workbench


@pytest.fixture
def fixed_point() -> tuple[TangleWorkbench, object]:
    """``(workbench, fp)`` with the k=10 saddle constructed and oriented."""
    case = build_k10(through="fixed_point")
    return case.workbench, case.fixed_point


@pytest.fixture
def initialized() -> tuple[TangleWorkbench, object]:
    """``(workbench, fp)`` with both fundamental segments initialized."""
    case = build_k10(through="seeded")
    return case.workbench, case.fixed_point


@pytest.fixture
def grown_unstable() -> tuple[TangleWorkbench, object, object]:
    """``(workbench, fp, manifold)`` with the unstable manifold grown 7 times.

    7 iterations develop a small tangle that crosses the stable manifold (a
    handful of crossings) while staying snappy.
    """
    case = build_k10(unstable_steps=7, through="grown_unstable")
    return case.workbench, case.fixed_point, _unstable_manifold(case.workbench, case.fixed_point)


@pytest.fixture
def grown_both() -> tuple[TangleWorkbench, object]:
    """``(workbench, fp)``: unstable grown 7 times, stable grown to turnaround."""
    case = build_k10(unstable_steps=7, through="grown")
    return case.workbench, case.fixed_point


@pytest.fixture
def small_tangle() -> tuple[TangleWorkbench, object]:
    """``(workbench, fp)`` after the 7-step tangle's crossings are computed."""
    case = build_k10(unstable_steps=7, through="intersected")
    return case.workbench, case.fixed_point


@pytest.fixture
def henon_tangle_with_bridges() -> tuple[TangleWorkbench, object]:
    """``(workbench, fp)``: 9 unstable steps, crossings, trimmed stable manifold,
    bridges and the iterate table (enough crossings for a full pseudoneighbor
    reference window)."""
    case = build_k10(through="bridges")
    return case.workbench, case.fixed_point


# --------------------------------------------------------------------------- #
# k = 10: sessions
# --------------------------------------------------------------------------- #
@pytest.fixture
def k10_session() -> tuple[TangleSession, object]:
    """``(session, fp)``: the k=10 tangle with bridges and the iterate table."""
    case = build_k10(through="bridges")
    return case.session, case.fixed_point


@pytest.fixture
def k10_partitioned() -> tuple[TangleSession, object]:
    """``(session, fp)`` with the k=10 trellis classified, punched and partitioned."""
    case = build_k10()
    assert case.session.trellis(case.fixed_point).stable_partitions
    return case.session, case.fixed_point


# --------------------------------------------------------------------------- #
# The k = 10 inversion saddle (k_value == 2 * period, det J = +1)
# --------------------------------------------------------------------------- #
@pytest.fixture
def henon_inversion_initialized() -> tuple[TangleWorkbench, object]:
    """``(workbench, fp)`` with the inversion saddle's fundamental segments only."""
    case = build_inversion(through="seeded")
    return case.workbench, case.fixed_point


@pytest.fixture
def henon_inversion() -> tuple[TangleWorkbench, object]:
    """``(workbench, fp)``: the inversion saddle grown 6 unstable / 5 stable
    steps on both branches, crossings computed."""
    case = build_inversion(through="intersected")
    return case.workbench, case.fixed_point


# --------------------------------------------------------------------------- #
# The nested map (2, 1)
# --------------------------------------------------------------------------- #
@pytest.fixture
def henon_p3_session() -> tuple[TangleSession, object, object, object]:
    """``(session, fp3, fp1, zone)``: the unblasted nested tangle and its zones.

    The outer saddle is grown at ``area_cutoff = 1e-4``; the build stops right
    after ``add_resonance_zones`` (no repin, no partition). ``zone`` is the
    LARGEST resonance zone (the historical fixture's choice).
    """
    case = build_nested(outer_blasts=0, p1_area_cutoff=1e-4, through="zones")
    fp1, fp3 = case.fixed_points
    zone = max(case.session.resonance_zones.values(), key=lambda z: z.area)
    return case.session, fp3, fp1, zone


@pytest.fixture
def p3_partitioned() -> tuple[TangleSession, object, object]:
    """``(session, fp3, fp1)``: the unblasted nested tangle, every trellis partitioned."""
    case = build_nested(outer_blasts=0, p1_area_cutoff=1e-4)
    fp1, fp3 = case.fixed_points
    assert case.session.trellis(fp3).stable_partitions
    return case.session, fp3, fp1


# --------------------------------------------------------------------------- #
# k = 2.8, blasted
# --------------------------------------------------------------------------- #
@pytest.fixture
def k28_partitioned() -> tuple[TangleSession, object]:
    """``(session, fp)``: the k=2.8 zone (trimmed at ``f(q0)``) blasted ONCE.

    The blast registers the forward images of one reference pair as a
    blast-child bridge: the fixture for "holes propagate backward only".
    """
    case = build_k28(blasts=1)
    assert case.session.trellis(case.fixed_point).stable_partitions
    return case.session, case.fixed_point


@pytest.fixture
def k28_two_blasts_partitioned() -> tuple[TangleSession, object]:
    """``(session, fp)``: the k=2.8 zone blasted TWICE, the symbolic-dynamics
    fixture with a non-trivial transition matrix."""
    case = build_k28(blasts=2)
    assert case.session.trellis(case.fixed_point).stable_partitions
    return case.session, case.fixed_point

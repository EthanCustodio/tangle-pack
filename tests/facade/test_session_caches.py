"""The :class:`~tanglepack.TangleSession` cache contract, as one literal table.

Every cached product of the session is checked against every event that may
or may not invalidate it (author decision 6, 2026-10-05). The expected outcome
of each cell is written out in :data:`EXPECTED`: ``HIT`` means the second call
hands back the very object of the first, ``MISS`` means it builds a new one.

Products (each asked for with the default, all-fixed-points selection):
``trellis``, ``arrangement``, ``bridge_classes``, ``minimal_trellis``,
``iterated_partition``, ``dual_graph``, ``symbolic_dynamics``.

Events, each applied to a fresh ``build_k10()`` between the two calls:

* ``hit`` -- nothing but reads (bridges, registry orderings, partitions).
* ``rebuild`` -- the second call passes ``rebuild=True``. Only the product's own
  rebuild is asserted; what it does to the products it reads is not part of
  the contract.
* ``generation_bump`` -- one ``iterate_bridge`` (a cheap workbench mutation that
  registers no new crossing), then the partition is redone through the session
  fan-outs. The partition signature and the strong pip come back unchanged, so
  the generation alone is what the cache sees.
* ``partition_signature`` -- a re-partition at the SAME generation (the
  trellis's results cleared and redone against another strong pip, then the
  original pip restored), so only the partition signature moves.
* ``pip_change`` -- ``set_strong_pip`` to another candidate, nothing else.

Extra rows pin that every workbench mutation path drops the cached trellis
(growth, ``iterate_bridge``, ``rebuild_bridges``, ``add_resonance_zones``, a
zone's ``restore`` and ``compute_intersections`` on the two-fixed-point nested
session), and that the cache is kept per fixed-point selection.

Whether a session result EQUALS a direct build is
``test_session_equivalence.py``; nothing here compares contents.
"""

from __future__ import annotations

from typing import Any, Callable

import pytest

from cases import Case, build_k10, build_nested
from tanglepack import TangleSession

HIT, MISS = True, False

PRODUCTS: tuple[str, ...] = (
    "trellis",
    "arrangement",
    "bridge_classes",
    "minimal_trellis",
    "iterated_partition",
    "dual_graph",
    "symbolic_dynamics",
)
EVENTS: tuple[str, ...] = (
    "hit",
    "rebuild",
    "generation_bump",
    "partition_signature",
    "pip_change",
)

#: product -> (hit, rebuild, generation_bump, partition_signature, pip_change).
#: The trellis and its arrangement are keyed on the generation alone; the
#: classes, minimal trellis and iterated partition on (generation, partition
#: signature); the dual graph and the dynamics on (generation, signature,
#: gathered strong pips).
EXPECTED: dict[str, tuple[bool, bool, bool, bool, bool]] = {
    "trellis":            (HIT, MISS, MISS, HIT,  HIT),
    "arrangement":        (HIT, MISS, MISS, HIT,  HIT),
    "bridge_classes":     (HIT, MISS, MISS, MISS, HIT),
    "minimal_trellis":    (HIT, MISS, MISS, MISS, HIT),
    "iterated_partition": (HIT, MISS, MISS, MISS, HIT),
    "dual_graph":         (HIT, MISS, MISS, MISS, MISS),
    "symbolic_dynamics":  (HIT, MISS, MISS, MISS, MISS),
}


# --------------------------------------------------------------------------- #
# helpers
# --------------------------------------------------------------------------- #
def _repartition(session: TangleSession) -> None:
    """Redo the pips, pseudoneighbors, holes and partition through the fan-outs."""
    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()


def _alternate_pip(session: TangleSession, fp: Any) -> int:
    """A strong-pip candidate of ``fp`` other than the chosen one."""
    trellis = session.trellis(fp)
    alternates = [c for c in trellis.strong_pip_candidates if c != trellis.strong_pip]
    assert alternates, "the k10 case must offer more than one strong-pip candidate"
    return alternates[0]


def _reads(case: Case) -> None:
    """Read-only accesses that must never invalidate a cache."""
    session, fp = case.session, case.fixed_point
    session.workbench.bridges
    case.registry.by_unstable_cdist
    session.trellis(fp).stable_partitions
    session.homotopy_partition()


def _bump_generation(case: Case) -> None:
    """One ``iterate_bridge``, then the partition redone; signature and pip kept."""
    session, fp = case.session, case.fixed_point
    generation = session.workbench.generation
    signature = session.homotopy_partition().signature()
    pip = session.strong_pip(fp)

    session.workbench.iterate_bridge(session.workbench.uniiterated_bridges[0])
    _repartition(session)

    assert session.workbench.generation != generation
    assert session.homotopy_partition().signature() == signature
    assert session.strong_pip(fp) == pip


def _change_partition_signature(case: Case) -> None:
    """Re-partition against another pip at the same generation, then restore the pip."""
    session, fp = case.session, case.fixed_point
    generation = session.workbench.generation
    signature = session.homotopy_partition().signature()
    trellis = session.trellis(fp)
    original = trellis.strong_pip

    trellis.clear_results()
    trellis.classify_strong_pips()
    trellis.set_strong_pip(_alternate_pip(session, fp))
    session.compute_pseudoneighbors(fp)
    session.punch_holes(fp)
    session.partition_stable_manifold(fp)
    trellis.set_strong_pip(original)

    assert session.workbench.generation == generation
    assert session.homotopy_partition().signature() != signature
    assert session.strong_pip(fp) == original


def _change_pip(case: Case) -> None:
    """Choose another strong pip; generation and partition stay put."""
    session, fp = case.session, case.fixed_point
    generation = session.workbench.generation
    signature = session.homotopy_partition().signature()

    session.trellis(fp).set_strong_pip(_alternate_pip(session, fp))

    assert session.workbench.generation == generation
    assert session.homotopy_partition().signature() == signature


_EVENT_ACTIONS: dict[str, Callable[[Case], None]] = {
    "hit": _reads,
    "generation_bump": _bump_generation,
    "partition_signature": _change_partition_signature,
    "pip_change": _change_pip,
}

_CELLS = [
    pytest.param(product, event, expected, id=f"{product}-{event}")
    for product, row in EXPECTED.items()
    for event, expected in zip(EVENTS, row)
]


# --------------------------------------------------------------------------- #
# the table
# --------------------------------------------------------------------------- #
def test_expected_table_is_complete() -> None:
    """Every product has one literal outcome per event."""
    assert tuple(EXPECTED) == PRODUCTS
    assert all(len(row) == len(EVENTS) for row in EXPECTED.values())


@pytest.mark.parametrize(("product", "event", "expected_hit"), _CELLS)
def test_session_cache_table(product: str, event: str, expected_hit: bool) -> None:
    """The second call hits or misses exactly as :data:`EXPECTED` says."""
    case = build_k10()
    accessor = getattr(case.session, product)
    first = accessor()

    if event == "rebuild":
        second = accessor(rebuild=True)
    else:
        _EVENT_ACTIONS[event](case)
        second = accessor()

    assert (second is first) is expected_hit
    if expected_hit:
        assert accessor() is first


# --------------------------------------------------------------------------- #
# every workbench mutation path drops the cached trellis
# --------------------------------------------------------------------------- #
def _define_zone(case: Case) -> Any:
    """A resonance zone trimmed at a strong-pip candidate inside the stable extent."""
    session, fp = case.session, case.fixed_point
    registry = case.registry
    trellis = session.trellis(fp)
    trellis.classify_strong_pips()
    outermost = max(registry[i].stable_cdist for i in registry.all_ids())
    inner = [c for c in trellis.strong_pip_candidates if registry[c].stable_cdist < outermost]
    assert inner, "expected a strong-pip candidate inside the stable extent"
    pip = max(inner, key=lambda c: registry[c].stable_cdist)
    trellis.set_strong_pip(pip)
    session.add_resonance_zones([pip])
    return session.resonance_zones[(fp, 0)]


def _grow(case: Case) -> None:
    case.session.grow_n_times(case.fixed_point, "unstable", num_iterations=1)


def _iterate_bridge(case: Case) -> None:
    workbench = case.session.workbench
    workbench.iterate_bridge(workbench.uniiterated_bridges[0])


def _rebuild_bridges(case: Case) -> None:
    case.session.workbench.rebuild_bridges(case.fixed_point)


def _add_resonance_zones(case: Case) -> None:
    case.session.add_resonance_zones([case.registry.all_ids()[0]])


def _restore(case: Case) -> None:
    case.zones[-1].restore(case.session.workbench)


_MUTATIONS: dict[str, Callable[[Case], None]] = {
    "growth": _grow,
    "iterate_bridge": _iterate_bridge,
    "rebuild_bridges": _rebuild_bridges,
    "add_resonance_zones": _add_resonance_zones,
    "restore": _restore,
}
#: Mutations that need a resonance zone made before the trellis is cached.
_NEEDS_ZONE = frozenset({"add_resonance_zones", "restore"})


@pytest.mark.parametrize("mutation", list(_MUTATIONS))
def test_trellis_misses_after_every_mutation_path(mutation: str) -> None:
    """Growth, a re-cut, an iterated bridge, a zone and its restore all miss.

    ``restore`` is the regression guard: it is the gap the manual
    ``invalidate_trellises`` calls used to leave open (stale trellises after a
    blast or restore, 2026-07).
    """
    case = build_k10(through="bridges")
    if mutation in _NEEDS_ZONE:
        case.zones.append(_define_zone(case))
    session, fp = case.session, case.fixed_point
    first = session.trellis(fp)

    _MUTATIONS[mutation](case)

    fresh = session.trellis(fp)
    assert fresh is not first
    assert session.trellis(fp) is fresh


def test_trellis_misses_after_a_recompute_on_the_nested_session() -> None:
    """``compute_intersections`` swaps the registry; the cached trellis follows it."""
    case = build_nested(outer_blasts=0, p1_area_cutoff=1e-4, through="zones")
    session = case.session
    fp1, fp3 = case.fixed_points
    stale = session.trellis(fp3)

    session.compute_intersections([fp3, fp1])

    fresh = session.trellis(fp3)
    assert fresh is not stale
    assert fresh.registry is session.workbench.intersection_registry


# --------------------------------------------------------------------------- #
# selection keys and the alphabet
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("product", PRODUCTS)
def test_cache_is_kept_per_fixed_point_selection(product: str) -> None:
    """The all-fixed-points and the single-fixed-point selections are two entries."""
    case = build_k10()
    accessor = getattr(case.session, product)

    whole = accessor()
    single = accessor(case.fixed_point)

    assert single is not whole
    assert accessor() is whole
    assert accessor(case.fixed_point) is single


def test_bridge_class_letters_survive_rebuilds_and_repartitions() -> None:
    """The alphabet outlives every cache entry: a rebuild keeps each letter, and
    after a re-partition every lettered class carries the alphabet's letter."""
    case = build_k10()
    session = case.session
    first = session.bridge_classes()
    alphabet_size = len(session.bridge_alphabet)

    rebuilt = session.bridge_classes(rebuild=True)
    assert rebuilt == first
    assert [e.letter for e in rebuilt] == [e.letter for e in first]
    assert len(session.bridge_alphabet) == alphabet_size

    _change_partition_signature(case)
    for entry in session.bridge_classes():
        if entry.letter is not None:
            assert session.bridge_alphabet.assigned[entry.bridge_class] == entry.letter
    assert len(session.bridge_alphabet) >= alphabet_size

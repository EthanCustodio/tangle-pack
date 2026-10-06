"""The :class:`~tanglepack.TangleSession` caches never serve a stale product.

Every cached product of the session is checked against every event that may
change what it should be (author decision 6, 2026-10-05; rewritten the same
day to stop pinning hit/miss policy). The contract is FRESHNESS: after any
event, the session's product equals -- structurally, through the letter-free
views of ``helpers/direct.py``, never by object identity -- a direct,
cache-free build from the current state. Whether a given event is served from
the cache or rebuilt is the session's business, with two exceptions that ARE
the contract:

* ``no_change`` -- nothing but reads between two calls: the second call hands
  back the very object of the first (the cache is a cache).
* ``rebuild`` -- ``rebuild=True`` hands back a NEW object of that product
  (what it does to the products it reads is not part of the contract).

Products (each asked for with the default, all-fixed-points selection):
``trellis``, ``arrangement``, ``bridge_classes``, ``minimal_trellis``,
``iterated_partition``, ``dual_graph``, ``symbolic_dynamics``.

Events, each applied to a fresh ``build_k10()`` after EVERY product has been
asked for once (so every cache holds an entry that the event may make stale):

* ``no_change`` -- reads only (bridges, registry orderings, partitions).
* ``rebuild`` -- the product asked for again with ``rebuild=True``.
* ``generation_bump`` -- one ``iterate_bridge`` (a cheap workbench mutation),
  then the partition redone through the session fan-outs. On ``k10`` it
  registers no new crossing, so every view stays put: the cell checks that a
  generation move leaves the session consistent, not that it is noticed.
* ``partition_signature`` -- a re-partition at the SAME generation (the
  trellis's results cleared and redone against another strong pip, then the
  original pip restored). It changes every product from the bridge classes on.
* ``pip_change`` -- ``set_strong_pip`` to another candidate, nothing else. It
  changes the dual graph and the symbolic dynamics.

Further tests check the trellis is fresh after every workbench mutation path
(growth, ``iterate_bridge``, ``rebuild_bridges``, ``add_resonance_zones``, a
zone's ``restore``, and ``compute_intersections`` on the two-fixed-point nested
session), that each fixed-point selection gets its own fresh product, and that
bridge-class letters outlive rebuilds and re-partitions.
"""

from __future__ import annotations

from typing import Any, Callable

import pytest

from cases import Case, build_k10, build_nested
from helpers.direct import direct_view, trellis_view, view_of
from helpers.zones import define_inner_zone
from tanglepack import TangleSession
from tanglepack.topology.Trellis import Trellis

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
    "no_change",
    "rebuild",
    "generation_bump",
    "partition_signature",
    "pip_change",
)


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
    """The anchor-nearest strong-pip candidate of ``fp`` other than the chosen one.

    Chosen by stable cdist, never by list position: the candidates come in
    registry-id order, which is not reproducible between builds.
    """
    trellis = session.trellis(fp)
    alternates = [c for c in trellis.strong_pip_candidates if c != trellis.strong_pip]
    assert alternates, "the k10 case must offer more than one strong-pip candidate"
    return min(alternates, key=lambda c: float(trellis.intersection(c).stable_cdist))


def _reads(case: Case) -> None:
    """Read-only accesses."""
    session, fp = case.session, case.fixed_point
    session.workbench.bridges
    case.registry.by_unstable_cdist
    session.trellis(fp).stable_partitions
    session.homotopy_partition()


def _bump_generation(case: Case) -> None:
    """One ``iterate_bridge``, then the partition redone."""
    session = case.session
    generation = session.workbench.generation

    session.workbench.iterate_bridge(session.workbench.uniiterated_bridges[0])
    _repartition(session)

    assert session.workbench.generation != generation


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
    """Choose another strong pip; nothing else moves."""
    session, fp = case.session, case.fixed_point
    pip = session.strong_pip(fp)

    session.trellis(fp).set_strong_pip(_alternate_pip(session, fp))

    assert session.strong_pip(fp) != pip


_EVENT_ACTIONS: dict[str, Callable[[Case], None]] = {
    "no_change": _reads,
    "generation_bump": _bump_generation,
    "partition_signature": _change_partition_signature,
    "pip_change": _change_pip,
}


def _assert_fresh(case: Case, product: str, result: Any, fixed_points: Any = None) -> None:
    """``result`` equals a direct, cache-free build of ``product`` from the current state."""
    assert view_of(product, result) == direct_view(case, product, fixed_points), (
        f"the session's {product} is stale"
    )


# --------------------------------------------------------------------------- #
# product x event
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("event", EVENTS)
@pytest.mark.parametrize("product", PRODUCTS)
def test_session_product_is_fresh_after_every_event(product: str, event: str) -> None:
    """After the event the product equals a direct build; it is the same object
    only after ``no_change`` and a new one only after ``rebuild=True``."""
    case = build_k10()
    for name in PRODUCTS:
        getattr(case.session, name)()
    accessor = getattr(case.session, product)
    first = accessor()

    if event == "rebuild":
        second = accessor(rebuild=True)
        assert second is not first
    else:
        _EVENT_ACTIONS[event](case)
        second = accessor()
        if event == "no_change":
            assert second is first

    _assert_fresh(case, product, second)


# --------------------------------------------------------------------------- #
# every workbench mutation path leaves the trellis fresh
# --------------------------------------------------------------------------- #
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
def test_trellis_is_fresh_after_every_mutation_path(mutation: str) -> None:
    """Growth, a re-cut, an iterated bridge, a zone and its restore all leave the
    session's trellis equal to a direct build of the workbench.

    ``restore`` is the regression guard: it is the gap the manual
    ``invalidate_trellises`` calls used to leave open (stale trellises after a
    blast or restore, 2026-07); a zone and its restore change the crossings, so
    a stale trellis would show. Growth, ``iterate_bridge`` and
    ``rebuild_bridges`` register no crossing on ``k10`` and are checked for
    consistency only.
    """
    case = build_k10(through="bridges")
    if mutation in _NEEDS_ZONE:
        case.zones.append(define_inner_zone(case.session, case.fixed_point))
    session, fp = case.session, case.fixed_point
    session.trellis(fp)

    _MUTATIONS[mutation](case)

    direct = Trellis.from_workbench(session.workbench, fp)
    assert trellis_view(session.trellis(fp)) == trellis_view(direct)


def test_trellis_is_fresh_after_a_recompute_on_the_nested_session() -> None:
    """``compute_intersections`` swaps the registry; the session's trellis follows it."""
    case = build_nested(outer_blasts=0, p1_area_cutoff=1e-4, through="zones")
    session = case.session
    fp1, fp3 = case.fixed_points
    session.trellis(fp3)

    session.compute_intersections([fp3, fp1])

    fresh = session.trellis(fp3)
    assert fresh.registry is session.workbench.intersection_registry
    assert trellis_view(fresh) == trellis_view(Trellis.from_workbench(session.workbench, fp3))


# --------------------------------------------------------------------------- #
# selection keys and the alphabet
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("product", PRODUCTS)
def test_every_fixed_point_selection_is_fresh(product: str) -> None:
    """Interleaved asks for the all-fixed-points and the single-fixed-point
    selections each return a direct build of THAT selection."""
    case = build_k10()
    accessor = getattr(case.session, product)

    accessor()
    accessor(case.fixed_point)

    _assert_fresh(case, product, accessor())
    _assert_fresh(case, product, accessor(case.fixed_point), case.fixed_point)


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

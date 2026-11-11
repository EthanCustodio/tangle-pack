"""A heteroclinic crossing as a strong pip, and the zone a ring of them bounds (2026-11-10).

The hydrogen atom in parallel fields at ``B = 3.5, E = 1`` has two period-1
saddles ``A = (v*, 0)`` and ``B = (-v*, 0)`` swapped by ``-I``, joined by a
heteroclinic tangle. ``q0 in W^u(B) x W^s(A)`` is A's pip and ``q0' in W^u(A) x
W^s(B)`` is B's; together they close ONE ring ``A -u-> q0' -s-> B -u-> q0 -s->
A``, the "diamond" around the origin. The pips are chosen here, test-side, as
the balanced crossing (``min |log(c_s / c_u)|``) of the lowest-action
heteroclinic orbit among the strong-pip candidates (the rule of
``scripts/hydrogen_pipeline.py``).

One build per module (``area_cutoff = 1e-5``, 8 unstable and 7 stable steps)
is saved once and loaded fresh by every test, so no test sees another's
mutations. No word, letter, matrix or count is pinned.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import TangleSession
from tanglepack.examples.hydrogen.hydrogen_map import HydrogenSystem
from tanglepack.loom.ResonanceZone import ring_closed, ring_cycles

ORIENTATION_A = {"unstable": np.array([-1.0, -1.5]), "stable": np.array([-1.0, 1.5])}
ORIENTATION_B = {stability: -direction for stability, direction in ORIENTATION_A.items()}
AREA_CUTOFF = 1e-5
UNSTABLE_STEPS, STABLE_STEPS = 8, 7


@pytest.fixture(scope="module")
def saved(tmp_path_factory) -> tuple[HydrogenSystem, str]:
    """The heteroclinic tangle up to its bridges, saved once: ``(system, path)``."""
    system = HydrogenSystem(3.5, 1.0)
    v_star = system.find_fixed_point()
    session = TangleSession(system.map, system.map_inverse, system.jacobian,
                            area_cutoff=AREA_CUTOFF, workers=1)
    fps = [session.construct_fixed_point(v_star), session.construct_fixed_point(-v_star)]
    for fp, orientation in zip(fps, (ORIENTATION_A, ORIENTATION_B)):
        session.orient_eigenvectors(fp, orientation)
        session.initialize_both_manifolds(fp)
        session.grow_n_times(fp, "unstable", num_iterations=UNSTABLE_STEPS)
        session.grow_n_times(fp, "stable", num_iterations=STABLE_STEPS)
    session.compute_intersections(fps)
    for fp in fps:
        session.trim_stable_manifolds(fp)
    for fp in fps:
        session.create_bridges(fp)
    path = str(tmp_path_factory.mktemp("hydrogen") / "heteroclinic")
    session.save(path)
    return system, path


@pytest.fixture
def tangle(saved) -> tuple[TangleSession, object, object, int, int]:
    """``(session, A, B, pA, pB)``: a fresh load, ``pA`` on ``W^s(A)`` and ``pB`` on ``W^s(B)``."""
    system, path = saved
    session = TangleSession.load(path, system.map, system.map_inverse, system.jacobian)
    session.workers = 1
    A, B = session.fixed_points
    return session, A, B, _balanced_pip(session, B, A), _balanced_pip(session, A, B)


def _balanced_pip(session: TangleSession, unstable_fp, stable_fp) -> int:
    """The balanced candidate of the lowest-action ``W^u(unstable_fp) x W^s(stable_fp)`` orbit."""
    registry = session.workbench.intersection_registry
    candidates = [
        iid for iid in session.trellis(session.fixed_points).strong_pip_candidates
        if registry[iid].manifold_a_key[0] is unstable_fp and registry[iid].manifold_b_key[0] is stable_fp
    ]
    action = {iid: registry[iid].stable_cdist * registry[iid].unstable_cdist for iid in candidates}
    primary = [iid for iid in candidates if action[iid] < min(action.values()) * (1 + 1e-3)]
    return min(primary, key=lambda iid: abs(np.log(registry[iid].stable_cdist / registry[iid].unstable_cdist)))


def _depth(session: TangleSession, fp) -> int:
    """``m`` of the registered orbit ``q0, f(q0), ..., f^m(q0)`` of ``fp``'s chosen pip."""
    trellis = session.trellis(fp)
    orbit = [trellis.registry.find(session.chosen_strong_pip(fp))]
    while (image := trellis.iterate(orbit[-1], 1)) is not None:
        orbit.append(image)
    return len(orbit) - 1


def test_a_heteroclinic_pip_is_the_pip_of_its_stable_fixed_point(tangle):
    session, A, B, pA, pB = tangle
    for fp, pip in ((A, pA), (B, pB)):
        assert session.set_strong_pip(fp, pip) == pip
        assert session.trellis(fp).pips == {fp: pip}
    with pytest.raises(ValueError, match="stable branch"):
        session.set_strong_pip(A, pB)


def test_a_ring_of_heteroclinic_pips_bounds_one_zone(tangle):
    session, A, B, pA, pB = tangle
    registry = session.workbench.intersection_registry
    corners = [registry[pip].coords for pip in (pA, pB)]

    session.add_resonance_zones([pA, pB])

    (zone,) = session.resonance_zones.values()
    assert set(map(id, zone.fixed_points)) == {id(A), id(B)}
    cuts = zone.boundary_intersections
    assert ring_closed(cuts) and [sorted(cycle) for cycle in ring_cycles(cuts)] == [list(range(len(cuts)))]
    midpoint = np.mean([A.coordinates[0], B.coordinates[0]], axis=0)
    assert all(zone.contains_point(point) for point in (midpoint, *corners))
    assert {fp: registry.find(session.chosen_strong_pip(fp)) for fp in (A, B)} == {A: pA, B: pB}


def test_a_lone_heteroclinic_pip_closes_no_ring(tangle):
    session, _A, _B, pA, _pB = tangle
    with pytest.raises(ValueError, match="ring"):
        session.add_resonance_zones([pA])


def test_chosen_heteroclinic_pips_survive_a_recompute_and_cut_the_topology(tangle):
    session, A, B, pA, pB = tangle
    session.set_strong_pip(A, pA)
    session.set_strong_pip(B, pB)

    session.workbench.compute_intersections([A, B], preserve_ids=True)
    session.workbench.rebuild_bridges()

    registry = session.workbench.intersection_registry
    pips = {fp: registry.find(session.chosen_strong_pip(fp)) for fp in (A, B)}
    assert pips == {A: pA, B: pB}
    assert {fp: session.strong_pip(fp) for fp in (A, B)} == pips
    assert session.topology().trellis.pips == pips


def test_the_ring_zone_topology_is_reliable_and_its_active_classes_in_the_zone(tangle):
    """Each class lies wholly inside the ring zone or wholly outside it (an escape
    lobe), and every active class -- B's bridges included -- inside: a bridge of
    either saddle belongs to the zone whose ``fixed_points`` hold both."""
    session, _A, _B, pA, pB = tangle
    session.add_resonance_zones([pA, pB])
    (zone,) = session.resonance_zones.values()

    topology = session.topology()

    assert topology.dynamics.is_reliable
    for entry in topology.classes:
        zones = {session.classify_bridge(topology.trellis.bridge_between(*member.bridge_id)) is zone
                 for member in entry.members}
        assert len(zones) == 1
        assert entry.zone_key is (zone.key if zones == {True} else None)
    assert any(not entry.inert for entry in topology.classes)
    assert all(entry.zone_key is zone.key for entry in topology.classes if not entry.inert)


def test_blasting_the_ring_keeps_the_pip_depth(tangle):
    session, A, B, pA, pB = tangle
    session.add_resonance_zones([pA, pB])
    (zone,) = session.resonance_zones.values()
    before = {fp: _depth(session, fp) for fp in (A, B)}

    result = session.blast_zone(zone, 1, fixed_point=[A, B], min_separation=1e-5)

    assert all(result.pip_depth[fp] >= before[fp] for fp in (A, B))

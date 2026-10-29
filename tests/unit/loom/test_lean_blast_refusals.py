"""A lean workbench refuses whatever would re-read the curve it freed.

A lean blast frees the interiors of the bridges nothing reads again
(``TangleWorkbench.release_curves``). From then on a read of a freed bridge,
and every operation that would re-read freed curve, raises ``CurveReleased``;
the topology and the drawing of the kept bridges still work, and so do the
joined arcs a trellis holds.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pytest

from cases import build_k28
from tanglepack.loom.Blast import _pieces
from tanglepack.numerics.Bridge import Bridge, CurveReleased
from tanglepack.topology.plotting import plot_bridges_by_class


@pytest.fixture
def lean():
    """The k28 zone, blasted lean (width guard) until a bridge is released."""
    case = build_k28(blasts=0, through="zones")
    session, zone = case.session, case.zones[0]
    for _ in range(8):
        session.blast_zone(zone, 1, fixed_point=[zone.fixed_point], min_separation="width", lean=True)
        if any(bridge.released for bridge in session.workbench.bridges):
            break
    assert session.workbench.lean
    return session, zone, zone.fixed_point


_REFUSALS = {
    "compute_intersections": lambda s, z, fp: s.compute_intersections([fp]),
    "create_bridges": lambda s, z, fp: s.create_bridges(fp),
    "rebuild_bridges": lambda s, z, fp: s.rebuild_bridges(),
    "clear_bridges": lambda s, z, fp: s.clear_bridges(),
    "grow_n_times": lambda s, z, fp: s.grow_n_times(fp, "unstable", 1),
    "grow_until_turnaround": lambda s, z, fp: s.grow_until_turnaround(fp, "stable"),
    "grow_until_bounded": lambda s, z, fp: s.grow_until_bounded(fp, "unstable", 10.0),
    "grow_until_arclength": lambda s, z, fp: s.grow_until_arclength(fp, "unstable", 10.0),
    "grow_until_intersection": lambda s, z, fp: s.grow_until_intersection(fp, "unstable"),
    "grow_until": lambda s, z, fp: s.grow_until(fp, lambda wb: True),
    "grow_until_iterates_closed": lambda s, z, fp: s.grow_until_iterates_closed(fp),
    "grow_until_faces_closed": lambda s, z, fp: s.grow_until_faces_closed(fp, [0]),
    "trim_stable_manifolds": lambda s, z, fp: s.trim_stable_manifolds(fp),
    "trim_stable_at": lambda s, z, fp: s.trim_stable_at((fp, "stable", 0, 0), 1e-3),
    "resonance_zone": lambda s, z, fp: s.resonance_zone(z.resolve_boundary_intersection_id(s.workbench)),
    "restore": lambda s, z, fp: z.restore(s.workbench),
    "session save": lambda s, z, fp: s.save("unused"),
    "workbench save": lambda s, z, fp: s.workbench.save("unused"),
    "check_resolution": lambda s, z, fp: s.check_resolution(),
    "arrangement": lambda s, z, fp: s.arrangement(),
    "blast not lean": lambda s, z, fp: s.blast_zone(z, 1, fixed_point=[fp], diagnose=False),
}


@pytest.mark.parametrize("operation", list(_REFUSALS))
def test_a_lean_workbench_refuses(lean, operation, tmp_path, monkeypatch):
    """Each operation that would re-read freed curve raises CurveReleased, writing nothing."""
    monkeypatch.chdir(tmp_path)
    with pytest.raises(CurveReleased, match="lean"):
        _REFUSALS[operation](*lean)
    assert not list(tmp_path.iterdir())


def test_same_stability_crossings_refuses(lean):
    """The same-stability reading reads every curve."""
    from tanglepack.numerics.diagnostics import same_stability_crossings

    with pytest.raises(CurveReleased):
        same_stability_crossings(lean[0].workbench)


def test_a_released_bridge_refuses_its_points(lean):
    """A released bridge keeps its id and ends, and refuses every row read."""
    session = lean[0]
    released = [b for b in session.workbench.bridges if b.released]
    assert released
    bridge = released[0]
    assert bridge.root is not None and bridge.tail is not None
    assert bridge.test_point is not None and bridge.point_count > 0
    for read in (bridge.get_point_array, bridge.get_cdist_array, lambda: bridge.ids):
        with pytest.raises(CurveReleased, match="without lean"):
            read()


def test_the_topology_and_the_kept_bridges_still_draw(lean):
    """topology() recomputes on the lean workbench, and its minimal trellis draws."""
    session, zone, fp = lean
    topology = session.topology()
    assert topology.dynamics.classes
    fig, ax = plt.subplots()
    try:
        plot_bridges_by_class(
            topology.trellis, topology.dynamics, ax, bridge_ids=topology.minimal.kept_bridge_ids
        )
    finally:
        plt.close(fig)
    session.blast_zone(zone, 1, fixed_point=[fp], min_separation="width", lean=True)


def test_a_dropped_bridge_a_lean_blast_released_is_not_re_admitted():
    """min_separation=None re-admits a width drop, but not one a lean blast
    released: the blast refuses before it maps anything."""
    case = build_k28(blasts=0, through="zones")
    session, zone = case.session, case.zones[0]
    workbench = session.workbench
    for _ in range(8):
        session.blast_zone(zone, 1, fixed_point=[zone.fixed_point], min_separation="width", lean=True)
        if any(b.dropped_at is not None and b.released for b in workbench.bridges):
            break
    else:
        pytest.skip("no width drop was released")
    generation = workbench.generation
    with pytest.raises(CurveReleased, match="cannot re-admit"):
        session.blast_zone(zone, 1, fixed_point=[zone.fixed_point], min_separation=None, lean=True)
    assert workbench.generation == generation


def test_a_joined_arcs_pieces_survive_a_release():
    """A joined arc over blast curve reads as before once its pieces are kept
    (``Blast._pieces``), however much else is freed."""
    case = build_k28(blasts=0, through="zones")
    session, zone = case.session, case.zones[0]
    session.blast_zone(zone, 2, fixed_point=[zone.fixed_point], diagnose=False)
    workbench = session.workbench
    branches = {m._span(stop_at_final=False)[1].index for m in workbench.manifolds.values()}
    off_branch = {
        b.id[0]: b for b in workbench.bridges
        if b.id is not None and b._span()[1].index not in branches
    }
    a, b = next(
        (a, off_branch[a.id[1]]) for a in off_branch.values()
        if a.id[1] in off_branch and off_branch[a.id[1]].manifold_key == a.manifold_key
    )
    joined = Bridge.joined([a, b])
    before = joined.get_point_array()

    workbench.release_curves(_pieces(workbench, joined))

    assert not (a.released or b.released)
    np.testing.assert_array_equal(joined.get_point_array(), before)

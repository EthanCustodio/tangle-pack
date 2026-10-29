"""A lean workbench refuses whatever would re-read the curve it freed.

A lean blast frees the interiors of the bridges nothing reads again
(``TangleWorkbench.release_curves``). From then on a read of a freed bridge,
and every operation that would re-read freed curve, raises ``CurveReleased``;
the topology and the drawing of the kept bridges still work.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import pytest

from cases import build_k28
from tanglepack.numerics.Bridge import CurveReleased
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

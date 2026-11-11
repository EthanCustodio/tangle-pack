"""Plotting tier, smoke: the trellis figure package and the loom cartoon on zones.

The package (``scripts/trellis_figures.py``) writes four panels and the
combined figure, each as PNG, SVG and pickle; the loom convenience lays a
nested topology's cartoon out on its zones (enclosures included). No style,
position or count of a particular build is pinned.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from cases import build_k10, build_nested  # noqa: E402
from tanglepack.loom import cartoon_zones, plot_trellis_cartoon  # noqa: E402
from tanglepack.topology import TrellisCartoon  # noqa: E402
from tanglepack.topology.plotting import cartoon_enclosures  # noqa: E402


def test_write_trellis_figures_writes_every_panel(tmp_path: Path) -> None:
    """Five figures, three files each."""
    from trellis_figures import PANELS, write_trellis_figures

    case = build_k10()
    folder = write_trellis_figures(case.session, tmp_path / "k10", "k10",
                                   cartoon_kw={"iterates": 1})
    expected = {f"{name}{suffix}" for name in (*PANELS, "combined")
                for suffix in (".png", ".svg", ".fig.pickle")}
    assert {path.name for path in folder.iterdir()} == expected
    assert plt.get_fignums() == []


def test_loom_cartoon_lays_nested_circles_out_on_the_zones() -> None:
    """The nested tangle's circle sits inside its parent's; every drawn end is a ranked point."""
    case = build_nested()
    session = case.session
    outer, inner = case.fixed_points
    fig, ax = plt.subplots()
    try:
        layout = plot_trellis_cartoon(
            session.topology(), session.resonance_zones.values(), session.workbench,
            ax=ax, show="all", holes=True,
        )
    finally:
        plt.close(fig)
    assert layout.zone.parents == {id(inner): id(outer)}
    for bridge_id in layout.bridges:
        for crossing in bridge_id:
            assert crossing in layout.points or any(
                key[0] == crossing for key in layout.points if isinstance(key, tuple)
            )



def test_loom_cartoon_routes_the_enclosing_lobes_the_long_way() -> None:
    """The wrapper's enclosures reach the drawing: every enclosing lobe is drawn as
    with ``cartoon_enclosures`` passed by hand, and some differ from the short way."""
    session = build_nested().session
    topology, zones = session.topology(), list(session.resonance_zones.values())
    trellis = topology.trellis
    interior, parents, order = cartoon_zones(zones, session.workbench)
    enclosures = cartoon_enclosures(
        trellis, [b.id for b in trellis.bridges if b.id is not None], parents, session.fixed_points
    )
    assert enclosures
    fig, ax = plt.subplots()
    try:
        wrapped = plot_trellis_cartoon(topology, zones, session.workbench, ax=ax, show="all")
        long_way, short_way = (
            TrellisCartoon.plot_trellis_cartoon(
                trellis, topology.classes, topology.dynamics, minimal=topology.minimal, ax=ax,
                show="all", interior_side=interior, parents=parents, arc_order=order,
                enclosures=given,
            )
            for given in (enclosures, None)
        )
    finally:
        plt.close(fig)
    for bridge_id in enclosures:
        np.testing.assert_array_equal(wrapped.paths[bridge_id], long_way.paths[bridge_id])
    assert any(
        not np.allclose(long_way.paths[bridge_id], short_way.paths[bridge_id])
        for bridge_id in enclosures
    )

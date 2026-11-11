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

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

from cases import build_k10, build_nested  # noqa: E402
from tanglepack.loom import plot_trellis_cartoon  # noqa: E402


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

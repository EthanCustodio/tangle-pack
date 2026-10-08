"""The exact iterate rule (WP2): bracket match, one batched map call, map gate.

The image of a crossing is the unique crossing on the advanced branches whose
unstable and stable brackets meet the source's carried one map step; the real
map only confirms it, in ONE batched call per inference pass.
"""

from __future__ import annotations

import logging

import numpy as np

from helpers.logs import assert_logged
from tanglepack.numerics.IterateTable import IterateTable


def _detected_links(registry) -> set[tuple[int, int]]:
    """The ``+1`` links of the detected (non-anchor) crossings."""
    return {
        (source, target)
        for source, n, target in registry.iterate_table.items()
        if n == 1 and not registry[source].is_synthetic
    }


def _gate(workbench, source, image) -> float:
    """The map gate of a link: the image of x's unstable segment plus y's own segment."""
    p0, p1 = (node.get_point() for node in source.unstable_segment)
    q0, q1 = (node.get_point() for node in image.unstable_segment)
    f0, f1 = workbench.dynamical_system.map_batch(np.vstack([p0, p1]))
    return float(np.linalg.norm(f1 - f0) + np.linalg.norm(np.ravel(q1) - np.ravel(q0)))


def test_one_map_call_per_inference_pass(henon_tangle_with_bridges, monkeypatch):
    """From an empty table, one pass rebuilds every link with ONE map call."""
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry
    links = set(registry.iterate_table.items())
    registry.iterate_table = IterateTable()

    calls = []
    system = workbench.dynamical_system
    original = system.map
    monkeypatch.setattr(system, "map", lambda point: calls.append(1) or original(point))

    assert workbench.infer_iterates() > 0
    assert len(calls) == 1
    assert set(registry.iterate_table.items()) == links


def test_map_gate_refuses_a_displaced_image(henon_tangle_with_bridges, caplog):
    """An image whose coordinates the map contradicts is not linked; nothing else moves."""
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry
    links = _detected_links(registry)
    # the end of a chain, so moving it touches no link of its own
    source, target = min(
        (source, target) for source, target in links if registry.iterate_table[target, 1] is None
    )
    image = registry[target]
    radius = _gate(workbench, registry[source], image)
    image.coords = tuple(np.ravel(image.coords) + np.array([10.0 * radius, 0.0]))

    registry.iterate_table = IterateTable()
    with caplog.at_level(logging.WARNING, logger="tanglepack.numerics.IterateInference"):
        workbench.infer_iterates()

    assert registry.iterate_table[source, 1] is None
    assert _detected_links(registry) == links - {(source, target)}
    assert_logged(caplog, logging.WARNING, "tanglepack.numerics.IterateInference", count=1)

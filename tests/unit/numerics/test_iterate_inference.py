"""The exact iterate rule (WP2): bracket match, one batched map call, map gate.

The image of a crossing is the unique crossing on the advanced branches whose
unstable and stable brackets meet the source's carried one map step; the real
map only confirms it, in ONE batched call per inference pass.
"""

from __future__ import annotations

import logging

import numpy as np
import pytest

from helpers.logs import assert_logged
from tanglepack.numerics import IterateInference
from tanglepack.numerics.IterateInference import infer_iterates
from tanglepack.numerics.IterateTable import IterateTable
from tanglepack.numerics.precise import dd_div, two_prod


def _detected_links(registry) -> set[tuple[int, int]]:
    """The ``+1`` links of the detected (non-anchor) crossings."""
    return {
        (source, target)
        for source, n, target in registry.iterate_table.items()
        if n == 1 and not registry[source].is_synthetic
    }


def _gate(workbench, source, image) -> float:
    """The map gate of a link: the image of x's unstable segment plus y's own segment."""
    p0, p1 = source.segment_ends("unstable")[0]
    q0, q1 = image.segment_ends("unstable")[0]
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


def test_an_image_claimed_across_passes_links_neither_claimant(henon_tangle_with_bridges, caplog):
    """B8: a pass whose match ``x -> y`` meets a ``y`` already linked to ``z``
    unlinks ``z`` and links neither: both are recorded as contesting ``y``."""
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry
    x, y = min(_detected_links(registry))
    z = min(i for i, ix in registry if not ix.is_synthetic and i not in (x, y))
    registry.iterate_table = IterateTable()
    registry.register_iterate(z, 1, y)
    assert registry.iterate_table[z, 1] == y

    with caplog.at_level(logging.WARNING, logger="tanglepack.numerics.IterateInference"):
        infer_iterates(workbench, [x])

    table = registry.iterate_table
    assert table[x, 1] is None and table[z, 1] is None and table[y, -1] is None
    assert registry.unresolved_links == {x: (y,), z: (y,)}
    (record,) = assert_logged(caplog, logging.WARNING, "tanglepack.numerics.IterateInference", count=1)
    assert record.args[0] == y


def test_an_image_claimed_three_times_in_one_pass_links_none(henon_tangle_with_bridges, monkeypatch):
    """B8: the claimants of one image are grouped before any is linked, so no
    order of the sources lets a third one link after the first two collide."""
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry
    image, *claimants = sorted(i for i, ix in registry if not ix.is_synthetic)[:4]
    monkeypatch.setattr(
        IterateInference, "_bracketed_image", lambda wb, reg, i, *rest: image if i in claimants else None
    )
    monkeypatch.setattr(IterateInference, "_map_gated", lambda wb, reg, matches: matches)
    registry.iterate_table = IterateTable()

    infer_iterates(workbench, claimants)

    assert registry.iterate_table[image, -1] is None
    assert all(registry.iterate_table[i, 1] is None for i in claimants)
    assert registry.unresolved_links == dict.fromkeys(claimants, (image,))


def test_a_multi_match_is_recorded_and_a_second_pass_changes_nothing(
    henon_tangle_with_bridges, monkeypatch
):
    """B8: a source whose brackets meet several images is recorded with them, in
    stable order, unlinked; re-deciding the same entries leaves the generation alone."""
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry
    # Brackets wide enough that every crossing on the advanced branches meets them.
    monkeypatch.setattr(IterateInference, "_brackets", lambda crossing: (np.array([0.0, 1e100, 0.0, 1e100]), np.zeros(4)))
    registry.iterate_table = IterateTable()
    detected = [i for i, ix in registry if not ix.is_synthetic]
    start = registry.generation

    infer_iterates(workbench, detected)  # no anchor: recording is the pass's only mutation
    generation = registry.generation
    infer_iterates(workbench, detected)

    assert generation > start
    assert set(registry.unresolved_links) == set(detected)
    for source, candidates in registry.unresolved_links.items():
        assert len(candidates) > 1 and set(candidates) <= set(detected)
        stable_cdists = [registry[j].stable_cdist for j in candidates]
        assert stable_cdists == sorted(stable_cdists)
        assert registry.iterate_table[source, 1] is None
    assert registry.generation == generation


def test_an_image_whose_stable_bracket_misses_the_carried_one_is_no_image(
    henon_tangle_with_bridges, monkeypatch
):
    """Both brackets must meet: a real link ``x -> y`` is not found once x's stable
    bracket, carried one step, sits just beside y's, though the stable window still
    holds y and the unstable brackets still meet."""
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry
    x, y = min(_detected_links(registry))
    source, image = registry[x], registry[y]
    stretch_s = workbench.manifolds[source.manifold_b_key].stretch_param
    brackets = IterateInference._brackets
    hi, lo = brackets(source)
    s_lo, s_hi = hi[2:]
    y_lo, y_hi = brackets(image)[0][2:]
    # Just beside y's bracket, on the side away from y's own cdist. The window
    # (the carried bracket widened by the branch's widest, about twice y's here)
    # still holds y: dropping the stable comparisons links x again.
    gap = 1e-3 * (y_hi - y_lo)
    if image.stable_cdist > 0.5 * (y_lo + y_hi):
        carried = (y_lo - gap - (s_hi - s_lo) / stretch_s, y_lo - gap)
    else:
        carried = (y_hi + gap, y_hi + gap + (s_hi - s_lo) / stretch_s)
    beside = (np.concatenate((hi[:2], [carried[0] * stretch_s, carried[1] * stretch_s])),
              np.concatenate((lo[:2], [0.0, 0.0])))
    monkeypatch.setattr(
        IterateInference, "_brackets", lambda crossing: beside if crossing is source else brackets(crossing)
    )
    registry.iterate_table = IterateTable()

    infer_iterates(workbench, [x])

    assert registry.iterate_table[x, 1] is None and registry.iterate_table[y, -1] is None
    assert x not in registry.unresolved_links


@pytest.mark.parametrize(("ulps", "linked"), [(-0.25, True), (0.25, False)], ids=["meets", "misses"])
def test_the_bracket_meet_reads_the_residuals(henon_tangle_with_bridges, monkeypatch, ulps, linked):
    """x's stable bracket, carried one step, starts at y's padded upper end as a
    plain float: a quarter ULP below it (in the residual) the two meet and x -> y
    is linked, a quarter ULP above it they do not."""
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry
    x, y = min(_detected_links(registry))
    source = registry[x]
    stretch_s = workbench.manifolds[source.manifold_b_key].stretch_param
    image_branch = source.manifold_b_key[0].advance_key(source.manifold_b_key, 1)
    _cdists, ids, near_hi, near_lo = registry.stable_brackets(image_branch)
    row = ids.tolist().index(y)
    padded_hi, padded_lo = IterateInference._padded(near_hi[row], near_lo[row])
    end_hi, end_lo = float(padded_hi[3]), float(padded_lo[3]) + ulps * np.spacing(padded_hi[3])
    # A source stable end whose carry is (end_hi, end_lo): dd_div's hi is fl(h / s).
    guess = end_hi * stretch_s
    h = next(c for c in (guess, np.nextafter(guess, 0.0), np.nextafter(guess, np.inf)) if c / stretch_s == end_hi)
    p, e = two_prod(end_hi, stretch_s)
    lo = end_lo * stretch_s - ((h - p) - e)
    assert dd_div(h, lo, stretch_s)[0] == end_hi
    hi, lo_row = IterateInference._brackets(source)
    width = (hi[3] - hi[2]) / stretch_s
    beside = (np.concatenate((hi[:2], [h, (end_hi + width) * stretch_s])), np.concatenate((lo_row[:2], [lo, 0.0])))
    brackets = IterateInference._brackets
    monkeypatch.setattr(
        IterateInference, "_brackets", lambda crossing: beside if crossing is source else brackets(crossing)
    )
    registry.iterate_table = IterateTable()

    infer_iterates(workbench, [x])

    assert (registry.iterate_table[x, 1] == y) is linked
    assert x not in registry.unresolved_links


def test_remove_iterate_forgets_the_link_and_its_twin(henon_tangle_with_bridges):
    """B8: removing ``x -> y`` leaves neither ``(x, 1)`` nor ``(y, -1)``."""
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry
    x, y = min(_detected_links(registry))
    assert registry.iterate_table[x, 1] == y
    generation = registry.generation

    registry.remove_iterate(x, 1)

    assert registry.iterate_table[x, 1] is None and registry.iterate_table[y, -1] is None
    assert registry.generation > generation


def test_a_later_pass_does_not_relink_a_contested_image(henon_tangle_with_bridges):
    """B8: a pass over one claimant re-decides only its own entry; the other
    claimant, left contesting ``y`` by an earlier pass, still contests it."""
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry
    x, y = min(_detected_links(registry))
    z = min(i for i, ix in registry if not ix.is_synthetic and i not in (x, y))
    registry.iterate_table = IterateTable()
    registry.register_iterate(z, 1, y)
    infer_iterates(workbench, [x])
    generation = registry.generation

    infer_iterates(workbench, [x])

    assert registry.iterate_table[x, 1] is None and registry.iterate_table[y, -1] is None
    assert registry.unresolved_links == {x: (y,), z: (y,)}
    assert registry.generation == generation


def test_a_pass_re_decides_its_own_sources(henon_tangle_with_bridges, monkeypatch):
    """B8: once the brackets single out one image again, the next pass drops the
    source's entry and links it: no entry outlives the pass that re-decides it."""
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench.intersection_registry
    links = set(registry.iterate_table.items())
    with monkeypatch.context() as patch:
        patch.setattr(IterateInference, "_brackets", lambda crossing: (np.array([0.0, 1e100, 0.0, 1e100]), np.zeros(4)))
        registry.iterate_table = IterateTable()
        infer_iterates(workbench)
    assert registry.unresolved_links

    infer_iterates(workbench)

    assert registry.unresolved_links == {}
    assert set(registry.iterate_table.items()) == links

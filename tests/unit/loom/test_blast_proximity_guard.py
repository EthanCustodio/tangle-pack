"""The blast proximity guard drops near-coincident sibling bridges.

Near the fixed point, successive bridge images pile up (homoclinic accumulation)
until distinct unstable curves run within machine precision and merge into an
artifact. ``min_separation`` stops the blast on the well-resolved side of that
limit by dropping a child whose interior comes too close to curve already in
the tangle (the bridges present when the blast started, every child an
earlier step produced and the siblings kept so far). A child ending on the chosen pip's orbit is exempt (B3, N-a): the
guard never decides the depth of ``T(q0)``.
"""

from __future__ import annotations

import importlib

import numpy as np
import pytest
from scipy.spatial import cKDTree

from tanglepack.loom.Blast import RESOLUTION_SEPARATION
from tanglepack.numerics import DEFAULT_AREA_CUTOFF

import cases

MIN_SEPARATION = 0.02


def _interior(bridge) -> np.ndarray:
    """A bridge's polyline with the outer 10 % trimmed off each end.

    Adjacent bridges legitimately meet at their shared crossings, so the guard
    compares interiors; this is the documented contract of ``min_separation``
    (the bodies of two bridges, not their touching ends).

    Args:
        bridge: The bridge.

    Returns:
        The ``(N, 2)`` interior points (all points for a short bridge).
    """
    pts = np.asarray(bridge.get_point_array())
    if len(pts) < 5:
        return pts
    margin = max(1, len(pts) // 10)
    return pts[margin:-margin]


def _distance(points: np.ndarray, cloud: list[np.ndarray]) -> float:
    """Smallest distance from ``points`` to the union of ``cloud``.

    Args:
        points: The query points.
        cloud: Reference point arrays.

    Returns:
        The distance, ``inf`` when either side is empty.
    """
    cloud = [c for c in cloud if len(c)]
    if not len(points) or not cloud:
        return float("inf")
    return float(cKDTree(np.vstack(cloud)).query(points)[0].min())


def _pip_orbit(session, fp) -> set[int]:
    """The chosen pip of ``fp`` and its registered forward iterates (empty without one)."""
    chosen = session.chosen_strong_pip(fp)
    if chosen is None:
        return set()
    trellis = session.trellis(fp)
    orbit = [trellis.registry.find(chosen)]
    while (image := trellis.iterate(orbit[-1], 1)) is not None:
        orbit.append(image)
    return set(orbit)


@pytest.mark.parametrize("pip_chosen", [True, False], ids=["pip", "no_pip"])
def test_min_separation_drops_close_bridges(henon_p3_session, pip_chosen):
    """Kept children keep their distance unless they end on the pip's orbit;
    dropped ones really were too close.

    Guards the 2026 "zig-zag" fix: every kept child's interior that does not
    end on the chosen pip's orbit stays at least ``min_separation`` from the
    curve accumulated before it, and every child the guard dropped came within
    ``min_separation`` of it. With the pip forgotten nothing is exempt, so the
    guard demonstrably fires (the fixture's one too-close child ends on the
    orbit).
    """
    session, _fp3, fp1, zone = henon_p3_session
    if not pip_chosen:
        session.forget_strong_pip(fp1)
    reference = [_interior(b) for b in session.workbench.bridges]

    guarded = session.blast_zone(
        zone, num_iterations=3, fixed_point=fp1, min_separation=MIN_SEPARATION
    )
    orbit = _pip_orbit(session, fp1)

    if not pip_chosen:
        assert guarded.too_close > 0
    assert sum(len(s.discarded_too_close) for s in guarded.steps) == guarded.too_close
    assert guarded.completed_iterations == 3

    checked = 0
    for step in guarded.steps:
        # kept order is the acceptance order: each is checked against the
        # accumulated curve plus the siblings accepted before it
        accepted: list[np.ndarray] = []
        for child in step.kept_interior:
            interior = _interior(child)
            if not orbit.intersection(child.id or ()):
                assert _distance(interior, reference + accepted) >= MIN_SEPARATION
            accepted.append(interior)
            checked += 1
        for child in step.discarded_too_close:
            assert _distance(_interior(child), reference + accepted) < MIN_SEPARATION
            checked += 1
        reference += [_interior(c) for c in step.produced_children]
    assert checked > 0


def test_min_separation_none_is_unchanged_default(henon_p3_session):
    """Default (None) keeps the original behavior: nothing dropped for proximity."""
    session, _fp3, fp1, zone = henon_p3_session
    result = session.blast_zone(zone, num_iterations=1, fixed_point=fp1)
    assert result.too_close == 0


def test_resolution_separation_follows_the_area_cutoff(henon_p3_session):
    """``"resolution"`` drops at ``RESOLUTION_SEPARATION`` scaled as ``area_cutoff ** (2/3)``.

    The fixture's outer saddle is grown at ``1e-4``; with its pip forgotten a
    six-iteration blast drops a child, which records the derived separation.
    """
    session, _fp3, fp1, zone = henon_p3_session
    session.forget_strong_pip(fp1)
    expected = RESOLUTION_SEPARATION * (session.area_cutoff / DEFAULT_AREA_CUTOFF) ** (2 / 3)
    result = session.blast_zone(
        zone, num_iterations=6, fixed_point=fp1, min_separation="resolution"
    )
    dropped = [c for s in result.steps for c in s.discarded_too_close]
    assert dropped
    assert all(c.dropped_at == pytest.approx(expected) for c in dropped)


def _blast_view(result) -> list:
    """Every step's children by verdict, each child as its id and its points."""

    def children(bridges):
        return [(b.id, np.asarray(b.get_point_array()).tobytes()) for b in bridges]

    return [
        (children(step.kept_interior), children(step.discarded_too_close),
         children(step.discarded_exterior), children(step.already_known))
        for step in result.steps
    ] + [result.points_created, result.too_close, result.stop_reason]


def test_workers_give_the_same_blast(monkeypatch):
    """A ``cases.SERIES`` blast with ``min_separation`` keeps and drops the same
    children, creates the same points and registers the same crossings at
    workers 1 and 4 (every map batch split, the separation queried on 4 cores)."""
    monkeypatch.setattr(importlib.import_module("tanglepack.numerics.DynamicalSystem"),
                        "MIN_PARALLEL_POINTS", 16)
    series = cases.SERIES["k28_series"]
    runs = []
    for workers in (1, 4):
        case = series.build()
        case.session.workers = workers
        zone = max(case.zones, key=lambda zone: zone.area)
        results = [
            case.session.blast_zone(zone, num_iterations=1, fixed_point=[zone.fixed_point],
                                    min_separation=series.min_separation, diagnose=True)
            for _ in range(series.blasts)
        ]
        registry = case.workbench.intersection_registry
        runs.append((
            [_blast_view(result) for result in results],
            sorted((iid, ix.unstable_cdist, ix.stable_cdist) for iid, ix in registry),
            sorted(registry.iterate_table.items()),
        ))
    serial, split = runs
    assert any(views[-2] for views in serial[0]), "some child must be dropped"
    assert split == serial

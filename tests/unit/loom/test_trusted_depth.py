"""``BlastResult.trusted_depth``: how deep the registered pip orbit contracts like the real one.

Inside ``TRUST_LINEAR_RADIUS`` of its orbit point each period of the pip orbit
must divide the distance by the stable multiplier to within ``TRUST_BAND``; the
first period that does not ends the trusted orbit (2026-10-29, the test of
``scripts/henon_k25_errors_figures.py``). It is a flag: at most ``pip_depth``,
and nothing reads it. The kernel is checked on hand-placed orbits (synthetic
crossings at chosen distances), the blast on a real build by rule only.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import numpy as np
import pytest

from helpers.fakes import bare_fixed_point
from tanglepack.loom import Blast
from tanglepack.loom.Blast import TRUST_BAND, TRUST_LINEAR_RADIUS, _trusted_depth
from tanglepack.numerics.IntersectionRegistry import IntersectionRegistry

BETA = 0.25


def _orbit(distances, period=1):
    """A session stand-in whose registry holds crossings at ``distances`` from the orbit point, and their ids."""
    fp = bare_fixed_point(period, beta=BETA, coordinates=[(0.0, 0.0)] * period)
    registry = IntersectionRegistry()
    ids = [
        registry.add_synthetic(
            (d, 0.0), unstable_cdist=1.0 + j, stable_cdist=d,
            manifold_a_key=(fp, "unstable", j % period, 0),
            manifold_b_key=(fp, "stable", j % period, 0),
        )
        for j, d in enumerate(distances)
    ]
    return SimpleNamespace(workbench=SimpleNamespace(intersection_registry=registry)), fp, ids


def test_an_orbit_that_contracts_exactly_is_trusted_to_its_end():
    session, fp, ids = _orbit([0.5 * TRUST_LINEAR_RADIUS * BETA**j for j in range(8)])
    assert _trusted_depth(session, fp, ids) == len(ids) - 1


@pytest.mark.parametrize("off_band", [TRUST_BAND[0] * 0.9, TRUST_BAND[1] * 1.1])
def test_the_first_off_band_step_ends_it(off_band):
    distances = [0.5 * TRUST_LINEAR_RADIUS * BETA**j for j in range(8)]
    distances[5:] = [d / off_band for d in distances[5:]]  # only step 4 -> 5 is off band
    session, fp, ids = _orbit(distances)
    assert _trusted_depth(session, fp, ids) == 4


def test_outside_the_linear_radius_nothing_is_tested():
    # Steps that start beyond the radius are free; the first one inside is tested.
    distances = [4 * TRUST_LINEAR_RADIUS, 3 * TRUST_LINEAR_RADIUS, 0.5 * TRUST_LINEAR_RADIUS]
    distances += [distances[-1] * BETA**j for j in range(1, 4)]
    session, fp, ids = _orbit(distances)
    assert _trusted_depth(session, fp, ids) == len(ids) - 1


def test_a_period_divides_the_distance_by_the_multiplier():
    # Period 2: the orbit alternates orbit points; only every second point is compared.
    period = 2
    distances = [0.5 * TRUST_LINEAR_RADIUS * BETA ** (j - j % period) for j in range(8)]
    session, fp, ids = _orbit(distances, period=period)
    assert _trusted_depth(session, fp, ids) == len(ids) - 1


def _composed_orbit(session, fp) -> list[int]:
    """``[q0, ..., f^m(q0)]`` of ``fp``'s chosen pip, composed through ``trellis.iterate``."""
    trellis = session.trellis(fp)
    orbit = [trellis.registry.find(session.chosen_strong_pip(fp))]
    while (image := trellis.iterate(orbit[-1], 1)) is not None:
        orbit.append(image)
    return orbit


def test_a_blast_trusts_its_own_pip_orbit(henon_p3_session):
    session, _fp3, fp1, zone = henon_p3_session
    result = session.blast_zone(zone, num_iterations=2, fixed_point=fp1, min_separation=0.02)
    orbit = _composed_orbit(session, fp1)
    assert result.trusted_depth.keys() == result.pip_depth.keys() == {fp1}
    assert result.pip_depth[fp1] == len(orbit) - 1
    assert result.trusted_depth[fp1] == _trusted_depth(session, fp1, orbit)


def test_a_blast_logs_a_trust_shallower_than_its_depth(henon_p3_session, monkeypatch, caplog):
    # No step is in an empty band, so the trust ends at the first period.
    monkeypatch.setattr(Blast, "TRUST_BAND", (1.0, 1.0))
    monkeypatch.setattr(Blast, "TRUST_LINEAR_RADIUS", np.inf)
    session, _fp3, fp1, zone = henon_p3_session
    with caplog.at_level(logging.INFO, logger=Blast.__name__):
        result = session.blast_zone(zone, num_iterations=2, fixed_point=fp1, min_separation=0.02)

    assert result.trusted_depth[fp1] == fp1.period - 1 < result.pip_depth[fp1]
    logged = [r.args for r in caplog.records if r.name == Blast.__name__ and r.levelno == logging.INFO]
    assert (fp1.label, result.pip_depth[fp1], result.trusted_depth[fp1]) in logged

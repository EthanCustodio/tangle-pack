"""``BlastResult.trusted_depth``: how deep the registered pip orbit contracts like the real one.

Inside ``TRUST_LINEAR_RADIUS`` of its orbit point each period of the pip orbit
must divide the distance by the stable multiplier to within ``TRUST_BAND``; the
first period that does not ends the trusted orbit (2026-10-29, the test of
``scripts/henon_k25_errors_figures.py``). It is a flag: at most ``pip_depth``,
and nothing reads it. The kernel is checked on hand-placed orbits (synthetic
crossings at chosen distances), the blast on a real build by rule only.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from helpers.fakes import bare_fixed_point
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


def test_a_blast_trusts_no_deeper_than_it_registers(henon_p3_session):
    session, _fp3, fp1, zone = henon_p3_session
    result = session.blast_zone(zone, num_iterations=2, fixed_point=fp1, min_separation=0.02)
    assert result.trusted_depth.keys() == result.pip_depth.keys()
    assert all(0 <= result.trusted_depth[fp] <= result.pip_depth[fp] for fp in result.pip_depth)

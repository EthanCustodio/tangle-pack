"""The hydrogen example: an area-preserving, invertible, batch-capable Poincare map.

``examples/hydrogen/hydrogen_map.py`` is the heteroclinic example (two period-1
saddles swapped by ``-I``). Its map must be what the library assumes of any map:
area-preserving, with an inverse that undoes it (the reversor claim ``M^-1 = R M
R``), and batch-capable over the points of a ``(2, N)`` array. The symmetry
``-I`` must commute with it, the identification ``z ~ -z`` must preserve area and
conjugate the map to ``map_identified``, the finite-difference Jacobian must be
consistent, and ``find_fixed_point`` must return a saddle.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack.examples.hydrogen.hydrogen_map import HydrogenSystem

#: Points in the allowed region at ``B = 3.5, E = 1``: inside the diamond, beyond
#: each saddle and near one.
POINTS = np.array([[0.3, 0.2], [-0.8, 0.4], [1.0, -0.3], [0.98, 0.02]])


@pytest.fixture(scope="module")
def system() -> HydrogenSystem:
    """The parameters of the heteroclinic zone test."""
    return HydrogenSystem(3.5, 1.0)


def test_inverse_undoes_the_map(system):
    """``map_inverse(map(x))`` and ``map(map_inverse(x))`` return ``x`` to 1e-12."""
    for point in POINTS:
        assert np.abs(system.map_inverse(system.map(point)) - point).max() < 1e-12
        assert np.abs(system.map(system.map_inverse(point)) - point).max() < 1e-12


def test_batch_equals_per_point(system):
    """A ``(2, N)`` batch maps to the per-point images exactly, both directions."""
    for fn in (system.map, system.map_inverse):
        assert np.array_equal(fn(POINTS.T).T, np.vstack([fn(point) for point in POINTS]))


def test_minus_identity_commutes_with_the_map(system):
    """``M(-x) = -M(x)``: the potential is even in ``v``."""
    assert np.abs(system.map(-POINTS.T) + system.map(POINTS.T)).max() < 1e-12


def test_the_map_preserves_area_and_its_jacobian_is_consistent(system):
    """``det J = 1`` to 1e-7, and the default step agrees with a tenfold one to 1e-4."""
    for point in POINTS:
        jacobian = system.jacobian(point)
        assert abs(np.linalg.det(jacobian) - 1.0) < 1e-7
        coarse = HydrogenSystem.jacobian_of(system.map, eps=1e-4)(point)
        assert np.abs(jacobian - coarse).max() < 1e-4


def test_identification_preserves_area_and_conjugates_the_map(system):
    """``identify`` has ``det = 1``, ``unidentify`` is a right inverse, and
    ``map_identified(identify(x)) = identify(map(x))`` (either preimage works)."""
    images = HydrogenSystem.identify(POINTS.T)
    assert np.abs(HydrogenSystem.identify(HydrogenSystem.unidentify(images)) - images).max() < 1e-12
    for point in POINTS:
        jacobian = HydrogenSystem.jacobian_of(HydrogenSystem.identify, eps=1e-6)(point)
        assert abs(np.linalg.det(jacobian) - 1.0) < 1e-8
    expected = HydrogenSystem.identify(system.map(POINTS.T))
    assert np.abs(system.map_identified(images) - expected).max() < 1e-10


def test_the_fixed_point_is_a_saddle(system):
    """``(v*, 0)`` is fixed to 1e-12 and its eigenvalues are ``lambda > 1 > 1/lambda > 0``."""
    fixed_point = system.find_fixed_point()
    assert fixed_point[1] == 0.0
    assert np.abs(system.map(fixed_point) - fixed_point).max() < 1e-12
    small, large = sorted(np.linalg.eigvals(system.jacobian(fixed_point)).real)
    assert 0.0 < small < 1.0 < large
    assert abs(small * large - 1.0) < 1e-7

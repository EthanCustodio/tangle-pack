import pytest
import numpy as np

from tanglepack import ManifoldInitializer, Point
from tanglepack import DynamicalSystem
from tanglepack import FixedPointSolver
from tanglepack.examples import (
    HENON_K10,
    HENON_P3,
    henon_jacobian,
    henon_map as _henon_map_factory,
    henon_map_inverse as _henon_map_inverse_factory,
    saddle_guesses,
)

# The binary-horseshoe Hénon map and its inverse.
henon_map = _henon_map_factory(*HENON_K10)
henon_map_inverse = _henon_map_inverse_factory(*HENON_K10)


# --------------------------------------------------------------------------- #
# Plan 1.6 -- ManifoldInitializer called a nonexistent FixedPoint.has_inversion()
# --------------------------------------------------------------------------- #
def _k10_fixed_point():
    """The k=10 binary-horseshoe saddle at [4, -4], oriented, with k_value set."""
    henon = DynamicalSystem(henon_map, henon_map_inverse)
    fp_solver = FixedPointSolver(henon)
    fixed_point = fp_solver.construct_fixed_point(saddle_guesses(*HENON_K10)["saddle"])
    man_maker = ManifoldInitializer(henon)
    man_maker.orient_manifolds(
        fixed_point,
        {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])},
    )
    return man_maker, fixed_point


def test_two_branch_initialization_grows_opposite_directions():
    """Branch 1 is the eigenvector direction negated, so the two first points
    straddle the fixed point. Previously an AttributeError (``has_inversion``)."""
    man_maker, fixed_point = _k10_fixed_point()

    segments = man_maker.construct_kevin_way(fixed_point, "unstable", num_branches=2)

    assert set(segments) == {(0, 0), (0, 1)}

    origin = np.asarray(fixed_point.coordinates[0]).reshape(-1)
    first_0 = np.asarray(
        man_maker.get_first_point(fixed_point, 0, 0, "unstable")
    ).reshape(-1)
    first_1 = np.asarray(
        man_maker.get_first_point(fixed_point, 0, 1, "unstable")
    ).reshape(-1)

    assert np.allclose(first_0 - origin, -(first_1 - origin))


# --------------------------------------------------------------------------- #
# Plan 1.7 -- construct_kevin_way walks the chain with FixedPoint.advance_key.
# The dict it returns must be unchanged on period 1 and period 3.
# --------------------------------------------------------------------------- #
_p3_map = _henon_map_factory(*HENON_P3)
_p3_map_inverse = _henon_map_inverse_factory(*HENON_P3)
_p3_jacobian = henon_jacobian(*HENON_P3)


def _p3_fixed_point():
    system = DynamicalSystem(_p3_map, _p3_map_inverse, _p3_jacobian)
    fp_solver = FixedPointSolver(system)
    fixed_point = fp_solver.construct_fixed_point(saddle_guesses(*HENON_P3)["period_3"])
    man_maker = ManifoldInitializer(system)
    man_maker.orient_manifolds(
        fixed_point,
        {"unstable": np.array([0, -1]), "stable": np.array([-1, -1])},
    )
    return man_maker, fixed_point


@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_kevin_way_period_one_keys_and_points(stability):
    man_maker, fixed_point = _k10_fixed_point()

    segments = man_maker.construct_kevin_way(fixed_point, stability)

    assert list(segments) == [(0, 0)]
    assert segments[(0, 0)].root is fixed_point.branch_points[0]
    assert len(segments[(0, 0)].get_point_array()) == 3


@pytest.mark.parametrize(
    ("stability", "chain"),
    [("unstable", [0, 1, 2]), ("stable", [0, 2, 1])],
    ids=["unstable", "stable"],
)
def test_kevin_way_period_three_orbit_chain(stability, chain):
    """Chain order: unstable walks 0, 1, 2; stable walks 0, 2, 1 (the inverse map
    grows the stable fundamental segment).

    The chain is read off the iterate links: from the innermost seeded point of
    orbit 0, the growth-direction link (``next_iterate`` for unstable,
    ``prev_iterate`` for stable) visits the other orbit points' segments in
    chain order, at strictly increasing canonical distance.
    """
    man_maker, fixed_point = _p3_fixed_point()

    segments = man_maker.construct_kevin_way(fixed_point, stability)

    assert set(segments) == {(0, 0), (1, 0), (2, 0)}
    owner = {
        id(node): key
        for key, manifold in segments.items()
        for node in manifold.get_point_array(return_nodes=True)
    }
    growth_link = "next_iterate" if stability == "unstable" else "prev_iterate"
    node = min(
        (
            n
            for n in segments[(0, 0)].get_point_array(return_nodes=True)
            if isinstance(n, Point) and n.cdist is not None
        ),
        key=lambda n: n.cdist,
    )
    visited, cdists = [], []
    while node is not None and id(node) in owner and len(visited) < len(chain):
        visited.append(owner[id(node)][0])
        cdists.append(node.cdist)
        node = getattr(node, growth_link)
    assert visited == chain
    assert cdists == sorted(cdists) and len(set(cdists)) == len(cdists)
    for (orbit_index, _branch_index), manifold in segments.items():
        assert manifold.root is fixed_point.branch_points[orbit_index]
        # every segment sits on the orbit point it is keyed to
        first = manifold.walk_fwd(None, manifold.root, 0)
        assert first is not None
        assert np.linalg.norm(
            np.asarray(first.get_point()).reshape(-1)
            - np.asarray(fixed_point.coordinates[orbit_index]).reshape(-1)
        ) < np.linalg.norm(
            np.asarray(first.get_point()).reshape(-1)
            - np.asarray(
                fixed_point.coordinates[(orbit_index + 1) % fixed_point.period]
            ).reshape(-1)
        )

import numpy as np

from tanglepack import ManifoldMachine
from tanglepack import DynamicalSystem


def henon_map(point):
    """Hénon map, binary-horseshoe parameters."""
    k, b = (10, 1)
    x, y = point[0], point[1]
    return np.array([y - k + x**2, -b * x])


def henon_map_inverse(point):
    """Inverse Hénon map."""
    k, b = (10, 1)
    x, y = point[0], point[1]
    return np.array([-y / b, x + k - (y**2) / (b**2)])


def test_machine_initialization():
    system = DynamicalSystem(henon_map, henon_map_inverse)
    machine = ManifoldMachine(system)
    assert machine.system is system


# --------------------------------------------------------------------------- #
# Plan 1.8 -- new_grow_manifold rewritten as a single advance_key chain walk.
#
# The rewrite must be behaviour-preserving for every non-inversion fixed point
# (k_value == period there, so the old nested loop degenerated to one pass).
# ``_old_new_grow_manifold`` below is the pre-rewrite implementation copied
# verbatim; the tests grow two identical workbenches, one with each, and compare
# the resulting point arrays exactly.
# --------------------------------------------------------------------------- #
import types

from tanglepack import TangleWorkbench
from tanglepack.numerics.BaseManifold import BaseManifold
from tanglepack.numerics.BranchPoint import BranchPoint


def _old_new_grow_manifold(self, fixed_point, stability, branch_index=None):
    """Verbatim copy of ManifoldMachine.new_grow_manifold before plan 1.8."""
    orbit_indices = fixed_point.get_iterable_array(stability, shift=1)

    if fixed_point.check_inversion():
        branches_to_grow = fixed_point.get_branch_array()
    else:
        branches_to_grow = (
            list(range(fixed_point.num_branches))
            if branch_index is None
            else [branch_index]
        )

    for b in branches_to_grow:
        current_manifold = BaseManifold(
            fixed_point.branch_points[orbit_indices[0]],
            stability,
            stretch_param=1,
            fixed_point=fixed_point,
            branch_index=b,
        )

        temp_root = current_manifold.root
        if isinstance(current_manifold.root, BranchPoint):
            current_manifold.root = current_manifold.walk_fwd(None, temp_root)

        if current_manifold.root is None:
            continue  # branch not initialized

        current_manifold.stretch_param = current_manifold.root.stretch_param

        for i in range(fixed_point.period):

            for bi in (
                fixed_point.get_branch_array()
                if fixed_point.check_inversion()
                else [b]
            ):

                self.iterate_manifold(current_manifold)

                next_index = (i + 1) % fixed_point.period
                next_orbit_idx = orbit_indices[next_index]

                next_manifold = BaseManifold(
                    root=fixed_point.branch_points[next_orbit_idx],
                    stability=stability,
                    stretch_param=current_manifold.stretch_param,
                    fixed_point=fixed_point,
                    branch_index=bi,
                )

                temp_root = next_manifold.root
                if isinstance(next_manifold.root, BranchPoint):
                    next_manifold.root = next_manifold.walk_fwd(None, temp_root)

                current_manifold = next_manifold


def _p3_map(point):
    k, b = 2, 1
    x, y = point
    return np.stack([y - k + x**2, -b * x], axis=0)


def _p3_map_inverse(point):
    k, b = 2, 1
    x, y = point
    return np.stack([-y / b, x + k - (y**2) / (b**2)], axis=0)


def _p3_jacobian(point):
    k, b = 2, 1
    x, y = point
    return np.array([[2 * x, 1], [-b, 0]])


def _batched_henon(point):
    k, b = 10, 1
    x, y = point
    return np.stack([y - k + x**2, -b * x], axis=0)


def _batched_henon_inverse(point):
    k, b = 10, 1
    x, y = point
    return np.stack([-y / b, x + k - (y**2) / (b**2)], axis=0)


def _build_k10(use_old: bool):
    wb = TangleWorkbench(_batched_henon, _batched_henon_inverse)
    if use_old:
        wb._man_machine.new_grow_manifold = types.MethodType(
            _old_new_grow_manifold, wb._man_machine
        )
    fp = wb.construct_fixed_point([4, -4])
    wb.orient_eigenvectors(
        fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
    )
    wb.initialize_both_manifolds(fp)
    wb.grow_n_times(fp, "unstable", num_iterations=2)
    wb.grow_n_times(fp, "stable", num_iterations=2)
    return wb, fp


def _build_p3(use_old: bool):
    wb = TangleWorkbench(_p3_map, _p3_map_inverse, _p3_jacobian)
    if use_old:
        wb._man_machine.new_grow_manifold = types.MethodType(
            _old_new_grow_manifold, wb._man_machine
        )
    fp = wb.construct_fixed_point([[0, 1], [-1, 0], [-1, 1]])
    wb.orient_eigenvectors(
        fp, {"unstable": np.array([0, -1]), "stable": np.array([-1, -1])}
    )
    wb.initialize_both_manifolds(fp)
    wb.grow_n_times(fp, "unstable", num_iterations=4)
    wb.grow_n_times(fp, "stable", num_iterations=3)
    return wb, fp


def _point_arrays(workbench):
    return {
        (key[1], key[2], key[3]): manifold.get_point_array()
        for key, manifold in workbench.manifolds.items()
    }


def _assert_same_growth(builder):
    old_wb, _old_fp = builder(use_old=True)
    new_wb, _new_fp = builder(use_old=False)

    old_points = _point_arrays(old_wb)
    new_points = _point_arrays(new_wb)

    assert set(old_points) == set(new_points)
    for key in old_points:
        assert old_points[key].shape == new_points[key].shape, key
        assert np.array_equal(old_points[key], new_points[key]), key


def test_new_grow_manifold_matches_old_period_one():
    _assert_same_growth(_build_k10)


def test_new_grow_manifold_matches_old_period_three():
    _assert_same_growth(_build_p3)

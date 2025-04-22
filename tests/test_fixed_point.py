import pytest
import numpy as np
from tanglepack.BranchPoint import BranchPoint
from tanglepack.FixedPoint import FixedPoint


def test_create_fixed_point():

    period = 3
    num_branches = 2
    p = FixedPoint(period, num_branches)

    assert len(p.coordinates) == period
    assert p.branch_points[0].num_branches == num_branches
    assert len(p.stable_eigenvalues) == period
    assert len(p.stable_eigenvectors) == period
    assert len(p.unstable_eigenvalues) == period
    assert len(p.unstable_eigenvectors) == period
    assert len(p.jacobians) == period

    assert np.shape(p.stable_eigenvectors[0]) == (2, 1)
    assert np.shape(p.unstable_eigenvectors[0]) == (2, 1)

    assert np.shape(p.jacobians[0]) == (2, 2)





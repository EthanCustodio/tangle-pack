import pytest
import numpy as np
from tanglepack import BranchPoint
from tanglepack import FixedPoint


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




# --------------------------------------------------------------------------- #
# Plan 1.7 / 2.5 -- FixedPoint.advance_key: the single (orbit, branch) advance rule
# --------------------------------------------------------------------------- #
def _bare_fixed_point(period: int, *, inversion: bool = False) -> FixedPoint:
    """A FixedPoint carrying only the data ``advance_key`` needs (period, k_value)."""
    fp = FixedPoint(period, 2)
    if inversion:
        fp.unstable_eigenvalues = [-2.0] * period
    else:
        fp.unstable_eigenvalues = [2.0] * period
    fp.set_k_value()
    return fp


def test_advance_key_period_one_is_the_identity():
    fp = _bare_fixed_point(1)
    key = (fp, "unstable", 0, 0)

    assert fp.advance_key(key, 1) == key
    assert fp.advance_key(key, -1) == key
    assert fp.advance_key(key, 7) == key


def test_advance_key_period_three_cycles_the_orbit_index():
    fp = _bare_fixed_point(3)
    assert fp.check_inversion() is False

    key = (fp, "unstable", 0, 0)
    forward = [key]
    for _ in range(3):
        forward.append(fp.advance_key(forward[-1], 1))

    assert [k[2] for k in forward] == [0, 1, 2, 0]
    # branch and stability pass through untouched
    assert all(k[1] == "unstable" and k[3] == 0 for k in forward)

    backward = [key]
    for _ in range(3):
        backward.append(fp.advance_key(backward[-1], -1))
    assert [k[2] for k in backward] == [0, 2, 1, 0]

    # stability is carried, never interpreted
    stable_key = (fp, "stable", 1, 0)
    assert fp.advance_key(stable_key, 1) == (fp, "stable", 2, 0)


def test_advance_key_inversion_flips_branch_on_wrap():
    fp = _bare_fixed_point(2, inversion=True)
    assert fp.check_inversion() is True
    assert fp.k_value == 4

    key = (fp, "unstable", 0, 0)
    forward = []
    current = key
    for _ in range(4):
        current = fp.advance_key(current, 1)
        forward.append((current[2], current[3]))

    assert forward == [(1, 0), (0, 1), (1, 1), (0, 0)]

    backward = []
    current = key
    for _ in range(4):
        current = fp.advance_key(current, -1)
        backward.append((current[2], current[3]))

    assert backward == [(1, 1), (0, 1), (1, 0), (0, 0)]


def test_advance_key_is_a_group_action_with_period_k_value():
    for fp in (
        _bare_fixed_point(1),
        _bare_fixed_point(3),
        _bare_fixed_point(2, inversion=True),
        _bare_fixed_point(3, inversion=True),
    ):
        for branch in range(2):
            key = (fp, "unstable", 0, branch)
            assert fp.advance_key(key, fp.k_value) == key
            assert fp.advance_key(key, 0) == key
            for a in range(-5, 6):
                for b in range(-5, 6):
                    assert fp.advance_key(fp.advance_key(key, a), b) == fp.advance_key(
                        key, a + b
                    )


def test_advance_key_rejects_a_foreign_fixed_point():
    fp = _bare_fixed_point(3)
    other = _bare_fixed_point(3)
    with pytest.raises(ValueError):
        fp.advance_key((other, "unstable", 0, 0), 1)


@pytest.mark.parametrize("inversion", [False, True])
def test_advance_key_rejects_an_out_of_range_branch(inversion):
    """Both paths reject a branch index the fixed point does not have."""
    fp = _bare_fixed_point(2, inversion=inversion)
    assert fp.num_branches == 2

    for bad_branch in (-1, 2):
        with pytest.raises(ValueError, match="branch_index"):
            fp.advance_key((fp, "unstable", 0, bad_branch), 1)


def test_advance_key_rejects_a_branch_beyond_a_single_branch_point():
    fp = FixedPoint(3, 1)
    fp.unstable_eigenvalues = [2.0] * 3
    fp.set_k_value()

    assert fp.advance_key((fp, "unstable", 0, 0), 1) == (fp, "unstable", 1, 0)
    with pytest.raises(ValueError, match="branch_index"):
        fp.advance_key((fp, "unstable", 0, 1), 1)

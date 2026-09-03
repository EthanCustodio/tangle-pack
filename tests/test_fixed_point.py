import pytest
import numpy as np
from tanglepack import BranchPoint
from tanglepack import FixedPoint


def test_create_fixed_point():

    period = 3
    p = FixedPoint(period)

    assert len(p.coordinates) == period
    # Both eigendirection slots are allocated regardless of inversion; how many
    # of them the k-chain visits is FixedPoint.num_branches (plan 2.8).
    assert p.branch_points[0].num_branches == 2
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
    fp = FixedPoint(period)
    sign = -1.0 if inversion else 1.0
    fp.unstable_eigenvalues = [sign * 2.0] * period
    fp.stable_eigenvalues = [sign * 0.5] * period
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

    for bad_branch in (-1, 2):
        with pytest.raises(ValueError, match="branch_index"):
            fp.advance_key((fp, "unstable", 0, bad_branch), 1)


def test_advance_key_carries_the_second_eigendirection_of_a_simple_saddle():
    """A saddle without inversion has one branch per eigendirection, each its own
    invariant chain: the branch index is carried, never flipped, and both are
    valid keys even though ``num_branches`` (the chain length in branches) is 1."""
    fp = _bare_fixed_point(3)
    assert fp.num_branches == 1

    assert fp.advance_key((fp, "unstable", 0, 0), 1) == (fp, "unstable", 1, 0)
    assert fp.advance_key((fp, "unstable", 0, 1), 1) == (fp, "unstable", 1, 1)


# --------------------------------------------------------------------------- #
# Plan 2.8 -- num_branches is DERIVED from inversion, never passed in
# --------------------------------------------------------------------------- #
def test_num_branches_is_derived_from_inversion():
    assert _bare_fixed_point(3).num_branches == 1
    assert _bare_fixed_point(3, inversion=True).num_branches == 2
    assert _bare_fixed_point(1).num_branches == 1
    assert _bare_fixed_point(1, inversion=True).num_branches == 2


def test_num_branches_needs_the_k_value():
    fp = FixedPoint(3)
    with pytest.raises(ValueError, match="k_value"):
        fp.num_branches


def test_num_branches_agrees_with_get_branch_array():
    for fp in (
        _bare_fixed_point(1),
        _bare_fixed_point(3),
        _bare_fixed_point(2, inversion=True),
    ):
        assert fp.get_branch_array() == list(range(fp.num_branches))


def test_fixed_point_takes_no_branch_count():
    """The constructor derives everything from the period."""
    with pytest.raises(TypeError):
        FixedPoint(3, 2)


def test_both_eigendirection_slots_are_always_allocated():
    """A planar saddle has two eigendirections whether or not the chain uses both.

    ``num_branches`` counts the branches the k-chain visits; the BranchPoint still
    carries a slot for each eigendirection, so a simple saddle can be initialized
    on both sides of its fixed point.
    """
    fp = FixedPoint(3)
    for branch_point in fp.branch_points:
        assert branch_point.num_branches == 2
        assert len(branch_point.forward_branches) == 2
        assert len(branch_point.backward_branches) == 2


# --------------------------------------------------------------------------- #
# Plan 2.5 -- set_k_value guards the eigenvalue signs
# --------------------------------------------------------------------------- #
def test_set_k_value_rejects_disagreeing_eigenvalue_signs():
    """An orientation-reversing map (det J < 0) cannot be modelled by one k_value."""
    fp = FixedPoint(1)
    fp.unstable_eigenvalues = [6.48]
    fp.stable_eigenvalues = [-0.154]
    with pytest.raises(ValueError, match="sign"):
        fp.set_k_value()


def test_set_k_value_tolerates_unset_stable_eigenvalues():
    """Eigenvalues left at their 0.0 placeholder carry no sign to disagree with."""
    fp = FixedPoint(2)
    fp.unstable_eigenvalues = [-2.0, -2.0]
    fp.set_k_value()
    assert fp.k_value == 4

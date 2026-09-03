"""The inversion saddle: the ``k_value == 2 * period`` path, exercised for real.

Every fixed point in the suite before this file has positive eigenvalues, so
``k_value == period`` and the branch index never moves. The ``henon_inversion``
fixture (see ``tests/conftest.py``) is the k=10, b=1 Hénon map's OTHER saddle, at
``(-2.3166, 2.3166)``, whose two eigenvalues are both negative (-4.406 and
-0.227, product 1): a genuine orientation-preserving inversion point on the same
map the rest of the suite uses.

The chain of manifold pieces is therefore ``(0, 0) -> (0, 1) -> (0, 0)``: one map
step lands on the opposite branch, and only two steps return.
"""

from __future__ import annotations

import numpy as np
import pytest

from invariants import (
    assert_cdist_monotonic,
    assert_no_geometric_spikes,
    assert_one_to_one,
)
from tanglepack import TangleWorkbench


# --------------------------------------------------------------------------- #
# (a) advance_key round trips with the branch flip
# --------------------------------------------------------------------------- #
def test_the_fixture_really_is_an_inversion_point(henon_inversion):
    _workbench, fp = henon_inversion

    assert fp.period == 1
    assert fp.check_inversion() is True
    assert fp.k_value == 2
    assert fp.num_branches == 2
    assert fp.unstable_eigenvalues[0] < 0
    assert fp.stable_eigenvalues[0] < 0
    # orientation preserving: det J = lambda_u * lambda_s = 1
    lambda_u = float(np.asarray(fp.unstable_eigenvalues[0]).ravel()[0])
    lambda_s = float(np.asarray(fp.stable_eigenvalues[0]).ravel()[0])
    assert lambda_u * lambda_s == pytest.approx(1.0)


@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_advance_key_flips_the_branch_and_returns_after_k_value_steps(
    henon_inversion, stability
):
    _workbench, fp = henon_inversion
    key = (fp, stability, 0, 0)

    assert fp.advance_key(key, 1) == (fp, stability, 0, 1)
    assert fp.advance_key(key, 2) == key
    assert fp.advance_key(key, -1) == (fp, stability, 0, 1)
    assert fp.branch_cycle(stability) == [key, (fp, stability, 0, 1)]


# --------------------------------------------------------------------------- #
# (b) construct_kevin_way builds both branches, in opposite directions
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_kevin_way_builds_both_branches_in_opposite_directions(
    henon_inversion_initialized, stability
):
    """The second piece is the image of the first, so it leaves along -v."""
    workbench, fp = henon_inversion_initialized

    keys = {k for k in workbench.manifolds if k[1] == stability}
    assert keys == {(fp, stability, 0, 0), (fp, stability, 0, 1)}

    anchor = np.asarray(fp.coordinates[0], dtype=float).ravel()
    directions = []
    for branch_index in (0, 1):
        manifold = workbench.manifolds[(fp, stability, 0, branch_index)]
        assert manifold.manifold_key == (fp, stability, 0, branch_index)
        assert manifold.branch_index == branch_index
        points = manifold.get_point_array()
        step = np.asarray(points[1], dtype=float).ravel() - anchor
        directions.append(step / np.linalg.norm(step))

    eigenvector = np.asarray(
        fp.unstable_eigenvectors[0]
        if stability == "unstable"
        else fp.stable_eigenvectors[0],
        dtype=float,
    ).ravel()
    eigenvector = eigenvector / np.linalg.norm(eigenvector)

    assert abs(float(np.dot(directions[0], eigenvector))) == pytest.approx(1.0, abs=1e-3)
    assert float(np.dot(directions[0], directions[1])) == pytest.approx(-1.0, abs=1e-3)


# --------------------------------------------------------------------------- #
# (c) growth keeps the geometric invariants on BOTH branches
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("stability", ["unstable", "stable"])
@pytest.mark.parametrize("branch_index", [0, 1])
def test_growth_preserves_the_invariants_on_both_branches(
    henon_inversion, stability, branch_index
):
    workbench, fp = henon_inversion
    manifold = workbench.manifolds[(fp, stability, 0, branch_index)]

    assert len(manifold.get_point_array()) > 5, "the branch never grew"
    assert_cdist_monotonic(manifold)
    assert_no_geometric_spikes(manifold)
    assert_one_to_one(manifold)


def _walk(manifold, branch_index):
    """Every node of one branch, root first."""
    nodes = []
    previous, current = manifold.root, manifold.walk_fwd(
        None, manifold.root, branch_index=branch_index
    )
    while current is not None:
        nodes.append(current)
        previous, current = current, manifold.walk_fwd(previous, current)
    return nodes


@pytest.mark.parametrize("branch_index", [0, 1])
def test_one_map_step_scales_cdist_by_per_step_beta(henon_inversion, branch_index):
    """The measured per-iterate cdist ratio IS ``per_step_beta``.

    This is the empirical statement of the whole map-step row: the eigenvalue is
    that of DM^period, so one map step multiplies the unstable canonical distance
    by ``|lambda_u| ** (1 / period)``. Taking the ``k_value``-th root instead
    (the pre-fix formula) is off by a factor of two in the exponent here --
    2.099 against the 4.406 the points actually move -- and no iterate chain in
    the registry matches it.
    """
    workbench, fp = henon_inversion
    manifold = workbench.manifolds[(fp, "unstable", 0, branch_index)]
    beta = fp.per_step_beta("unstable")

    ratios = [
        node.next_iterate.cdist / node.cdist
        for node in _walk(manifold, branch_index)
        if getattr(node, "next_iterate", None) is not None
        and node.cdist
        and node.next_iterate.cdist
    ]

    assert len(ratios) > 10, "no iterate pairs to measure on this branch"

    # The points carry the MEASURED stretch (per_step_factor of the ratio over
    # k_value real map applications), which trails the analytic eigenvalue by the
    # linearization error of the seed step -- ~5e-6 relative here. The wrong
    # k_value-th root is off by 110%, so this tolerance still pins the formula.
    wrong = float(np.abs(np.asarray(fp.unstable_eigenvalues[0]).ravel()[0])) ** (
        1.0 / fp.k_value
    )
    for ratio in ratios:
        assert ratio == pytest.approx(beta, rel=1e-4)
        assert ratio != pytest.approx(wrong, rel=0.5)


def test_forward_iterates_are_registered_on_the_inversion_tangle(henon_inversion):
    """The crossings' iterate chains are found, at the per-map-step factor.

    This is the downstream half of the map-step fix: the registry predicts an
    image at ``(u * beta, s / beta)`` and matches on both canonical distances, so
    a wrong beta finds nothing at all on an inversion tangle.
    """
    workbench, fp = henon_inversion
    registry = workbench._intersection_registry
    beta = fp.per_step_beta("unstable")

    transverse = [
        (source_id, target_id)
        for source_id, n, target_id in registry.iterate_table.items()
        if n == 1 and registry[source_id].unstable_cdist > 0
    ]
    assert len(transverse) >= 2, (
        "no forward iterate of a transverse crossing was registered on the "
        "inversion tangle"
    )
    for source_id, target_id in transverse:
        source, target = registry[source_id], registry[target_id]
        assert target.unstable_cdist / source.unstable_cdist == pytest.approx(
            beta, rel=1e-3
        )
        assert target.stable_cdist / source.stable_cdist == pytest.approx(
            1.0 / beta, rel=1e-3
        )


# --------------------------------------------------------------------------- #
# (d) crossings carry both branch indices
# --------------------------------------------------------------------------- #
def test_intersections_carry_keys_on_both_branches(henon_inversion):
    workbench, fp = henon_inversion
    registry = workbench._intersection_registry

    assert len(registry) > 0, "the inversion tangle produced no crossings"

    # The anchors (cdist 0 on both manifolds) are synthetic crossings that exist
    # on all four (unstable branch, stable branch) pairs by construction, so they
    # would satisfy this vacuously; only real transverse crossings count.
    transverse = [
        ix for _ix_id, ix in registry if ix.unstable_cdist > 0 and ix.stable_cdist > 0
    ]
    assert len(transverse) >= 4, "no transverse crossings in the inversion tangle"

    seen_unstable, seen_stable = set(), set()
    for ix in transverse:
        assert ix.manifold_b_key is not None
        for key, stability in (
            (ix.manifold_a_key, "unstable"),
            (ix.manifold_b_key, "stable"),
        ):
            if key is None:
                continue
            assert key[0] is fp
            assert key[1] == stability
            assert key[3] in (0, 1)
        if ix.manifold_a_key is not None:
            seen_unstable.add(ix.manifold_a_key[3])
        seen_stable.add(ix.manifold_b_key[3])

    assert seen_unstable == {0, 1}
    assert seen_stable == {0, 1}


# --------------------------------------------------------------------------- #
# (e) an orientation-REVERSING map is rejected, not silently mismodelled
# --------------------------------------------------------------------------- #
def _henon_b_negative(point):
    k, b = 10, -1
    x, y = point
    return np.stack([y - k + x**2, -b * x], axis=0)


def _henon_b_negative_inverse(point):
    k, b = 10, -1
    x, y = point
    return np.stack([-y / b, x + k - (y**2) / (b**2)], axis=0)


def _henon_b_negative_jacobian(point):
    k, b = 10, -1
    x, y = point
    return np.array([[2 * x, 1], [-b, 0]])


def test_mixed_eigenvalue_signs_are_rejected():
    """det J = b = -1: one eigenvalue negative, one positive.

    The unstable and stable manifolds then have DIFFERENT inversion status, which
    a single ``k_value`` cannot represent, so ``set_k_value`` must refuse rather
    than silently model the point as if only the unstable side inverted.
    """
    workbench = TangleWorkbench(
        _henon_b_negative, _henon_b_negative_inverse, _henon_b_negative_jacobian
    )
    with pytest.raises(ValueError, match="sign"):
        workbench.construct_fixed_point([3.1623, 3.1623])

"""The inversion saddle: the ``k_value == 2 * period`` path, exercised for real.

Every fixed point in the suite before this file has positive eigenvalues, so
``k_value == period`` and the branch index never moves. The ``henon_inversion``
fixture (see ``tests/conftest.py``) is the k=10, b=1 Hénon map's OTHER saddle, at
``(-2.3166, 2.3166)``, whose two eigenvalues are both negative (-4.406 and
-0.227, product 1): a genuine orientation-preserving inversion point on the same
map the rest of the suite uses.

The chain of manifold pieces is therefore ``(0, 0) -> (0, 1) -> (0, 0)``: one map
step lands on the opposite branch, and only two steps return.

The growth invariants on both branches, the registered forward iterates, the
crossing keys on both branches and the fixture's inversion facts run as laws on
the ``inversion`` case of ``tests/invariants/``.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack.numerics.ManifoldInitializer import _MIN_SEED_STEP


# --------------------------------------------------------------------------- #
# (a) advance_key round trips with the branch flip
# --------------------------------------------------------------------------- #
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
# (c) one map step scales cdist by per_step_beta on BOTH branches
# --------------------------------------------------------------------------- #
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
    # linearization error of the seed step. The seed step is
    # max(accuracy, _MIN_SEED_STEP), so that error is of order _MIN_SEED_STEP
    # relative (~5e-6 here); the tolerance allows 20x that. The wrong k_value-th
    # root is off by 110%, so this tolerance still pins the formula.
    rel = 20 * _MIN_SEED_STEP
    wrong = float(np.abs(np.asarray(fp.unstable_eigenvalues[0]).ravel()[0])) ** (
        1.0 / fp.k_value
    )
    for ratio in ratios:
        assert ratio == pytest.approx(beta, rel=rel)
        assert ratio != pytest.approx(wrong, rel=0.5)



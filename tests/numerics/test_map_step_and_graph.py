"""Plan 2.5 / 2.7 -- one map step on FixedPoint, and the per-branch intersection graph.

2.5: ``FixedPoint`` owns the map step. ``advance_key`` (Phase 1.7) advances a
manifold key, ``per_step_beta`` is the per-map-step canonical-distance factor and
``branch_cycle`` is the chain of manifold pieces one map step visits in turn.
Nothing else may re-derive any of the three.

2.7: adjacency in the intersection graph is per BRANCH. Canonical distances are
measured from each branch's own anchor, so a global cdist sort is meaningless for
a period > 1 orbit -- it interleaves branches and joins crossings that are not
neighbours on any curve.
"""

from __future__ import annotations

import numpy as np
import pytest

from tanglepack import FixedPoint, TangleWorkbench
from tanglepack.numerics.IterateTable import IterateTable
from tanglepack.topology.Pseudoneighbor import forward_unstable_branch_cycle
from tanglepack.topology.StrongPip import forward_stable_branch_cycle


def _bare_fixed_point(
    period: int,
    *,
    inversion: bool = False,
    lambda_u: float = 4.0,
    lambda_s: float | None = None,
) -> FixedPoint:
    """A FixedPoint carrying only the eigen/period data the map step needs."""
    fp = FixedPoint(period)
    sign = -1.0 if inversion else 1.0
    if lambda_s is None:
        lambda_s = 1.0 / lambda_u
    fp.unstable_eigenvalues = [sign * lambda_u] * period
    fp.stable_eigenvalues = [sign * lambda_s] * period
    fp.set_k_value()
    return fp


# --------------------------------------------------------------------------- #
# 2.5 -- per_step_beta
# --------------------------------------------------------------------------- #
def test_per_step_beta_unstable_is_the_period_th_root_of_the_eigenvalue():
    """The stored eigenvalue is that of DM^period, so one map step is its
    ``period``-th root -- NOT its ``k_value``-th root.

    On an inversion point M^period already carries a point onto the opposite
    branch with the full eigenvalue's growth, so the two roots differ by a
    factor of two in the exponent and only this one matches the dynamics.
    """
    for fp in (
        _bare_fixed_point(1),
        _bare_fixed_point(3),
        _bare_fixed_point(3, inversion=True),
    ):
        beta = fp.per_step_beta("unstable")
        assert beta == pytest.approx(4.0 ** (1.0 / fp.period), rel=0, abs=0)
        assert beta > 1.0
        # a branch return is k_value steps and costs lambda ** num_branches
        assert beta ** fp.k_value == pytest.approx(4.0 ** fp.num_branches)


def test_per_step_beta_is_the_k_th_root_only_without_inversion():
    """The old (wrong) k_value-th root and the right one agree iff k == period."""
    plain = _bare_fixed_point(3)
    assert plain.per_step_beta("unstable") == pytest.approx(
        4.0 ** (1.0 / plain.k_value)
    )

    inverted = _bare_fixed_point(3, inversion=True)
    assert inverted.per_step_beta("unstable") != pytest.approx(
        4.0 ** (1.0 / inverted.k_value)
    )


def test_per_step_beta_stable_is_the_reciprocal_on_an_area_preserving_map():
    """det J = 1 means lambda_u * lambda_s = 1, so the two factors are reciprocal.

    Every call site before this row wrote the stable contraction as
    ``stable_cdist / beta_unstable``; the method must reproduce that exactly.
    """
    fp = _bare_fixed_point(3, lambda_u=4.0, lambda_s=0.25)
    assert fp.per_step_beta("stable") == 1.0 / fp.per_step_beta("unstable")
    assert 0.0 < fp.per_step_beta("stable") < 1.0


def test_per_step_beta_stable_uses_the_stable_eigenvalue_when_the_product_is_not_one():
    """A dissipative map (|det J| != 1) has no reciprocal relation to lean on."""
    fp = _bare_fixed_point(2, lambda_u=4.0, lambda_s=0.1)
    assert fp.unstable_eigenvalues[0] * fp.stable_eigenvalues[0] != pytest.approx(1.0)
    assert fp.per_step_beta("stable") == pytest.approx(0.1 ** (1.0 / fp.period))
    assert fp.per_step_beta("stable") != 1.0 / fp.per_step_beta("unstable")


def test_per_step_beta_rejects_an_unknown_stability():
    fp = _bare_fixed_point(1)
    with pytest.raises(ValueError, match="stability"):
        fp.per_step_beta("sideways")


def test_per_step_beta_matches_the_expression_it_replaces(fixed_point):
    """The k=10 saddle has no inversion (k_value == period == 1), so this is bit
    for bit the ``lambda_u ** (1 / k_value)`` the call sites used to write."""
    _workbench, fp = fixed_point
    assert fp.k_value == fp.period
    lambda_u = float(np.abs(np.asarray(fp.unstable_eigenvalues[0]).ravel()[0]))
    assert fp.per_step_beta("unstable") == lambda_u ** (1.0 / fp.k_value)


# --------------------------------------------------------------------------- #
# 2.5 -- branch_cycle
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize(
    "fp",
    [
        _bare_fixed_point(1),
        _bare_fixed_point(3),
        _bare_fixed_point(2, inversion=True),
        _bare_fixed_point(3, inversion=True),
    ],
)
@pytest.mark.parametrize("stability", ["unstable", "stable"])
def test_branch_cycle_is_the_advance_key_chain(fp, stability):
    cycle = fp.branch_cycle(stability)

    assert len(cycle) == fp.k_value
    assert len(set(cycle)) == fp.k_value
    assert cycle[0] == (fp, stability, 0, 0)
    for a, b in zip(cycle, cycle[1:]):
        assert fp.advance_key(a, 1) == b
    # the cycle closes
    assert fp.advance_key(cycle[-1], 1) == cycle[0]


def test_branch_cycle_reproduces_the_topology_branch_orders():
    """The two hand-rolled cycles in topology/ are exactly this one method."""
    for fp in (
        _bare_fixed_point(1),
        _bare_fixed_point(3),
        _bare_fixed_point(2, inversion=True),
    ):
        assert forward_unstable_branch_cycle(fp) == fp.branch_cycle("unstable")
        assert forward_stable_branch_cycle(fp) == fp.branch_cycle("stable")


def test_branch_cycle_rejects_an_unknown_stability():
    fp = _bare_fixed_point(3)
    with pytest.raises(ValueError, match="stability"):
        fp.branch_cycle("sideways")


# --------------------------------------------------------------------------- #
# 2.5 -- one advance rule, one scaling unit
# --------------------------------------------------------------------------- #
def test_workbench_has_no_private_key_advance():
    """``_advance_key_forward`` is gone; ``FixedPoint.advance_key`` is the rule."""
    assert not hasattr(TangleWorkbench, "_advance_key_forward")
    assert not hasattr(TangleWorkbench, "_per_step_beta")


def test_scale_cdist_is_in_map_steps(henon_p3_session):
    """``Trellis.scale_cdist(c, n, ...)`` scales by ONE map step per n."""
    session, fp3, _fp1, _zone = henon_p3_session
    trellis = session.trellis(fp3)

    beta = fp3.per_step_beta("unstable")
    assert trellis.scale_cdist(1.0, 1, "unstable", fp3) == pytest.approx(beta)
    assert trellis.scale_cdist(1.0, 1, "stable", fp3) == pytest.approx(1.0 / beta)

    # k_value map steps is one full branch return, which costs
    # |lambda| ** num_branches -- the eigenvalue itself on this non-inversion orbit.
    lambda_u = trellis.lambda_u(fp3)
    assert fp3.num_branches == 1
    assert trellis.scale_cdist(1.0, fp3.k_value, "unstable", fp3) == pytest.approx(
        lambda_u**fp3.num_branches
    )
    assert trellis.scale_cdist(2.0, 0, "unstable", fp3) == pytest.approx(2.0)


# --------------------------------------------------------------------------- #
# 2.7 -- IterateTable.items()
# --------------------------------------------------------------------------- #
def test_iterate_table_items_yields_every_relation_once():
    table = IterateTable()
    table.register_iterate(3, 1, 7)
    table.register_iterate(7, 2, 11)
    table[11, -1] = 9  # recorded as the forward relation 9 -> 11

    assert sorted(table.items()) == [(3, 1, 7), (7, 2, 11), (9, 1, 11)]
    for source_id, n, target_id in table.items():
        assert table[source_id, n] == target_id


def test_iterate_table_items_is_empty_for_a_fresh_table():
    assert list(IterateTable().items()) == []


# --------------------------------------------------------------------------- #
# 2.7 -- graph adjacency is per branch
# --------------------------------------------------------------------------- #
def _adjacency(graph, stability):
    return [
        (u, v)
        for u, v, d in graph.edges(data=True)
        if d.get("type") == "adjacency" and d.get("stability") == stability
    ]


def test_k10_intersection_graph_is_unchanged(henon_tangle_with_bridges):
    """Pinned before the rebuild: one stable branch, so nothing may move."""
    workbench, _fp = henon_tangle_with_bridges
    graph = workbench.build_intersection_graph()

    assert graph.number_of_nodes() == 8
    assert graph.number_of_edges() == 17
    assert len(_adjacency(graph, "stable")) == 7
    assert len(_adjacency(graph, "unstable")) == 7
    assert (
        sum(1 for _u, _v, d in graph.edges(data=True) if d.get("type") == "iterate")
        == 3
    )


def test_stable_edges_join_consecutive_crossings_on_one_branch(henon_p3_session):
    """A stable edge is a stable ARC: same branch, adjacent in that branch's order."""
    session, _fp3, _fp1, _zone = henon_p3_session
    workbench = session.workbench
    registry = workbench._intersection_registry
    graph = workbench.build_intersection_graph()

    per_branch: dict[tuple, list[int]] = {}
    for ix_id, ix in registry:
        if ix.manifold_b_key is not None:
            per_branch.setdefault(ix.manifold_b_key, []).append(ix_id)
    for ids in per_branch.values():
        ids.sort(key=lambda i: registry[i].stable_cdist)
    rank = {i: r for ids in per_branch.values() for r, i in enumerate(ids)}

    edges = _adjacency(graph, "stable")
    assert edges, "the period-3 session must have stable adjacency edges"
    for u, v in edges:
        assert registry[u].manifold_b_key == registry[v].manifold_b_key, (
            f"stable edge {u}->{v} joins two different stable branches: "
            f"{registry[u].manifold_b_key[1:]} vs {registry[v].manifold_b_key[1:]}"
        )
        assert abs(rank[u] - rank[v]) == 1, (
            f"stable edge {u}->{v} joins non-consecutive crossings on the branch "
            f"(ranks {rank[u]} and {rank[v]})"
        )

    expected = sum(len(ids) - 1 for ids in per_branch.values())
    assert len(edges) == expected


def test_unstable_edges_are_exactly_the_bridges(henon_p3_session):
    session, _fp3, _fp1, _zone = henon_p3_session
    workbench = session.workbench
    graph = workbench.build_intersection_graph()

    bridge_pairs = sorted(
        (b.first_intersection, b.second_intersection)
        for b in workbench.bridges
        if b.first_intersection is not None
        and b.second_intersection is not None
        and b.first_intersection != b.second_intersection
    )
    assert sorted(_adjacency(graph, "unstable")) == bridge_pairs


def test_registry_graph_falls_back_to_per_branch_unstable_order(henon_p3_session):
    """Without bridges the registry orders the unstable side by its own branch."""
    session, _fp3, _fp1, _zone = henon_p3_session
    registry = session.workbench._intersection_registry
    graph = registry.graph()

    per_branch: dict[tuple, list[int]] = {}
    for ix_id, ix in registry:
        if ix.manifold_a_key is not None:
            per_branch.setdefault(ix.manifold_a_key, []).append(ix_id)
    for ids in per_branch.values():
        ids.sort(key=lambda i: registry[i].unstable_cdist)

    expected = sorted(
        (ids[i], ids[i + 1]) for ids in per_branch.values() for i in range(len(ids) - 1)
    )
    assert sorted(_adjacency(graph, "unstable")) == expected


def test_registry_graph_takes_the_bridges_it_is_given(henon_p3_session):
    session, _fp3, _fp1, _zone = henon_p3_session
    registry = session.workbench._intersection_registry
    ids = registry.all_ids()[:4]
    bridges = [(ids[0], ids[1]), (ids[2], ids[3])]

    graph = registry.graph(bridges=bridges)
    assert sorted(_adjacency(graph, "unstable")) == sorted(bridges)


def test_build_intersection_graph_does_not_mutate_the_registry_graph(
    henon_tangle_with_bridges,
):
    """The workbench decorates a COPY: its extra attributes stay out of the registry."""
    workbench, _fp = henon_tangle_with_bridges
    registry = workbench._intersection_registry

    decorated = workbench.build_intersection_graph()
    decorated.add_node("sentinel")
    decorated.add_edge("sentinel", "sentinel", type="adjacency", stability="stable")

    assert "sentinel" not in registry.graph()

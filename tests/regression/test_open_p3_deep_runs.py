"""The deep period-3 runs: a known open problem, pinned where it is firm.

At 15 or 16 unstable steps, or after 6 blasts, the period-3 tangle's symbolic
dynamics is unreliable: one class is unreachable by the dual-graph walk
(``A:p3@2.0/R#2 <-> A:p3@2.0/R#3``) and a virtual ``new1`` appears (Phase 0
probe P3; ``KNOWN_ISSUES`` key ``open_p3_deep_runs``). What IS firm there:

* invariant I1 -- all holes of one origin share their side of their bridge.
  A bridge's unstable branch is carried by ``Bridge.manifold_key``; inferring it
  from the endpoint intersections instead was unreliable on a period > 1 orbit
  (every branch's anchor sits at unstable cdist 0), so the backward propagation
  picked a container on the wrong branch and an orbit's holes flipped
  ``bridge_side`` partway down the chain. The deep runs are where a propagated
  hole's image sub-arc and the nearest vertex of its whole containing bridge sit
  on different folds (fixed 2026-10-02). The symbolic dynamics must still run.
* virtual symbols are sinks of the transition graph.

That every class resolves without a virtual symbol is ``xfail(strict=True)``.
"""

from __future__ import annotations

import pytest

from cases import BUILDERS, Case, issue_marks
from tanglepack.topology import check_holes_share_bridge_side

#: The deep runs (``tests/cases.py``): 15 and 16 unstable steps, 6 blasts.
DEEP_RUNS: tuple[str, ...] = ("p3_15", "p3_16", "p3_6_blasts")


def _build(name: str) -> Case:
    """A fresh build of one deep run."""
    return BUILDERS[name]()


@pytest.mark.parametrize("name", DEEP_RUNS)
def test_deep_p3_holes_share_bridge_side(name: str) -> None:
    """I1 holds deep into the tangle, and the symbolic dynamics runs."""
    case = _build(name)
    (fp,) = case.fixed_points
    trellis = case.session.trellis(fp)
    assert trellis.holes, f"{name} must punch holes"
    check_holes_share_bridge_side(
        trellis.holes, orientation_preserving=trellis.orientation_preserving
    )
    assert case.session.symbolic_dynamics().classes


@pytest.mark.parametrize(
    "name",
    [
        pytest.param(name, marks=issue_marks("open_p3_deep_runs", name), id=name)
        for name in DEEP_RUNS
    ],
)
def test_deep_p3_every_class_resolves_without_virtual_symbols(name: str) -> None:
    """Every class gets an itinerary and no pair falls outside the class table."""
    dyn = _build(name).session.symbolic_dynamics()
    assert not dyn.virtual_classes
    assert not [cd for cd in dyn.classes.values() if cd.unresolved_reason]


def test_virtual_symbols_are_transition_sinks() -> None:
    """A virtual symbol has no outgoing transition and a zero matrix row (real data:
    the 15-step run, unreliable as it is)."""
    dyn = _build("p3_15").session.symbolic_dynamics()
    graph = dyn.transition_graph()
    virtual = [node for node, data in graph.nodes(data=True) if data["kind"] == "virtual"]
    if not virtual:
        pytest.skip("the 15-step run produced no virtual symbol (the open problem is gone)")
    names, matrix = dyn.transition_matrix()
    for node in virtual:
        assert graph.out_degree(node) == 0
        assert graph.in_degree(node) > 0
        if node in names:
            assert not matrix[names.index(node)].any()

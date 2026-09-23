"""
TangleSession.symbolic_dynamics / describe_symbolic_dynamics and the plot
delegates: results, caching, invalidation, and the two-blast k=2.8 pin.

Mirrors ``test_session_dual_graph.py``. Nothing here pins a registry id, a
bridge count or a class letter: the two-blast structure is located through
the homotopy names of each class's elements.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from minimal_helpers import build_pieces
from tanglepack.topology import plotting
from tanglepack.topology.DualGraph import DualGraph
from tanglepack.topology.SymbolicDynamics import SymbolicDynamics, symbolic_dynamics


def _rule_texts(dyn: SymbolicDynamics) -> dict[str, str]:
    """``{letter: word}`` over unrefined symbols, as text."""
    return {letter: dyn.word(letter, refined=False) for letter in dyn.rules}


def _refined_rule_texts(dyn: SymbolicDynamics) -> dict[str, str]:
    return {name: dyn.word(name) for name in dyn.refined_rules}


def _same_dynamics(a: SymbolicDynamics, b: SymbolicDynamics) -> bool:
    """Word-level equivalence: the same rules, refined rules and matrix."""
    names_a, matrix_a = a.transition_matrix()
    names_b, matrix_b = b.transition_matrix()
    return (
        _rule_texts(a) == _rule_texts(b)
        and _refined_rule_texts(a) == _refined_rule_texts(b)
        and names_a == names_b
        and np.array_equal(matrix_a, matrix_b)
        and a.is_reliable == b.is_reliable
    )


def _homotopy_names(dyn: SymbolicDynamics, cd) -> frozenset[str]:
    """The homotopy names of a class's two elements, e.g. ``{"R_1", "R_3"}``."""
    return frozenset(
        dyn.naming.homotopy_name(ref).text
        for ref in (cd.bridge_class.source, cd.bridge_class.target)
    )


def _class_named(dyn: SymbolicDynamics, *names: str):
    """The active class whose elements carry exactly these homotopy names."""
    wanted = frozenset(names)
    matches = [cd for cd in dyn.classes.values() if _homotopy_names(dyn, cd) == wanted]
    assert len(matches) == 1, f"expected one class {sorted(wanted)}, found {matches}"
    return matches[0]


# --------------------------------------------------------------------------- #
# k=10: results and caching
# --------------------------------------------------------------------------- #
def test_k10_symbolic_dynamics_matches_a_direct_build(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    dual = DualGraph(pieces.minimal, pieces.iterated, strong_pips=pieces.strong_pips)
    direct = symbolic_dynamics(dual, session.bridge_classes([fp]))
    dyn = session.symbolic_dynamics([fp])
    assert isinstance(dyn, SymbolicDynamics)
    assert _same_dynamics(dyn, direct)
    assert dyn.rules, "the k=10 fixture has at least one resolved class"


def test_k10_symbolic_dynamics_is_built_over_the_cached_pieces(k10_partitioned):
    session, fp = k10_partitioned
    dyn = session.symbolic_dynamics([fp])
    assert dyn.table is session.bridge_classes([fp])
    assert dyn.naming.iterated is session.iterated_partition([fp])
    active = [cd for cd in dyn.classes.values() if cd.kind == "active"]
    assert len(active) == 1
    letter = active[0].letter
    assert active[0].word == f"{letter} u^-1 {letter}^-1"
    assert dyn.is_reliable


def test_cache_hits_return_the_same_object(k10_partitioned):
    session, fp = k10_partitioned
    assert session.symbolic_dynamics() is session.symbolic_dynamics()
    assert session.symbolic_dynamics(fp) is session.symbolic_dynamics(fp)
    assert session.symbolic_dynamics([fp]) is session.symbolic_dynamics([fp])


def test_rebuild_flag_forces_a_fresh_equivalent_object(k10_partitioned):
    session, fp = k10_partitioned
    before = session.symbolic_dynamics()
    dual_before = session.dual_graph()
    after = session.symbolic_dynamics(rebuild=True)
    assert after is not before and _same_dynamics(before, after)
    # A rebuild rebuilds the dual graph and the table it reads too.
    assert session.dual_graph() is not dual_before
    assert after.table is session.bridge_classes()
    assert session.symbolic_dynamics() is after


def test_workbench_mutation_invalidates_the_cache(k10_partitioned):
    session, fp = k10_partitioned
    before = session.symbolic_dynamics()
    # Bump the workbench generation and restore a fully partitioned state, so
    # the post-mutation call recomputes successfully rather than raising.
    session.grow_n_times(fp, "unstable", num_iterations=1)
    session.compute_intersections([fp], preserve_ids=True)
    session.create_bridges(fp)
    session.classify_strong_pips()
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()
    after = session.symbolic_dynamics()
    assert after is not before
    assert after.table is session.bridge_classes()
    for cd in after.classes.values():
        if cd.itinerary is not None:
            assert len(cd.itinerary) % 2 == 0


def test_pip_change_at_the_same_generation_invalidates_the_cache(k10_partitioned):
    session, fp = k10_partitioned
    trellis = session.trellis(fp)
    alternatives = [c for c in trellis.strong_pip_candidates if c != trellis.strong_pip]
    if not alternatives:
        pytest.skip("the fixture has a single strong-pip candidate")
    generation = session.workbench.generation
    before = session.symbolic_dynamics()
    table_before = session.bridge_classes()
    trellis.set_strong_pip(alternatives[0])
    assert session.workbench.generation == generation
    after = session.symbolic_dynamics()
    assert after is not before
    # The classes do not depend on the pip choice; the dual graph does.
    assert session.bridge_classes() is table_before
    assert after.table is table_before


def test_describe_symbolic_dynamics(k10_partitioned):
    session, fp = k10_partitioned
    text = session.describe_symbolic_dynamics([fp])
    assert text and isinstance(text, str)
    assert text == session.symbolic_dynamics([fp]).describe()
    active = next(cd for cd in session.symbolic_dynamics([fp]).classes.values()
                  if cd.kind == "active")
    letter = active.letter
    assert f"{letter} -> {letter} u^-1 {letter}^-1" in text
    assert "transition matrix" in text


# --------------------------------------------------------------------------- #
# Plot delegates
# --------------------------------------------------------------------------- #
_PLOT_NAMES = (
    "plot_walk",
    "plot_bridges_by_class",
    "plot_transition_graph",
    "plot_itinerary_table",
    "walk_legend_handles",
)
_plots_missing = [name for name in _PLOT_NAMES if not hasattr(plotting, name)]
needs_plots = pytest.mark.skipif(
    bool(_plots_missing), reason=f"plotting is missing {_plots_missing}"
)


@needs_plots
def test_session_plot_delegates_draw_on_the_given_axes(k10_partitioned):
    """The delegates run on Agg and hand back what the plotting functions return."""
    import networkx
    from matplotlib.lines import Line2D
    from matplotlib.table import Table

    session, fp = k10_partitioned
    dyn = session.symbolic_dynamics([fp])
    active = next(cd for cd in dyn.classes.values() if cd.kind == "active")
    assert active.search is not None and active.search.walk is not None
    fig, (ax_graph, ax_table, ax_plane) = plt.subplots(1, 3)
    try:
        graph = session.plot_transition_graph(ax=ax_graph, fixed_points=[fp])
        assert isinstance(graph, networkx.DiGraph)
        assert isinstance(
            session.plot_transition_graph(dyn, ax=ax_graph, refined=False), networkx.DiGraph
        )
        table = session.plot_itinerary_table(ax=ax_table, fixed_points=[fp])
        assert isinstance(table, Table) and table.axes is ax_table
        assert isinstance(session.plot_itinerary_table(dyn, ax=ax_table), Table)
        colours = session.plot_bridges_by_class(ax=ax_plane, fixed_points=[fp])
        assert isinstance(colours, dict) and colours
        assert session.plot_bridges_by_class(dyn, ax=ax_plane, refined=False)
        assert session.plot_dual_graph(session.dual_graph([fp]), ax=ax_plane) is ax_plane
        line = session.plot_walk(active.search.walk, ax=ax_plane, fixed_points=[fp])
        assert isinstance(line, Line2D) and line.axes is ax_plane
        assert ax_plane.lines, "the bridges and the walk are drawn on the plane axes"
    finally:
        plt.close(fig)


@needs_plots
def test_session_plot_delegates_forward_to_plotting(k10_partitioned, monkeypatch):
    session, fp = k10_partitioned
    dyn = session.symbolic_dynamics()
    active = next(cd for cd in dyn.classes.values() if cd.kind == "active")
    walk = active.search.walk
    calls = []
    monkeypatch.setattr(
        plotting, "plot_transition_graph",
        lambda d, ax=None, **kw: calls.append(("graph", d, ax, kw)) or "graph",
    )
    monkeypatch.setattr(
        plotting, "plot_walk",
        lambda g, w, ax=None, **kw: calls.append(("walk", (g, w), ax, kw)) or "walk",
    )
    monkeypatch.setattr(
        plotting, "plot_bridges_by_class",
        lambda t, d, ax=None, **kw: calls.append(("bridges", (t, d), ax, kw)) or "bridges",
    )
    monkeypatch.setattr(
        plotting, "plot_itinerary_table",
        lambda d, ax=None, **kw: calls.append(("table", d, ax, kw)) or "table",
    )
    monkeypatch.setattr(
        plotting, "plot_dual_graph_cartoon",
        lambda g, d, ax=None, **kw: calls.append(("cartoon", (g, d), ax, kw)) or "cartoon",
    )
    assert session.plot_transition_graph(refined=False) == "graph"
    assert session.plot_walk(walk, color="k") == "walk"
    assert session.plot_bridges_by_class(show_inert=False) == "bridges"
    assert session.plot_itinerary_table(refined=True) == "table"
    assert calls[0][1] is dyn and calls[0][3] == {"refined": False}
    assert calls[1][1] == (session.dual_graph(), walk) and calls[1][3] == {"color": "k"}
    assert calls[2][1] == (session.trellis(), dyn) and calls[2][3] == {"show_inert": False}
    assert calls[3][1] is dyn and calls[3][3] == {"refined": True}
    assert session.plot_dual_graph_cartoon(show_walks=False) == "cartoon"
    assert calls[4][1] == (session.dual_graph(), dyn) and calls[4][3] == {"show_walks": False}


# --------------------------------------------------------------------------- #
# k=2.8, two blasts: the three-class structure, pinned by element names
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_k28_two_blasts_structure(k28_two_blasts_partitioned):
    session, fp = k28_two_blasts_partitioned
    dyn = session.symbolic_dynamics([fp])
    assert dyn is session.symbolic_dynamics([fp])
    assert dyn.is_reliable, dyn.describe()

    active = [cd for cd in dyn.classes.values() if cd.kind == "active"]
    assert len(active) == 3, dyn.describe()
    for cd in dyn.classes.values():
        assert cd.itinerary is not None, cd.unresolved_reason
        assert len(cd.itinerary) % 2 == 0
        assert cd.source == "walk"
        assert cd.verified is True
    inert = [cd for cd in dyn.classes.values() if cd.kind == "inert"]
    assert [_homotopy_names(dyn, cd) for cd in inert] == [frozenset({"L_1", "L_3"})]
    inert_letter = inert[0].letter
    assert inert[0].symbols == []

    # Letters follow the smallest member cdist: the anchor class {R_1, R_5} is
    # a, then {R_3, R_5} is b and {R_1, R_3} is c; the inert {L_1, L_3} is u.
    a = _class_named(dyn, "R_1", "R_5")
    b = _class_named(dyn, "R_3", "R_5")
    c = _class_named(dyn, "R_1", "R_3")
    assert (a.letter, b.letter, c.letter, inert_letter) == ("a", "b", "c", "u")
    table = session.bridge_classes([fp])
    assert table.entries[0].bridge_class == a.bridge_class
    assert table.entries[0].min_unstable_cdist == 0.0
    assert [entry.min_unstable_cdist for entry in table] == sorted(
        entry.min_unstable_cdist for entry in table
    )

    # c -> a u^-1 a^-1: three symbols, X then inert^-1 then X^-1 with X = {R_1, R_5}.
    assert [(s.bridge_class, s.direction) for s in c.symbols] == [
        (a.bridge_class, +1),
        (inert[0].bridge_class, -1),
        (a.bridge_class, -1),
    ]
    assert c.word == f"{a.letter} {inert_letter}^-1 {a.letter}^-1"
    # a -> a u^-1 b^-1.
    assert [(s.bridge_class, s.direction) for s in a.symbols] == [
        (a.bridge_class, +1),
        (inert[0].bridge_class, -1),
        (b.bridge_class, -1),
    ]
    # b -> c: the single pair (R_1^3, R_3).
    assert [(s.bridge_class, s.direction) for s in b.symbols] == [(c.bridge_class, +1)]
    assert b.word == c.letter
    assert len(b.itinerary) == 2

    # Only {R_1, R_5} refines, into exactly two children (anchor outward), and
    # under the footprint cut rule every member of it matches a child.
    assert set(dyn.refined) == {a.bridge_class}
    children = dyn.refined[a.bridge_class]
    assert [child.name for child in children] == [f"{a.letter}_1", f"{a.letter}_2"]
    assert dyn.word(f"{a.letter}_1") == dyn.word(f"{a.letter}_2")
    assert dyn.word(c.letter) == f"{a.letter}_1 {inert_letter}^-1 {a.letter}_2^-1"
    assert dyn.unmatched_members == {}
    assert {bid for child in children for bid in child.members} == {
        member.bridge_id for member in a.entry.members if not member.is_loop
    }

    # Unrefined matrix over (a, b, c) = [[1,1,0],[0,0,1],[2,0,0]].
    names, matrix = dyn.transition_matrix(refined=False)
    order = [names.index(cd.letter) for cd in (a, b, c)]
    restricted = matrix[np.ix_(order, order)]
    assert restricted.tolist() == [[1, 1, 0], [0, 0, 1], [2, 0, 0]]
    assert matrix[names.index(inert_letter)].tolist() == [0] * len(names)

    assert session.describe_symbolic_dynamics([fp]) == dyn.describe()


# --------------------------------------------------------------------------- #
# nested period-3: smoke only (no pinned words; the trellis is not grown far
# enough for every class to resolve)
# --------------------------------------------------------------------------- #


@pytest.mark.slow
def test_p3_symbolic_dynamics_smoke(p3_partitioned):
    """Over both fixed points the dynamics builds and every itinerary is well formed."""
    session, _fp3, _fp1 = p3_partitioned
    dyn = session.symbolic_dynamics()
    assert isinstance(dyn, SymbolicDynamics)
    assert dyn.naming.tagged  # more than one stable branch is partitioned
    assert len(dyn.classes) > 0
    for cd in dyn.classes.values():
        if cd.itinerary is None:
            assert cd.unresolved_reason is not None
            continue
        assert len(cd.itinerary) % 2 == 0
        for first, second in zip(cd.itinerary[0::2], cd.itinerary[1::2]):
            assert first.side == second.side
        assert not any(symbol.cross_side for symbol in cd.symbols)
    names, matrix = dyn.transition_matrix()
    assert matrix.shape == (len(names), len(names))
    assert dyn.describe()

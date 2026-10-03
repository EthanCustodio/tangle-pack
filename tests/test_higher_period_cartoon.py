"""
Higher-period and nested tangles: the shared case builders, the cross-branch
fixes they need, and the circular-zone dual-graph cartoon.

The period-3 and nested sessions come from
``tanglepack.examples.henon_cases`` -- the same builders the figure scripts
use -- and are session-scoped (the tests only read them).
"""

from __future__ import annotations

import importlib
import logging
from types import SimpleNamespace

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pytest

from minimal_helpers import build_pieces
from tanglepack.examples.henon_cases import build_nested, build_period3
from tanglepack.topology import plotting
from tanglepack.topology.DualGraph import DualGraph
from tanglepack.topology.ElementNaming import ElementNaming
from tanglepack.topology.SymbolicDynamics import symbolic_dynamics
from tanglepack.topology.TopologyResults import PartitionInterval

family_module = importlib.import_module("tanglepack.topology.PartitionFamily")


# --------------------------------------------------------------------------- #
# Fixtures
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="session")
def p3_built():
    """The period-3 orbit alone (read only)."""
    return build_period3()


@pytest.fixture(scope="session")
def nested_built():
    """The nested period-1 + period-3 tangle, outer zone blasted twice (read only)."""
    return build_nested()


def _k10_dynamics(k10_partitioned):
    """``(pieces, dual graph, symbolic dynamics)`` of the k=10 fixture."""
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    if any(entry.letter is None for entry in pieces.table.active):
        pieces.table = session.bridge_classes([fp])
    dual = DualGraph(pieces.minimal, pieces.iterated, strong_pips=pieces.strong_pips)
    return pieces, dual, symbolic_dynamics(dual, pieces.table)


# --------------------------------------------------------------------------- #
# Chord pairing across stable branches (synthetic)
# --------------------------------------------------------------------------- #
class _FakeTrellis:
    """Just enough trellis for ``_empty_stretches``: a crossing's stable branch."""

    def __init__(self, branches: dict[int, tuple]):
        self._branches = branches

    def intersection(self, iid: int):
        return SimpleNamespace(manifold_b_key=self._branches[iid])


def _hole_homotopy(lo_id: int, hi_id: int, branch) -> SimpleNamespace:
    """A homotopy family with one open-open interval, the hole ``(lo, hi)``."""
    hole = PartitionInterval(
        lo_id=lo_id, hi_id=hi_id, lo_cdist=1.0, hi_cdist=2.0,
        closed_lo=False, closed_hi=False,
    )
    return SimpleNamespace(results={(branch, "right"): SimpleNamespace(intervals=[hole])})


def test_chords_pair_only_the_base_branchs_crossings(monkeypatch, caplog):
    own, foreign = ("fp", "stable", 1, 0), ("fp", "stable", 2, 0)
    branches = {1: own, 2: own, 10: own, 11: own, 12: own, 13: own, 98: foreign, 99: foreign}
    minimal = SimpleNamespace(trellis=_FakeTrellis(branches), hole_bridge_ids=[(1, 2)])
    # The image of the hole bridge crosses the foreign branch twice between
    # its two folds on its own branch.
    chain = [(10, 11), (11, 99), (99, 98), (98, 12), (12, 13)]
    monkeypatch.setattr(family_module, "_image_chain", lambda trellis, bid: chain)
    with caplog.at_level(logging.INFO, logger="tanglepack.topology.PartitionFamily"):
        empty = family_module.IteratedHomotopyPartition._empty_stretches(
            minimal, _hole_homotopy(1, 2, own)
        )
    assert (13, "base") in empty[10]
    assert (12, "chord") in empty[11]
    assert 99 not in empty and 98 not in empty
    assert any("foreign stable crossing" in r.message for r in caplog.records)


def test_a_lobe_ending_on_another_branch_marks_nothing(monkeypatch, caplog):
    own, foreign = ("fp", "stable", 1, 0), ("fp", "stable", 2, 0)
    branches = {1: own, 2: own, 10: own, 11: own, 12: foreign}
    minimal = SimpleNamespace(trellis=_FakeTrellis(branches), hole_bridge_ids=[(1, 2)])
    monkeypatch.setattr(family_module, "_image_chain", lambda trellis, bid: [(10, 11), (11, 12)])
    with caplog.at_level(logging.WARNING, logger="tanglepack.topology.PartitionFamily"):
        empty = family_module.IteratedHomotopyPartition._empty_stretches(
            minimal, _hole_homotopy(1, 2, own)
        )
    assert set(empty) == {1, 2}  # only the hole itself
    assert any("ends on another stable branch" in r.message for r in caplog.records)


# --------------------------------------------------------------------------- #
# Circular-zone cartoon, k=10 (one branch)
# --------------------------------------------------------------------------- #
def test_circle_layout_puts_every_node_on_its_sides_circle(k10_partitioned):
    _pieces, dual, dynamics = _k10_dynamics(k10_partitioned)
    fig, ax = plt.subplots()
    try:
        layout = plotting.plot_dual_graph_cartoon(dual, dynamics, ax=ax, shape="circle")
        assert isinstance(layout, plotting.ZoneLayout)
        ((fp_id, (cx, cy, radius)),) = layout.circles.items()
        (branch_key,) = layout.arcs
        _fp, alpha, sweep = layout.arcs[branch_key]
        assert alpha == 0.0 and sweep == pytest.approx(np.pi)  # period 1: upper half
        assert set(layout.nodes) == set(layout.normals) == set(layout.segments)
        for ref, (x, y) in layout.nodes.items():
            sigma = layout.side_sign(branch_key, ref.side)
            offset = (
                plotting.CARTOON_SEGMENT_OFFSET if ref in layout.singletons
                else plotting.CARTOON_NODE_OFFSET
            )
            assert np.hypot(x - cx, y - cy) == pytest.approx(radius + sigma * offset)
            normal = np.array(layout.normals[ref])
            radial = np.array([x - cx, y - cy]) / np.hypot(x - cx, y - cy)
            assert normal @ radial == pytest.approx(sigma)
        # Anchorward ends sit at smaller angles: the anchor is the clockwise end.
        for ref, (t_lo, t_hi, _r) in layout.segments.items():
            assert t_lo <= t_hi
            assert alpha <= t_lo and t_hi <= alpha + sweep
        # The same bridges and the default sides (right inside) as without zones.
        bridge_gids = {p.get_gid() for p in ax.patches if (p.get_gid() or "").startswith("bridge:")}
        assert len(bridge_gids) == layout.bridges_drawn == len(dual.minimal.kept_bridge_ids)
        assert layout.outside[branch_key] == "left"
        assert "0.0" in {t.get_text() for t in ax.texts}
    finally:
        plt.close(fig)


def test_circle_layout_honours_interior_side_and_rejects_bad_shapes(k10_partitioned):
    _pieces, dual, dynamics = _k10_dynamics(k10_partitioned)
    (branch_key,) = dual.partition.branch_keys
    layout = plotting.dual_graph_zone_layout(dual, interior_side={branch_key: "left"})
    assert layout.outside[branch_key] == "right"
    fig, ax = plt.subplots()
    try:
        with pytest.raises(ValueError, match="shape"):
            plotting.plot_dual_graph_cartoon(dual, dynamics, ax=ax, shape="hexagon")
        with pytest.raises(ValueError, match="label_position"):
            plotting.plot_dual_graph_cartoon(
                dual, dynamics, ax=ax, shape="circle", label_position="nowhere"
            )
    finally:
        plt.close(fig)
    labels = [h.get_label() for h in plotting.dual_graph_cartoon_legend_handles(
        walks=False, empty_side=True
    )]
    assert labels[-1] == "empty (unstable) side"


def test_line_cartoon_labels_every_row_with_its_branch_code(k10_partitioned):
    _pieces, dual, dynamics = _k10_dynamics(k10_partitioned)
    fig, ax = plt.subplots()
    try:
        plotting.plot_dual_graph_cartoon(dual, dynamics, ax=ax)
        assert "anchor  0.0" in {t.get_text() for t in ax.texts}
    finally:
        plt.close(fig)


# --------------------------------------------------------------------------- #
# Period 3
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_p3_words_are_equivariant_under_the_orbit_shift(p3_built):
    """The three active classes map one to the next, oriented alike (A2)."""
    dynamics = p3_built.session.symbolic_dynamics()
    active = [cd for cd in dynamics.classes.values() if cd.kind == "active"]
    assert len(active) == 3
    # Every active class starts at its anchor element, on whichever branch.
    for cd in active:
        assert cd.bridge_class.source.element_id == 0
        assert cd.bridge_class.target.element_id > 0
    words = {cd.letter: cd.word for cd in active}
    a, b, c = sorted(words)
    assert words[a] == b and words[b] == c
    assert words[c].split()[0] == a  # no orbit-shift inverse on the anchor bridge
    # Everything resolves and no walk is unreachable.
    for cd in dynamics.classes.values():
        assert cd.itinerary is not None, cd.unresolved_reason
        assert cd.search is None or cd.search.status != "unreachable"
    assert dynamics.is_reliable


@pytest.mark.slow
def test_p3_names_carry_orbit_codes_without_a_letter(p3_built):
    naming = p3_built.session.symbolic_dynamics().naming
    assert not naming.letters
    codes = {name.branch_code for name in naming.names}
    assert codes == {"0.0", "1.0", "2.0"}


@pytest.mark.slow
def test_p3_circle_follows_the_zone_boundary(p3_built):
    """Ring order, direction and anchor positions are the zone's own."""
    session = p3_built.session
    (fp,) = p3_built.fixed_points
    interior, parents, arc_order = session.cartoon_zones()
    assert parents == {}
    assert set(interior) == set(fp.branch_cycle("stable"))
    order = arc_order[id(fp)]
    # The map sends z0 -> z1 -> z2 clockwise in the plane, and the boundary
    # runs pip_j -> z_j -> pip_{j+1} clockwise (its interior on the right).
    assert order.clockwise
    assert set(interior.values()) == {"right"}
    ring = [key[2] for key in order.order]
    assert ring in ([0, 1, 2], [1, 2, 0], [2, 0, 1])
    fig, ax = plt.subplots()
    try:
        layout = session.plot_dual_graph_cartoon(ax=ax, shape="circle")
        sweep = np.pi / 3
        keys = order.order
        for this, following in zip(keys, keys[1:] + keys[:1]):
            _fp, alpha, signed = layout.arcs[this]
            _fp, alpha_next, signed_next = layout.arcs[following]
            # Elements run counter-clockwise from the anchor (the clockwise
            # end), and going clockwise the empty arc leaves the anchor for
            # the NEXT branch's pip end.
            assert signed == pytest.approx(sweep)
            pip_end_next = alpha_next + signed_next
            gap = (alpha - pip_end_next) % (2 * np.pi)
            assert gap == pytest.approx(sweep)
            # Each anchor sits within half a slot of its orbit point.
            real = order.anchor_angles[this]
            miss = abs((alpha - real + np.pi) % (2 * np.pi) - np.pi)
            assert miss < sweep
        assert {"0.0", "1.0", "2.0"} <= {t.get_text() for t in ax.texts}
    finally:
        plt.close(fig)


# --------------------------------------------------------------------------- #
# Nested
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_nested_resolves_everything_with_real_walks(nested_built):
    session = nested_built.session
    dynamics = session.symbolic_dynamics()
    assert all(cd.itinerary is not None for cd in dynamics.classes.values())
    statuses = {cd.search.status for cd in dynamics.classes.values() if cd.search is not None}
    assert "unique" in statuses and "unreachable" not in statuses
    assert session.minimal_trellis().image_bridge_ids


@pytest.mark.slow
def test_nested_names_carry_fixed_point_letters(nested_built):
    outer, inner = nested_built.fixed_points
    naming = nested_built.session.symbolic_dynamics().naming
    assert naming.letters == {id(outer): outer.label, id(inner): inner.label}
    assert {outer.label, inner.label} == {"A", "B"}
    for ref, name in naming.items():
        assert name.fixed_point_letter == ref.fixed_point.label
        assert name.mathtext.startswith(f"${{}}^{{{ref.fixed_point.label}}}")
        assert ref.label.startswith(f"{ref.fixed_point.label}:")
    assert isinstance(naming, ElementNaming)


@pytest.mark.slow
def test_nested_inner_circle_sits_inside_the_outer_one(nested_built):
    session = nested_built.session
    outer, inner = nested_built.fixed_points
    _interior, parents, _arc_order = session.cartoon_zones()
    assert parents == {id(inner): id(outer)}
    fig, ax = plt.subplots()
    try:
        layout = session.plot_dual_graph_cartoon(ax=ax, shape="circle")
        ox, oy, o_radius = layout.circles[id(outer)]
        ix, iy, i_radius = layout.circles[id(inner)]
        assert (ix, iy) == (ox, oy)
        assert i_radius < o_radius
        assert layout.parents == {id(inner): id(outer)}
        texts = {t.get_text() for t in ax.texts}
        assert {f"{inner.label} 0.0", f"{inner.label} 1.0", f"{outer.label} 0.0"} <= texts
    finally:
        plt.close(fig)


# --------------------------------------------------------------------------- #
# The trellis panel draws the minimal trellis only
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_bridges_by_class_can_be_restricted_to_the_minimal_trellis(nested_built):
    session = nested_built.session
    dynamics = session.symbolic_dynamics()
    minimal = session.minimal_trellis()
    members = {
        member.bridge_id
        for cd in dynamics.classes.values()
        for member in cd.entry.members
    }
    kept = set(minimal.kept_bridge_ids)
    assert members - kept, "the blasts add class members outside the minimal trellis"
    trellis = session.trellis()
    fig, ax = plt.subplots()
    try:
        before = len(ax.lines)
        plotting.plot_bridges_by_class(trellis, dynamics, ax=ax, bridge_ids=kept)
        drawn = len(ax.lines) - before
        expected = sum(
            1 for bid in members & kept if trellis.bridge_between(*bid) is not None
        )
        assert drawn == expected
    finally:
        plt.close(fig)


@pytest.mark.slow
def test_nested_classes_are_lettered_and_coloured_tangle_by_tangle(nested_built):
    session = nested_built.session
    dynamics = session.symbolic_dynamics()
    table = session.bridge_classes()
    groups = [entry.tangle for entry in table]
    ranks = [len(table.entries) if g is None else g for g in groups]
    assert ranks == sorted(ranks), "one tangle's classes, then the next, then connecting ones"
    assert len({g for g in groups if g is not None}) == 2
    # Within a tangle, by smallest member cdist.
    for tangle in {g for g in groups}:
        cdists = [e.min_unstable_cdist for e in table if e.tangle == tangle]
        assert cdists == sorted(cdists)
    # Active letters run a, b, c, ... tangle by tangle.
    active = [e.letter for e in table if e.letter is not None]
    assert active == sorted(active, key=lambda letter: (len(letter), letter))
    for refined in (False, True):
        palette = plotting.class_colors(dynamics, refined=refined)
        assert len(set(palette.values())) == len(palette), "no colour used twice"
    palette = plotting.class_colors(dynamics, refined=False)
    by_tangle = {}
    for cd in dynamics.classes.values():
        by_tangle.setdefault(cd.entry.tangle, set()).add(palette[cd.letter])
    families = [set(family) for family in plotting.TANGLE_COLOR_FAMILIES]
    for tangle, colours in by_tangle.items():
        if tangle is not None:
            assert any(colours <= family for family in families)


@pytest.mark.slow
def test_nested_outer_arcs_go_around_the_inner_circle(nested_built):
    """No inward arc of the outer circle enters the inner one; a bridge whose
    lobe encloses the inner tangle goes around its far side."""
    session = nested_built.session
    outer, inner = nested_built.fixed_points
    _interior, parents, _arcs = session.cartoon_zones()
    minimal = session.minimal_trellis()
    enclosures = session.cartoon_enclosures(minimal.kept_bridge_ids, parents)
    assert enclosures, "the outer anchor lobes enclose the period-3 tangle"
    assert all(ids == frozenset({id(inner)}) for ids in enclosures.values())
    fig, ax = plt.subplots()
    try:
        layout = session.plot_dual_graph_cartoon(ax=ax, shape="circle")
        cx, cy, _r = layout.circles[id(outer)]
        keepout = layout.keepout[id(outer)]
        assert keepout > layout.circles[id(inner)][2]
        routed = 0
        for patch in ax.patches:
            gid = patch.get_gid() or ""
            if not gid.startswith("bridge:"):
                continue
            first, second = (int(v) for v in gid.split(":")[1].split("-"))
            registry = session.workbench.intersection_registry
            if registry[first].manifold_b_key[0] is not outer:
                continue
            vertices = np.asarray(patch.get_path().vertices)
            radii = np.hypot(vertices[:, 0] - cx, vertices[:, 1] - cy)
            if len(vertices) > 4:  # a routed (polyline) arc
                routed += 1
                assert radii.min() >= keepout - 1e-9
                # The long way round sweeps more than half the circle.
                angles = np.unwrap(np.arctan2(vertices[:, 1] - cy, vertices[:, 0] - cx))
                swept = abs(angles[-1] - angles[0])
                assert (swept > np.pi) == ((first, second) in enclosures)
        assert routed >= len(enclosures)
    finally:
        plt.close(fig)


# --------------------------------------------------------------------------- #
# Nested blasts: each zone blasts only its own tangle
# --------------------------------------------------------------------------- #
def _tangle_words(build, fixed_point) -> dict:
    """Every class of one tangle mapped to its word, letters spelled out.

    Letters differ between sessions, so a class is named by its two homotopy
    elements (orbit code and short name, no fixed-point letter) and each token
    by its class's name and direction.
    """
    dynamics = build.session.symbolic_dynamics()

    def spelled(bridge_class):
        ends = (bridge_class.source, bridge_class.target)
        names = (dynamics.naming.homotopy_name(ref) for ref in ends)
        return tuple((name.orbit_code, name.short_text) for name in names)

    return {
        spelled(cd.bridge_class): [(spelled(s.bridge_class), s.direction) for s in cd.symbols]
        for cd in dynamics.classes.values()
        if cd.bridge_class.source.branch_key[0] is fixed_point
    }


@pytest.mark.slow
@pytest.mark.regression
def test_nested_inner_words_are_the_period3_words():
    """The outer blasts leave the inner tangle exactly as it is alone."""
    alone = build_period3(blasts=4)
    nested = build_nested(outer_blasts=2, inner_blasts=4)
    _outer, inner = nested.fixed_points
    words = _tangle_words(alone, alone.fixed_points[0])
    assert words and _tangle_words(nested, inner) == words


@pytest.mark.slow
@pytest.mark.regression
def test_nested_blast_order_does_not_matter(caplog):
    with caplog.at_level(logging.WARNING, logger="tanglepack.loom.Blast"):
        inner_first = build_nested(outer_blasts=4, inner_blasts=4)
    assert not [r for r in caplog.records if "skipping bridge" in r.message]

    outer_first = build_nested(outer_blasts=4)
    session = outer_first.session
    inner_zone = min(session.resonance_zones.values(), key=lambda zone: zone.area)
    for _ in range(4):
        session.blast_zone(
            inner_zone, 1, fixed_point=[inner_zone.fixed_point], min_separation=1e-4
        )
        session.classify_strong_pips()
        for fp, pip in zip(outer_first.fixed_points, outer_first.pips):
            session.set_strong_pip(fp, pip)
    session.compute_pseudoneighbors()
    session.punch_holes()
    session.partition_stable_manifold()

    for one, other in zip(inner_first.fixed_points, outer_first.fixed_points):
        words = _tangle_words(inner_first, one)
        assert words and _tangle_words(outer_first, other) == words
        assert len(inner_first.session.trellis(one).bridges) == len(
            session.trellis(other).bridges
        )


_P3_WORDS = {"a": "b", "b": "c", "c": "a u^-1 w^-1"}


def _active_words(build) -> dict:
    dynamics = build.session.symbolic_dynamics()
    return {cd.letter: cd.word for cd in dynamics.classes.values() if cd.kind == "active"}


@pytest.mark.slow
def test_p3_words_survive_four_blasts():
    """Blasting the closed period-3 zone changes nothing topologically."""
    assert _active_words(build_period3(blasts=4)) == _P3_WORDS


@pytest.mark.slow
def test_nested_default_words_are_pinned(nested_built):
    """The nested words at the defaults, as at 2315204 (2026-10-02)."""
    outer_words = {"d": "d uu^-1 e^-1", "e": "f", "f": "d uu^-1 d^-1"}
    assert _active_words(nested_built) == {**_P3_WORDS, **outer_words}


@pytest.mark.slow
@pytest.mark.regression
def test_nested_outer_blasts_leave_the_inner_bridges_alone(p3_built, nested_built):
    """The outer zone's blasts no longer iterate the period-3 bridges (2026-10-02)."""
    _outer, inner = nested_built.fixed_points
    (fp3,) = p3_built.fixed_points
    inner_bridges = nested_built.session.trellis(inner).bridges
    assert len(inner_bridges) == len(p3_built.session.trellis(fp3).bridges)

"""The trellis cartoon: what it draws is the trellis's combinatorics, in order.

Drawing properties over the law cases (each built fresh, function scope):
one representative per class, ranks in stable-cdist order, the route kind
read off the rows, the drawn selection, the iterates as image chains, and the
planarity law on the single-branch cases (two unstable curves never cross, so
neither do their cartoon curves). No law case has a hairpin (a one-branch
bridge with two rows), so the fold is held to the same law on synthetic nested
ends, and a flipped row makes one. No case fact is pinned.
"""

from __future__ import annotations

import itertools

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from cases import BUILDERS, LAW_CASES  # noqa: E402
from tanglepack.topology import TrellisCartoon, class_representatives, plot_trellis_cartoon  # noqa: E402
from tanglepack.topology.MinimalTrellis import image_chain  # noqa: E402
from tanglepack.topology.StablePartition import rows_of_bridge  # noqa: E402


@pytest.fixture(params=LAW_CASES)
def topology(request):
    """A fresh build of one law case, its topology recomputed."""
    return BUILDERS[request.param]().session.topology()


@pytest.fixture
def ax():
    """A fresh axes, closed after the test."""
    fig, axes = plt.subplots()
    yield axes
    plt.close(fig)


def _draw(topology, ax, **kwargs):
    """The cartoon of a topology, no zones (the layout's defaults)."""
    return plot_trellis_cartoon(
        topology.trellis, topology.classes, topology.dynamics, minimal=topology.minimal,
        ax=ax, **kwargs,
    )


def _names(topology) -> dict:
    """``{BridgeClass: letter}`` as the cartoon names classes (inert letters included)."""
    return {entry.bridge_class: topology.dynamics.classes[entry.bridge_class].letter
            for entry in topology.classes}


def test_one_representative_per_class_from_its_members(topology) -> None:
    """Every spanned class gets one of its own members, a hole bridge when it has one."""
    minimal = topology.minimal
    representatives = class_representatives(topology.classes, minimal, topology.dynamics)
    names = _names(topology)
    spanned = [entry for entry in topology.classes if entry.spanned]
    assert set(representatives) == {names[entry.bridge_class] for entry in spanned}
    holes, images = set(minimal.hole_bridge_ids), set(minimal.image_bridge_ids)
    for entry in spanned:
        chosen = representatives[names[entry.bridge_class]]
        assert chosen in entry.bridge_ids
        if holes & set(entry.bridge_ids):
            assert chosen in holes
        elif images & set(entry.bridge_ids):
            assert chosen in images


def test_ranks_follow_the_stable_cdist_from_the_anchor(topology, ax) -> None:
    """On every branch the anchor is rank 0 and the ranks increase strictly in stable cdist."""
    trellis = topology.trellis
    for compact in (True, False):
        layout = _draw(topology, ax, show="all", compact=compact, holes=True)
        assert set(layout.order) == {branch.key for branch in trellis.stable_branches}
        for branch_key, keys in layout.order.items():
            anchor, *rest = keys
            assert anchor == (anchor[0], branch_key) and trellis.intersection(anchor[0]).is_anchor
            cdists = [trellis.intersection(key).cdist_key("stable") for key in rest]
            assert all(not isinstance(key, tuple) for key in rest)
            assert cdists == sorted(cdists) and len(set(cdists)) == len(cdists)
            assert cdists[:1] != [(0.0, 0.0)]


def test_the_route_kind_is_read_off_the_rows(topology, ax) -> None:
    """``u`` iff one branch and one row, ``hairpin`` iff one branch and two rows, else ``cross``."""
    trellis = topology.trellis
    layout = _draw(topology, ax, show="all")
    assert layout.bridges
    for bridge_id, kind in layout.bridges.items():
        unstable_key = trellis.unstable_key_of(bridge_id)
        branches = {trellis.stable_key_at(i, unstable_key) for i in bridge_id}
        rows = set(rows_of_bridge(trellis, bridge_id))
        expected = "cross" if len(branches) == 2 else "u" if len(rows) == 1 else "hairpin"
        assert kind == expected, bridge_id


def test_show_selects_every_bridge_or_the_named_classes(topology, ax) -> None:
    """``all`` draws every trellis bridge; ``letters`` restricts to those classes' members."""
    trellis = topology.trellis
    every = _draw(topology, ax, show="all")
    assert set(every.bridges) == {bridge.id for bridge in trellis.bridges if bridge.id is not None}
    representatives = _draw(topology, ax)
    assert set(representatives.bridges) == set(representatives.representatives.values())
    names = _names(topology)
    chosen = sorted({names[entry.bridge_class] for entry in topology.classes if entry.spanned})[:2]
    members = {
        bridge_id
        for entry in topology.classes
        if names[entry.bridge_class] in chosen
        for bridge_id in entry.bridge_ids
    }
    assert set(_draw(topology, ax, show="all", letters=chosen).bridges) == members
    picked = set(_draw(topology, ax, letters=chosen).bridges)
    assert picked <= members and len(picked) == len(chosen)
    one = next(iter(members))
    assert set(_draw(topology, ax, bridge_ids=[one]).bridges) == {one}


def test_the_drawn_iterates_are_image_chains(topology, ax) -> None:
    """With ``iterates=1`` every drawn image pair is ``_image_chain`` of its bridge."""
    trellis = topology.trellis
    layout = _draw(topology, ax, show="all", iterates=1)
    for bridge_id in layout.bridges:
        chain = image_chain(trellis, bridge_id)
        assert layout.iterates[bridge_id] == ([chain] if chain else [])
    deeper = _draw(topology, ax, iterates=2)
    for bridge_id, steps in deeper.iterates.items():
        assert len(steps) <= 2
        for before, after in itertools.pairwise(steps):
            assert after == [pair for source in before for pair in (image_chain(trellis, source) or [])]


def _crossings(a: np.ndarray, b: np.ndarray) -> int:
    """How many segments of polyline ``a`` properly cross segments of polyline ``b``."""
    p, r = a[:-1, None, :], (a[1:] - a[:-1])[:, None, :]
    q, s = b[None, :-1, :], (b[1:] - b[:-1])[None, :, :]
    cross = r[..., 0] * s[..., 1] - r[..., 1] * s[..., 0]
    qp = q - p
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (qp[..., 0] * s[..., 1] - qp[..., 1] * s[..., 0]) / cross
        u = (qp[..., 0] * r[..., 1] - qp[..., 1] * r[..., 0]) / cross
    return int(np.count_nonzero((cross != 0) & (t > 0) & (t < 1) & (u > 0) & (u < 1)))


@pytest.mark.parametrize("name", ["k10", "k28_one_blast", "k28_two_blasts"])
def test_single_branch_cartoon_is_planar(name: str, ax) -> None:
    """No two drawn bridges of a one-branch tangle cross (two unstable curves never do)."""
    topology = BUILDERS[name]().session.topology()
    layout = _draw(topology, ax, show="all")
    assert len(layout.order) == 1 and len(layout.paths) > 1
    for (id_a, a), (id_b, b) in itertools.combinations(layout.paths.items(), 2):
        # Curves sharing a crossing leave it side by side: the shared ends are
        # touching, not crossing, so their first samples are left out.
        trim_a = slice(3 if id_a[0] in id_b else 0, -3 if id_a[1] in id_b else None)
        trim_b = slice(3 if id_b[0] in id_a else 0, -3 if id_b[1] in id_a else None)
        assert _crossings(a[trim_a], b[trim_b]) == 0, (id_a, id_b)


def _random_nesting(count: int, rng: np.random.Generator) -> list[tuple[int, int]]:
    """A random non-crossing perfect matching of ``0..count-1`` (a random Dyck word)."""
    pairs, stack, opened = [], [], 0
    for i in range(count):
        if stack and (opened == count // 2 or rng.random() < 0.5):
            pairs.append((stack.pop(), i))
        else:
            stack.append(i)
            opened += 1
    return pairs


def _past_the_pip_end(layout, curve: np.ndarray) -> np.ndarray:
    """The angular progress of a curve along its branch's arc: 0 at the anchor, 1 at the pip end."""
    (branch_key,) = layout.order
    fp_id, alpha, sweep = layout.zone.arcs[branch_key]
    cx, cy, _radius = layout.zone.circles[fp_id]
    return np.unwrap(np.arctan2(curve[:, 1] - cy, curve[:, 0] - cx) - alpha) / sweep


@pytest.mark.parametrize("seed", range(6))
def test_folded_hairpins_and_u_arcs_of_one_branch_never_cross(seed: int, ax) -> None:
    """Arcs nesting on the unfolded line give curves that never cross, hairpins included.

    No law case has a one-branch bridge with two rows, so the fold is held to
    the planarity law on synthetic ends: every rank on both rows of a branch,
    matched at random without crossing on the unfolded line (``w = 0`` the pip
    end). A matched pair straddling ``w = 0`` is a hairpin; it must run past
    the pip end and never wrap the anchor.
    """
    layout = _draw(BUILDERS["k10"]().session.topology(), ax)
    (branch_key,) = layout.order
    chart = TrellisCartoon._Chart(layout.zone, branch_key, 9)
    rows = sorted(chart.sigma)
    ends = sorted((chart.unfolded(rank, row), rank, row) for rank in range(1, 9) for row in rows)
    # The outermost pair joins the two rows' last ranks: one hairpin at least.
    inner = _random_nesting(len(ends) - 2, np.random.default_rng(seed))
    pairs = [(0, len(ends) - 1)] + [(i + 1, j + 1) for i, j in inner]
    curves = {(i, j): chart.arc(ends[i][0], ends[j][0]) for i, j in pairs}
    hairpins = [(i, j) for i, j in pairs if ends[i][2] != ends[j][2]]
    assert hairpins

    for i, j in hairpins:
        progress = _past_the_pip_end(layout, curves[i, j])
        assert progress.max() > 1.0 and progress.min() >= 0.0, (ends[i], ends[j])
    for (a, curve_a), (b, curve_b) in itertools.combinations(curves.items(), 2):
        # Two ends at one rank on the two rows are one crossing: touching there, not crossing.
        assert _crossings(curve_a, curve_b) == 0, (a, b)


def test_a_one_branch_bridge_with_two_rows_is_a_hairpin(monkeypatch, ax) -> None:
    """Flipping one end's row turns a bridge's U-arc into a hairpin round its pip end."""
    topology = BUILDERS["k10"]().session.topology()
    layout = _draw(topology, ax)
    bridge_id, kind = next(iter(layout.bridges.items()))
    assert kind == "u"
    flip = {"left": "right", "right": "left"}
    signed_row = TrellisCartoon.row_of_end

    def row_of_end(trellis, pair, endpoint):
        row = signed_row(trellis, pair, endpoint)
        return flip[row] if pair == bridge_id and endpoint == "second" else row

    monkeypatch.setattr(TrellisCartoon, "row_of_end", row_of_end)
    layout = _draw(topology, ax, bridge_ids=[bridge_id])
    assert layout.bridges == {bridge_id: "hairpin"}
    progress = _past_the_pip_end(layout, layout.paths[bridge_id])
    assert progress.max() > 1.0 and progress.min() >= 0.0


def test_representatives_need_the_minimal_trellis(ax) -> None:
    """Without ``minimal`` only an explicit bridge list or ``show="all"`` can be drawn."""
    topology = BUILDERS["k10"]().session.topology()
    trellis, classes = topology.trellis, topology.classes
    with pytest.raises(ValueError, match="minimal"):
        plot_trellis_cartoon(trellis, classes, ax=ax)
    with pytest.raises(ValueError, match="show"):
        plot_trellis_cartoon(trellis, classes, ax=ax, show="some", minimal=topology.minimal)
    with pytest.raises(ValueError, match="no class"):
        plot_trellis_cartoon(trellis, classes, ax=ax, show="all", letters=["no such class"])
    layout = plot_trellis_cartoon(trellis, classes, ax=ax, show="all")
    assert set(layout.bridges) == {b.id for b in trellis.bridges if b.id is not None}
    assert layout.representatives == {}

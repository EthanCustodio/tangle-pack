"""
The three rules that make a cone point's tangle readable (Mitchell 2009, Fig. 8):

* an anchor with one stable branch sits on the row its first arc arrives on
  (``Trellis.settle_one_branch_anchors``);
* a pair of elements whose member lobes hold different markers is several
  classes, one per lobe (``LobeSplit``, ``BridgeClass.lobe``);
* of tied walks, the one whose end chords keep the order of the arcs leaving
  that element survives (``SymbolicDynamics.fan_order``).

Each is tested on hand-built data shaped like the hydrogen trellis; the rules
themselves are the library's.
"""

from __future__ import annotations

from types import SimpleNamespace

import pytest

from helpers.fakes import bare_fixed_point, ref, stable_key
import importlib

from tanglepack.loom import persist
from tanglepack.topology.BridgeClass import BridgeClass, oriented_class
from tanglepack.topology.DualGraph import face_side
from tanglepack.topology.Trellis import Trellis

# The package re-exports classes under their modules' names: import the modules.
lobe_split = importlib.import_module("tanglepack.topology.LobeSplit")
symbolic = importlib.import_module("tanglepack.topology.SymbolicDynamics")

HOLE_A = ("hole", (254, 255), -2)
HOLE_B = ("hole", (254, 255), -1)
HOLE_C = ("hole", (31, 32), 0)
ORBIT = ("orbit", 1, 0)


@pytest.fixture
def inner():
    """One stable branch of one fixed point, and three of its right-row elements."""
    fp = bare_fixed_point()
    key = stable_key(fp)
    return fp, key, [ref(key, "right", i) for i in range(3)]


# --------------------------------------------------------------------------- #
# BridgeClass.lobe
# --------------------------------------------------------------------------- #
def test_a_lobe_tells_two_classes_of_one_pair_apart(inner):
    fp, _key, (r0, _r1, r2) = inner
    pair, direction = oriented_class(r0, r2, [fp])
    assert direction == 1 and pair.lobe == frozenset()
    g, b = pair.with_lobe({HOLE_A, HOLE_C}), pair.with_lobe({HOLE_B})
    assert g != b and g != pair and len({g, b, pair}) == 3
    assert (g.source, g.target) == (b.source, b.target) == (pair.source, pair.target)
    assert pair.with_lobe(None) is pair
    assert "hole 254.255^-2" in g.label and "hole 31.32^0" in g.label
    assert g.label != b.label


def test_a_lobe_survives_a_save():
    item = ("lobe", frozenset({HOLE_A, ORBIT}))
    saved = persist._encode_identity_item(item, key=None)
    assert persist._decode_identity_item(saved, key=None) == item


# --------------------------------------------------------------------------- #
# LobeSplit
# --------------------------------------------------------------------------- #
def _table(pair, bridge_ids):
    members = [SimpleNamespace(bridge_id=bid, is_loop=False) for bid in bridge_ids]
    return [SimpleNamespace(bridge_class=pair, members=members)]


def test_members_whose_lobes_disagree_split_the_pair(inner, monkeypatch):
    """The paper's g, the bridges going round the island like it, and b."""
    fp, _key, (r0, _r1, r2) = inner
    pair, _ = oriented_class(r0, r2, [fp])
    lobes = {
        (1, 236): {HOLE_A: True, HOLE_C: True, HOLE_B: False},  # g: y lies on it
        (321, 322): {HOLE_A: True, HOLE_C: True, HOLE_B: False, ORBIT: True},
        (237, 246): {HOLE_A: False, HOLE_C: False, HOLE_B: True, ORBIT: False},
        (439, 440): {HOLE_A: False, HOLE_C: False, HOLE_B: True, ORBIT: False},
    }
    monkeypatch.setattr(lobe_split, "bridge_lobe", lambda dual, bid: lobes[tuple(bid)])
    splits, tags = lobe_split.split_by_lobes(None, _table(pair, list(lobes)))
    assert list(splits) == [pair]
    g, b = frozenset({HOLE_A, HOLE_C, ORBIT}), frozenset({HOLE_B})
    assert tags == {(1, 236): g, (321, 322): g, (237, 246): b, (439, 440): b}
    split = splits[pair]
    # A walk piece is named by its own lobe; one agreeing with both is not named.
    assert split.lobe_for({HOLE_A: True, HOLE_C: True, ORBIT: True}) == g
    assert split.lobe_for({HOLE_B: True}) == b
    assert split.lobe_for({}) is None
    assert split.lobe_for(None) is None


def test_members_whose_lobes_agree_stay_one_class(inner, monkeypatch):
    """A marker on one lobe's curve (None, dropped) decides nothing."""
    fp, _key, (r0, _r1, r2) = inner
    pair, _ = oriented_class(r0, r2, [fp])
    lobes = {(1, 2): {HOLE_C: True}, (3, 4): {HOLE_C: True, ORBIT: False}, (5, 6): {}}
    monkeypatch.setattr(lobe_split, "bridge_lobe", lambda dual, bid: lobes[tuple(bid)])
    assert lobe_split.split_by_lobes(None, _table(pair, list(lobes))) == ({}, {})


def test_a_pair_across_two_branches_is_never_split(monkeypatch):
    fp = bare_fixed_point()
    a, b = ref(stable_key(fp, branch=0), "right", 0), ref(stable_key(fp, branch=1), "right", 0)
    pair, _ = oriented_class(a, b, [fp])
    monkeypatch.setattr(lobe_split, "bridge_lobe", lambda dual, bid: pytest.fail("read a lobe"))
    assert lobe_split.split_by_lobes(None, _table(pair, [(1, 2), (3, 4)])) == ({}, {})


# --------------------------------------------------------------------------- #
# fan_order
# --------------------------------------------------------------------------- #
def _arc(edge_key, kind, tail, head, reverse=False):
    return SimpleNamespace(edge_key=edge_key, kind=kind, tail_id=tail, head_id=head, reverse=reverse)


def _fan_fixture(inner, monkeypatch, d2_feet):
    """
    The start face of M(d2) in Fig. 8: going round from the start edge's end at
    the fixed point, g, the inner arc near r0, d1, the wall above a2's hole, a2,
    the wall below it, d3, and the start edge (B from 0.141 down to y).
    """
    fp, key, (r0, _r1, r2) = inner
    cdist = {1: 0.0, 262: 0.1408, 236: 1.019, 87: 0.915, 89: 1.395, 35: 1.186, 44: 1.036, 255: 0.988}
    arcs = [
        _arc("g", "unstable", 1, 236),
        _arc("B near r0", "stable", 236, 87),
        _arc("d1", "unstable", 87, 89),
        _arc("A above a2", "stable", 89, 35),
        _arc("a2", "unstable", 35, 44),
        _arc("A below a2", "stable", 44, 255),
        _arc("d3", "unstable", 255, 262),
        _arc("start", "stable", 262, 1),
    ]
    face = SimpleNamespace(faces=[SimpleNamespace(arcs=arcs)])
    end = {i: (arc.edge_key, face_side(arc)) for i, arc in enumerate(arcs)}
    trellis = SimpleNamespace(
        intersection=lambda i: SimpleNamespace(stable_cdist=cdist[i]),
        fixed_points=[fp],
    )
    dual = SimpleNamespace(trellis=trellis)
    other = ref(key, "right", 9)

    def walk(ring_position):
        return SimpleNamespace(multiplicity=1, far=end[ring_position])

    def record(name, feet, ring_position, search=None):
        cls = BridgeClass(r0, other if name != "d2" else r2)
        cls = cls.with_lobe({name}) if name != "d2" else cls
        rec = SimpleNamespace(
            bridge_class=cls, feet=feet, letter=name,
            search=search or SimpleNamespace(status="resolved", walk=walk(ring_position)),
        )
        return rec

    neighbours = [record("g", (0.0, 0.0, False), 1), record("d3", (0.1408, 0.2598, False), 3),
                  record("b", (0.3789, 0.3789, False), 5)]
    tied = [walk(1), walk(5)]
    d2 = record("d2", d2_feet, None, SimpleNamespace(status="ambiguous", walks=tied, walk=tied[0]))
    records = {rec.bridge_class: rec for rec in neighbours + [d2]}
    monkeypatch.setattr(
        symbolic, "_end_chord",
        lambda w, which: (face, end[7], w.far) if which == "start" else None,
    )
    monkeypatch.setattr(symbolic, "_class_feet", lambda tr, cd, el, recs, naming: cd.feet)
    return dual, d2, records, tied


def test_the_tied_walk_crossing_a_neighbour_is_dropped(inner, monkeypatch):
    """d2 leaves nearer the fixed point than d3, so its image cannot reach the
    far side of a2, past the image of d3: of g^-1 x d5 and d3 v a1, the first."""
    dual, d2, records, tied = _fan_fixture(inner, monkeypatch, (0.0, 0.1408, True))
    assert symbolic.fan_order(dual, d2, records, naming=None) is tied[0]


def test_feet_with_no_order_decide_nothing(inner, monkeypatch):
    """Feet interleaving every neighbour's give no order: the tie stands."""
    dual, d2, records, tied = _fan_fixture(inner, monkeypatch, (0.0, 0.5, False))
    assert symbolic.fan_order(dual, d2, records, naming=None) is None


def test_feet_beyond_every_neighbour_keep_the_other_walk(inner, monkeypatch):
    """The order cuts both ways: feet past b's admit only the walk ending past b's."""
    dual, d2, records, tied = _fan_fixture(inner, monkeypatch, (0.5, 0.6, False))
    assert symbolic.fan_order(dual, d2, records, naming=None) is tied[1]


def test_before_reads_strict_bounds():
    assert symbolic._before((0.0, 0.1, False), (0.2, 0.3, False))
    assert not symbolic._before((0.0, 0.2, False), (0.2, 0.3, False))
    assert symbolic._before((0.0, 0.2, True), (0.2, 0.3, False))
    assert not symbolic._before((0.2, 0.3, False), (0.0, 0.1, False))


# --------------------------------------------------------------------------- #
# Trellis.settle_one_branch_anchors
# --------------------------------------------------------------------------- #
def _anchor_trellis(anchor_sign, arrival_sign, *, second_branch=False, arrival_on_own=True):
    fp = bare_fixed_point()
    s_key, u_key = stable_key(fp), (fp, "unstable", 0, 0)
    other_s = (fp, "stable", 0, 1)
    foreign = stable_key(bare_fixed_point())
    anchor = SimpleNamespace(is_anchor=True, crossing_sign=anchor_sign,
                             manifold_a_key=u_key, manifold_b_key=s_key)
    arrival = SimpleNamespace(is_anchor=False, crossing_sign=arrival_sign)
    crossings = {1: anchor, 236: arrival}
    fake = SimpleNamespace(
        registry=list(crossings.items()),
        branches={s_key: object(), **({other_s: object()} if second_branch else {})},
        bridges=[SimpleNamespace(id=(1, 236))],
        unstable_key_of=lambda bid: u_key,
        stable_key_at=lambda cid, uk: s_key if (cid == 1 or arrival_on_own) else foreign,
        intersection=lambda cid: crossings[cid],
    )
    fake.__contains__ = lambda cid: cid in crossings
    return fake, anchor


class _Contains(SimpleNamespace):
    def __contains__(self, item):
        return self.__dict__["__contains__"](item)


def _settle(fake):
    return Trellis.settle_one_branch_anchors(_Contains(**vars(fake)))


def test_a_one_branch_anchor_takes_its_first_arcs_row():
    fake, anchor = _anchor_trellis(+1, +1)
    assert _settle(fake) == [1]
    assert anchor.crossing_sign == -1
    assert _settle(fake) == []  # once: the sign now agrees


@pytest.mark.parametrize(
    "case",
    [
        dict(anchor_sign=+1, arrival_sign=-1),  # already on the arc's row
        dict(anchor_sign=+1, arrival_sign=+1, second_branch=True),  # an ordinary saddle
        dict(anchor_sign=+1, arrival_sign=+1, arrival_on_own=False),  # arrives elsewhere
    ],
)
def test_every_other_anchor_is_left_alone(case):
    fake, anchor = _anchor_trellis(
        case["anchor_sign"], case["arrival_sign"],
        second_branch=case.get("second_branch", False),
        arrival_on_own=case.get("arrival_on_own", True),
    )
    assert _settle(fake) == []
    assert anchor.crossing_sign == case["anchor_sign"]

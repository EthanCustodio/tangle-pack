"""
Tests for ``tanglepack.topology.ElementNaming``: the ``L_i`` / ``R_i^j`` names
of the homotopy and iterated partition elements.

The synthetic tests build partition families by hand (a fake fixed point is
all an ``ElementRef`` asks for); the fixture tests pin the k=10 combinatorics
of the plan BY NAME TEXT, never by registry id.
"""

from __future__ import annotations

import pytest

from minimal_helpers import build_pieces
from tanglepack.topology.ElementNaming import ElementName, ElementNaming
from tanglepack.topology.PartitionFamily import (
    HomotopyPartition,
    IteratedHomotopyPartition,
)
from tanglepack.topology.TopologyResults import (
    ElementRef,
    PartitionInterval,
    StablePartitionResult,
)


# --------------------------------------------------------------------------- #
# Synthetic builders
# --------------------------------------------------------------------------- #
class _FakeFixedPoint:
    """The only thing a name asks of a fixed point is its period."""

    def __init__(self, period: int = 1) -> None:
        self.period = period


def _key(fp, orbit: int = 0, branch: int = 0):
    return (fp, "stable", orbit, branch)


def _result(branch_key, side, bounds, parents=None) -> StablePartitionResult:
    """A result with consecutive elements over ``bounds`` (cdists), stamped
    anchor outward; ``parents`` (per element) sets ``parent_element_id``."""
    intervals = []
    for element_id, (lo, hi) in enumerate(zip(bounds, bounds[1:])):
        interval = PartitionInterval(
            lo_id=None if element_id == 0 else element_id,
            hi_id=element_id + 1,
            lo_cdist=float(lo),
            hi_cdist=float(hi),
            closed_lo=True,
            closed_hi=False,
            element_id=element_id,
            branch_key=branch_key,
            side=side,
        )
        if parents is not None:
            interval.parent_element_id = parents[element_id]
        intervals.append(interval)
    return StablePartitionResult(branch_key=branch_key, side=side, intervals=intervals)


def _one_branch_families():
    """One branch, both sides: homotopy L has 3 elements, R has 3; the
    iterated family splits L_1 in two and R_3 in three."""
    fp = _FakeFixedPoint()
    key = _key(fp)
    homotopy = HomotopyPartition.from_results(
        [
            _result(key, "left", [0, 1, 2, 3]),
            _result(key, "right", [0, 1, 2, 3]),
        ]
    )
    iterated = IteratedHomotopyPartition(
        [
            _result(key, "left", [0, 0.5, 1, 2, 3], parents=[0, 0, 1, 2]),
            _result(key, "right", [0, 1, 2, 2.3, 2.6, 3], parents=[0, 1, 2, 2, 2]),
        ],
        homotopy=homotopy,
    )
    return fp, key, homotopy, iterated


# --------------------------------------------------------------------------- #
# ElementName
# --------------------------------------------------------------------------- #
def test_element_name_text_and_mathtext():
    key = _key(_FakeFixedPoint(3), orbit=1, branch=0)
    split = ElementName(key, "right", 1, 2, 3, False)
    assert split.side_letter == "R"
    assert split.is_split
    assert split.text == "R_1^2"
    assert str(split) == "R_1^2"
    assert split.mathtext == "$R_{1}^{2}$"
    assert split.tag == "p3@1.0"

    plain = ElementName(key, "left", 3, 1, 1, False)
    assert plain.side_letter == "L"
    assert not plain.is_split
    assert plain.text == "L_3"
    assert plain.mathtext == "$L_{3}$"

    tagged = ElementName(key, "right", 2, 1, 2, True)
    assert tagged.text == "p3@1.0:R_2^1"
    assert tagged.mathtext == "p3@1.0:$R_{2}^{1}$"
    tagged_plain = ElementName(key, "right", 2, 1, 1, True)
    assert tagged_plain.text == "p3@1.0:R_2"


def test_element_name_is_hashable_and_compares_by_value():
    key = _key(_FakeFixedPoint())
    a = ElementName(key, "right", 1, 2, 3, False)
    b = ElementName(key, "right", 1, 2, 3, False)
    assert a == b and hash(a) == hash(b)
    assert len({a, b, ElementName(key, "right", 1, 3, 3, False)}) == 2
    assert repr(a) == "ElementName('R_1^2')"


# --------------------------------------------------------------------------- #
# ElementNaming, synthetic
# --------------------------------------------------------------------------- #
def test_homotopy_names_count_anchor_outward_per_side():
    _fp, key, homotopy, iterated = _one_branch_families()
    naming = ElementNaming(homotopy, iterated)
    assert not naming.tagged
    assert [n.text for n in naming.homotopy_names] == [
        "L_1", "L_2", "L_3", "R_1", "R_2", "R_3",
    ]
    for side in ("left", "right"):
        for element_id in range(3):
            name = naming.homotopy_name(ElementRef(key, side, element_id))
            assert name.subscript == element_id + 1
            assert name.superscript == 1 and name.siblings == 1


def test_superscripts_follow_parent_groups_anchor_outward():
    _fp, key, homotopy, iterated = _one_branch_families()
    naming = ElementNaming(homotopy, iterated)
    assert [n.text for n in naming.names] == [
        "L_1^1", "L_1^2", "L_2", "L_3",
        "R_1", "R_2", "R_3^1", "R_3^2", "R_3^3",
    ]
    assert len(naming) == 9
    r33 = naming.name(ElementRef(key, "right", 4))
    assert (r33.subscript, r33.superscript, r33.siblings) == (3, 3, 3)
    assert r33.is_split


def test_unsplit_elements_print_plain():
    _fp, key, homotopy, iterated = _one_branch_families()
    naming = ElementNaming(homotopy, iterated)
    r1 = naming.name(ElementRef(key, "right", 0))
    assert not r1.is_split
    assert r1.text == "R_1" and r1.mathtext == "$R_{1}$"
    assert r1.superscript == 1 and r1.siblings == 1


def test_ref_of_name_round_trips_and_lookup_by_text():
    _fp, key, homotopy, iterated = _one_branch_families()
    naming = ElementNaming(homotopy, iterated)
    for ref in naming.refs:
        name = naming.name(ref)
        assert naming.ref_of(name) == ref
        assert naming.lookup(name.text) == ref
    assert naming.lookup("R_3^2") == ElementRef(key, "right", 3)
    # A name assembled by hand with the same fields resolves too.
    assert naming.ref_of(ElementName(key, "left", 1, 2, 2, False)) == ElementRef(key, "left", 1)
    with pytest.raises(KeyError, match="R_9"):
        naming.lookup("R_9")
    with pytest.raises(KeyError):
        naming.ref_of(ElementName(key, "left", 1, 3, 3, False))


def test_parent_of_and_children_of_round_trip():
    _fp, key, homotopy, iterated = _one_branch_families()
    naming = ElementNaming(homotopy, iterated)
    for parent in naming.homotopy_refs:
        children = naming.children_of(parent)
        assert children == sorted(children, key=lambda c: c.element_id)
        for position, child in enumerate(children, start=1):
            assert naming.parent_of(child) == parent
            assert naming.name(child).superscript == position
            assert naming.name(child).siblings == len(children)
    # Every iterated element has exactly one parent, and the groups tile the family.
    assert sorted(
        (c for p in naming.homotopy_refs for c in naming.children_of(p)),
        key=lambda r: (r.side, r.element_id),
    ) == sorted(naming.refs, key=lambda r: (r.side, r.element_id))
    assert naming.parent_of(ElementRef(key, "right", 4)) == ElementRef(key, "right", 2)
    assert naming.children_of(ElementRef(key, "left", 0)) == [
        ElementRef(key, "left", 0), ElementRef(key, "left", 1),
    ]


def test_the_two_families_are_never_mixed():
    """Refs collide by value: iterated L#3 is L_3, homotopy L#2 is L_3 too."""
    _fp, key, homotopy, iterated = _one_branch_families()
    naming = ElementNaming(homotopy, iterated)
    assert naming.homotopy_name(ElementRef(key, "left", 2)).text == "L_3"
    assert naming.name(ElementRef(key, "left", 2)).text == "L_2"
    assert naming.name(ElementRef(key, "left", 3)).text == "L_3"
    with pytest.raises(KeyError):
        naming.homotopy_name(ElementRef(key, "left", 3))  # only the iterated family has #3
    with pytest.raises(KeyError):
        naming.name(ElementRef(key, "right", 5))
    with pytest.raises(KeyError):
        naming.parent_of(ElementRef(key, "right", 5))
    with pytest.raises(KeyError):
        naming.children_of(ElementRef(key, "right", 3))
    assert ElementRef(key, "right", 5) not in naming
    assert ElementRef(key, "right", 4) in naming


def test_homotopy_only_naming_is_its_own_parent():
    _fp, key, homotopy, _iterated = _one_branch_families()
    naming = ElementNaming(homotopy)
    assert not naming.has_iterated
    assert [n.text for n in naming.names] == ["L_1", "L_2", "L_3", "R_1", "R_2", "R_3"]
    for ref in naming.refs:
        assert naming.name(ref) == naming.homotopy_name(ref)
        assert not naming.name(ref).is_split
        assert naming.parent_of(ref) == ref
        assert naming.children_of(ref) == [ref]
        assert naming.ref_of(naming.name(ref)) == ref
    assert len(naming) == 6
    assert "L_2  (#1)" in naming.describe()


def test_tagged_only_when_more_than_one_branch_is_partitioned():
    fp = _FakeFixedPoint(period=3)
    key_a = _key(fp, orbit=0)
    key_b = _key(fp, orbit=1)
    single = ElementNaming(
        HomotopyPartition.from_results(
            [_result(key_a, "left", [0, 1, 2]), _result(key_a, "right", [0, 1, 2])]
        )
    )
    assert not single.tagged
    assert [n.text for n in single.names] == ["L_1", "L_2", "R_1", "R_2"]

    two = ElementNaming(
        HomotopyPartition.from_results(
            [_result(key_a, "left", [0, 1, 2]), _result(key_b, "left", [0, 1])]
        )
    )
    assert two.tagged
    assert [n.text for n in two.names] == ["p3@0.0:L_1", "p3@0.0:L_2", "p3@1.0:L_1"]
    assert two.names[0].mathtext == "p3@0.0:$L_{1}$"
    assert two.lookup("p3@1.0:L_1") == ElementRef(key_b, "left", 0)
    # The tag matches the branch text of ElementRef.label.
    assert ElementRef(key_b, "left", 0).label.startswith(two.names[2].tag)
    assert "branch-tagged" in two.describe()


def test_missing_parent_id_reads_as_own_parent():
    fp = _FakeFixedPoint()
    key = _key(fp)
    homotopy = HomotopyPartition.from_results([_result(key, "left", [0, 1, 2])])
    iterated = IteratedHomotopyPartition(
        [_result(key, "left", [0, 1, 2], parents=[None, None])], homotopy=homotopy
    )
    naming = ElementNaming(homotopy, iterated)
    assert [n.text for n in naming.names] == ["L_1", "L_2"]
    assert naming.parent_of(ElementRef(key, "left", 1)) == ElementRef(key, "left", 1)


def test_bad_parent_raises_value_error():
    fp = _FakeFixedPoint()
    key = _key(fp)
    homotopy = HomotopyPartition.from_results([_result(key, "left", [0, 1, 2])])
    bad_parent = IteratedHomotopyPartition(
        [_result(key, "left", [0, 0.5, 1, 2], parents=[0, 0, 7])], homotopy=homotopy
    )
    with pytest.raises(ValueError, match="parent element #7"):
        ElementNaming(homotopy, bad_parent)

    other_side = IteratedHomotopyPartition(
        [_result(key, "right", [0, 1, 2], parents=[0, 1])], homotopy=homotopy
    )
    with pytest.raises(ValueError, match="homotopy partition does not"):
        ElementNaming(homotopy, other_side)


def test_names_are_deterministic_across_builds():
    _fp, _key_, homotopy, iterated = _one_branch_families()
    first = ElementNaming(homotopy, iterated)
    second = ElementNaming(homotopy, iterated)
    assert first.names == second.names
    assert first.refs == second.refs
    assert first.describe() == second.describe()


# --------------------------------------------------------------------------- #
# Fixture tests: the k=10 combinatorics of the plan, pinned by name text
# --------------------------------------------------------------------------- #
@pytest.fixture
def k10_naming(k10_partitioned):
    session, fp = k10_partitioned
    pieces = build_pieces(session, [fp])
    return pieces, ElementNaming(pieces.homotopy, pieces.iterated)


def test_k10_homotopy_names(k10_naming):
    _pieces, naming = k10_naming
    assert not naming.tagged
    assert {n.text for n in naming.homotopy_names} == {
        "L_1", "L_2", "L_3", "R_1", "R_2", "R_3",
    }


def test_k10_iterated_names(k10_naming):
    _pieces, naming = k10_naming
    assert {n.text for n in naming.names} == {
        "L_1^1", "L_1^2", "L_2", "L_3",
        "R_1^1", "R_1^2", "R_1^3", "R_2", "R_3^1", "R_3^2", "R_3^3",
    }
    assert len(naming) == 11
    # Split parents and their children, by name.
    assert [naming.name(c).text for c in naming.children_of(_homotopy(naming, "L_1"))] == ["L_1^1", "L_1^2"]
    assert naming.parent_of(naming.lookup("L_1^2")) == _homotopy(naming, "L_1")
    assert [naming.name(c).text for c in naming.children_of(_homotopy(naming, "R_1"))] == ["R_1^1", "R_1^2", "R_1^3"]
    assert [naming.name(c).text for c in naming.children_of(_homotopy(naming, "R_3"))] == ["R_3^1", "R_3^2", "R_3^3"]
    for text in ("L_2", "L_3", "R_2"):
        assert [naming.name(c).text for c in naming.children_of(_homotopy(naming, text))] == [text]


def _homotopy(naming: ElementNaming, text: str) -> ElementRef:
    """The homotopy ref printing ``text``."""
    for ref, name in zip(naming.homotopy_refs, naming.homotopy_names):
        if name.text == text:
            return ref
    raise KeyError(text)


def test_k10_names_agree_with_the_partition_structure(k10_naming):
    pieces, naming = k10_naming
    for result in pieces.iterated.as_list():
        for interval in result.intervals:
            ref = result.ref(interval.element_id)
            name = naming.name(ref)
            parent = naming.parent_of(ref)
            assert parent.element_id == interval.parent_element_id
            assert parent.branch_key == ref.branch_key and parent.side == ref.side
            assert name.subscript == interval.parent_element_id + 1
            assert name.side_letter == result.side[0].upper()
            assert naming.ref_of(name) == ref
            assert naming.lookup(name.text) == ref
            assert naming.homotopy_name(parent).subscript == name.subscript
    # Every iterated element lies within its parent's span.
    for parent in naming.homotopy_refs:
        span = pieces.homotopy.element(parent)
        for child in naming.children_of(parent):
            child_iv = pieces.iterated.element(child)
            assert child_iv.lo_cdist >= span.lo_cdist - pieces.homotopy.tol
            assert child_iv.hi_cdist <= span.hi_cdist + pieces.homotopy.tol


def test_k10_describe_lists_every_parent_with_its_children(k10_naming):
    _pieces, naming = k10_naming
    report = naming.describe()
    assert "R_1  (#0) -> R_1^1, R_1^2, R_1^3" in report
    assert "L_3  (#2) -> L_3" in report
    assert "11 iterated element(s)" in report


@pytest.mark.slow
def test_p3_names_are_branch_tagged(p3_partitioned):
    """More than one stable branch is partitioned, so every name carries its tag."""
    session, fp3, fp1 = p3_partitioned
    pieces = build_pieces(session, [fp3, fp1])
    naming = ElementNaming(pieces.homotopy, pieces.iterated)
    assert naming.tagged
    assert len(pieces.homotopy.branch_keys) > 1
    tags = {name.tag for name in naming.names}
    assert len(tags) == len(pieces.homotopy.branch_keys)
    for ref, name in naming.items():
        assert name.text.startswith(f"{name.tag}:")
        assert ref.label.startswith(name.tag)
        assert naming.ref_of(name) == ref
        assert naming.lookup(name.text) == ref
    assert len({n.text for n in naming.names}) == len(naming)

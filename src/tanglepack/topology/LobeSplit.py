"""
Classes one pair of elements does not name, told apart by what their lobes hold.

A bridge with both ends on one stable branch closes with the stable segment
between its ends into a Jordan curve, its LOBE (the arc meets no stable arc
inside, the segment no part of the arc). Two such arcs between the same two
elements are homotopic exactly when their lobes hold the same markers (holes
and periodic points): what lies between them is the difference of their lobes.
On an ordinary trellis every arc of a pair holds the same markers, so the pair
is the class. Where a stable branch ENDS at its fixed point (a cone point of a
quotient, e.g. the hydrogen map's identified plane), an arc can pass either
side of the island that branch bounds and still join the same two elements;
its lobe then holds the island or not (2026-11-10, hydrogen E = 0.5,
B = 2.141, Mitchell 2009 Fig. 8: the paper's b and g).

:func:`split_by_lobes` reads every member's lobe on a dual graph built from
the pair-named classes (the kept bridges, so the arrangement, do not depend on
how the others are classed) and splits a pair whose members fall into several
groups of agreeing lobes; :meth:`LobeSplit.lobe_for` names a walk piece of
such a pair by its own lobe. Both are combinatorial: :meth:`DualGraph.enclosed`
on the arrangement, no geometry.

A marker on the curve itself (the fixed point an arc starts from), one the
model cannot place, or one inside an island of a face the arc crosses (the
arrangement does not say which side of the island the arc passes) decides
nothing: two lobes agree when no marker is held by
one and left out by the other.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional, TYPE_CHECKING

from .BridgeClass import BridgeClass, BridgeClassTable
from .DualWalk import _half_toward, _stable_edges_between
from .StablePartition import row_of_end

if TYPE_CHECKING:
    from ..numerics.Bridge import BridgeId
    from .DualGraph import Chord, DualGraph, FaceNode, Marker

logger = logging.getLogger(__name__)

#: ``{marker: held}`` for the markers a lobe decides.
Held = dict


def _decided(routings) -> Optional[Held]:
    """The markers every routing agrees on. A face with an island (a boundary
    component apart from its ring) lets a chord pass either side of it; the
    arrangement does not say which, so such a marker decides nothing."""
    if not routings:
        return None
    first, *rest = routings
    return {
        marker: held
        for marker, held in first.items()
        if held is not None and all(other.get(marker) == held for other in rest)
    }


def _node_by_edge(dual: "DualGraph") -> dict:
    index = dual.__dict__.get("_lobe_node_index")
    if index is None:
        index = {node.edge_key: node for node in dual.stable_nodes.values()}
        dual.__dict__["_lobe_node_index"] = index
    return index


def chord_lobe(dual: "DualGraph", face: "FaceNode", chord: "Chord") -> Optional[Held]:
    """
    The markers held by a chord closed with the stable segment between its ends.

    Args:
        dual: The dual graph.
        face: The face the chord runs through.
        chord: ``(edge key, side, edge key, side)``, its two ends.

    Returns:
        ``{marker: held}`` for the markers the lobe decides; an empty dict for
        a chord with both ends on one edge; None when its ends are on two
        stable branches or the curve is beyond the model of
        :meth:`~.DualGraph.DualGraph.enclosed`.
    """
    nodes = _node_by_edge(dual)
    a, b = nodes.get(chord[0]), nodes.get(chord[2])
    if a is None or b is None or a.branch_key != b.branch_key:
        return None
    if a is b:
        return {}
    lo, hi = sorted((a.lo_cdist, a.hi_cdist, b.lo_cdist, b.hi_cdist))[1:3]
    between = _stable_edges_between(dual, a.branch_key, lo, hi) - {a.edge_key, b.edge_key}
    return _decided(dual.enclosed(between, [(face, chord)], {_half_toward(a, b), _half_toward(b, a)}))


def chord_of(dual: "DualGraph", bridge_id: "BridgeId") -> Optional[tuple["FaceNode", "Chord"]]:
    """
    A bridge the minimal trellis does not keep, as a chord of the face it lies in.

    Each end is placed by registry order alone: the stable edge its crossing
    lies strictly inside, on its row. None when an end is a vertex of the
    arrangement or on no stable edge, or the two ends read in two faces.
    """
    trellis = dual.trellis
    unstable_key = trellis.unstable_key_of(bridge_id)
    ends = []
    for endpoint, crossing_id in zip(("first", "second"), bridge_id):
        cdist = trellis.intersection(crossing_id).stable_cdist
        row = row_of_end(trellis, bridge_id, endpoint)
        ends += [
            (node.face_on(row), (node.edge_key, row))
            for node in dual.stable_nodes_on(trellis.stable_key_at(crossing_id, unstable_key))
            if node.lo_cdist < cdist < node.hi_cdist and row in node.sides
        ]
    if len(ends) != 2 or ends[0][0] is not ends[1][0]:
        return None
    return ends[0][0], ends[0][1] + ends[1][1]


def bridge_lobe(dual: "DualGraph", bridge_id: "BridgeId") -> Optional[Held]:
    """
    The markers a bridge's lobe holds.

    A kept bridge is an edge of the arrangement and closes with the stable
    edges between its ends; any other is a chord of its face
    (:func:`chord_of`). None for a bridge across two stable branches, or one
    the model cannot place.
    """
    trellis = dual.trellis
    bridge_id = tuple(bridge_id)
    unstable_key = trellis.unstable_key_of(bridge_id)
    branches = {trellis.stable_key_at(end, unstable_key) for end in bridge_id}
    if len(branches) != 1:
        return None
    edge = dual.unstable_edges.get(bridge_id)
    if edge is not None:
        lo, hi = sorted(float(trellis.intersection(end).stable_cdist) for end in bridge_id)
        along = _stable_edges_between(dual, branches.pop(), lo, hi) | {edge}
        return _decided(dual.enclosed(along))
    placed = chord_of(dual, bridge_id)
    return None if placed is None else chord_lobe(dual, *placed)


def _agree(one: Held, other: Held) -> bool:
    return all(other.get(marker, held) == held for marker, held in one.items())


@dataclass(frozen=True)
class LobeSplit:
    """
    One pair of elements holding several classes.

    Attributes:
        pair: The pair, with no lobe.
        groups: ``(lobe, held)`` per class: the markers its members' lobes
            hold, and every marker they decide.
    """

    pair: BridgeClass
    groups: tuple[tuple[frozenset, tuple], ...]

    def lobe_for(self, held: Optional[Held]) -> Optional[frozenset]:
        """The lobe of the one class whose members agree with ``held``; None
        when ``held`` is unread or agrees with none or several."""
        if held is None:
            return None
        fits = [lobe for lobe, decided in self.groups if _agree(held, dict(decided))]
        return fits[0] if len(fits) == 1 else None


def split_by_lobes(
    dual: "DualGraph", table: BridgeClassTable
) -> tuple[dict[BridgeClass, LobeSplit], dict["BridgeId", frozenset]]:
    """
    Split every pair whose members' lobes disagree.

    Members are grouped in table order: one joins the first group it agrees
    with (it agrees with every member there on every marker both decide), and
    starts a new group otherwise. A pair with one group is one class. A
    member whose lobe cannot be read joins the first group, with a WARNING.

    Args:
        dual: A dual graph over the pair-named classes.
        table: Their evidence table.

    Returns:
        ``({pair: LobeSplit}, {bridge id: lobe})`` for the pairs that split.
    """
    splits: dict[BridgeClass, LobeSplit] = {}
    lobes: dict["BridgeId", frozenset] = {}
    for entry in table:
        pair = entry.bridge_class
        if pair.source.branch_key != pair.target.branch_key:
            continue
        members = [member.bridge_id for member in entry.members if not member.is_loop]
        if len(members) < 2:
            continue
        groups: list[tuple[Held, list]] = []
        unread = []
        for bridge_id in members:
            held = bridge_lobe(dual, bridge_id)
            if held is None:
                unread.append(bridge_id)
                continue
            fits = [group for group in groups if _agree(held, group[0])]
            if len(fits) > 1:
                logger.warning(
                    "the lobe of bridge %s agrees with several groups of class %s; "
                    "it joins the first",
                    bridge_id, pair.label,
                )
            if fits:
                fits[0][0].update(held)
                fits[0][1].append(bridge_id)
            else:
                groups.append((dict(held), [bridge_id]))
        if len(groups) < 2:
            continue
        names = [frozenset(marker for marker, inside in held.items() if inside) for held, _ in groups]
        if len(set(names)) < len(names):
            logger.warning(
                "class %s: member lobes disagree but two groups hold the same markers; "
                "not split", pair.label,
            )
            continue
        if unread:
            logger.warning(
                "class %s splits by its lobes, but the lobes of %s cannot be read; "
                "they join the first group", pair.label, unread,
            )
            groups[0][1].extend(unread)
        split = LobeSplit(
            pair, tuple((name, tuple(sorted(held.items(), key=str))) for name, (held, _) in zip(names, groups))
        )
        splits[pair] = split
        for name, (_held, ids) in zip(names, groups):
            for bridge_id in ids:
                lobes[tuple(bridge_id)] = name
        logger.info(
            "class %s holds %d classes, told apart by their lobes: %s",
            pair.label, len(groups),
            "; ".join(
                f"{pair.with_lobe(name).label} ({len(ids)} bridges)"
                for name, (_held, ids) in zip(names, groups)
            ),
        )
    return splits, lobes

"""
Saving and loading a session: the workbench, and beside it what the session remembers.

A session save is the workbench's save (:mod:`tanglepack.numerics.persist`, the
same two files) plus a ``"session"`` section: the chosen strong pips, the class
letters, the resonance zones with their frozen boundary polygons, and the last
blast growth factor per zone. These are history, not functions of the
workbench: a pip is chosen, a letter is assigned in the order classes were
seen, and ``blast_zone`` tests membership against the polygon captured when the
zone was defined. The topology and the session's caches are not saved; the
first ``topology()`` after a load recomputes them.
"""

from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np

from ..numerics.BranchPoint import BranchPoint
from ..numerics.PointStore import End
from ..numerics.persist import (
    Decoder,
    Encoder,
    encode_workbench,
    floats_of,
    read,
    restore_workbench,
    typed_floats,
    write,
)
from ..topology.TopologyResults import Arc
from .ResonanceZone import ResonanceZone

if TYPE_CHECKING:
    from .TangleSession import TangleSession


def save_session(session: TangleSession, path: str | os.PathLike) -> None:
    """
    Write a session to ``<path>.npz`` and ``<path>.json``.

    Args:
        session: The session to save.
        path: The stem of the two files.
    """
    encoder = Encoder(session.workbench)
    meta, arrays = encode_workbench(session.workbench, encoder)
    fixed_point, key, crossing, end = encoder.fixed_point, encoder.key, encoder.crossing, encoder.end
    zones = list(session.resonance_zones.values())
    meta["session"] = {
        "pips": [[fixed_point(fp), crossing(pip)] for fp, pip in session._strong_pips.items()],
        # an identity is a frozenset of (branch key, side, lo id, hi id), plus ("lobe",
        # markers) for a class its pair does not name alone: sorted, to save one order
        "letters": [
            [sorted((_encode_identity_item(item, key) for item in identity), key=str), letter]
            for identity, letter in session.bridge_alphabet.assigned.items()
        ],
        "zones": [{
            "fixed_point": fixed_point(zone.fixed_point),
            "stable_branch_key": key(zone.stable_branch_key),
            "boundary_intersection": crossing(zone.boundary_intersection),
            "inside": None if zone.inside is None else fixed_point(zone.inside),
            "cut_intersections": [crossing(cut) for cut in zone.cut_intersections],
            "intersection_ids": list(zone.intersection_ids),
            "previous_tails": [
                [key(branch_key), end(tail if isinstance(tail, BranchPoint) else tail.id)]
                for branch_key, tail in zone.previous_tails.items()
            ],
            "boundary_arcs": [
                [arc.kind, arc.lo_id, arc.hi_id, key(arc.branch_key), arc.bridge_id, arc.reverse]
                for arc in zone.boundary_arcs
            ],
            "has_vertices": zone.boundary_vertices is not None,
        } for zone in zones],
        "blast_stats": [fixed_point(fp) for fp in session._blast_stats],
    }
    for i, zone in enumerate(zones):
        if zone.boundary_vertices is not None:
            arrays[f"zone{i}_vertices"] = np.asarray(zone.boundary_vertices)
    arrays["blast_factors"], meta["session"]["blast_factor_types"] = typed_floats(
        session._blast_stats.values())
    write(path, "session", meta, arrays, encoder)


def load_session(session: TangleSession, path: str | os.PathLike) -> None:
    """
    Put a saved session into a freshly built one.

    Args:
        session: A session just built around the save's map.
        path: The stem of the two files.

    Raises:
        ValueError: The save is of another format, holds an object array, is
            a workbench save (load that with ``TangleWorkbench.load``), or the
            session's map disagrees with the saved points.
    """
    meta, arrays = read(path)
    if meta["kind"] != "session":
        raise ValueError(f"{os.fspath(path)} is a {meta['kind']} save: load it with TangleWorkbench.load")
    decoder = restore_workbench(session.workbench, meta, arrays)
    saved = meta["session"]
    fixed_points, key, crossings = decoder.fixed_points, decoder.key, decoder.crossings
    session._strong_pips = {fixed_points[fp]: crossings[pip] for fp, pip in saved["pips"]}
    session.bridge_alphabet.assigned = {
        frozenset(_decode_identity_item(item, key) for item in identity): letter
        for identity, letter in saved["letters"]
    }
    for i, zone in enumerate(saved["zones"]):
        fixed_point = fixed_points[zone["fixed_point"]]
        session.resonance_zones[fixed_point] = ResonanceZone(
            fixed_point=fixed_point,
            stable_branch_key=key(zone["stable_branch_key"]),
            boundary_intersection=crossings[zone["boundary_intersection"]],
            cut_intersections=[crossings[cut] for cut in zone["cut_intersections"]],
            intersection_ids=zone["intersection_ids"],
            previous_tails={key(branch): _tail(decoder, *end) for branch, end in zone["previous_tails"]},
            boundary_arcs=[
                Arc(kind=kind, lo_id=lo, hi_id=hi, branch_key=key(branch_key),
                    bridge_id=None if bridge_id is None else tuple(bridge_id), reverse=reverse)
                for kind, lo, hi, branch_key, bridge_id, reverse in zone["boundary_arcs"]
            ],
            boundary_vertices=arrays[f"zone{i}_vertices"] if zone["has_vertices"] else None,
            # A save before format 4 declared no nesting.
            inside=None if zone.get("inside") is None else fixed_points[zone["inside"]],
        )
    factors = floats_of(arrays["blast_factors"], saved["blast_factor_types"])
    session._blast_stats = {fixed_points[fp]: factor for fp, factor in zip(saved["blast_stats"], factors)}


def _tail(decoder: Decoder, row: int, fixed_point: int, orbit_index: int) -> End | BranchPoint:
    """A zone's saved tail as the zone held it: the row's ``End``, or the periodic point."""
    end = decoder.end(row, fixed_point, orbit_index)
    return decoder.store.end(end) if isinstance(end, int) else end


def _encode_identity_item(item: tuple, key) -> list:
    """One item of a class identity: an element ``(branch, side, lo, hi)`` or
    ``("lobe", markers)`` (``Topology.identity``)."""
    if item[0] == "lobe":
        return ["lobe", sorted((_encode_marker(marker) for marker in item[1]), key=str)]
    branch, side, lo, hi = item
    return [key(branch), side, lo, hi]


def _decode_identity_item(item: list, key) -> tuple:
    """The inverse of :func:`_encode_identity_item`."""
    if item[0] == "lobe":
        return ("lobe", frozenset(_decode_marker(marker) for marker in item[1]))
    branch_key, side, lo, hi = item
    return (key(branch_key), side, lo, hi)


def _encode_marker(marker: tuple) -> list:
    """A dual-graph marker: ``("hole", origin ids, iterate)`` or ``("orbit", position, index)``."""
    if marker[0] == "hole":
        return ["hole", list(marker[1]), marker[2]]
    return list(marker)


def _decode_marker(marker: list) -> tuple:
    if marker[0] == "hole":
        return ("hole", tuple(marker[1]), marker[2])
    return tuple(marker)

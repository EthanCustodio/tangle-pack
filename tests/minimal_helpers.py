"""Shared builders for the minimal-trellis / partition-family / dual-graph tests."""

from __future__ import annotations

from dataclasses import dataclass

from tanglepack.topology.BridgeClass import BridgeClassTable
from tanglepack.topology.MinimalTrellis import MinimalTrellis, minimal_trellis
from tanglepack.topology.PartitionFamily import (
    HomotopyPartition,
    IteratedHomotopyPartition,
)
from tanglepack.topology.Trellis import Trellis


@dataclass
class Pieces:
    """Everything the new layer is built from, for one partitioned session."""

    session: object
    fixed_points: list
    full: Trellis
    holes: list
    table: BridgeClassTable
    homotopy: HomotopyPartition
    minimal: MinimalTrellis
    iterated: IteratedHomotopyPartition

    @property
    def strong_pips(self) -> list[int]:
        return [
            self.session.trellis(fp).strong_pip
            for fp in self.fixed_points
            if self.session.trellis(fp).strong_pip is not None
        ]


def build_pieces(session, fixed_points) -> Pieces:
    """Build the minimal trellis and both partition families directly (no session caches)."""
    fixed_points = list(fixed_points)
    full = session.trellis()
    holes = [hole for fp in fixed_points for hole in session.trellis(fp).holes]
    table = session.bridge_classes()
    homotopy = HomotopyPartition.from_results(
        session._gathered_partitions(), trellis=full
    )
    minimal = minimal_trellis(full, holes, table)
    iterated = IteratedHomotopyPartition.from_minimal(minimal, homotopy)
    return Pieces(session, fixed_points, full, holes, table, homotopy, minimal, iterated)

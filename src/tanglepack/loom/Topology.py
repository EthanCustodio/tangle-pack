"""
One recompute of the topology: stages 3-9 of the chalkboard for one numerical
state and one strong-pip choice.

Dev Notes:

* Why it exists (SIMP-8, WP8). Topology is a pure function of the numerical
  state and the pip choice (D1), so the session holds ONE result per
  ``(generation, fixed points, pips)``
  (:meth:`~tanglepack.loom.TangleSession.TangleSession.topology`) instead of a
  hand-rolled cache per stage keyed on a partition signature. Stages 3-6
  (pseudoneighbors, holes, stable partition) are recomputed from scratch on
  the selection's trellis when the object is built; stages 7-9 (classes,
  minimal trellis, iterated partition, dual graph, words) are lazy and
  computed once.
* Partitions, holes and pips come from the SELECTED fixed points only
  (API-2), and every stable branch of the selection is partitioned, so the
  old "branch missing a partition" WARNING is gone: a heteroclinic bridge's
  class needs both its fixed points selected, and ``bridge_classes`` raises
  naming the end it cannot resolve otherwise (``henon_cases.build_heteroclinic``
  is the first such case).
* The selection's trellis is cut at every selected pip's orbit (D10, D12) by
  ``TangleSession.trellis`` itself: a per-fixed-point trellis when its pip is
  set, a selection of several at its fixed points' pips on every read. Stage 3
  (the pip) runs first on each fixed point's own trellis, so a default pip is
  chosen before the selection is read. A class member may be a joined arc of
  the cut, which the workbench does not hold, so bridges are resolved through
  the trellis (``bridge_between``).
* Stages 4-6 run on the SELECTION's trellis (WP9): a heteroclinic pair of
  ``W^u(A)`` crossings on B's stable branch is a pair of B (D7) only on a
  trellis holding both tangles, and its backward chain runs into A's bridges.
  So heteroclinic pairs and their holes exist only when both tangles are
  selected, ``topology()`` and ``topology(fp)`` differ exactly where they do,
  and after ``topology()`` of several fixed points the per-fixed-point
  trellises carry no stage output (``topology(fp)`` fills one). For one fixed
  point the selection's trellis IS its trellis.
* Letters follow :meth:`Topology.identity`, the pair of the bounding crossings
  of a class's two elements, so a class keeps its letter across recomputes
  while its elements keep their bounds (C2, BCG-7). Element NAMES
  (``R_(o.b;i)``) stay positional.
"""

from __future__ import annotations

import logging
from functools import cached_property
from typing import TYPE_CHECKING, Optional, Sequence

from ..topology.BridgeClass import (
    BridgeClass,
    BridgeClassEntry,
    BridgeClassTable,
    bridge_classes,
    zone_label,
)
from ..topology.DualGraph import DualGraph
from ..topology.MinimalTrellis import MinimalTrellis, minimal_trellis
from ..topology.PartitionFamily import HomotopyPartition, IteratedHomotopyPartition
from ..topology.SymbolicDynamics import SymbolicDynamics, symbolic_dynamics

if TYPE_CHECKING:
    from ..numerics.FixedPoint import FixedPoint
    from ..topology.TopologyResults import StablePartitionResult
    from ..topology.Trellis import Trellis
    from .TangleSession import TangleSession

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())


class Topology:
    """
    Stages 3-9 of the chalkboard for one numerical state and pip choice (D1).

    Attributes:
        generation: The workbench generation it was computed at.
        fixed_points: The selected fixed points, in workbench order.
        trellis: The selection's trellis, cut at the pip orbits (its bridges
            are the ones classed).
        strong_pips: The strong pip of each selected fixed point (None where
            a fixed point has no candidate), in the same order.
        holes: Every hole of the selection (a heteroclinic reference's chain
            holes both tangles).
        partitions: Their stable partitions, one per ``(branch_key, side)``.
    """

    def __init__(self, session: "TangleSession", fixed_points: Sequence["FixedPoint"]) -> None:
        """
        Choose every selected fixed point's pip, then recompute stages 4-6 on
        the selection's trellis.

        Args:
            session: The session whose workbench, pip choice, zones and
                alphabet the topology reads.
            fixed_points: The selection, in workbench order.

        Raises:
            ValueError: A fixed point's chosen strong pip is no longer a
                candidate (see ``TangleSession.set_strong_pip``).
        """
        self._session = session
        self.generation = session.workbench.generation
        self.fixed_points = tuple(fixed_points)
        self.strong_pips = tuple(self._choose_pip(fp) for fp in self.fixed_points)
        self.trellis: "Trellis" = session.trellis(self.fixed_points)  # cut at strong_pips
        self.trellis.compute_pseudoneighbors()
        self.trellis.punch_holes()
        self.trellis.partition_stable_manifold()
        # New lists: a later manual stage call clears the trellis slots in place.
        self.holes = list(self.trellis.holes)
        self.partitions: list["StablePartitionResult"] = list(self.trellis.stable_partitions)

    def _choose_pip(self, fixed_point: "FixedPoint") -> Optional[int]:
        """Stage 3 on the fixed point's own trellis: the pip the session carries
        (the default, remembered, when none was chosen)."""
        session = self._session
        trellis = session.trellis(fixed_point)
        if trellis.strong_pip is None:
            session.classify_strong_pips(fixed_point)
        if trellis.strong_pip is not None:
            session.set_strong_pip(fixed_point, trellis.strong_pip)  # one set on the trellis is adopted
        elif session.chosen_strong_pip(fixed_point) is not None:
            raise ValueError(
                f"the chosen strong pip of {fixed_point!r} is no longer a "
                "strong-pip candidate; choose another with set_strong_pip()"
            )
        return trellis.strong_pip

    @cached_property
    def homotopy(self) -> HomotopyPartition:
        """The stable partitions as a family (stage 6)."""
        return HomotopyPartition.from_results(self.partitions, trellis=self.trellis)

    @cached_property
    def minimal(self) -> MinimalTrellis:
        """The minimal trellis over every bridge's class, the evidence table (stages 7-8)."""
        return minimal_trellis(
            self.trellis, self.holes, bridge_classes(self.trellis, self.partitions)
        )

    @cached_property
    def classes(self) -> BridgeClassTable:
        """The ALPHABET (D2): the minimal trellis's classes, active ones lettered
        by :meth:`identity` from the session alphabet, every one zoned."""
        table = self.minimal.alphabet
        for entry in table:
            if entry.active:
                entry.letter = self.letter_for(entry.bridge_class)
            entry.zone_key = self._zone_key(entry)
        return table

    def _zone_key(self, entry: BridgeClassEntry) -> Optional["FixedPoint"]:
        """The resonance zone every member bridge of a class lies in, or None
        (with a WARNING when the members disagree)."""
        session = self._session
        keys = set()
        for member in entry.members:
            zone = session.classify_bridge(self.trellis.bridge_between(*member.bridge_id))
            keys.add(None if zone is None else zone.key)
        if len(keys) > 1:
            logger.warning(
                "bridge class %s straddles resonance zones: members fall in %s",
                entry.bridge_class.label,
                sorted("none" if key is None else zone_label(key) for key in keys),
            )
            return None
        return next(iter(keys), None)

    @cached_property
    def iterated(self) -> IteratedHomotopyPartition:
        """The homotopy partition cut by the minimal trellis's image bridges."""
        return IteratedHomotopyPartition.from_minimal(self.minimal, self.homotopy)

    @cached_property
    def dual_graph(self) -> DualGraph:
        """The dual graph, unified on each pip's fundamental segment."""
        pips = [pip for pip in self.strong_pips if pip is not None]
        return DualGraph(self.minimal, self.iterated, strong_pips=pips)

    @cached_property
    def dynamics(self) -> SymbolicDynamics:
        """The words, closed under the word map, the closure lettered like the alphabet."""
        return symbolic_dynamics(self.dual_graph, self.classes, letter_for=self.letter_for)

    def identity(self, bridge_class: BridgeClass) -> frozenset:
        """
        A class across recomputes: the identities of its two elements (C2).

        Args:
            bridge_class: A class of this topology (its refs name homotopy
                elements).

        Returns:
            The unordered pair of
            :attr:`~tanglepack.topology.TopologyResults.PartitionInterval.identity`.
        """
        return frozenset(
            self.homotopy.element(ref).identity
            for ref in (bridge_class.source, bridge_class.target)
        )

    def letter_for(self, bridge_class: BridgeClass) -> str:
        """The session alphabet's letter for a class, by :meth:`identity`."""
        return self._session.bridge_alphabet.letter_for(self.identity(bridge_class))


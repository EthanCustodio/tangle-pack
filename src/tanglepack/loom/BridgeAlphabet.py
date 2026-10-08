"""
The persistent alphabet that names bridge classes.

A :class:`~tanglepack.topology.BridgeClass.BridgeClass` is a symbol; an
alphabet gives it a letter. The alphabet lives on the session and outlives any
one recompute: a class seen before keeps its letter, a class never seen before
takes the next unused one, so growing, blasting or re-partitioning "moves
through the alphabet" rather than reshuffling it.

Dev Notes:

* A class is named by its IDENTITY, the pair of its two elements'
  :attr:`~tanglepack.topology.TopologyResults.PartitionInterval.identity`
  (their bounding crossings; ``loom.Topology.identity``), never by the
  ``BridgeClass`` object: an ``ElementRef``'s ``element_id`` is positional,
  so a new hole renumbered every element after it and re-lettered classes
  that had not changed (BCG-7). A class keeps its letter while its two
  elements keep their bounding crossings (C2, WP8). Any hashable works.
* Letters run ``a``..``t`` then ``aa``, ``ab``, ... (the active series of
  :func:`~tanglepack.topology.SymbolicDynamics.letter`), so ``a`` and ``aa``
  never collide with a refined symbol such as ``a_1``.
* This alphabet letters ACTIVE classes only, and only the alphabet's (the
  minimal trellis's classes, D2). Inert classes take the disjoint series
  ``u..z, uu, ...`` from ``topology/SymbolicDynamics.inert_letters``: the two
  series share one generator over disjoint digits, so they never collide and
  nothing is skipped (SIMP-6, WP6).
"""

from __future__ import annotations

from typing import Hashable

from ..topology.SymbolicDynamics import letter


class BridgeAlphabet:
    """
    Letters for bridge classes, assigned once and remembered.

    Attributes:
        assigned: The letter of every class identity named so far, in
            assignment order.
    """

    def __init__(self) -> None:
        """Start with no letters assigned."""
        self.assigned: dict[Hashable, str] = {}

    def __len__(self) -> int:
        return len(self.assigned)

    def __contains__(self, identity: object) -> bool:
        return identity in self.assigned

    def letter_for(self, identity: Hashable) -> str:
        """
        The letter of a class identity, assigning the next unused one on first sight.

        Args:
            identity: The class to name, by its identity (the session passes
                ``Topology.identity(bridge_class)``).

        Returns:
            Its letter.
        """
        existing = self.assigned.get(identity)
        if existing is not None:
            return existing
        new = letter(len(self.assigned))
        self.assigned[identity] = new
        return new

    def class_of(self, letter_: str) -> Hashable:
        """
        The class identity a letter names.

        Args:
            letter_: A letter previously handed out.

        Returns:
            The identity.

        Raises:
            KeyError: If no class carries that letter.
        """
        for identity, assigned in self.assigned.items():
            if assigned == letter_:
                return identity
        raise KeyError(f"no bridge class is lettered {letter_!r}")

    def reset(self) -> None:
        """Forget every assignment; the next class named gets ``a`` again."""
        self.assigned.clear()

"""
The persistent alphabet that names bridge classes.

A :class:`~tanglepack.topology.BridgeClass.BridgeClass` is a symbol; an
alphabet gives it a letter. The alphabet lives on the session and outlives any
one classing: a class seen before keeps its letter, a class never seen before
takes the next unused one, so growing, blasting or re-partitioning "moves
through the alphabet" rather than reshuffling it.

Dev Notes:

* Identity is the class object (its two ``ElementRef`` s, compared componentwise
  with the fixed point by identity). A re-partition that renumbers element ids
  changes the pair, so the "same" homotopy class comes back as a new class with
  a new letter. Carrying identity across a refinement is the refined-symbols
  work, not this module's.
* Letters run ``a``..``z`` then ``aa``, ``ab``, ... (bijective base 26), so
  ``a`` and ``aa`` never collide with a future refined symbol such as ``a_1``.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..topology.BridgeClass import BridgeClass


def letter(index: int) -> str:
    """
    The ``index``-th letter of the alphabet: ``a``..``z``, ``aa``, ``ab``, ...

    Args:
        index: Zero-based position.

    Returns:
        The letter string.

    Raises:
        ValueError: If ``index`` is negative.
    """
    if index < 0:
        raise ValueError(f"letter index must be non-negative, got {index}")
    chars = []
    index += 1
    while index > 0:
        index, remainder = divmod(index - 1, 26)
        chars.append(chr(ord("a") + remainder))
    return "".join(reversed(chars))


class BridgeAlphabet:
    """
    Letters for bridge classes, assigned once and remembered.

    Attributes:
        assigned: The letter of every class named so far, in assignment order.
    """

    def __init__(self) -> None:
        """Start with no letters assigned."""
        self.assigned: dict["BridgeClass", str] = {}

    def __len__(self) -> int:
        return len(self.assigned)

    def __contains__(self, bridge_class: object) -> bool:
        return bridge_class in self.assigned

    def letter_for(self, bridge_class: "BridgeClass") -> str:
        """
        The class's letter, assigning the next unused one on first sight.

        Args:
            bridge_class: The class to name.

        Returns:
            Its letter.
        """
        existing = self.assigned.get(bridge_class)
        if existing is not None:
            return existing
        new = letter(len(self.assigned))
        self.assigned[bridge_class] = new
        return new

    def class_of(self, letter_: str) -> "BridgeClass":
        """
        The class a letter names.

        Args:
            letter_: A letter previously handed out.

        Returns:
            The class.

        Raises:
            KeyError: If no class carries that letter.
        """
        for bridge_class, assigned in self.assigned.items():
            if assigned == letter_:
                return bridge_class
        raise KeyError(f"no bridge class is lettered {letter_!r}")

    def reset(self) -> None:
        """Forget every assignment; the next class named gets ``a`` again."""
        self.assigned.clear()

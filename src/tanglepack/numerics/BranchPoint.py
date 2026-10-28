"""
Branch points: the points of a periodic orbit, where manifold branches start.

A :class:`BranchPoint` is one orbit point of a :class:`~.FixedPoint.FixedPoint`
seen from the manifolds: the root of up to four branches (one per direction
along each eigenvector) and the anchor of the chains that hang off it. It is a
VIEW over its fixed point and stores nothing of its own:

* its coordinates are the fixed point's ``coordinates[orbit_index]``;
* its canonical distance is 0.0 on both manifolds (the periodic point is where
  every one of its branches starts);
* the chain on each slot lives in the fixed point's slot table, keyed by
  ``(stability, orbit_index, branch)``;
* its image and preimage are the next and previous orbit points.

It is never a row of the :class:`~.PointStore.PointStore` and no chain holds
it: a chain hanging off it names it as its ``anchor``, the left neighbour of
position 0.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import numpy as np
from numpy.typing import NDArray

from .Intersection import Stability

if TYPE_CHECKING:
    from .FixedPoint import FixedPoint
    from .PointStore import Chain


class BranchPoint:
    """
    One orbit point of a periodic orbit, as the root of its manifold branches.

    Attributes:
        fixed_point (FixedPoint): The periodic orbit this point belongs to.
        orbit_index (int): Which orbit point it is.
        num_branches (int): The eigendirection SLOTS at this point, always two
            (one per direction along each eigenvector). Not
            :attr:`FixedPoint.num_branches`, which counts the branches one
            chain of the orbit visits (1, or 2 under inversion).
    """

    num_branches = 2

    def __init__(self, fixed_point: FixedPoint, orbit_index: int) -> None:
        """
        Args:
            fixed_point: The periodic orbit.
            orbit_index: Which of its orbit points this is.
        """
        self.fixed_point = fixed_point
        self.orbit_index = orbit_index

    @property
    def x(self) -> float:
        """The point's x-coordinate."""
        return float(self.fixed_point.coordinates[self.orbit_index][0])

    @property
    def y(self) -> float:
        """The point's y-coordinate."""
        return float(self.fixed_point.coordinates[self.orbit_index][1])

    def get_point(self) -> NDArray[np.float64]:
        """
        The point's coordinates.

        Returns:
            np.ndarray: A read-only ``(2,)`` view of the fixed point's orbit point.
        """
        point = np.asarray(self.fixed_point.coordinates[self.orbit_index], dtype=np.float64).view()
        point.flags.writeable = False
        return point

    def get_cdist(self, stability: Stability) -> float:
        """
        The point's canonical distance along a manifold of the given stability:
        0.0 on both, since every branch starts here.

        Args:
            stability (Stability): ``"unstable"`` or ``"stable"``.

        Raises:
            ValueError: ``stability`` is neither.
        """
        if stability not in ("unstable", "stable"):
            raise ValueError(f"stability is 'unstable' or 'stable', not {stability!r}")
        return 0.0

    def chain(self, stability: Stability, branch: int) -> Optional[Chain]:
        """The chain hanging on one slot (None for an empty slot)."""
        return self.fixed_point.slot_chain(stability, self.orbit_index, branch)

    def attach(self, chain: Chain, stability: Stability, branch: int) -> None:
        """
        Hang a branch's chain on one slot: the slot holds the chain and the
        chain names this point as its anchor.

        Raises:
            ValueError: The slot already holds another chain.
            IndexError: ``branch`` is not 0 or 1.
        """
        self.fixed_point.attach_chain(chain, stability, self.orbit_index, branch)
        chain.anchor = self

    @property
    def forward_branches(self) -> tuple[Optional[Chain], ...]:
        """The chains of the two unstable slots (read-only)."""
        return tuple(self.chain("unstable", b) for b in range(self.num_branches))

    @property
    def backward_branches(self) -> tuple[Optional[Chain], ...]:
        """The chains of the two stable slots (read-only)."""
        return tuple(self.chain("stable", b) for b in range(self.num_branches))

    def iterate(self, n: int = 1) -> BranchPoint:
        """The orbit point ``n`` map steps on (``n < 0``: back)."""
        orbit = self.fixed_point.branch_points
        return orbit[(self.orbit_index + n) % len(orbit)]

    @property
    def next_iterate(self) -> BranchPoint:
        """The point's image: the next orbit point."""
        return self.iterate(1)

    @property
    def prev_iterate(self) -> BranchPoint:
        """The point's preimage: the previous orbit point."""
        return self.iterate(-1)

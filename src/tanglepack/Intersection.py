from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Literal, TYPE_CHECKING

if TYPE_CHECKING:
    from .FixedPoint import FixedPoint

# ManifoldKey = (fixed_point, stability, orbit_index, branch_index)
# Identical to the key type used in TangleWorkbench.manifolds.
ManifoldKey = tuple["FixedPoint", Literal["unstable", "stable"], int, int]

import numpy as np
from numpy.typing import NDArray


class Intersection:
    """
    Represents a single crossing between the stable and unstable manifolds.

    Two instances at the same coordinates are distinct objects (eq=False keeps
    identity-based equality and hash, so Intersections can live in sets/dicts).

    Attributes:
        coords: Geometric (x, y) of the crossing.
        unstable_cdist: Position along the unstable manifold at this crossing.
        stable_cdist: Position along the stable manifold at this crossing.
        seg_ids: The pair of R-tree segment IDs that produced this crossing.
            None for synthetic intersections.
        label: Optional human-readable name.
    """

    def __init__(
        self,
        coords: tuple[float, float] = None,
        unstable_cdist: float = None,
        stable_cdist: float = None,
        seg_ids: Optional[frozenset[int]] = None,
        id: Optional[int] = None,
        label: Optional[str] = None,
        manifold_a_key: Optional[ManifoldKey] = None,
        manifold_b_key: Optional[ManifoldKey] = None,
    ):
        self.coords = coords
        self.unstable_cdist = unstable_cdist
        self.stable_cdist = stable_cdist
        self.seg_ids = seg_ids
        self.id = id
        self.label = label
        self.manifold_a_key = manifold_a_key
        self.manifold_b_key = manifold_b_key

    @classmethod
    def from_segments(
        cls,
        coords: tuple[float, float],
        unstable_cdist: float,
        stable_cdist: float,
        seg1_id: int,
        seg2_id: int,
        manifold_a_key: Optional[ManifoldKey] = None,
        manifold_b_key: Optional[ManifoldKey] = None,
        label: Optional[str] = None,
    ) -> Intersection:
        """Create an Intersection backed by two R-tree segment IDs."""
        return cls(
            coords=coords,
            unstable_cdist=unstable_cdist,
            stable_cdist=stable_cdist,
            seg_ids=frozenset({seg1_id, seg2_id}),
            label=label,
            manifold_a_key=manifold_a_key,
            manifold_b_key=manifold_b_key,
        )

    @property
    def fixed_points(self) -> tuple:
        """Return the distinct FixedPoint objects involved in this intersection."""
        fps = []
        if self.manifold_a_key is not None:
            fps.append(self.manifold_a_key[0])
        if self.manifold_b_key is not None and self.manifold_b_key[0] is not fps[0]:
            fps.append(self.manifold_b_key[0])
        return tuple(fps)

    @classmethod
    def synthetic(
        cls,
        coords: tuple[float, float],
        unstable_cdist: float,
        stable_cdist: float,
        label: Optional[str] = None,
    ) -> Intersection:
        """
        Create an Intersection not backed by a detected segment crossing.

        Use this for:
        - The fixed point itself
        - Manually specified turning points
        - Any crossing you want to declare programmatically
        """
        return cls(coords, unstable_cdist, stable_cdist, None, label)

    # --- helpers ---
    @property
    def is_synthetic(self) -> bool:
        """True if this intersection was not detected from crossing segments."""
        return self.seg_ids is None

    def get_point(self) -> NDArray[np.float64]:
        """Return coords as a (2,) array, consistent with Point.get_point()."""
        return np.array(self.coords, dtype=np.float64)

from __future__ import annotations

from typing import Optional, Literal, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from .BaseManifold import BaseManifold
    from .FixedPoint import FixedPoint
    from .Point import Point

# ManifoldKey = (fixed_point, stability, orbit_index, branch_index)
# Identical to the key type used in TangleWorkbench.manifolds.
ManifoldKey = tuple["FixedPoint", Literal["unstable", "stable"], int, int]


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
        unstable_manifold: The unstable curve the crossing was detected on --
            an indexed manifold or an iterated bridge. ``Tangle.create_bridges``
            groups crossings by it, because a bridge and the manifold it was cut
            out of share a ``manifold_key`` but are distinct polylines.
        unstable_segment: The two adjacent points of ``unstable_manifold`` that
            bracket the crossing, captured at resolve time. These are the head
            and tail candidates a cut at this crossing uses; capturing them here
            keeps them stable even after a separator point is spliced into the
            same segment. None for a synthetic crossing.
    """

    def __init__(
        self,
        coords: Optional[tuple[float, float]] = None,
        unstable_cdist: Optional[float] = None,
        stable_cdist: Optional[float] = None,
        seg_ids: Optional[frozenset[int]] = None,
        id: Optional[int] = None,
        label: Optional[str] = None,
        manifold_a_key: Optional[ManifoldKey] = None,
        manifold_b_key: Optional[ManifoldKey] = None,
        unstable_manifold: Optional["BaseManifold"] = None,
        unstable_segment: Optional[tuple["Point", "Point"]] = None,
    ):
        self.coords = coords
        self.unstable_cdist = unstable_cdist
        self.stable_cdist = stable_cdist
        self.seg_ids = seg_ids
        self.id = id
        self.label = label
        self.manifold_a_key = manifold_a_key
        self.manifold_b_key = manifold_b_key
        self.unstable_manifold = unstable_manifold
        self.unstable_segment = unstable_segment

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
        unstable_manifold: Optional["BaseManifold"] = None,
        unstable_segment: Optional[tuple["Point", "Point"]] = None,
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
            unstable_manifold=unstable_manifold,
            unstable_segment=unstable_segment,
        )

    @property
    def fixed_points(self) -> tuple["FixedPoint", ...]:
        """
        Return the distinct FixedPoint objects involved in this intersection.

        Either key may be absent: a crossing born on an iterated bridge carries
        only ``manifold_b_key`` until the bridge's unstable key is propagated,
        and a synthetic crossing may carry neither.

        Returns:
            The distinct fixed points in (a, b) order: one element for a
            homoclinic crossing, two for a heteroclinic one, and the empty
            tuple when no manifold key is set.
        """
        fps: list["FixedPoint"] = []
        for key in (self.manifold_a_key, self.manifold_b_key):
            if key is None:
                continue
            if not any(key[0] is seen for seen in fps):
                fps.append(key[0])
        return tuple(fps)

    @classmethod
    def synthetic(
        cls,
        coords: tuple[float, float],
        unstable_cdist: float,
        stable_cdist: float,
        label: Optional[str] = None,
        manifold_a_key: Optional[ManifoldKey] = None,
        manifold_b_key: Optional[ManifoldKey] = None,
    ) -> Intersection:
        """
        Create an Intersection not backed by a detected segment crossing.

        Use this for:
        - The periodic point itself (the anchor at cdist (0, 0) of a branch pair)
        - Manually specified turning points
        - Any crossing you want to declare programmatically

        Args:
            coords: Geometric (x, y) of the crossing.
            unstable_cdist: Position along the unstable manifold.
            stable_cdist: Position along the stable manifold.
            label: Optional human-readable name.
            manifold_a_key: Key of the unstable branch this crossing sits on.
            manifold_b_key: Key of the stable branch this crossing sits on.

        Returns:
            An Intersection with ``seg_ids`` and ``id`` both None; the id is
            assigned by :meth:`IntersectionRegistry.add`.

        Note:
            The keys matter for anchors: every branch anchor shares cdist
            (0, 0), so the branch keys are what keeps them apart in the
            registry's collision test.
        """
        return cls(
            coords=coords,
            unstable_cdist=unstable_cdist,
            stable_cdist=stable_cdist,
            seg_ids=None,
            label=label,
            manifold_a_key=manifold_a_key,
            manifold_b_key=manifold_b_key,
        )

    # --- helpers ---
    @property
    def is_synthetic(self) -> bool:
        """True if this intersection was not detected from crossing segments."""
        return self.seg_ids is None

    def get_point(self) -> NDArray[np.float64]:
        """Return coords as a (2,) array, consistent with Point.get_point()."""
        return np.array(self.coords, dtype=np.float64)

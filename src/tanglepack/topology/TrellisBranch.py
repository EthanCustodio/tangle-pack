from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional, TYPE_CHECKING

import numpy as np
from numpy.typing import NDArray

if TYPE_CHECKING:
    from ..FixedPoint import FixedPoint
    from ..Intersection import ManifoldKey

"""
Dev Notes:

A TrellisBranch is the topological view of one branch of one manifold of one
periodic-orbit point. It is intentionally lightweight: it stores the *ordering*
of intersections along the branch (encoded as the order of intersection_ids) and
the eigen/orbit metadata the algorithms need, but it does not hold a reference to
the registry. Anything that needs the underlying Intersection objects or their
cdists goes through the owning Trellis.

intersection_ids is always kept sorted by canonical distance ascending — i.e.
from the anchoring periodic point outward. "Toward the anchor" therefore means
walking toward index 0.

Open question — return_period: for a period-q orbit a single branch maps onto
itself under M^k_value (k_value = q, or 2q with inversion). We expose that as a
convenience, but the precise branch period p used by the Pseudoneighbor Algorithm
should be confirmed against the orbit/inversion bookkeeping when that algorithm
is implemented.
"""


@dataclass
class TrellisBranch:
    """
    One branch of a stable or unstable manifold within a finite trellis.

    A branch corresponds exactly to one ``TangleWorkbench.manifolds`` key. It
    carries the intersections that lie on it, ordered by canonical distance from
    the anchoring periodic point outward, plus the eigen/orbit metadata the
    topological algorithms read.

    Attributes:
        key: The manifold key (fixed_point, stability, orbit_index, branch_index)
            identifying this branch. Identical to the key in TangleWorkbench.manifolds.
        fixed_point: The fixed point this branch is anchored to.
        stability: "unstable" or "stable".
        orbit_index: Index of the anchoring point within the periodic orbit.
        branch_index: 0 or 1 (1 only for fixed points with inversion).
        intersection_ids: Registry IDs of the intersections on this branch,
            sorted ascending by canonical distance (anchor outward).
    """

    key: "ManifoldKey"
    fixed_point: "FixedPoint"
    stability: Literal["unstable", "stable"]
    orbit_index: int
    branch_index: int
    intersection_ids: list[int] = field(default_factory=list)

    # ── eigen / orbit metadata ──────────────────────────────────────────────

    @property
    def eigenvalue(self) -> Optional[float]:
        """
        Unstable eigenvalue magnitude (lambda_u) governing cdist scaling.

        Canonical distance scales by lambda_u per forward iterate on the unstable
        side and by 1/lambda_u on the stable side; a single magnitude is used for
        both (the map is area-preserving). Returns None if eigenvalues are unset.
        """
        evals = getattr(self.fixed_point, "unstable_eigenvalues", None)
        if not evals:
            return None
        return float(abs(np.asarray(evals[0]).ravel()[0]))

    @property
    def return_period(self) -> int:
        """
        Number of applications of the map M that send this branch onto itself.

        Equal to the fixed point's k_value (period, doubled under inversion).
        Falls back to the bare period if k_value has not been set.
        """
        return int(getattr(self.fixed_point, "k_value", self.fixed_point.period))

    @property
    def anchor_coords(self) -> NDArray[np.float64]:
        """Phase-space coordinates of the anchoring periodic point as (2,)."""
        return np.asarray(self.fixed_point.coordinates[self.orbit_index]).ravel()

    # ── ordering along the branch ───────────────────────────────────────────

    def ordered_ids(self, toward_anchor: bool = False) -> list[int]:
        """
        Return the intersection IDs in canonical-distance order.

        Args:
            toward_anchor: If False (default) order runs from the anchor outward
                (ascending cdist). If True, order runs toward the anchor
                (descending cdist).

        Returns:
            List of registry IDs.
        """
        return list(reversed(self.intersection_ids)) if toward_anchor else list(
            self.intersection_ids
        )

    def neighbor(self, intersection_id: int, toward_anchor: bool = True) -> Optional[int]:
        """
        Return the adjacent intersection on this branch.

        Walking "toward the anchor" moves to the next-lower canonical distance,
        matching the next-intersection primitive used by the Pseudoneighbor
        Algorithm.

        Args:
            intersection_id: Registry ID of the reference intersection.
            toward_anchor: Direction to step. True walks toward the anchoring
                periodic point (decreasing cdist); False walks outward.

        Returns:
            Registry ID of the neighbor, or None if there is none in that
            direction or intersection_id is not on this branch.
        """
        try:
            idx = self.intersection_ids.index(intersection_id)
        except ValueError:
            return None
        target = idx - 1 if toward_anchor else idx + 1
        if 0 <= target < len(self.intersection_ids):
            return self.intersection_ids[target]
        return None

    def __len__(self) -> int:
        return len(self.intersection_ids)

    def __repr__(self) -> str:
        return (
            f"TrellisBranch({self.stability}, orbit={self.orbit_index}, "
            f"branch={self.branch_index}, n_intersections={len(self.intersection_ids)})"
        )

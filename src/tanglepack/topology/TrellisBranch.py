"""
One branch of one manifold of one orbit point, as an ordered list of ids.

:class:`TrellisBranch` is the topological unit the algorithms walk: the
crossings on a single branch, ordered by canonical distance outward from the
anchoring periodic point, plus the orbit metadata they need.

Dev Notes:

A TrellisBranch is the topological view of one branch of one manifold of one
periodic-orbit point. It is intentionally lightweight: it stores the *ordering*
of intersections along the branch (encoded as the order of intersection_ids) and
the orbit metadata the algorithms need, but it does not hold a reference to
the registry. Anything that needs the underlying Intersection objects or their
cdists goes through the owning Trellis.

intersection_ids is always kept sorted by canonical distance ascending — i.e.
from the anchoring periodic point outward. "Toward the anchor" therefore means
walking toward index 0.

Return period: a branch maps onto itself under M^k_value
(k_value = q for a period-q orbit, 2q with inversion). The pseudoneighbor window
and lap use ``return_period = k_value``, and a return scales a canonical distance
by ``beta ** k_value`` (``FixedPoint.per_step_beta``); no branch reads an
eigenvalue.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal, Optional, TYPE_CHECKING


if TYPE_CHECKING:
    from ..numerics.FixedPoint import FixedPoint
    from ..numerics.Intersection import ManifoldKey


@dataclass
class TrellisBranch:
    """
    One branch of a stable or unstable manifold within a finite trellis.

    A branch corresponds exactly to one ``TangleWorkbench.manifolds`` key. It
    carries the intersections that lie on it, ordered by canonical distance from
    the anchoring periodic point outward, plus the orbit metadata the
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

    # ── orbit metadata ──────────────────────────────────────────────

    @property
    def return_period(self) -> int:
        """
        Number of applications of the map M that send this branch onto itself.

        Equal to the fixed point's k_value (period, doubled under inversion).
        Falls back to the bare period if k_value has not been set.
        """
        return int(getattr(self.fixed_point, "k_value", self.fixed_point.period))

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

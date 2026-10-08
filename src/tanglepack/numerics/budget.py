"""The point budget of one growth step.

Past the escape step a manifold arm's point count can jump by orders of
magnitude in ONE refinement step (the stress campaign's E3/E7 runs went from
~1e4 to ~1e6 points and then out of memory inside a single call), so the only
place that can stop the killer step is inside it.
:class:`~.ManifoldMachine.ManifoldMachine` counts the points each top-level
call creates and raises :class:`PointBudgetExceeded` BEFORE the map evaluation
that would take it over its budget.

The budget is a resource guard, not an invariant: it is an exception type,
never an ``assert``, and it says nothing about the tangle's topology.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

#: Points one top-level machine call may create (one growth step of a whole
#: ``k_value`` chain, or one bridge image). ~2 GB at the measured ~1.06 KB per
#: point; the largest fixture step chain (inversion at 1e-7) creates ~540k.
DEFAULT_MAX_POINTS_PER_STEP = 2_000_000


class PointBudgetExceeded(RuntimeError):
    """
    One growth step (or one bridge image) would create more points than allowed.

    Raised before the map evaluation that would cross the budget, so
    ``points_created <= budget`` always holds. The state left behind is NOT
    resumable: a growth chain may be partly mapped and under-refined, and a
    blasted bridge's points may already be linked to a half-refined image while
    the bridge is still un-iterated. Rebuild the session with fewer steps, a
    stop condition (``grow_until_*``), a coarser ``area_cutoff`` or a larger
    ``max_points_per_step`` (``None`` disables the budget).

    Attributes:
        budget: The per-call point budget that was exceeded.
        points_created: Points the call had created when it stopped.
        stage: Where it stopped: mapping the curve (``"map"``) or refining it
            (``"refine"``).
        max_abs_coord: Largest ``|coordinate|`` among the points created so far
            in the call -- a huge value means the arm is escaping to infinity.
        manifold_key: Key of the manifold being mapped, when it has one.
    """

    def __init__(
        self,
        budget: int,
        points_created: int,
        stage: Literal["map", "refine"],
        max_abs_coord: float,
        manifold_key: Optional[tuple] = None,
    ) -> None:
        """
        Record the overrun and build the message.

        Args:
            budget: The exceeded per-call budget.
            points_created: Points created before stopping.
            stage: ``"map"`` or ``"refine"``.
            max_abs_coord: Largest ``|coordinate|`` created so far.
            manifold_key: Key of the manifold being mapped, if any.
        """
        self.budget = budget
        self.points_created = points_created
        self.stage = stage
        self.max_abs_coord = max_abs_coord
        self.manifold_key = manifold_key
        super().__init__(
            f"one growth step or bridge image would create more than {budget:,} "
            f"points (stopped while {'mapping' if stage == 'map' else 'refining'} "
            f"after {points_created:,}; largest |coordinate| {max_abs_coord:.3g}). "
            "A huge coordinate means an escaping arm: grow with a stop condition "
            "(grow_until_bounded and the other grow_until_* drivers); otherwise "
            "use fewer steps or blasts, a coarser area_cutoff, or a larger "
            "max_points_per_step (None disables it). The session is not "
            "resumable; rebuild it."
        )


@dataclass(frozen=True)
class GrowthStats:
    """
    What the last growth step of one manifold cost (``workbench.growth_stats``).

    Attributes:
        points_created: Points the step created over its whole ``k_value`` chain.
        growth_factor: ``points_created`` over the previous step's, None after
            the first step. Its product with ``points_created`` is the hint
            ``grow_n_times`` compares with the budget -- an UNDERestimate past
            escape, where the factor itself jumps 5-8x.
        max_abs_coord: Largest ``|coordinate|`` the step mapped or created.
    """

    points_created: int
    growth_factor: Optional[float]
    max_abs_coord: float

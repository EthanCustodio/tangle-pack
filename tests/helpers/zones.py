"""Resonance zones for the loom and facade tests.

The one recipe for "a zone whose trim really shortens the stable branch",
shared by the blast / restore tests and the session cache table (it replaces
the two local ``_define_zone`` copies, 2026-10-05 refactor Phase 7b).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from tanglepack import TangleSession
    from tanglepack.loom.ResonanceZone import ResonanceZone
    from tanglepack.numerics.FixedPoint import FixedPoint


def define_inner_zone(session: "TangleSession", fixed_point: "FixedPoint") -> "ResonanceZone":
    """
    Trim at the outermost strong-pip candidate strictly inside the stable extent.

    The default strong pip of the k=10 flow sits at the largest stable
    canonical distance, so trimming there is a geometric no-op. The outermost
    candidate strictly inside the current extent is chosen instead, so the
    trim removes crossings and a later restore has something to rebuild.

    Args:
        session: A session at (or past) the bridges stage.
        fixed_point: The fixed point whose zone to define.

    Returns:
        The zone the session registered for ``(fixed_point, 0)``.

    Raises:
        AssertionError: No candidate lies inside the stable extent.
    """
    registry = session.workbench.intersection_registry
    trellis = session.trellis(fixed_point)
    trellis.classify_strong_pips()
    outermost = max(registry[i].stable_cdist for i in registry.all_ids())
    inner = [c for c in trellis.strong_pip_candidates if registry[c].stable_cdist < outermost]
    assert inner, "expected a strong-pip candidate inside the stable extent"
    pip = max(inner, key=lambda c: registry[c].stable_cdist)
    trellis.set_strong_pip(pip)
    session.add_resonance_zones([pip])
    return session.resonance_zones[(fixed_point, 0)]

"""Law tier: crossings.

Only an unstable manifold can cross a stable one (CLAUDE.md, "Fundamental
invariant"): every registered crossing is one unstable key times one stable
key, and no two unstable (or two stable) curves cross geometrically, within
one branch, across branches or across fixed points. A crossing's cdists are
defined and bracketed by the segment it was found on, its sign is the cross
product of the two increasing-cdist directions, and every registered ``+1``
link preserves the canonical area and scales the cdists by
``per_step_beta``. Recomputing the crossings keeps ids and is idempotent
(the ``recompute`` layer: one fresh build per case, since these two laws
mutate the session).
"""

from __future__ import annotations

from helpers.law_tier import layer_test

test_crossings_laws = layer_test("crossings")
test_recompute_laws = layer_test("recompute")

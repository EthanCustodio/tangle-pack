"""Law tier: crossings.

Only an unstable manifold can cross a stable one (CLAUDE.md, "Fundamental
invariant"): every registered crossing is one unstable key times one stable
key, and no two unstable (or two stable) curves cross geometrically, within
one branch, across branches or across fixed points. A crossing's cdists are
defined and bracketed by the segment it was found on, its sign is the cross
product of the two increasing-cdist directions, and every registered ``+1``
link scales each cdist by the ``per_step_beta`` of its own saddle (so a
homoclinic link preserves the canonical area and a heteroclinic one scales it
by ``beta_A / beta_B``, NUM-20). Recomputing the crossings keeps the topology (D1, NUM-12),
keeps every id and is idempotent (the ``recompute`` layer: one fresh build per
case, since these laws mutate the session).
"""

from __future__ import annotations

from helpers.law_tier import layer_tests

test_crossings_laws, test_crossings_set_aside = layer_tests("crossings")
test_recompute_laws, test_recompute_set_aside = layer_tests("recompute")

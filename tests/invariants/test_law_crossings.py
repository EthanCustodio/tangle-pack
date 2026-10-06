"""Law tier: crossings.

Only an unstable manifold can cross a stable one (CLAUDE.md, "Fundamental
invariant"): every registered crossing is one unstable key times one stable
key, and no two unstable (or two stable) curves cross geometrically, within
one branch, across branches or across fixed points. A crossing's cdists are
defined and bracketed by the segment it was found on, its sign is the cross
product of the two increasing-cdist directions, and every registered ``+1``
link preserves the canonical area and scales the cdists by
``per_step_beta``. Recomputing the crossings keeps ids and is idempotent
(fresh builds: these two laws mutate the session).
"""

from __future__ import annotations

from helpers import laws
from helpers.law_tier import law_test

test_crossings_unstable_by_stable = law_test(laws.check_crossings_unstable_by_stable)
test_no_same_stability_crossing = law_test(laws.check_no_same_stability_crossing)
test_crossings_cdists_defined = law_test(laws.check_crossings_cdists_defined)
test_crossing_cdist_bracketed = law_test(laws.check_crossing_cdist_bracketed)
test_crossing_sign_is_cross_product = law_test(laws.check_crossing_sign_is_cross_product)
test_area_along_iterate_links = law_test(laws.check_area_along_iterate_links)
test_iterate_link_scales_by_beta = law_test(laws.check_iterate_link_scales_by_beta)
test_recompute_preserves_ids = law_test(
    laws.check_recompute_preserves_ids, fixture="fresh_law_case"
)
test_recompute_is_idempotent = law_test(
    laws.check_recompute_is_idempotent, fixture="fresh_law_case"
)

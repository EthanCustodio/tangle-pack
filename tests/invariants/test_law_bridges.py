"""Law tier: bridges and their identity.

A bridge is partial exactly when it has no ``BridgeId``; its id is its two
crossings in increasing unstable cdist, both on the bridge's OWN unstable
branch (the period-3 regression of the 2026-09 codebase audit) and inside its
span; ``bridges_at`` indexes exactly the two ends; each id has one registered
copy; two bridges of one branch only touch; and the derived image genealogy
lands on the advanced branch and round-trips through the preimage.
"""

from __future__ import annotations

from helpers import laws
from helpers.law_tier import law_test

test_bridge_partial_iff_no_id = law_test(laws.check_bridge_partial_iff_no_id)
test_bridge_id_in_unstable_order = law_test(laws.check_bridge_id_in_unstable_order)
test_bridge_endpoints_on_own_branch = law_test(laws.check_bridge_endpoints_on_own_branch)
test_bridges_at_exact = law_test(laws.check_bridges_at_exact)
test_bridge_single_copy = law_test(laws.check_bridge_single_copy)
test_bridges_do_not_overlap = law_test(laws.check_bridges_do_not_overlap)
test_bridge_image_round_trip = law_test(laws.check_bridge_image_round_trip)

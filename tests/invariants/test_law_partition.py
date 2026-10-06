"""Law tier: holes and the stable (homotopy) partition.

The partition covers each branch on each side in contiguous positional
elements, every crossing has exactly one owner per side, pinched singletons
own exactly their point. Every hole is classified; the holes of one origin
share their bridge side (I1); a direct hole's side is the side of its bridge
its coordinates sit on. The OPENINGS LAW (Phase 0 probe P1): a direct hole
opens its inward pair on its own bridge's row; every opening sits on its
bridge's row at that bound (a backward hole's row may differ from its
origin's when its image lobe's stable side was never computed); at a bound
that IS the registered iterate of the origin's bound the opening is the
origin's; a bound opens nothing only at an anchor or the branch end.
Propagated holes land on the branch the cycle predicts and stop at
periodicity; nothing at iterate ``>= k_value`` is punched (the firm half of
holes-backward-only; the ``+1 .. +(k-1)`` exemption is provisional and
untested); reference pairs are consecutive on both branches inside the
reference window.
"""

from __future__ import annotations

from helpers import laws
from helpers.law_tier import law_test

test_partition_covers_branch = law_test(laws.check_partition_covers_branch)
test_partition_unique_owner = law_test(laws.check_partition_unique_owner)
test_partition_singletons = law_test(laws.check_partition_singletons)
test_holes_are_classified = law_test(laws.check_holes_are_classified)
test_holes_share_bridge_side = law_test(laws.check_holes_share_bridge_side)
test_direct_hole_side_is_coordinate_side = law_test(
    laws.check_direct_hole_side_is_coordinate_side
)
test_direct_hole_opens_inward_pair = law_test(laws.check_direct_hole_opens_inward_pair)
test_openings_on_own_bridge_row = law_test(laws.check_openings_on_own_bridge_row)
test_openings_linked_bound = law_test(laws.check_openings_linked_bound)
test_openings_missing_only_at_anchor_or_tail = law_test(
    laws.check_openings_missing_only_at_anchor_or_tail
)
test_propagated_holes_land_on_predicted_branch = law_test(
    laws.check_propagated_holes_land_on_predicted_branch
)
test_propagation_terminates = law_test(laws.check_propagation_terminates)
test_no_direct_hole_beyond_fundamental = law_test(laws.check_no_direct_hole_beyond_fundamental)
test_reference_pairs_valid = law_test(laws.check_reference_pairs_valid)

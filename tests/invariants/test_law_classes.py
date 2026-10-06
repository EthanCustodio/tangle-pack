"""Law tier: bridge classes.

Every identified bridge lands in exactly one class; a class is oriented
anchor outward (``anchor_outward_key``) and a member's direction is ``+1``
exactly when its first end is the source; the combinatorial row of a bridge
end (crossing sign alone) is the side the bridge geometrically approaches
from; a same-branch bridge has one row at both ends (I2); a tangle's anchor
bridge is a ``+1`` member of an active class that leads its group; the table groups
by tangle (connecting classes last) and orders by smallest member cdist; a
class mixes fixed points only as a connecting class; the classes of a fixed
point reach every stable branch of its cycle; only active classes are
lettered.
"""

from __future__ import annotations

from helpers import laws
from helpers.law_tier import law_test

test_every_bridge_in_one_class = law_test(laws.check_every_bridge_in_one_class)
test_class_orientation_anchor_outward = law_test(laws.check_class_orientation_anchor_outward)
test_row_of_end_is_geometry = law_test(laws.check_row_of_end_is_geometry)
test_bridge_rows_consistent = law_test(laws.check_bridge_rows_consistent)
test_anchor_bridge_class_leads_its_tangle = law_test(
    laws.check_anchor_bridge_class_leads_its_tangle
)
test_class_table_order = law_test(laws.check_class_table_order)
test_classes_do_not_mix_tangles = law_test(laws.check_classes_do_not_mix_tangles)
test_classes_use_every_orbit_branch = law_test(laws.check_classes_use_every_orbit_branch)
test_only_active_classes_lettered = law_test(laws.check_only_active_classes_lettered)

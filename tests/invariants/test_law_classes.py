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

from helpers.law_tier import known_issue_test, layer_test

test_classes_laws = layer_test("classes")
test_classes_known_issue = known_issue_test("classes")

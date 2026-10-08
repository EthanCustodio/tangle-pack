"""Law tier: bridge classes.

The classing laws read the EVIDENCE table (``minimal_trellis().classes``, every
bridge's class); the lettering laws read the ALPHABET (``bridge_classes()``,
the classes with a member the minimal trellis keeps, D2). Every identified
bridge lands in exactly one class or is an orphan loop left unclassed; a
class is never a loop and is oriented anchor outward (``anchor_outward_key``),
a member's direction is ``+1`` exactly when its first end is the source; a
loop folds FORWARD up its image-chain parents into a non-loop member of its
class (D2); the combinatorial row of a bridge
end (crossing sign alone) is the side the bridge geometrically approaches
from; a same-branch bridge has one row at both ends (I2); a tangle's anchor
bridge is a ``+1`` member of an active class that leads its group (in a tangle
whose every stable branch carries a hole, D9); the table groups
by tangle (connecting classes last) and orders by smallest member cdist; a
class mixes fixed points only as a connecting class; the classes of a fixed
point reach every stable branch of its cycle; only active classes are
lettered; and the alphabet is exactly the evidence classes with a member the
minimal trellis keeps, the others never lettered (D2).
"""

from __future__ import annotations

from helpers.law_tier import layer_tests

test_classes_laws, test_classes_set_aside = layer_tests("classes")

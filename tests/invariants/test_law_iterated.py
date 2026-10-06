"""Law tier: the iterated homotopy partition, element names and the minimal trellis.

An iterated element lies inside its homotopy parent; every homotopy boundary
survives the cut with its closedness and every new boundary is an image
bridge endpoint on that row; every crossing (and every stable arc midpoint of
the minimal arrangement) has one owner per side; an element under a
same-branch image bridge records the bridge covering it. Names agree with the
structure (side letter, parent subscript, round trips). The minimal trellis
keeps the hole bridges and the images of the active ones and nothing else,
keeps the right nodes, and its sparse arrangement carries exactly the kept
bridges and closes up.
"""

from __future__ import annotations

from helpers.law_tier import layer_test

test_iterated_laws = layer_test("iterated")

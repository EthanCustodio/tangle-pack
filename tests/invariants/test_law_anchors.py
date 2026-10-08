"""Law tier: the anchors (the periodic points registered as crossings at cdist (0, 0)).

Each periodic point is ONE anchor, the first crossing of every branch attached
to it (D4: two stable and two unstable branches at an inversion point; this
replaced the 2026-10-05 "one anchor per unstable branch" law in WP3), the
arrangement's real arcs leaving it are exactly those branches (one degree-4
node, never four coincident ones), and an anchor's crossing sign is the
handedness of its home pair's oriented eigendirections.
"""

from __future__ import annotations

from helpers.law_tier import layer_tests

test_anchors_laws, test_anchors_set_aside = layer_tests("anchors")

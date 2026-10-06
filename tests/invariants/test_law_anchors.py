"""Law tier: the anchors (the periodic points registered as crossings at cdist (0, 0)).

Every unstable branch carries exactly one anchor (author, 2026-10-05; the
inversion saddle's two anchors per branch are a known issue), and an
anchor's crossing sign is the handedness of the two oriented eigendirections.
"""

from __future__ import annotations

from helpers.law_tier import known_issue_test, layer_test

test_anchors_laws = layer_test("anchors")
test_anchors_known_issue = known_issue_test("anchors")

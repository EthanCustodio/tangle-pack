"""Law tier: the anchors (the periodic points registered as crossings at cdist (0, 0)).

Every unstable branch carries exactly one anchor (author, 2026-10-05; the
inversion saddle's two anchors per branch are a known issue), and an
anchor's crossing sign is the handedness of the two oriented eigendirections.
"""

from __future__ import annotations

from helpers import laws
from helpers.law_tier import law_test

test_one_anchor_per_unstable_branch = law_test(laws.check_one_anchor_per_unstable_branch)
test_anchor_sign_matches_eigendirections = law_test(
    laws.check_anchor_sign_matches_eigendirections
)

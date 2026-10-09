# Decisions

What the author has decided. Binding on every later run. Newest last.

## 2026-10-02 — blasting nested tangles (slug `hole-side-and-own-blast`, base `2315204`)

- **Said** (plan verdict, `plan_ok_planner_human.json`): approved with no notes; the plan's
  defaults stand.
- **Consequence 1 — direct-hole side is combinatorial.** A directly punched hole's
  `Hole.bridge_side` comes from the crossing signs, in the style of `row_of_end`:
  `first` = lower unstable cdist end, `outward = +1 if second.stable_cdist > first.stable_cdist
  else -1`, `left iff first.crossing_sign * outward > 0`. Not the lobe-polygon area, not a
  nearest-vertex measurement.
- **Consequence 2 — a propagated hole's side is read on its image sub-arc** (the `span` inside
  the containing bridge), never against the whole bridge polyline.
- **Consequence 3 — each nested zone blasts only its own fixed point**
  (`fixed_point=[zone.fixed_point]` in `build_nested`); blast order no longer matters.
- **Consequence 4 — CLAUDE.md** carries one sentence: the period-3 side trip at 15-16 steps and
  6+ blasts is fixed, but those runs are unreliable (an unreachable class, a virtual `new1`;
  open). No other CLAUDE.md change.
- Left open by the plan's "Not doing" list (not decided, not to be assumed): p3 >= 6 blasts
  unreliable; `ManifoldMachine.iterate_manifold -> merge_manifolds` non-contiguous crash;
  `_hole_openings` whole-polyline nearest vertex; pip choice lost on generation bump (`_repin`);
  `Blast.py` tidy-up; constant nested slug; noisy inert-word WARNING; blast speed.

## 2026-10-02 — outcome verdict on `hole-side-and-own-blast` (first outcome)

- **Said** (`debrief_ok_auditor_human.json`): "figures", no further notes.
- **Consequence:** proof figures were made on the branch (commit `f8a6d73`, three scripts under
  `scripts/refactor_proof/hole-side-and-own-blast/`; PNGs untracked under
  `figures/refactor_proof/hole-side-and-own-blast/` in the branch worktree) and the auditor
  re-audited (`outcome_auditor_human_2.md`, second pass). The proof scripts are outside the
  plan's five-file list by instruction; rubric 5 is judged against `3b33baa`.

## 2026-10-07 — open-issues round (catalog `docs/open_issues.md`, base `684d70b`)

- **Said** (planning session answers, 2026-10-07), recorded as D1-D8 in `docs/open_issues.md`:
  - **D1 full recompute.** "We only start doing the topology after all the numerics is done …
    any growth that happens requires entirely recomputing the topology." Planned loop: blast,
    compute the topology, blast again, recompute, look for topological changes. The class
    count "may increase multiple times, strictly increasing, never decreasing. Linear growth
    over many iterates is a clear bug though."
  - **D2** the alphabet is the minimal trellis's bridges only.
  - **D3** "the shortest walk on the dual graph gives the word, that is the algorithm" and a
    walk must always exist; the virtual `new1` class is "a very clear bug"; f(bridge) makes the
    word for singleton bridges that cannot walk. The period-3 path "seems buggy".
  - **D4** "Each anchor point should have two stable branches and two unstable branches …
    Each of those individual branches of course have a left side and right side." The left of
    one inversion branch is the right of the other; zones break the space up and must be kept
    in mind; cross-side is a flag only for the SAME branch, never for a bridge between two
    branches.
  - **D5** keep the period-k +1..+(k-1) iterate rule, fix the gap where table links are missing.
  - **D6** signs decide hole openings and rows; geometry only checks.
  - **D7** heteroclinic pairs give holes in B's partition.
  - **D8** simplify freely ("make as many choices as you'd like to make the design simpler …
    robust and professional"); public-API alignment is a later, separate plan.
- **Consequence:** WP0-WP12 in `docs/open_issues.md`, each run analyst -> planner ->
  test-engineer -> coder -> code-reviewer on its own branch. Claude's calls C1-C6 under D8
  (anchor-only branch = one element, element identity by bounding holes, det<0 out of scope,
  innermost-zone blasting, new-image-only same-stability count, per-call `max_points`) are
  revisable by the author.

## 2026-10-08 — DEC-19, the period-k growth stagger (raised by WP5)

- **Said:** "Accept the asymmetry" (offered: trim the trellis to complete orbits; warn and grow more; accept).
- **Consequence (D9):** a stable branch whose reference pair's +j images lie beyond the computed manifold has no fundamental-segment hole; that is correct for the finite trellis, not a bug. WP5's WARNING naming the branch stays. Laws that assume every branch has a window hole (e.g. the inversion anchor class leading its tangle) must be scoped to branches that have one, not left as known issues.

## 2026-10-08 — DEC-20, DEC-21, DEC-22 (raised by WP7 and WP8)

- **DEC-20, said:** "Cut S at the pip orbit." **Consequence (D10):** the partition and the dual graph use the stable manifold cut at the strong pip's orbit (chalkboard stage 3, as pseudoneighbors already do); this replaces the 2026-09-14 "the pip moves the unified set without trimming" rule.
- **DEC-21, said:** "I'm not sure this is an issue? We do not need computed bridges to find words, we strictly use the dual graph unless it is a singleton. That allows us to potentially find words even though we haven't computed the actual iterates. A word can include the class itself and any other classes." Confirmed on follow-up: "Yes, but mark them." **Consequence (D11):** a walk pair no computed bridge spans is a legitimate class: lettered, walked from its elements' landings, NOT a reliability failure, and marked "no computed bridge yet" (informational). Supersedes C9's "a pair no bridge spans is never a class".
- **DEC-22, said:** "Balanced trellis." **Consequence (D12):** cut the crossing set at an unstable depth matched to the stable cut at the pip, so blasting cannot reshape the fundamental segment at depth; the class count then never drops under blasting.

## 2026-10-08 — WP10 (numerics and performance leftovers, slice 1), Claude's calls under D3 / D8

Revisable by the author; recorded so later runs know they are calls, not the author's words.

- **C11 (under D3) — incomplete evidence gives no verdict.** A registered image chain whose
  word parts from the walk because the image passes the end of the computed stable branch
  (diagnosis `"stable end"`: two consecutive chain crossings on one stable branch with one
  sign) misses its crossings with the trimmed-away part, so it is incomplete: `word_agrees`
  is None, the diagnosis is kept and logged at INFO, and it never makes `verified` False or
  the dynamics unreliable. A `"walk"` disagreement stays a WARNING and a reliability failure.
  Effect: k28 at 6-7 blasts is reliable; at 8+ blasts the DEC-23 ties keep it unreliable.
- **C12 (under D8) — the law tier builds each case once per RUN.** `law_case` is session
  scoped (was module scoped: every one of ~12 law modules rebuilt all 11 cases, 40 s of
  setup); the fingerprint guard still fails a test that mutates the shared build.

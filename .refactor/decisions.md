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

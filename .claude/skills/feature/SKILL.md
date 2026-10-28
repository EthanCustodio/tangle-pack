---
name: feature
description: Add new behaviour to tanglepack on a worktree branch, with tests for the new behaviour, plan approval up front and a reviewed diff at the end.
argument-hint: <the feature>
disable-model-invocation: true
---

# /feature $ARGUMENTS

Same flow as `/refactor` (read `.claude/skills/refactor/SKILL.md` and follow its steps), with
these differences:

- **Branch** `feature/<slug>`, worktree `../tangle-pack-<slug>`.
- **Plan** additionally states the behaviour being added, where it sits in the pipeline (the stage
  in CLAUDE.md), the API surface (names, signatures), and which tests will show it works: a law
  in `tests/invariants/` when the feature is a physical or combinatorial law, otherwise
  `tests/unit/...`. The test policy in `tests/conftest.py` holds: no pinned case-specific word,
  letter, matrix or count.
- **Workflow** `tp-implement` with `mode: "feature"`: the implementer writes those tests too.
- **Report** also lists what changed in any output and why, and the CLAUDE.md / design-note lines
  that need updating (the implementer updates them on the branch when the plan says so).

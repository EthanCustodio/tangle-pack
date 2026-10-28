---
name: house-style
description: The author's code style for tanglepack (compact, human, pythonic, no needless machinery) and the branch convention every workflow follows.
user-invocable: false
---

# House style

The author wants code a person enjoys reading. Good, working, testable code, and no more of it
than the problem needs. CLAUDE.md's Coding Style section (Google docstrings, type hints, logging,
Dev Notes) still applies; this is the taste on top.

**Write**

- Plain functions and small classes. Reach for a class only when there is real state to hold.
- Comprehensions, unpacking, `enumerate`, `zip`, dataclasses, early returns.
- A docstring where the *why* is not obvious. Comments explain intent, never restate code.
- Keep the module's existing idiom. Read the surrounding file before adding to it.

**Do not write**

- Helpers called once, wrapper layers, abstract base classes with one subclass, config objects
  for two parameters, or "manager"/"handler" classes.
- Defensive code for states the invariants rule out. CLAUDE.md's invariants are physics: two
  unstable manifolds never intersect, maps preserve area, holes propagate backward only. Code may
  *assert* them; it must never model their violation.
- Speculative generality, feature flags, or backwards-compatibility shims nobody asked for.
- Silent behaviour changes. A refactor that changes an output says so.

**Humanising** means a reader who knows the maths can follow the code without a debugger:

- One idea per function, in the order the mathematics happens; the top of a module reads as its
  outline.
- Names from the paper, not from the plumbing: `unstable_cdist`, not `val2`.
- Delete narration comments and commented-out code. Keep the comment that says *why*: an author's
  rule, an invariant, the date it was decided.
- Collapse indirection a reader has to chase: a wrapper that only forwards, a flag one caller sets.

**The test of a change:** it is shorter or clearer than what it replaced, the suite passes, and
someone who knows the maths can read it top to bottom.

# Branches and running things

- Python: `MPLBACKEND=Agg env/bin/python`; tests: `MPLBACKEND=Agg env/bin/python -m pytest -q`.
  In a worktree use the main repo's `env/bin/python`: pytest imports the worktree's own `src`
  (`pythonpath` in pyproject), a script needs `PYTHONPATH=src` or it imports the main repo's code.
- Library changes (`/refactor`, `/feature`, `/perf`) live on `<kind>/<slug>` in the worktree
  `../tangle-pack-<slug>`, made from the author's `HEAD`. Never commit to the author's branch from
  a worktree job and never touch their uncommitted changes.
- Scripts, figures and quick edits (`/experiment`, `/figure`, `/quick`) stay on the current branch
  and are committed only when the author asks.
- The author's binding decisions are in `.refactor/decisions.md`. Read the entries that touch what
  you change; append one (dated, one paragraph) when the author states a new rule.

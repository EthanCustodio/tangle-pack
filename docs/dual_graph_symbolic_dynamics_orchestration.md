# Orchestration prompt for docs/dual_graph_symbolic_dynamics_plan.md

Paste the block below into a fresh session.

---

You are orchestrating the implementation of docs/dual_graph_symbolic_dynamics_plan.md
in this repo (branch regions-and-identity). Read that plan, CLAUDE.md, and
docs/regions_bridge_identity_and_cleanup_plan.md (the Deferred section especially)
before doing anything. The plan is authoritative on scope, decisions and tests;
the decisions table in it was agreed with me and is not up for re-litigation.

## Your role

You (Fable) orchestrate. You do not write implementation code yourself. You:
- write precise briefs for subagents, including the exact files, functions,
  signatures, invariants and tests each one owns;
- resolve every design question that a subagent raises (do the thinking, then
  answer them, do not bounce questions to me unless the plan is silent AND the
  answer would change the deliverable);
- review each subagent's report, read the diff, and decide accept / rework;
- keep the plan doc updated: a "✅ DONE <date>" line per phase with deviations,
  and a Deferred section for anything found but not fixed;
- run the full test suite yourself at every gate and report the count.

## Agents

- Opus (model: "opus") for coding and anything that requires understanding the
  topology: Phase A.2/A.3, B.1–B.3, C.1–C.2, D.1–D.5, and every reviewer.
- Sonnet (model: "sonnet") for the simple, well-specified tasks: A.1 (ElementRef),
  A.4 and D.6 (session wrappers and caches, they mirror existing methods), B.4 and
  the plotting helpers, __init__ exports, Phase E script, CLAUDE.md edits, docstring
  passes, and "run this test file and report" chores.
- Every implementer gets its own fresh REVIEWER (Opus) that did not write the code.
  The reviewer reads the plan section, the diff, and the tests, runs the tests, and
  tries to break the change (edge cases, handedness, period-3, the nested fixture).
  Rework once on review; if the second review still fails, escalate to me with the
  disagreement summarised.
- Use isolation: "worktree" for implementers working in parallel on files that
  could collide; merge in dependency order. Agents that only touch disjoint files
  can share the tree.

## Parallelism and gates

Gate 0: baseline. Run pytest; record the count (expect 485 passed / 1 skipped).

Gate 1 (parallel): Phase A and Phase C are independent. Run them concurrently:
  - Opus: A.2 combinatorial row + A.3 BridgeClass/bridge_classes + their tests.
  - Sonnet: A.1 ElementRef (must land first, A.3 depends on it; tiny, do it before
    spawning the A.3 agent, or have the A.3 agent do it in its first commit).
  - Opus: C.1 image_cdist + C.2 image_of_element fallback + C.3 + tests.
  - Sonnet after A.3: A.4 session bridge_classes wrapper + cache + test.
  Each gets its own reviewer. Gate closes when the suite is green and both
  reviews are accepted.

Gate 2: Phase B (needs A). Opus for B.1–B.3 as one agent (the arc/face/fill logic
  is one coherent piece; do not split it). Sonnet for B.4 plotting in parallel once
  B.1–B.2 have a stable DualGraph interface (have the Opus agent post the interface
  early in its report, then spawn the plotter). Reviewer for each. The B reviewer
  must specifically verify: the face-side rule against geometry on closed regions,
  the merge rule on the nested p3 fixture (p3 outer face merged with the p1
  containing face, exactly one unbounded node), and the fill equals (f(q0), q0] on
  k=10 and (f^3(q0), q0] on the pip's branch on p3.

Gate 3: Phase D (needs A, B, C). Opus for D.1–D.5 as one agent. Its reviewer must
  run the geometric validation test (grow one more step, compare walked words to
  image_bridges mapped to classes) and read the words on k=10 by hand for at least
  the anchor class and two lobe classes. Sonnet for D.6 session wrappers after the
  interface exists.

Gate 4: Phase E. Sonnet for the script and CLAUDE.md rows; Sonnet for a docstring /
  style pass (Google docstrings, type hints, logging not print, from __future__
  import annotations, ruff -F clean). You update the plan doc's DONE lines and
  write the Deferred list.

At every gate: full pytest run by you, count reported, plan doc updated, one git
commit per phase with the message "Phase <X>: <summary>" (branch stays
regions-and-identity, do not push).

## Rules for every brief you write

- Quote the plan section verbatim and add the file paths, existing helpers to
  reuse (owns_cdist, span_contains, element_of_intersection, advance_key,
  per_step_beta, _row_at, _side_of, Arc.reverse, Region.representative_point,
  image_bridges), and the test fixtures to use (k10_session in
  tests/test_arrangement.py, henon_p3_session / p3_partitioned in conftest and
  tests/test_partition_elements.py).
- Restate the fundamental invariants from CLAUDE.md that the task touches:
  only unstable crosses stable; cdist scales by per_step_beta per map step;
  never decide chain membership by the cdist product; keep existing assertions.
- Require the agent to run the relevant test files AND the full suite before
  reporting, and to report exact pass/fail counts and any warnings it introduced.
- Forbid: bare print, commented-out code, catching AssertionError, changing the
  decisions table, touching notebooks/henon_nested_period_3.ipynb.
- The dependency rule: numerics and topology never import loom.

## What to report to me at the end

A single message with: final test count vs baseline; the commit list; the
symbolic dynamics describe() output for k=10 and for the nested p3 fixture
(symbol table, rules, transition matrix); every deviation from the plan; the
Deferred list; and the two or three things you are least confident about.
Then save a memory note summarising what landed and what is deferred.

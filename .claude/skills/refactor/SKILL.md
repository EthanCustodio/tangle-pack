---
name: refactor
description: Restructure tanglepack code without changing what it computes, on a worktree branch, with plan approval up front and a reviewed diff at the end.
argument-hint: <what to refactor and why>
disable-model-invocation: true
---

# /refactor $ARGUMENTS

A refactor changes how the code reads, never what it computes. If an output must change, stop
and say so: that is a `/feature`. Follow the `house-style` skill.

1. **Scout.** Read the code in question and the `.refactor/decisions.md` entries that touch it.
   If a `/campaign` file names this run, read its section. Find the tests that cover it.
2. **Plan.** Enter plan mode. The plan states: the slug, the goal in one sentence, the file list
   (the scope; nothing outside it changes), the steps, and *what must not change* (outputs,
   public API, numbers). Ask at most a couple of questions with `AskUserQuestion`, only for real
   choices. `ExitPlanMode` is the first checkpoint.
3. **Set up.**
   `git worktree add ../tangle-pack-<slug> -b refactor/<slug> HEAD`, record the base sha, write the
   approved plan to `../tangle-pack-<slug>/.plan.md`, and record the baseline in that worktree
   BEFORE anything is edited there (cd into it; `env/bin/python` is the main repo's):
   `MPLBACKEND=Agg <repo>/env/bin/python -m pytest -q -p no:cacheprovider --cov=tanglepack --cov-report=json:.base-cov.json`.
   Keep the pass/xfail counts. Both files stay untracked (add them to the worktree's
   `.git/info/exclude` or simply never `git add` them).
4. **Implement.** Run the saved workflow `tp-implement` with
   `args: {mode: "refactor", slug, repo: "<abs main repo>", worktree: "<abs path>", plan: "<abs>/.plan.md", base: "<sha>", base_cov: "<abs>/.base-cov.json", base_suite: "<counts>"}`.
   A characterization test is written only when the touched code has no coverage.
5. **Report** (the second checkpoint). Show: `git diff --stat base..`, net line count, suite
   before/after, the coverage guard result, the confirmed review findings and what was done with
   them, and any leftovers. Ask to merge.
6. **Merge.** On approval: `git merge --ff-only refactor/<slug>` on the author's branch (rebase
   the branch first if `HEAD` moved), `git worktree remove ../tangle-pack-<slug>`, delete the branch.
   If the author stated a rule during the run, append it to `.refactor/decisions.md`. If a
   campaign file names this run, tick it off.

---
name: review
description: Multi-lens review of a tanglepack branch or commit range (physics invariants, correctness, house style, tests), each finding adversarially verified.
argument-hint: "[branch | base..head | worktree path]  (default: current branch vs main, else the working tree vs HEAD)"
disable-model-invocation: true
---

# /review $ARGUMENTS

1. Resolve the target to a `base` sha and a `worktree` path (the repo itself unless a worktree is
   named). With no argument: the merge-base of `HEAD` and `main` if they differ, else `HEAD`
   (reviewing the uncommitted working tree). Say in one line what is being reviewed and its
   `git diff --stat`.
2. Run the saved workflow `tp-review` with `args: {base, worktree}`.
3. Report the confirmed findings with `ReportFindings`, most severe first. Fix nothing unless the
   author asks.

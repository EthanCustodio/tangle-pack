---
name: perf
description: Profile a slow or memory-heavy path in tanglepack, try several optimizations in parallel worktrees, measure each, and merge the winner.
argument-hint: <what is slow, and the benchmark if known>
disable-model-invocation: true
---

# /perf $ARGUMENTS

Optimizations never change what is computed: results must be equal (byte-equal registries where
the code promises it, the suite green). Follow the `house-style` skill.

1. **Scout.** Pick or write the benchmark: `stress/profile/profile_numerics.py` (see
   `stress/profile/README.md`) or a small timing script in the scratchpad that runs the slow path
   on one case with an explicit `area_cutoff`. Read the hot code.
2. **Plan.** Enter plan mode with the benchmark command, the current numbers if already known,
   and 2-4 candidate optimizations (one sentence each, the files each touches). `ExitPlanMode`.
3. **Measure.** Run the saved workflow `tp-perf` with
   `args: {target, benchmark: "<command>", candidates: [{name, idea, files}], base: "<sha>", repo: "<abs main repo>"}`.
   Each candidate runs in its own worktree and comes back with timing, memory and suite status.
4. **Choose.** Show the table. Ask the author which to keep (none is an answer).
5. **Merge.** The chosen candidate's worktree branch goes through steps 5-6 of `/refactor`
   (report, approval, ff merge); delete the other branches and worktrees. Record the
   before/after numbers in the commit message.

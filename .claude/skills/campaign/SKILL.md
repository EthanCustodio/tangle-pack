---
name: campaign
description: Plan a larger tanglepack effort as a short sequence of runs (at most 4), each a /refactor, /feature or /perf call, tracked in one campaign file.
argument-hint: <the goal of the campaign>
disable-model-invocation: true
---

# /campaign $ARGUMENTS

The author's rule: a campaign is **at most 4 runs, fewer is better**, and this first call only
plans it. Each run is then its own `/refactor`, `/feature` or `/perf` call, typed by the author.

1. **Scout** the goal: read the code involved, measure if the goal is about speed or memory, and
   read the relevant `.refactor/decisions.md` entries. Use `/investigate`-style probes as needed.
2. **Plan** in plan mode: the goal, why it needs more than one run, and the runs in order. Each
   run gets a slug, its skill, a goal sentence, the files it touches, and how its success is
   measured. Merge runs that could be one. `ExitPlanMode`.
3. **Write** `.refactor/campaigns/<YYYY-MM-DD>-<slug>.md`: a short header (goal, date, base sha),
   then one `## Run N: <slug> (/skill)` section per run with a `- [ ] done` line. This file is the
   only campaign state.
4. **Tell** the author the command for run 1, e.g. `/refactor run 1 of .refactor/campaigns/<file>`.

Each run's skill reads its section, and on merge ticks `- [x] done (<merge sha>, <one-line
result>)`. When the last run merges, add a one-paragraph outcome at the top of the file.

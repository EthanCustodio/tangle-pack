---
name: investigate
description: Answer one question about tanglepack (why something happens, whether a claim holds, where a bug comes from) read-only, with file:line evidence and measured numbers.
argument-hint: <the question>
disable-model-invocation: true
---

# /investigate $ARGUMENTS

Read-only: no edits to the repo. Probe scripts go in the scratchpad.

1. **Restate** the question as something falsifiable, in one line.
2. **Look.** Read the code path yourself when it is a few files; spawn up to 3 `Explore` agents in
   parallel when it spans subsystems. Check `.refactor/decisions.md`, `docs/design_note.md` and
   `docs/open_issues.md` for prior answers.
3. **Measure** when the answer is a number or a behaviour: a probe script that builds the case
   (`tests/cases.py` builders or `examples/henon_cases.py`, explicit `area_cutoff`) and prints
   what decides the question. Run it; do not guess what it would print.
4. **Answer** in a few paragraphs: the answer first, then the evidence (`file:line`, the probe's
   numbers), then what remains uncertain. Distinguish a numerical artifact (under-resolution, a
   same-stability "crossing") from a real behaviour, per CLAUDE.md's invariants.
5. **Offer** the next step in one line: a `.refactor/decisions.md` entry, a memory, a `/refactor`
   or `/feature`, or nothing.

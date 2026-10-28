---
name: quick
description: A small, well-understood edit to tanglepack (a fix, a rename, a docstring, a parameter) made directly on the current branch and checked with the relevant tests.
argument-hint: <the edit>
disable-model-invocation: true
---

# /quick $ARGUMENTS

1. Read the code, make the edit (follow the `house-style` skill), and keep it to what was asked.
2. Run the test tier that covers it (`MPLBACKEND=Agg env/bin/python -m pytest -q tests/<tier>`, or
   the single test file); the whole suite if the edit touches `numerics/` or `topology/` logic.
3. Show the diff summary and the test result. Commit only when asked.

If the edit turns out to touch more than ~3 files, change behaviour the request did not mention,
or need a design choice, stop and suggest `/refactor` or `/feature` instead.

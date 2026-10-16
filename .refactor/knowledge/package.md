# package — knowledge file (what crosses the layers: docs, CLAUDE.md, Sphinx, campaign references)

Created 2026-10-15 by the B12 analyst at `HEAD` `94fca6b` (read-only run). Updated 2026-10-16 by
the B12a auditor: branch `refactor/one-design-note` @ `b39c7e2` (NOT yet merged; main is still
`94fca6b`). Where the two differ, "B12a" below is the branch state.

## B12a state (branch `b39c7e2`)

- `docs/design_note.md` (722 lines) is THE design: `## A.0` global rules, `A.1`-`A.8`, `A.8b`
  minimal trellis, `A.8c` word, `A.8d` iterated partition / refined level, `A.9`-`A.11`, then
  `H` history (slices B0-B12a with merge hashes, WPs, superseded rules, the 9 audit divergences
  closed without a chalkboard row), `S` stress findings that bind, `O` open items, `I` id index.
  Each stage: Definition / anchors "(D13)" with the statement folded in / "Implemented by" names.
  Ids resolve in `.refactor/decisions.md` and `docs/open_issues.md`.
- `src/`: 0 decision ids except 7 inside log strings (`TangleWorkbench.py:1529` NUM-3,
  `BridgeClass.py:618` D2, `MinimalTrellis.py:355` D9, `StablePartition.py:338` D9, `:430` NUM-5,
  `topology/diagnostics.py:713,721` D6) and 1 date (`PartitionFamily.py:812`); 0 `B*/Q*/E*/F7/N-a`;
  15 "design note §A.x" pointers. Old-code narration outside the history-word list survives:
  `BridgeIterator.py:20-21,44` ("plan row 1.1/2.2 removed"), `BaseManifold.py:16-17`,
  `ManifoldMachine.py:632,700-701`, `Tangle.py:593,1017-1018`, `DualGraph.py:50-51`,
  `geometry.py:9-11`.
- `CLAUDE.md`: 185 lines / 62,737 bytes; `## The pipeline, stage by stage` (`:28-61`, one
  paragraph per stage + "Global rules", each ending "(design note §A.x)") replaces "Bridge classes
  and the dual graph"; `:1-27` and `## Coding Style`..EOF byte-identical to the pre-B12a copy;
  table has `cartoon.py` and `topology/diagnostics.py` rows (numerics `diagnostics` is in the
  "Intersection detection" prose because `rowcheck.py` keys rows by basename). namecheck: 23
  "misses", all builtins / protected lines / panel names, none a dead name.
- Both the note (A.11) and CLAUDE.md (`:60`) say the session has "no fan-outs": imprecise,
  `TangleSession.classify_strong_pips(fixed_point=None)` still loops over fixed points (B9a kept it
  for pip memory); the comment above it (`TangleSession.py:540-546`) is accurate.
- MERGE HAZARD (believed, not run): main's working tree holds the same paths as commit 0 as
  modified / untracked files, so `git merge --ff-only` in main refuses until they are cleared;
  commit 0 equals main's copies except the post-copy auditor edits to `.refactor/`.
- Proof tools (scratchpad `b12/`): `stripped_ast.py BASE BRANCH` (docstring-free AST equality per
  file; comments never reach the AST), `namecheck.py ROOT FILE`, `rowcheck.py`, `cites_src.txt`
  (the 313 base citations, 86 ids). Reuse them for any comments-only slice.

## Pre-B12a state (main `94fca6b`, still main's HEAD)

## Map

- `CLAUDE.md` (working copy 173 lines; `HEAD` 171 + the two sweep-script lines). Sections:
  preamble `:1-10`, physics `:11-21` "Fundamental Invariant" and `:22-27` "Area Preservation"
  (author's words; NOTE they carry campaign references: `:20` "Since WP10 ... since 2026-10-07",
  `:26` "(NUM-20, WP10b ...)"), "Bridge classes and the dual graph" `:28-50` (the campaign's
  running log, dated author rules inline), Commands `:51-82`, architecture + table `:83-162`,
  Coding Style `:164-173`.
- `docs/alignment_plan_2026-10-09.md` §A `:27-432` (the design-note draft), §B slices, §C `:1129`,
  §E Q1-Q16. `docs/pipeline_chalkboard.md` (stages 0-10; divergence table `:81-103`, 17 rows, dated
  2026-10-07). `docs/open_issues.md` holds every D/C/NUM/INT/PK/BCG/NEST/API/SIMP/PERF/INV statement
  (Author decisions table `:1239-1279`, WP table `:1281-1301`); C11-C16, C21 also
  `.refactor/decisions.md:78-129, :194`; DEC-25 `.refactor/decisions.md:214`.
- `docs/source/*.rst`: automodule list = 41 directives, 40 distinct targets (`tanglepack` twice) for 51 `.py` files; `tanglepack.numerics.graphviz`
  (`tanglepack.numerics.rst:183-189`) no longer exists; 12 modules have no automodule:
  `examples.henon_cases`, `loom.{BridgeAlphabet, cartoon, Topology}`,
  `topology.{BridgeClass, diagnostics, DualGraph, DualWalk, ElementNaming, MinimalTrellis,
  PartitionFamily, SymbolicDynamics}`. `conf.py` adds `../src` (resolves to `docs/src`, which does
  not exist; autodoc works only because the package is installed) and needs `myst_parser` + `furo`
  (not installed in `env/`).

## Invariants / facts

- Decision-id citations (regex `D|C|WP|NUM|INT|SIMP|BCG|PK|NEST|API|DEC|PERF|INV`): 313 in `src/`
  (30 files, 86 distinct ids), 381 in `tests/`, 176 in `CLAUDE.md`.
- CLAUDE.md name check (scripts `scratchpad/b12/namecheck.py`, `rowcheck.py`): 23 hard misses on
  12 lines (`_resolve_inertness`, `image_classes`, `image_loops`, `graph(bridges=)`, `graphviz`,
  `element_for`, `_anchor_branch`, `measured_side` / `side_source` on `Hole`, `EvidenceDiagnosis` /
  `ImageEvidence` / `registered_evidence` under SymbolicDynamics, `TangleSession.cartoon_enclosures`,
  five session plot delegates + `cartoon_zones()`, `terminated_early`,
  `tests/test_higher_period_cartoon.py`); no row for `loom/cartoon.py` or `topology/diagnostics.py`.

## Smells

- Campaign history in docstrings (361 src lines carry an id or a date).
- `CLAUDE.md` "Invariants the code asserts" lists two non-asserts: `image_of(region)` preserves
  area (an acceptance filter, `Arrangement.image_of`) and "refined classes never change
  `BridgeClass` identity" (structural, no assert).

## Features

- Sphinx docs: unbuildable as configured (missing extensions); 1,206 warnings when forced.

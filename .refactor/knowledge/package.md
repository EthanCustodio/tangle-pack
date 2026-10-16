# package — knowledge file (what crosses the layers: docs, CLAUDE.md, Sphinx, campaign references)

Created 2026-10-15 by the B12 analyst at `HEAD` `94fca6b` (read-only run).

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

# Request: WP10b — Remaining code issues before the final verification

*From the orchestrator, 2026-10-09. Follows WP9c; precedes WP12.*

## The question

**Close the code issues still open in `docs/open_issues.md` that are not WP12's rerun/doc work.** Re-check each against current main first.

1. **NUM-12 (WP10 slice 2):** growth after a blast loses every crossing born on blast image curves, because `compute_intersections` re-indexes manifolds only; it also recuts the bridge set and clears the `iterated` flags (blast genealogy lost). WP10's proposed fix: keep the image curves, re-index them query-only and let `create_bridges` cut them. Under D1 the cycle "grow / blast / recompute" must be possible in any order without losing computed crossings. Add a law/regression: blast, then grow one step, then recompute — every crossing registered before growth that still lies on a computed curve is still registered (by position), and the class count does not decrease.
2. **NUM-20:** on the p1(u) x p3(s) recipes (nested 13u/15s, 13u/16s) heteroclinic iterate links fail `area_along_iterate_links` and `iterate_link_scales_by_beta` (e.g. link 7 -> 215). Root-cause `IterateInference` for A(u) x B(s) crossings (the bracket rule and per-step scaling must use each crossing's OWN unstable and stable fixed points) and fix.
3. **PK-4 (WP7 slice 2):** foreign-branch chords of a lobe excursion ("abuts no hole" at crossing 44 on p3_6_blasts). WP7 proposed: pair chords by STABLE order (Jordan alternation), extending the author's empty-stretch rule; it changes only refinement labels, never walks or words. Implement if it is the clean rule (a Claude's call under D8, record it as a C-number).
4. **WP11 review nits:** delete dead `ManifoldMachine._parabolic_fit` and its pinning test; drop the always-zero `tol` from `CartoonLayout`/`ZoneLayout`/`_cluster`/`_rank_of`; fix the stale `span_contains`/`owns_cdist` docstrings; strengthen `test_stable_trim` with an exact-node case.
5. **NUM-19 (possibly stale):** verify each sub-item and close or fix.

## What "done" means

Each item closed in `docs/open_issues.md` with evidence, or explicitly left with a reason. No behaviour change on law cases except where an item says so (prove with id-free signatures). KNOWN_ISSUES stays at 1. Suite near ~132 s.

## Binding context

- `docs/open_issues.md` (D1-D14, C1-C18) and `docs/pipeline_chalkboard.md`; `.refactor/decisions.md` is binding.
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins.

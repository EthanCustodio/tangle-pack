# Request: WP9b — Break walk ties with the non-crossing constraint (D13)

*From the author, 2026-10-09. Follows WP10.*

## The question

**When two or more shortest dual-graph walks of equal length read different words, which one is the image?** The author's rule (D13, DEC-23): f(bridge) cannot cross f(W^u) — the images of the unstable arcs that bound a lobe are walls for the walk. Enrich the walk with that non-crossing constraint so the side is picked combinatorially. It stays a walk (D3); no registered image is needed and topology never calls the map.

Known ties:
- `heteroclinic` (WP9): with D7 the heteroclinic lobe's hole splits B's window element; two equal walks pass the lobe on either side (they differ by a loop around that hole); 3 classes ambiguous; strict xfail `(every_active_class_walks, heteroclinic)` tagged WP9.
- k28 at >= 8 blasts (WP8b/WP10): ambiguous closure walks of the same kind.

## What "done" means

- Derive the constraint precisely (which unstable arcs' images bound the walk, how they enter the dual graph: e.g. as walls between faces, or as a winding/side constraint on walks), state it on the chalkboard (stage 8), and implement it in `DualWalk` / `DualGraph`.
- The heteroclinic KNOWN_ISSUE flips; k28 at 8+ blasts has no ambiguous walk; no law case changes its words except where it was ambiguous (prove with id-free signatures).
- A law: wherever a member's registered image exists, the walk chosen agrees with it at word level (evidence check), and no class stays `ambiguous` on any law case.
- Docs: CLAUDE.md symbolic-dynamics bullet, chalkboard, knowledge, `docs/open_issues.md` (DEC-23 implemented, NEST-1 closed if done).

## Binding context

- `docs/open_issues.md` (D1-D13, C1-C12) and `docs/pipeline_chalkboard.md`; `.refactor/decisions.md` is binding.
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins; KNOWN_ISSUES only shrinks. Keep the suite near its current ~136 s.

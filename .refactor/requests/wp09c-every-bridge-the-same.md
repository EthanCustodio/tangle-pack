# Request: WP9c — Every bridge is the same (D14)

*From the author, 2026-10-09. Follows WP11.*

## The question

**Where does the heteroclinic case's member-less connecting class `bn` = {A:R_1, B:L_(2.0;1)} come from, and what is the uniform rule that removes the special case?** The author (D14, quoted in `.refactor/decisions.md`): "If there is no computed bridge that connects the tangles then there is no bridge class which connects the tangles either. Treat every bridge the same. They just connect two intersection points. Those two intersections lie in stable partition elements. Side matters for those partition elements. A bridge connecting two partition elements on the same tangle is no different than a bridge connecting two partition elements on different tangles."

Read D11 together with D14: a pair a walk names is a class (D11), but a pair across tangles that no bridge connects means the walk (or landing, or partition) that named it is wrong.

## What "done" means

- Trace `bn` to the walk that names it (`named_by`), find why that walk crosses from A's element to B's element where no bridge connects them, and fix the construction (landing, partition element across tangles, unified window, walk admissibility, ...). Then the strict xfail `(every_active_class_walks, heteroclinic)` flips with no re-scoping of the law.
- Remove every special case for connecting / cross-tangle classes that the uniform rule makes unnecessary (in classes, walks, marker law, plotting ordering only if it is not the author's ordering rule).
- A law: every class in the closed alphabet whose two elements lie on different fixed points has a computed bridge member (or the case is a proven construction error).
- No behaviour change on any non-heteroclinic law case or deep run (id-free signatures).
- Docs: CLAUDE.md, chalkboard (the ★ several-saddles note), knowledge, `docs/open_issues.md` (DEC-24 implemented, NEST-1).

## Binding context

- `docs/open_issues.md` (D1-D14, C1-C16) and `docs/pipeline_chalkboard.md`; `.refactor/decisions.md` is binding.
- `CLAUDE.md` invariants are physics. Topology never calls the map. `numerics` and `topology` never import `loom`.
- No golden pins; KNOWN_ISSUES only shrinks. Suite near ~132 s.

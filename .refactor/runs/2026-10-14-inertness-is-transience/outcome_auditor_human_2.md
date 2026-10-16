# Outcome — B5 `inertness-is-transience`

Branch `refactor/inertness-is-transience` @ `3fb310a` (worktree `../tangle-pack-inertness-is-transience`), base `7d3a294`. The plan was approved at round 2, and review passed in round 1 of coding. There is no debrief verdict yet. Nothing is merged; the merge is yours. The main tree still carries uncommitted docs and `.refactor/` changes; the branch does not touch them.

## What "inert" means now

Inertness is now decided in two places.

**1. Before the dual graph, your zeroth-hole rule.** A bridge that carries a zeroth hole is inert. A zeroth hole is the one punched on a pseudoneighbor pair (`hole.pair is not None`). The minimal trellis now maps every hole bridge except these, and it decides this from `holes` alone, without reading the class table.

The theorem behind the rule:
- the lobe `b ∪ σ` of a zeroth hole holds exactly one marker, the hole itself;
- holes propagate backward only, so that hole's image is no marker;
- so `f(b) ∪ f(σ)` encloses no marker, `f(b)` is homotopic to the stable segment `f(σ)`, and the class word is empty.

It is proven under two conditions: no selected orbit lies inside the lobe, and no other hole sits beside `b`. It holds on all 15 builds.

**2. After the words, transience.** Every class gets its verdict here (Q3, I-a). A class is inert iff its word is empty or names only inert classes. An unresolved class is active, and so is every class whose word reaches one. Transience has to stay the definition because the zeroth-hole rule is only sufficient: one hole-free class, heteroclinic `u`, is inert.

## The answer to your question now

**What makes a class inert?** No cycle (and no unresolved class) can be reached from it in the class word graph. Registered image evidence plays no part any more.

Verified on the 15 builds (the 11 law cases, `p3_6_blasts`, `p3_16`, k28 at 8 blasts, nested with 6 inner blasts):
- **Unchanged:** words, minimal trellis, `P1`, refined classes and `is_reliable`.
- **Changed:** five classes go from active to inert. Base `7d3a294` had called each of them active only because it had no registered image of them.

## What changed

The two hunks that matter most.

**`minimal_trellis` (`MinimalTrellis.py:300-335`).** It decides from the holes, not the class table.
```python
# before                                  # after
entry = classes.entry_of(bid)             classes.entry_of(bid)        # only for the unclassed WARNING
...                                       ...
if not entry.active:                      if frozenset(bid) in zeroth: # hole.pair is not None
    inert.append(bid); continue               inert.append(bid); continue
```

**The inertness fixed point.** The old `BridgeClass._resolve_inertness` (24 lines over `image_classes` and `has_image_evidence`, fed by the 40-line `_collect_image_evidence`) is replaced by `SymbolicDynamics.inert_classes` (`:499`):
```python
inert: set[BridgeClass] = set()
while grown := {cls for cls, word in words.items()
                if cls not in inert and word is not None and all(token in inert for token in word)}:
    inert |= grown
return inert
```
It is called once (`:1471`), after the closure loop has recorded each class's class-level tokens (`:1465`) and before translation. Then:
- `entry.inert` is set on every entry;
- letters follow in one pass, active classes through `letter_for` and inert ones through the `u, v, …` series.

`Topology.classes` (`loom/Topology.py:152`) runs `self.dynamics` and returns the zoned `_alphabet`. Its lettering loop is gone. PK-8's WARNING is deleted.

## What got shorter or clearer

- **Source:** 170 lines added, 257 deleted; the four files go from 3,055 to 2,968 lines. `BridgeClass.py` drops from 920 to 791. `SymbolicDynamics.py` grows from 1,490 to 1,516, almost all of it docstrings plus the 6-line fixed point. `MinimalTrellis.py` grows from 438 to 449, all of it the new Dev Note giving your argument.
- **Deleted:** `_resolve_inertness`, `_collect_image_evidence`, `has_image_evidence`, the `image_classes` / `image_loops` / `unresolved` slots, and the image section of `describe()`. The grep for them over `src` and `tests` is empty.
- **Order:** the verdict now has to exist before any letter is assigned. Before, lettering in `Topology.classes` used a verdict computed from the evidence.

## Every behaviour change

1. **Words, minimal trellis, `P1`, refined classes, `is_reliable`:** identical on all 15 builds. The test engineer's comparison script also checked the itineraries, the iterated element names and the refined matrix on the remaining nodes.
2. **Five verdicts go from active to inert**, each with an empty word and a `trivial` walk:
   - inversion_two_blasts `g`;
   - heteroclinic `d` and `bj`;
   - p3_6_blasts `ab`;
   - nested-inner-6 `ae`.

   Each takes an inert letter and leaves the transition graph and the matrix. In a fresh session, the active letters after it shift down by one. No class went from inert to active.
3. **`MinimalTrellis.unmapped` shrinks:** by 1 on inversion_two_blasts, p3_6_blasts and nested-inner-6, and by 2 on heteroclinic. Those bridges now appear in `inert_hole_bridge_ids`, which now means "carries a zeroth hole".
4. **`session.bridge_classes()` / `Topology.classes` now runs the whole dynamics.** This costs more for scripts that only want the table.
5. **One WARNING is reworded.** The "class %s unresolved" WARNING now names the class by its label (`A:p1@0.0/L#0 <-> A:p1@0.0/L#2`) instead of a letter, because letters do not exist yet during the walk. It fires on no build, only on the `p3_partitioned` fixture. Otherwise no new WARNING appears, and PK-8's is gone (it never fired).
6. **API removals:** `BridgeClassEntry.image_classes`, `.image_loops`, `.unresolved` and `.has_image_evidence` are gone, and `describe()` takes no `names` argument. `entry.inert` is now `Optional[bool]`: `None` until the dynamics runs, and on evidence entries outside the alphabet. Before the dynamics, `BridgeClassTable.active` returns every entry. `inert_classes` is new and public.

## What was verified, and how

- **Suite on the branch (I ran it):** 1193 passed, 38 skipped, 5 deselected, 1 xfailed in 200 s. `-rxX` lists only `orientation_reversing` (C3, KNOWN_ISSUES). The base count was also 1193: −3 deleted tests, +2 `inert_classes` cases, +1 test-engineer test.
- **Diff (I read all of it):** `git diff --stat 7d3a294..3fb310a` touches exactly the nine planned files.
- **inversion_two_blasts, base vs branch (my own script, classes keyed by element identity):** every class's word is identical. The only differences are `{0.0/left[23,2] 0.0/left[3,24]}` going from `g` active to `v` inert, with an empty word, and `unmapped [(23,24)] -> []`. `is_reliable` is True on both.
- **`p3_partitioned` fixture, base vs branch (my own script):**
  - bridge `(5,6)` goes from `c` active to `u` inert;
  - bridge `(4,5)` stays `b` active and unresolved;
  - `unmapped [(4,5),(5,6)] -> []`;
  - `is_reliable` is False on both;
  - the WARNING carries the label on the branch.
- **Not re-run by me (from the test engineer's and reviewer's logs):**
  - the full 15-build comparison;
  - the base-vs-branch runs of the new laws (on the base with the branch's `laws.py`, `inert_is_transient` fails on exactly `g`, `d` and `bj`, and the empty-word half passes everywhere);
  - the break-it mutations.
- **Callers outside the diff:** `scripts/henon_bridge_classes.py` and `stress/numerics/e11_blast_depth.py` read `entry.inert` after `topology.classes`. That property now runs the dynamics, so they get the new verdicts. I did not run them.

## What did not go to plan

- **The plan took two rounds.** Round 1 mapped every hole bridge. That would have made `P1` finer on 8 builds and added three WARNINGs for skipped image pairs (D9). You sent it back with the zeroth-hole rule, and round 2 kept the minimal trellis exactly as at base.
- **The renamed minimal-trellis test had to be rewritten.** The planner had expected `assert minimal.unmapped` to keep passing. It could not: on `p3_partitioned` every bridge that was unmapped at base carries a zeroth hole, so `unmapped` is now empty. The test now strips `pair` from the holes and re-runs `minimal_trellis`, so the unmapped path still runs on real data.
- **More unit tests needed adapting than the plan listed:**
  - two PK-8 tests were deleted; one of them assumed an inert class whose word is itself, which cannot exist under transience;
  - two closure tests were rewritten so the verdict comes from the word;
  - four `test_bridge_class.py` tests lost reads of the deleted evidence.
- **One rubric check fails a different way than the plan described.** Rubric 10's mutation (`word is not None` → `True`) makes the test fail with a `TypeError` before its assertion is reached. The test engineer's semantic mutation (`word or ()`) does fail on the assertion.
- **The test engineer added a test, and it closes a real gap.** `test_the_session_alphabet_carries_the_dynamics_verdicts` snapshots `session.bridge_classes()` before the dynamics is read. It catches break 5 (dropping `self.dynamics` from `Topology.classes`), which the law tier misses: an earlier law runs the dynamics on the shared case first.

## Where I am least sure this matches what you meant

1. **Your rule is a sufficient condition, not the definition.** You wrote "every bridge that has an initial hole punched in it should be inert". The branch uses that to decide what the minimal trellis maps. The verdict, though, is transience, because heteroclinic `u` holds no hole and is still inert. If you meant inert to be the zeroth-hole property itself, `u` would turn active.
2. **One zeroth-hole class stays active.** On the unblasted nested fixture `p3_partitioned`, the zeroth-hole bridge `(4,5)` is not mapped, yet its class is active, because it is unresolved (a singleton landing with no registered chain). This is outside the theorem's hypothesis, which needs a resolvable walk. It is reported, not pinned, and no build of the 15 has such a class. Should an unresolved zeroth-hole class be called inert by your rule, or stay active by convention?
3. **Two loose ends:**
   - `inert_classes` is public but not re-exported from `tanglepack.topology`, because `__init__.py` was outside the file list.
   - `CLAUDE.md` is now stale on inertness. It still says "Inert = at least one resolvable image…", "Inert hole bridges are kept but not mapped" and "`session.bridge_classes()` IS that alphabet" with evidence-based verdicts. This is deferred to B12.

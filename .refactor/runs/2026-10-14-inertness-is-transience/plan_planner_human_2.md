# Plan — B5 `inertness-is-transience`, round 2 (base `7d3a294`)

**Your rule holds, and it now decides what the minimal trellis maps.** A bridge holding a zeroth
hole (the hole punched on a pseudoneighbor pair, iterate 0) is inert. The reason is short:
- the pair's bridge `b` and the stable segment `σ` between its ends hold no point of X, so no arc of
  the trellis enters the lobe `b ∪ σ`, and the lobe holds exactly one marker, its own hole;
- holes propagate backward only, so that hole's image is no marker, and `f(b) ∪ f(σ)` encloses none;
- so `f(b)` is homotopic to the stable segment `f(σ)`, its ends lie in one element, and the word is
  empty.

The analyst checked this on all 15 builds: every zeroth-hole class has an empty word (a `trivial`
walk). It is a theorem under two conditions that hold on every build: no selected orbit sits inside
the lobe, and no other hole sits beside `b`. It needs only `trellis.holes`, so the minimal trellis
can decide it before the dual graph exists. The rule is sufficient but it is not the definition.
Heteroclinic `u` holds no hole and is still inert (its one member is a piece of the anchor bridge's
image, and that member's image is a loop). Transience over the words therefore remains the verdict
for every class (Q3, I-a), and a new law checks your rule against the words. Your note is recorded
in `.refactor/decisions.md` (2026-10-14).

## 1. The answer

**What changes:** inertness moves from the registered-evidence fixed point to two places:
- your zeroth-hole rule, which decides which hole bridges the minimal trellis maps;
- transience in the word graph, which gives every class its verdict after the words.

**What stays:** with this split, the minimal trellis, `P1`, the image bridges, the refined classes,
the words and `is_reliable` are identical to `HEAD` on all 15 builds (analyst prototype R2). The
15 builds are the 11 law cases, `p3_6_blasts`, `p3_16`, k28 at 8 blasts, and nested with 6 inner
blasts.

**The verdict changes:** five classes go from active to inert:
- inversion_two_blasts `g`;
- heteroclinic `d` and `bj`;
- p3_6_blasts `ab`;
- nested-inner-6 `ae`.

All five are zeroth-hole classes that `HEAD` called active only because it had no registered image
of them. No class moves the other way.

## 2. Recommendation

Make the change. Round 1's plan stands except for step 1, and it gains one law.

## 3. Scope

Slug **`inertness-is-transience`**. The branch is `refactor/inertness-is-transience` in
`../tangle-pack-inertness-is-transience`, made from `7d3a294`.

Files (nine, as in round 1; see question 1):

1. `src/tanglepack/topology/BridgeClass.py`
2. `src/tanglepack/topology/MinimalTrellis.py`
3. `src/tanglepack/topology/SymbolicDynamics.py`
4. `src/tanglepack/loom/Topology.py`
5. `tests/helpers/laws.py`
6. `tests/unit/topology/test_bridge_class.py`
7. `tests/unit/topology/test_symbolic_dynamics.py`
8. `tests/unit/topology/test_minimal_trellis.py` (the assertion at line 94)
9. `tests/helpers/direct.py` (its pre-lettering loop, lines 107-125)

Size: roughly 250 changed lines, mostly deletions in `BridgeClass.py`.

## 4. Steps

1. **`MinimalTrellis.py`: map every hole bridge except one carrying a zeroth hole.**
   - In the hole-bridge loop, note which bridges carry a zeroth hole. A zeroth hole is
     `hole.pair is not None`.

     ```python
     hole_ids.setdefault(frozenset(bridge.id), bridge.id)
     if hole.pair is not None:   # a zeroth hole: its lobe's image holds no marker, the word is empty
         zeroth.add(frozenset(bridge.id))
     ```

   - In the mapping loop, the unclassed check stays first. Then
     `if frozenset(bid) in zeroth: inert.append(bid); continue` replaces `if not entry.active`
     (`:326`). The class table is no longer read for this decision.
   - Keep the field `inert_hole_bridge_ids`, but re-document it: "the hole bridges carrying a
     zeroth hole: their image is homotopic to a stable segment (the class word is empty), so they are
     kept and never mapped". Its summary and describe lines stay.
   - Update three docstrings: the class docstring ("images of the active ones" becomes "images of
     every hole bridge without a zeroth hole"), the `classes` argument (it no longer decides mapping),
     and the `unmapped` attribute ("the mapped-eligible hole bridges with no registered image").
   - Fix the DEBUG text at `:331`: "hole bridge %s has no registered image".

2. **`BridgeClass.py`: delete the evidence machinery.** Unchanged from round 1.
   - Delete `_resolve_inertness` (`:826-849`), `_collect_image_evidence` (`:783-823`) and their call
     sites (`:712-714`).
   - Delete `has_image_evidence` (`:250-253`) and the `image_classes`, `image_loops` and `unresolved`
     slots (`:230-233`), with their docstring lines and their part of `describe()`.
   - `_image_chain` stays: loop folding and the minimal trellis read it.
   - `inert` stays a slot, `None` until the dynamics sets it.
   - The Dev Notes "INERT = MAPS TO NO ACTIVE CLASS" (`:51-71`) become one paragraph. Inert is
     transience in the word graph (Q3, I-a). Your zeroth-hole rule is the sufficient condition the
     minimal trellis uses (2026-10-14). The 2026-09-16 evidence form gets one line of history.

3. **`SymbolicDynamics.py`: the verdict, then the letters.** Unchanged from round 1.
   - A public `inert_classes(words)` returns the fixed point "a class is inert iff its word is empty
     or names only inert classes". An unresolved class (word `None`) is active, and so is every class
     whose word reaches it.

     ```python
     def inert_classes(words: dict[BridgeClass, Optional[list[BridgeClass]]]) -> set[BridgeClass]:
         """Q3 (I-a): inert iff no cycle is reachable; an unresolved class is active."""
         inert: set[BridgeClass] = set()
         while grown := {cls for cls, word in words.items()
                         if cls not in inert and word is not None
                         and all(token in inert for token in word)}:
             inert |= grown
         return inert
     ```

   - The closure walks every class and records its class-level tokens (`None` when it has no
     itinerary). Its `add` assigns no letter and copies entries with `letter=None, inert=None`.
   - After the closure:
     - set `entry.inert` from `inert_classes(words)`;
     - letter in one pass in table order: active classes through `letter_for` (else `letter(i)`),
       inert ones through `inert_letters`;
     - `entry.letter` is set on active entries only;
     - translate as today.
   - The initial lettering (`:1383-1386`) goes. So does PK-8's WARNING (`:919-925`): under I-a an
     inert word names only inert classes by definition.

4. **`Topology.py`: verdicts before anyone reads the alphabet.** Unchanged from round 1.
   - A `_alphabet` cached property zones `self.minimal.alphabet`, and `dynamics` passes it.
   - `classes` touches `self.dynamics` and then returns `_alphabet`.
   - The lettering loop (`:148-150`) goes. `session.bridge_classes()` now runs the dynamics.

5. **`tests/helpers/laws.py`.**
   - **New law `inert_is_transient`** (layer `symbolic`), as in round 1:
     - build the graph over ALL classes of the closed table from `cd.symbols` (loops excluded);
     - assert `kind == "inert"` iff neither a cycle nor an unresolved class is reachable from the
       class;
     - compute reachability by its own DFS, never by calling `inert_classes`.
   - **New law `zeroth_hole_classes_have_empty_words`** (layer `symbolic`). For every hole with
     `pair is not None` in `_holes(case)`, take the class of its bridge
     (`minimal.classes.entry_of(bridge.id)`). Its `ClassDynamics` must have an itinerary, no symbol
     after loops are dropped, and `kind == "inert"`. Register it beside `inert_classes_outside_transitions`.
   - **`minimal_trellis_bridges`:**
     - `{frozenset(b) for b in minimal.inert_hole_bridge_ids}` equals the set of bridges carrying a
       zeroth hole, computed from `_holes(case)`;
     - no bridge in `image_chains` or `unmapped` carries a zeroth hole;
     - every other hole bridge is unclassed, mapped along its registered chain, or in `unmapped`;
     - the three `entry_of(bid).active` reads (`:2312`, `:2314`) go.
   - `anchor_bridge_class_leads_its_tangle` reads the verdict from the dynamics
     (`dyn.classes[cls].kind == "active"`), as in round 1. `inert_classes_outside_transitions` and
     `only_active_classes_lettered` are unchanged.

6. **Unit tests.**
   - `test_bridge_class.py`: delete `test_resolve_inertness_fixed_point` and its `_entry` helper.
   - `test_symbolic_dynamics.py`: two cases on `inert_classes`:
     - `a -> b, b -> (empty), c -> c`: `a` and `b` inert, `c` active;
     - `a -> b, b unresolved`: both active.
   - `test_minimal_trellis.py:94`: `assert table.entry_of(bid).active` becomes "unmapped ⇒ no
     zeroth hole": no hole of `trellis.holes` with `pair is not None` has `bid` as its bridge. The
     test is renamed `test_a_hole_bridge_without_a_registered_image_is_unmapped`.
   - `direct.py`: delete the pre-lettering loop. Its `classes()` returns the unlettered alphabet.

## 5. Behaviour

- **Unchanged on every build** (all 15): class-level words, the minimal trellis (kept, image and
  dropped bridges), `P1` (iterated elements and names), refined classes and their matrix, the closed
  alphabet's class set, and `is_reliable`.
- **Inert verdicts** (reported, not pinned): the five classes in section 1 go from active to inert.
  - Each takes an inert letter. In a fresh session, every active letter after it shifts down by one.
  - Each leaves the transition graph and matrix. The matrix on the remaining nodes is unchanged.
- **`MinimalTrellis.unmapped` shrinks:** −1 on inversion_two_blasts, p3_6_blasts and
  nested-inner-6, −2 on heteroclinic. These are the five classes' zeroth-hole bridges, which now
  move to `inert_hole_bridge_ids`.
- **Logs:** no new WARNING. PK-8's WARNING is gone; it fired on no build.
- **API:**
  - `BridgeClassEntry.image_classes`, `.image_loops`, `.unresolved` and `.has_image_evidence` no
    longer exist;
  - `entry.inert` is `None` on evidence entries outside the alphabet;
  - `MinimalTrellis.inert_hole_bridge_ids` keeps its name and now means "carries a zeroth hole";
  - `session.bridge_classes()` computes the dynamics.

## 6. Not doing

- **A pre-dual-graph rule for hole-free classes.** "A lobe holding no marker whose image is a
  marker" would catch heteroclinic `u`, but it needs `DualGraph.enclosed`; transience covers it.
- **Proving propagated-hole classes are always active.** They are active on the whole corpus, but
  this is unproven.
- **Seeding `entry.inert` in `bridge_classes` for zeroth-hole classes.** Nothing but the minimal
  trellis would read the seed.
- **Other slices:** refinement rules (B7), the inert letter series itself, any quotient.
- **Docs:** `CLAUDE.md`'s inertness text ("Inert hole bridges are kept but not mapped" becomes
  "zeroth-hole bridges"; "Inert = at least one resolvable image…" becomes transience) and the
  alignment plan's B5 wording are left to B12 / docs. The tree's `CLAUDE.md` is dirty and not the
  branch's.

## 7. Rubric

1. `MPLBACKEND=Agg env/bin/python -m pytest -q` passes on the branch. `-rxX` lists only
   `cases.KNOWN_ISSUES`.
2. **Pinning tests on both commits.**
   - The base commit is a detached worktree at `7d3a294` with the branch's `laws.py` copied in. On
     it, the empty-word half of `zeroth_hole_classes_have_empty_words` (itinerary present, no symbol)
     passes on every law case.
   - On the base, `inert_is_transient` and the `kind == "inert"` half of
     `zeroth_hole_classes_have_empty_words` FAIL on exactly inversion_two_blasts (`g`) and
     heteroclinic (`d`, `bj`), and pass on the other nine law cases. This shows both laws have
     teeth.
   - On the branch, both laws pass on every law case, and both new `inert_classes` unit cases pass.
   - The report quotes the failing classes, matched by element identity.
3. `git diff --stat 7d3a294..` touches exactly the nine files in section 3 (eight if question 1
   takes the alternative).
4. **No behaviour change beyond section 5.** One script runs base and branch on all 15 builds,
   matching classes by element identity, and its output goes in the report. Equal on base and
   branch:
   - letter-free class-level `rules`;
   - `minimal.kept_bridge_ids`, `image_bridge_ids` and `skipped_pairs`;
   - the iterated partition's element intervals;
   - the refined classes' iterated pairs;
   - `is_reliable`.

   The only differences allowed are the five verdicts and the `unmapped` deltas of section 5. No
   WARNING record appears on the branch that is absent on the base.
5. `grep -rn "has_image_evidence\|_resolve_inertness\|_collect_image_evidence\|image_classes\|image_loops" src tests`
   is empty.
6. `minimal_trellis` reads neither `entry.active` nor `entry.inert`
   (`grep -n "\.active\|\.inert" src/tanglepack/topology/MinimalTrellis.py` has no hit inside the
   function). On every law case, `inert_hole_bridge_ids` is exactly the set of bridges carrying a
   zeroth hole: the updated `minimal_trellis_bridges` passes.
7. In `test_minimal_trellis.py`, the old assertion `table.entry_of(bid).active` is gone. The renamed
   test asserts that no unmapped bridge carries a zeroth hole, and it passes.
8. `Topology.classes` contains no `letter_for` call. On a fresh session, `session.bridge_classes()`
   returns entries whose `inert` is never `None`. `only_active_classes_lettered` and
   `test_session_bridge_classes.py` pass unchanged.
9. `SymbolicDynamics.py` has no PK-8 WARNING. `inert_classes` is called exactly once, after the
   closure loop and before `translate_itinerary`. The closure's `add` assigns no letter.
10. The `a -> b, b unresolved` unit case marks both classes active, and it fails if `word is not None`
    is replaced by `True` in `inert_classes`.
11. The report lists the changed verdicts per build (expected: the five of section 1) and the
    `unmapped` deltas. Any difference from section 5 is explained, not pinned.

## 8. Questions

1. **Nine files or eight?**
   - `test_minimal_trellis.py:94` is edited because it states the old rule ("unmapped ⇒ active").
     Under the new rule the assertion would most likely still pass, but it would be true for a reason
     the code no longer uses.
   - `direct.py` pre-letters classes the way the deleted `Topology.classes` loop did, so it would
     letter classes the dynamics then calls inert.

   **Default (yes = this):** nine files, both edited as in step 6.
   **Alternative:** leave `direct.py` alone. That makes eight files, and `direct.py` still
   pre-letters classes the dynamics calls inert. Its views are letter-free, so its tests likely
   still pass.

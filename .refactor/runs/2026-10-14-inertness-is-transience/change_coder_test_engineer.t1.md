# Change — B5 `inertness-is-transience` (coder, iteration 1)

**Base commit:** `7d3a294` (author's `HEAD`, main).
**Branch:** `refactor/inertness-is-transience`, worktree `../tangle-pack-inertness-is-transience`, head `b9145d0`.
**Base worktree:** `../tangle-pack-inertness-is-transience-base`, detached at `7d3a294`. It is left clean: rubric 2's copy of `laws.py` was restored with `git checkout`.

## Commits (`git log --oneline 7d3a294..`)

```
b9145d0 Laws: inert_is_transient, zeroth_hole_classes_have_empty_words; minimal_trellis_bridges never maps a zeroth-hole bridge (B5)
ec33d31 Bridge classes: delete the registered-evidence inertness (B5)
8d22685 Symbolic dynamics: inert = transient in the word graph, lettered after the verdict (B5, Q3 I-a)
fa40cf5 Minimal trellis: map every hole bridge except one carrying a zeroth hole (B5, author 2026-10-14)
```

## `git diff --stat 7d3a294`

```
 src/tanglepack/loom/Topology.py               |  19 ++-
 src/tanglepack/topology/BridgeClass.py        | 187 ++++----------------------
 src/tanglepack/topology/MinimalTrellis.py     |  77 ++++++-----
 src/tanglepack/topology/SymbolicDynamics.py   | 144 ++++++++++++--------
 tests/helpers/direct.py                       |  20 +--
 tests/helpers/laws.py                         |  71 ++++++++--
 tests/unit/topology/test_bridge_class.py      |  79 +----------
 tests/unit/topology/test_minimal_trellis.py   |  39 ++++--
 tests/unit/topology/test_symbolic_dynamics.py |  66 ++++-----
 9 files changed, 298 insertions(+), 404 deletions(-)
```

These are exactly the plan's nine files (question 1's default).

## Suite

`MPLBACKEND=Agg env/bin/python -m pytest -q -rxX` on the branch:

```
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11: det J < 0 is out of scope (C3); FixedPoint.set_k_value rejects a saddle with exactly one negative eigenvalue
1192 passed, 38 skipped, 5 deselected, 1 xfailed in 203.75s (0:03:23)
```

The analyst's baseline at `HEAD` was 1193 passed, 38 skipped, 5 deselected, 1 xfailed. The net −1 is accounted for test by test:
- removed `test_resolve_inertness_fixed_point` (−1);
- removed the two PK-8 tests (−2);
- added the two `inert_classes` cases (+2).

The new laws ride the existing per-(layer, case) law tests, so they add no test items. The only xfail is the `cases.KNOWN_ISSUES` entry.

## What was done, step by step

1. **`MinimalTrellis.py`**
   - The hole loop records `zeroth` (`hole.pair is not None`).
   - The mapping loop keeps the unclassed check first, then `if frozenset(bid) in zeroth: inert.append(bid); continue`. The class table is consulted only for the unclassed WARNING; its entry is no longer read.
   - `inert_hole_bridge_ids` is re-documented as "carries a zeroth hole". Updated: the module text and a new Dev Note (the author's argument, 2026-10-14), the class, `classes` and `unmapped` docstrings, and the DEBUG text ("hole bridge %s has no registered image; not mapped").
2. **`BridgeClass.py`**
   - Deleted: `_resolve_inertness`, `_collect_image_evidence` and their call site; `has_image_evidence`; the `image_classes`, `image_loops` and `unresolved` slots, with their docstring lines and the images part of `describe()` (`describe()` lost its `names` argument, now unused).
   - `inert: Optional[bool] = None`.
   - The Dev Note is now one paragraph: transience (Q3, I-a); the zeroth-hole rule as the minimal trellis's sufficient condition (2026-10-14); one line of history for the 2026-09-16 evidence form.
   - `_image_chain` stays.
3. **`SymbolicDynamics.py`**
   - New public `inert_classes(words)`: the plan's fixed point, verbatim in substance.
   - The closure's `add` no longer letters: it copies evidence entries with `letter=None, inert=None` and tracks membership in a `listed` set.
   - The loop records `words[cls]`: the class-level tokens, loops dropped, `None` when unresolved. While walking, a `ClassDynamics` is created with `letter = cls.label` and `kind = "active"`.
   - After the loop, in this order:
     - `inert_classes(words)`, called once;
     - `entry.inert` set on every entry of the closed table;
     - `inert_letters(table)` plus `_active_letters(table, letter_for)`, one pass in table order;
     - `entry.letter` set on active entries only;
     - `cd.letter` and `cd.kind` set;
     - the closure INFO log;
     - translation as before.
   - The initial lettering is gone, and so are PK-8's WARNING in `_edges` and the now-unused `itertools.count`. The Dev Notes were updated.
4. **`Topology.py`**
   - A cached `_alphabet` zones `minimal.alphabet`.
   - `dynamics` passes `_alphabet`.
   - `classes` (now a plain property) touches `self.dynamics` and returns `_alphabet`.
   - The lettering loop is gone.
5. **`laws.py`**
   - New law `inert_is_transient`: its own DFS over all classes' `cd.symbols`; a class is persistent iff some node reachable from it, itself included, is unresolved or reaches itself. It also checks `entry.inert is (kind == "inert")`.
   - New law `zeroth_hole_classes_have_empty_words`: itinerary present, no symbol, `kind == "inert"`, asserted in that order.
   - Both are registered after `inert_classes_outside_transitions`.
   - `minimal_trellis_bridges`:
     - `inert_hole_bridge_ids` equals the zeroth set;
     - nothing in `image_chains` or `unmapped` carries a zeroth hole;
     - `hole_set − zeroth == unclassed ∪ mapped ∪ unmapped`;
     - the `.active` reads are gone.
   - `anchor_bridge_class_leads_its_tangle` reads `dyn.classes[cls].kind == "active"`.
6. **Unit tests and `direct.py`**: see "Test edits beyond the plan's letter" below.

## Rubric 4 — base vs branch on all 15 builds

Script: `scratchpad/b5cmp/dump.py`, run in each worktree with `PYTHONPATH=src:tests`; `compare.py` diffs the dumps. The scratchpad is `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/`.

Classes are keyed by element identity: `{fp:o.b/side[lo_id,hi_id], …}`, homotopy bounds.

Compared, letter-free:
- the class-level `rules`, as `(identity, direction)` tokens;
- `kept_bridge_ids`, `image_bridge_ids` and `skipped_pairs`;
- the iterated partition's `signature()` (every interval's ids and closedness);
- each refined child's iterated pair, by iterated interval identity;
- `is_reliable`;
- the class set;
- WARNING records (`tanglepack` logger) as a multiset.

Verbatim output (`compare_output.txt`):

```
== k10: 2 classes, inert 1 -> 1; image bridges 3 -> 3; iterated elements 11 -> 11; unmapped 0 -> 0 [] -> []; inert hole bridges 1 -> 1; skipped 0 -> 0; reliable True -> True; warnings 0 -> 0
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
== k28_one_blast: 3 classes, inert 2 -> 2; image bridges 0 -> 0; iterated elements 8 -> 8; unmapped 0 -> 0 [] -> []; inert hole bridges 2 -> 2; skipped 0 -> 0; reliable True -> True; warnings 0 -> 0
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
== k28_two_blasts: 4 classes, inert 1 -> 1; image bridges 3 -> 3; iterated elements 13 -> 13; unmapped 0 -> 0 [] -> []; inert hole bridges 1 -> 1; skipped 0 -> 0; reliable True -> True; warnings 0 -> 0
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
== p3: 5 classes, inert 2 -> 2; image bridges 0 -> 0; iterated elements 16 -> 16; unmapped 0 -> 0 [] -> []; inert hole bridges 2 -> 2; skipped 0 -> 0; reliable True -> True; warnings 0 -> 0
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
== nested: 9 classes, inert 3 -> 3; image bridges 3 -> 3; iterated elements 29 -> 29; unmapped 0 -> 0 [] -> []; inert hole bridges 3 -> 3; skipped 0 -> 0; reliable True -> True; warnings 0 -> 0
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
== nested_inner_blasts: 9 classes, inert 3 -> 3; image bridges 3 -> 3; iterated elements 29 -> 29; unmapped 0 -> 0 [] -> []; inert hole bridges 3 -> 3; skipped 0 -> 0; reliable True -> True; warnings 0 -> 0
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
== inversion: 4 classes, inert 2 -> 2; image bridges 0 -> 0; iterated elements 14 -> 14; unmapped 0 -> 0 [] -> []; inert hole bridges 2 -> 2; skipped 0 -> 0; reliable True -> True; warnings 0 -> 0
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
== inversion_one_blast: 7 classes, inert 1 -> 1; image bridges 3 -> 3; iterated elements 21 -> 21; unmapped 0 -> 0 [] -> []; inert hole bridges 1 -> 1; skipped 0 -> 0; reliable True -> True; warnings 0 -> 0
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
== inversion_two_blasts: 8 classes, inert 1 -> 2; image bridges 3 -> 3; iterated elements 22 -> 22; unmapped 1 -> 0 [[23, 24]] -> []; inert hole bridges 1 -> 2; skipped 0 -> 0; reliable True -> True; warnings 0 -> 0
    verdict {A:0.0/left[23,2], A:0.0/left[3,24]}: g active -> v inert (zeroth-hole class: True, word trivial)
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
== inversion_period2: 6 classes, inert 2 -> 2; image bridges 0 -> 0; iterated elements 24 -> 24; unmapped 0 -> 0 [] -> []; inert hole bridges 2 -> 2; skipped 0 -> 0; reliable True -> True; warnings 0 -> 0
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
== heteroclinic: 75 classes, inert 3 -> 5; image bridges 8 -> 8; iterated elements 158 -> 158; unmapped 2 -> 0 [[92, 93], [113, 114]] -> []; inert hole bridges 2 -> 4; skipped 0 -> 0; reliable True -> True; warnings 0 -> 0
    verdict {A:0.0/left[114,4], A:0.0/left[92,113]}: d active -> v inert (zeroth-hole class: True, word trivial)
    verdict {A:0.0/left[92,113], A:0.0/left[None,93]}: bj active -> y inert (zeroth-hole class: True, word trivial)
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
== p3_6_blasts: 43 classes, inert 1 -> 2; image bridges 4 -> 4; iterated elements 92 -> 92; unmapped 1 -> 0 [[42, 43]] -> []; inert hole bridges 1 -> 2; skipped 0 -> 0; reliable True -> True; warnings 0 -> 0
    verdict {A:1.0/right[4,43], A:1.0/right[42,3]}: ab active -> v inert (zeroth-hole class: True, word trivial)
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
== p3_16: 43 classes, inert 2 -> 2; image bridges 4 -> 4; iterated elements 92 -> 92; unmapped 0 -> 0 [] -> []; inert hole bridges 2 -> 2; skipped 0 -> 0; reliable True -> True; warnings 0 -> 0
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
== k28_8_blasts: 4 classes, inert 1 -> 1; image bridges 3 -> 3; iterated elements 13 -> 13; unmapped 0 -> 0 [] -> []; inert hole bridges 1 -> 1; skipped 0 -> 0; reliable True -> True; warnings 4 -> 4
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
== nested_inner6: 47 classes, inert 2 -> 3; image bridges 7 -> 7; iterated elements 105 -> 105; unmapped 1 -> 0 [[46, 47]] -> []; inert hole bridges 2 -> 3; skipped 0 -> 0; reliable True -> True; warnings 0 -> 0
    verdict {B:1.0/right[46,7], B:1.0/right[8,47]}: ae active -> w inert (zeroth-hole class: True, word trivial)
    equal: rules, kept/image/skipped, iterated partition, refined pairs, is_reliable
```

**Reading:**
- On all 15 builds the class sets, letter-free rules, minimal trellis, iterated partition, refined pairs and `is_reliable` are equal.
- The only differences are the five verdicts (`g`; `d`, `bj`; `ab`; `ae`), all zeroth-hole classes with a trivial walk.
- The `unmapped` deltas match the plan: −1 on inversion_two_blasts, p3_6_blasts and nested_inner6, −2 on heteroclinic. The moved bridges are exactly those that entered `inert_hole_bridge_ids`.
- No WARNING appears on the branch that is absent on the base. k28 x8 logs its same four on both.

## Rubric 2 — the two new laws on base and branch

Script: `scratchpad/b5cmp/laws2.py`. In the base worktree it ran with the branch's `laws.py` copied in, restored afterwards. It runs both laws on all 11 law cases and also judges each half separately, every hole and every class:
- empty-word half: an itinerary present and no symbol;
- kind half: `kind == "inert"`;
- transience: its own DFS.

**Base** (`laws_base.txt`):
- Nine law cases pass both laws and all three halves: k10, k28_one_blast, k28_two_blasts, p3, nested, nested_inner_blasts, inversion, inversion_one_blast and inversion_period2.
- inversion_two_blasts:
  ```
  inert_is_transient FAIL ('g', 'active', 'A:p1@0.0/L#3 <-> A:p1@0.0/L#5'); zeroth_hole_classes_have_empty_words FAIL ('g', 'A:p1@0.0/L#3 <-> A:p1@0.0/L#5')
  empty-word half failures: none
  kind half failures: [('g', '{A:0.0/left[23,2], A:0.0/left[3,24]}')]
  transience failures: [('g', 'active', '{A:0.0/left[23,2], A:0.0/left[3,24]}')]
  ```
- heteroclinic:
  ```
  inert_is_transient FAIL ('d', 'active', 'A:p1@0.0/L#2 <-> A:p1@0.0/L#4'); zeroth_hole_classes_have_empty_words FAIL ('d', 'A:p1@0.0/L#2 <-> A:p1@0.0/L#4')
  empty-word half failures: none
  kind half failures: [('d', '{A:0.0/left[114,4], A:0.0/left[92,113]}'), ('bj', '{A:0.0/left[92,113], A:0.0/left[None,93]}')]
  transience failures: [('d', 'active', '{A:0.0/left[114,4], A:0.0/left[92,113]}'), ('bj', 'active', '{A:0.0/left[92,113], A:0.0/left[None,93]}')]
  ```
- So on the base the empty-word half passes on every law case, and both the kind half and transience fail on exactly inversion_two_blasts (`g`) and heteroclinic (`d`, `bj`).

**Branch** (`laws_branch.txt`): both laws pass on all 11 law cases, with no failures in any half. The returned counts:

| Law case | `inert_is_transient` | `zeroth_hole_classes_have_empty_words` |
|---|---|---|
| k10 | 2 | 1 |
| k28_one_blast | 3 | 2 |
| k28_two_blasts | 4 | 1 |
| p3 | 5 | 2 |
| nested | 9 | 3 |
| nested_inner_blasts | 9 | 3 |
| inversion | 4 | 2 |
| inversion_one_blast | 7 | 1 |
| inversion_two_blasts | 8 | 2 |
| inversion_period2 | 6 | 2 |
| heteroclinic | 75 | 4 |

The law tier also runs them in the suite: `tests/invariants/test_law_symbolic.py` passes. Both `inert_classes` unit cases pass.

## Behaviour changes

On the 15 builds, only those in section 5:
- **Five verdicts** go from active to inert: inversion_two_blasts `g`; heteroclinic `d` and `bj`; p3_6_blasts `ab`; nested_inner6 `ae`.
  - In a fresh session they take inert letters, respectively `v`; `v` and `y`; `v`; `w`. Active letters assigned after them shift down.
  - They leave the transition graph and matrix.
- **`unmapped` shrinks** as listed above, and those bridges move to `inert_hole_bridge_ids`.

Elsewhere, all reported, none pinned:
1. **The `p3_partitioned` fixture** (unblasted nested, `p1_area_cutoff=1e-4`; not one of the 15 builds). It changes in the same way:
   - `unmapped` [(4,5), (5,6)] becomes [];
   - both bridges carry zeroth holes and move to `inert_hole_bridge_ids` (2 → 4);
   - the class of (5,6), trivial word, goes from active `c` to inert `u`;
   - the class of (4,5) is UNRESOLVED (its singleton landing has no registered chain), so it stays active by the convention, letter `b`;
   - later active letters shift (`d,e,f` → `c,d,e`).

   This fixture is the one class on any build where a zeroth-hole class has no word. The new law runs on law cases only, and none has this, so this is not a counterexample to the theorem. It is a case its hypothesis (a resolvable walk) excludes.
2. **Log text.** While walking, a class has no letter yet, so `_resolve_itinerary`'s "class %s unresolved" WARNINGs now name the class by its label (`A:p3@1.0/R#2 <-> …`) instead of its letter. They are the same records, worded differently. They fire on no law case or deep build; they do fire on the `p3_partitioned` fixture.
3. **PK-8's WARNING is gone.** It fired on no build.
4. **API.**
   - `BridgeClassEntry.image_classes`, `.image_loops`, `.unresolved` and `.has_image_evidence` are gone.
   - `entry.inert` is `None` until the dynamics sets it, and on evidence entries outside the alphabet.
   - `BridgeClassEntry.describe()` takes no `names` argument, prints no images, and prints `undecided` for `inert is None`.
   - `BridgeClassTable.active` is every entry before the dynamics runs.
   - `MinimalTrellis.inert_hole_bridge_ids` now means "carries a zeroth hole".
   - `session.bridge_classes()` / `Topology.classes` runs the dynamics.
   - `symbolic_dynamics` now sets `inert` and `letter` on the alphabet's entries.
   - `inert_classes` is new and public in `SymbolicDynamics`. It is not re-exported from `tanglepack.topology`, because `__init__.py` is not in the file list.
5. **One detail beyond the plan's wording, in `_active_letters`.** An active entry that already carries a letter keeps it (`entry.letter or letter_for(...) or letter(i)`), as `_active_letters` did at the base. Nothing in the session pre-letters any more. Hand-built test tables do pre-letter, and they rely on it.

## Test edits beyond the plan's letter (all inside the nine files)

- **`test_minimal_trellis.py`**
  - The renamed test could not keep `assert minimal.unmapped`. On its fixture, `p3_partitioned`, every base `unmapped` bridge carried a zeroth hole, so `unmapped` is now empty. The planner's "would most likely still pass" did not hold.
  - The test now asserts:
    1. no unmapped bridge carries a zeroth hole (rubric 7);
    2. the fixture does have hole bridges with no registered image;
    3. re-running `minimal_trellis` on the same holes with `pair` stripped (`dataclasses.replace`, i.e. read as propagated) lists exactly those bridges as `unmapped` and none as inert. This keeps the "unmapped" path exercised on real data.
  - The test is renamed `test_a_hole_bridge_without_a_registered_image_is_unmapped`, and the module docstring is updated.
- **`test_symbolic_dynamics.py`**
  - Deleted `test_inert_class_naming_an_active_class_warns_and_stays_out`, which tested the deleted PK-8 WARNING.
  - Deleted `test_inert_class_mapping_onto_inert_classes_is_silent`. Its premise (an inert class whose word is itself) is impossible under transience: a self-cycle is active.
  - Rewrote `test_a_word_naming_a_spanned_class_closes_the_alphabet`. The closure class's verdict now comes from its own word, not the evidence flag, so the test runs once with an empty word (inert `u`, `letter_for` not asked) and once with a word naming the self-mapping `a` (active, lettered `q` by `letter_for`).
  - `test_a_pair_no_bridge_spans_is_a_flagged_class` now walks the member-less class onto `a`, so it stays active and still exercises `letter_for`. With the old walk `L_3 L_3` it would be inert.
  - Added the plan's two cases.
- **`test_bridge_class.py`**
  - Besides deleting `test_resolve_inertness_fixed_point` with `_entry` (and its now-unused `_ref`, `_class`, `_INERT_FP` and the `BridgeClassEntry` import), three tests read the evidence verdict and were adapted:
    - the two loop tests select entries with loops instead of `table.inert`;
    - the orphan test drops its "former class is now active" line;
    - the lookups test takes `table.active[0]` only, the same entry as before.
  - The fold-born test drops its `image_loops` assertion (rubric 5).
- **`direct.py`**: the pre-lettering loop and the `BridgeAlphabet` import are deleted; `classes()` returns the unlettered alphabet.

## Rubric, my reading

1. **Pass.** 1192 passed, 38 skipped, 5 deselected, 1 xfailed. `-rxX` lists only the orientation_reversing KNOWN_ISSUE.
2. **Pass.** On the base, with the branch's `laws.py` copied in, the empty-word half passes on all 11 law cases. `inert_is_transient` and the kind half fail on exactly inversion_two_blasts (`g` = `{A:0.0/left[23,2], A:0.0/left[3,24]}`) and heteroclinic (`d` = `{A:0.0/left[114,4], A:0.0/left[92,113]}`, `bj` = `{A:0.0/left[92,113], A:0.0/left[None,93]}`). Everything passes on the branch, and both unit cases pass.
3. **Pass.** Exactly the nine files.
4. **Pass.** See the output above: equal everywhere except the five verdicts and the `unmapped` deltas, with no new WARNING.
5. **Pass.** `grep -rn "has_image_evidence\|_resolve_inertness\|_collect_image_evidence\|image_classes\|image_loops" src tests` is empty.
6. **Pass.** `grep -n "\.active\|\.inert" src/tanglepack/topology/MinimalTrellis.py` hits only lines 224 and 238 (`self.inert_hole_bridge_ids` in `summary`/`describe`), not `minimal_trellis`. The updated `minimal_trellis_bridges` passes on every law case.
7. **Pass.** The `table.entry_of(bid).active` assertion is gone. The renamed test asserts no unmapped bridge carries a zeroth hole and passes. The fixture change is explained above.
8. **Pass.**
   - `Topology.classes` has no `letter_for` call; the only call is in `dynamics`.
   - Fresh sessions (`scratchpad/r8.py`): `session.bridge_classes()` has no `inert is None` entry on k10 (2), inversion_two_blasts (8) or heteroclinic (53), and letters exactly the active ones.
   - `only_active_classes_lettered` and `test_session_bridge_classes.py` pass unchanged.
9. **Pass.**
   - No PK-8 WARNING.
   - `inert_classes(` is called once, at `SymbolicDynamics.py:1471`: after the closure loop, before translation.
   - `add` assigns no letter.
10. **Pass, with a caveat on how it fails.** The `a -> b, b unresolved` case asserts both are active. With `word is not None` replaced by `True` (checked by mutating the file and restoring it), the case fails: `TypeError` at `SymbolicDynamics.py:516`, because `all(...)` iterates the `None` word. The other case still passes. So it does not fail on the assertion; the mutated code raises before the assertion is reached.
11. **Pass.** The changed verdicts per build are the five of section 1. The `unmapped` deltas match section 5. The differences off the 15 builds (the `p3_partitioned` fixture, the log wording, the letter retention in `_active_letters`) are explained above, not pinned.

## Artifacts of the runs

All under `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/`:
- `b5cmp/dump.py` and `b5cmp/compare.py`: rubric 4. Outputs are in `b5cmp/compare_output.txt`, with per-build JSON in `b5cmp/{base,branch}/*.json`.
- `b5cmp/laws2.py`: rubric 2. Outputs are in `b5cmp/laws_base.txt` and `b5cmp/laws_branch.txt`.
- `r8.py`: rubric 8.
- `suite_branch.txt`: the full suite run.

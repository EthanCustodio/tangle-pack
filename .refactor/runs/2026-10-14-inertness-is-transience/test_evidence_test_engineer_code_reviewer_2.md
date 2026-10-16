# Test evidence: B5 `inertness-is-transience` (test engineer, round 2)

- **Plan:** `artifacts/plan_planner_human_2.md`.
- **Coder report:** `artifacts/change_coder_test_engineer.t1.md`.
- **Base:** `7d3a294`, detached worktree `../tangle-pack-inertness-is-transience-base`. It is clean at the end; every copy I made into it was restored with `git checkout`.
- **Branch:** `refactor/inertness-is-transience`, worktree `../tangle-pack-inertness-is-transience`. The coder's head was `b9145d0`. I added one commit, so the head is now **`3fb310a`**. Every break-it mutation was restored with `git checkout`, and the tree is clean.
- **Scripts and outputs:** all are mine, written independently of the coder's `b5cmp/` scripts. They live in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/te_b5/`:
  - `dump.py`, `run.sh`, `compare.py` → `compare_output.txt` (rubric 4);
  - `laws_run.py`, `laws.sh` → `laws_output.txt` (rubric 2);
  - `r8.py` → `r8_output.txt` (rubric 8);
  - `mut1/`, `mut2/` with `mut{1,2}_compare.txt` (break-it checks 1 and 2);
  - `suite_branch.txt` (the full suite).

## Verdict summary

| Rubric | Result | How I checked it |
|---|---|---|
| 1 full suite | **PASS**: 1193 passed, 38 skipped, 5 deselected, 1 xfailed. `-rxX` lists only `orientation_reversing` (C3, `cases.KNOWN_ISSUES`). | I re-ran it on the branch, my new test included. |
| 2 pinning on both commits | **PASS** | I re-ran it with the branch's `laws.py` copied into the base, through my own script and through the real law tier. |
| 3 nine files | **PASS**: `git diff --name-only 7d3a294..` lists exactly the nine files, my commit included (`test_symbolic_dynamics.py` is one of them). | I re-ran it. |
| 4 no behaviour change beyond §5 | **PASS** on all 15 builds | My own dump and compare, stricter than the rubric (see below). |
| 5 grep for the deleted API | **PASS**: no output (grep exit 1). | I re-ran it. |
| 6 `minimal_trellis` reads no verdict; the inert set equals the zeroth set | **PASS** | grep, the law tier, and the 15-build dump. |
| 7 test_minimal_trellis | **PASS**: the old assertion is gone, and the renamed test passes in the suite. It also fails under break 1. | Suite plus break 1. |
| 8 fresh-session verdicts | **PASS** on all 11 law cases. The facade test passes and is unchanged. **New test added.** | `r8.py`, plus the new unit test. |
| 9 one `inert_classes` call; `add` letters nothing | **PASS** | grep and reading the code. |
| 10 unresolved case | **PASS**, and the coder's caveat is settled (see break 3). | Mutation. |
| 11 changed verdicts and `unmapped` deltas | **PASS**: exactly the five classes and the deltas of §5. | `compare_output.txt` |

## Rubric 4: base vs branch on all 15 builds

The command, run once per (worktree, build) pair with `xargs -P 16`:

```
cd <worktree> && PYTHONPATH=src:tests MPLBACKEND=Agg env/bin/python te_b5/dump.py <build> te_b5/<base|branch>/<build>.json
env/bin/python te_b5/compare.py
```

The 15 builds:
- the 11 `cases.LAW_CASES`;
- `p3_6_blasts` and `p3_16` (`cases.BUILDERS`);
- `k28_8_blasts` = `cases.build_k28(blasts=8)`;
- `nested_inner6` = `cases.build_nested(inner_blasts=6)`.

**Class identity.** A class is keyed by the sorted pair of its two homotopy elements `{fp:o.b/side[lo_id,hi_id]}`. Letters appear only in the VERDICT lines.

**Compared, letter-free.** This goes further than the rubric:
- the class set;
- each class's class-level word, as `(identity, direction)` tokens;
- each class's full iterated itinerary, by iterated-interval identity;
- the walk status and the itinerary source;
- the refined children's iterated pairs;
- `kept_bridge_ids`, `image_bridge_ids` and `skipped_pairs`;
- every iterated element, with its interval identity AND its name (`P1` with names), plus `iterated.signature()` (FixedPoint addresses normalised);
- the refined transition matrix on letter-free nodes `(class identity, refined index)`. Branch nodes must equal base nodes minus the flipped classes' nodes, and the entries must be equal on the remaining nodes;
- `is_reliable`;
- every WARNING-or-higher record on any logger, as a multiset.

Verbatim `compare_output.txt`:

```
== k10: 2 classes; inert 1 -> 1; unmapped [] -> []; inert_hole_bridge_ids [[5, 6]] -> [[5, 6]] (zeroth-hole bridges [[5, 6]]); matrix nodes 2 -> 2; reliable True -> True; WARNINGs 0 -> 0
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
== k28_one_blast: 3 classes; inert 2 -> 2; unmapped [] -> []; inert_hole_bridge_ids [[3, 4], [4, 5]] -> [[3, 4], [4, 5]] (zeroth-hole bridges [[3, 4], [4, 5]]); matrix nodes 1 -> 1; reliable True -> True; WARNINGs 0 -> 0
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
== k28_two_blasts: 4 classes; inert 1 -> 1; unmapped [] -> []; inert_hole_bridge_ids [[9, 10]] -> [[9, 10]] (zeroth-hole bridges [[9, 10]]); matrix nodes 4 -> 4; reliable True -> True; WARNINGs 0 -> 0
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
== p3: 5 classes; inert 2 -> 2; unmapped [] -> []; inert_hole_bridge_ids [[3, 4], [4, 5]] -> [[3, 4], [4, 5]] (zeroth-hole bridges [[3, 4], [4, 5]]); matrix nodes 3 -> 3; reliable True -> True; WARNINGs 0 -> 0
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
== nested: 9 classes; inert 3 -> 3; unmapped [] -> []; inert_hole_bridge_ids [[7, 8], [8, 9], [37, 38]] -> [[7, 8], [8, 9], [37, 38]] (zeroth-hole bridges [[7, 8], [8, 9], [37, 38]]); matrix nodes 7 -> 7; reliable True -> True; WARNINGs 0 -> 0
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
== nested_inner_blasts: 9 classes; inert 3 -> 3; unmapped [] -> []; inert_hole_bridge_ids [[7, 8], [8, 9], [41, 42]] -> [[7, 8], [8, 9], [41, 42]] (zeroth-hole bridges [[7, 8], [8, 9], [41, 42]]); matrix nodes 7 -> 7; reliable True -> True; WARNINGs 0 -> 0
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
== inversion: 4 classes; inert 2 -> 2; unmapped [] -> []; inert_hole_bridge_ids [[1, 2], [2, 3]] -> [[1, 2], [2, 3]] (zeroth-hole bridges [[1, 2], [2, 3]]); matrix nodes 2 -> 2; reliable True -> True; WARNINGs 0 -> 0
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
== inversion_one_blast: 7 classes; inert 1 -> 1; unmapped [] -> []; inert_hole_bridge_ids [[26, 27]] -> [[26, 27]] (zeroth-hole bridges [[26, 27]]); matrix nodes 7 -> 7; reliable True -> True; WARNINGs 0 -> 0
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
== inversion_two_blasts: 8 classes; inert 1 -> 2; unmapped [[23, 24]] -> []; inert_hole_bridge_ids [[26, 27]] -> [[23, 24], [26, 27]] (zeroth-hole bridges [[23, 24], [26, 27]]); matrix nodes 8 -> 7; reliable True -> True; WARNINGs 0 -> 0
    VERDICT {A:0.0/left[23,2] A:0.0/left[3,24]}: g active -> v inert (zeroth-hole class True, walk trivial, word [])
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
== inversion_period2: 6 classes; inert 2 -> 2; unmapped [] -> []; inert_hole_bridge_ids [[13, 14], [14, 15]] -> [[13, 14], [14, 15]] (zeroth-hole bridges [[13, 14], [14, 15]]); matrix nodes 4 -> 4; reliable True -> True; WARNINGs 0 -> 0
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
== heteroclinic: 75 classes; inert 3 -> 5; unmapped [[92, 93], [113, 114]] -> []; inert_hole_bridge_ids [[77, 78], [78, 79]] -> [[77, 78], [78, 79], [92, 93], [113, 114]] (zeroth-hole bridges [[77, 78], [78, 79], [92, 93], [113, 114]]); matrix nodes 76 -> 74; reliable True -> True; WARNINGs 0 -> 0
    VERDICT {A:0.0/left[114,4] A:0.0/left[92,113]}: d active -> v inert (zeroth-hole class True, walk trivial, word [])
    VERDICT {A:0.0/left[92,113] A:0.0/left[None,93]}: bj active -> y inert (zeroth-hole class True, walk trivial, word [])
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
== p3_6_blasts: 43 classes; inert 1 -> 2; unmapped [[42, 43]] -> []; inert_hole_bridge_ids [[41, 42]] -> [[41, 42], [42, 43]] (zeroth-hole bridges [[41, 42], [42, 43]]); matrix nodes 43 -> 42; reliable True -> True; WARNINGs 0 -> 0
    VERDICT {A:1.0/right[4,43] A:1.0/right[42,3]}: ab active -> v inert (zeroth-hole class True, walk trivial, word [])
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
== p3_16: 43 classes; inert 2 -> 2; unmapped [] -> []; inert_hole_bridge_ids [[53, 54], [54, 55]] -> [[53, 54], [54, 55]] (zeroth-hole bridges [[53, 54], [54, 55]]); matrix nodes 42 -> 42; reliable True -> True; WARNINGs 0 -> 0
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
== k28_8_blasts: 4 classes; inert 1 -> 1; unmapped [] -> []; inert_hole_bridge_ids [[9, 10]] -> [[9, 10]] (zeroth-hole bridges [[9, 10]]); matrix nodes 4 -> 4; reliable True -> True; WARNINGs 4 -> 4
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
== nested_inner6: 47 classes; inert 2 -> 3; unmapped [[46, 47]] -> []; inert_hole_bridge_ids [[45, 46], [53, 54]] -> [[45, 46], [46, 47], [53, 54]] (zeroth-hole bridges [[45, 46], [46, 47], [53, 54]]); matrix nodes 47 -> 46; reliable True -> True; WARNINGs 0 -> 0
    VERDICT {B:1.0/right[46,7] B:1.0/right[8,47]}: ae active -> w inert (zeroth-hole class True, walk trivial, word [])
    equal: class set, words, itineraries, refined pairs, walk status, kept/image/skipped, iterated elements+names, matrix on remaining nodes, is_reliable
```

**Reading.**
- All 15 builds are equal on everything compared. The only differences are the five verdicts and the `unmapped` deltas:
  - inversion_two_blasts `g`;
  - heteroclinic `d` and `bj`;
  - p3_6_blasts `ab`;
  - nested_inner6 `ae`.
- Each of the five is a zeroth-hole class with a `trivial` walk and an empty word.
- `unmapped` drops by 1 on inversion_two_blasts, p3_6_blasts and nested_inner6, and by 2 on heteroclinic. Each moved bridge reappears in `inert_hole_bridge_ids`.
- On every build the branch's `inert_hole_bridge_ids` equals the zeroth-hole bridge set (rubric 6 at runtime, all 15 builds).
- The matrix loses exactly the flipped classes' nodes (8→7, 76→74, 43→42, 47→46), and every remaining entry is unchanged.
- No WARNING is new. k28_8_blasts logs the same four records on both commits: two numerics resolution flags, the same-stability count, and one cross-side pair flag.

## Rubric 2: the two new laws on the base and the branch

**Base.** I copied the branch's `tests/helpers/laws.py` into the base worktree (`cmp` identical), ran the checks there, then ran `git checkout tests/helpers/laws.py`; the base was clean again. `laws_run.py` calls `laws.check_inert_is_transient` and `laws.check_zeroth_hole_classes_have_empty_words`, then judges each half separately with its own DFS.

```
cd <worktree> && PYTHONPATH=src:tests MPLBACKEND=Agg env/bin/python te_b5/laws_run.py <law case>
```

Base, the failing cases verbatim. All other law cases pass both laws and all three halves:
- k10, k28_one_blast, k28_two_blasts, p3, nested, nested_inner_blasts;
- inversion, inversion_one_blast, inversion_period2.

```
== inversion_two_blasts
  check_inert_is_transient FAIL ('g', 'active', 'A:p1@0.0/L#3 <-> A:p1@0.0/L#5')
  check_zeroth_hole_classes_have_empty_words FAIL ('g', 'A:p1@0.0/L#3 <-> A:p1@0.0/L#5')
  empty-word half failures: none
  kind half failures: [('g', '{A:0.0/left[23,2] A:0.0/left[3,24]}')]
  transience failures (own DFS): [('g', 'active', '{A:0.0/left[23,2] A:0.0/left[3,24]}')]
== heteroclinic
  check_inert_is_transient FAIL ('d', 'active', 'A:p1@0.0/L#2 <-> A:p1@0.0/L#4')
  check_zeroth_hole_classes_have_empty_words FAIL ('d', 'A:p1@0.0/L#2 <-> A:p1@0.0/L#4')
  empty-word half failures: none
  kind half failures: [('d', '{A:0.0/left[114,4] A:0.0/left[92,113]}'), ('bj', '{A:0.0/left[92,113] A:0.0/left[None,93]}')]
  transience failures (own DFS): [('d', 'active', '{A:0.0/left[114,4] A:0.0/left[92,113]}'), ('bj', 'active', '{A:0.0/left[92,113] A:0.0/left[None,93]}')]
```

On the base:
- The empty-word half passes on all 11 law cases.
- `inert_is_transient` and the kind half fail on exactly inversion_two_blasts (`g`) and heteroclinic (`d`, `bj`).
- Matched by element identity, these are the same classes the rubric-4 compare flips.

The same check through the real law tier on the base, with the branch's `laws.py` copied in:
```
env/bin/python -m pytest -q tests/invariants/test_law_symbolic.py tests/invariants/test_law_iterated.py tests/invariants/test_law_classes.py
FAILED tests/invariants/test_law_symbolic.py::test_symbolic_laws[inversion_two_blasts]   (inert_is_transient, zeroth_hole_classes_have_empty_words)
FAILED tests/invariants/test_law_iterated.py::test_iterated_laws[inversion_two_blasts]   (minimal_trellis_bridges: inert set != zeroth set)
FAILED tests/invariants/test_law_symbolic.py::test_symbolic_laws[heteroclinic]
FAILED tests/invariants/test_law_iterated.py::test_iterated_laws[heteroclinic]
4 failed, 29 passed, 32 skipped
```

The rewritten `anchor_bridge_class_leads_its_tangle` (classes layer) passes on the base.

**Branch** (`laws_output.txt`): both laws and every half pass on all 11 law cases. The counts returned match the coder's table, from k10 2/1 to heteroclinic 75/4. Both `inert_classes` unit cases pass.

## Rubric 6

- `grep -n "\.active\|\.inert" src/tanglepack/topology/MinimalTrellis.py` hits only `:224` and `:238`, both inside `MinimalTrellis.summary` / `describe`. `minimal_trellis` starts at `:268`, so there is no hit inside the function.
- The function calls `classes.entry_of(bid)` only for the unclassed WARNING and discards the result.
- `minimal_trellis_bridges` passes on the 11 law cases in the suite.
- On all 15 builds, `inert_hole_bridge_ids == zeroth-hole bridges` (see the compare above).

## Rubric 8

`r8.py` builds each law case fresh and asserts that `"dynamics" not in vars(session.topology())`. It then calls `session.bridge_classes()` and reads the entries.

```
heteroclinic: dynamics cached before the call: False; 53 entries; inert None: none; letter iff active: True; inert == dynamics kind: True
inversion: ... 4 entries; inert None: none; letter iff active: True; inert == dynamics kind: True
inversion_one_blast: ... 7 entries; inert None: none; ...True; ...True
inversion_period2: ... 6 entries; inert None: none; ...True; ...True
inversion_two_blasts: ... 8 entries; inert None: none; ...True; ...True
k10: ... 2 entries; inert None: none; ...True; ...True
k28_one_blast: ... 3; k28_two_blasts: ... 4; nested: ... 9; nested_inner_blasts: ... 9; p3: ... 5  (all: inert None: none, letter iff active, inert == kind)
```

(The full lines are in `r8_output.txt`.)

- `tests/facade/test_session_bridge_classes.py`: 9 passed. `git diff 7d3a294 -- <it>` is empty.
- `Topology.classes` is `self.dynamics; return self._alphabet`. The only `letter_for` call is in `dynamics` (`:190`).

## Rubric 9

- `grep -rn "inert_classes(" src/` finds one call site, `SymbolicDynamics.py:1471`. It sits after the `for entry in entries` closure loop (`:1437-1468`) and before the `translate_itinerary(` call (`:1494`).
- `grep "PK-8\|names a non-inert"` in `SymbolicDynamics.py` finds nothing.
- `add` (`:1413-1429`) copies the entry with `letter=None, inert=None` or builds a member-less one. It touches only `listed` and `closure`, and assigns no letter.

## New test (committed, `3fb310a`)

`tests/unit/topology/test_symbolic_dynamics.py::test_the_session_alphabet_carries_the_dynamics_verdicts(k10_partitioned)`, in the existing "Real data" section.

**What it does.** On a fresh session, where the dynamics is not cached, it calls `session.bridge_classes()` and snapshots each entry's `(class, inert, letter)` at once. Only then does it read `session.symbolic_dynamics()`. It asserts:
- some class is inert;
- every snapshot verdict equals the dynamics' `kind`;
- `letter is None` exactly when the class is inert.

**What it guards.** Rubric 8, the order "verdict first, letters after", without the law tier's order dependence. Under break 5 below, only `test_session_bridge_classes.py::test_active_class_lies_in_the_zone...` failed in the facade + law + unit + regression tiers. The law tier missed it, because an earlier law (`anchor_bridge_class_leads_its_tangle`) runs the dynamics on the shared case before `only_active_classes_lettered` reads it.

**Base.** I copied the file into the base worktree and ran only this test: **1 passed**. The verdicts were already present at alphabet time on the base; the file was restored.

**Branch.** Passes.

**Teeth.** It fails under breaks 4 and 5 below. A first draft asserted on the live entries after `symbolic_dynamics()` and passed under break 5, because the dynamics mutates the same entry objects. I fixed it with the snapshot before committing.

## Break-it checks (on the branch, each restored with `git checkout <file>`)

| # | Mutation | Result |
|---|---|---|
| 1 | `MinimalTrellis.py`: drop `if frozenset(bid) in zeroth: inert.append(bid); continue` (zeroth-hole bridges get mapped too) | **Caught.** `test_law_iterated`: `minimal_trellis_bridges` fails on all 11 law cases (`assert set() == {frozenset({5, 6})}` etc.). Four `iterated_cut_provenance` set-asides (inversion, inversion_period2, k28_one_blast, p3) fail because the law now applies. The renamed `test_a_hole_bridge_without_a_registered_image_is_unmapped` fails: 16 failed. **P1 comparison** (`mut1_compare.txt`, base vs mutant): kept/image/iterated partition/itineraries DIFFER on k28_one_blast, p3 and heteroclinic (heteroclinic's refined pairs and matrix too). skipped_pairs DIFFER on inversion_two_blasts, with a new WARNING on it and on k28_one_blast. On k10 and k28_two_blasts the zeroth bridge has no registered image, so only `unmapped` moves there. |
| 2 | `MinimalTrellis.py`: replace the zeroth test with an unconditional `inert.append(bid); continue` (skip ALL hole bridges) | **Caught.** The image chains vanish: image bridges 3→0 (k10, k28_two_blasts, nested, inversion_one_blast) and 8→0 (heteroclinic). Words change on 1/2, 1/4, 1/9, 1/7 and 21/75 classes. `is_reliable` goes True→False on all five; heteroclinic gains a class (76). Suite slice: 21 failed (iterated and symbolic laws on 11 cases, three `test_minimal_trellis.py` tests). |
| 3 | `SymbolicDynamics.inert_classes`: `word is not None and all(... for token in word)` → `all(... for token in word or ())` (an unresolved word read as empty) | **Caught on the assertion.** `test_an_unresolved_class_is_active_and_so_is_every_class_reaching_it` fails with `AssertionError: assert {BridgeClass(...), BridgeClass(...)} == set()`; the other case passes. |
| 3' | The coder's mutation, `word is not None` → `True` | Confirmed as the coder said: the same test fails with `TypeError: 'NoneType' object is not iterable` at `SymbolicDynamics.py:516`, before the assertion. **This settles the caveat.** That mutation cannot reach the assertion, because the fixed point iterates the `None` word. The semantically meaningful break, mutation 3, fails on the assertion. |
| 4 | `Topology._alphabet`: re-add `if entry.active: entry.letter = self.letter_for(entry.bridge_class)` (letter before the verdict) | **Caught.** `only_active_classes_lettered` fails on all 11 law cases (`assert ('b' is not None) == False`: the inert class kept the letter `b`). `test_session_bridge_classes.py::test_active_class_lies_in_the_zone_after_trimming_at_the_image_pip` fails. The new test fails (`('b' is None) is True`). Slice: 12 failed. |
| 5 | `Topology.classes`: drop `self.dynamics` (the alphabet is returned before any verdict) | Extra check. **Caught** by the facade test above and by the new test (`AssertionError: k=10 has an inert class`). The law tier did not catch it (see "New test"). |

## Full suite on the branch (`3fb310a`)

```
MPLBACKEND=Agg env/bin/python -m pytest -q -rxX
XFAIL tests/invariants/test_law_case_sanity.py::test_orientation_reversing_case_builds[orientation_reversing] - C3: INV-11: det J < 0 is out of scope (C3); FixedPoint.set_k_value rejects a saddle with exactly one negative eigenvalue
1193 passed, 38 skipped, 5 deselected, 1 xfailed in 198.65s (0:03:18)
```

The coder had 1192; the +1 is the new test. No existing test was weakened or deleted by me.

## What I did not test, and why

- **The `p3_partitioned` fixture** (unblasted nested, not one of the 15 builds). The coder reports a zeroth-hole class there that stays active because it is UNRESOLVED. That is the one place where the zeroth-hole rule (map-skip) and the verdict disagree: the bridge is not mapped, yet its class is active by the unresolved convention. I did not pin it: plan §5 says report, not pin, and the law runs on law cases only. The reviewer should note it as the hypothesis boundary of the author's theorem (a resolvable walk).
- **Letter shifts in a fresh session** (active letters after a flipped class move down by one) are not pinned. The house rule is no letter pins; the compare shows them in the VERDICT lines only.
- **Rubric 4's WARNING comparison** covers the build plus `topology()` and the dynamics. It does not cover plotting or the scripts.

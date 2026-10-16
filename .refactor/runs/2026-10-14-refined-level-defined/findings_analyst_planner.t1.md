# Findings: B7 `refined-level-defined` (analyst, iteration 1)

**The tree is dirty.** `main` @ `3fb310a`. These files have uncommitted changes: `CLAUDE.md`, `.refactor/decisions.md`, `.refactor/knowledge/*.md`, `docs/alignment_plan_2026-10-09.md` and `docs/code_fidelity_audit_2026-10-09.md`. There is also the staged rename of `artifacts/request_origin_analyst.md`, and these paths are untracked: `artifacts/`, `.refactor/runs/*`, `figures/henon_period3_*` and `scripts/henon_period3_*`. The change branches from `HEAD` and must not touch any of them. The plan text I compare against below is the working-tree `docs/alignment_plan_2026-10-09.md`.

**Kind of request.** This is a definitional refactor plus two new laws. The arithmetic must not change. The targets are:

- `SymbolicDynamics.refined_rules` and its documentation;
- the odd-count WARNING in `IteratedHomotopyPartition._empty_stretches`;
- two new laws in `tests/helpers/laws.py` (layer `symbolic`) and a series run.

I did not re-plan the slice.

## The answer, up front

1. **Theorem T fails on one series as the plan states it, so Rubric 5 is triggered.** The run is `k32_series` at blasts 6 and 7. In class `A:p1@0.0/R#7 <-> A:p1@0.0/R#9`, hole-bridge member `(6, 7)` has a different `P1` itinerary from its siblings `(54, 55)`, `(64, 65)` and `(96, 97)`. The difference is at the occurrence level; the class-level word is the same. Following Q15, this goes to the author before step 4 (the series half of `members_share_the_P1_itinerary`) is planned.
   - The cause is the **chain-order chord pairing on the base branch** (`PartitionFamily.py:957`).
   - On that lobe, chain-order pairing disagrees with the plane. Point-in-polygon on the lobe image shows that the stretches inside `f(L)` are `(11,10)` and `(9,12)`. Chain-order pairing declares `(9,10)` and `(11,12)` empty instead.
   - I ran a scratch monkeypatch that pairs in stable order on the base branch as well. With it, T holds on every build and step, and `M`, `M_1`, node names, class counts and `is_reliable` are unchanged everywhere.
   - That is a change to `from_minimal`'s cuts, which B7's not-doing list forbids. Details are in §2b.
2. **A second, literal "failure" is incomplete evidence, not a counterexample.** It occurs on k28 x8 and on `k28_series` steps 6–8. One member's registered image chain runs past the computed stable branch: two consecutive crossings sit on one branch with one sign, which is the existing C11 "stable end" criterion. The law needs to exclude such chains, as `_compare_evidence` already does for word verdicts. With that exclusion, those four builds and steps pass.
3. **On the eleven law cases, T holds with 0 disagreements.** It is applicable on k10 (1 class), inversion_two_blasts (2) and heteroclinic (7), which is 10 classes, not the plan's 5. Period-k classes are tested: heteroclinic's 7 and p3_16's 4 are all on the period-3 tangle, as are `p3_series` blasts 7–9. Inversion classes (two branches) are tested too.
4. **The in-splitting identities hold on all 62 builds and series steps** (15 builds plus 47 series steps), with 0 violations.
5. **The odd-count WARNING never fires** on any build or step at `HEAD`. `incomplete_footprints` will therefore be empty on the whole corpus; only the hand-built unit test populates it.

## 1. The code at HEAD (verified by reading)

| What | Where at `HEAD` | Matches the request? |
|---|---|---|
| Module Dev Notes bullet "every refined child inherits the parent's word" | `SymbolicDynamics.py:47-58` (bullet starts at 47; "inherits" is on :50) | request said :48-52; close |
| `translate_itinerary` | `SymbolicDynamics.py:525` (loops dropped at :577-579; `occurrence` swapped for direction -1 at :594) | yes |
| `refine` | `:606` (occurrences from ALL classes' words, inert ones included, :625-630) | yes |
| `_refined_symbol` | `:656` (an occurrence matching no child returns the bare symbol, :667-669) | yes |
| `SymbolicDynamics` docstring of `refined_rules` | `:758-760` ("all with the same word") | — |
| `refined_rules` construction | `:796-809` (`refined_word` :802, the per-child copy :805-807) | yes |
| `match_members` | `:813` | yes |
| `_edges` | `:940`; it skips inert tokens and uses `symbol.base` as the target when refined | — |
| `transition_graph` | `:960` (asserts every target is a node) | yes |
| `transition_matrix` | `:995` (silently drops a target that is not a node, `:1012`) | yes |
| `trellis_itinerary` | `DualWalk.py:1301` | yes |
| odd-count WARNING | `PartitionFamily.py:959-966` (`if len(on_branch) % 2:` at :959, `logger.warning(` :960-966); inside the staticmethod `_empty_stretches` `:892`, called once from `from_minimal` `:748` | yes |
| base-branch chain order | `PartitionFamily.py:957` (`if branch_key != base_branch: on_branch.sort(...)`) | — |
| `describe()` | inherited from `PartitionFamily.describe` `:413`; `IteratedHomotopyPartition` (`:633`) does NOT override it | — |
| `check_refined_children_inherit_word` | `tests/helpers/laws.py:2977` | yes |
| `check_member_matching_consistent` | `:2995` | — |
| `check_matrix_is_token_counts` | `:3071` | — |
| `check_k10_is_the_full_shift` | `:3092`; its docstring ALREADY reads "refined matrix its in-splitting `[[1, 1], [1, 1]]` (A.8d: two occurrences, two children)" | step 6 is already done |
| `LAYERS["symbolic"]` | `:3218-3236` (17 checks) | yes |
| `LAW_CASES_ONLY` | `:3242` = `("reference_pairs_clear_of_X", "k10_is_the_full_shift")` | — |
| `_NO_REFINEMENT` set-asides | `tests/cases.py:1375-1379`: `refined_children_inherit_word` and `member_matching_consistent` on k28_one_blast, p3, inversion, inversion_period2 | — |

Baseline: `tests/invariants/test_law_symbolic.py`, `tests/unit/topology/test_symbolic_dynamics.py` and `tests/unit/topology/test_partition_family.py` give **55 passed, 28 skipped** at `HEAD` (verified, 9.6 s).

## 2. Theorem T at HEAD

**How it was measured (verified).** The script is `scratchpad/b7/probe.py`. For every class of `session.symbolic_dynamics().classes`, for every non-loop member whose `_image_chain(topology().trellis, member)` is non-empty, it computes:

```
trellis_itinerary(trellis, session.iterated_partition(), chain, member.direction)
translate_itinerary(it, dyn.naming, dyn.table, letters, fixed_points)
```

The token sequence it compares is `(bridge_class, direction, occurrence)`. Loops are already dropped by `translate_itinerary`. Each sequence is then freely reduced (a token next to its inverse with the same occurrence cancels).

- A class is "tested" when two or more members have chains.
- As a cross-check, the per-member tokens equal `cd.evidence[*].symbols` on every build and step: `evidence_mismatch = 0` everywhere.
- **Free reduction never fired** (`reduction_needed = 0` everywhere).

The series runs replicate `cases.blast_series` (the same `_blast` call per step, `session.symbolic_dynamics()` after each blast). Classes are shown as `letter:kind<members>`.

### 2a. Law cases and deep runs

| Build | Classes tested (members with chains) | Disagreements |
|---|---|---|
| k10 | 1: `a` active (2) | 0 |
| inversion_two_blasts | 2: `a` active (3), `u` inert (2) | 0 |
| heteroclinic | 7, all on the period-3 tangle: `bc`(4) `w` inert(4) `bi`(4) `bk`(4) `bm`(2) `bn`(2) `bo`(2) | 0 |
| k28_one_blast, k28_two_blasts, p3, nested, nested_inner_blasts, inversion, inversion_one_blast, inversion_period2 | 0 (not applicable) | — |
| k28 x8 (`build_k28(blasts=8)`) | 4: `a`(23) `u`(17) `b`(12) `c`(7) | **1** (stable end, below) |
| p3_16 | 4, period 3: `u` inert(3) `aa`(2) `ab`(2) `ac`(2) | 0 |
| p3_6_blasts | 0 | — |
| nested inner 6 (`build_nested(inner_blasts=6)`) | 0 | — |

### 2b. Series (47 steps, 141 class-tests in total)

| Series | Classes tested per step | Disagreements |
|---|---|---|
| k28_series 0..8 | 0,0,0,2,3,4,4,4,4 | steps 6, 7, 8: 1 each (stable end) |
| p3_series 0..9 (**period 3**) | 0,0,0,0,0,0,0,2,3,3 | 0 |
| inversion_series 0..4 | 0,0,2,3,4 | 0 |
| k32_series 0..10 | 0,0,0,2,3,4,9,13,13,13,13 | **step 6: 1, step 7: 3** |
| nested_outer_series 0..11 | 0,0,0,2,3,4,4,4,4,4,4,4 | 0 |

### Disagreements, verbatim

Letters are the session's at that step. Each token is followed by its occurrence pair.

**k28 x8, class `A:p1@0.0/R#0 <-> A:p1@0.0/R#4` (`a`).** The same disagreement appears on `k28_series` steps 6–8, where the letter is `c` and the tokens are otherwise identical.

```
(0, 3):   a(R_(0.0;1)^1,R_(0.0;5)^3) u^-1(L_(0.0;1)^2,L_(0.0;3)) b^-1(R_(0.0;3),R_(0.0;5)^1)
(44, 45): a(R_(0.0;1)^1,R_(0.0;5)^3) A:p1@0.0/L#2 <-> A:p1@0.0/R#2(L_(0.0;3),R_(0.0;3))
```

- The chain of `(44, 45)` is `[(74, 78), (78, 79)]`. Crossings 74 and 78 lie on `s(0,0)` and both have sign -1, so the image passes the end of the computed stable branch.
- `cd.evidence` already marks this member `diagnosis.cause = "stable end"`, `word_agrees = None` (verified).
- This is incomplete evidence, not a T counterexample. With such chains excluded (any two consecutive chain crossings on one branch with one sign, the `_diagnose` criterion at `SymbolicDynamics.py:1262-1270`), the four disagreements vanish: one member is skipped per build or step and the tested class counts are unchanged (verified, `SKIP_STABLE_END=1` run).

**k32_series step 6, class `A:p1@0.0/R#7 <-> A:p1@0.0/R#9` (letter `k`, active, unsplit).**

```
(6, 7):   n(R_(0.0;4),R_(0.0;20)) v^-1(L_(0.0;1)^2,L_(0.0;3)^3) m(R_(0.0;17),R_(0.0;19)^1) v(L_(0.0;1)^4,L_(0.0;3)^2) l^-1(R_(0.0;5)^2,R_(0.0;19)^4)
(54, 55): n(R_(0.0;4),R_(0.0;20)) v^-1(L_(0.0;1)^2,L_(0.0;3)^3) m(R_(0.0;17),R_(0.0;19)^1) v(L_(0.0;1)^5,L_(0.0;3)^1) l^-1(R_(0.0;5)^2,R_(0.0;19)^3)
```

- At step 7 the same `(6, 7)` line disagrees with `(54, 55)`, `(64, 65)` and `(96, 97)`, which all read like `(54, 55)` above.
- The walked word equals the `(54, 55)` reading. So the walk agrees with the majority, and `(6, 7)` is the odd one out.
- Word level agrees (`verified = True`). Today this shows only as the INFO "word verified but the registered itinerary differs" (`SymbolicDynamics.py:1352-1358`).
- At step 8 the trellis changes (18 → 45 classes) and the disagreement disappears.

**Why k32 step 6 disagrees.** I verified this with `scratchpad/b7/dig.py` and `jordan.py`.

- Both members are hole bridges with nested lobes. In stable order on the right row, `54 < 7 < 6 < 55`. Their images keep that order: `90 < 13 < 8 < 95`, so this is not an order flip.
- The image chain of `(6, 7)` is `8 → 9 → 10 → 11 → 12 → 13`, all on `s(0,0)`. Its stable order is `13, 8, 11, 10, 9, 12`: the image spirals, and its two left excursions `(11,12) ⊃ (9,10)` nest.
- Chain-order pairing of the interior crossings gives the chords `(9,10)` and `(11,12)`.
- Point-in-polygon of each stable gap's midpoint against `f(b) ∪ σ'` (the image lobe's polygon) gives:
  - `(8,11)` out, `(11,10)` **in**, `(10,9)` **out**, `(9,12)` **in**;
  - for the sibling image `90..95`: `(92,93)` in, `(93,94)` out, `(94,91)` in.
- So the stretches the lobe image really covers are `(11,10)` (already hole `(10,11)`'s stretch) and `(9,12)`. These are the pairs in **stable** order, which is the C19 rule already used on foreign branches.
- Chain order instead marks `(10,9)` empty, and the sibling's image legitimately crosses there (at 93 and 94, inside a stretch outside `f(L)`).
- With the chain-order chord, crossing 10 is pinched into a singleton `L_1^4` and 9 is put in `L_3^2`. The sibling's `(93, 94)` reads `(L_1^5, L_3^1)`. Under stable-order pairing, both members' pairs land in `([10,93], [94,9])`. This last claim is a hand computation, *believed*; the patch run below agrees.

**Patch experiment (verified, scratch only).** I monkeypatched `IteratedHomotopyPartition._empty_stretches` so it sorts by stable cdist on the base branch too (one line: drop the `branch_key != base_branch` test). I reran all 62 builds and steps. Results:

- T: 0 disagreements on k32 steps 6 and 7. The k28 stable-end case still disagrees until it is excluded.
- `M`, `M_1`, node names, class counts, active counts, `is_reliable`: identical on all 62.
- In-splitting: 0 violations.
- What does change is the number of split classes, all of them **inert**, which sit outside both matrices. It drops by 1 on 15 builds and steps: heteroclinic, p3_16, p3_6_blasts, nested inner 6, p3_series 6–9 and k32_series 2–5 and 8–10.
- Iterated element names change wherever chain and stable order differ. That can only happen when a base branch carries 4 or more interior crossings, so the k=2.8 fixture's single chord `(9,10)` is unaffected.

**Period-k coverage.** T is tested on the period-3 tangle in heteroclinic (7 classes), p3_16 (4) and p3_series blasts 7–9 (2, 3, 3), and on the inversion saddle (two branches) in inversion_two_blasts and inversion_series 2–4. All agree.

Note that the plan's "(p3 or k32 at depth)" treats k32 as period-k. It is **period 1** (the probe reads `periods = [1]` on all its classes).

## 3. In-splitting identities (verified on all 62)

The script reads `transition_matrix(refined=True)`, `transition_matrix(refined=False)` and `symbol_nodes(refined=True)`; a node's `parent` is the class letter. It checks three things:

1. sibling rows of `M_1` are equal;
2. for every row `r` of `M_1` and every parent column `d`, `Σ_k M_1[r][d_k] = M[parent(r)][d]`;
3. the node set of `M_1` equals {children of the split active classes} ∪ {letters of the unsplit active ones}, and the set of parents equals the node set of `M`.

**0 violations on every build and step.**

The unit test the plan asks for has teeth: on the existing `_Layout` (`test_symbolic_dynamics.py:65`), setting `refined_rules["a_2"]` to its first token changes `M_1` from `[[1,1],[1,1]]` to `[[1,1],[1,0]]` (verified). `transition_matrix` reads `refined_rules` live, so the sibling-row check fails as required.

There is one trap for the law author. Because `_refined_symbol` (`:667-669`) returns an unmatched token unchanged and `transition_matrix` drops a target that is not a node (`:1012`), a token whose occurrence matched no child would silently lose a matrix entry. The column-sum identity catches exactly that; the sibling-row check alone would not.

## 4. Odd-count case (verified)

- Today it is `logger.warning` at `PartitionFamily.py:960-966`, counted over the INTERIOR crossings `crossings[1:-1]` per stable branch. On the base branch an odd interior count is the same as an odd total.
- **It fired on none of the 62 builds and steps.** I captured the logger, and the capture is proven by the INFO "iterated homotopy partition: …" being caught.
- Structural facts for step 2:
  - `_empty_stretches` is a `@staticmethod` that returns only the stretch dict, so `incomplete_footprints` has to be returned alongside it or collected in `from_minimal` (`:748`) and passed to `__init__` (`:645`), which has no such parameter today.
  - `describe()` is the base class's (`:413`), so printing the footprints needs an override or a hook.
  - The WARNING sits inside a per-hole-bridge loop. "One INFO per build naming them" means collecting first and logging once after the loop.

## 5. What differs from the plan section (B7, A.8d) and the request

1. **T is not a law on every build until the author decides on k32_series steps 6 and 7** (Rubric 5, Q15). The plan's probe (2026-10-09, before B1, B6a and B5) ran on no k32 step; its series was never probed for T.
2. **Stable-end chains must be excluded from T**, or the law fails on k28 x8 and `k28_series` 6–8. A.8d itself says an odd chain has "left the computed stable manifold". The plan's law text does not exclude such chains; the C11 criterion does.
3. **Applicable law cases** are k10 (1), inversion_two_blasts (2) and heteroclinic (**7**, plan: 2), which is 10 classes, not 5. The other 8 law cases need a `NOT_APPLICABLE` entry. k28 x8 has 4 classes (plan: 3). p3_16 (4) was not in the plan's list.
4. **File list.** The set-asides live in `tests/cases.py` (`NOT_APPLICABLE`, `:1356`), and the series records are built by `cases.blast_series` (`:1151`) into `SeriesStep` (`:1085`). That function keeps no session, so a per-step T count in `test_class_count_under_blasting.py` needs a new `SeriesStep` field (or its own blast loop). `tests/cases.py` is therefore a **seventh file** the request does not list. Renaming `check_refined_children_inherit_word` would also touch `tests/cases.py:1377` and the docstring of `tests/invariants/test_law_symbolic.py:12`.
5. **Step 6 is already done**: the docstring at `laws.py:3093-3096` already says "in-splitting" and cites A.8d.
6. **B5 did not make `P1` finer.** This is confirmed: the class counts, `M` and `M_1` I measured match the B5 analyst's reading. A.8d's "finer on 6 of 11 cases" is moot, as the request says.
7. The PartitionFamily Dev Notes (`:94-106`) claim that the base branch keeps chain order and that this "changes refinement labels only, never walks or words". The labels claim is now contradicted at the occurrence level on k32 step 6, where member `(6, 7)`'s refined occurrence differs from the walk's.

## Recommendation (the smallest change)

- Steps 1, 2, 3, 5 and 6 can be planned as written. They change no number; the in-splitting identities hold everywhere and the odd-count case never fires.
- For step 4, the planner should:
  - (a) write the T law to skip members whose image chain has a "stable end" run (the existing C11 criterion). This is needed in any case.
  - (b) put the k32 question to the author. Q15 says stop, and the author must pick one of:
    - **(i)** change the base-branch chord pairing to stable order. The evidence here supports this: the geometry on k32 step 6, and the C19 Jordan argument, which reads a branch from its anchor. It changes only inert splits and iterated element names on the corpus. It is a `from_minimal` cut change, which B7's not-doing list and the author's 2026-09-23 chain-order rule both exclude, so it is either a separate slice before B7's series law or a widening of B7.
    - **(ii)** keep chain order, register `(members_share_the_P1_itinerary, k32_series)` as a `KNOWN_ISSUE`, and name the work package that fixes it.

  I lean towards (i). I have not proved that the Jordan reading holds on the base branch in general; it is verified on this one lobe.

## Out of scope (noticed, not for B7)

- `refine` counts occurrences from inert classes' words too (`:625`), so an active class can be split by an occurrence that only an inert word names. The in-splitting identities still hold, but such a child has no incoming edge in `M_1`. I measured no instance and did not search for one.
- `inversion_series` and `nested_outer_series` show T tested on 2–4 classes per step from blast 2–3 on. These are the only period-1 evidence beyond k28 and k32.
- `test_foreign_chords.py` checks the stable-order pairing geometrically on foreign branches only. A base-branch analogue would have caught the k32 lobe.

Scratch scripts are in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b7/`: `probe.py` (T, in-splitting and odd-count, with the `SKIP_STABLE_END` env var), `probe_patched.py` with `patched_empty.py` (the stable-order monkeypatch), `dig.py` and `jordan.py` (the k32 point-in-polygon). Raw results are in `out/`, `outp/` and `outs/` (`*.jsonl`).

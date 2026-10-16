# Findings — B5 `inertness-is-transience` (analyst, HEAD `7d3a294`, 2026-10-14)

*(The analyst's Write tool refused the artifact; the orchestrator saved this text verbatim.)*

**Dirty tree.** `git status` at `HEAD` `7d3a294` (main) shows uncommitted documentation only:
- modified: `.refactor/decisions.md`, `.refactor/knowledge/{loom,numerics,testing,topology}.md`, `CLAUDE.md`, `docs/alignment_plan_2026-10-09.md`, `docs/code_fidelity_audit_2026-10-09.md`;
- a staged rename of the previous `artifacts/request_origin_analyst.md` into `.refactor/runs/2026-10-12-oracle-laws/`;
- untracked: `artifacts/`, `.refactor/runs/2026-10-1{2,3,4}-*`, `figures/henon_period3_*`, `scripts/henon_period3_*_sweep.py`.

No `src/` or `tests/` file is dirty. The branch starts from `HEAD`, and none of these files may go into it.

Baseline suite at `HEAD` (*verified*, ran it): **1193 passed, 38 skipped, 5 deselected, 1 xfailed** (222 s).

**Kind of request.** This is a refactor that moves one rule. Inertness moves from the evidence table (`BridgeClass.py`) to the word graph (`SymbolicDynamics.py`), and the lettering moves with it. The minimal trellis stops filtering by activity. I read the request as: verify the code spans, re-run the plan's probe at HEAD, and report what differs from the plan section. I do not re-plan the slice.

## Answer in brief

- **The plan holds at HEAD.** Class-level, letter-free words are **identical on all 15 builds**: the 11 law cases, `p3_6_blasts`, `p3_16`, k28 x8 and nested with 6 inner blasts. `is_reliable` is unchanged (True on all 15). No class goes from inert to active.
- **5 classes change verdict**, all from active to inert, and all were active only for lack of evidence: `inversion_two_blasts` g; `heteroclinic` d and bj; `p3_6_blasts` ab; nested inner 6 ae.
- **No unresolved class exists anywhere on the corpus**, in either run. So no word names an unresolved class, and "unresolved is active" is never exercised.
- **The plan's step 3 and its step-4 law disagree on exactly that case** (§3). This is the one correction the planner must carry.
- **The plan's seven files are not enough as written.** `session.bridge_classes()` letters before the dynamics runs. At least `tests/facade/test_session_bridge_classes.py` and `tests/helpers/direct.py` read `entry.active` / `entry.letter` on that table before any dynamics exists. The smallest fix (§5) stays inside `Topology.py`.

## 1. The inertness code at `HEAD` (*verified*, read at `7d3a294`)

| What | Where |
|---|---|
| Dev Notes "INERT = MAPS TO NO ACTIVE CLASS" | `src/tanglepack/topology/BridgeClass.py:51-71` |
| Slots `inert`, `image_classes`, `image_loops`, `unresolved` (docstring :201-213) | `BridgeClass.py:230-233` |
| `active` property | `BridgeClass.py:239-242` |
| `has_image_evidence` | `BridgeClass.py:250-253` |
| `describe()` prints the evidence | `BridgeClass.py:312-353` |
| Call sites in `bridge_classes` | `BridgeClass.py:712-714` |
| `_collect_image_evidence` | `BridgeClass.py:783-823` |
| `_resolve_inertness` | `BridgeClass.py:826-849` |
| `_image_chain` (stays: loop folding, minimal trellis, trellis itinerary) | `BridgeClass.py:725` |
| `minimal_trellis` mapping loop | `MinimalTrellis.py:310-334` |
| — its `entry.active` test | `MinimalTrellis.py:326` |
| `inert_hole_bridge_ids` (field, summary, describe, constructor) | `MinimalTrellis.py:91/117, 216, 229-233, 425` |
| `Topology.classes` | `loom/Topology.py:143-152` |
| — its lettering loop | `loom/Topology.py:148-150` (`zone_key` at :151) |
| `Topology.dynamics` passes `letter_for` | `loom/Topology.py:182-185` |
| `BridgeAlphabet.letter_for` (unchanged by B5) | `loom/BridgeAlphabet.py:55` |
| `inert_letters` (reads `table.inert`) | `SymbolicDynamics.py:466-478` |
| `_active_letters` | `SymbolicDynamics.py:481-483` |
| `translate_itinerary` sets kind from `entry.inert` | `SymbolicDynamics.py:551` |
| `symbol_nodes` (active only) | `SymbolicDynamics.py:871-899` |
| `_edges` | `SymbolicDynamics.py:904-933` |
| — PK-8's WARNING | `SymbolicDynamics.py:919-925` |
| Driver: initial letters | `SymbolicDynamics.py:1383-1386` |
| Driver: closure `add` | `SymbolicDynamics.py:1389-1415` |
| Driver: kind | `SymbolicDynamics.py:1426` |
| Driver: class-level tokens `named` | `SymbolicDynamics.py:1446-1449` |
| Driver: translate loop | `SymbolicDynamics.py:1462-1466` |
| Law `check_minimal_trellis_bridges` | `tests/helpers/laws.py:2295-2336` |
| Law `check_inert_classes_outside_transitions` | `tests/helpers/laws.py:3006-3015` |

The plan's citation `Topology.py:148` is still correct.

**PK-8's WARNING text differs from the plan's quote.** The code says "inert class %s has the word %r, which names a non-inert class". It fires only when an inert word names a non-inert symbol.
- It fired on **no build** in the probe (*verified*: WARNINGs were captured during both runs).
- Under I-a it cannot fire at all, so deleting it removes dead code.

**Ordering fact the planner needs** (*verified* by reading).
- The driver's `named` list at :1446-1449 already holds every class-level token.
- It is computed by the same `oriented_class(parent_of(x), parent_of(y))` that `translate_itinerary` uses at :537-538. A `direction == 0` token is a loop.
- So the transience fixed point can run between the closure loop and the translate loop, on `named`, with no new walking.
- The letters and the `kind` that `translate_itinerary` needs come after it.

## 2. The probe re-run at `HEAD`

**Method** (*verified*). The script is `probe.py`, with `analyze.py` and `loops.py` beside it, all in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b5/`.
- **Prototype:** it monkeypatches `tanglepack.loom.Topology.bridge_classes` to mark every evidence entry active. As a result `minimal_trellis` maps every hole bridge, and the dynamics walks every class.
- **Runs:** each case is built once with the `tests/cases.py` builders. Then `session.topology(rebuild=True)` runs once at base and once under the patch.
- **I-a verdict:** computed from the patched run's class-level words. A class is inert iff its word is resolved and every token in it is inert; an empty word counts as inert. Unresolved classes count as active.
- **Matching:** classes are matched by element identity (fixed point, branch, side, `lo_id`, `hi_id`), never by letter.

| Case | Words equal | Image bridges | Iterated elements | Inert (base → I-a) | Verdict changes | `is_reliable` |
|---|---|---|---|---|---|---|
| k10 | yes (2 classes) | 3 → 3 | 11 → 11 | 1 → 1 | — | T → T |
| k28_one_blast | yes (3) | 0 → **1** | 8 → **11** | 2 → 2 | — | T → T |
| k28_two_blasts | yes (4) | 3 → 3 | 13 → 13 | 1 → 1 | — | T → T |
| p3 | yes (5) | 0 → **2** | 16 → **21** | 2 → 2 | — | T → T |
| nested | yes (9) | 3 → **5** | 29 → **34** | 3 → 3 | — | T → T |
| nested_inner_blasts | yes (9) | 3 → **5** | 29 → **34** | 3 → 3 | — | T → T |
| inversion | yes (4) | 0 → **2** | 14 → **19** | 2 → 2 | — | T → T |
| inversion_one_blast | yes (7) | 3 → 3 | 21 → 21 | 1 → 1 | — | T → T |
| inversion_two_blasts | yes (8) | 3 → 3 | 22 → 22 | 1 → **2** | **g** → inert | T → T |
| inversion_period2 | yes (6) | 0 → **2** | 24 → **29** | 2 → 2 | — | T → T |
| heteroclinic | yes (75) | 8 → **10** | 158 → **164** | 3 → **5** | **d, bj** → inert | T → T |
| p3_6_blasts | yes (43) | 4 → 4 | 92 → 92 | 1 → **2** | **ab** → inert | T → T |
| p3_16 | yes (43) | 4 → **6** | 92 → **98** | 2 → 2 | — | T → T |
| k28 x8 | yes (4) | 3 → 3 | 13 → 13 | 1 → 1 | — | T → T |
| nested, 6 inner | yes (47) | 7 → 7 | 105 → 105 | 2 → **3** | **ae** → inert | T → T |

Letters are HEAD's. Every changed class was active at HEAD with `has_image_evidence` False.

**Further results** (*verified*):
- **Closed alphabet:** the closed-alphabet class set (by identity) is identical base vs I-a on every build. The closure adds and drops nothing.
- **Class matrix:** equal on the surviving nodes everywhere. The five moved classes simply leave the node set:
  - inversion_two_blasts 7 → 6
  - heteroclinic 72 → 70
  - p3_6_blasts 42 → 41
  - nested inner 6 45 → 44

  They are transient, so dropping them removes only rows that no cycle passes through.
- **Refined active nodes:** new in-splits occur only on **heteroclinic `bk` (1 → 2 children)** and **p3_16 `ab` (1 → 2)**, besides the dropped nodes of the five moved classes. On k28_one_blast, p3, nested, inversion x0 and inversion_period2, the finer `P1` changes element names and itineraries but adds **no** active refined node.
- **Verdicts from HEAD's own words:** transience computed from HEAD's own words gives the same inert set as the I-a run on all 15 builds. The verdict does not depend on the extra image bridges.
- **Loops:** every alphabet class holding a loop member is inert under I-a on all 15 builds (`loops.py`). The facade test's "a class holding a loop is inert" therefore still holds empirically. *Believed*, not proven in general.
- **New `unmapped` entries:** hole bridges with no registered image now land in `MinimalTrellis.unmapped`:
  - +1 each on k10, k28_two, nested, nested_inner, inversion_one and p3_6;
  - +2 on nested inner 6.

  The DEBUG text at `MinimalTrellis.py:331-333` still says "active hole bridge".
- **New WARNING:** one new "image pair … of hole bridge … has no Bridge object … skipped" (D9) appears on k28_one_blast, inversion_two_blasts and k28 x8. Two of these are law cases. No test asserts a quiet log there (*verified*: `skipped_pairs` is read only in `laws.py` and `test_minimal_trellis.py:40`).

**Against the plan's probe (2026-10-09, pre-B1):**
- "Words unchanged on 14 builds" holds at HEAD on 15.
- "Image bridges +0..2" holds.
- **Builds where `P1` gets finer:** the plan names 5 (k28_one_blast, p3, nested, inversion x0, heteroclinic) and says "`P1` finer on 6". HEAD gives **8**: the 5 named plus **nested_inner_blasts, inversion_period2, p3_16**.
- "Iterated-element deltas +0..6" holds (+3, +5, +6).
- **Refined in-splitting:** the plan's "the refined matrix in-splits finer there" holds only on heteroclinic and p3_16. Elsewhere only the names change.
- **Inert-verdict changes:** the plan names none; HEAD has the five above.
- **k28 x8:** 4 classes, 1 inert, consistent with "4/4/4 since B1". Nothing changes there.

## 3. Unresolved classes and "unresolved is active"

**On the corpus** (*verified*: every class has a word, both runs). None of B6a's three unresolved reasons occurs on any of the 15 builds:
- a truncated search;
- a landing span meeting no stable node;
- member lobes holding different markers.

So no class's word names an unresolved class, and the clause changes nothing on the corpus.

**How the clause propagates under step 3** (*verified* by construction of the fixed point).
- An unresolved class never enters the inert set.
- A class with an unresolved token can never have "every token inert".
- So **a class is inert iff no cycle AND no unresolved class is reachable from it**.

**Where the plan's step 4 disagrees.** Its law reads "`inert` equals no cycle reachable in the transition graph over ALL classes". An unresolved class has no rule, so in that graph it is a sink. The law would call it, and every class reaching only it, **inert**, which contradicts step 3.
- On the corpus the two readings agree (0 disagreements, *verified*), but only because nothing is unresolved.
- The law must read "neither a cycle nor an unresolved class is reachable".
- It cannot use `transition_graph()`, which drops inert classes and their tokens. It must build its own all-class graph from the class-level words (`cd.symbols` → `symbol.bridge_class`).
- A unit case should pin this, since the corpus never reaches it: `a -> b` with `b` unresolved, both active. The plan's `a -> b, b -> (empty), c -> c` test does not cover it.

## 4. Consumers outside the seven files (*verified* by grep at HEAD)

`session.bridge_classes()` returns `Topology.classes` (`TangleSession.py:339`). Today that table letters active classes at :148-150 without running the dynamics. Under the plan, `entry.inert` and `entry.letter` stay None until the dynamics runs, and `active = not inert` is True when `inert` is None.

**Breaks or changes meaning:**
- **`tests/facade/test_session_bridge_classes.py:73-80`:** asserts `(entry.letter is not None) == entry.active` and "a class holding a loop is inert" on `session.bridge_classes()`, and never calls the dynamics. **It fails as planned.**
- **`tests/helpers/direct.py:107-125`:** letters `entry.active` entries from a fresh `BridgeAlphabet` before calling `symbolic_dynamics` with no `letter_for`.
  - Its views are letter-free, so the equivalence tests probably still pass (*believed*).
  - But it would pre-letter classes the dynamics then calls inert, and the loop mirrors the one B5 deletes.
- **`tests/unit/topology/test_minimal_trellis.py:94`:** asserts every `unmapped` hole bridge is `active`. Under I-a, `unmapped` gains inert-class hole bridges, so this breaks wherever its fixture has one.

**In scope, affected:**
- **`laws.py:2116-2125` `check_only_active_classes_lettered`:** reads `session.bridge_classes()`. With session-scoped law cases, its result depends on whether an earlier law already ran the dynamics.
- **`laws.py:2024` `check_anchor_bridge_class_leads_its_tangle`:** reads `entry.active` on the EVIDENCE table. For entries outside the alphabet, `inert` stays None forever, so `active` is vacuously True.

**Unaffected:**
- `tests/cases.py:1192-1193`: the SeriesStep counts are read after `session.symbolic_dynamics()`, on the same entry objects.
- `tests/regression/test_blast_children_are_evidence_only.py:53`: reads `dynamics.table`.

**Not tests:**
- `scripts/henon_bridge_classes.py:275-318` and `stress/numerics/e11_blast_depth.py:359-393` read `entry.inert` from `session.bridge_classes()`.
- CLAUDE.md's description of inertness (`_resolve_inertness`, `image_classes`, "Inert = at least one resolvable image…") goes stale. CLAUDE.md is dirty in the tree and is not the branch's business.

**Letter shift** (the letter-free rubric is unaffected). Once a class moves to the inert series, every active letter assigned after it in a FRESH session shifts by one, because `BridgeAlphabet.letter_for` uses `letter(len(assigned))`. This affects inversion_two_blasts, heteroclinic, p3_6_blasts and nested inner 6.

## 5. Recommendation (smallest change)

Keep the plan's five steps, with three adjustments that do not re-plan the slice.

1. **Make the step-4 law match step 3.** `inert_is_transient` should read: inert iff neither a cycle nor an unresolved class is reachable in the graph over ALL classes' class-level words. Add a unit case `a -> b, b unresolved` (both active) next to the plan's three-class table.
2. **Keep `session.bridge_classes()` annotated, inside `Topology.py`.**
   - The dynamics reads the zoned, unlettered alphabet (`self.minimal.alphabet` plus `zone_key`, as a private cached stage).
   - `Topology.classes` returns that same table after touching `self.dynamics`, so `letter` and `inert` are set on the shared entry objects.
   - This keeps `test_session_bridge_classes.py` and `check_only_active_classes_lettered` meaningful without touching them.
   - Cost: `bridge_classes()` now computes the dynamics. Nothing at HEAD calls it where the dynamics fails.
   - Alternative: accept `test_session_bridge_classes.py` and `direct.py` as extra files.
3. **Also edit two test files:** `tests/unit/topology/test_minimal_trellis.py:94` (drop the `active` assertion on `unmapped`) and `tests/helpers/direct.py`'s lettering loop (dead under B5). That makes 9 files instead of 7, all tests-only. If the author holds rubric 1 to seven files, the planner should ask rather than leave `direct.py` pre-lettering inert classes.

**Implementation hint** (not a plan).
- Run the transience fixed point on the driver's `named` tokens after the closure loop (`SymbolicDynamics.py:1449`) and before letters and translation.
- The closure's `add` (:1389-1415) then letters nothing. Lettering happens once, in table order, after the verdicts.
- This keeps today's letter order for the classes that stay active.

## Out of scope (noticed, not recommended here)

- `BridgeClassEntry.inert = None` makes `active` (`not inert`) True for every evidence entry outside the alphabet. No reader needs it once the laws read the dynamics.
- The new D9 "no Bridge object" WARNING on three builds, two of them law cases, is true: the image of an inert hole bridge reaches an arc no blast computed. It may read as log noise.
- `MinimalTrellis.py:331`: the DEBUG text still says "active hole bridge".
- Plan wording: "P1 finer on 6 builds" against the five named; HEAD gives 8 (§2).

## Verified / believed

- *Verified* (ran):
  - the baseline suite at HEAD;
  - the probe on all 15 builds: word equality, verdicts, the image-bridge / iterated-element / unmapped / skipped deltas, class matrices on surviving nodes, refined active nodes, `is_reliable`, captured warnings, and loops ⇒ inert;
  - every line span in §1;
  - every consumer in §4, by grep and reading.
- *Believed*:
  - the `direct.py` equivalence tests still pass under B5 (its views are letter-free);
  - "holds a loop ⇒ inert" holds in general;
  - the monkeypatch prototype equals the planned implementation on words. It walks the same classes with the same walls, and inertness feeds only `kind` and letters, which the comparison ignores.

**Knowledge file updated:** `.refactor/knowledge/topology.md`, new closing section "B5 `inertness-is-transience` analyst" (code map with line spans, probe numbers, out-of-scope consumers), and its "Last updated" line.

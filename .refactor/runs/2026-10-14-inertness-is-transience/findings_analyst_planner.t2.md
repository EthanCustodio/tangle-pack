# Findings, round 2: B5 `inertness-is-transience` (analyst, HEAD `7d3a294`, 2026-10-14)

*(The analyst's Write tool refused the artifact; the orchestrator saved this text verbatim.)*

**The working tree is dirty, but only in documentation and run files.** `git status` at `7d3a294` (main) shows the same state as round 1:
- modified: `.refactor/decisions.md`, `.refactor/knowledge/{loom,numerics,testing,topology}.md`, `CLAUDE.md` and two `docs/` files;
- a staged rename of the old `artifacts/request_origin_analyst.md`;
- untracked: `artifacts/`, `.refactor/runs/2026-10-1{2,3,4}-*`, `figures/henon_period3_*`, `scripts/henon_period3_*_sweep.py`.

No file under `src/` or `tests/` is dirty. The branch starts from `HEAD` and none of these files go into it.

## The author's note

> "I believe you can determine if a bridge is inert before constructing the dual graph. Every bridge that has an initial hole punched in it should be inert. The zeroth holes."

**Answer: yes, and it holds on every build.** Every bridge carrying a zeroth hole lies in an inert class, on all 15 builds. Its class word is not just inert but **empty**: the bridge's image lands in one element.

This is a theorem, given two conditions stated in §3, and not a corpus accident. It needs only `trellis.holes` and the class table, so it can be decided before the minimal trellis and the dual graph.

**It does not decide every class.** One hole-free class is inert: heteroclinic `u`, a piece of the anchor bridge's image whose own image is a loop. So the author's rule is a sufficient condition, not the definition.

**What I recommend:**
- The author's rule decides what the minimal trellis maps. Every hole bridge is mapped except the zeroth-hole bridges, as before the plan.
- Transience over the words stays the verdict for every class.
- A new law checks the theorem: every zeroth-hole class has an empty word.

With this rule the minimal trellis, `P1`, the image bridges, the refined classes and the words are **identical to HEAD on all 15 builds**. The only change is the verdicts of the five classes round 1 found. The plan's "P1 finer on 8 builds", its new `unmapped` entries and its three new D9 WARNINGs all disappear. Details in §5.

**Kind of request.** This is a revision of a refactor plan: one rule's home and its pre-dual-graph use. The targets are:
- `minimal_trellis` (`src/tanglepack/topology/MinimalTrellis.py:310-334`, the `entry.active` test at `:326`);
- `bridge_classes` / `_resolve_inertness` (`src/tanglepack/topology/BridgeClass.py:712-714, 826-849`);
- the planned transience step in `symbolic_dynamics`.

## 1. Method (*verified*, run at `HEAD`)

The scripts are in `/tmp/claude-1000/-home-dezhu-Development-TanglePack-tangle-pack/c679944b-fa37-474d-8e9a-1a75fda47087/scratchpad/b5r2/`. Each case is built with the `tests/cases.py` builders: the 11 law cases, `p3_6_blasts`, `p3_16`, `build_k28(blasts=8)`, and `build_nested(inner_blasts=6)`. Classes are matched by element identity (fixed point, branch, side, `lo_id`, `hi_id`), never by letter. Letters below are HEAD's.

- **`holes.py`** classifies every hole in `topo.holes`.
  - A hole is a *zeroth* (reference) hole when `Hole.pair is not None`, which is the same as `iterate == 0`. Holes with `iterate < 0` are propagated.
  - Its bridge is `trellis.bridge_between(*hole.bounding_ids)`, and its class is `topo.minimal.classes.entry_of(bridge.id)` in the evidence table.
  - The script builds five sets:
    - **A**: classes holding a zeroth-hole bridge.
    - **B**: classes holding any hole bridge.
    - **C**: HEAD's inert set.
    - **D**: the I-a inert set from round 1's prototype.
    - **E**: classes holding a loop member.
- **`probe2.py`** prototypes the author's rule. It monkeypatches `loom.Topology.bridge_classes` in two ways:
  - **R1**: `inert` = holds a zeroth-hole bridge, alone.
  - **R2**: zeroth-hole OR HEAD's evidence verdict.

  It runs `topology(rebuild=True)` and compares the result against HEAD: words, image bridges, iterated elements, `unmapped`/`skipped`, refined classes, verdicts against D, `is_reliable`, and new WARNINGs.
- **`pre.py`** checks two things: whether a zeroth-hole bridge also carries a propagated hole, and the `-1` iterates of zeroth-hole bridge ends.

## 2. The comparison, per build (*verified*)

| Build | zeroth / propagated holes | A (zeroth-hole classes) | C (HEAD inert) | D (I-a inert) | E (loop classes) | A vs D |
|---|---|---|---|---|---|---|
| k10 | 1 / 2 | u | u | u | u | A = D |
| k28_one_blast | 2 / 2 | u, v | u, v | u, v | v | A = D |
| k28_two_blasts | 1 / 3 | u | u | u | — | A = D |
| p3 | 2 / 6 | u, v | u, v | u, v | u, v | A = D |
| nested | 3 / 9 | u, v, w | u, v, w | u, v, w | v, w | A = D |
| nested_inner_blasts | 3 / 9 | u, v, w | u, v, w | u, v, w | v, w | A = D |
| inversion | 2 / 4 | u, v | u, v | u, v | u, v | A = D |
| inversion_one_blast | 1 / 6 | u | u | u | u | A = D |
| inversion_two_blasts | 2 / 6 | **g**, u | u | **g**, u | u | A = D |
| inversion_period2 | 2 / 8 | u, v | u, v | u, v | u, v | A = D |
| heteroclinic | 4 / 69 | **d**, **bj**, v, w | u, v, w | **u**, **d**, **bj**, v, w | v, w | **A ⊊ D**: D \ A = {u} |
| p3_6_blasts | 2 / 40 | u, **ab** | u | u, **ab** | u | A = D |
| p3_16 | 2 / 40 | u, v | u, v | u, v | u, v | A = D |
| k28 x8 | 1 / 3 | u | u | u | — | A = D |
| nested, 6 inner | 3 / 43 | u, v, **ae** | u, v | u, v, **ae** | v | A = D |

What the table shows:

- **A ⊆ D on all 15 builds.** A = D on 14 of them. The exception is heteroclinic, where D \ A = {`u`}.
- **The five classes that change verdict are exactly the zeroth-hole classes HEAD calls active:** `g`, `d`, `bj`, `ab`, `ae`. In every case HEAD calls them active only because it has no registered image of them.
- **C ∩ B ⊆ A everywhere.** Every hole class HEAD already calls inert holds a zeroth hole. So **D = A ∪ C on all 15 builds.**
- **B ∩ D = A everywhere.** No class holding only propagated holes is inert, on any build.
- **E ⊆ A ⊆ D everywhere.** Every loop-holding class is a zeroth-hole class.
  - The loops are the folded images of the zeroth-hole bridge. Heteroclinic `v`, for example, holds 24 members, 21 of them loops folded from `(77,78)`, `(79,80)`, `(83,84)` and `(89,90)`.
- **The words.** Every class in A has an **empty** word.
  - Its walk is `trivial`: the itinerary is one pair `(e, e)` whose homotopy parents coincide, so the pair is a loop. This was read on heteroclinic `d`, `bj`, `v`, `w`.
  - **No inert class on any build has a non-empty word**, and no active class has an empty one. On the corpus, inert simply means "the word is empty".
  - None of the four alternatives the brief asked about (a word naming only inert classes, a word naming an active class) occurs.
- **No bridge carries both a zeroth and a propagated hole** (`pre.py`, all 15 builds).

**The one difference: heteroclinic `u`.**
- `u` = {`A:0.0/left[114,4]`, `A:0.0/left[None,93]`}. It has one member, the bridge `(4,5)`, which holds no hole.
- `(4,5)` is the middle piece of the anchor bridge's image chain `(0,4) -> [(0,4), (4,5), (5,6)]`. `(0,4)` carries the propagated hole at iterate −3 of the zeroth hole `(113,114)`.
- `(4,5)`'s own image is the loop `(6,110)` in `A:0.0/left` element 0. Its word is empty.
- So `u` is a **hole-free inert class**, inert under HEAD's evidence rule and under I-a alike. Under the author's rule alone (R1) it would be **active**, with an empty word.

## 3. The mathematics

**The zeroth lobe holds exactly one marker.** Let `b` be a zeroth-hole bridge, with ends `x, x'`. Let `σ` be the stable segment between them.
- `σ` lies in the window `S \ f^k(S)`, written `S \ f(S)` below. Because `(x, x')` is a pseudoneighbor pair, both open arcs `b` and `σ` hold no point of X (B2 tests both).
- The lobe `L = b ∪ σ` therefore contains **no arc of the trellis T**:
  - an unstable arc cannot cross `b` (u×u), and could enter `L` only by putting a crossing on `σ`;
  - a stable arc cannot cross `σ` (s×s), and could enter only by putting a crossing on `b`.
- So `L` holds exactly one marker, its own zeroth hole, provided two things hold:
  - (i) no selected periodic orbit lies inside `L` with its whole computed tangle;
  - (ii) no other hole sits beside `b`.
- (ii) follows from clarity: a propagated hole sits on the container of an image sub-arc whose ends are X points. Inside `b` such ends would contradict clarity, and the sub-arc is never `b` itself, since a trajectory has one reference. The corpus confirms it: no bridge carries both kinds of hole.

**Its image holds none.**
- Holes propagate backward only, so the image of the zeroth hole, `(origin, +1)`, is no hole. `DualGraph.image_marker` returns None for it (`DualGraph.py:1172-1173`), and `lobe_image_markers` drops it (`DualWalk.py:826-830`).
- `f` is a homeomorphism, so `f(b) ∪ f(σ)` encloses no marker. `f(b)` is therefore homotopic, relative to the punctures (the markers), to the stable segment `f(σ)`.
- `f(σ)` contains no crossing of T in its interior. If it held one at `y`, then `f⁻¹(y)` would be an X point on the open `σ`.
- Its ends `f(x), f(x')` therefore lie in one homotopy element.
- So the class's image is a loop and its **word is empty**. That makes the class inert under every rule:
  - HEAD's evidence rule, once the image is registered: the loop is inert evidence;
  - I-a: an empty word is transient.

`f(b)` can still fold through the window: its chords cut the ITERATED partition, which is why I-a made `P1` finer. At class level, though, every fold encloses no hole and reads as a loop.

**What this says about the forward images and the backward chain.**
- The zeroth hole is the forward-most hole of its orbit. Its lobe is the last one a hole marks.
- The propagated holes are its backward images. For each propagated hole, the image of its bridge `c_{-n}` must enclose the hole at `-n+1`, so its word is non-empty and names the class of the bridge beside that hole.
- `f⁻¹(b)` is never a bridge of T. The open window's preimage lies outside T, and at most one end of `b` sits on the pip orbit: `pre.py` finds one such end on 5 builds, never two. So `c_{-1}` strictly contains `f⁻¹(b)`, and `f(c_{-1})` runs beyond `b`.
- The corpus finds every propagated-hole class active. This is ***believed***, not proven: the argument shows a non-empty word, not a reachable cycle.

**Is "holds a zeroth hole ⇒ inert" a theorem?**
- **It is a theorem under conditions (i) and (ii) above**, given that D13's walk computes the homotopy class faithfully. That last part is the walk machinery's own claim, and it is *verified* on all 15 builds: every zeroth-hole class's walk is `trivial`.
- **It is true on the whole corpus.**
- **It is false nowhere I looked.** Condition (i) could fail only if a whole inner tangle sat inside a window lobe without crossing it. No build does this.

**Does it decide every class?** No, only the classes that hold a zeroth hole.
- **Hole-free classes can be inert** (heteroclinic `u`). A class can also be inert with no zeroth hole whenever its lobe holds no marker whose image is a marker.
- **Classes holding only propagated holes** are active on the whole corpus, but nothing proves they must be.
- So the converse "inert ⇒ holds a zeroth hole" is **false** on a law case.

## 4. What B5 should compute inertness from

| Option | Verdicts vs D (I-a) | Minimal trellis / `P1` / words vs HEAD | Verdict |
|---|---|---|---|
| Author's rule alone (R1) | equal on 14; heteroclinic `u` becomes **active** with an empty word: a zero row in the transition graph, contradicting both "maps to no active class" (2026-09-16) and I-a | identical on all 15 | insufficient as a definition |
| Transience alone, every hole bridge mapped (the sent-back plan) | = D | `P1` finer on 8 builds, image bridges +1..2, new `unmapped` and 3 new D9 WARNINGs (round 1) | works, but changes more than it needs to |
| **Author's rule before the dual graph + transience as the verdict** (R2 behaviour) | **= D on all 15** | **identical on all 15**: image bridges, iterated elements, refined classes, words, `is_reliable`; `unmapped` −1 on inversion_two_blasts, p3_6_blasts, nested inner 6, −2 on heteroclinic; 0 new WARNINGs | **recommended** |

The R2 numbers are *verified*, from `probe2.py`. The prototype used HEAD's evidence for the hole-free part. Transience gives the same set, D, because D = A ∪ C on every build.

So:
- **The author's rule is the pre-dual-graph decision.** It is the only thing the minimal trellis needs.
- **Transience is the definition for every class.** It agrees with the author's rule wherever that rule applies and extends it to hole-free classes.
- **The theorem becomes a law that checks the rule.**

The author's rule can be the definition only for the minimal trellis's question, "which hole bridges to map", and not for the alphabet's letters.

## 5. Recommendation (the smallest change to the sent-back plan)

Keep steps 2-5 of `artifacts/plan_planner_human_2.md`: delete the evidence machinery, decide the transience verdict in `symbolic_dynamics`, letter after the verdict, and add the law `inert_is_transient` with the unresolved fix from round 1. Change only these:

1. **Step 1, `MinimalTrellis.py`.** Map every hole bridge **except one carrying a zeroth hole**, decided from the `holes` argument `minimal_trellis` already takes. A bridge carries a zeroth hole when some `hole.pair is not None` with `bridge_between(*hole.bounding_ids)` equal to it.
   - Drop the `entry.active` test at `:326`. The class table is no longer read for this decision.
   - Keep the field that records the unmapped zeroth-hole bridges (today `inert_hole_bridge_ids`), renamed or re-documented as "zeroth-hole bridges: their image is a loop (the class word is empty), never mapped".
   - The plan's "maps every hole bridge" and its whole `P1`-finer behaviour paragraph go. `P1` is HEAD's on all 15 builds, which also satisfies Q5 ("P1 finer is a defined consequence") vacuously.
2. **Step 4, `tests/helpers/laws.py`.**
   - Add the law `zeroth_hole_classes_have_empty_words` (layer `symbolic`): every class holding a bridge with a zeroth hole has an empty word, and is therefore inert under `inert_is_transient`. It is the theorem of §3, and it checks the pre-dual-graph rule against the words.
   - `minimal_trellis_bridges` then reads: "every hole bridge without a zeroth hole, with a registered image, contributes its chain".
3. **Files.** With the rule in `minimal_trellis` (bridge-level, as the author worded it), `tests/unit/topology/test_minimal_trellis.py:94` changes from "unmapped ⇒ active" to "unmapped ⇒ no zeroth hole". That line was already in the plan's nine files. The file list does not grow.

Leave the alternative alone: seeding `entry.inert = True` in `bridge_classes` for zeroth-hole classes. Nothing but the minimal trellis would read that seed, and lettering stays after the dynamics in either case, because heteroclinic `u` is decided only there (round 1, §5 point 2, still holds).

**Behaviour under this recommendation** (*verified* by R2 on 15 builds):
- Words: unchanged.
- Minimal trellis, `P1` and refined classes: unchanged.
- Verdicts: the same five classes go from active to inert, and no class goes the other way.
- `is_reliable`: True everywhere.
- `unmapped` shrinks where the five classes had unmapped hole bridges.

## Out of scope

- A sharper pre-dual-graph rule, "a lobe holding no marker whose image is a marker gives an empty word", would also catch heteroclinic `u`. It needs `DualGraph.enclosed`, so it is not pre-dual-graph. Transience covers it.
- Activity of the propagated-hole classes is unproven (§3).
- `CLAUDE.md`'s inertness text ("Inert = at least one resolvable image…", "Inert hole bridges are kept but not mapped") becomes "zeroth-hole bridges are kept but not mapped; inert = transient". The tree's `CLAUDE.md` is dirty and is not the branch's concern.

## Verified / believed

**Verified (ran it):**
- the sets A to E and their relations on all 15 builds;
- empty words for every class in A, and the `trivial` walks on heteroclinic;
- heteroclinic `u`'s provenance (the chain of `(0,4)`, and the image loop `(6,110)`);
- no bridge carries both kinds of hole;
- R1 and R2 against HEAD and D (words, image bridges, iterated elements, refined classes, `unmapped`/`skipped`, `is_reliable`, WARNINGs);
- the `-1` iterates of zeroth-hole ends;
- `image_marker` and `lobe_image_markers` dropping the zeroth hole's image (read).

**Believed:**
- the proof sketch of §3, beyond its corpus check, including conditions (i) and (ii) in general;
- that propagated-hole classes are always active;
- that the planned implementation equals the monkeypatch prototypes.

**Knowledge file updated:** `.refactor/knowledge/topology.md`. The B5 section now has a round-2 subsection on the zeroth-hole rule, and its "Last updated" line is updated.

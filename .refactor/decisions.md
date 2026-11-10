# Decisions

What the author has decided. Binding on every later run. Newest last.

## 2026-10-02 — blasting nested tangles (slug `hole-side-and-own-blast`, base `2315204`)

- **Said** (plan verdict, `plan_ok_planner_human.json`): approved with no notes; the plan's
  defaults stand.
- **Consequence 1 — direct-hole side is combinatorial.** A directly punched hole's
  `Hole.bridge_side` comes from the crossing signs, in the style of `row_of_end`:
  `first` = lower unstable cdist end, `outward = +1 if second.stable_cdist > first.stable_cdist
  else -1`, `left iff first.crossing_sign * outward > 0`. Not the lobe-polygon area, not a
  nearest-vertex measurement.
- **Consequence 2 — a propagated hole's side is read on its image sub-arc** (the `span` inside
  the containing bridge), never against the whole bridge polyline.
- **Consequence 3 — each nested zone blasts only its own fixed point**
  (`fixed_point=[zone.fixed_point]` in `build_nested`); blast order no longer matters.
- **Consequence 4 — CLAUDE.md** carries one sentence: the period-3 side trip at 15-16 steps and
  6+ blasts is fixed, but those runs are unreliable (an unreachable class, a virtual `new1`;
  open). No other CLAUDE.md change.
- Left open by the plan's "Not doing" list (not decided, not to be assumed): p3 >= 6 blasts
  unreliable; `ManifoldMachine.iterate_manifold -> merge_manifolds` non-contiguous crash;
  `_hole_openings` whole-polyline nearest vertex; pip choice lost on generation bump (`_repin`);
  `Blast.py` tidy-up; constant nested slug; noisy inert-word WARNING; blast speed.

## 2026-10-02 — outcome verdict on `hole-side-and-own-blast` (first outcome)

- **Said** (`debrief_ok_auditor_human.json`): "figures", no further notes.
- **Consequence:** proof figures were made on the branch (commit `f8a6d73`, three scripts under
  `scripts/refactor_proof/hole-side-and-own-blast/`; PNGs untracked under
  `figures/refactor_proof/hole-side-and-own-blast/` in the branch worktree) and the auditor
  re-audited (`outcome_auditor_human_2.md`, second pass). The proof scripts are outside the
  plan's five-file list by instruction; rubric 5 is judged against `3b33baa`.

## 2026-10-07 — open-issues round (catalog `docs/open_issues.md`, base `684d70b`)

- **Said** (planning session answers, 2026-10-07), recorded as D1-D8 in `docs/open_issues.md`:
  - **D1 full recompute.** "We only start doing the topology after all the numerics is done …
    any growth that happens requires entirely recomputing the topology." Planned loop: blast,
    compute the topology, blast again, recompute, look for topological changes. The class
    count "may increase multiple times, strictly increasing, never decreasing. Linear growth
    over many iterates is a clear bug though."
  - **D2** the alphabet is the minimal trellis's bridges only.
  - **D3** "the shortest walk on the dual graph gives the word, that is the algorithm" and a
    walk must always exist; the virtual `new1` class is "a very clear bug"; f(bridge) makes the
    word for singleton bridges that cannot walk. The period-3 path "seems buggy".
  - **D4** "Each anchor point should have two stable branches and two unstable branches …
    Each of those individual branches of course have a left side and right side." The left of
    one inversion branch is the right of the other; zones break the space up and must be kept
    in mind; cross-side is a flag only for the SAME branch, never for a bridge between two
    branches.
  - **D5** keep the period-k +1..+(k-1) iterate rule, fix the gap where table links are missing.
  - **D6** signs decide hole openings and rows; geometry only checks.
  - **D7** heteroclinic pairs give holes in B's partition.
  - **D8** simplify freely ("make as many choices as you'd like to make the design simpler …
    robust and professional"); public-API alignment is a later, separate plan.
- **Consequence:** WP0-WP12 in `docs/open_issues.md`, each run analyst -> planner ->
  test-engineer -> coder -> code-reviewer on its own branch. Claude's calls C1-C6 under D8
  (anchor-only branch = one element, element identity by bounding holes, det<0 out of scope,
  innermost-zone blasting, new-image-only same-stability count, per-call `max_points`) are
  revisable by the author.

## 2026-10-08 — DEC-19, the period-k growth stagger (raised by WP5)

- **Said:** "Accept the asymmetry" (offered: trim the trellis to complete orbits; warn and grow more; accept).
- **Consequence (D9):** a stable branch whose reference pair's +j images lie beyond the computed manifold has no fundamental-segment hole; that is correct for the finite trellis, not a bug. WP5's WARNING naming the branch stays. Laws that assume every branch has a window hole (e.g. the inversion anchor class leading its tangle) must be scoped to branches that have one, not left as known issues.

## 2026-10-08 — DEC-20, DEC-21, DEC-22 (raised by WP7 and WP8)

- **DEC-20, said:** "Cut S at the pip orbit." **Consequence (D10):** the partition and the dual graph use the stable manifold cut at the strong pip's orbit (chalkboard stage 3, as pseudoneighbors already do); this replaces the 2026-09-14 "the pip moves the unified set without trimming" rule.
- **DEC-21, said:** "I'm not sure this is an issue? We do not need computed bridges to find words, we strictly use the dual graph unless it is a singleton. That allows us to potentially find words even though we haven't computed the actual iterates. A word can include the class itself and any other classes." Confirmed on follow-up: "Yes, but mark them." **Consequence (D11):** a walk pair no computed bridge spans is a legitimate class: lettered, walked from its elements' landings, NOT a reliability failure, and marked "no computed bridge yet" (informational). Supersedes C9's "a pair no bridge spans is never a class".
- **DEC-22, said:** "Balanced trellis." **Consequence (D12):** cut the crossing set at an unstable depth matched to the stable cut at the pip, so blasting cannot reshape the fundamental segment at depth; the class count then never drops under blasting.

## 2026-10-08 — WP10 (numerics and performance leftovers, slice 1), Claude's calls under D3 / D8

Revisable by the author; recorded so later runs know they are calls, not the author's words.

- **C11 (under D3) — incomplete evidence gives no verdict.** A registered image chain whose
  word parts from the walk because the image passes the end of the computed stable branch
  (diagnosis `"stable end"`: two consecutive chain crossings on one stable branch with one
  sign) misses its crossings with the trimmed-away part, so it is incomplete: `word_agrees`
  is None, the diagnosis is kept and logged at INFO, and it never makes `verified` False or
  the dynamics unreliable. A `"walk"` disagreement stays a WARNING and a reliability failure.
  Effect: k28 at 6-7 blasts is reliable; at 8+ blasts the DEC-23 ties keep it unreliable.
- **C12 (under D8) — the law tier builds each case once per RUN.** `law_case` is session
  scoped (was module scoped: every one of ~12 law modules rebuilt all 11 cases, 40 s of
  setup); the fingerprint guard still fails a test that mutates the shared build.

## 2026-10-09 — DEC-23, equal shortest walks (raised by WP9)

- **Said:** chose "Non-crossing constraint" (offered: non-crossing constraint; registered image decides; stay ambiguous).
- **Consequence (D13):** when equal-length shortest walks read different words, the walk is enriched with the constraint that f(bridge) cannot cross the images of the unstable arcs bounding the lobe; the side is picked combinatorially, with no registered image needed. Applies to the heteroclinic lobe tie and the same-kind ties on k28 at >= 8 blasts.

## 2026-10-09 — WP9b (the non-crossing rule, D13), Claude's calls

Revisable by the author; recorded so later runs know they are calls, not the author's words.

- **C13 (under D13) — the marker rule applies to every walk with a member lobe, not only to
  ties.** D13 is worded for equal walks, but its physics (f is a homeomorphism, so the image of a
  member's lobe holds exactly the images of the markers the lobe holds) binds every walk. On the
  law case `heteroclinic` the class `bh` has a UNIQUE shortest walk (trivial, empty word) that
  holds the wrong markers and that its registered image contradicts; a rule restricted to ties
  would leave that contradiction (and `is_reliable` False) on a law case. Effect: `bh` reads
  `bi bj bi^-1`, its registered word, through a two-crossing detour.
- **C14 (under D3 / D13) — the walk runs between the landing SPANS.** D3's "from the landing of
  one end to the landing of the other" is read exactly: the walk leaves from and arrives at the
  stable edges the image span `f(X)` meets, not every edge of the iterated element (the hole
  bridges bounding `X` map onto the kept image bridges at the ends of `f(X)`, which `f(b)` cannot
  cross). This alone removes every tie on k28 at 8-13 blasts and k32 at 8-10 blasts, and is a
  no-op on every law case.
- **DEC-24 raised** (open): a member-less connecting class whose equal walks differ by a loop
  around another tangle, with no heteroclinic iterate link in N (`heteroclinic`'s `bn`). See
  `docs/open_issues.md`.

## 2026-10-09 — WP11 (simplification sweep, slice 1), Claude's calls under D8 / C3

Revisable by the author; recorded so later runs know they are calls, not the author's words.

- **C15 (under D8) — the topology reads no absolute cdist slack.** Anchors are told by
  identity (`Intersection.is_anchor`: synthetic at exactly (0, 0)), partition ownership,
  containment, the pip-orbit cut, the fundamental-segment overlap and `interval_at` are
  exact, and the only slack left is the RELATIVE error a scaled cdist carries
  (`SCALING_RTOL` at a scaled landing end, `collision_rtol` in `Trellis.same_point`).
  Measured before the change: forcing the `cdist_tol` that topology and loom read to 0
  (and to 1e-12) moved no output on any law case, deep run or blast series, and the full
  suite passed. `cdist_tol` stays for geometry only: degenerate bridge cuts in numerics
  and clipping a polyline to its arc. The cartoon identifies partition boundaries exactly
  too (they are registry floats).
- **C16 (under C3) — orientation is a FixedPoint fact.** `FixedPoint.preserves_orientation`
  is the sign of the solver's single-step Jacobian (`partial_jacobians[0]`, any period;
  None before a solve), and `Trellis.orientation_preserving` reads it; the topology never
  evaluates the dynamical system's Jacobian. The six side flips are one rule,
  `topology.side_after(side, steps, preserving)`, kept although det J < 0 stays out of
  scope (C3).

## 2026-10-09 — DEC-24, a member-less connecting class (raised by WP9b)

- **Said:** "Take whatever path is the most logical with our foundational setup. It is not that tricky. 'One class connects the two tangles and has no computed bridge' is a bit nonsense. If there is no computed bridge that connects the tangles then there is no bridge class which connects the tangles either. Treat every bridge the same. They just connect two intersection points. Those two intersections lie in stable partition elements. Side matters for those partition elements. A bridge connecting two partition elements on the same tangle is no different than a bridge connecting two partition elements on different tangles."
- **Consequence (D14):** no special case for connecting classes anywhere. A member-less connecting class must be traced to whatever produced it (the walk or landing that named the pair) and fixed there; read together with D11, a walk-named pair is a class like any other, but a pair across tangles that no bridge connects signals a construction error to root out, not a class to model.

## 2026-10-09 — WP9c (every bridge is the same, D14), Claude's calls under D13 / D14

Revisable by the author; recorded so later runs know they are calls, not the author's words.

- **Flag for the author: D14's premise does not hold on the fixture.** `bn` = {A:R_1, B:L_(2.0;1)}
  is no construction error: `named_by` gives `bn <- bl <- bi`, and `bi` has two COMPUTED members,
  `(91, 92)` and `(93, 94)` (arcs of W^u(p3 2.0) from p3's stable branch to the heteroclinic
  crossings on p1's window). `bl = f(bi)` and `bn = f^2(bi)` are member-less only because the
  images of 92 and 93 lie past p3's unstable cut (D12): D11's situation, and at p3 17 unstable steps
  they are computed. D14's rule ("treat every bridge the same") is what fixes it: the gap was that
  D13 had a rule only for classes with a member LOBE.
- **C17 (under D13 / D14) — walls.** D13's own wording, "f(bridge) cannot cross the images of the
  unstable arcs", for every class: a class's arcs lie on its members' unstable branches, or one map
  step past its namer's (D11); arcs on different unstable branches are disjoint, and so are their
  images, so the decided walk (`unique` / `trivial`) of a class walked EARLIER whose branches are
  disjoint from this class's is a wall, and a tie loses the walks whose chords interleave with a
  wall's on a face ring (`shortest_walks(walls=)`, `walk_chords`, `DualGraph.chords_cross`). A tie
  filter only (a unique walk is never second-guessed); when every walk crosses a wall the tie
  stands. Only a walk of multiplicity 1 is a wall or is judged by one: a walk standing for
  several face paths names none of them (review fix). Effect: `heteroclinic`'s `bn` resolved with the word base printed first, the dynamics
  reliable; at 17u the connecting class `bm` resolved to its registered word; nothing else moves
  on any law case, deep run or series step.
- **C18 (under D14) — what is tested.** D14 is tested as "a class connects two tangles only where
  a computed bridge in the trellis connects them" (the clause in `alphabet_closed_under_words`),
  not as "every connecting class has a computed member", which contradicts D11 (`bl`, `bn`).

## 2026-10-09 — WP10b (leftovers), Claude's calls

Revisable by the author; recorded so later runs know they are calls, not the author's words.

- **C19 (under D8; PK-4) — foreign chords pair in STABLE order.** A hole's image chain can
  leave its base branch for an excursion across another stable branch (another orbit
  branch or fixed point). The lobe H (bridge plus the stable segment between its ends)
  holds no point of any stable branch of the trellis, so f(H) holds no anchor; walking any
  OTHER branch from its anchor, the branch enters and leaves f(H) at successive crossings
  with the image chain (Jordan alternation), so the empty stretches there are its crossings
  paired in stable order. The base branch keeps the author's chain-order pairing
  (2026-09-23). The base is read at a non-anchor end (an anchor lies on every branch, D4),
  and the two "cannot happen" branches (a base across two branches; a non-anchor crossing
  abutting a stretch on another branch) are assertions. Checked by point-in-polygon on
  p3_6_blasts, p3_15, p3_16 and heteroclinic; changes refinement labels only (class tables,
  holes, words, refined rules and matrices identical).
- **C20 (under D1; NUM-12) — blast image curves are computed manifold.** The workbench keeps
  every image curve a blast computes and re-indexes it (query-only) in every recompute until
  growth has merged it into a branch; `create_bridges` cuts it like any curve, and one
  crossed fewer than twice is held whole as a partial piece. A recompute, a growth round or
  a zone (re)definition after a blast keeps every crossing, bridge, partial piece and
  `iterated` flag, and blast / recompute / blast equals blasting straight.

## 2026-10-09 — WP12 (re-verify and document), Claude's calls

Revisable by the author; recorded so later runs know they are calls, not the author's words.

- **C21 (under D3) — an order flip at a member end is noise, no verdict.** The twin of C11. An
  evidence (a registered member image) whose word parts from the walk, where a member end and a
  bound of that end's own element have registered `+1` images in the REVERSE stable order on
  one branch, is diagnosed `"order flip"`: the map preserves the order along a branch, so the
  registered cdists contradict the dynamics there (k=2.8 at 9 blasts: end 8.5e-6 outward of its
  bound, images 1.2e-6 the other way, both crossings resolution ratio ~1e4: a deep lambda-lemma
  copy hugging the hole-bridge end). `word_agrees` None, INFO, never a reliability failure.
  Lookup only. Effect: k=2.8 at 9-13 blasts reliable; k=3.2 from 9 blasts still unreliable
  (INT-11 kind 2: no registered preimage to prove a flip, the resolution floor, NUM-3/4).
  Law cases unchanged.
- **C22 (under D1 / D8) — growth un-trims a zone; restore never shortens.** A stable trim is a
  view: `grow_n_times` re-finds every branch's physical end (keeping a trim across growth would
  hide every crossing the growth was meant to discover, D1), and the topology never depends on
  the trim (the trellis is cut at the pip orbit, D10). `ResonanceZone.restore` leaves a branch
  that growth has run past alone (it used to cut it back: k=2.8 two blasts 362.9 -> 42.25). A
  zone that re-trims itself after growth is D8 API work.
- **Reporting only — the unmarked lobe.** A crossing on the image lobe of a hole with no
  stretch of its own marks nothing by the author's 2026-09-23 rule; it is now skipped as
  `"unmarked lobe"` at INFO instead of the "grow or blast" WARNING (the fold IS registered).
  No partition, class or word changes.
- **DEC-25 raised (BCG-8, open): the inversion series loses classes at blast 6.** The
  references `(23,24), (26,27)` are replaced by `(56,57)` alone, a clean witness (0.16 / 0.26 in
  stable cdist), no noise; its `+1` pair is never spanned because the reference bridge lies
  outside the zone (NUM-21), so the second stable branch has no window hole (D9's stagger):
  holes 9 -> 5, seed 9 -> 6, two classes ambiguous and contradicted. D12 does not stop it (the
  pip-orbit laps keep deepening). It contradicts D1's "a decrease is a bug" while following the
  pseudoneighbor definition. Options: (a) accept (the count converges as N deepens; D1 holds up
  to reference switches that are pure discoveries); (b) compute the reference arcs' images even
  outside the zone (needs growth, the image is partial); (c) refine old reference pairs rather
  than replace them (DEC-1's open sub-question). Also for the author: does D1 bind the CLOSED
  alphabet across a switch? (k=3.2 at blast 9: closed 45 -> 39 while the seed stays at 34 and
  the holes grow 56 -> 58.)
- **DEC-14 recommendation (the author's call): the default blasting cutoff 1e-7 defines the
  symbolic result.** At 1e-5 the seed decreases on k28, k32, p3 and nested_inner (E12); the k28
  drop is a lambda-lemma copy 5.6e-5 from a pair end that the coarse polyline orders on the
  wrong side, a pair-ordering resolution limit the per-crossing ratio (183) does not see. Coarse
  blasting cutoffs are for exploration only.

## 2026-10-09 — D15, holes propagate backward only on every orbit (reverses D5)

- **Author:** looking at the period-3 itinerary figures, "it seems like we are iterating the
  holes forward for some reason. Whereas in all other cases we only iterate the original holes
  backwards. That is the same for higher period cases."
- **Decision:** no hole at any forward iterate, the `+1 .. +(k_value-1)` appearances on the
  other branches of a period-k (or inversion) orbit included. The other branches get their
  holes from the backward chain, which walks the whole unstable branch cycle.
  `_is_forward_beyond_fundamental` -> `_is_forward_iterate` (`iterate > 0`); the PK-3 gap
  WARNING and `_missing_appearance` are deleted (a missing `+j` costs nothing now); law
  `no_direct_hole_beyond_fundamental` + `window_holes_on_every_branch` -> `no_forward_hole`;
  regression `test_period_k_window_gaps.py` deleted; `("walks_avoid_disjoint_images", "p3")`
  NOT_APPLICABLE (one decided walk left).
- **Measured:** p3 0-5 blasts 12 -> 8 holes, 9 -> 5 classes (the 4 lost were inert; words
  `a -> b -> c -> a u^-1 v^-1`); p3 8 blasts and p3_16 25/27 -> 23 seed, 47/49 -> 43 closed;
  p3_15, p3_6_blasts, nested_inner6 unchanged (their `+j` were never registered); inversion 0
  blasts 8 -> 6 holes, 6 -> 4 classes, every stable branch still holed. Suite 1117 pass.
- **Open (probed, not changed):** past 5 blasts the p3 class count is driven by the reference
  switch (BCG-4): the references become blast-born crossings at unstable cdist ~2600-3100 whose
  backward chain runs ~20 iterates (42 holes). DEC-25's "no window hole on the second branch"
  is no longer a gap under D15.

## 2026-10-09 — the alignment plan (`docs/alignment_plan_2026-10-09.md`, section E, Q1-Q13)

- **Said** (via the orchestrator, verbatim): "Do not drop the refined matrix and iterated partition.
  They are very important. The rest is fine."
- **Consequence — Q5 (the refined level) is a GOAL.** The audit's finding (refined children inherit the
  parent's word; the matrix is "assumed") is to be SOLVED, not removed: plan A.8d defines the iterated
  partition `P1 = P ∨ f(∂P)` with the author's empty-stretch rule (2026-09-23), the refined symbols as
  occurrences, and the refined word; states Theorem T (every member of a class has the same
  `P1`-level freely reduced itinerary), under which the inherited word is the IN-SPLITTING of the class
  graph (conjugate; identical sibling rows; children's columns sum to the parent's); slice B7 makes T
  and the in-splitting identities laws on every build and derives k=10's `[[1,1],[1,1]]` from `[[2]]`.
  Nothing of the iterated partition or refinement is deleted.
- **Consequence — every other question took the plan's default:**
  - Q1 (S-a): the pseudoneighbor test checks BOTH open intervals against X (the stable one too),
    reversing the 2026-09-15 "known asymmetry, no change" (slice B2; the series re-measured before
    DEC-25).
  - Q2 (N-a): the blast never drops a child holding a pip-orbit point; `m` stays "deepest registered"
    (B3). N-b (an explicit `depth`) not chosen.
  - Q3 (I-a): the minimal trellis maps EVERY hole bridge; inertness = transience in the class
    transition graph, read after the words; the evidence fixed point goes (B5). Measured: no
    class-level word moves on 14 builds; the iterated partition gets finer on 6.
  - Q4 (W-a): D3 stands (the walk is the word, the registered chain an oracle); search bounds become
    honest — a walk found under a truncated budget, a span meeting no node, members disagreeing on
    markers are `unresolved` (B6a); the evidence comparison moves out of the algorithm (B6b); the word
    is FREELY reduced (loops dropped, adjacent `c c^-1` cancelled).
  - Q6: the spectral radius stays unshown (2026-09-23 rule).
  - Q7: the manifold data model (arrays) is a later plan (B11).
  - Q8: D1 restated — at fixed depth `m` the hole set and class count are non-decreasing under
    blasting; a deeper `m` is a discovery; D1 does not bind the closed alphabet.
  - Q9: a class exists iff a computed bridge spans it or the word of an existing class names it
    (closure from the spanned classes; D11 and D14 read together).
  - Q10: recipe constants kept, their reasons written (B12).
  - Q11: the suite may pin ONE textbook result, the k=10 full 2-shift (the single exception to the
    2026-10-05 no-fixture-fact policy; B0).
  - Q12: `build_inversion_period2` takes `stable_steps=10` (class `e` resolves, same 6 classes /
    10 holes; B0).
  - Q13 (a): a detection the polylines cannot resolve — `resolution_ratio` below the threshold, or a
    neighbour of the same two curves closer along the unstable polyline than the chord deviation (a
    near-coincident run is one unresolved stretch, whatever its parity) — is registered and reported
    with the remedy "refine" and is NOT a crossing of N on any path; the blast-path parity vote
    (`Tangle._collapse_noise_crossings`) goes. The margin stays 3 and the corpus minimum (heteroclinic
    3.04) is recorded in the rubric (B1). Reverses "flags never discard" for this one case; NUM-4 is
    thereby decided.
- **Claude's defaults pending the author's objection (plan Q14-Q16):** an incomplete footprint (odd
  chain crossings on a branch) is data, never a reliability failure; a build contradicting Theorem T
  stops slice B7 and comes back to the author; no second refinement level is built.

## 2026-10-12 — B0 `oracle-laws` (alignment plan slice B0; plan `artifacts/plan_planner_human_2.md`)

- **Said** (via the orchestrator, verbatim): "Go ahead and run B0 through the planner and coder
  agents", then "Do you need to do more? If you do just do it".
- **Consequence — B0 approved as planned** (`plan_ok_planner_human.json`, iteration 1), with the
  planner's corrections standing: the hole-containment law is CLOSED at an end the iterate table
  links to a container crossing and strict elsewhere; no check-2 `KNOWN_ISSUES` (none failed on
  base); `inversion_period2` joins `every_active_class_walks`' NOT_APPLICABLE.
- **Plan question 1 took its default (yes):** the Q11 exception (the k=10 full-2-shift pin) is
  recorded in one sentence of `tests/conftest.py`'s no-pin paragraph.
- **Scope widened by one file** (orchestrator addendum, under the second instruction):
  `tests/regression/test_open_p3_deep_runs.py`, because its loop reads `LAYERS` and failed the
  k10 pin as vacuous on the deep builds. Rule adopted: `laws.LAW_CASES_ONLY =
  ("reference_pairs_clear_of_X", "k10_is_the_full_shift")`, skipped by the deep-run loop, so the
  scaling-only pair oracle never runs at depth (A.5). Six files in all, none under `src/`.
- **Test-engineer and code-reviewer stages ran on the orchestrator's judgement** (the second
  instruction); review passed in one round (`cycle_ok_test_engineer_code_reviewer.json`).
- **Not changed, left for the author:** check 3 bounds no distance between a carried end and the
  container polyline (reviewer's note; measured, see `artifacts/outcome_auditor_human_2.md`).

## 2026-10-12 — B1 `resolved-crossings`, review round 1 (author decision at the revise)

- **Said** (plan assay): approved the plan as written; the three questions take their defaults.
- **Said** (after the reviewer's `revise` on rubric 1-2): "Rule stands; widen scope to fix the four tests."
- **Consequence.** Q13 (a) holds everywhere: `N` is the resolved set, the trellis of the pip is cut at the
  deepest RESOLVED point of the pip's orbit, and an unresolved orbit point (k28 x8's f^8(q0) at ratio
  2.991; k32 x9 at 1.48; nested_outer x9 at 1.78) shortens `T(q0)` by one iterate. The remedy is
  "refine". The deep-run numbers that follow (k28 x8 24/15/23 -> 4/4/4; k32 x9-10 holes 58 -> 56,
  closed 39 -> 45, reliable) are reported, never pinned. The plan's file list is widened by the four
  files the failures live in: `tests/helpers/laws.py` (the pip-orbit-cut oracle reads `N`),
  `tests/regression/test_coarse_blasts_rows_by_signs.py`, `tests/regression/test_evidence_diagnosis.py`,
  `tests/regression/test_pip_orbit_cut.py` (rework their preconditions or builds, e.g. a finer
  `area_cutoff` so the detection resolves; never delete the defect they guard). Rubric 1's "none of the
  tests named in section 6 is edited" is amended to these four; rubric 3 lists 18 paths.

## 2026-10-12 — B1 `resolved-crossings`, outcome (review passed in round 2; branch `refactor/resolved-crossings` @ `76d0710`, base `1029fcc`)

- **Said:** nothing new since the revise-round decision above (the rule stands; scope 18 paths).
  The review passed in round 2 (`cycle_ok_test_engineer_code_reviewer.json`, iteration 2); no
  debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): N = `resolved_ids(threshold)`
  on seven readers; law cases unchanged; every deep-run change is one mechanism -- the pip orbit's
  deepest detection is unresolved (k28 x8 `f^9(q0)` id 245, ratio 2.991 -- the entry above calls it
  f^8, counted from the pip it is f^9; k32 x9 `f^10(q0)` 523, 1.48; nested_outer x9 `f^10(q0)` 447,
  1.78), is never linked, and `T(q0)` stays at the previous blast's trellis under every further
  blast until refined (k32 x9-10 and nested_outer x9-11 equal x8; k28 x8 equals x7).
- **Open for the author (from the audit):** (1) freeze vs placing an unresolved orbit point by
  `image_cdist`; (2) the two `threshold = None` regression builds (keep / NOT_APPLICABLE / unit
  re-guard); (3) `CLAUDE.md`'s "flags never discard" and k=2.8 x8 fixture facts stay stale until
  B12; `Trellis.py:204` default for the `MinimalTrellis` sparse snapshot (fix at `MinimalTrellis.py:411`).

## 2026-10-12 — B6a `walk-bounds-honest`, outcome (review passed in round 1; branch `refactor/walk-bounds-honest` @ `8a1cdc3`, base `76d0710`)

- **Said** (plan assay): asked for a plain description of the issue first, then approved the plan
  as written; question 1 takes its default (`lobe_image_markers` raises `MembersDisagree`, caught
  in `_resolve_itinerary` before any walk). No debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): a walk is a word only when
  established. A truncated search (any status), a landing span overlapping no stable node (new status
  `span_meets_no_node`) and member lobes disagreeing on a marker each leave the class unresolved with
  its reason; an unresolved class lends no wall. `_nodes_meeting` keeps a second tier (overlap > 0) so a
  span narrower than its own scaling error keeps the one edge it lies in (k32 blasts 6/8). No rule,
  reliability or closed-size change on any law case, deep run or series step; heteroclinic at
  `max_walks=2` becomes unreliable (was reliable on a truncated `unique`).
- **Left open:** export `MembersDisagree` from `tanglepack.topology` (out of scope); `CLAUDE.md` still
  describes the old "falls back to all of them" span rule; whether `SCALING_RTOL * |c|` is too coarse
  for narrow landing spans (numerics, later).

## 2026-10-13 — B8 `unresolved-sets`, outcome (review passed in round 1; branch `refactor/unresolved-sets` @ `ae18859`, base `8a1cdc3`)

- **Said** (plan assay, 2026-10-13): approved the plan as written -- eight files, the order-free
  collision rule, `forget_strong_pip` resets the cached trellis, the candidacy test runs after the
  zone is defined; question 1 takes its default (test candidacy inside `resonance_zone` as is,
  `add_resonance_zones` unchanged). No debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): an image claimed by two or
  more crossings links none of them (held link removed by `IntersectionRegistry.remove_iterate`), and
  a multi-match is recorded; both go to `registry.unresolved_links`, which the law
  `iterate_table_complete` keeps disjoint from the `+1` sources. `TangleSession.forget_strong_pip`
  exists; `resonance_zone` adopts its pip only when it is a candidate on the trimmed trellis (else one
  WARNING, previous choice kept). No state change on any law case or series step (`unresolved_links`
  empty everywhere but `near_tangency`); the workbench generation counter reads one higher per zone.
- **Left open:** candidates that make `topology()` assert once set (k10 2, k28 x1/x2 4, inversion
  2/14); `registry._generation += 1` written from `IterateInference`; the law's new clauses guarded
  only by a review probe; `check_iterate_table_complete` scans the registry, not N.

## 2026-10-13 — B4a `holes-algorithm-vs-diagnostics`, outcome (review passed in round 1; branch `refactor/holes-algorithm-vs-diagnostics` @ `f166573`, base `ae18859`)

- **Said** (plan assay, 2026-10-13): approved the plan as written -- eleven files, three commits
  (verbatim move, measured side on demand, tests), no behaviour change; question 1 takes its default:
  the inverse-map placement of propagated holes (`_map_backward`) stays in `StablePartition` until
  B4b, and the "no map in `StablePartition`" rubric moves to B4b. No debrief verdict yet. Not merged:
  the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): the 16 geometric readings
  live in `topology/diagnostics.py` (byte-identical move); `Hole` has no `measured_side` /
  `side_source`; the measured side is `diagnostics.measured_sides(trellis, holes)`, built only on the
  flag path (`punch_holes(diagnose=True)`), and `side_source(hole)` is a function. I1/I2 stay
  assertions in `StablePartition` as wrappers over cores in `diagnostics`. Every hole, measured side,
  row diagnosis, flag and fingerprint is identical on the 11 law cases and `p3_coarse_blasts`.
  StablePartition.py 1,947 -> 1,354 lines.
- **Left open:** `_map_backward` (B4b); stale references outside the file list (`numerics/geometry.py:7`,
  `scripts/refactor_proof/hole-side-and-own-blast/hole_sides.py`, `CLAUDE.md` holes bullet and `Hole`
  row); `measured_sides` cannot replay an orbit started from the propagation fallback (B4b removes it).

## 2026-10-14 — B4b `holes-are-combinatorial`, outcome (review passed in round 1; branch `refactor/holes-are-combinatorial` @ `ba23a1f`, base `f166573`)

- **Said** (plan assay, 2026-10-13): approved the plan as written -- nine files plus one
  `TangleSession` docstring line; a hole is its five combinatorial fields; placement moves to
  `plotting.hole_point`; question 1 takes its default (a): `hole_point` carries the reference point
  back with `trellis.dynamical_system.map_inv`, so figures are identical and `StablePartition` calls
  no map; the new WARNING watches scaled span ends only. No debrief verdict yet. Not merged: the merge
  is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): `Hole` is
  `(near_intersection_id, pair, bridge_side, bounding_ids, openings, iterate, origin)`; `punch_holes`
  punches references only and asserts no backward pair; propagation starts at the reference, no
  dedups, no fallback; `_hole_openings` names the anchor branch by `stable_key_at` plus the one-curve
  identification; `StablePartition.py` 1,091 lines, no numpy, no map; drawing is
  `plotting.hole_point` (which calls `map_inv`), also read by `measured_sides` and three laws. Every
  hole's five fields, fingerprint, word and measured side identical on the 11 law cases and
  `p3_coarse_blasts`; drawn points identical (0.0). Logs: `near_tangency`'s 14 "cannot propagate"
  WARNINGs gone; one new NUM-5 WARNING (fires on `k28_deep_blast_sliver`, no law case).
- **Left open:** `hole_point` raises a bare `StopIteration` on a propagated hole with no reference
  (unreachable through `Trellis.punch_holes`); stale "seeds already-punched orbits and regions"
  comments (`Trellis.py:1250`, `StablePartition.py:306`); I1/I2 cores in diagnostics (B4a);
  `CLAUDE.md` holes bullet and `Hole` row (B12).

## 2026-10-14 — B2 `pseudoneighbor-both-intervals`, outcome (review passed in round 1; branch `refactor/pseudoneighbor-both-intervals` @ `9cd29d5`, base `ba23a1f`)

- **Said** (plan assay, 2026-10-14): approved the plan as written -- one
  `_interval_is_clear(..., stability, ...)` with a min/max m-range, `per_step_beta(stability)`,
  candidates `own_intersection_ids` with the D7 guard mirrored; a pair kept only when both open
  intervals are clear; three unit tests; e13 over the seven E12 cases at both cutoffs (14 cells); the
  one question takes its default: e13's results jsonl is committed as the fourth file. No debrief
  verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): the pseudoneighbor test now
  matches the PDF definition (both open intervals against X, one function). No behaviour change on
  any law case or any of e13's 14 cells (0 stable fires in 182 steps). Seven unit-test fixtures whose
  data contradicted the definition were changed in data only (probe crossings moved to the other
  stable branch, or nudged); the reviewer accepted this against rubric 2's literal wording. DEC-25
  gets no new input: the inversion series still drops 2 -> 1 references at blast 6, with the stable
  clause playing no part.
- **Left open:** the m-range `min/max` is guarded only by an uncommitted hand check; two older tests
  (`walk_stops_at_reference_window_end`, `heteroclinic_candidate_ignored`) have no teeth for their
  named rule; `CLAUDE.md` "a crossing past a cut is a point of X" wording (B12); the old function name
  in `docs/` and `docs/open_issues.md`.

## 2026-10-14 — B3 `depth-of-n`, outcome (review passed in round 1; branch `refactor/depth-of-n` @ `7d3a294`, base `9cd29d5`)

- **Said** (plan assay, 2026-10-14): approved the plan as written, question 1 = (a) -- N-a as
  decided and restated: a child ending on the chosen pip's composed `+1` orbit is exempt from
  `min_separation` (one INFO), `BlastResult.pip_depth` composed (never `forward_depth`), the
  `blast_zone` docstring, `SeriesStep.pip_depth` measured test-side, the series law
  "non-decreasing" (question 2's default), the proximity-guard test restated, a new unit-test module
  on a multi-iteration call. Five files. The series stalls -- (ii) unresolved `f^{m+1}(q0)`, (iii)
  undetected, (iv) the zone filter -- left for a later decision. No debrief verdict yet. Not merged:
  the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): `m` = the registered `+1`
  chain from the chosen pip; `min_separation` never decides it. No change on any law case or on any
  step of the five series (pip_depth k28 `[1..8, 8]`, k32 `[1..9, 9, 9]`, nested_outer
  `[1..9, 9, 9, 9]`, p3 `[14..23]`, inversion `[5..9]`). Inside one multi-iteration call a
  pip-bearing child is kept: nested_outer 11 iterations ends at m = 9 (was 7), and that call now
  fires the existing 361 same-stability and 1 unresolved-detection WARNINGs (the exemption keeps
  near-coincident curve). One INFO per exempted child (one per single-iteration blast with the
  recipes' 1e-5). The child filter runs after every parent of the iteration is mapped (no effect of
  its own, instrumented by the reviewer). Suite 1193 pass / 1 xfail (KNOWN_ISSUES).
- **Left open:** the stalls (ii) on k28 / k32 / nested_outer / inversion, (iii) on p3 /
  heteroclinic, (iv) the zone filter -- options: refine (B1), place by `image_cdist` (PK-6-like,
  reverses B1 for the pip orbit, drifts toward N-b), a finer image near the anchor for (iii); whether
  the guard should win over the exemption when the curve is truly under-resolved; the stale B3
  section and decision N in `docs/alignment_plan_2026-10-09.md` (`:149-154`, `:903-943`,
  `:1166-1167`); `scripts/henon_blast_period_3.py` dead on `main` (pinned pip 12 not a candidate);
  the test-side orbit composition written three times.

## 2026-10-14 — B5 `inertness-is-transience`, plan round 1 sent back (author note at the revise)

- **Said** (`artifacts/plan_ok_planner_human.json`, status `revise`): "I believe you can determine if
  a bridge is inert before constructing the dual graph. Every bridge that has an initial hole punched
  in it should be inert. The zeroth holes."
- **Consequence:** the analyst tested the rule on all 15 builds (`artifacts/findings_analyst_planner.t2.md`).
  It holds everywhere: every class holding a zeroth-hole (reference-hole) bridge has an EMPTY word,
  and it is a theorem when no selected orbit sits inside the zeroth lobe and no other hole sits beside
  the bridge (the lobe holds one marker, its image holds none). It is sufficient, not the definition:
  heteroclinic `u` is a hole-free inert class. Re-planned: the minimal trellis maps every hole bridge
  EXCEPT one carrying a zeroth hole, decided from `holes` before the dual graph (no class table read);
  transience over the words (Q3, I-a) stays the verdict for every class; a new law
  `zeroth_hole_classes_have_empty_words` checks the theorem. Minimal trellis, `P1`, image bridges,
  refined classes and words stay as at `HEAD` on all 15 builds (round 1's "P1 finer on 8 builds" and
  its three new D9 WARNINGs are gone); only five verdicts change.

## 2026-10-14 — B5 `inertness-is-transience`, outcome (plan approved at round 2; review passed in round 1 of coding; branch `refactor/inertness-is-transience` @ `3fb310a`, base `7d3a294`)

- **Said** (plan assay, round 2, 2026-10-14): approved the round-2 plan as written -- the author's
  zeroth-hole rule (a bridge carrying a reference hole, iterate 0, is inert: its lobe holds one
  marker whose image is no marker, so its word is empty) decides, before the dual graph, which hole
  bridges the minimal trellis maps (every hole bridge except the zeroth-hole ones, read from
  `holes`); transience over the words is the verdict for every class (Q3, I-a); new law
  `zeroth_hole_classes_have_empty_words`; nine files (question 1's default). No debrief verdict yet.
  Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): the registered-evidence
  inertness is deleted (`_resolve_inertness`, `_collect_image_evidence`, `has_image_evidence`, the
  `image_classes` / `image_loops` / `unresolved` slots); `entry.inert` is `None` until
  `symbolic_dynamics` sets it from `inert_classes(words)`, then letters (active via `letter_for`,
  inert series); `Topology.classes` runs the dynamics; PK-8's WARNING gone. Words, minimal
  trellis, `P1`, refined classes and `is_reliable` unchanged on all 15 builds; five verdicts active
  -> inert (inversion_two_blasts `g`, heteroclinic `d` / `bj`, p3_6_blasts `ab`, nested inner 6
  `ae`); `unmapped` shrinks by those bridges. Suite 1193 pass / 1 xfail (KNOWN_ISSUES).
- **Left open:** the zeroth-hole rule is a sufficient condition, not the definition (heteroclinic
  `u` is a hole-free inert class); the unblasted nested fixture `p3_partitioned` has a zeroth-hole
  class (bridge `(4,5)`) that is UNRESOLVED and therefore active while its bridge is not mapped
  (outside the theorem's hypothesis); `inert_classes` not re-exported from `tanglepack.topology`;
  `CLAUDE.md` inertness text stale (B12).

## 2026-10-14 — B7 `refined-level-defined`, outcome (plan approved round 1 with question 1's default, option (ii); review passed in round 1 of coding; branch `refactor/refined-level-defined` @ `a6a4f68`, base `3fb310a`)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved): "Question 1 takes the default, option
  (ii): keep the base-branch chain-order chord pairing in B7, register
  (members_share_the_P1_itinerary, k32_series) as a strict-xfail KNOWN_ISSUE naming B7b, and plan B7b
  base-chords-in-stable-order as the next slice. Seven files; no number changes." No debrief verdict
  yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): `P1 = P ∨ f(∂P)` and the
  empty-stretch rule (base chords chain order 2026-09-23, foreign chords stable order C19) stated in
  the `PartitionFamily` module docstring; the odd-fold-count WARNING is gone, replaced by
  `IteratedHomotopyPartition.incomplete_footprints` (+ one INFO per build, a `describe()` line),
  collected by a sibling `_incomplete_footprints` because `_empty_stretches`' signature is
  monkeypatched by out-of-scope tests; Theorem T and the in-splitting stated in `SymbolicDynamics`
  (docstrings only, AST-equal); laws `members_share_the_P1_itinerary` (law cases only) and
  `refined_matrix_is_an_in_splitting`; the series run of T under its own id
  `series_members_share_the_P1_itinerary`, strict-xfail on `k32_series` (wp `B7b`). No number changes
  on the 15 builds; suite 1200 passed / 46 skipped / 5 deselected / 2 xfailed.
- **Deviations the author did not approve in advance:** an EIGHTH file
  (`tests/invariants/test_law_case_sanity.py`, +2/-1: the series id in `_STANDALONE_ISSUES`, the wp
  regex accepts `B\d+[a-z]?`), ordered by the orchestrator because the seven-file plan could not carry
  the approved known issue past that guard; the known-issue key carries the `series_` prefix, not the
  verdict's literal `members_share_the_P1_itinerary`.
- **Left open (for B7b or the author):** the k32 class is `{R_(0.0;8), R_(0.0;10)}` in ElementNaming
  (the KnownIssue text's "R_7 <-> R_9" is the positional `R#7`/`R#9`); the walked word agrees with the
  siblings and only the hole-bridge member `(6, 7)` disagrees; the module docstring's "laws on every
  build" overstates T (law cases + series only); B7b reverses the author's 2026-09-23 base-branch
  chain-order rule; commit trailers say `Claude Opus 5.5` though the orchestrator asked for another name.

## 2026-10-15 — B7b `base-chords-in-stable-order`, outcome (plan approved round 1 with both defaults; review passed in round 1 of coding; branch `refactor/base-chords-in-stable-order` @ `63477e3`, base `a6a4f68`)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved 2026-10-15): "Approved as written. Both
  questions take their defaults: question 1 no fold of _incomplete_footprints into _empty_stretches;
  question 2 yes, the geometric reading goes in tests/helpers/laws.py with a series test in the
  regression tier. Seven files. The base-branch chords pair in stable order, reversing the 2026-09-23
  chain-order rule on the base branch." No debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): ONE chord rule on every stable
  branch -- a lobe image's interior crossings pair in STABLE order, base branch included (C19's Jordan
  argument; `_empty_stretches` drops the `branch_key != base_branch` test and the dead
  `(base_branch,) = ends`). The 2026-09-23 chain-order clause is REVERSED. Evidence: on all 47 base-lobe
  readings of the corpus (16 builds + 47 series steps) the gaps alternate out/in from the anchor (base
  inside) and the covered gaps are the stable-order pairs; the one lobe where the orders differ is
  `k32_series` hole `(6,7)` at steps 2-10; the stop rule did not fire. Behaviour changes on `k32_series`
  only: iterated elements (2-5 23->24, 6-7 40->39, 8-10 113->111), `unmatched_members` 1->0 and Theorem T
  1->0 disagreements at 6-7, `refined_rules` at 8-10 (inert `v` 4 children -> 3, one token `v_4 -> v_3`).
  Letter-free words, `M`, `M_1`, all class counts, split counts and `is_reliable` identical on all 63.
  The request's "split inert -1 on 15 builds" was an artefact of B7's buggy scratch patch (it lost the
  `empty.setdefault(b, ...)` line). KNOWN_ISSUES shrinks to C3 alone. New test reading
  `helpers.laws.base_chord_misreadings` (point-in-polygon, clauses a-d) on 4 builds + 5 series. Suite
  1211 passed / 47 skipped / 5 deselected / 1 xfailed.
- **Left open:** the CLAUDE.md "Partition families" bullet still states base chords in chain order
  (B12); `_incomplete_footprints`' INFO log, `describe()` line and the class attribute doc still say
  "last" crossing (now the outermost); the new build test reads 4 of the 11 build lobes (the other 7
  read clean by hand); outside the k32 lobe the base tests cannot catch a DROPPED chord (the one inside
  gap is also a hole stretch), only a wrong one; commit trailers say `Claude Fable 5.1`.

## 2026-10-15 — B6b `word-finished`, outcome (plan approved round 1 with all six defaults; review passed in round 1 of coding; branch `refactor/word-finished` @ `65edec4`, base `63477e3`)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved 2026-10-15): "Approved as written. All six
  questions take their defaults: (1) WalkSearch.status has three values resolved/ambiguous/unresolved
  and 'trivial' belongs only to a class (empty word); (2) a search cut short with candidates becomes
  unresolved with reason 'truncated', walks kept, the truncated field goes; (3) class-level free
  reduction per A.8c as one free_reduce applied to cd.symbols and to the closure/inertness words,
  translate_itinerary stays unreduced so the T law keeps its same-occurrence reading, the registered
  evidence word is reduced by the same function before comparison; (4) members_disagree stays on the
  class only (six walk reasons); (5) MembersDisagree exported from tanglepack.topology; (6) one slice,
  ten files." No debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`):
  - "trivial" now means only an empty class word. `WalkSearch.status` takes the values `resolved`,
    `ambiguous` and `unresolved`. `WalkSearch.reason` holds one of six failure values (`unreachable`,
    `no_start_node`, `no_goal_node`, `no_admissible_walk`, `truncated`, `span_meets_no_node`). The
    `WalkSearch.truncated` flag is gone, and a search cut short by the budget keeps its walks but is
    `unresolved`.
  - `SymbolicDynamics.free_reduce` cancels at class level. It has three readers: `cd.symbols`, the
    closure/inertness words and the registered evidence word. `translate_itinerary` and the T law
    are unreduced.
  - The evidence comparison now lives in `topology/diagnostics.py`, behind
    `verify_words(records, trellis, naming, table, letters)`. `_compare_evidence` was folded into
    `verify_words`. `_member_pair` stays in SD.
  - No number changes on the 63 corpus items. Free reduction fires 0 times there.
  - Suite: 1223 passed / 47 skipped / 5 deselected / 1 xfailed (C3).
  - Line counts: SD 1537 -> 1312, DualWalk 1370 -> 1383, diagnostics 713 -> 990 (exactly at its cap).
- **Deviations the author did not approve in advance:**
  - `WalkReason` and `free_reduce` were also exported.
  - The per-class DEBUG line `class %s: %s, word %r, verified %s` was dropped.
  - A capped search with several walks logs the max_walks WARNING (now ending "the search decides
    nothing") and no longer also the ambiguous WARNING.
  - `plotting._class_status` lost its `-` fallback and now returns `cd.status`.
  - The loops unit test still pins the unreduced `translate_itinerary` and adds `free_reduce == []`.
  - The plan's rubric 7 said the capped tie reads `ambiguous` on base. It actually reads `unique`
    with `truncated=True`.
- **Left open:**
  - The `CLAUDE.md` symbolic-dynamics bullet and the `DualWalk` table row still list the old
    statuses and `WalkSearch.truncated` (B12).
  - `diagnostics.py` has zero headroom under its 990 cap.
  - `free_reduce` runs twice per class (closure, then `cd.symbols`), so a cancellation's DEBUG line
    would print twice.
  - `plotting._class_status` is untested.
  - The ambiguous, failure and truncated rows of the status mapping are exercised by unit tests
    only, because the corpus has none.
  - Commit trailers say `Claude Fable 5.1`.

## 2026-10-15 — B10 `dead-api`, outcome (plan approved round 1 with all eight defaults; review passed in round 1 of coding; branch `refactor/dead-api` @ `517aa01`, base `65edec4`)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved 2026-10-15): "Author approved the B10
  dead-api plan 2026-10-15 with all eight defaults": (Q1) rewrite the two `*_new_intersection_structures.py`
  demo scripts with `filter`, dropping their graph and array sections; (Q2) delete
  `TangleSession.classify_bridges` and rewrite its script and test as a comprehension over
  `classify_bridge`; (Q3 (a)) change only `henon_blast_period_3.py:123` and report the script dead;
  (Q4) delete the `from_branch` test with the method; (Q5) drop the `"_get_iterate"` string in
  `e02_gpu_end_to_end.py`; (Q6) leave the notebooks; (Q7) delete `IterateTable.as_forward_array` /
  `as_backward_array` in C5; (Q8) trim `topology.__all__` to 48 names and keep every import line.
  No debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`):
  - Deleted: `grow_x_times`, `_get_refined_point`, `_get_iterate`, `trim_stable_at_intersection`,
    the registry's networkx graph, `graphviz.py`, `build_intersection_graph`,
    `visualize_intersection_graph`, `on_cdist_range`, `from_branch`, `nearest_by_unstable_cdist`,
    both `as_*_array` pairs, `classify_bridges`, `BlastResult.terminated_early`. Kept:
    `release_manifold` (live), every law-tier oracle, `refine_manifold`, `filter` /
    `from_fixed_point` / `iterate_orbit`.
  - `import tanglepack` no longer loads networkx. `topology.__all__` 104 -> 48.
  - No number changes on the 63 corpus items (probe byte-identical). Suite 1215 passed / 47 skipped /
    5 deselected / 1 xfailed (C3); 8 test ids gone as planned, 1 renamed.
  - Net −1013 lines over src + tests (−1034 with scripts and stress).
- **Deviations the author did not approve in advance:** an unused `TYPE_CHECKING` import dropped in
  C3; the `Iterable` import went in C4, not C5; `remove_iterate` no longer reads `target_id`; a few
  extra docstrings reworded so none mentions a registry graph; the comment in
  `henon_bridge_classification.py` moved (2 lines, not 1).
- **Left open:** CLAUDE.md rows `:133` (`graph(bridges=)`), `:136` (`graphviz`), `:154`
  (`terminated_early`) and `docs/source/tanglepack.numerics.rst:183-189` (automodule of the deleted
  `graphviz`) are stale (B12); the 56 F401s in `topology/__init__.py` (B12's import block);
  `scripts/manifold_cutting_test.ipynb` and `notebooks/archive/period_three.ipynb` still call
  `grow_x_times` (Q6); `henon_blast_period_3.py` still dead at `:95`;
  `test_map_step_and_graph.py` keeps "graph" in its name; commit trailers say `Claude Fable 5.1`.

## 2026-10-15 — B9a `session-three-concerns`, outcome (plan approved round 1 with every default Q0-Q6; review passed in round 1 of coding; branch `refactor/session-three-concerns` @ `5f0db5f`, base `517aa01`)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved 2026-10-15): "Author approved the B9
  plan 2026-10-15 with every default": (Q0) split B9 into B9a `session-three-concerns` now and B9b
  `one-element-at-end` next; (Q1) no thin plot delegate stays on the session; (Q2 (a))
  `cartoon_zones` moves verbatim to `loom/cartoon.py` and reads the STORED zones -- zones are NOT
  computed from `Trellis.pips` (that would change `inversion_period2`'s circle cartoon and pull the
  zone boundary builder and `ring_cycles` / `ring_order` into topology); (Q3) `plot_resonance_zones`
  goes to the same module, taking zones; (Q4) the describe fan-outs (`describe_holes`,
  `describe_stable_partitions`, `_describe`) are deleted; (Q5) the
  `test_a_pip_trimmed_away_is_not_replaced_by_plotting` test is deleted (`Trellis.plot_*` never
  touches the pip memory); (Q6) no new "element on any side" API. No debrief verdict yet. Not
  merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`):
  - `TangleSession.py` 1744 -> 843 lines; it holds no `plot_*`, no cartoon geometry, no fan-out,
    no `partition_element_for`, no describe fan-out. `cartoon_zones` / `plot_resonance_zones` live
    in `loom/cartoon.py` (exported from `tanglepack.loom`); `cartoon_enclosures(trellis, bridge_ids,
    parents, fixed_points)` in `topology/plotting.py` (not in `topology.__all__`, still 48 names).
  - Callers draw a circle cartoon in three explicit steps (`cartoon_zones` -> `cartoon_enclosures`
    when `parents` -> `plotting.plot_dual_graph_cartoon`); per-fixed-point plots go through
    `session.trellis(fp).plot_*`; manual stage calls through `session.topology(fp)` or the trellis.
  - No number, figure or script stdout changes (63-item probe byte-identical, 32 canvases, 42
    itinerary PNGs, 11 scripts). Suite 1208 passed / 47 skipped / 5 deselected / 1 xfailed (C3);
    7 test ids deleted, 3 renamed, as planned. Net -749 lines (src -613, tests -143, scripts +7).
- **Deviations the author did not approve in advance:** imports pruned in the commit that stops
  using them; the `plotting.py:6` docstring fix in C2, not C1; small docstring edits in listed
  scripts and `test_partition_elements.py`; a test-engineer commit strengthening
  `test_session_caches._repartition_by_hand` (pre-existing toothless assertion).
- **Left open:** B9b (the element-at-end fold); stale docstring citations of
  `TangleSession.partition_element_for` at `Trellis.py:1413`, `TopologyResults.py:427` (B9b) and of
  "the nested-zone classifier in TangleSession" at `ResonanceZone.py:295`; CLAUDE.md `TangleSession`
  row (B12); `tanglepack.loom.plot_resonance_zones` untested; `henon_blast_period_3.py` still dead
  at `:95`; commit trailers say `Claude Fable 5.1`.

## 2026-10-15 — B9b `one-element-at-end`, outcome (plan approved round 1 with every default Q1-Q5; review passed in round 1 of coding; branch `refactor/one-element-at-end` @ `94fca6b`, base `5f0db5f`)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved 2026-10-15): "Author approved the B9b
  one-element-at-end plan 2026-10-15 with all five defaults": (Q1 (B)) a function-level import of
  `HomotopyPartition` in `bridge_classes`, which calls `owner_at_end`; `_element_at_end` deleted
  (not option A's two shared helpers beside `row_of_end`); (Q2) the class rule's unowned-crossing
  AssertionError becomes `owner_of_intersection`'s ValueError on that unreachable path; (Q3) the
  three tests citing `element_for` are renamed; (Q4) no `plot_resonance_zones` smoke test here
  (B12); (Q5) `ResonanceZone.py` untouched (its two `TangleSession` citations are accurate). No
  debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`):
  - `PartitionFamily.owner_at_end` is the one reader of "which element owns this bridge end on its
    row", on every computed path; the class rule's `ends` calls it. `Trellis.element_for`,
    `StablePartition._elements_at_bridge`, the `StablePartitionResult.elements_at_bridge` field
    (with its two writers) and `BridgeClass._element_at_end` are deleted. The any-side reading is
    `owner_of_intersection` on the `stable_key_at` branch; no "any side" API added (B9a Q6).
  - No number changes: agreement probe 7,090 / 14,180 / 1,776 with readings byte-identical to
    base; 63-item outcome probe byte-identical. Suite 1207 passed / 47 skipped / 5 deselected / 1
    xfailed (C3). One test deleted, three renamed. Net -168 lines (src -149, tests -19).
- **Deviations the author did not approve in advance:** docstring wording beyond the analyst's
  candidate (the `bridge_classes` Note, the singular "it" in the `StablePartitionResult` Note);
  the missing-partition ValueError text is now `owner_of_intersection`'s (type unchanged, still
  names the branch); commit trailers say `Claude Opus 5.5` (the spawn asked for `Fable 5.1`).
- **Left open:** CLAUDE.md `Trellis` row still lists `element_for` (`:139` in the working tree,
  B12); `plot_resonance_zones` untested (B12); the unused `Intersection` import at
  `StablePartition.py:238`; `henon_blast_period_3.py` still dead at `:95`; four stale committed
  higher-period PNGs; the rewritten `test_partition_elements.py` does not catch a broken C7 anchor
  rule or a wrong row in `owner_at_end` (the law tier does).

## 2026-10-15 — B12a `one-design-note`, outcome (plan approved round 1 with the B12a/B12b split and every default Q1-Q12; review revise in round 1, pass in round 2; branch `refactor/one-design-note` @ `b39c7e2`, base `94fca6b`)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved 2026-10-15): "Author approved the B12
  plan 2026-10-15 with the split (Q11: B12a now = commit 0 + design note + src de-citation +
  CLAUDE.md; B12b next = laws, smoke test, __init__ trim, rst lists, chalkboard, Q10 reasons, PNG
  regeneration, plan Outcome) and all twelve defaults": (Q1) the 649 MB `figures/henon_period3_*`
  sweeps are NOT committed, their two scripts and CLAUDE.md lines are; (Q2) B12b deletes the 9
  TAUTOLOGY laws, keeps `dual_bipartite_degree`, every MIXED and ORACLE law (-4 ids); (Q3) B12b
  regenerates only the 3 content-stale PNG folders (`heteroclinic_p1_p3`, `nested_p1_p3_inner6`,
  `inversion_blast2`) in their own commit; (Q4) the note is `docs/design_note.md`, sections
  A.0-A.11 (+ H, S, O, I); (Q5) CLAUDE.md's physics sections stay byte-identical, their ids
  (`WP10`, `NUM-20`, `WP10b`, a date) exempt; (Q6) the 8 log strings (7 ids + 1 date) stay in
  `src/`; (Q7) B12b: `topology.__all__` -> 51 with `MembersDisagree`, `WalkReason`, `free_reduce`,
  the other 53 unused imports dropped; (Q8) the other id families (`B*`, `Q*`, `E*`, `F7`, `N-a`)
  and dates leave `src/` too, "design note §A.x" stays; (Q9) Sphinx: a list check only (+12
  modules, -graphviz, `conf.py` path `../../src`), no `-W` build; (Q10) the chalkboard's 17 rows
  closed in place in B12b, the 9 audit divergences closed without a row listed in the note's H now;
  (Q11) two branches, B12a first; (Q12) the recipe reason for zones at `f(q0)` is the measured
  equivalence (zone at `f(q0)` after n blasts = zone at `q0` after n-1, k28/k32, n <= 2), written
  in B12b. No debrief verdict yet. Not merged: the merge is the author's.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`):
  - `docs/design_note.md` (722 lines) is the design: one section per stage, every one of the 86
    ids `src/` cited at base anchored, plus H (slices, WPs, superseded rules), S, O, I.
    `src/` keeps only the 7 log-string ids and the 1 log-string date; 15 "design note §A.x"
    pointers. Stripped AST equal on 51 files; src -60 lines (+441/-501), 35 files (34 planned +
    `FixedPointSolver.py`, a recorded deviation the reviewer asked for).
  - CLAUDE.md: "Bridge classes and the dual graph" replaced by "The pipeline, stage by stage" (one
    paragraph per stage, each naming its functions and ending "design note §A.x"); table at HEAD;
    80,246 -> 62,737 bytes (173 -> 185 lines); `:1-27` and Coding Style byte-identical.
  - Suite 1207 passed / 47 skipped / 5 deselected / 1 xfailed (C3), collected ids identical.
- **Round-1 defects** (caught by reading, no script could): the note said `image_of` maps the
  boundary with the real map, listed `from_table` as deleted by B4b (B4b added it), said the
  in-splitting runs on the series (only T does), and that the minimal trellis ignores the class
  table; plus history sentences in `FixedPoint.py` / `FixedPointSolver.py`. All fixed in round 2.
- **Left open:** B12b (laws, smoke test, `__init__` trim, rst, chalkboard, Q10 reasons, PNGs, the
  plan's Outcome section); old-code narration outside the rubric's word list (`BridgeIterator.py`
  "plan row 1.1/2.2 removed", `BaseManifold.py:16-17`, `ManifoldMachine.py:632,700-701`,
  `Tangle.py:593,1017-1018`, `DualGraph.py:50-51`, `geometry.py:9-11`); the note and CLAUDE.md say
  "no fan-outs" though `TangleSession.classify_strong_pips` still fans out over fixed points (kept
  by B9a for pip memory); note cosmetics (`:506` unwrapped, `:377` broken line, `:695` "990-line
  cap" vs 988); commit trailers say `Claude Opus 5.5`. These notes were appended in the MAIN tree
  after commit 0 copied it: the branch's `.refactor/` is one entry behind main's working copy.

## 2026-10-16 — B12b `laws-and-loose-ends`, outcome (the campaign's LAST run; plan approved round 1 with all six defaults Q1-Q6; review revise in round 1, pass in round 2; branch `refactor/laws-and-loose-ends` @ `e163f7a`, base `a550a5d`, nine commits)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved 2026-10-16): "Author approved the B12b
  laws-and-loose-ends plan 2026-10-16 (the campaign's last run) with all six defaults": (Q1) the
  design note's "Verified against" line goes in a final docs-only commit naming the last commit
  that touches code, tests or figures; (Q2) the other 13 figure folders are NOT regenerated (all
  pixel-identical; three differ only in the PNG `Software` chunk); (Q3) no reproducible-SVG
  machinery: figures are judged by PNG bytes, normalised SVGs and pickles re-rendered to PNG;
  (Q4) `CLAUDE.md:162` and `:101` fixed as well as `:60`, one line each; (Q5) `outer_blasts=2`
  keeps the author's phrase "chosen to make the case closed" plus the measured distinction (one
  outer blast resolves every class through the trellis path, two make every outer class walk);
  (Q6) the coarse-blasts regression's call to the deleted `check_openings_follow_signs` is
  REPLACED by `check_partition_covers_branch`, not dropped.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`):
  - Nine TAUTOLOGY laws deleted whole with their `LAYERS` entries and two imports only they used
    (`laws.py` 106 -> 97 laws, -209 lines); `dual_bipartite_degree` and every MIXED/ORACLE law
    unchanged by AST. Collected ids 1260 -> 1257 (the four `refined_children_inherit_word`
    set-aside skips out, `test_plot_resonance_zones_shades_every_zone` in). Suite 1208 passed /
    43 skipped / 5 deselected / 1 xfailed (C3), no XPASS.
  - `topology.__all__` = 51 (adds `MembersDisagree`, `WalkReason`, `free_reduce`), 53 unused
    imports dropped (`__init__` 195 -> 133 lines); side effects: `tanglepack.topology.ElementNaming`
    / `.PartitionFamily` are now the submodules, `import tanglepack` no longer loads `StrongPip`.
  - Stripped AST equal on 50 of 51 `src` files (the 51st is `__init__`); 63-item probe equal.
  - rst automodule set == module set (51); `conf.py` path `../../src`; chalkboard closed in place
    (5th column); alignment plan's Outcome section (17 rows); note "Verified against `ad1d96a`".
  - Three figure folders regenerated (51 files, `partitions.png` byte-identical so untouched);
    PNGs byte-reproducible with system `python3` (18/18).
- **Round-1 defects** (caught by reading, not grep): `henon_cases.py` "real (``unique``) walks"
  (no such status since B6b; the word came from base through the plan's text) and `geometry.py`
  "pure functions on plain arrays ... one tolerance" (two of four take a manifold / bridge; four
  different tolerances); plus "pre-2.2" left dangling in `BridgeIterator.py`, the note's index
  still "B0-B12a", and the Outcome intro calling every hash a branch tip. All fixed in round 2 as
  two appended commits (`ad1d96a`, `e163f7a`), history not rewritten.
- **Left open after the campaign:** `topology/plotting.py:1947` stale ``unique`` status wording;
  the `plot_resonance_zones` smoke test pins the handle count only;
  `scripts/henon_symbolic_itineraries_period3.py:23-27` stale period-3 facts;
  `scripts/henon_blast_period_3.py:95` dead (standing rule: reported only); DEC-25; B11 deferred;
  `ManifoldMachine.py:627` unused `final_node`; `BaseManifold.py:11` stray to-do; ignored
  `src/tanglepack_gui/`, `src/tanglepack_webdash/` (only `__pycache__`); `conf.py:35` `_static`
  missing; no Sphinx build without `myst_parser` / `furo`; commit trailers mix Fable 5.1 / Opus
  5.5. No debrief verdict yet. Not merged: the merge is the author's.

## 2026-10-21 — `blast-steps-compose`, outcome (not a slice of the closed alignment plan; plan approved round 1 with every default except Q6; review passed in round 1 of coding; branch `refactor/blast-steps-compose` @ `a6d2ffe`, base `c1ca6d0`, six commits)

- **Said** (`artifacts/plan_ok_planner_human.json`, approved 2026-10-21): every default (Q1 drop
  remembered on the bridge as `Bridge.dropped_at`, carried by `rebuild_bridges`; Q2 re-admitted
  only at a strictly smaller separation or None; Q3 numbers may move as discovery; Q4 a
  composition test; Q5 one `BlastResult.skipped_dropped` count, no second list; Q7 depth stall
  diagnosed, not fixed; Q8 the k=2.5 rerun and legends after the merge; Q9 cloud rule B) EXCEPT
  Q6: "Q6 I think we should tie them together. If things are getting smaller and smaller scale we
  need to refine more to get good manifolds. I get what you're saying about it dropping more
  children though too. I don't imagine it will be a big issue though. If it is, then find a
  simple reasonable solution".
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): `n` single-iteration
  blasts equal one `n`-iteration blast (bridge ids, `iterated`, drops, `pip_depth`); no law-case
  number moved; only `k32_series` 9-10 and `nested_outer_series` 8-11 bridge / N counts fell;
  one-call blasts with drops slightly stricter (k28 n=11 557 -> 554). Q6 landed as an OPT-IN
  `min_separation="resolution"` = `1e-5 * (area_cutoff / 1e-7) ** (2/3)`; default stays None and
  no recipe uses it. The segment-distance fallback did not fire.
- **Open for the author:** whether Q6's tie stays opt-in (a), replaces the 1e-5 constants in the
  recipes (b), or becomes what every guard-wanting caller gets (c). Record the choice here when
  made. Not merged: the merge is the author's.

## 2026-10-23 — N0 `numerics-arrays-skeleton`, outcome (FIRST run of the new numerics-arrays campaign, not a slice of the closed alignment plan; plan approved round 2 of planning with every default; review passed in round 1 of coding; branch `refactor/numerics-arrays-skeleton` @ `41056f8`, base `a8bb8a0`, five commits)

- **Said** (request, 2026-10-23): D1 skeleton then split, one slice per run; D2 a tracked profiler
  run before and (unchanged) after; D3 byte-identical where it can be, every law always, a measured
  and explained difference allowed, moved topology listed and judged; D4 index-based public API but
  human readable (`bridge.head` / `bridge.tail`, `manifold.root` / `tail`, `len`, `.points`,
  `.cdists`); D5 `BranchPoint` stays a class with its API, touched only in its own slice, last or
  near last; D6 arrays throughout, an array namespace, serialisable state, no scheduler, GPU opt-in.
  Plan gate (`artifacts/plan_ok_planner_human.json`): "Approved as written; all Q1-Q9 defaults
  accepted (Q8: same-stability count joins the campaign as N2; Q6 order N1..N10)."
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): committed on the branch
  `stress/profile/{profile_numerics.py, profile_figures.py, identity_probe.py, README.md}`, the
  baseline `figures/numerics_profile/a8bb8a0{,_scale1e-9,_scale_traced}/` and
  `docs/numerics_arrays_plan_2026-10-23.md` (N1 slim-segment-index, N2 same-stability count, N3
  point store, N4 order arrays, N5 brackets on ids, N6 readers on views, N7 delete Point, N8
  BranchPoint alone, N9 array namespace, N10 save/load + after-profile). No `src/` or `tests/` change.
  The 45 GB peak that motivated the campaign is NOT on `HEAD` (blast-steps-compose fixed it): the
  k = 2.5 series plateaus at 2.9 GB (1e-7) / 6.5 GB (1e-9); the segment index's Python bookkeeping is
  57 % of held memory, `Point` objects 41 %, rtrees ~0; 82-87 % of blast time is the same-stability
  count, not linked-list work.
- **Carried forward (reviewer's non-grading notes):** the `identity_probe.py` Dev Note over-claims
  ("every read of a library object goes through one reader"; registry fields and topology
  `describe()` are read directly): reword at merge or in N1's diff. N10 must own giving
  `profile_figures.figure_stage_peaks` (and the by-component census panel) a baseline series: add
  `stress/profile/profile_figures.py` to N10's Files (doc fix at merge or in N10's plan).
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-24 — N1 `slim-segment-index`, outcome (slice 1 of the numerics-arrays campaign; plan approved round 1 of planning with every default; review REVISE in round 1 of coding, PASS in round 2; branch `refactor/slim-segment-index` @ `8ded634`, base `41056f8`, eleven commits)

- **Said** (plan gate, `artifacts/plan_ok_planner_human.json`): "Approved as written; all defaults
  accepted (Q1 (a) slim id()-keyed edge set; Q2 nothing beyond the curve; Q3 (i) same count; Q4 keep
  only paired records; Q5 branch-cut bridges keep claims; Q6 two scenarios in identity_probe.py; Q7
  rename to indexed_segments)." D1-D6 / Q1-Q9 of the campaign plan bound the run unchanged.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): a query-only curve keeps one
  int edge key (`Tangle._query_edges`) and a `_Segment` only when its query paired; census at the 1e-7
  plateau 1914.0 -> 990.6 MiB, segment-index band 1045.4 -> 122.0 MiB; blast time flat against the
  same-day base (399.05 -> 399.57 s); probe identical on 62 dumps; suite 1221 passed, C3 alone. Round
  1's indexed path re-indexed query-only edges on `reset=False` (16 vs 22 crossings on k28); fixed by
  one guard in `_insert_segment` (the base's extra claim dropped, no observable effect). Accepted
  behaviour changes: later sids shift when a separator splices into a query-only bridge (sid values
  decide nothing); `_query_edges` never forgets on `reset=False` (Q1's risk); last-bit differences
  after a `reset=False` pass; `add_manifold` returns `set[int]`; `populate_intersections(new_sids)`.
- **Carried forward to the merge:** two reviewer wording slips (the Outcome row's peak RSS is
  absolute, the doc's convention is over-baseline; "nor claims ... as before N1" at `Tangle.py:36`,
  `:666-667` and in the row is false about the base). **To later slices:** the stress report still
  reads e03's old `segments` key; the pre-N1 `reset=False` under-detection (re-query or retire).
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-24 — N2 `same-stability-count-new-only`, outcome (slice 2 of the numerics-arrays campaign; plan approved round 1 of planning with every default; review PASS in round 1 of coding with recommendations R1-R4; branch `refactor/same-stability-count-new-only` @ `6bd37cf`, base `0b1db28`, six commits)

- **Said** (plan gate, `artifacts/plan_ok_planner_human.json`): "Approved as written; all defaults
  accepted (Q1 (a)+(b): count item in identity_probe.py always on with diagnose=True, and the count
  in each profile_numerics.py scale row; Q2 new by curves, unchanged; Q3 orientation kept, no box
  pre-filter; Q4 randomised property test; Q5 keep _PAIR_BLOCK chunking; Q6 reword the cost
  docstring)." D1-D6 / Q1-Q9 of the campaign plan bound the run unchanged.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): `_candidate_pairs` with `new`
  pairs an old segment only with the new positions in its forward range (`pool`/`start`/`count`); the
  pairs equal the base's filtered pairs element for element (auditor 3,750 calls, TE 4,896, reviewer
  1,200); probe identical on 62 dumps incl. the new `same_stability_counts`; count fields equal on all
  32 scale blasts. Sequential same-day profile: `diagnose` 360.3 -> 53.7 s, blast 415.8 -> 109.2 s,
  peak RSS over baseline 2396.2 -> 2374.0 MiB, census 990.6 MiB unchanged. No library behaviour
  change. The probe now runs every blast with `diagnose=True` (inherited by N3-N10).
- **Carried forward to the merge (text only):** the Outcome row quotes the CONCURRENT profiler pair
  (2248.6 -> 2427.2 MiB, a rise that does not reproduce) and a false "-488 to +404 MiB" range (measured
  -727 to +404) and a misleading "traced peak 1552.0 -> 1544.8" (run-to-run; ~1.55 GB on both): swap
  in the sequential pair (reviewer R1). **To later slices:** run profiler pairs sequentially; peak RSS
  is set by the count's unblocked narrow phase (12.65 M candidates at blast 13, ~1.55 GB traced), not
  by held state, so census drops (N3) will not show in peak RSS; blocking the narrow phase is a
  candidate small slice outside the plan. N3 request drafted in the outcome §10.
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-24 — N3 `point-store`, outcome (slice 3 of the numerics-arrays campaign; plan approved round 1 of planning with every default; review PASS in round 1 of coding with notes (a)-(e); branch `refactor/point-store` @ `5ce5685`, base `138b34c`, seven commits)

- **Said** (plan gate, `artifacts/plan_ok_planner_human.json`): "Approved as written; all defaults
  accepted (Q1 one persistent Point per real row, id() sites untouched; Q2 BasePoint.__slots__ = ()
  and rewrite the walker test; Q3 no BranchPoint row, four [branch_point] tests rewritten
  BranchPoint<->BranchPoint; Q4 drop edist entirely; Q5 the workbench owns point_store, machine
  appends to the curve store, loose_store for hand-built points; Q6 merge-dropped duplicates and
  fictitious chain are ordinary rows; Q7 append-only until N10, census changes now; Q8 phantom-only
  tie unit test; Q9 get_point returns a read-only view)." D1-D6 / Q1-Q9 of the campaign plan bound the
  run unchanged.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): one append-only `PointStore`
  per workbench (`numerics/PointStore.py`: xy, cdist NaN-none, image / preimage int64 -1-none, x1.5
  growth, `nodes` row -> Point); `Point` slotted over its row (100 B object + row int, was 408 B);
  phantoms are rows with no Point (778,302 at the 1e-7 plateau, 285.0 MiB of objects gone);
  `_cache_preiterate` deleted; the tie rule reads `has_iterate_link()` on both heads. Census at the 1e-7
  plateau 990.6 -> 398.5 MiB (all point bytes 829.6 -> 237.4 MiB, ~151 B per real point all-in); peak
  RSS over baseline 2338.8 -> 1883.6 MiB (sequential pair; TE's pair 2483.2 -> 1891.4); blast +2.8 %
  (TE +3.4 %), `walk` +1.0 s, `infer_iterates` +2.1 s. Probe identical on 62 dumps (auditor re-dumped
  HEAD); ties record-equal (17 / 2 / 2, keeper head_1); suite 1234 passed / 43 skipped / C3 alone
  (auditor rerun). Accepted behaviour changes are attribute-level only: no `edist`, no `__dict__` /
  `_coords`, read-only x / y and view, `TypeError` on a cross-store or non-Point link, `loose_store()`,
  the fictitious tail kept (+2 rows), the "already exists" check before any link write, and (not in
  the plan's section 5) a stored NaN coordinate or cdist reads as None.
- **Carried forward to the merge (text and one-liners):** "preimage" for "preiterate" in three places
  (`PointStore.py:13`, `ManifoldMachine.py:777`, `Point.py:134`; on a stable curve the phantom is an
  image, 187 of 1457 on k10); `Point.py:11` "one Point per row" is false for phantom rows;
  `_get_preiterate` returns a writable view (`ManifoldMachine.py:814`: `store._xy` -> `store.xy`); the
  Outcome row's "`edist` and `_coords` are gone from `Point` (no reader)" (`_coords` had three reads).
  Optional: the "already exists" check above `store.append` (`:362` / `:370`); the PointStore Dev Note
  on which rows "die". **To later slices:** `profile_figures.MEMORY_GROUPS` has no `point_store`
  group (eight fixed slots); public names for the store columns (the machine and Point read `_xy`,
  `_image`, `_preimage`); the NaN-reads-as-None change; the per-step same-store guard in
  `_iterate_manifold` (reviewer (e)); CLAUDE.md's data-model paragraph and the absence of `PointStore`
  there (N4 makes it false). N4 request drafted in the outcome §10.
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-25 — N4 `curves-on-order-arrays`, outcome (slice 4 of the numerics-arrays campaign; plan approved round 1 of planning with every default Q1-Q13; review REVISE rounds 1-2 for prose only (engine accepted round 1), PASS round 3; branch `refactor/curves-on-order-arrays` @ `ab7bd20`, base `8fc3933`, nine commits)

- **Said** (plan gate, `artifacts/plan_ok_planner_human.json`): "Approved as written; all defaults
  Q1-Q13 accepted: one run; curves as views over chains; two int32 store columns _chain/_position;
  keeper rule reproduced exactly incl. branch (b); (f)-(h) raise AssertionError; BranchPoint
  insert_point_* set the chain anchor (D5: nothing else in the class); role 5
  (_branch_forward/_branch_backward) retired in N4 on the 0-through-walk evidence, departing from
  campaign Q9 which gave it to N8; public read-only store columns plus writers; NaN in arrays, None on
  Point; stretch_param per chain; merge-branch tests (a)-(e) on both trees; names kept; small
  pull-ins; MEMORY_GROUPS folded into eight slots." D1-D6 / Q1-Q9 of the campaign plan bound the run,
  except Q9's role-5 assignment (overridden by the author's Q6 above).
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): the order along a curve is a
  `PointStore.Chain` (int64 `ids`, `stretch_param`, `anchor`), rows located by two int32 store columns;
  every curve is a `(root, tail)` view re-resolved per read; `Point.forward` / `backward` /
  `insert_point_*`, the walkers, the through-walk toggle and `_insert_point_geometrically` are gone;
  merge = one `searchsorted` (different chains, no shared cdist) or the old zipper on positions;
  refine = one `chain.insert` per layer. Probe identical on 62 dumps (auditor cross-compared three
  dump pairs, and N3's pre-N3 base dump against N4's branch: identical); suite 1248 passed / 43
  skipped / C3 alone (auditor rerun). Sequential profiler pairs at 1e-7: blast 113.1 / 112.2 -> 85.0 s,
  `walk` 19.2 -> 1.3, `refine` 8.2 -> 1.0, `infer_iterates` unchanged 10.3 (N5's), census 398.5 ->
  394.6 MiB, `Point` 100 -> 76 B; peak RSS over baseline -70 MiB (coder, branch on a dirty tree) vs
  -30 MiB (TE, clean commit). Behaviour changes are API-level plus states no build reaches: (f)/(g)/(h)
  raise; refinement at the anchor reads NaN on its left (library path: base TypeError -> completes;
  explicit `branch_index`: 9 vs 7 / 25 vs 21 points); anchor-bridge images merge from the slot.
- **Carried forward to the merge (text and one-liners):** the Outcome row should quote the TE's
  clean-commit profiler pair (or both); "share no cdist" -> "lie on different chains and share no
  cdist" (row and `ManifoldMachine.py:525-526`); `BaseManifold.py:248-249` (and `:194`) "None reads to
  the end of the root's chain" is false (None = the tail; measured 176 vs 232). Optional (D5 fence):
  `BranchPoint.py:5-12` (`forward`/`backward` pair, `only_forward` Dev Note). **To later slices:** F3
  (getter test's oracle shares `_span`; compare against `chain.ids`, N6); stale "linked list" text in
  `design_note.md:82`, `BasePoint`, `gpu.py`, `Trellis`, `BranchPoint`, `tests/helpers/invariants.py`
  (N6-N9); `_image_chains` reads `store._chain` / `_position` privately (no public view);
  `Chain.version` unread; plan doc §H / N8 still list the through-walk for N8 (already gone). N5
  request drafted in the outcome §10 (BranchPoint ends of segments without a row is its hard question).
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-25 — N5a `brackets-on-rows`, outcome (first half of N5 `segments-and-brackets-on-ids`; plan approved round 1 of planning with every default Q1-Q10, Q7 splitting N5 into N5a / N5b; review REVISE round 1 for two Outcome-row sentences and one test comment, PASS round 2; branch `refactor/brackets-on-rows` @ `9b72985`, base `01e991a`, eight commits)

- **Said** (plan gate, `artifacts/plan_ok_planner_human.json`): "Approved as written; all defaults
  Q1-Q10 accepted: split N5, this run is N5a brackets-on-rows (N5b segments-and-edge-keys-on-rows
  next); a periodic-point bracket end keeps the BranchPoint object (no -1 sentinel);
  Intersection.store; readers outside numerics move now via segment_ends, no shim; _polyline_gap on
  row equality; populate_intersections left (measured 0.23 s); public column views once per batch;
  batched inference search in scope; the cut converts rows to nodes once; Outcome rows 5a/5b."
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): `Intersection.unstable_segment`
  / `stable_segment` are int store rows (the `BranchPoint` at a periodic-point end), with `store`,
  `segment_ends(stability)` and the `brackets` memo; the iterate search bisects the same key list and
  masks one `(m, 4)` bracket array per bucket (`stable_brackets`; `candidates_near_stable_cdist` gone);
  `_polyline_gap` on row equality; the cut turns rows into nodes once. Probe identical on 62 dumps
  (auditor: three dump pairs; id 503 of k32_series .b09/.b10 is the only 0.0 bracket end, equal on both
  trees); suite 1250 passed / 43 skipped / C3 alone (auditor rerun); two ids added (the trap-(a) test
  and the TE's stable-bracket test). Five sequential profiler pairs, ALL under a concurrent user game:
  blast -9.1 to -10.6 s, `infer_iterates` 10.5-11.7 -> 1.10-1.34 s, `resolved_ids` -0.62, census +0.83
  MiB (registry), `segment_index` unchanged, peak RSS noise (base alone spans 1805-1883 MiB). Auditor:
  the probe is identical under m4 (stable comparisons dropped) and m8 (candidates reversed) -- the
  stable half of the bracket rule decides nothing on the corpus; only the two new unit tests pin it.
  No library output changed; API: rows in brackets, new `store` / `segment_ends` / `brackets` /
  `stable_brackets`, two new one-store assertions.
- **Carried forward to the merge:** nothing required. Optional text: `Tangle.py:274` "bracketing
  points" -> rows; `IterateInference.py:241` elliptical window comment. **To later slices:** N5b
  (`_Segment` rows, row edge keys, batched query-only add, the cut's `_query_edges.update`; trap (b):
  sid counter, one sid per NEW query edge paired or not, `frozenset((new_sid, cand_id))` order,
  `tuple(seg_id_pair)` operand order); the `brackets` memo assumes no cdist write after detection
  (nothing asserts it; N8 must check); peak RSS should not be bounded at 1 % from one pair; N4's
  carry-forwards (F3, stale "linked list" text, private `_chain`/`_position` reads, `Chain.version`,
  §H's through-walk line) still open. N5b request drafted in the outcome §10.
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-26 — N5b `segments-and-edge-keys-on-rows`, outcome (second half of N5; plan approved round 1 of planning with every default Q1-Q8; review REVISE round 1 for one Outcome-row sentence (peak RSS), PASS round 2; branch `refactor/segments-and-edge-keys-on-rows` @ `c944cf6`, base `1c25afb`, five commits)

- **Said** (plan gate, `artifacts/plan_ok_planner_human.json`): "Author approved the N5b plan with all
  Q1-Q8 defaults (2026-10-26)." -- i.e. `_Segment` fields renamed `row0` / `row1` plus `store`; one
  row key `((lo + 1) << 32) | (hi + 1)` for `_edge_to_sid` and `_query_edges`, `_edge` keeps its name,
  `_edge_keys` asserts rows < 2**31; sids by `list(islice(Tangle._ids, n))` from the same counter;
  pairs built per edge exactly as before (`frozenset((new_sid, cand_id))`, then `|=`); the cut's
  re-keying vectorised and kept; one public view per accessor call on the resolve path and ONE
  `store.xy` view per batched add; one run; the cdist type pinned by a unit test, not the probe.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): `_Segment(id, manifold,
  store, row0, row1, in_rtree)` with `coords(xy=None)` / `cdists(stability)`; int row keys in both
  maps; `Tangle._store` via `_adopt_store` (one store asserted at every add; `clear_all` resets it);
  the query-only add is `_query_curve` (one `intersection_v` per curve, a record only for a pair,
  every new edge's sid returned); `tuple(seg_id_pair)` and the operand order untouched. Probe
  identical on 62 dumps (auditor: own branch dumps of 8 builds / 14 dumps against the coder's base
  dumps, identical); shadow trace equal on 19 targets; suite 1252 passed / 43 skipped / C3 alone
  (auditor rerun); two ids added (the trap-(b) pins, both passing on `1c25afb`, auditor-verified).
  Two sequential profiler pairs under a user game: blast 74.9 -> 65.2 s (-13 %), `index` 9.0 -> 0.63,
  `cut` 0.63 -> 0.08, census 395.2 -> 383.9 MiB (`segment_index` the only component that moves),
  peak RSS noise. Auditor mutants: dropping the cut's re-keying fails the separator test; an unsorted
  `_edge` fails `test_readding_a_bridge_is_idempotent`; dropping the BranchPoint mask in
  `_query_curve` (m9) fails NOTHING and leaves the dumps identical (the path is reached once in the
  suite, deciding nothing). No library output changed.
- **Carried forward to the merge:** nothing required. Optional text: row 5b "-12.8" -> "-12.9";
  `Tangle.py:355-357` capture comment says "nodes" / "never off the index" (rows now). **To later
  slices:** the BranchPoint root edge of `_query_curve` unpinned (m9); scalar `_edge` has no 31-bit
  assert; the one-store assert unpinned (m6); `loom/ResonanceZone.py:27-30` `segment.p0_seg1` prose
  (N6); `_insert_segment`'s dead `else seg.id`; `diagnose` is 82 % of blast (53 of 65 s) and no slice
  of N6-N10 targets it; float bits still hang on `tuple(frozenset(pair))` (making the operand order
  explicit would be a D3 exception); the periodic point now has four spellings (object, -1 key
  offset, `np.maximum(e0, 0)` + mask, `Chain.anchor`) for N8 to retire; N4/N5a carry-forwards as
  listed in the outcome §9. N6 request drafted in the outcome §10.
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-26 — N6 `readers-on-views`, outcome (slice 6 of the numerics-arrays campaign; plan approved round 1 of planning with every default Q1-Q8; review REVISE round 1 for one false Outcome-row clause plus two wording fixes, PASS round 2; branch `refactor/readers-on-views` @ `e25adf1`, base `9dc702b`, nine commits)

- **The author said** (`artifacts/plan_ok_planner_human.json`): "Author approved the N6 plan with all
  Q1-Q8 defaults, including Q7 (accept the peak-RSS rise as retained free heap, judged by the two
  allocator controls) (2026-10-26)." -- i.e. the bridge walk memo and `Chain.version` deleted outright;
  `_segment_rows` keyed by row + per-store offset and one negative number per periodic point, the float
  reads moved too, two key pins; `grow_until_turnaround` only its first-segment read moved (§H already
  fixed by N3/N4); `_image_chains` keeps its private `_chain` / `_position` reads; the `BranchPoint`
  root edge of `_query_curve` pinned here; one run; the default-allocator peak-RSS rise accepted as
  retained free heap if both controls hold; `StablePartition.py:835` moved here.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): `get_point_array` returns a
  fresh writable array per call; new `BaseManifold.points` / `cdists` (root to tail); the same-stability
  field, `arc_polyline` (`stability=` unread), the anchors, the turnaround, `trim_stable_at` and five
  topology readers read rows / arrays. Probe identical on 62 dumps (auditor: own branch dumps of 8 builds
  / 14 dumps, identical); reader shadow 19 targets 0 differences; suite 1255 / 43 / C3 alone (auditor
  rerun); three ids added, all passing on the base (auditor-verified). Blast 63.7 -> 53.2 s (-16 %),
  `diagnose` 52.3 -> 42.4 s, census 383.9 -> 350.7 MiB. Peak RSS over baseline +156..+217 MiB and held
  RSS +360..+400 MiB under glibc's default allocator; both controls hold (auditor reran both: trim
  313 -> 302 / 474 -> 440 MiB, pinned mmap 1808.6 -> 1791.5). Unlisted in the plan's §5: the trellis's
  cached bridge polyline is now a writable shared array, and `arc_polyline` / `cdists` raise or return
  an object array on a periodic point not built by `FixedPointSolver` (no library path).
- **Carried forward to the merge:** nothing required. Optional text in row 6 (the reviewer's comma; the
  `diagnose` parenthetical; "the branch holds 311 -> 301"). **To later slices:** stale memo prose outside
  the file list (`CLAUDE.md`, `TangleSession.py:743-745`, `Trellis.py:578-579`, `test_geometry.py:102`,
  two test names); the profiler's `walk_memo` field; the writable cached polyline; mC / mH, the
  trim-root and arc-reversal mutants unpinned; `diagnose` 80 % of blast and owned by no slice; N7
  request drafted in the outcome §10.
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-26 — N7a `engine-on-rows`, outcome (first half of N7 `delete-point`, slice 7a of the numerics-arrays campaign; plan approved with every default Q1-Q10; review REVISE round 1 (R1: `BaseManifold.__init__` set the ends before `branch_index`, so a periodic-point root with a row tail and no `store=` raised; R2: `PointStore.unlink` was a dead write and two prose phrases credited it), PASS round 2; branch `refactor/engine-on-rows` @ `88f1a47`, base `d279815`, eleven commits)

- **The author said** (`artifacts/plan_ok_planner_human.json`): "Author approved the N7 plan with all
  Q1-Q10 defaults (2026-10-26): split by risk, this run is N7a engine-on-rows; N7b delete-point follows as
  its own request; Q2-Q5 and Q9 decided now, applied in N7b; Q6 and Q10 applied in N7a." -- i.e. Q1 two
  runs cut by risk; Q2 ends stored as rows (or the `BranchPoint`) plus `_store`, the API record
  `End(id, x, y, cdist)` exported from `tanglepack.numerics` only, `Bridge.head` a read-only alias of
  `root`, `store=` required for a row root; Q3 `BranchPoint` inlines `BasePoint`'s public members
  verbatim, no private base, both cdist fields kept, each reader on its current field; Q4 delete
  `return_nodes`, `_walk_nodes`, `first_node` (no row-returning variant); Q5 drop
  `arc_polyline(stability=)` at its six tracked callers with a pin of the storage-order reversal; Q6
  the read-only guard in `Trellis.cached_bridge_polyline`; Q7 N7a fixes only its own false sentences, the
  rest of the memo / linked-list prose, the rst and `design_note.md:82-85` go in N7b's text commit, test
  names kept, mC / mH pinned in N8, `CLAUDE.md` is the author's (replacement text in the plan's Q7); Q8
  the census on rows in N7a, `BasePoint` / `store.nodes` out and `MEMORY_GROUPS` relabelled in N7b; Q9
  `loose_store` deleted in N7b, one module-level test store in `tests/helpers/polylines.py`; Q10
  `ManifoldView.root` / `.tail` deleted.
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): the engine (`ManifoldMachine`,
  `ManifoldInitializer`, `Tangle`, `trim_stable_at`) appends, links and reads rows and makes no `Point`;
  curves store `_root` / `_tail` as rows of `store` (or the `BranchPoint`); `root` / `tail` still return
  N3's lazy adapter `Point`; a `BranchPoint` slot holds its `Chain`; new `ids`, `store`, `Chain.xy_at`,
  `store=`; the cached bridge polyline is read-only again. Probe identical on 62 dumps; whole store
  byte-equal on 7 cases (auditor: own dump of 10 MORE cases, store columns, chains, anchors, every
  manifold / bridge end row, version, slots and generation all equal); suite 1258 / 43 / 5 / C3 alone
  (auditor rerun); +3 / -0 ids. Census 350.5 -> 249.3 MiB (`points_image` 101.06 -> 0.24), blast flat
  (-1 %); malloc_trim control 302 -> 236 / 439 -> 332 MiB (auditor rerun). src +142 lines net (the
  adapters stay until N7b).
- **Carried forward to the merge:** nothing required. **To N7b:** the plan's §6 list; R1 unpinned; the two
  F401s in `profile_numerics.py`; the census on a store-less tree reads `n_points` 0 silently; the
  `_image_chains` head index unpinnable without a hand-built state; `refine_manifold` /
  `_refine_layer` `branch_index` unused; `Point.py`'s module docstring is false since N7a (goes with the
  file); `CLAUDE.md:131` "its slot names the first point" is false since N7a (the author's file).
  N7b request drafted in the outcome §10.
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-27 — N7b `delete-point`, outcome (second half of N7, slice 7b of the numerics-arrays campaign; plan approved with every default Q1-Q10; review REVISE round 1 (three false `BranchPoint.py` docstrings: `has_iterate_link`'s tie-rule reader, `insert_*_iterate`'s reversed `Raises` and argument name, "lies on one chain"), PASS round 2; branch `refactor/delete-point` @ `647786b`, base `4c9c3cb`, six commits)

- **The author said** (`artifacts/plan_ok_planner_human.json`): "Approved with all defaults Q1-Q10: one
  run, tests-first commits; five pins; raise on bad ends; keep insert_point_* names to N8; drop unused
  refine params; census raises on old trees, F401s go; arc_polyline(stability=) dropped at all 10
  callers; Sphinx = rubric 16 command; hole_sides.py :74 only; tests read ends via public ids."
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): `Point` and `BasePoint` are
  deleted; `root` / `tail` / `Bridge.head` read `End(id, x, y, cdist)` (a `NamedTuple` built per read by
  `PointStore.end`, exported from `tanglepack.numerics` only) or the `BranchPoint`; `BranchPoint` has no
  base class (`BasePoint` inlined verbatim, both cdist fields); `_stored` raises `TypeError` /
  `ValueError`; `store.nodes`, `node`, `node_at`, `loose_store`, `_walk_nodes`, `first_node`,
  `return_nodes`, `_manifold_from_points`, `arc_polyline(stability=)` and the unused refine parameters
  are gone; `tanglepack.__all__` 38 -> 36. Probe identical on 62 dumps and whole store byte-equal on 7
  cases (auditor: own fresh dumps of both trees); suite 1253 / 43 / 5 / C3 alone (auditor rerun);
  1302 -> 1297 ids (-10 Q7 / +5 pins). Census 249.32 -> 231.85 MiB, blast flat (auditor pair 52.20 ->
  51.75 s). src -398 lines net.
- **Carried forward to the merge:** nothing required; optional: `numerics/geometry.py:34` now-unused
  `Stability` import (new F401). **To N8:** `_stored` accepts an `End` of another store or an
  out-of-range row silently (reads uninitialised capacity); `BranchPoint` prose ("a fixed point or a
  crossing", "doubly linked list", `only_forward` Dev Note) and its dead members (`edist`, `edists`,
  `*_stretch_params`, `has_iterate_link`); the cdist-field choice; the `-1` head placeholder in
  `_image_chains`; mC / mH pin; the census's `BasePoint` fallback. **The author's:** `CLAUDE.md:131`
  (replacement in the outcome §9); the three untracked scripts (`henon_k25_length_growth.py:100,113`,
  `henon_k25_period3_class_focus.py:99`, `henon_k25_period3_class_iterates.py:184`); Sphinx's 25
  `PointStore` warnings (duplicate `Chain` attributes, two targets for `End`), judged acceptable by the
  reviewer. N8 request drafted in the outcome §10.
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-27 — N8 `branch-point-storage`, outcome (slice 8 of the numerics-arrays campaign; plan approved with every default Q1-Q10; review REVISE round 1 (Outcome row 8: the Tests sentence counted 1296 ids, not the tip's 1299, and the mmap control listed 3 of 13 pairs), PASS round 2; branch `refactor/branch-point-storage` @ `91677b8`, base `d0d78c5`, eight commits)

- **The author said** (`artifacts/plan_ok_planner_human.json`): "Approved with all defaults Q1-Q10:
  BranchPoint becomes a view over its FixedPoint (no store row); slot table on FixedPoint with
  attach/chain; only get_cdist(stability); iterate cycle as a read-only orbit view; IndexError on
  out-of-store rows; census drops fallback and stops at Chain; mH pin as written; prose fixed here;
  zero-shadow on base; Outcome row 8 only."
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): `BranchPoint(fixed_point,
  orbit_index)` stores nothing of its own (coordinates from `fp.coordinates[o]`, cdist 0.0, slots in
  `FixedPoint._slots[(stability, o, branch)]` via `chain` / `attach`, image = next orbit point);
  `BranchPoint.py` 346 -> 133 lines; no store row is added (the request's literal "a store row" was
  declined by Q1). Gone: `set_x` / `set_y`, `.cdist`, `cdists`, `edist(s)`, the stretch lists, the
  stored iterate links and their seven methods, `has_iterate_link`, `_has_iterate`, `insert_point_*`
  (now `attach`), the census's `BasePoint` fallback. `_stored` raises `IndexError` outside the store.
  Probe identical on 62 dumps and whole store byte-equal on 7 cases (auditor: fresh dumps of both
  trees); suite 1255 / 43 / 5 / C3 alone (auditor rerun); ids 1297 -> 1299 (-5 / +7). Census
  231.85 MiB both (-313 B), blast flat (auditor pair 54.54 -> 54.55 s at load 4-5). src -203 lines net.
- **Rubric 14(a) (pinned-mmap peak, +25 MiB per pair, "a failure stops the run") failed literally on 5
  of 13 pairs and passed on the rest; the reviewer judged the bound below the sampler's noise
  (bimodal peaks on both trees, equal maxima, trim control equal) and the run continued.** The
  auditor's 14th pair reads 1714.0 -> 1685.9. Recorded for the author's grade, and for the form of the
  check in later slices (median or max of N >= 3 per tree).
- **Carried forward to the merge:** nothing required; optional: the `BaseManifold.py:20-21` wording
  ("the next orbit point" is the unstable direction only) and `BranchPoint.attach`'s unlisted
  `IndexError`. **To N9 / later:** `attach` accepts a misspelt stability silently (phantom slot);
  `attach_chain`'s `IndexError` unpinned; stale prose `IterateInference.py:37`, `gpu.py:27` (N9),
  `test_manifold_initializer.py:101-102`; N10 serialises `FixedPoint._slots`. **The author's:**
  `CLAUDE.md:131` (now also its `BranchPoint` sentence; replacement in the outcome §9) and the
  `FixedPoint` row of the key-classes table; the three untracked scripts (unchanged since N7b).
  N9 request drafted in the outcome §10.
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-27 — N9 `array-namespace`, outcome (slice 9 of the numerics-arrays campaign; plan approved with every default Q1-Q9; review REVISE round 1 (three false sentences: Outcome row 9's Tests count after the test engineer's added id, the design note's "every other stage is host-only", and `gpu.py`'s reasoned list holding the store and registry without a reason), PASS round 2; branch `refactor/array-namespace` @ `4dd587f`, base `a2a2c4b`, five commits)

- **The author said** (`artifacts/plan_ok_planner_human.json`): "Approved with all defaults Q1-Q9: the seam
  only (array_namespace + the two geometry kernels; everything else named host-only in gpu.py); namespace read
  off the inputs; new numerics/namespace.py not exported; kernels never transfer, mixed raises; enable_gpu
  unchanged; fake-cupy plus device-parity tests; gpu.py Measured result restated with 2026-10-27 numbers; memory
  rubric = median and max of >=3 interleaved mmap runs per tree, trim control decisive; Outcome row 9 only."
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): §C N9's five-stage list is NARROWED to one
  seam (the plan doc's §C N9 and D6 stay as written history; row 9 records the narrowing). New
  `numerics/namespace.py` (`array_namespace`, not exported, never imports CuPy); `geometry.curvature_area_batch` /
  `chord_deviation_batch` compute in their inputs' namespace; no library caller passes device arrays, so nothing
  the library computes changed. `enable_gpu` / `disable_gpu` unchanged (AST equal). `gpu.py`'s docstring rewritten
  ("What runs where", re-measured result). Probe identical on 62 dumps and whole store byte-equal on 7 cases
  (auditor: fresh dumps of both trees); kernels byte-equal on numpy with and without CuPy imported (auditor);
  suite 1261 / 43 / 5 / C3 alone (auditor rerun); ids 1299 -> 1305 (+6, none removed). Census 231.85 MiB both,
  blast flat. src +120 / -29 (mostly docstring); `namespace.py` 39 lines.
- **The memory rubric's new form (Q8) held** on both runs that measured it (coder, test engineer): median and
  maximum inside +25 MiB, trim control within 0.5 MiB. Proposed for N10's rubric.
- **Carried forward to the merge:** nothing required; optional: `gpu.py:34` name the bulk segment build
  (`Tangle.add_manifolds`) beside `_query_curve`; the design note's "the user's map runs on the device" -> "the
  user's batched map calls". **The author's:** `CLAUDE.md`'s GPU row and data-model paragraph (replacements in the
  outcome §9); the three untracked scripts (unchanged since N7b); `stress/numerics/harness.py:60` puts its own
  tree's `src` first (a copied script silently loads main's `tanglepack`). N10 request drafted in the outcome §10.
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-28 — N10a `save-load`, outcome (slice 10a of the numerics-arrays campaign; N10 split by Q1; plan approved with every default Q1-Q11; review REVISE round 1 (a loaded `FixedPoint.coordinates` came back a list, three false sentences on set order / pair order / row 10a's counts), PASS round 2; branch `refactor/save-load` @ `90aba08`, base `ebc9806`, eight commits)

- **The author said** (`artifacts/plan_ok_planner_human.json`): "Approved with all defaults Q1-Q11: split, this run
  N10a save-load (N10b after-profile next); per-index sid counter as commit 1; session state saved beside the
  workbench; segment index saved exactly; PointBudgetExceeded state saved as-is with a Note; numerics/persist.py +
  loom/persist.py; npz+JSON format 1; round-trip tests on 11 law cases + continuation + one subprocess; attach_chain
  stability check pinned; N10b both cutoffs; Outcome rows 10a/10b."
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): `Tangle._ids` (process-wide) is replaced by a
  per-index `Tangle._next_sid`; a second build in one process is now byte-equal to a fresh one (no fresh-process
  output moved: probe identical 62 dumps, store byte-equal on 7 cases). `TangleWorkbench.save` / `load` and
  `TangleSession.save` / `load` (one npz + one JSON, format 1, `allow_pickle=False`, no float in JSON); a loaded
  session re-saves byte-equal on the 11 law cases in a fresh process and continues byte-identically (59 series
  steps, growth, blasts 14-16 at scale). Sets are saved sorted (a deviation, accepted). `attach_chain` refuses a
  misspelt stability. Suite 1380 / 43 / 5 / C3 alone; ids 1305 -> 1424. Census 231.85 MiB, blast within 0.4 %.
- **Rule from now on:** once merged, any change to what a save holds bumps `FORMAT_VERSION` (reviewer's rule).
- **Carried forward:** `CLAUDE.md:174` ("`Tangle._ids` is a process-wide counter"), the key-classes rows for save /
  load, the data-model paragraph and GPU row (the author's); `docs/open_issues.md:579` (WP8); `FixedPoint.py:81`;
  `stress/numerics/harness.py:60`; the three untracked scripts. N10b request drafted in the outcome §10.
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-28 — N10b `after-profile`, outcome (the LAST slice of the numerics-arrays campaign; plan approved with every default Q1-Q6; review PASS round 1; branch `refactor/after-profile` @ `9aa741e`, base `08b1d12`, three commits)

- **The author said** (`artifacts/plan_ok_planner_human.json`): "Approved with all defaults Q1-Q6: branch
  refactor/after-profile; skip a case whose stage calls differ with a note; rerun the traced folder; README past tense
  + row 10b + Closed line, CLAUDE.md drafted only; only stage_peak_rss.png changes; fix the README probe numbers."
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): the profiler ran unchanged at `dc1bd0c`
  (`src` / `tests` = `08b1d12`'s trees) into `figures/numerics_profile/dc1bd0c{,_scale1e-9,_scale_traced}/`, drawn
  with `--compare` against the `a8bb8a0*` folders (read, never re-measured). `figure_stage_peaks(runs, base, out)`
  draws each law case's baseline peak dashed and writes nothing for a scale-only folder; every other png, and every
  png without `--compare`, byte-equal. The campaign, `a8bb8a0` -> `dc1bd0c`: census 1914.0 -> 231.9 MiB, peak
  3034.7 -> 2015.0 MiB, blast 458.4 -> 52.2 s at 1e-7; 4683.2 -> 502.3 MiB, 6632.6 -> 1930.8 MiB, 541.4 -> 67.7 s
  at 1e-9; nothing computed moved. Suite 1380 / 43 / 5 / C3 alone; probe identical 62 / 6024 / 744.
- **The campaign plan is closed** ("Closed 2026-10-28 by N10b" under the plan doc's title). No further request drafted.
- **Carried forward (the author's):** every open item of N0-N10 is gathered in the outcome §9: the `CLAUDE.md` lines
  (`:131`, `:140`, `:174`, the key-classes rows, the drafted profiler line after `:90`), the Merge column of Outcome
  rows 1-10b, the three untracked scripts, `docs/open_issues.md:579`, `stress/numerics/harness.py:60`, the D6
  wording, a load-time map check and a `vars()` law; next-campaign candidates `diagnose`'s narrow phase (81 % of
  blast, and the peak's transient) and the segment index (48 % of the census).
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-28 — D1 `diagnose-on-new-segments`, outcome (the ONE run of the diagnose campaign, by Q1; plan approved with every default Q1-Q7; review PASS round 1; branch `refactor/diagnose-on-new-segments` @ `7d9d009`, base `88c62ea`, five commits)

- **The author said** (`artifacts/plan_ok_planner_human.json`): "Approved with all defaults Q1-Q7: one run; v2 strip
  sweep in numerics/diagnostics.py; oracle test through public functions; unresolved half out of scope;
  diagnose=True default kept; new profile folders + README paragraph; _query_curve intersection_v crash carried
  forward, not in this campaign."
- **Consequence, as audited** (`artifacts/outcome_auditor_human_2.md`): the same-stability sweep runs inside strips
  across the sweep axis, only where a new segment lies (`diagnostics.py` 481 -> 537 lines, stateless, nothing else
  in `src`). Nothing reported moved: probe identical 62 / 6024 / 744 (re-run by the auditor), every
  `SameStabilityCount` field equal per blast at 1e-7 (32) and 1e-9 (29), listings equal in order and every float,
  warning text byte-identical. `diagnose` 42.0 -> 22.4 s and blast 52.2 -> 32.5 s at 1e-7; 52.6 -> 30.0 s and
  67.7 -> 45.2 s at 1e-9. Peak RSS falls ~700-960 MiB at 1e-7 but RISES 60-95 MiB (about 4 %) at 1e-9, reproduced.
  `same_stability_crossings`' Returns corrected to "sweep order" (it said curve order: false); listing order, pair
  orientation and the axis-0 tie are now pinned by tests. Suite 1401 / 43 / 5 / C3 alone.
- **The campaign closes with this run** (Q1), one run of the four allowed.
- **Carried forward (the author's):** the `Tangle._query_curve` `intersection_v` latent crash (Q7; a small
  `query-curve-guard` run if wanted); the remaining whole-field ~10 s / ~20 s per series (persistent cache rejected,
  Q2); `separation` as the next blast cost; the 1e-9 peak rise.
- **Open for the author:** the grade; the merge (the author's). Not merged at writing.

## 2026-10-29 — Run 1 `cpu-workers` of the speedup campaign (/feature; branch `feature/cpu-workers`, base `5ed9940`)

- **The author's rulings (the approved plan):** `workers` is OPT-IN, default 1 (nothing runs on another core
  unless asked); `map_executor="process"` uses the FORK start method (the initializer installs the maps in the
  child, so a closure is never pickled; `ValueError` where fork is missing); the "lean blast" was rejected
  (the blast keeps every step it has; only the cores change). `workers` is a runtime resource: `save` / `load`
  neither store nor read it (format 1 unchanged).
- **Consequence:** for an elementwise map nothing computed moves at any worker count (identity probe identical
  62 / 6024 / 744 base vs branch at 1 worker and branch at 1 vs 8 workers with every batch split).
  Deviations from the plan, all measured: two split thresholds (`MIN_PARALLEL_POINTS` 32768 for a vectorised
  map in threads, `MIN_PARALLEL_POINTS_EXPENSIVE` 16 for a per-point map or a process pool); the separation
  guard is one logarithmic cloud of KD-trees queried with `distance_upper_bound` (the plan's per-sibling trees
  measured slower, 8.5 s against 4.1 s), and its queries take no workers (slower at 16); the Duffing example
  replaces the campaign's pendulum.
- **Measured** (`bench_workers.py`, E-cores, median of 3): the Duffing `solve_ivp` growth 3.58 -> 0.98 s at 8
  workers (3.65x); the leapfrog map does not speed up (too cheap); the 1e-7 scale blast 56.9 -> 37.6 s at 8
  (`diagnose` 1.76x), the 1e-9 blast 72.5 -> 58.3 s (1.24x); the sweep's serial whole-field part (19.3 of
  44.3 s at 1e-7, 39.0 of 56.3 s at 1e-9) is the ceiling. Suite 1421 / 43 / 5 / C3 alone.
- **Open for the author:** the merge; whether the sweep's whole-field prefix is the next target.

## 2026-10-29 — Run 2 `gpu-device` of the speedup campaign (/perf; branch `perf/gpu-device`, base `5dec003`)

- **The author's rulings:** keep `device-sweep` (the whole sweep on the device: edge dedup, sweep sort, strips,
  candidate pairs and a fused pre-test with loosened margins; the host re-tests the survivors exactly) and
  `device-map` (Duffing's leapfrog as one `--fmad=false` ElementwiseKernel in numpy's operation order,
  bit-identical; `gpu.py` documents how to write one); `device-pairs` dropped. The sweep must be usable with the
  map left on the host, so `enable_gpu(..., route_map=False)` exists: a Hénon run gets the sweep speedup with a
  byte-equal state (a routed Hénon map's CuPy last bits change the state, base behaviour, and by 32 blasts the
  same-stability counts too).
- **Measured** (`bench_gpu.py`, cores 16-31, RTX 5060 Ti, median of 3): the 1e-7 scale series, 32 blasts,
  CPU 50.1 s / `diagnose` 38.7 s -> `route_map=False` 17.4 s / 5.9 s (2.9x), digest and counts equal to the CPU;
  routed 18.0 s. Duffing growth 200 steps 0.27 s CPU / 2.33 s base GPU / 0.35 s fused; 2000 steps 1.93 / 19.6 /
  0.56 s, the state byte-equal CPU vs GPU. `min_batch_points` default 64 -> 16 (a fused map breaks even at ~10).
- **Open for the author:** `gpu_enabled` is not re-exported from `tanglepack`; `device-pairs`' 1.6x is recorded
  only here.

## 2026-10-29 — Run 3 `sweep-runner` of the speedup campaign (/feature; branch `feature/sweep-runner`, base `3a35d74`)

- **The author's rulings (plan):** the sweep runner is a script, `scripts/sweep_pool.py` (`Job`, `run_pool`), so the
  library still has no scheduler (D6); `processes x workers <= cores` or `ValueError`; resume is per finished job
  (the zone sweep: per finished `k`), never a mid-job workbench save; a job ended by a memory or time cap is not
  relaunched unless `retry`. The zone sweep's `--processes` replaces `--chunk` / `--merge`; the in-process path is
  unchanged (its `zones.csv` byte-equal to the base's). The author's nine k25 entropy scripts are tracked (commit 1
  verbatim, SHA-256 equal to the untracked copies on main at the cut), and only `henon_k25_entropy.py` (`workers`,
  the summary CSV) and `henon_k25_entropy_post_array_refactor.py` (`pool`) are edited.
- **Deviations, measured or forced:** a k25 sweep that reaches `--max-blasts` records no `stopped` reason, so a
  job is done when its JSON records a stop OR a row at its last blast; the zone sweep's seed `k` stays in the
  first chunk of the upward run (its row is written by that job), and each run of `sweep_order` gets its own chunks
  (`--processes 1` is one job per side of the seed); `--max-rss-gb` / `--max-seconds` of `pool` are the sweep's
  own soft stop, and the pool's hard cap is the same memory but TWICE the time (`HARD_TIME_FACTOR`: the sweep reads
  its clock only between blasts, so an equal cap always killed the job before the soft time stop fired; the memory
  stop is predictive, so an equal memory cap is a true backstop); a zone-sweep `k` whose build fails leaves
  `rows/<k>.skipped` (its error) and counts as finished, so a deterministic failure no longer relaunches its job on
  every rerun (delete the file to rebuild it).

## 2026-10-29 — Run 1 `lean-blast` of the deeper-blasts campaign (/feature; branch `feature/lean-blast`, base `8ce14b6`)

- **The author's rulings (plan):** a lean blast releases at the END of each iteration, never inside one (nothing
  can be freed before that iteration's topology is known); lean takes `min_separation` None or `"width"` and never
  diagnoses; a lean workbench refuses recomputes, cuts, growth, trims, zones, save, the same-stability readings,
  the dense arrangement and a blast that is not lean (`CurveReleased`); the default path is unchanged (identity
  probe: 62 builds identical).
- **Deviations, measured or forced:** the plan's keep rule let a later topology read a freed bridge (k32_series:
  a hole bridge released steps earlier), so `release_curves` also keeps the LINEAGE of the kept frontier, its
  ancestors (preimage rows) and every bridge their images tile (ancestors alone fail too: a mapped hole bridge's
  image). Unguarded, nothing dies, so lean frees about 1% or less (scale series 0.5% at 1e-7); with the width
  guard the store ends 3-5x smaller. The straddle rule keeps a frontier's existing images, the run of their chain
  and its phantoms, finer than "the neighbour's whole image chain". The width contracts by `per_step_beta("stable")`
  (< 1) per step, not `w / beta_s` (that would grow it), of the saddle whose stable branch holds the seed ring
  (`Bridge.width_saddle`; review round 2: a heteroclinic lobe of `W^u(A)` on `W^s(B)` read A's factor, NUM-20). `WIDTH_SEPARATION_FACTOR` stays 1.0 (0.01 matches the
  1e-5 guard's onset on the 1e-7 scale series but drops nothing in the k28 test builds). Author at merge (2026-10-29): the width guard stays OFF by default (`min_separation=None`) and the factor uncalibrated; whoever opts in sets `Blast.WIDTH_SEPARATION_FACTOR`. Widths
  are saved as an optional column (format 1 unchanged for a save without one). No series has a joined arc on a
  released chain: the joined-pieces rule is unexercised by any build, so a test hands a trellis a joined arc of
  two blast children a lean blast would release, and the keep rule's clauses are asserted directly
  (`_lean_keep` names every reader); the all-fixed-points minimal trellis is never released without its clause
  on any suite build, so only the direct assertion depends on it. A lean blast's released drops stay out: `min_separation=None` refuses to re-admit
  one (`CurveReleased`) instead of keeping every drop's rows.

## 2026-10-29 — Run 2 `anchor-depth` of the deeper-blasts campaign (/feature; branch `feature/anchor-depth`, base `d4c068d`)

- **The author's rulings (plan):** clause 3 of N, the near-anchor guard, at `c = 5 = 1/(1 - 0.8)` (the trust band) and ON by
  default (`anchor_resolution_factor`, None off); the graded mesh near the orbit is opt-in (`near_orbit_radius`, None
  default); trusted depth is a flag, never a decision.
- **Deviations, measured or forced:** (1) the stable-neighbour half of clause 3 is NOT in N: at c = 5 it took real,
  well-resolved crossings (ratio 4000-6000) out of `heteroclinic` (three of one fold) and, on the k28 / k32 / nested series,
  a crossing of the original trellis whose stable neighbour is a deep blast crossing of another curve; an arc then loses an
  odd number of crossings and the arrangement, the heteroclinic case sanity and the hole-sign assertion fail (76 tests).
  Restricting it to one unstable curve still broke `heteroclinic`. Author at merge (2026-10-30): the pair rule stays out. (2) Clause 3 floors
  `h` at the xy resolution `64 eps max(1, |P|)` (`geometry.xy_resolution`): with the mesh `h` falls below what float64
  resolves and the orbit was registered to 1 ulp of P (depth 23, trusted 21 at 1e-7). (3) The pip-end zone rule takes only
  the bridge ARRIVING at `f^j(q0)`, `j >= 1` (second crossing on the orbit of the zone's own corner pip, `boundary_intersection`, never the session's chosen one): the plan's "ends on the orbit"
  iterated the turnstile's exit lobe leaving `q0` and changed every blasted law case; the arriving arc is inside by
  orientation (design note A.10). (4) The mesh reads `r` as the CHORD's distance to the orbit point (a deep image curve
  passes P between nodes 0.04 apart), and grades the seed once at initialisation (`ManifoldMachine.grade_seed`: the first
  fundamental domain carried inward on the seed's ray, linked as iterates; mapping a stable seed forward amplified its
  error 31x per step). (5) Five regression builds that reproduce pre-guard defects at deep or coarse blasts run with the
  guard off, as they already run with clause 1 off.
- **Measured (k=2.5 scale series, 1e-5 separation, to the plateau, one core each):** pip / trusted depth base 9 / 7 and
  12 / 9 (1e-7 / 1e-9); guard 7 / 7 and 9 / 9; guard + pip end the same (the guard binds before the zone filter); guard +
  pip end + mesh (r0 = 1e-3) 18 / 18 at both cutoffs, ended by the xy-resolution floor (d = 2.8e-13), blast 51 s / 103 s
  against base 76 s / 98 s, plateau points 2.10M / 5.33M against 2.18M / 5.10M. Identity probe: knobs off identical on all
  62 builds; defaults differ only on `k32_series` from blast 7 and `nested_outer_series` from blast 8 (shallower pip orbit).

## 2026-10-30 — Run 3 `cdist-precision` of the deeper-blasts campaign (/refactor; branch `refactor/cdist-precision`, base `da963fb`)

- **The author's choice (plan):** the residual column. The float `cdist` column stays, computed exactly as before;
  a float64 `cdist_lo` beside it holds what the rounding dropped, so the precise value is the double-double
  `cdist + cdist_lo` (`numerics/precise.py`, each kernel's `hi` the old float expression bit for bit). Plain cdists,
  merges and the public API (`cdist`, `get_cdist_array`) are unchanged; 8 bytes a row; save format 2, format 1
  loading with residual 0.0.
- **Readers left on the plain float, measured or forced:** (1) the merge, its run cuts and the chain check: a precise
  merge keeps the 17 plain ties the zipper drops on the k=2.8 one-blast build (cdist 1 to 36, the two points 4e-5 to
  9e-3 apart) and reorders neighbours one ULP apart (4 on k=2.8, 1 on the inversion builds), which changes the law
  cases; switching it is the author's decision. (2) The crossing sign: a precise reading flips one sign of the k=2.8
  series at blast 7, on a stable segment whose ends are one ULP apart. (3) `BRACKET_PAD` stays 1e-12: 1e-24 left the
  identity probe unchanged but linked three crossings fewer on the k=2.8 deep-blast sliver build (422 -> 418
  crossings). (4) Partition ownership and the dual-graph walks: their intervals and nodes hold plain floats
  (`intervals_from_marks`, `PartitionFamily`, `DualGraph`, outside the run's files).
- **Measured:** identity probe identical on all 62 builds (the residuals dumped apart, `<name>.lo.npz`). Unguarded k=2.5
  scale series at 1e-7 (`near_orbit_radius=1e-3`, `min_separation=None`, stopped at 9.1M rows, blast 15): base and branch
  compute the same rows, crossings, links and depth (pip / trusted 16 / 16); neighbours tie as plain floats from blast 10
  (8 ties; 773 at blast 15, `dc/c` down to 1.1e-16) and the residuals tell every one apart, none against the chain; no
  merge drop, no crossing-order tie or flip and no ambiguous link on either. Cost: +17.3 MiB at the 1e-7 guarded
  plateau (2.10M rows, store 104 -> 121 MiB); blasts about 15% slower deep in the series (13.2 -> 15.2 s at blast 14),
  in the registry's stable buckets and the per-source precise meet test of `infer_iterates`.
- **Review round 1: the order key is lexicographic.** `Intersection.cdist_key` is `(cdist, residual)` compared as a
  tuple, not the normalised precise pair: the plain float decides and the residual only breaks a plain tie. A precise
  order could contradict the plain-merged chain (the one-ULP neighbours above), cutting a bridge backwards and failing
  the sparse `Arrangement`'s plain-order assertion; the key cannot. `Tangle.create_bridges` asserts its sorted
  brackets never step back along the chain. Dedup and the iterate brackets still compare the precise value. Left
  open: two crossings on ONE segment whose ends the residuals reverse, tied as plain floats, are still ordered by
  their residuals (against `t`); it needs the merge on precise values (the author's decision above).
- **At merge (2026-10-30), the author accepted** the reworded bridge invariant (a `BridgeId` runs in chain order,
  the unstable-cdist order up to the ULP a lerped cdist steps back by at a fold of tied nodes) and merged with the
  merge, the crossing sign, `BRACKET_PAD`, ownership and the dual walks left plain. A precise merge would change
  computed results: a `/feature` of its own if wanted.

## 2026-10-30 — Run 4 `class-representatives` of the deeper-blasts campaign (/feature; branch `feature/class-representatives`, base `b4a5582`)

- **The author's rulings (plan):** `blast_zone(representatives=...)` opt-in, policies `fewest_points` (default in the
  tests and script; ties by smallest id), `smallest_id`, `largest_id`; held bridges are ALWAYS re-offered, nothing is
  remembered on a bridge (no `numerics/Bridge.py` change); pip-orbit, minimal-mapped and unclassed candidates are always
  iterated; a bucket is the refined child, else the unmatched pair, else the class, a loop by class and element.
- **Deviations, forced:** (1) each step links the iterate table (`infer_iterates`) before it reads the topology, so a
  step inside one call sees what a fresh call would and `n` one-step calls equal one `n`-step call. (2) No regression
  file of its own: the D1 series tests of `tests/regression/test_class_count_under_blasting.py` run under both
  `representatives=None` and `"fewest_points"` (a fixture), Theorem T and the B7b base-chord law included, so each law
  is stated once; `cases.blast_series` and `cases._blast` take `representatives=` (tests-only, outside the plan's list)
  so the D1 records are built by the one loop. The plan's "the alphabet stays closed" is no test: `symbolic_dynamics`
  adds every class a word names to its table, so it holds by construction. (3) "The series laws" are those D1 laws: the full `LAYERS` checks fail on the
  FULL series too at depth (k28 blast 6+, k32 7+, nested_outer 8+: Euler, rows, same-stability, table completeness), so
  they are no law of a deep blast.
- **Measured (`scripts/blast_representatives.py`, `figures/blast_representatives/`):** identity probe identical on all
  62 builds with `representatives=None`. Equivalence against the full series on every `cases.SERIES`,
  `nested_inner_blasts` and `heteroclinic`, three policies: alphabet, closed alphabet, words, refined matrix up to
  relabelling, holes, references, `pip_depth`, `trusted_depth` equal at every step, Theorem T never failing on either,
  except (a) k32_series blasts 6-9 (`fewest_points`; 6 `smallest_id`; 6-10 `largest_id`): the full blast loses its
  reference `(9, 10)` to crossing 93, born on the image chain of bridge `(34, 35)`, which the representative blast holds
  (its `(24, 25)`), so the reference switch, and with it holes 12 -> 13 and the alphabet 12 -> 13, arrives when the
  held bridge is iterated; (b) from blast 4 on k28, inversion, k32 and nested_outer, 1-4 trajectory pseudoneighbor pairs
  of the full blast are missing (a thinner X), never a reference. Cost on the k=2.5 scale series at 1e-7, blast 12:
  848 -> 433 bridges, 1.11M -> 0.58M points (1e-5 guard), 1361 -> 529 and 1.80M -> 0.71M (unguarded): 1.9-2.6x, not the
  10x target. Cause: the pip depth stops at 7 (the near-anchor guard) by blast 6, and from blast 8 the frontier lies past
  the cut trellis, is unclassed and always iterated (92 candidates at blast 12: 87 unclassed, 5 bucketed). Lean +
  representatives unguarded, 1e-7: blast 18, 11708 bridges, 15.3M rows held, max unstable cdist 5.8e15, stopped by the
  6M-point per-call budget at blast 19 (169 s for blast 18), depth 7 / 7 throughout. Prune question (k=2.5 scale,
  representatives, guard on, 1e-9 / 1e-10 / 1e-11 / 1e-12, `min_separation` 1e-5 / "resolution" / None, capped by 4M
  points per call, 15M held and 40 min): every run reaches pip = trusted depth 9 / 10 / 11 / 12 (the guard's), and at that
  common trusted depth and at the last common blast (16-18) the three separations agree on every item (4 classes); the
  prune hides nothing at trusted depth there. Only the 1e-5 runs at 1e-9 / 1e-10 exhausted their frontier (blasts 34 /
  36; 1e-10 then has 15 classes, beyond what the other two reached within budget).
- **Open, the author's call:** the representatives do not pay while the pip depth is stalled; holding unclassed bridges
  past the cut (or the graded mesh, which lifts the depth to 18) is the lever, not part of this run.

## 2026-10-31 — D12 reversed (/feature; branch `feature/drop-unstable-pip-cut`, base `f75fde4`)

- **The author's rule:** use every computed manifold. A chosen pip cuts the STABLE branches only, branch `j` at
  `c_s(f^j q0)`, `j < k_value`: `T(q0) = (W^u computed, W^s[P, q0])`, and every computed unstable crossing whose
  stable end lies inside that cut is a crossing of N. The class count may fall under blasting; the entropy may not.
  `m` (`BlastResult.pip_depth`, `trusted_depth`) stays a report; the blast's pip-orbit `min_separation` exemption
  stays; no public API changes.
- **Evidence (the /investigate of 2026-10-31):** D12's unstable cut hid real crossings: the class drops it prevented
  are refutations of reference pairs by crossings well resolved (ratio 128-9650), 37-43% inside the pair, persistent at
  10x / 100x finer `area_cutoff`. It froze the topology whenever the pip orbit stalled: k28_series blast 8 gives 15
  classes uncut against 4 cut; k32_series blasts 8-10 give 34 / 34 / 30 classes uncut (13 cut) with h 0.540465 ->
  0.545622 -> 0.548464. Every law case but the three inversion ones is unchanged (crossings, holes, classes, words,
  matrix), and so is `p3_series`; the heteroclinic case's selection trellis now holds all four heteroclinic crossings.
- **Laws:** `series_class_count_never_decreases` is replaced by `series_entropy_never_decreases` (`ln rho` of the
  refined transition matrix, 0 for an empty or nilpotent one, per `SeriesStep.entropy`); the cut law states no
  crossing past a stable cut and every resolved crossing inside the stable cuts held.
- **New KNOWN_ISSUES, against the "may only shrink" note (INV-12, `wp="author: inversion investigation"`):** the
  unblasted `inversion` case keeps one reference `(23, 24)`, one hole, one class and no active class, so ten laws
  (`anchor_bridge_class_leads_its_tangle`, `registered_images_keep_markers`, `partition_singletons`,
  `backward_holes_keep_bridge_side`, `anchor_opening_faces_its_hole`, `propagated_holes_land_on_predicted_branch`,
  `propagation_terminates`, `propagated_hole_inside_container`, `refined_matrix_is_an_in_splitting`,
  `matrix_is_token_counts`) and `base_chords_are_inside_the_lobe_image` check nothing on it; `inversion_series`
  blast 4 has two ambiguous walks and its entropy falls 0.264 -> 0 (reference switch to `(56, 57)`), under the full
  and the representative blast alike, so `series_entropy_never_decreases` and `series_no_walk_is_ambiguous` are strict
  xfails on it without a policy-specific id. `members_share_the_P1_itinerary` x `inversion` is out of NOT_APPLICABLE:
  its one class `u` now has two members with readable image chains, so the law tests one class there and passes
  (the set-aside test would fail).

## 2026-10-31 — pip-depth machinery removed (/feature; branch `feature/drop-pip-depth-machinery`, base `c58f6f5`)

- **The author's rule:** with D12's unstable cut reversed the pip depth `m` limits nothing, so the three pieces built
  only to deepen it go: the blast's pip-orbit `min_separation` exemption (B3), `BlastResult.trusted_depth` (with
  `TRUST_BAND` / `TRUST_LINEAR_RADIUS`) and the opt-in graded mesh (`near_orbit_radius`, `ManifoldMachine.grade_seed`).
  Clause 3 of N (the near-anchor guard, its factor 5 = 1/(1 - 0.8) and its `xy_resolution` floor) stays, `pip_depth`
  stays as the one depth report, and the pip-end frontier rule (bridges arriving at `f^j(q0)`) stays. Refinement is
  the old `near_orbit_radius=None` path byte for byte.
- **Measured (each series to its depth, main with the exemption against the branch; last step shown):**

  | series | blasts | main m / N / bridges | branch m / N / bridges | holes | classes | h |
  |---|---|---|---|---|---|---|
  | k28_series | 8 | 8 / 239 / 168 | same | 24 | 15 | 0.528049 |
  | k32_series | 10 | 7 / 717 / 507 | same | 50 | 30 | 0.548464 |
  | p3_series | 9 | 23 / 67 / 49 | same | 42 | 23 | 0.101522 |
  | nested_outer_series | 11 | 8 / 493 / 349 | 7 / 486 / 344 | 12 | 9 | 0.528049 |
  | inversion_series | 4 | 9 / 47 / 33 | same | 5 | 6 | 0 |

  Holes, classes, `h` and `is_reliable` are equal at every step of all five; only nested_outer moves, from blast 7 on
  (m 8 -> 7, 4-7 fewer crossings of N and 3-5 fewer bridges). The exemption had kept near-coincident curves near `P`
  (361 same-stability WARNINGs on nested_outer; 37k-420k per blast at k=2.5, 1e-13).
- **Save format 3:** its settings hold no mesh radius. A format-1 or format-2 save loads when its stored radius is
  None (or absent) and raises `ValueError` naming the graded mesh otherwise; format 1 still loads its residuals as 0.0.
- **Tests:** the exemption and trusted-depth tests are gone; the width and cloud guards are tested on pip-orbit
  children like any other (`test_pip_bearing_child_is_dropped`, `test_the_pip_orbit_child_is_not_exempt`); the series
  law compares each blast's `pip_depth` with the depth composed test-side; the lean joined-arc test blasts k28 five
  times, not six (at six no joined pair is released any more).

## 2026-11-01 — Zone nesting is declared, never inferred (the author's rule; /feature; branch `feature/declared-zone-nesting`, base `7732025`)

- **Rule:** "There really shouldn't be any identification of nestedness in the actual code. The user who adds the different zones together should specify that information." `add_resonance_zones(nesting={inner_pip: outer_pip})` sets `ResonanceZone.inside`; zones not named are not nested, however they overlap. This amends C4 (DEC-11, NEST-4): an outer blast still leaves a declared inner zone's bridges alone. Save format 4 holds `inside`; formats 1-3 load with no nesting.
- **Why:** the hydrogen heteroclinic pair are mirror-image zones under `-I`, with areas equal up to rounding. `other.area < zone.area` in `blast_zone` called the one smaller by 1e-12 inner, so the other lost every bridge in their overlap: `{A: 0, B: 2}` or `{A: 2, B: 0}` by energy, rho 2.699 or 2.0 on the lower plateau where the transition graph gives 2.95231. All 23 energies from 0.27 to 0.325 followed the sign of the area difference.
- **Removed:** the area comparison in `blast_zone`, the smallest-area choice in `classify_bridge` (a declared inner zone before its host, then the bridge's own fixed point's zone, then the first stored), and the area and containment inference of `cartoon_zones.parents`. `ResonanceZone.area` stays as a measurement.

## 2026-11-10 — A class is read off its face ring (the author approved building; /feature; branch `feature/ring-classes`, base `0695321`)

- **Rule:** the pair of elements names one homotopy class only in a simply connected zone. In a disk face of the minimal arrangement an arc is determined by where its ends sit on the ring: a bridge or walk piece cutting exactly one kept bridge off the ring, with no hole between, is that bridge's class (`DualGraph.alongside`, `Trellis.parallel_bridges`, `DualWalk.walk_itinerary`), and a tie of a class a word named is decided by the product of the decided walks of the kept bridges on a hole-free side of its naming piece (`SymbolicDynamics._product_walk`; the word map is a homomorphism). Combinatorial throughout; silent off a disk face, where the pair stands as before.
- **Why:** hydrogen E = 0.5, B = 2.141 (Mitchell 2009, Fig. 8). The outer zone is an annulus around the inner zone and the inner fixed point is a cone point, so arcs pass it on its free side: the bridges alongside the inner boundary arc `g` join the same two elements as `b`, and the class `d2` had two tied walks, both words of real classes on one pair.
- **Replaces:** `diagnostics.wrapped_end` (the sign of a lobe polygon's area, a geometric reading that decided), the orbit-marker re-reading in `walk_itinerary`, and `shortest_walks(namer_walls=)` / `_crosses_namer` (commits `62c64c9`, `0695321`). `diagnostics.rounds_the_stable_end` and its row exemption stay.
- **Not built:** a class cutting SEVERAL kept bridges off its ring (a product) keyed by its ring position. Such a class is still its pair; when its members sit in several ring positions it is one WARNING (`diagnostics.classes_in_several_ring_positions`), as on the same map with the inner zone cut at an ordinary-sign pip, where the run is not reliable.

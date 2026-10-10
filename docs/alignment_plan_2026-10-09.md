# Aligning tanglepack to its mathematics — the plan (final draft, 2026-10-09)

Written against `docs/code_fidelity_audit_2026-10-09.md` (the audit), `docs/pipeline_chalkboard.md`,
the three algorithm notes in `docs/`, `CLAUDE.md`, `.refactor/decisions.md` (D1-D15, C1-C22) and the
working tree with the uncommitted D15 change. Every audit claim this plan builds on was re-read in the
code; where the plan reads the audit differently it says so in A.0.

The author's goal, in their words: *"strong simple fundamental algorithms, and to solve bugs by
clarifying those algorithms rather than piecemeal."* The plan has one shape: **state each stage as a
definition, then make the code that definition and nothing else.** Where a definition is missing it is
written here as a candidate, marked `DECISION`, with a recommended default and the consequence of each
option; section E collects those for the author. Slices (section B) are ordered so that decisions come
before the deletions that depend on them and every slice leaves the suite green or changes a law on
purpose.

Final draft (awaiting the author's answers to E): A and E complete (E gained Q11-Q13); the slices B0, B1, B4a, B4b, B6a, B8, B10 are
in house-planner form (B1 waits on Q13, the others are decision-free); B2, B3, B5, B6b, B7, B9, B11,
B12 wait on the author and stay as skeletons.
C and D cover the expanded slices in full and the rest as estimates. Probe results quoted below were
run by the orchestrator (`scratchpad/plan/probe_q123.py`, `probe_q2b.py`, `probe_p1to5.py`,
`probe_q3b.py`) or here (`probe_p2.py`).

---

## A. The target model

One definition per stage. "Code today" says how far the code is from it; "Change" names the slice.

### A.0 Where this plan reads the audit differently (agreed with the orchestrator)

- **Stage 5, "lookup first is a fiction".** Measured correctly (`from_table` False on every call on
  k10/k28/inversion/heteroclinic), but it is not a defect of the algorithm: a reference pair lies in the
  window `S \ f(S)`, whose preimage leaves `S`, so its backward images are NEVER registered crossings,
  by construction. The PDF's own mechanism is "simple division". The defect is the documentation
  ("table first") and the silently discarded flag at `StablePartition.py:513`. The fix is to state the
  rule (containment is decided in the canonical coordinate, exact up to each cdist's own noise), drop the
  flag, and add ONE WARNING where the container can be wrong: a scaled span end within `SCALING_RTOL`
  of a crossing of the container's branch (NUM-5's mechanism). That is B4b.
- **"D14 contradicts D11".** Narrower than stated. `.refactor/decisions.md` (WP9c) records that D14's
  premise did not hold on the fixture: `bn` has a computed ANCESTOR under the word map
  (`bn <- bl <- bi`, `bi` spanned). The consistent reading, accepted as the definition (Q9): *a class
  exists iff a computed bridge spans it or the word of an existing class names it* (closure from the
  spanned classes). That is what the code does; it is stated in A.8c.
- **"The marker rule is not a complete invariant".** True of the marker rule alone; with the walls
  (images of disjoint unstable arcs are disjoint) and the window rule the constraint set is the full set
  of facts the real `f(b)` satisfies. Completeness is still not a theorem, so a tie that survives them is
  honestly `ambiguous`. The plan keeps that status and removes only the tie-breaks that are not
  theorems (first-in-order after a cap).
- **k=10 "never checked end to end".** Confirmed, and the real build passes today
  (`a -> a u^-1 a^-1`, refined matrix `[[1,1],[1,1]]`, reliable, 3 holes, 8 crossings in T(q0)).
  The oracle law can be added first and passes on the base (B0) — but it is a case-specific pin, which
  the author's 2026-10-05 test policy forbids; it needs Q11.
- **"Dead or test-only API".** Four of the names (`euler_characteristic`, `Arrangement.preimage_of`,
  `StableNode.other_face`, `preimage_bridges`) are read by the law tier as ORACLES (`laws.py:1633,
  1697, 2214, 2250`, `test_bridge_identity.py:225`). They are test-only, not dead; B10 keeps them or
  moves them into `tests/helpers`, it does not delete them.

### A.1 Stage 0 — manifolds

**Definition.** For each branch, a polyline approximation of the manifold with canonical distance
`c` non-decreasing along it; the fundamental segment `[x, M^k x]` mapped with refinement so that every
chord deviates from the image curve by less than `area_cutoff` (chord-deviation rule `h = 1.5 A / L`).
`FixedPoint` is the one source of branch bookkeeping (`advance_key`, `per_step_beta`, `branch_cycle`,
`preserves_orientation`).

**Code today.** Faithful in shape. The data model (two linked lists in `Point`, `BranchPoint`, walk
caches) is a representation choice, not a definition; `edist` and the stretch slots are dead.

**Change.** None to the definition. The array representation is slice B11 (large, last, optional).

### A.2 Stage 1 — the crossing set N

**Definition.** `N` = the RESOLVED transversal crossings of an unstable polyline with a stable
polyline, each `(u_key, c_u, s_key, c_s, sign, brackets, resolution)`, found by ONE procedure on every
path (first compute, recompute, blast image): query each unstable segment against the stable index,
resolve each hit once, dedup by keys and cdists (`DEDUP_RTOL`), bracket on the current polyline.
Every detection is registered and reported; a detection is a crossing of `N` iff the polylines
RESOLVE it: its separation from its brackets exceeds the chord deviation of both polylines there
(`resolution_ratio ≥ threshold`, today 3, derived bound 1) AND no other detection of the same two
curves lies closer to it along the unstable polyline than that chord deviation (a near-coincident
RUN is one unresolved stretch, whatever its parity). An unresolved detection is not a point of `N` and
not of `X`: no bridge ends at it, no iterate link starts or ends at it, the arc across it is one arc,
and `report_resolution` names it with the remedy — refine (NUM-3/4) — exactly like a crossing past the
pip cut is excluded without being deleted (D10). The anchor is one synthetic crossing per periodic
point (D4).

`DECISION Q13 — a resolution criterion inside the definition of N.` Today every detection is fed to
the topology ("flags never discard", NUM-4 deferred) EXCEPT on the blast path, where
`Tangle._collapse_noise_crossings` (`Tangle.py:569-627`, called only from
`populate_intersections_for_manifold` ← `BridgeIterator.py:363`) keeps the median of an odd
near-coincident run and drops an even one. The orchestrator's probe shows the vote is load-bearing:
k28 at 8 blasts, collapse disabled at build — crossings 248 → 252, holes 24 → 40, classes 15 → 23,
closed 23 → 31, words diverging after `a -> a u^-1 b^-1`, `b -> c`; BOTH results report `reliable`.
The two runs are ODD runs of three, `(75, 1, 2)` and `(125, 250, 251)`: the build kept the median of
each (75, 125) and the recompute gains the two outer detections. Along the unstable polyline the
members sit 2.9e-7 .. 8.5e-5 apart at `c_u` 8e6 .. 5e7 (relative 1e-13 .. 2e-12), and every one of
them is well resolved by the PER-CROSSING ratio (1: 2297, 2: 4642, 250: 624, 251: 633; the kept medians
75 and 125 (12.0) are the LEAST resolved of their runs). So clause 1 of the definition holds none of
them out; clause 2 (a run closer than the chord deviation is one unresolved stretch) does all the work,
and "keep the best-resolved detection of a run" would keep the wrong one — a run is unresolved as a
whole. Treated as crossings the four spawn 16 holes and 8 classes, and 251 steals the `+1` link of
the real crossing 125 (the B8 collision). Two OTHER detections of that build, 148 and 249, are below
the ratio threshold today and fed to the topology (`unresolved_crossings(3)`); on every law case that
set is empty (minimum ratios: k10 522, k28 x1 5780, x2 2450, p3 / nested / nested_inner 8.7, inversion
124 / 19.5 / 12.2, period2 33.7, **heteroclinic 3.04**). Across the whole corpus (11 law cases, p3 to 9,
inversion to 4, k32 to 10, nested_outer to 11 blasts) the collapse fires ONCE, on k28 x8. Options:
- **(a, recommended default)** the definition above: an unresolved detection is registered, reported,
  and NOT in `N` until refinement resolves it, on every path. Recompute commutes; on k28 x8 six
  detections are held out (the two runs' outer members 1, 2, 250, 251 by clause 2, and 148, 249 by
  clause 1), so its topology is RE-DERIVED by the slice, not assumed to be the base's 24 / 15 / 23.
  The threshold enters the definition of `N` as the polyline's own accuracy: `resolution_ratio =
  separation / chord deviation`, so ratio 1 is the derived bound (the crossing's displacement equals
  the polygon error) and 3 is an undeclared margin that the heteroclinic minimum 3.04 sits 1 % above —
  with (a) the margin is load-bearing and must be chosen (Q13's second part), not tuned to silence.
- **(b)** refine locally around every unresolved detection until each grazing run vanishes or becomes
  resolved crossings (NUM-3/4). The right remedy, and (a)'s WARNING points to it; alone it costs growth
  on every blast and still needs (a)'s rule for the detections it has not yet refined.
- **(c)** keep a collapse but run it on every path. Consistent, still a vote by parity; not recommended.

**Change.** Slice B1 (after Q13).

### A.3 Stage 2 — the iterate table

**Definition.** A partial bijection `N → N`, `y = M(x)` iff `y` lies on the advanced branches and both
its brackets meet `x`'s carried one step, uniquely, and the map confirms it. A crossing with several
candidates, or an image claimed by two sources, is an UNRESOLVED LINK: recorded as data
(`registry.unresolved_links`), never linked, never resolved by registry order.

**Code today.** Bracket rule and map gate as stated. A second source for an image "logs and the first
in registry order wins" (`IterateInference.py:125-132`): order-dependent. A multi-match is logged at
INFO and forgotten (`:227`).

**Change.** Slice B8.

### A.4 Stage 3 — the strong pip and the trellis of the pip

**Definition.** `q0` is a strong pip iff no crossing of `q0`'s fixed point, mapped onto `q0`'s stable
branch by the minimal number of backward steps in the branch cycle (table first, else β-scaling),
lands strictly inside the open box `(0, c_s(q0)) × (0, c_u(q0))`; `q0`'s own orbit is recognised by id
(table) or by a 2-D collision (`collision_rtol`). Choosing `q0` cuts the trellis to
`T_m(q0) = (W^u[P, M^m q0], W^s[P, q0])`: stable branch `j` at `c_s(M^j q0)`, `j < k`; the unstable
branch of `M^i q0` at `c_u(M^i q0)` for the `k` deepest `i ≤ m`. **`m` is part of the definition.**

`DECISION N — what fixes m.` Today `m` = the deepest forward iterate of `q0` the table holds, and that
depth is set by blast policy: `min_separation` (library default `None`, but passed as 1e-5 by every
recipe: `cases.py:301` `build_k28`, `henon_cases.py:183`) drops a child "too close" and never feeds it
again (k32 blast 9: `f^10(q0)` dropped, the pip orbit stalls), and the zone-interior filter never
iterates a bridge outside the zone (the BCG-8 reference). Options:
- **(N-a, recommended default)** `m` stays "deepest registered", and the blast NEVER drops a child
  holding a pip-orbit point (`min_separation` exempts it; the zone filter is the definition of a blast
  and stays). Consequence: the pip orbit deepens with every blast as long as its bridge is in the zone;
  D12's "class count never drops" stays a measured law. ~20 lines in `Blast.py`.
- **(N-b)** `m` is an explicit parameter (`set_strong_pip(fp, id, depth=)`, default N-a's value), so two
  blasts can be compared at ONE `T_m`. Consequence: D1 can be stated sharply — at fixed `m` the only
  thing a blast can do inside `T_m` is add crossings, so the hole set and class count are non-decreasing
  by construction (A.10). ~60 lines (`Trellis._pip_orbit_cuts`, session plumbing, one law).
- **(N-c)** leave as is: `m` decided by `min_separation` and the zone polygon. Not recommended: the
  definition of `N` then depends on a resolution guard.

**Code today.** `is_strong_pip` is the PDF algorithm (one representative per crossing suffices because
the box lies inside the computed trellis). `cut_at_pip_orbit` is clean; `m` is decided as above.

**Change.** Slice B3 (N-a now; N-b if the author wants the sharp D1).

### A.5 Stage 4 — pseudoneighbor pairs

**Definition (PDF).** Consecutive crossings `x, x'` on a stable branch of `T_m(q0)` are a pair iff
`W^s(x, x') ∩ X = W^u(x, x') ∩ X = ∅`, `X = ∪_n M^n(N)`. Both open intervals are tested the same way:
for each crossing `r` and each `n` that lands `M^n(r)` on the interval's branch (closed form over the
branch cycle), the landing is the table's crossing when linked (exact, by id), else `r` read against the
pair's registered images in `r`'s own frame, else placed by β-scaling with the 2-D collision guard for
the pair's own orbit. Reference pairs are those in the window `(c_s(M^k q0), c_s(q0)]` on the pip's
branch; trajectories are their table iterates while inside `T_m(q0)`; no forward iterate gets a hole
(D15). A pair touching the anchor is no pair (the anchor is a point of `N` on the arc).

`DECISION S — the stable interval.` Today only the unstable interval is tested against `X`; the stable
interval is "consecutive on the computed branch", which excludes registered crossings only
(`Pseudoneighbor.py:329-470`; its Dev Notes and NUM-18 say so; the author chose "no change" on
2026-09-15, before the audit framed it as half a definition).
- **(S-a, recommended default)** test both intervals with ONE function
  `_interval_is_clear(trellis, x0, x1, stability)` (the stable version swaps the roles: landings on the
  STABLE branch, cdist scaled by `1/β`, the branch cycle of the stable side's fixed point). ~40 lines
  net, the Dev-Notes paragraphs on the asymmetry deleted. Consequence: NUM-18's orphan window pairs
  (p3, inversion) die; references may switch at other depths; BCG-4 / BCG-8 / DEC-25 are re-measured
  AFTER this (the audit's §6.2). A pair that dies is a discovery, not a regression (author, 2026-10-07).
- **(S-b)** keep the asymmetry and WRITE it as the definition ("N-consecutive on S, X-clear on U").
  Consequence: the code matches a definition that is not the PDF's; the orphan pairs stay.

**Change.** Slice B2.

### A.6 Stage 5 — holes

**Definition.** Each reference pair's bridge gets a hole beside it, on its `bridge_side` read off the
signs. The reference bridge is mapped BACKWARD: both ends move by `c ↦ c β^{-1}` on `advance_key(key, -1)`
(the table never holds these points, A.0), the bridge whose consecutive crossings bracket the span is
the container and gets a hole with the carried side (`side_after`), until the container was visited
(periodic) or none exists (left the trellis). Containment is decided in the canonical coordinate; its
only uncertainty is each cdist's own noise, flagged when a span end lies within `SCALING_RTOL` of a
container-branch crossing. A hole is combinatorial data: `(origin, iterate, bounding_ids, bridge_side,
openings)`. Where it is DRAWN is plotting.

**Code today.** `propagate_reference_holes` is this. Around it: the `from_table` flag read and dropped;
a direct punch of backward-iterate pairs plus two dedups and a fallback (`punch_holes` else-branch,
`seen`, `punched_regions`, `fallbacks`) for pairs that never exist (orchestrator probe: `iterate < 0`
pairs = 0 on all 11 law cases, k28 x8, p3 x6 and every step of the k28 / p3 / inversion series); hole
COORDINATES computed in the algorithm (`epsilon`, `_step_into`, `_stable_arc_midpoint`,
`_nearest_arc_point`, `_image_arc_midpoint`) though read by `plotting.py` only (`:285, :298`); the
reference hole carried back by the REAL inverse map (`_map_backward`, topology calling `map_inv`) to
fill `Hole.measured_side`, a flag the law `backward_holes_keep_bridge_side` already re-measures
independently; `_hole_openings` carries an inline copy of `Trellis.stable_key_at` plus an inversion
flip stated nowhere else. ~990 code lines for ~250 of algorithm; ~340 of geometric diagnosis
(`RowDiagnosis`, `SideDiagnosis`, `geometry_flags`, `bridge_row_diagnosis`, `side_diagnoses`), ~170 of
placement.

**Change.** Slices B4a (diagnostics out) and B4b (placement out, dead path out, the one WARNING).

### A.7 Stage 6 — the stable partition

**Definition (PDF + D6).** Each stable branch is cut at every hole bridge's endpoints; an interval is
open at an end iff a hole lies on that half, where the half is `outward ⇔ (sign(x) > 0) == (side == left)`;
a point excluded on both sides is a singleton; a branch with no crossing past its anchor is one element
(C1). Rows: `row(first) = left ⇔ sign > 0`, `row(second) = left ⇔ sign < 0`.

**Code today.** `_hole_marks_on_branch` + `intervals_from_marks`: ~75 lines, faithful.

**Change.** None to the definition; the four implementations of "element at a bridge end"
(`Trellis.element_for`, `BridgeClass._element_at_end`, `PartitionFamily.owner_at_end`,
`StablePartition._elements_at_bridge`) become one in B7 (two of them live in code B7 removes).

### A.8 Stage 7 — bridge classes and inertness

**Definition.** The class of a bridge is the unordered pair `{X, Y}` of partition elements its two ends
sit in on their rows, oriented anchor-outward (the homotopy class of arcs between `X` and `Y` relative
to the stable manifold and the holes). A bridge with both ends in one element is a loop and belongs to
the class of the bridge whose image chain holds it (D2); a loop no chain reaches is an orphan. The
alphabet is the classes with a member the minimal trellis keeps, closed under the word map (A.8c).

`DECISION I — inertness.` The chalkboard says "inert iff every image is a loop or inert"; the code
computes it as a fixed point over REGISTERED image evidence (`_resolve_inertness`: no evidence = active),
so it is monotone in evidence, not a property of the dynamics, and it is load-bearing: an inert hole
bridge is kept but NOT mapped in the minimal trellis. Orchestrator probe (every class forced active,
every hole bridge mapped, 14 builds): image bridges +0..2 (k28 x1 0→1, p3 0→2, nested 3→5,
inversion x0 0→2, heteroclinic 8→10), iterated elements +0..6 on 6 builds, **class-level words changed
on ZERO classes, reliability unchanged on every build.** Nothing in the record says why those cuts may
be omitted.
- **(I-a, recommended default)** *Inertness is a property of the word graph, read AFTER the words:* a
  class is inert iff no infinite forward path leaves it in the class transition graph (every path
  reaches the empty word), i.e. it is transient. The minimal trellis maps EVERY hole bridge. This is the
  chalkboard's definition made exact (induction on "every image is a loop or inert" is exactly "no
  infinite path"), it needs no evidence rule, and it removes inertness from the CONSTRUCTION: the extra
  image bridges are legitimate walls (`f(b)` cannot cross the image of a disjoint hole bridge, inert or
  not), so adding them can only remove wrong walks — which the probe confirms (no word moved).
  Consequence: iterated elements change on 6/11 cases (moot under R-a; refined names change under R-b);
  the inert letter series and the matrix's exclusion of inert classes stay as display rules.
  ~−60 lines in `BridgeClass.py`, ~−30 in `MinimalTrellis.py`, +25 in `SymbolicDynamics.py`.
- **(I-b)** keep the evidence fixed point and WRITE the reason inert images may be omitted. Nobody has.
- **(I-c)** drop inertness altogether (every class active, lettered `a..`). Harmless mathematically but
  loses the author's 2026-09-16/23 display rules for no gain over I-a.

**Change.** Slice B5.

### A.8b Stage 8 — the minimal trellis

**Definition.** The minimal trellis keeps every hole bridge and the bridges tiling the first forward
image of every hole bridge (I-a) — the arcs the word of any class can meet — with the stable manifold
whole; its arrangement is sparse (an absent unstable ray is an absent wall). An image pair no `Bridge`
spans is uncomputed manifold (D9/INT-9): skipped with a WARNING, never invented.

**Code today.** As stated except for the inert exclusion. Clean.

### A.8c Stage 8 — the word of a class

**Definition (the object).** Let `b` be any member of class `C = {X, Y}`. `f(b)` is an unstable arc from
`f(X)` to `f(Y)`; it meets the stable manifold `S` only at its ends and inside the window `S \ f(S)` (its
interior crossings are images of crossings of `b` with `f^{-1}(S) \ S`). The **itinerary** of `f(b)`
lists, at each crossing, the element on the entry row then the element on the exit row; consecutive
pairs are classes or their inverses (`a`, `a^-1`); pairs in one element are loops and drop. The
**word** `w(C)` is that reduced itinerary. It is a class invariant: two members are homotopic relative
to `S ∪ holes`, so are their images, and the reduced itinerary is the homotopy invariant of an arc in
the complement. (This is the statement the audit says is missing; both readings compute it.)

**Definition (the algorithm, D3).** `w(C)` is read off the dual graph of the minimal trellis: the walk
from the stable edges the landing span `f(X)` meets to those `f(Y)` meets (C14), crossing `S` only at
unified nodes (the window), never crossing an unstable edge, whose closed curve holds exactly the images
of the markers the member lobe holds (D13), and whose chords interleave with no decided walk of a class
on disjoint unstable branches (walls, C17). Among the walks satisfying these constraints the word is the
one of minimal length; several of minimal length = `ambiguous`; none within the search bound, or a bound
reached before uniqueness is established = `unresolved`. Singletons cannot walk: their word is the
registered chain (`trellis_itinerary`). The alphabet is closed under `w` (C9/D11): a pair named by a
word is a class, spanned or not, flagged "no computed bridge yet"; a class exists iff spanned or named.

`DECISION W — what the registered chain is.` D3: the walk is the answer, the chain evidence. The audit
measured 100 % agreement wherever a chain gives a verdict (except the already-flagged inversion blast 6).
- **(W-a, recommended default: keep D3, finish it)** The chain stays an oracle, in TESTS (the law
  `walks_agree_with_registered_images` on every case with a clean chain) and as a verification record;
  the search bounds become honest (B6a, decision-free): a walk found under a truncated budget is
  `unresolved` (today `is_reliable` ignores it: `max_walks=1` yields "unique" and reliable on
  heteroclinic), a span meeting no stable node is `unresolved` (today "the walk may leave from any of
  them"), members disagreeing on a marker is `unresolved` (today the first member is read). The
  evidence comparison (`registered_evidence`, `EvidenceDiagnosis`, C11/C21 noise detectors, ~300 lines)
  moves out of the algorithm into `topology/diagnostics.py` (B6b).
- **(W-b)** a surviving tie is decided by a member's clean chain when one exists (it IS `f(b)`).
  Changes D3; needs the C11/C21 detectors INSIDE the algorithm. Less simple.
- **(W-c)** the word is the registered chain of one member; D11 classes get no word until a bridge is
  computed. Deletes ~1,500 lines but reverses D3 and D11. Listed because the audit lists it.

**Change.** B6a (decision-free), B6b (after Q4 and B7).

### A.8d Stage 8 — the iterated partition

`DECISION R (shared with A.9) — is the refined level a goal?` The iterated partition exists to give
refined symbols `a_1, a_2, …`. Today:
- `IteratedHomotopyPartition.from_minimal` is an ENDPOINT rule (at each image-bridge crossing that ends
  an empty stretch — hole, lobe base, chord — open the flank facing it), rewritten three times against
  fixture expectations; its odd foreign-chord count (C19) is a parity violation treated as a WARNING;
  no mathematical object is named.
- every refined child inherits the parent's word verbatim (`SymbolicDynamics.py:762`), so the refined
  matrix is NOT computed: its rows per parent are identical and its spectral radius equals the
  unrefined matrix's. For k=10 `[[1,1],[1,1]]` is right because the horseshoe is the full shift.
- Orchestrator probe: landing spans at the HOMOTOPY level straddle nothing (296 landings over 14 builds,
  0 straddles, 0 unresolved), so under R-a `contained` is an assertion, not a status.

Candidates:
- **(R-a, recommended default)** Words in CLASSES; the transition matrix is over the active classes
  (`M[C][D]` = multiplicity of `D^{±1}` in `w(C)`; k=10: `[[2]]`, radius 2, the full 2-shift). The
  iterated partition and refinement are REMOVED: the dual graph's stable nodes carry the homotopy
  element of their edge (edges still end at image-bridge crossings), the walk reads homotopy elements,
  `land_element` returns the span only. Deletes `IteratedHomotopyPartition`, `Cut`, `_empty_stretches`,
  `refine`, `RefinedClass`, `match_members`, `refined_rules`, the iterated naming (~1,000 lines) and the
  `R_1^2`-style names. Figures lose the refined labels. The level-2 object, if wanted later, is (R-c).
- **(R-b)** keep the iterated partition as LABELS only, matrix over classes. Honest, but keeps ~1,000
  lines whose definition is still the endpoint rule (and I-a changes those labels on 6 builds).
- **(R-c)** define the refined level: the iterated partition is `P ∨ f(P)` — `S` cut at every crossing
  of the images of the hole bridges, a cut point belonging to the side away from the empty footprint
  (`∪ holes ∪ f(hole lobes) ∩ S`) — and `f(a_i)` is computed by walking at level 2. The audit's
  set-theoretic candidate reproduces the endpoint rule on the 11 shallow cases and cuts MORE at depth.
  New mathematics and a new stage; its own plan after R-a, if the author wants the refined matrix.

**Change.** Slice B7.

### A.9 Stage 9 — the transition matrix

**Definition.** Nodes = the active classes (inert classes are transient, A.8); `M[C][D]` = the number of
tokens `D` or `D^{-1}` in `w(C)`. No spectral radius shown (author 2026-09-23) until its meaning is
settled; with R-a it is well defined (the entropy bound of the class shift) and can return as a
documented number if the author wants (Q6).

### A.10 Stage 10 — recompute

**Definition (D1).** The topology is a pure function of `(registry, bridges, pips, m)`. Grow, blast and
recompute commute (A.2). At fixed `m` (N-b) a blast can only add crossings to `T_m(q0)`, so the hole set
and the class count are non-decreasing; a change of `m` is a discovery and may change anything
(author rule 2026-10-07). Zones, trims and letters feed no topology output (audit probe).

**Code today.** True except on the collapse path (A.2) and for the depth coupling (A.4).

### A.11 Session and API (no mathematics)

The session's three concerns: hold the workbench, cache one `Topology` per `(generation, selection,
pips[, depth])`, remember pips and letters. Everything else (243 lines of plot delegates, 182 of
cartoon geometry, 264 of fan-outs that feed no product, zone stitching) moves to `plotting` or goes.
`resonance_zone()` adopts a pip only if it is a candidate; a pip that stops being one is forgotten
through `forget_strong_pip`, never a `topology()` that raises with no way out.

---

## B. Slices

Order of execution. Decision-free slices (full form below): **B0 → B6a → B8 → B4a → B4b → B10.**
B1 is in full form too but waits on **Q13** (it puts a resolution criterion inside the definition of
N; it goes first once answered). Waiting on the author (skeleton): B2 (Q1), B3 (Q2), B5 (Q3), B7 (Q5),
B6b (Q4, after B7 under R-a), B9, B11 (Q7), B12 (last).

### B0 — `oracle-laws` (tests only)

**The answer.** The suite has 98 laws, about a third of which can only fail if the code reports its
own failure, and none re-checks the PDF definitions (pip box, pair intervals, hole containment) on a
real build by an independent computation, nor the one textbook answer (k=10 = the full 2-shift). Four
laws close that; three of them are rules and need no decision; the k=10 oracle is a case-specific pin
and needs Q11. One case, `inversion_period2`, is a law case with an unresolved class
(`e`: "a landing is a singleton element and no member has a registered image chain"); one more growth
step on EITHER side resolves it with identical counts (probe: 6 classes / 10 holes before and after) —
Q12.

**Scope.** Slug `oracle-laws`. Files: `tests/helpers/laws.py`, `tests/cases.py`,
`tests/invariants/test_law_partition.py` (docstring), `tests/invariants/test_law_symbolic.py`
(docstring). No `src/` file.

**Steps.**
1. `check_chosen_pip_is_strong(case)` in layer `partition`: for each selected fixed point, re-derive
   the box test from the registry and `FixedPoint` alone (no `StrongPip` import): for every crossing
   `r` of the fixed point, `d = (pos(r.s_key) - pos(q0.s_key)) % k`, `c_s' = c_s(r) β^d`,
   `c_u' = c_u(r) β^{-d}`; skip `r` when `(c_u', c_s')` collides with `(c_u(q0), c_s(q0))` within
   `SCALING_RTOL` on both (the one identity assumption, stated); assert no `r` has
   `0 < c_s' < c_s(q0)` and `0 < c_u' < c_u(q0)`. Count the crossings tested (non-vacuous).
2. `check_reference_pairs_clear_of_X(case)` in layer `partition`: for each reference pair and BOTH
   open intervals (unstable and stable), for every crossing `r` of the arc's fixed point (D7: the
   unstable interval's cycle is the unstable side's; the stable interval's is the stable side's),
   enumerate the `n` whose landing branch is the interval's branch and whose scaled cdist falls in the
   open interval (finite by `β > 1`), skipping the table-linked landings that ARE `x0` / `x1` and the
   scaled ones colliding with an endpoint on both cdists; assert none remains. Scaling only, on the LAW
   cases only (shallow; the library's "own frame" reading exists for depth and is not reproduced).
   On base the STABLE half fails on `p3` and `inversion` (NUM-18's orphan pairs):
   `KNOWN_ISSUES[("reference_pairs_clear_of_X", name)] = KnownIssue("NUM-18 …", wp="B2")` for those,
   strict xfail, so B2 is falsifiable by the entry's removal.
3. `check_propagated_hole_inside_container(case)` in layer `partition`: for each hole with
   `iterate = n < 0`, carry the reference bridge's two END crossings back `|n|` times with the REAL
   `map_inv` (test side), project each onto the container bridge's polyline (nearest segment), and
   assert both projected cdists lie strictly between the container's two crossing cdists. Independent
   of `image_cdist` and `_containing_bridge`.
4. `check_k10_is_the_full_shift(case)` in layer `symbolic` (gated by Q11): letter-free — exactly one
   ACTIVE class; its word's active tokens are two occurrences of itself with opposite directions; the
   class matrix is `[[2]]`; while refinement exists, the refined matrix is `[[1,1],[1,1]]`.
   `NOT_APPLICABLE[(law, name)] = "the textbook answer is known for k=10 only"` for every other law
   case (the set-aside test keeps it visible).
5. (Q12, default yes) `build_inversion_period2(stable_steps=10)`; the docstring's "10 runs the stable
   manifold past the pip orbit" remark is moot since D10 and goes.

**Behaviour.** None in the library. The suite gains 4 laws × 11 cases; two strict xfails
(step 2 on `p3`, `inversion`) until B2.

**Not doing.** No deletion of tautological laws (B12 revisits them once the code they mirror is
rewritten); no deep-series oracle (scaling-only is wrong at depth by design, A.5).

**Rubric.**
1. Suite passes on the branch; `git diff --stat` touches only the four files.
2. Each new check imports nothing from `StrongPip`, `Pseudoneighbor`, `StablePartition`
   (`grep -n "from tanglepack.topology.StrongPip\|Pseudoneighbor import\|_containing_bridge"` on the
   new functions: none).
3. Teeth: with `is_strong_pip` monkeypatched to accept everything and the default pip forced to the
   second candidate, check 1 fails on `k10`; with `_unstable_interval_is_clear` patched to `True`,
   check 2 fails on `k28_two_blasts`; with `_containing_bridge` patched to return the NEXT bridge on the
   branch, check 3 fails on `k10`.
4. Check 2's stable half fails on base for `p3` and `inversion` and nowhere else (the KNOWN_ISSUES
   entries are exactly those two).
5. Check 4 passes on `k10` on base; `-rxX` lists its 10 not-applicable set-asides.
6. `inversion_period2` at `stable_steps=10`: `topology().dynamics.is_reliable` True, 6 classes, 10 holes
   (the probe's numbers; not pinned in a test, read by the reviewer).

**Questions.** Q11, Q12 (section E).

### B1 — `resolved-crossings` (waits on Q13)

**The answer.** The blast path votes on near-coincident detections (`_collapse_noise_crossings`:
odd run → median, even run → dropped) and every other path does not, so a plain recompute after deep
blasting changes `N` (k28 x8: 248 → 252). Deleting the vote is not enough: with it gone, k28 x8 goes
from 24 holes / 15 classes to 40 / 23 with diverging words, both "reliable" — the pipeline has NO rule
for a detection the polylines cannot resolve once it is registered. A.2 states the rule (Q13 (a)): a
detection is in `N` iff the polylines resolve it; an unresolved one is registered, reported with the
remedy "refine", and excluded on every path the way a crossing past the pip cut is. The topology then
commutes with recompute by construction, and k28 x8's answer is the collapsed one, now for a stated
reason instead of a parity vote.

**Scope.** Slug `resolved-crossings`. Files: `src/tanglepack/numerics/Tangle.py`,
`src/tanglepack/numerics/Intersection.py` (`resolved`), `src/tanglepack/numerics/TangleWorkbench.py`
(`compute_intersections`, `create_bridges`, `report_resolution`), `src/tanglepack/numerics/IterateInference.py`
(sources and images from `N`), `src/tanglepack/topology/Trellis.py` (branch orderings from `N`),
`tests/unit/numerics/test_noise_crossing_collapse.py` (deleted) → `test_resolved_crossings.py`,
`tests/cases.py`, `tests/regression/test_class_count_under_blasting.py`. Eight files — flagged: one
rule applied at the four places that read the registry as `N`.

**Steps.**
0. Probe (test-side script kept under `tests/_tools/`): for every law case and series step, list the
   same-curve detection pairs closer along the unstable polyline than the chord deviation there, and the
   minimum per-crossing ratio; attach the output to the report (rubric 4).
1. `Intersection.resolution_ratio` stays as computed; add `Intersection.run_separation` = the distance
   along the unstable polyline to the nearest other detection of the same `(u_key, s_key)` on the same or
   an adjacent unstable segment (None when alone), set by `Tangle.resolve_crossings` and
   `populate_intersections_for_manifold` alike (one helper, called on both paths).
2. `IntersectionRegistry.resolved_ids(threshold)` = ids with `is_synthetic` or (`resolution_ratio ≥
   threshold` and `run_separation` None or `> max(h_u, h_s)`); `unresolved_ids(threshold)` the
   complement. `N` is `resolved_ids` wherever the registry is read as the crossing set:
   `create_bridges` (cuts at resolved crossings only — an unresolved detection is inside one bridge),
   `infer_iterates` (sources and candidate images), `Trellis.from_workbench` (branch orderings;
   `own_intersection_ids`). Four call sites, one name.
3. Delete `_collapse_noise_crossings` and `_purge_crossing` (its only caller, `Tangle.py:626`);
   `populate_intersections_for_manifold` returns every detection.
4. `report_resolution`: the unresolved set is ONE WARNING per call naming the detections (ids,
   coordinates, ratio, separation) and the remedy ("refine: lower `area_cutoff` or grow the stretch;
   NUM-3/4"); it stays a flag — the exclusion is the definition, not the flag.
5. Delete `test_noise_crossing_collapse.py`; `test_resolved_crossings.py` builds a hand-made polyline
   zig-zag across a stable segment (3 detections within the chord deviation) and asserts: all three
   registered, none in `resolved_ids`, one bridge across them, no iterate link to them, the WARNING
   once; and a clean transversal crossing with ratio ≥ 3 is resolved.
6. `cases.blast_series` records `recompute_crossings` per step (one `compute_intersections(…,
   preserve_ids=True)` + `rebuild_bridges()` per step) and `test_recompute_keeps_the_crossing_set(name)`
   over `cases.SERIES` asserts it equals the step's `len(N)` (resolved crossings), at every step.
7. `resolution_ratio_threshold` keeps its value 3 and gains a docstring: the derived bound is 1
   (separation equals the polygon error), 3 is the margin; it is part of the definition of `N` (A.2).

**Behaviour.** Law cases: clause 1 holds nothing out (`unresolved_crossings(3)` is empty at build on all
11; the collapse fired 0 times there); clause 2 is NOT yet measured on the law cases (no per-detection
chord deviation without recomputing it), so step 0 of the slice is a probe that prints, per law case,
every pair of same-curve detections closer along the unstable polyline than the chord deviation; the
expected answer is none, and any hit is a law-case behaviour change the report names. Series: p3 (9),
inversion (4), k32 (10), nested_outer (11): recompute already exact, no runs. **k28 at 8 blasts is the
one case the rule decides**, by clause 2: the outer members 1, 2, 250, 251 of the two odd runs are held
out (their ratios 624 .. 4642 are high; clause 1 does nothing to them), AND 148 and 249 are held out by
clause 1 (below threshold today, fed to the topology by the base). Six exclusions, not the base's two:
the slice RE-DERIVES k28 x8's holes / seed / closed and prints them beside the base's 24 / 15 / 23 and
the uncollapsed 40 / 23 / 31; the claim is "the topology of the resolved crossings", whatever the
numbers, and the report says which of the base's holes depend on 148 or 249. The B8 collision on
k28 x8 (125 and 251 claiming 213) disappears with 251. Threshold: unchanged at 3; heteroclinic's
minimum ratio 3.04 is recorded in the rubric so a drift under the margin is caught, and the margin's
value is Q13's second part.

**Not doing.** No local refinement (option (b), NUM-3/4: its own slice if the author wants the remedy
automated); no change to `DEDUP_RTOL`; no change to the same-stability flags (they stay flags: a u×u
hit is never a crossing anyway).

**Rubric.**
1. Suite passes on the branch; `git diff --stat` touches only the eight files.
2. `grep -rn "_collapse_noise_crossings\|_purge_crossing\|Collapsed" src tests` → nothing; `grep -rn
   "resolved_ids" src` → exactly the four readers in step 2 plus the registry.
3. `test_recompute_keeps_the_crossing_set[k28_series]` FAILS on base (248 ≠ 252 at blast 8) and passes
   on the branch; passes on both for the other four series.
4. Fingerprint of every law case equal on base and branch; `unresolved_ids(3)` is EMPTY on every law
   case under BOTH clauses (the step-0 probe's output is attached to the report); the minimum
   per-crossing ratio on `heteroclinic` is printed (3.04 on base) — a value under 3 is a behaviour change
   the report must own.
5. k28 at 8 blasts on the branch: `unresolved_ids` holds exactly `{1, 2, 250, 251, 148, 249}` (the
   first four by clause 2 alone — asserted by evaluating clause 1 on them: all ≥ 3 — the last two by
   clause 1); recompute leaves `N` unchanged; holes / seed / closed and letter-free words are printed
   beside the base's (24 / 15 / 23) and the uncollapsed build's (40 / 23 / 31) with the holes that
   depended on 148 / 249 named; the WARNING fires once per `compute_intersections` /
   `blast_zone(diagnose=True)` naming the six.
6. The hand-built zig-zag test (step 5) fails on base (base registers and cuts at all three).
7. `Tangle.py` shorter by ≥ 60 lines; no parity test (`% 2`) anywhere in `Tangle.py`.

**Questions.** Q13 (section E), both parts: the criterion, and the margin (3 today, derived bound 1,
heteroclinic at 3.04).

### B6a — `walk-bounds-honest`

**The answer.** Three places let the walk return a word it has not established: (i) a budget that ran
out AFTER a walk was found leaves status `unique` with `truncated=True`, and `is_reliable` reads only
`ambiguous`/`unresolved` (`SymbolicDynamics.py:839-848`; `max_walks=1` makes heteroclinic "reliable");
(ii) a landing span meeting no stable node "may leave from any of them" (`DualWalk.py:729-757`); (iii)
two members holding different markers are a WARNING and the first is read (`:815-828`). Under A.8c
each is `unresolved` with its reason. This is W-a's decision-free half: it changes no definition, it
stops claiming what is unproven.

**Scope.** Slug `walk-bounds-honest`. Files: `src/tanglepack/topology/DualWalk.py`,
`src/tanglepack/topology/SymbolicDynamics.py`, `tests/unit/topology/test_dual_walk.py`,
`tests/unit/topology/test_symbolic_dynamics.py`.

**Steps.**
1. `_nodes_meeting` returns `[]` when nothing meets; `shortest_walks` returns
   `WalkSearch(start, goal, "span_meets_no_node")` (one WARNING kept, worded as a resolution flag).
2. `lobe_image_markers` returns `None` and sets a flag the search records (`members_disagree`) when two
   member lobes hold different marker images; `_resolve_itinerary` makes the class unresolved with
   reason "members hold different markers: the class is not one homotopy class (resolution)".
3. A `WalkSearch` with `truncated=True` is never `unique`: `_resolve_itinerary` treats
   `search.truncated` as unresolved ("search budget reached before uniqueness was established"). The
   default `max_walks=64` is kept (nothing truncates on any case or series; audit: `max_walks=8`
   changes nothing).
4. `is_reliable` docstring: unchanged formula (it reads `itinerary is None`), now covering the three.

**Behaviour.** None on any law case, deep run or series step (none of the three fires: the WARNING
strings "meets none of its", "hold different markers" and the truncated path appear in no recorded
log; rubric 3). On a hand-built graph: statuses change as stated.

**Not doing.** Statuses are not renamed (B6b collapses the enumeration); the evidence comparison is not
moved (B6b); an UNPLACED marker image (`image_marker` None) stays informational — it weakens the
constraint but never picks a wrong walk; no change to `MAX_DETOUR`.

**Rubric.**
1. Suite passes on the branch; `git diff --stat` touches only the four files.
2. `symbolic_dynamics(..., max_walks=1)` on `heteroclinic` yields `is_reliable` False on the branch and
   True on base (the audit's probe, re-run by the reviewer).
3. Letter-free `rules` of every law case and of `p3_6_blasts`, `p3_16`, `k28` at 8 blasts equal on base
   and branch.
4. Three new unit tests on hand-built dual graphs: a span meeting no node → `unresolved`; two members
   disagreeing on a marker → `unresolved`; a budget of 1 with two candidate walks → `unresolved`, each
   failing on base.
5. `grep -n "may leave from any of them" src` → nothing.

### B8 — `unresolved-sets`

**The answer.** Two warn-and-continue sites decide something the definition leaves undecided: an
iterate-table collision is resolved by registry order (`IterateInference.py:125-132`), and a zone's
pip is adopted without a candidacy test (`TangleSession.py:1509`) while a pip that stops being a
candidate makes `topology()` raise with no API to forget it (`Topology.py:_choose_pip`). Both become
data: `registry.unresolved_links` and `forget_strong_pip`.

**Scope.** Slug `unresolved-sets`. Files: `src/tanglepack/numerics/IterateInference.py`,
`src/tanglepack/numerics/IntersectionRegistry.py`, `src/tanglepack/loom/TangleSession.py`,
`tests/unit/numerics/test_iterate_inference.py` (or the existing inference test module),
`tests/facade/test_session_topology.py`.

**Steps.**
1. `IntersectionRegistry.unresolved_links: dict[int, tuple[int, ...]]` (source id → the candidate
   image ids), cleared at the start of each `infer_iterates` pass over its sources, bumping
   `_generation` when it changes.
2. `_bracketed_image`: a multi-match records `unresolved_links[src] = tuple(images)` (INFO kept).
3. The collision: when `(j, -1)` is already in the table, record BOTH sources under `j`'s candidates
   (`unresolved_links[i] = (j,)`, `unresolved_links[table[j, -1]] = (j,)`) and UNLINK the earlier one
   (`register_iterate` gains `remove_iterate(i, 1)`), so the table never holds a link decided by order.
   One WARNING naming the three ids.
4. `TangleSession.forget_strong_pip(fp)`: drops the remembered pip; the next `topology()` chooses the
   default again. `resonance_zone()` remembers the zone's pip only when `trellis._is_candidate(pip)`
   and WARNS otherwise (the zone is still defined). `Topology._choose_pip`'s `ValueError` names
   `forget_strong_pip`.
5. Law (layer `crossings`, in the existing iterate-table check): `unresolved_links` is disjoint from the
   linked sources, and every unresolved source has ≥ 1 candidate.

**Behaviour.** Law cases: none (no collision fires; no recipe adopts a non-candidate pip). k28 x8
recompute BEFORE B1: the one collision (sources 125 and 251 claiming image 213, whose preimage is 250;
none of the three has another link — orchestrator probe, no cascade) leaves 2 fewer `+1` links than the
base's 144 and 2 entries in `unresolved_links`. AFTER B1 the collision disappears with the unresolved
detections 250/251, and no collision fires anywhere in the corpus: B8's collision rule is then a safety
net the definition keeps honest, not the fix. `resonance_zone()` on a non-candidate: a WARNING and no
pip adoption instead of a later `ValueError` from `topology()`.

**Not doing.** Odd foreign-chord counts and I2 on image bridges (`PartitionFamily.py:760, 960-967`)
are left to B7, which deletes or keeps that code by Q5; marker disagreement is B6a.

**Rubric.**
1. Suite passes on the branch; `git diff --stat` touches only the five files.
2. `grep -n "already has the preimage" src` → the one WARNING, followed by `remove_iterate` /
   recording, no `continue` that leaves the earlier link in place.
3. Unit test: two sources whose brackets name one image → neither linked, both in `unresolved_links`;
   fails on base (base links the first).
4. Facade test: `resonance_zone(id)` with `id` a non-candidate → `chosen_strong_pip(fp)` unchanged, one
   WARNING, `topology()` does not raise; fails on base (base raises in `topology()`).
5. Facade test: `set_strong_pip`, then `forget_strong_pip`, then `topology()` → the default pip.
6. Fingerprints of every law case equal on base and branch.

### B4a — `holes-algorithm-vs-diagnostics`

**The answer.** `StablePartition.py` (1,947 lines, ~990 of code) holds ~250 lines of algorithm (A.6,
A.7) and ~340 of geometric self-diagnosis: `_side_of`, `_anchorward_look`, `_approach_displacement`,
`_oriented_bridge_polyline`, `_build_oriented_bridge_polyline`, `_bridge_side_of`, `RowEnd`,
`RowDiagnosis`, `bridge_row_diagnosis`, `bridge_row_violation`, `check_bridge_rows_consistent`,
`_orbit_side`, `SideDiagnosis`, `side_diagnoses`, `bridge_side_violations`,
`check_holes_share_bridge_side`, `geometry_flags`, `_warn_geometry_flags`, plus `_map_backward` (the
topology calling the real inverse map to fill `Hole.measured_side`, a flag the law
`backward_holes_keep_bridge_side` re-measures independently). The algorithm module should read as the
PDF; the readings go to `topology/diagnostics.py`, mirroring `numerics/diagnostics.py` (flags, never
discards).

**Scope.** Slug `holes-algorithm-vs-diagnostics`. Files: `src/tanglepack/topology/StablePartition.py`,
`src/tanglepack/topology/diagnostics.py` (new), `src/tanglepack/topology/Trellis.py`
(`punch_holes(diagnose=)` imports), `src/tanglepack/topology/TopologyResults.py` (`Hole`),
`src/tanglepack/topology/__init__.py`, and the test imports: `tests/helpers/laws.py`
(`_side_of`, `bridge_row_violation`), `tests/unit/topology/test_stable_partition.py`
(`_bridge_side_of`, `_side_of`), `tests/plotting/test_plot_topology.py` (`_side_of`). Eight files —
flagged: a MOVE with import updates, no logic change.

**Steps.**
1. Create `topology/diagnostics.py` with the 18 names above moved verbatim (one commit "move", then
   one commit "tidy" so the reviewer can diff the move as a rename).
2. `Hole` loses `measured_side` and `side_source` (the measured side is returned by
   `diagnostics.side_diagnoses(trellis)` on demand, computed there by `_map_backward`; `side_source`
   is `iterate == 0 → signs, else carried`, a function not a field).
3. `Trellis.punch_holes(diagnose=True)` calls `diagnostics.warn_geometry_flags(self)` after the
   algorithm; the two invariant checks (`check_bridge_rows_consistent`, `check_holes_share_bridge_side`)
   stay assertions and stay in the algorithm module (they are I1/I2, not geometry): keep them in
   `StablePartition.py`, ~45 lines.
4. `topology/__init__.py` re-exports the moved names from `diagnostics` (no caller outside tests
   imports them by module path except `Trellis.py`).

**Behaviour.** None. `Hole` loses two fields (no library reader outside the diagnostics;
`plotting.py` reads `coords`, `bridge_side`, `openings` only — verify by grep in the rubric).

**Not doing.** Placement, the direct-pair path and the WARNING (B4b); any change to `row_of_end`,
`_hole_openings`, `intervals_from_marks`.

**Rubric.**
1. Suite passes on the branch; `git diff --stat` touches only the eight files.
2. `StablePartition.py` ≤ 1,100 lines; `topology/diagnostics.py` holds the 18 names and imports
   nothing from `StablePartition` except `row_of_end`, `rows_of_bridge`, `owns_cdist`.
3. `grep -rn "measured_side\|side_source" src` → only `topology/diagnostics.py`.
4. `grep -rn "map_inv\|dynamical_system" src/tanglepack/topology/StablePartition.py` → nothing (the
   algorithm never calls the map).
5. Fingerprints of every law case equal on base and branch; `geometry_flags` output identical on the
   11 law cases (the reviewer runs `side_diagnoses` on both).

### B4b — `holes-are-combinatorial`

**The answer.** With the diagnostics out, what remains in `StablePartition.py` beyond the PDF is
(i) hole placement in the plane (`epsilon`, `_step_into`, `_stable_arc_midpoint`, `_near_far`,
`_nearest_arc_point`, `_image_arc_midpoint`, the `coords` arithmetic in `punch_holes` and
`_punch_in_bridge`, ~170 lines), read by `plotting.py` only; (ii) the direct punch of backward-iterate
pairs and its two dedups and fallback (`punch_holes` else-branch, `seen`, `punched_regions`,
`fallbacks`, ~60 lines) for pairs that never exist; (iii) the `from_table` flag read and discarded at
`:513`; (iv) `_hole_openings`' inline copy of `Trellis.stable_key_at`. A hole is `(origin, iterate,
bounding_ids, bridge_side, openings)`; where to draw it is `plotting.hole_point(trellis, hole)`.

**Scope.** Slug `holes-are-combinatorial`. Files: `src/tanglepack/topology/StablePartition.py`,
`src/tanglepack/topology/plotting.py`, `src/tanglepack/topology/TopologyResults.py` (`Hole.coords`
removed), `src/tanglepack/topology/Trellis.py` (`punch_holes(epsilon=)` gone),
`tests/unit/topology/test_stable_partition.py` (`_near_far` tests move to a plotting test or go with
the helper). `TangleSession.punch_holes` forwards `**kwargs` and names `epsilon` only in its docstring
(`TangleSession.py:1312`): one docstring line, counted in the five files.

**Steps.**
1. `punch_holes(trellis, pairs=None)`: punches the REFERENCE pairs only (`iterate == 0`); assert no
   recorded pair has `iterate < 0` (the orchestrator's probe: none on any build or series step; the
   assertion makes the deletion falsifiable); forward pairs skipped as today (D15). Delete the
   else-branch placement, `_near_far` where only placement used it, `seen`, `punched_regions`,
   `fallbacks` and the "has no punched bridge; propagating from its iterate" path.
2. `propagate_reference_holes`: the span is `tuple(sorted(image_cdist(...)[1] for end in ids))` with a
   one-line comment stating A.0's fact (never in the table, by construction); after
   `_containing_bridge` returns, WARN once per branch when either span end lies within
   `SCALING_RTOL * |c|` of a crossing of the container's branch (the container could be its
   neighbour: NUM-5's mechanism). `_punch_in_bridge` loses `carried` and `span` placement; `Hole` is
   built from ids, side and openings alone.
3. `_hole_openings` at an anchor end is the general D6 formula plus ONE identification, stated as such
   (orchestrator probe: non-inversion cases take `s(o,b) = stable_key_at` at every anchor end; the
   inversion cases take `s(o,1-b)` at 2/4, 2/2, 2/2 and 4/8 of theirs): *the two stable branches of an
   inversion point are one smooth curve through it (stage 0), so the anchorward half of `s(o,b)` at the
   anchor IS the outward half of `s(o,1-b)`.* The block names the branch by `stable_key_at` and applies
   that identification when the opening is the anchorward half; no "flip", no second rule.
4. `plotting.hole_point(trellis, hole, epsilon=0.05)` reproduces today's placement (reference: chord
   midpoint nudged toward the bridge; propagated: image sub-arc midpoint nudged toward the chord) using
   the moved helpers; `plot_holes` calls it (`plotting.py` reads `Hole.origin`, `Hole.iterate` and the
   point, nothing else: `:268-298`).
5. `Hole.coords` removed; `Trellis.punch_holes(epsilon=)` removed.

**Behaviour.** None in any topology output (holes are compared by `(origin, iterate, bounding_ids,
bridge_side, openings)` in every law and fingerprint). Figures: identical (same arithmetic, moved).
Logs: one new WARNING kind, firing on no law case (rubric 5).

**Not doing.** "One element at a bridge end" (B7); the `_containing_bridge` bisect itself; any change
to which bridge gets the hole.

**Rubric.**
1. Suite passes on the branch; `git diff --stat` touches only the listed files.
2. `StablePartition.py` ≤ 700 lines, and `grep -n "epsilon\|coords\|_step_into\|np.linalg" src/tanglepack/topology/StablePartition.py` → nothing.
3. The assertion in step 1 is present and a unit test with a hand-built `iterate=-1` pair hits it.
4. Fingerprints of every law case equal on base and branch (holes compared by their five fields).
5. The new WARNING appears in no law-case log; a unit test with a hand-built container whose neighbour
   crossing sits at `0.5 * SCALING_RTOL` from the span end triggers it.
6. `plot_holes` on `k28_two_blasts` renders hole markers at the same coordinates as base (the reviewer
   compares `hole_point` against the base's `Hole.coords` to 1e-12 on every hole).

### B10 — `dead-api`

**The answer.** Zero-user API is deleted; test-only API that serves as an ORACLE stays or moves into
`tests/helpers`. Dispositions (callers verified by grep):

| Name | Where | Callers | Disposition |
|---|---|---|---|
| `grow_x_times` | `ManifoldMachine.py:490` | scripts only (10) | delete; scripts use `grow_n_times` |
| `release_manifold` | `Tangle.py:711` | `TangleWorkbench` only (internal) | delete if the workbench path is itself unused (verify), else keep private |
| `trim_stable_at_intersection` | `ResonanceZone.py:473`, exported twice | none | delete + exports |
| `_get_iterate` | `ManifoldMachine.py:848` | none | delete |
| `_get_refined_point` | `ManifoldMachine.py:815` | `test_refinement.py` (one kernel test) | delete; the test reads the vectorised kernel's output for the same pair |
| `visualize_intersection_graph` | `graphviz.py:62`, workbench `:2042` | none | delete |
| `build_intersection_graph` + `IntersectionRegistry.graph` / `_rebuild_adjacency_edges` / `_consecutive_on_branch` | `graphviz.py:38`, workbench `:2031`, registry `:416-545` | `test_map_step_and_graph.py` | delete the registry's networkx graph (~230 lines) and `graphviz.py`'s graph builder; keep `graphviz.py` only if a plotting view survives (likely delete, 358 lines) |
| `on_cdist_range`, `from_branch`, `nearest_by_unstable_cdist`, `as_forward_array`/`as_backward_array` | registry | tests only / none | delete |
| `classify_bridges` | `TangleSession.py:1591` | `test_resonance_zone.py:115` | delete with its one test |
| `terminated_early` | `Blast.py:172` | `test_blast_budget.py` | delete; tests read `stop_reason` |
| `partition_element_for` | `TangleSession.py:1353` | `test_partition_elements.py` | B9 (session shrink), not here |
| `euler_characteristic`, `Arrangement.preimage_of`, `StableNode.other_face`, `preimage_bridges` | arrangement, dual graph, bridge iterator | law tier oracles | KEEP (oracles) |
| `refine_manifold` | `ManifoldMachine.py:627` | 27 script users | KEEP |
| `topology.__all__` (101 names) | `topology/__init__.py` | — | trim to the stage exports + plotting; B12 finishes |

**Scope.** Slug `dead-api`. Files: `src/tanglepack/numerics/ManifoldMachine.py`,
`src/tanglepack/numerics/IntersectionRegistry.py`, `src/tanglepack/numerics/TangleWorkbench.py`,
`src/tanglepack/numerics/graphviz.py` (deleted), `src/tanglepack/loom/ResonanceZone.py`,
`src/tanglepack/loom/Blast.py`, `src/tanglepack/loom/TangleSession.py`, the `__init__` exports, and
the tests named in the table. Flagged larger (deletions across ~10 files, no logic change).

**Steps.** One commit per table row, each deleting the name, its exports, its tests and its docs
mention (`CLAUDE.md` table rows updated in B12, not here).

**Behaviour.** None.

**Not doing.** Session API (B9); the plotting split (SIMP-11); anything the law tier reads.

**Rubric.**
1. Suite passes; `git diff --stat` touches only the listed files; net ≤ −900 lines.
2. For every deleted name `grep -rn <name> src tests scripts stress docs` → nothing (scripts updated).
3. `import tanglepack` and `from tanglepack.topology import *` succeed; `len(topology.__all__)` ≤ 70.
4. The law tier is unchanged (no `laws.py` edit in the diff).

### Skeletons (wait on the author)

| # | Slug | Goal | Implements | Depends on |
|---|---|---|---|---|
| B2 | `pseudoneighbor-both-intervals` | One `_interval_is_clear(stability)`; the stable interval tested against X; Dev Notes shortened; B0's KNOWN_ISSUES entries removed; re-measure BCG-4/BCG-8 series and report before DEC-25. | A.5 (S-a) | B0, Q1 |
| B3 | `depth-of-N` | `min_separation` never drops a pip-orbit child (N-a); optionally `depth=` (N-b) and the sharp D1 law at fixed depth. | A.4 | Q2 |
| B5 | `inertness-is-transience` | Minimal trellis maps every hole bridge; inert read off the word graph after the words; evidence bookkeeping deleted. Law: class-level words unchanged on every law case (the probe's result, made a rubric). | A.8 (I-a) | Q3, B4b |
| B7 | `words-in-classes` | Remove the iterated partition and refinement (R-a); matrix over classes; dual-graph nodes carry homotopy elements; `contained` an assertion; one "element at a bridge end"; figures follow. LARGE (~1,000 lines deleted). | A.8d, A.9 | Q5, B5 |
| B6b | `word-finished` | Statuses → `resolved / trivial / ambiguous / unresolved(reason)`; evidence comparison out of the algorithm into `topology/diagnostics.py`; the definition of the word in the module docstring; walls and markers kept as theorems. | A.8c (W-a) | Q4, B7 |
| B9 | `session-three-concerns` | Plot delegates and cartoon geometry to `plotting`; fan-outs and `partition_element_for` deleted; zones computed from `Trellis.pips`. LARGE (~700 lines moved/deleted). | A.11 | B7, B8 |
| B11 | `manifold-arrays` | Arrays `xy, c, pre_xy, mapped` per branch; same refinement kernel; bridges as index ranges. LARGE, numerics-only; its own plan. | A.1 | Q7, B1 |
| B12 | `one-design-note` | The 303 decision ids leave `src/` for one design note; CLAUDE.md rewritten as section A (one paragraph per stage); the tautological laws replaced or deleted now that the code they mirrored is gone; chalkboard divergence table closed. | all | last |

---

## C. What is deleted overall

Decision-free slices (firm):

| Slice | What | Lines |
|---|---|---|
| B1 (Q13) | `_collapse_noise_crossings`, `_purge_crossing`, `test_noise_crossing_collapse.py`; `resolved_ids` added at four readers | −70 src +60 src, −80 tests +80 tests |
| B6a | fallbacks replaced by statuses | ±20 |
| B8 | first-wins branch; `unresolved_links`, `forget_strong_pip` added | +40 |
| B4a | geometric diagnostics MOVED to `topology/diagnostics.py`; `Hole.measured_side`, `side_source` | −340 moved, −25 |
| B4b | placement MOVED to `plotting.py`; direct backward-pair path, dedups, fallback, `from_table` flag, `Hole.coords`, `epsilon` | −170 moved, −90 |
| B10 | dead API, registry networkx graph, `graphviz.py` | −900 (src + tests) |

Pending the author (estimates):

| Slice | What | Lines | Condition |
|---|---|---|---|
| B2 | the asymmetry's Dev Notes, one function for both intervals | −40 | Q1 |
| B5 | inertness evidence bookkeeping | −90 | Q3 |
| B7 | `IteratedHomotopyPartition`, `Cut`, `_empty_stretches`, `refine`, `RefinedClass`, `match_members`, `refined_rules`, iterated naming, three of four "element at a bridge end" | −1,000 | Q5 (R-a) |
| B6b | walk statuses, 3-valued logic; evidence comparison moved out | −250, −300 moved | Q4 |
| B9 | session plot delegates, cartoon geometry, fan-outs, `partition_element_for` | −700 (≈250 moved) | — |
| B11 | linked-list data model → arrays | ~1,200 → ~500 | Q7 |
| B12 | 303 decision-id citations in `src/`; tautological laws | −300 comments, tests ±0 | — |
| not planned | walk machinery (only under W-c) | −1,500 | not the default |

---

## D. Risks — where a slice could change the mathematics, and how it proves it did not

- **B0** changes nothing; its risk is a FALSE oracle. Proof: rubric 3 (teeth: each check fails when
  the library's stage is sabotaged) and rubric 4 (the stable-interval half fails exactly where NUM-18
  predicts).
- **B1** changes the DEFINITION of `N` (Q13): an unresolved detection is out until refined, on every
  path. The mathematics it could change: a genuine transversal crossing with `ratio < 3` or a
  neighbour closer than the chord deviation would be held out of `N` until refinement. Proof that it
  did not on the corpus: `unresolved_ids` empty on every law case (rubric 4) and on p3 / inversion /
  k32 / nested_outer at every series step (the collapse never fired there); on k28 x8 the held-out set
  is exactly the four detections the base's vote dropped, and the topology equals the base's (rubric 5).
  Proof that the change is intended where it bites: the two uncollapsed runs are even (a Z-kink of the
  polyline across the stable curve), and the base already excluded them — by parity; the branch excludes
  them by resolution and SAYS so. **The threshold becomes load-bearing:** `ratio = separation / chord
  deviation`, so 1 is the derived bound; 3 is a margin that was tuned to keep the fixtures silent
  (`TangleWorkbench.py:231`) and heteroclinic's minimum (3.04) sits 1 % above it. The rubric records
  that minimum; the margin's value is the author's (Q13, second part). Clause 2 has no free constant
  (the chord deviation is the polyline's own).
- **B6a** can only turn a claimed word into `unresolved`. Proof: letter-free `rules` equal on every
  law case, deep run and series step (rubric 3); the three new statuses are exercised on hand-built
  graphs only.
- **B8** changes which links the table holds when a collision occurs (now neither, was the first).
  Proof: law cases' fingerprints equal (no collision there); the k28 x8 recompute has 2 fewer `+1`
  links and 2 entries in `unresolved_links`; the bijection law passes on both.
- **B4a/B4b** are moves. Proof: fingerprints equal on every law case (holes by their five fields);
  `side_diagnoses` identical; `hole_point` equals the base's `Hole.coords` to 1e-12. The one new
  WARNING fires on no law case and is exercised by a hand-built container.
- **B10** is deletion of unused code. Proof: grep-clean; the law tier untouched.
- **B2** (pending) changes which pairs survive at depth: intended, a discovery under the PDF
  definition. Proof: B0(iii) passes on every law case; the series are re-measured and each lost or
  gained reference is listed with its witness (`reference_witness`) before DEC-25 is decided.
- **B3** (pending) changes `m`: intended. Proof: the pip orbit never stalls on k32 (`f^10(q0)`
  registered at blast 9); under N-b the class count is non-decreasing at fixed depth by the sharp law.
- **B5** (pending) changes iterated elements on 6/11 cases; class-level words must not move. Proof: the
  orchestrator's probe made a rubric (letter-free `rules` equal on all 14 builds).
- **B7/B6b** (pending) must not change any class-level word. Proof: letter-free `rules` equal on every
  law case, deep run and series step; the k=10 oracle; the registered-chain oracle on every clean chain.
- **B11** (pending) is the only slice that could change the polylines. Proof: identical crossing sets
  and cdists to `DEDUP_RTOL` on every case.

---

## E. Open questions for the author (each with a default so "yes to all defaults" is a usable answer)

1. **Q1 (DECISION S) — Test the open STABLE interval of a pseudoneighbor pair against X, exactly as the
   unstable one is tested?** Default: **yes** (S-a). You chose "no change" on 2026-09-15 as a known
   asymmetry; the audit shows it is half the PDF definition and that the orphan pairs it lets through
   (NUM-18) are no pairs by your own definition. Consequence: some reference pairs die at depth; the
   series are re-measured before DEC-25.

2. **Q2 (DECISION N) — Who fixes the depth `m` of `T_m(q0)`?** Default: **N-a** (the blast never drops a
   child holding a pip-orbit point; `m` = deepest registered). Add **N-b** (an explicit `depth`) if you
   want D1 as a theorem at fixed depth rather than a measured law.

3. **Q3 (DECISION I) — Inertness as transience in the word graph, computed AFTER the words, with every
   hole bridge mapped in the minimal trellis?** Default: **yes** (I-a). Measured: no class-level word
   moves on 14 builds. Alternative: write down why an inert class's image bridges may be left out of the
   partition (nobody has).

4. **Q4 (DECISION W) — Keep D3 (the walk is the word; the registered chain is an oracle)?** Default:
   **yes** (W-a). W-b (a clean chain decides a tie) changes D3; W-c (the chain is the word) reverses D3
   and D11 and deletes the walk. The decision-free half (honest bounds, B6a) proceeds either way.

5. **Q5 (DECISION R) — Is the REFINED transition matrix a goal?** Today it is inherited, not computed.
   Default: **R-a**: words in classes, matrix over active classes (k=10 → `[[2]]`), iterated partition
   and refinement removed (~1,000 lines), with the level-2 object (`P ∨ f(P)`, `f(a_i)` walked at level
   2) as a later plan if you want it. R-b keeps the refined names as labels only.

6. **Q6 — The spectral radius.** With R-a the class matrix's radius is well defined (entropy bound of
   the class shift). Default: **still not shown** (your 2026-09-23 rule) until you say otherwise.

7. **Q7 — The manifold data model (arrays instead of linked lists, B11).** Default: **later, its own
   plan**, after the definitional slices; it changes no definition and is the riskiest code change.

8. **Q8 — D1 restated.** Default: "at fixed depth `m` the hole set and class count are non-decreasing
   under blasting; a deeper `m` is a discovery and may change anything" (covers BCG-8's blast-6 drop:
   `m` went to 11 there). Does D1 also bind the CLOSED alphabet? Default: **no** (k32 45 → 39 at blast 9
   is the refined words naming fewer unspanned pairs; under R-a the question dissolves).

9. **Q9 — D11 and D14 together.** Accepted by the orchestrator as the definition: "a class exists iff a
   computed bridge spans it or the word of an existing class names it". Recorded here for your veto.

10. **Q10 — Recipe constants with no stated reason** (zones at `f(q0)` not `q0` for k28/k32/inversion;
    heteroclinic 16/10; `outer_blasts=2`). Default: **keep the numbers, write the reason** in
    `henon_cases.py` where one exists, else mark "chosen to make the case closed". No slice depends on it.

11. **Q11 — May the suite pin ONE textbook result?** Your 2026-10-05 policy forbids case-specific pins
    ("arbitrary and irrelevant"). The k=10 horseshoe is not arbitrary: its symbolic dynamics is the full
    2-shift, and the real build gives it today. Default: **yes, as the single exception**, stated in the
    policy docstring as "a textbook answer is an oracle, not a fixture fact". If no, B0 drops step 4.

12. **Q12 — `inversion_period2` at one more stable step?** The law case has an unresolved class (`e`,
    a singleton landing with no registered chain); `stable_steps=10` resolves it with the same 6 classes
    and 10 holes. Default: **yes** (the 9/9 choice predates D10's cut, whose remark the docstring still
    carries). If no, the oracle laws tolerate an unresolved class (they are stage 3-5 checks and the
    k=10 law), and nothing else changes.

13. **Q13 (DECISION, A.2) — May a resolution criterion enter the definition of N?** Today "flags never
    discard" (NUM-4 deferred), except that the blast path votes by parity on near-coincident detections.
    The vote decides k28 at 8 blasts (24 vs 40 holes, 15 vs 23 classes, both "reliable"); nothing else
    in the corpus is touched. Default: **(a)** — a detection the polylines cannot resolve (ratio below
    threshold, or a neighbour of the same two curves closer than the chord deviation) is registered and
    reported with the remedy "refine", and is NOT a crossing of `N` on any path, as a crossing past the
    pip cut is not. Alternatives: (b) refine automatically until it resolves (the remedy; its own slice,
    and it still needs (a) for what it has not refined yet); (c) keep a parity vote on every path (not
    recommended). **Second part — the margin.** `resolution_ratio = separation / chord deviation`; the
    derived bound is 1 (the crossing's displacement equals the polygon error); 3 is an undeclared margin
    tuned so the fixtures stay silent, and the heteroclinic law case's least-resolved crossing is at
    3.04. Under (a) the margin decides membership of `N`. Default: **keep 3** and record the corpus
    minimum in the rubric (a drift under it is a visible behaviour change), OR say the value you want
    (1 = the derived bound; 2 = a factor-two margin). Note clause 2 (a run closer than the chord
    deviation) has no free constant and is the clause that decides k28 x8.

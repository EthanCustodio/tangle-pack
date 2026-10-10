# The pipeline at the chalkboard

This page states the whole tanglepack pipeline the way it is drawn at a chalkboard: one rule per stage, its input and output, and the invariants it must keep. The code is judged against it. When the code does something this page cannot explain, either the page is missing a rule or the code is carrying complexity it should not have. Issues found that way go in [`open_issues.md`](open_issues.md).

The first draft was written on 2026-10-07 by the audit workflow's chalkboard review. It has been adjusted to the author's decisions of the same day (D1-D8 in `open_issues.md`), of 2026-10-08 (D9-D12) and of 2026-10-09 (D13, D14); the divergence table below is final for the open-issues round (WP12).

## Notation

- M is an area-preserving planar map. A saddle orbit is z₀..z_{p−1}, with p its minimal period.
- λ_u is the unstable eigenvalue of DM^p, and β = |λ_u|^{1/p} is the per-step factor on canonical distance.
- k = `k_value` = p without inversion, or 2p with inversion (λ_u < 0).
- The canonical distance c is the linearising arclength coordinate on a branch, 0 at the periodic point.
- **Every piece of branch bookkeeping is one rule:**
  - `advance_key(key, n)` says which branch M^n lands on.
  - `c ↦ c·β^{±n}` says where on it.
  - A full return to the same branch is k steps, so it scales c by β^k = |λ_u|^{num_branches}.

Inversion, period-k and several saddles are **no special case** below, except at the points marked ★.

## Stages

0. **Manifolds.** Input: a saddle and its eigendata. Output: a polyline for each grown stable and unstable branch, with c monotone along it. Rule: the fundamental segment [x, M^k x] is mapped forward (unstable) or backward (stable). Invariant: c is non-decreasing along the curve.
   - ★ Inversion: the two branches of one stability are the two halves of ONE smooth curve through z, and M^p swaps them.

1. **Trellis.** Input: finite pieces of the manifolds. Output: the finite crossing set N, each crossing = (unstable branch, c_u, stable branch, c_s, sign = sgn(u+ × s+)), plus the anchors.
   - Law: u×u, s×s and self-crossings cannot happen. When one shows up it is a numerical artifact, never a crossing. The index never pairs two segments of one stability (one rtree per stability, WP10); a transversal same-stability polyline crossing is only a resolution flag.
   - The rotation at every node is fixed by `sign` alone.
   - ★ Inversion: the anchor is ONE point with two stable and two unstable branches attached. Each branch has its own left and right side; the left of stable branch 0.0 is the same region as the right of 0.1. A bridge of u(o,b) leaving the anchor meets s(o,b) there (C7), whose sign is the anchor's own.

2. **Iterate table.** M(crossing) is another crossing, or absent.
   - Topology only looks crossings up in the table. Where a link is absent, it predicts with `advance_key` and β-scaling; it never calls the map.
   - Numerics may call the map to make the table exact.
   - Invariant: table links form a partial bijection. Compare c_u and c_s individually, never their product.
   - Topology compares crossings by identity and by registry floats; the only slack is the relative error of a scaled cdist (WP11, C15).

3. **Strong pip q₀** (`Strong_Pip_Algorithm.pdf`). q₀ is a strong pip iff no crossing, mapped onto q₀'s stable branch, lands in the open box (0, c_s(q₀)) × (0, c_u(q₀)).
   - Output: the pip orbit (one cut point per stable branch), the resonance zone, the fundamental segment F = (c_s(M^k q₀), c_s(q₀)] on the pip's branch, and **the trellis of the pip T(q₀) = (W^u[P, M^m q₀], W^s[P, q₀])** (D10, D12; `T(q0)`, `Trellis.cut_at_pip_orbit`): stable branch j is cut at c_s(M^j q₀), j < k, and the unstable branch of M^i q₀ at c_u(M^i q₀) for the k deepest registered i (m = the last forward iterate of q₀ the table holds). M⁻¹ maps its unstable half into itself and M its stable half, so it is a genuine trellis. Every later stage reads its crossing set N and its arcs; a crossing past a cut is a point of X, not of N, and the arc across it is one arc of T(q₀). A missing orbit point is placed by `image_cdist` and flagged.
   - ★ Several saddles: a trellis holds a crossing only when both its branches are selected, so a tangle alone holds no heteroclinic crossing (its arcs run through them). The selection's trellis, cut at every pip, is T(q₀_A) ∪ T(q₀_B): each branch at its own pip orbit, a heteroclinic crossing in N iff it is inside A's unstable cut and B's stable cut. M⁻¹ maps the union of the unstable halves into itself and M the union of the stable halves, so it is a genuine trellis (WP9).

4. **Pseudoneighbors** (`Pseudoneighbor_Algorithm.pdf`). N is T(q₀)'s (stage 3). Consecutive crossings x, x′ on a stable branch are a pair iff the unstable arc between them contains no point of X = ∪_n Mⁿ(N). A trajectory is extended only while its images are crossings of T(q₀).
   - Reference pairs are the pairs in F. Their +1..+(k−1) iterates on the other branches of a period-k orbit are recorded but get no hole (D15, 2026-10-09: holes propagate backward only, on every orbit).
   - A landing is compared on ONE polyline: the registered landing's own, else r's own against the pair's registered images there (WP8); β-scaling onto the pair is the fallback only when neither image is registered.
   - A missing pip cut is placed by `image_cdist` (a missing `+j` appearance costs nothing since D15).
   - ★ Heteroclinic (D7): a pair of heteroclinic crossings on B's stable branch can be a pair of B's tangle. Its unstable arc is A's, so X is read with A's branch cycle, k and β, over the crossings whose unstable side is A's. It needs both tangles selected: stages 4-6 run on the selection's trellis (WP9).

5. **Holes** (`Stable_Manifold_Partition_Algorithm`). Each reference pair puts one hole in its lobe. A reference bridge is mapped BACKWARD until it maps into a bridge already visited; each backward image lies inside one existing bridge, and that bridge gets the hole. Forward iterates (≥ k) get no hole.
   - A hole's side is its side of its OWN bridge, in the dynamical orientation, and M preserves it (det J > 0).
   - ★ A heteroclinic reference's backward chain runs along A's unstable cycle into A's bridges, so it holes A's partition too.
   - Invariants:
     - I1: one origin ⇒ one side.
     - I2: a bridge whose two ends are on one stable branch has the same row at both ends (its signs differ).

6. **Stable partition.** Each stable branch is cut at the endpoints of every hole bridge. An interval is open at an end iff a hole lies on that half. A point excluded on both sides is a singleton.
   - Everything here is sign arithmetic (D6):
     - `row(first) = left ⇔ sign > 0`, `row(second) = left ⇔ sign < 0`.
     - The half on a given side of the bridge at crossing x is `outward ⇔ (sign(x) > 0) == (side == left)`.
   - Geometry is only a cross-check, and it flags under-resolution.

7. **Bridge classes.** The class of a bridge is the unordered pair {X, Y} of elements its two ends sit against on their rows: its homotopy class relative to the holes, oriented anchor-outward.
   - If both ends are in one element, the bridge is a loop. It is folded FORWARD into the class of the bridge whose image contains it (D2).
   - A class is inert iff every image is a loop or inert.
   - The class set is finite and grouped by tangle.
   - The alphabet uses the minimal trellis's bridges only (D2), CLOSED under the word map (stage 8).

8. **Minimal trellis, dual graph, word.** The minimal trellis keeps the hole bridges plus the images of the active ones. The iterated partition is the homotopy partition cut by those images. The dual graph has faces as nodes, plus one stable node per edge side, unified only on the fundamental segment.
   - **The word of a class is the shortest walk** in the dual graph from the landing SPAN f(X) of one end to the landing span f(Y) of the other (D3): it leaves from and arrives at the stable edges the spans meet, since the hole bridges bounding X map onto the kept image bridges at the ends of f(X), which f(b) cannot cross (D13, C14). The walk always exists; an unreachable walk or a virtual class is a construction bug.
   - **The non-crossing rule (D13).** f is a homeomorphism and W^u never crosses W^u. A member bridge b with both ends on one stable branch and the stable segment σ_b between them bound a disk D_b (the lobe); f(b) ∪ f(σ_b) bounds f(D_b), which holds f(m) iff D_b holds m, for every **marker** m: a periodic point (its anchor vertex; f moves it to the next orbit point) or a hole (beside its OWN bridge on its bridge side; f moves hole n of an orbit to hole n + 1, and a reference hole's image sits beside the image chain of its bridge). The word is the shortest walk whose curve, closed by σ' along S between its end edges, holds exactly {f(m) : m ∈ D_b}, read combinatorially (a flood of the faces the curve cuts; any routing past a face's other boundary components) and looked for up to four crossings past the shortest walk. A walk crossing σ', or whose chords cross in a face, is never f(b) (f(b) is simple and meets f(σ_b) only at its ends). It applies to every class with a member lobe, tie or not (C13). **Walls** (D13's own words: f(b) cannot cross the images of the unstable arcs; C17): every class's arcs lie on known unstable branches (its members', or one map step past its namer's, D11); arcs on different unstable branches are disjoint and so are their images, so the decided walk of an earlier class on disjoint branches is a wall, and a tie loses every walk whose chord interleaves with a wall's on a face ring (only a walk of one face path is a wall or is judged by one). A class with no member lobe is no special case: a connecting class's tie is decided by the images of the heteroclinic bridges it is disjoint from (D14). A tie no rule decides stays flagged.
   - The iterated partition opens the stretches a hole's image lobe f(H) leaves empty: its base, and its chords, paired in chain order on the base branch and in stable order on any other branch (f(H) holds no anchor, so a branch walked from its anchor enters and leaves it in turn; C19).
   - Every pair of a walk is the class of one piece of f(B). A word naming a class outside the alphabet adds it, and that class is walked in turn (C9), until no word names a new pair (element pairs are finite). A pair no computed bridge spans is a class all the same (D11): its word is walked from its elements' landings, and it is flagged "no computed bridge yet", which never makes the words unreliable.
   - The walk may cross S only where an image can: at the images of its ends and in S \ f(S). S stops at the pip orbit (stage 3, D10), so that is the pip window, and the other k − 1 branches are walls.
   - At every crossing the walk reads the entry element, then the exit element. Disjoint consecutive pairs are classes or their inverses, and loops drop out.
   - Singletons cannot walk; their word is f(bridge) read off the registered image chain. Everywhere else that chain is the cross-check.
   - ★ Cross-side (D4) is a flag only for a pair on the SAME branch.

9. **Refinement and transition matrix.** A class whose symbol occurs with two or more distinct image pairs splits into a₁, a₂, …. The nodes are the active refined symbols, and the edge weights are token multiplicities. Inert classes are not nodes.

10. **Recompute (D1).** Topology is a pure function of the numerical state. Any growth or blast recomputes stages 1-9 from scratch. A blast's image curves are computed manifold: every recompute re-indexes them, so grow, blast and recompute commute (NUM-12, C20). Successive recomputes are compared by element identity (the bounding crossings of each element: hole-bridge endpoints, whose ids persist in a session; C10), and the class count may only go up: N is T(q₀)'s (D12), so a lap the pip orbit has not completed, a lambda-lemma copy deeper than M^m q₀, is not in N and cannot reshape the fundamental segment. Measured (E12, WP12, cutoff 1e-7): the seed count is non-decreasing and plateaus on k28, k32, p3, nested_outer, nested_inner and heteroclinic; inversion drops once, at blast 6, by a genuine switch (BCG-8, DEC-25); at the coarse cutoff 1e-5 several series drop (NUM-5: near-parallel copies ordered wrongly across the window), so 1e-7 is the recommended cutoff (DEC-14). While the hole set stays the same, so does the alphabet closed under the word map (D11).

★ **Several saddles.** A heteroclinic crossing is a node of both tangles, and of neither alone: the topology of several saddles is computed on ONE trellis, the selection's (stage 3). Nesting of tangles is the rotation system's outer face of each component (the face holding due east at its easternmost vertex, never an area, WP10) plus ONE containment test per component, and a zone's blast iterates the bridges whose innermost zone it is (C4). A heteroclinic bridge is a bridge like any other (D14): it is classed, walked and a wall by the same rules; a connecting class exists only where a computed bridge connects the tangles (C18), and its images are D11 classes until growth computes them (at p3 17u they are). Nothing else changes.

## Where the code diverged on 2026-10-07

From the audit; the details and evidence are in `open_issues.md`.

| Stage | Divergence | Issue | WP |
|---|---|---|---|
| 1 | Four coincident anchors at an inversion point | INV-1 | fixed (WP3) |
| 2 | Iterate table filled by a 5% cdist heuristic, never after a blast | INT-6, BCG-3 | fixed (WP2) |
| 4 | Lap factor λ_u instead of β^k | INV-2 | fixed (WP2) |
| 4 | Missing table links leave branches without holes | PK-3 | accepted (D9): warned (WP5), nothing invented |
| 5-6 | Holes and rows read from polyline geometry, not signs | INT-3 | fixed (WP4) |
| 7 | Every bridge classed; loops folded backward | BCG-1, BCG-2 | fixed (WP6) |
| 8 | Unreachable walks and virtual classes treated as states | PK-1, INT-9 | no longer states (WP7); a walk pair no computed bridge spans is a flagged class (D11, WP8b) |
| 3, 8 | The partition and dual graph do not cut S at the pip orbit (unified window ≠ S \ f(S) when S runs past it) | PK-5 | fixed (WP8b, D10): the trellis is cut at the pip orbit |
| 10 | Element identity is a positional index | BCG-7 | fixed (WP8): letters by bounding crossings |
| 4 | Heteroclinic pairs rejected (A's arc tested against B's branch cycle); stages 4-6 per tangle | NEST-2 | fixed (WP9): the arc's own cycle, stages on the selection's trellis |
| — | An outer blast iterated the bridges inside an inner zone | NEST-4 | fixed (WP9, C4): innermost zone |
| 4 | Class-count decreases at depth: references lost to witnesses inside the cdist placement error bar (noise not ruled out; one clear lambda-lemma copy) | BCG-5 | fixed (WP8b, D12): every witness lay deeper than the pip's own orbit; N is cut there |
| 8 | Equal shortest walks picked by order; the walk started on any edge of the iterated element | DEC-23 | WP9b: span ends + the marker law (D13); WP9c: walls (D13, D14) |
| 8 | Chords of a lobe image on a foreign stable branch left unpaired ("abuts no hole") | PK-4 | fixed (WP10b, C19) |
| 10 | A recompute after a blast dropped every crossing on its image curves | NUM-12 | fixed (WP10b) |
| 8 | Registered evidence contradicts the walk at depth through cdist noise (an image order the map cannot produce) | INT-11 | fixed (WP12, C21: an order flip gives no verdict); kind 2 (no registered preimage to prove it) is the resolution floor, NUM-3/4, deferred |
| 4, 10 | The inversion series loses classes at blast 6: a genuine reference switch whose `+1` pair is never spanned (the reference arc lies outside the zone) | BCG-8 | open: DEC-25 |

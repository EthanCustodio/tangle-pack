# Error analysis of the k=2.5 entropy sweeps (2026-10-29)

`errors.png` (+ .svg, .fig.pickle) has twelve panels showing what limits the pip depth, plotted against `area_cutoff`. It is drawn by `scripts/henon_k25_errors_figures.py` from `data/`.

The data come from `scripts/henon_k25_errors.py`. That script reruns every pruned sweep: `min_separation` 1e-5, one blast at a time, until the frontier is exhausted. Instead of the topology, it records numerical errors at every blast. Every sweep reproduces the original run's pip depth exactly.

## Findings

1. **The registered pip orbit stops being the real orbit a few levels before the reported depth** (panels a, b, l and e, hollow markers).
   - The real orbit contracts toward P by β at each step: 5.56 for p1, and λ per period for p3. The registered orbit does so down to about |f^j − P| ≈ 1e-9 … 1e-10. That is the image curves' chord deviation at that cutoff (panel f).
   - Past that point the registered crossings lie closer to P than the polylines can place them. Their positions, and so their order on the stable branch that cuts T(q0), are noise. Their resolution ratios are still around 1e4, because the ratio never compares a crossing's distance to the anchor against h.
   - Trustworthy depth: p1 7–13, and p3 16–33 (1e-7 to 1e-14), against 9–22 and 17–40 reported.

2. **Most of what looked like progress at fine cutoffs comes from beyond the trustworthy depth.**

   | case | trustworthy h | h from beyond the trustworthy depth |
   |---|---|---|
   | p1 alone | 0.5280489 (x³ − x² − 2) at every cutoff | 0.5281139 and 0.5281272 at 1e-13 |
   | p3 alone | 0.1877778 (1e-13 and 1e-14, depth ≤ 33) | 0.1885180 and 0.1892423 |
   | nested | 0.5281517, reached at a trustworthy depth from 1e-8 on (p1 9, p3 16) | — |

   The nested value stands.

3. **The depth stalls for two reasons** (panels c, d).
   - **1e-7, and p3 at 1e-8:** the image curve reaches f^{m+1}(q0), but the crossing's resolution ratio is 0.7–2.1, below 3. The unstable chord deviation (h ≈ 1.5·A/L on long, nearly straight segments) is larger than the short stable segment near P.
   - **1e-9 and finer:** the bridge that ends at f^m(q0) is never iterated. `classify_bridge` tests its index-middle node against the frozen zone polygon. The bridge is a long arc (15.6 to 53.5 in length) whose middle lies outside the zone, so `zone=None`. Confirmed on:
     - p1 at 1e-10, 1e-12 and 1e-13;
     - p3 at 1e-12.

     At 1e-7, the same 15.6-long bridge was classified inside the zone and iterated.

4. **The canonical-distance precision runs out at the deepest blasts** (panel h). The smallest relative cdist step reaches float64 eps at p1 1e-12 and 1e-13 (blasts 13–20).

5. **Same-stability polyline crossings jump at the deep blasts** (panel i). p1 at 1e-13 gives 37k–420k per blast at blasts 14–21. These are near-coincident unstable curves close to P.

6. **Prune comparison** (panel j). Against no `min_separation`:
   - the unpruned series give more classes and a different h at p1 1e-11 (depth 14) and p3 1e-12 (depth 33);
   - both differences lie beyond the trustworthy depth (12 and 29), so this gives no evidence yet that the prune hides real topology;
   - at 1e-9 the two series are identical.

Rerun: `PYTHONPATH=src env/bin/python scripts/henon_k25_errors.py --cases p1 --cutoffs 1e-7 1e-8` (it redraws at the end), or `scripts/henon_k25_errors_figures.py` alone.

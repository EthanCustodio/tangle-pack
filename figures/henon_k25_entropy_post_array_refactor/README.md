# Hénon k=2.5 entropy, rerun after the numerics arrays campaign (main 5ed9940, 2026-10-28/29)

Made by `scripts/henon_k25_entropy_post_array_refactor.py` (the k=2.5 entropy scripts pointed at this folder).
Same figures as `figures/henon_k25_entropy/`. Pruned sweeps use `min_separation` 1e-5 as before and run until
the frontier is exhausted (the first blast that keeps no new interior bridge in any zone), or until they hit
their time or memory cap.
`components/` is at 1e-9, as before. `components/finest/` reads each case at its finest exhausted sweep.
The unpruned series (`components/unpruned/`) and the length growth stop on memory or time.
The p3 length-growth curve is copied from the pre-refactor run. That curve is NumPy-only past the build,
its initial curve is identical, and p1's 22 recomputed curve steps match the old run bit for bit.

| sweep | blasts | stop | pip depth | h(M1) | active | peak GB | s |
|---|---|---|---|---|---|---|---|
| nested 1e-07 | 91 | frontier exhausted at blast 91 | {'p3': 17, 'p1': 9} | 0.5281517 | 59 | 1.3 | 304 |
| nested 1e-08 | 60 | frontier exhausted at blast 60 | {'p3': 20, 'p1': 11} | 0.5281517 | 77 | 1.5 | 209 |
| nested 1e-09 | 60 | frontier exhausted at blast 60 | {'p3': 27, 'p1': 12} | 0.5281517 | 117 | 2.5 | 256 |
| nested 1e-10 | 60 | frontier exhausted at blast 60 | {'p3': 30, 'p1': 13} | 0.5281517 | 136 | 4.3 | 371 |
| nested 1e-11 | 60 | frontier exhausted at blast 60 | {'p3': 31, 'p1': 14} | 0.5281517 | 136 | 8.4 | 509 |
| nested 1e-12 | 60 | frontier exhausted at blast 60 | {'p3': 33, 'p1': 15} | 0.5281517 | 149 | 17.5 | 991 |
| p1 1e-07 | 30 | frontier exhausted at blast 30 | {'p1': 9} | 0.5280489 | 19 | 1.2 | 122 |
| p1 1e-08 | 27 | frontier exhausted at blast 27 | {'p1': 11} | 0.5280489 | 19 | 1.5 | 92 |
| p1 1e-09 | 27 | frontier exhausted at blast 27 | {'p1': 12} | 0.5280489 | 19 | 2.2 | 131 |
| p1 1e-10 | 27 | frontier exhausted at blast 27 | {'p1': 13} | 0.5280489 | 19 | 4.0 | 221 |
| p1 1e-11 | 27 | frontier exhausted at blast 27 | {'p1': 14} | 0.5280489 | 19 | 7.8 | 355 |
| p1 1e-12 | 27 | frontier exhausted at blast 27 | {'p1': 15} | 0.5280489 | 44 | 16.1 | 522 |
| p1 1e-13 | 27 | frontier exhausted at blast 27 | {'p1': 22} | 0.5281272 | 66 | 34.8 | 1230 |
| p1 1e-14 | 10 | memory limit (50 GB): blast 11 predicted to peak at 51.3 GB | {'p1': 11} | 0.5280489 | 19 | 37.5 | 287 |
| p3 1e-07 | 91 | frontier exhausted at blast 91 | {'p3': 17} | 0.1535942 | 26 | 0.7 | 153 |
| p3 1e-08 | 60 | frontier exhausted at blast 60 | {'p3': 20} | 0.1773338 | 44 | 0.6 | 63 |
| p3 1e-09 | 60 | frontier exhausted at blast 60 | {'p3': 27} | 0.1859099 | 84 | 0.8 | 84 |
| p3 1e-10 | 60 | frontier exhausted at blast 60 | {'p3': 30} | 0.1877778 | 103 | 1.1 | 137 |
| p3 1e-11 | 60 | frontier exhausted at blast 60 | {'p3': 31} | 0.1877778 | 103 | 2.1 | 187 |
| p3 1e-12 | 60 | frontier exhausted at blast 60 | {'p3': 33} | 0.1877778 | 103 | 4.1 | 318 |
| p3 1e-13 | 60 | frontier exhausted at blast 60 | {'p3': 37} | 0.1885180 | 127 | 8.3 | 588 |
| p3 1e-14 | 60 | frontier exhausted at blast 60 | {'p3': 40} | 0.1892423 | 129 | 17.3 | 1072 |
| p3 1e-15 | 33 | memory limit (30 GB): blast 34 predicted to peak at 31.2 GB | {'p3': 41} | 0.1892423 | 129 | 30.2 | 1261 |
| unpruned nested 1e-09 | 17 | memory limit (40 GB): blast 18 predicted to peak at 58.8 GB | {'p3': 25, 'p1': 12} | 0.5281517 | 86 | 25.7 | 1890 |
| unpruned p1 1e-09 | 18 | memory limit (45 GB): blast 19 predicted to peak at 56.0 GB | {'p1': 12} | 0.5280489 | 19 | 25.6 | 1789 |
| unpruned p3 1e-09 | 45 | time limit (6000 s) before blast 46 | {'p3': 27} | 0.1859099 | 84 | 7.5 | 7985 |

## h against the pip depth

`convergence_vs_depth` plots h against the pip depth m, one line per cutoff, with the y axis zoomed to the plateau. Blasts that leave the depth unchanged collapse onto one point. Panels: (a) p1 alone, (b) p3 alone, (c) nested against the p1 depth, (d) nested against the p3 depth. Each dotted line is the other case's final h at its finest sweep that ran to completion. Nested settles at 0.5281517 at p1 depth 9 / p3 depth 15. p1 alone leaves 0.5280489 only at 1e-13, where depths 18 and 21 bring it to 0.5281272. `entropy_vs_depth` shows the same points on the full 0 to ln 2 axis.

## Caveat from the error analysis (`errors/`)

The deepest pip-orbit points of every sweep lie closer to P than the curves' chord deviation, so their positions are noise. Only the depth at which the orbit still contracts by β is trustworthy:
- **p1 alone:** 0.5280489 is the trustworthy value. The 1e-13 rise to 0.5281272 comes from depths 15–22, which are beyond it.
- **p3 alone:** 0.1877778 is trustworthy. The values above it come from beyond the trustworthy depth.
- **nested:** 0.5281517 stands.

See `errors/README.md`.

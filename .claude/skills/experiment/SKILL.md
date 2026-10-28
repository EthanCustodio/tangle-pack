---
name: experiment
description: Write and run a tanglepack experiment or data-generation script (sweeps, measurements, series) in scripts/, writing figures and data, and summarize the result.
argument-hint: <the question the experiment answers>
disable-model-invocation: true
---

# /experiment $ARGUMENTS

Stays on the current branch; commit only when asked. No plan mode unless the run is expensive
(over ~10 minutes or ~10 GB): then state the parameter grid, cost estimate and outputs and ask.

**The script** `scripts/<name>.py` (look at `scripts/henon_period3_k_sweep.py` for the shape):

- A module docstring that states the question, the method, what is (and is not) computed, the
  outputs, and the growth/budget policy. This is the record of the experiment; keep it accurate.
- `argparse` for every knob, with the defaults that answer the question. Frozen parameters are
  passed explicitly: always an explicit `area_cutoff`; build cases with
  `tanglepack.examples.henon_cases` / `henon.py` rather than copying recipes.
- Outputs to `figures/<name>/` (figures) or `data/<name>/` (large arrays, npz), plus a small CSV of
  the measured numbers and a `summary.png` when there is a trend to see. Use `--chunk`/`--merge`
  only when a sweep is long enough to parallelize.
- `logging` for progress, not print, except the final summary line.
- Run it as `PYTHONPATH=src MPLBACKEND=Agg env/bin/python scripts/<name>.py ...`.

**Run** a quick smoke setting first (a coarse cutoff or a few frames), then the real run in the
background. Read the CSV and look at the summary figure.

**Report** the answer to the question in a few lines with the numbers, the output paths, and
anything surprising. Topology changes under deeper growth or blasting are discovery, not failure.
If the result settles a question the author cares about, offer to record it (decisions.md entry
or memory).

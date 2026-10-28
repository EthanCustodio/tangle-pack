---
name: figure
description: Design and generate a tanglepack figure (trellis, partitions, cartoons, sweeps, report plots), iterating on the rendered image until it reads well.
argument-hint: <what the figure should show, and for whom>
disable-model-invocation: true
---

# /figure $ARGUMENTS

Stays on the current branch; commit only when asked.

1. **Reuse first.** Topology drawing lives in `src/tanglepack/topology/plotting.py` (bridges by
   class, cartoons, walks, itinerary tables, legends) and `loom/cartoon.py`; panel layout in
   `scripts/symbolic_figures.py` (`draw_panels`). Use them; add a plotting routine to the library
   only if the figure will be reused (then it gets a smoke test in `tests/plotting/`).
2. **Design.** Load the `dataviz` skill for colour, marks and layout choices on charts. For
   manifold pictures keep the library's conventions: `class_colors`, inert dashed, anchor on the
   right in cartoons. Decide the one thing the reader must see and make it the most salient mark.
3. **Generate** from a script (an existing one with a new flag, or `scripts/<name>.py`), never
   from an interactive session, so the figure can be regenerated. Save `.png`, `.svg` and
   `.fig.pickle` into `figures/<name>/` (freeze networkx curved labels before pickling, as
   `symbolic_figures.py` does).
4. **Look.** Read the PNG. Check: labels legible and not overlapping, nothing clipped, the key
   feature visible at a glance, the legend matches the marks, colours distinguishable. Fix and
   regenerate; usually 2-3 rounds.
5. **Show** the final PNG path(s) and a line on what each panel shows. Offer
   `scripts/show_figure.py <pickle>` for a live look.

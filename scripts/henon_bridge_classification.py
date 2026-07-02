"""
Classify every bridge of the nested Henon tangle by resonance zone.

Builds the exact same manifolds as ``henon_resonance_zone.py`` (outer period-1 fixed
point + inner period-3 orbit), defines one resonance zone per fixed point, then walks
every bridge and reports which zone it falls in. The summary is a table of bridge
counts: one row per zone (plus an "outside all zones" row), one column per fixed
point, so you can read off how many of each fixed point's bridges sit in each zone.

A bridge is attributed to the innermost zone whose boundary contains its midpoint, so
the inner period-3 zone's own boundary bridges count toward that inner zone.

Run:  PYTHONPATH=src python scripts/henon_bridge_classification.py
"""

import logging

logging.basicConfig(level=logging.WARNING)

import numpy as np

from tanglepack import TangleSession

import matplotlib.pyplot as plt


def henon_map(point):
    k, b = 2, 1
    x, y = point
    return np.array([y - k + x**2, -b * x])


def henon_map_inverse(point):
    k, b = 2, 1
    x, y = point
    return np.array([-y / b, x + k - (y**2) / (b**2)])


def henon_jacobian(point):
    k, b = 2, 1
    x, y = point
    return np.array([[2 * x, 1], [-b, 0]])


# --- build the nested tangle (identical to henon_resonance_zone.py) -------------
session = TangleSession(henon_map, henon_map_inverse, henon_jacobian)

# inner period-3 orbit
session.workbench._man_machine.area_cutoff = 1e-7
fp3 = session.construct_fixed_point([[0, 1], [-1, 0], [-1, 1]])
session.orient_eigenvectors(
    fp3, {"unstable": np.array([0, -1]), "stable": np.array([-1, -1])}
)
session.initialize_both_manifolds(fp3)
session.grow_n_times(fp3, "unstable", num_iterations=13)
session.grow_n_times(fp3, "stable", num_iterations=9)

# outer period-1 fixed point
session.workbench._man_machine.area_cutoff = 1e-4
fp1 = session.construct_fixed_point([4, -4])
session.orient_eigenvectors(
    fp1, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
)
session.initialize_both_manifolds(fp1)
session.grow_n_times(fp1, "unstable", num_iterations=11)
session.grow_until_turnaround(fp1, "stable")

# co-compute the nested tangle and cut bridges
session.compute_intersections([fp3, fp1])
session.trim_stable_manifolds(fp3)
session.trim_stable_manifolds(fp1)
session.create_bridges(fp3)
session.create_bridges(fp1)
session.infer_iterate_table()

# pick a strong pip on each tangle
T1 = session.trellis(fp1)
T1.classify_strong_pips()
T3 = session.trellis(fp3)
T3.classify_strong_pips()
# The pip previously pinned here (id 7) predates the fsolve solver recalibration;
# intersection ids shifted when the growth counts were rebalanced. Use the default
# pip classify_strong_pips() selects (smallest unstable cdist); re-pin a specific
# candidate with T1.set_strong_pip(id) if a different zone boundary is wanted.

# one resonance zone per fixed point (single final recompute keeps the nesting)
session.add_resonance_zones([T1.strong_pip, T3.strong_pip])


# --- classify every bridge ------------------------------------------------------
classification = session.classify_bridges()  # {bridge: ResonanceZone | None}


# --- report ---------------------------------------------------------------------
def fp_label(fp):
    """Short readable name for a fixed point, e.g. 'fp(period 3)'."""
    return f"fp(period {fp.period})"


# Stable column order: the workbench's fixed points.
fixed_points = list(session.fixed_points)
fp_cols = [fp_label(fp) for fp in fixed_points]

# Stable row order: each defined zone (insertion order), then the catch-all.
zone_rows = []  # (row_label, zone_key) ; zone_key is None for the "outside" row
for (fp, branch_index), rz in session.resonance_zones.items():
    label = f"Zone: {fp_label(fp)} branch {branch_index}"
    zone_rows.append((label, rz.key))
zone_rows.append(("Outside all zones", None))

# counts[zone_key][fp] -> int
counts = {key: {fp: 0 for fp in fixed_points} for _, key in zone_rows}
for bridge, zone in classification.items():
    key = zone.key if zone is not None else None
    counts[key][bridge.fixed_point] += 1

# Pretty-print an aligned table.
row_header_w = max(len(label) for label, _ in zone_rows + [("Total", None)])
col_w = max(8, max(len(c) for c in fp_cols) + 1)

print("\nNested Henon tangle — bridge classification by resonance zone")
print("=" * 70)
print(f"Total bridges classified: {len(classification)}")
print()

# header
header = " " * row_header_w + "  " + "".join(c.rjust(col_w) for c in fp_cols)
header += "Total".rjust(col_w)
print(header)
print("-" * len(header))

# zone rows
col_totals = {fp: 0 for fp in fixed_points}
for label, key in zone_rows:
    row_counts = counts[key]
    row_total = sum(row_counts.values())
    line = label.ljust(row_header_w) + "  "
    line += "".join(str(row_counts[fp]).rjust(col_w) for fp in fixed_points)
    line += str(row_total).rjust(col_w)
    print(line)
    for fp in fixed_points:
        col_totals[fp] += row_counts[fp]

# totals row
print("-" * len(header))
total_line = "Total".ljust(row_header_w) + "  "
total_line += "".join(str(col_totals[fp]).rjust(col_w) for fp in fixed_points)
total_line += str(sum(col_totals.values())).rjust(col_w)
print(total_line)
print()


plt.figure(figsize=(8, 8))
session.plot_resonance_zones(alpha=0.3)

for fp in (fp3, fp1):
    session.workbench.plot_tangle(fp, "stable", color="r", linewidth=0.8)
    session.workbench.plot_tangle(fp, "unstable", color="b", linewidth=0.8)

# Rebuild trellises (the recompute rebuilt the registry) before plotting pips.
T1 = session.trellis(fp1, rebuild=True)
T1.classify_strong_pips()
T3 = session.trellis(fp3, rebuild=True)
T3.classify_strong_pips()
T1.plot_strong_pip(s=40)
T1.plot_strong_pip_candidates()
T3.plot_strong_pip(s=40)
T3.plot_strong_pip_candidates()

session.workbench.plot_intersections(fp1, show_ids=True)
session.workbench.plot_all_bridges()

plt.xlim([-8, 8])
plt.ylim([-8, 8])
plt.title("Nested resonance zones (outer period-1, inner period-3)")
plt.show()

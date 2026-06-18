import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)

import tanglepack, numpy as np
import matplotlib.pyplot as plt
from tanglepack import Trellis


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


wb = tanglepack.TangleWorkbench(henon_map, henon_map_inverse, henon_jacobian)

wb._man_machine.area_cutoff = 1e-7

# Period-3 inner fixed point
fp3 = wb.construct_fixed_point([[0, 1], [-1, 0], [-1, 1]])
# wb.orient_eigenvectors(fp3, {"unstable": np.array([0, 1]), "stable": np.array([1, 1])})
wb.orient_eigenvectors(
    fp3, {"unstable": np.array([0, -1]), "stable": np.array([-1, -1])}
)
wb.initialize_both_manifolds(fp3)
# wb.grow_n_times(fp3, "unstable", num_iterations=12)
# wb.grow_n_times(fp3, "stable", num_iterations=12)
wb.grow_n_times(fp3, "unstable", num_iterations=10)
wb.grow_n_times(fp3, "stable", num_iterations=6)


wb._man_machine.area_cutoff = 1e-4
fp1 = wb.construct_fixed_point([4, -4])
wb.orient_eigenvectors(fp1, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])})
wb.initialize_both_manifolds(fp1)
wb.grow_n_times(fp1, "unstable", num_iterations=7)
wb.grow_until_turnaround(fp1, "stable")


# Nested tangles: co-index BOTH fixed points into one Tangle/registry so the inner
# (period-3) and outer (period-1) crossings coexist. Computing them separately would
# wipe the first tangle, since each compute_intersections() rebuilds the registry.
wb.compute_intersections([fp3, fp1])
wb.trim_stable_manifolds(fp3)
wb.trim_stable_manifolds(fp1)
wb.create_bridges(fp3)
wb.create_bridges(fp1)
wb.infer_iterate_table()


T1 = Trellis.from_workbench(wb, fp1)
T1.classify_strong_pips()
print(f"Strong pip candidates T1: {T1.strong_pip_candidates}")
print(f"Strong pip {T1.strong_pip}")

T3 = Trellis.from_workbench(wb, fp3)

T3.classify_strong_pips()
print(f"Strong pip candidates: {T3.strong_pip_candidates}")
print(f"Strong pip {T3.strong_pip}")

plt.figure()
wb.plot_tangle(fp3, "stable", color="r", linewidth=1)
wb.plot_tangle(fp3, "unstable", color="b", linewidth=1)

wb.plot_intersections(fp3)
# wb.plot_all_bridges()
T3.plot_strong_pip_candidates()
T3.plot_strong_pip(s=12)

wb.plot_tangle(fp1, "stable", color="r", linewidth=1)
wb.plot_tangle(fp1, "unstable", color="b", linewidth=1)

wb.plot_intersections(fp1, show_ids=True)
# wb.plot_all_bridges()
T1.plot_strong_pip_candidates()
T1.plot_strong_pip(s=12)


plt.xlim([-15, 15])
plt.ylim([-15, 15])

plt.show()

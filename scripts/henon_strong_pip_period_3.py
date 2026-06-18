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
wb.grow_n_times(fp3, "unstable", num_iterations=8)
wb.grow_n_times(fp3, "stable", num_iterations=6)


wb.compute_intersections(fp3)
wb.trim_stable_manifolds(fp3)
bridges = wb.create_bridges(fp3)
wb.infer_iterate_table()


T = Trellis.from_workbench(wb, fp3)

T.classify_strong_pips()
print(f"Strong pip candidates: {T.strong_pip_candidates}")
print(f"Strong pip {T.strong_pip}")

plt.figure()
wb.plot_tangle(fp3, "stable", color="r", linewidth=1)
wb.plot_tangle(fp3, "unstable", color="b", linewidth=1)

wb.plot_intersections(fp3)
# wb.plot_all_bridges()
T.plot_strong_pip_candidates()
T.plot_strong_pip(s=12)

plt.xlim([-1.2, 0])
plt.ylim([-0.1, 1.2])

plt.show()

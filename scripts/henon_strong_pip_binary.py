import tanglepack, numpy as np
import matplotlib.pyplot as plt
from tanglepack import Trellis


def henon_map(point):
    k, b = 10, 1
    x, y = point
    return np.array([y - k + x**2, -b * x])


def henon_map_inverse(point):
    k, b = 10, 1
    x, y = point
    return np.array([-y / b, x + k - (y**2) / (b**2)])


wb = tanglepack.TangleWorkbench(henon_map, henon_map_inverse)
fp = wb.construct_fixed_point([4, -4])
wb.orient_eigenvectors(fp, {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])})
wb.initialize_both_manifolds(fp)
wb.grow_n_times(fp, "unstable", num_iterations=6)
wb.grow_until_turnaround(fp, "stable")

wb.compute_intersections(fp)
wb.trim_stable_manifolds(fp)
bridges = wb.create_bridges(fp)

new_bridges = wb.iterate_bridge(bridges[2])
wb.iterate_bridge(new_bridges[0])

wb.infer_iterate_table()


T = Trellis.from_workbench(wb, fp)

T.classify_strong_pips()
print(f"Strong pip candidates: {T.strong_pip_candidates}")
print(f"Strong pip {T.strong_pip}")

plt.figure()
wb.plot_tangle(fp, "stable", color="r")
wb.plot_intersections(fp)
wb.plot_all_bridges()
T.plot_strong_pip_candidates()
T.plot_strong_pip()

plt.xlim([-15, 15])
plt.ylim([-15, 15])

plt.show()

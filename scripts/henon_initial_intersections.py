import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)

import tanglepack
import numpy as np
import matplotlib.pyplot as plt


def henon_map(point):
    """
    Defines the henon map for binary horshoe parameters to test basic functionality
    """

    k, b = (10, 1)

    x = point[0]
    y = point[1]

    return np.array([y - k + x**2, -b * x])


def henon_map_inverse(point):
    """Defines the inverse henon map for"""

    k, b = (10, 1)

    x = point[0]
    y = point[1]

    return np.array([-y / b, x + k - (y**2) / (b**2)])


henon = tanglepack.DynamicalSystem(henon_map, henon_map_inverse)
fp_solver = tanglepack.FixedPointSolver(henon)
man_maker = tanglepack.ManifoldInitializer(henon)
man_machine = tanglepack.ManifoldMachine(henon)

tangle = tanglepack.Tangle()

initial_guess = [4, -4]

fixed_point = fp_solver.construct_fixed_point(initial_guess, 2)

print(f"The fixed point is: {fixed_point.coordinates[0]}")

approx_dirs = {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}

man_maker.orient_manifolds(fixed_point, approx_dirs)

# initial_unstable_segment = man_maker.get_initial_fundamental_segment(
#     fixed_point, 0, 0, "unstable"
# )
initial_stable_segment = man_maker.construct_kevin_way(fixed_point, "stable")

initial_unstable_segment = man_maker.construct_kevin_way(fixed_point, "unstable")

unstable_manifold = initial_unstable_segment[(0, 0)]
stable_manifold = initial_stable_segment[(0, 0)]

num_iterations = 6

# grow unstable manfold
man_machine.grow_x_times(fixed_point, "unstable", num_iterations)
unstable_manifold._find_tail()

# grow stable manifold
man_machine.grow_x_times(fixed_point, "stable", num_iterations)
stable_manifold._find_tail()


def intersections():
    tangle.add_manifold(unstable_manifold)
    tangle.add_manifold(stable_manifold)

    unstable_ids = tangle._manifold_segs[unstable_manifold]
    stable_ids = tangle._manifold_segs[stable_manifold]

    hits = []  # (segA, segB) pairs

    for sid in unstable_ids:
        segA = tangle._seg_lookup[sid]
        for segB in tangle.intersections_for_segment(segA):
            hits.append((segA, segB))

    # print(f"Found {len(hits)} intersections")
    # print(f"Intersections: {tangle._intersecting_segments}")


intersections()
tangle.populate_intersection_dict()
print(f"Intersections: {tangle._intersecting_coords.values()}")

fig = plt.figure()

show_points = False
unstable_manifold.plot("blue", show_points=show_points)
stable_manifold.plot("red", show_points=show_points)
plt.scatter(*fixed_point.coordinates[0], c="k", s=7)

for point in tangle._intersecting_coords.values():
    plt.scatter(*point, c="k", s=7, zorder=10)

plt.xlim([-15, 15])
plt.ylim([-15, 15])

# fig.savefig(
#     "tangle_plot.png",  # pdf/svg/eps/etc. all work
#     dpi=300,  # print-quality resolution
#     bbox_inches="tight",
# )  # trim extra whitespace

plt.show()

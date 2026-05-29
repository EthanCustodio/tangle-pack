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

    k, b = (2, 1)

    x = point[0]
    y = point[1]

    return np.array([y - k + x**2, -b * x])


def henon_map_inverse(point):
    """Defines the inverse henon map for"""

    k, b = (2, 1)

    x = point[0]
    y = point[1]

    return np.array([-y / b, x + k - (y**2) / (b**2)])


def henon_jacobian(point):
    """defines the jacobian for the henon map"""

    k, b = (2, 1)

    x = point[0]
    y = point[1]

    return np.array([[2 * x, 1], [-b, 0]])


henon = tanglepack.DynamicalSystem(henon_map, henon_map_inverse, henon_jacobian)
fp_solver = tanglepack.FixedPointSolver(henon)
man_maker = tanglepack.ManifoldInitializer(henon)
man_machine = tanglepack.ManifoldMachine(henon)

tangle = tanglepack.Tangle()

initial_guess = [[0, 1], [-1, 0], [-1, 1]]
initial_guess_zero = [4, -4]
# initial_guess = [6.104, 0]

fixed_point = fp_solver.construct_fixed_point(initial_guess, 2)
fixed_point_zero = fp_solver.construct_fixed_point(initial_guess_zero, 2)


print(f"fixed point: {fixed_point.coordinates}")

orbit_index = 0

# initial_unstable_segment = man_maker.get_initial_fundamental_segment(
#     fixed_point, orbit_index, 0, "unstable"
# )
# initial_stable_segment = man_maker.get_initial_fundamental_segment(
#     fixed_point, orbit_index, 0, "stable"
# )

approx_dir = {"unstable": [0, 1], "stable": [1, 1]}

man_maker.orient_manifolds(fixed_point, approx_dir)

approx_dir = {"unstable": [-1, 0], "stable": [0, 1]}

man_maker.orient_manifolds(fixed_point_zero, approx_dir)

# initial_unstable_segments = man_maker.get_all_initial_segments(fixed_point, "unstable")
# initial_stable_segments = man_maker.get_all_initial_segments(fixed_point, "stable")

initial_unstable_segments = man_maker.construct_kevin_way(fixed_point, "unstable")
initial_stable_segments = man_maker.construct_kevin_way(fixed_point, "stable")

initial_zero_unstable = man_maker.construct_kevin_way(fixed_point_zero, "unstable")
initial_zero_stable = man_maker.construct_kevin_way(fixed_point_zero, "stable")

print(f"points! {initial_unstable_segments[(0, 0)].get_point_array()}")

# unstable_manifold = initial_unstable_segments[(0, 0)]
# stable_manifold = initial_stable_segments[(0, 0)]

# grow unstable manifold
num_iterations = 7
man_machine.grow_x_times(fixed_point, "unstable", num_iterations)
man_machine.grow_x_times(fixed_point, "stable", num_iterations)

# for manifold in initial_unstable_segments:
#     manifold._find_tail()

# unstable_manifold._find_tail()

# grow stable manifold
# num_iterations = 7
num_iterations = 7
man_machine.grow_x_times(fixed_point_zero, "stable", num_iterations)
man_machine.grow_x_times(fixed_point_zero, "unstable", num_iterations)


# approx_dir = {"unstable": [0, -1], "stable": [-1, -1]}

# man_maker.orient_manifolds(fixed_point, approx_dir)

# initial_unstable_segments_new = man_maker.construct_kevin_way(fixed_point, "unstable")
# initial_stable_segments_new = man_maker.construct_kevin_way(fixed_point, "stable")

# num_iterations = 7
# man_machine.grow_x_times(fixed_point, "unstable", num_iterations)
# man_machine.grow_x_times(fixed_point, "stable", num_iterations)

# plt.figure()
# initial_unstable_segment.plot(show_points=True)
# plt.show()

# unstable_manifold = initial_unstable_segments[(0, 0)]
# unstable_manifold._find_tail()

# stable_manifold = initial_stable_segments[(0, 0)]
# stable_manifold._find_tail()


fig = plt.figure()

unstable_manifold = initial_zero_unstable[(0, 0)]
unstable_manifold._find_tail()
stable_manifold = initial_zero_stable[(0, 0)]
stable_manifold._find_tail()

unstable_manifold.plot(color="blue")
stable_manifold.plot(color="red")
plt.scatter(*fixed_point_zero.coordinates[0], c="k", s=20)

for i in range(fixed_point.period):

    unstable_manifold = initial_unstable_segments[(i, 0)]
    unstable_manifold._find_tail()

    stable_manifold = initial_stable_segments[(i, 0)]
    stable_manifold._find_tail()

    unstable_manifold.plot(color="blue")
    stable_manifold.plot(color="red")
    plt.scatter(*fixed_point.coordinates[i], c="k", s=12)

# unstable_manifold.plot(show_points=True, color="blue")
# stable_manifold.plot(show_points=True, color="red")
# plt.scatter(*fixed_point.coordinates[0], c="k", s=7)

# for point in tangle._intersecting_coords.values():
#     plt.scatter(*point, c="k", s=7, zorder=10)


plt.xlim([-6, 6])
plt.ylim([-6, 6])

plt.title("k=2, b=1 Period 3 Nested Tangle")

# fig.savefig(
#     "period_3_tangle_plot.png",  # pdf/svg/eps/etc. all work
#     dpi=600,  # print-quality resolution
#     bbox_inches="tight",
# )  # trim extra whitespace

plt.show()

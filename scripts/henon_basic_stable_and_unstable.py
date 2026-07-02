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

initial_guess = [4, -4]

fixed_point = fp_solver.construct_fixed_point(initial_guess, 1)

print(f"The fixed point is: {fixed_point.coordinates[0]}")

initial_unstable_segment = man_maker.get_initial_fundamental_segment(
    fixed_point, 0, 0, "unstable"
)
initial_stable_segment = man_maker.get_initial_fundamental_segment(
    fixed_point, 0, 0, "stable"
)

unstable_manifold = initial_unstable_segment
stable_manifold = initial_stable_segment

print(f"unstable stretchiness: {unstable_manifold.stretch_param}")
print(f"stable stretchiness: {stable_manifold.stretch_param}")

# Grow both manifolds. The old counts (12 unstable / 16 stable) predate the
# fsolve fixed-point fix; at k=10 the horseshoe folding makes the point count
# explode exponentially past ~8 iterations, so those counts now run for hours.
num_iterations = 8

man_machine.grow_x_times(fixed_point, "unstable", num_iterations)
unstable_manifold._find_tail()

man_machine.grow_x_times(fixed_point, "stable", num_iterations)
stable_manifold._find_tail()


unstable_points = unstable_manifold.get_point_array()
stable_points = stable_manifold.get_point_array()

plt.figure()

unstable_manifold.plot("blue")
stable_manifold.plot("red")
plt.scatter(*fixed_point.coordinates[0], c="k", s=7)

plt.xlim([-15, 15])
plt.ylim([-15, 15])
plt.show()

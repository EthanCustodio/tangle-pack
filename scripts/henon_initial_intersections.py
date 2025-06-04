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
# initial_guess = [6.104, 0]

fixed_point = fp_solver.construct_fixed_point(initial_guess, 2)

print(f"The fixed point is: {fixed_point.coordinates[0]}")

initial_unstable_segment = man_maker.get_initial_fundamental_segment(
    fixed_point, 0, 0, "unstable"
)
initial_stable_segment = man_maker.get_initial_fundamental_segment(
    fixed_point, 0, 0, "stable"
)

unstable_manifold = initial_unstable_segment
stable_manifold = initial_stable_segment

# grow unstable manifold
num_iterations = 12
for i in range(num_iterations):
    unstable_manifold = man_machine.grow_manifold(unstable_manifold)

# grow stable manifold
num_iterations = 16
for i in range(num_iterations):
    stable_manifold = man_machine.grow_manifold(stable_manifold)


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

plt.figure()

unstable_manifold.plot("blue")
stable_manifold.plot("red")
plt.scatter(*fixed_point.coordinates[0], c="k", s=7)

for point in tangle._intersecting_coords.values():
    plt.scatter(*point, c="k", s=7, zorder=10)

plt.xlim([-15, 15])
plt.ylim([-15, 15])
plt.show()

import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)

import tanglepack
from tanglepack.examples import (
    HENON_K10,
    henon_map as _henon_map_factory,
    henon_map_inverse as _henon_map_inverse_factory,
)
import numpy as np
import matplotlib.pyplot as plt


henon_map = _henon_map_factory(*HENON_K10)
henon_map_inverse = _henon_map_inverse_factory(*HENON_K10)


henon = tanglepack.DynamicalSystem(henon_map, henon_map_inverse)
fp_solver = tanglepack.FixedPointSolver(henon)
man_maker = tanglepack.ManifoldInitializer(henon)
man_machine = tanglepack.ManifoldMachine(henon)

tangle = tanglepack.Tangle()

initial_guess = [4, -4]

fixed_point = fp_solver.construct_fixed_point(initial_guess)

print(f"The fixed point is: {fixed_point.coordinates[0]}")

approx_dirs = {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}

man_maker.orient_manifolds(fixed_point, approx_dirs)

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
    """Index both manifolds and resolve every crossing between them.

    ``add_manifolds`` bulk-loads the whole segment set into the rtree in one pass
    and records the candidate crossing pairs; ``resolve_crossings`` then resolves
    each pair into an ``Intersection`` and hands them back -- the Tangle keeps only
    the segment index, never the crossings themselves.
    """
    tangle.add_manifolds([unstable_manifold, stable_manifold])
    return tangle.resolve_crossings()


crossings = intersections()
print(f"Intersections: {[ix.coords for ix in crossings]}")

fig = plt.figure()

show_points = False
unstable_manifold.plot("blue", show_points=show_points)
stable_manifold.plot("red", show_points=show_points)
plt.scatter(*fixed_point.coordinates[0], c="k", s=7)

for crossing in crossings:
    plt.scatter(*crossing.coords, c="k", s=7, zorder=10)

plt.xlim([-15, 15])
plt.ylim([-15, 15])

plt.show()

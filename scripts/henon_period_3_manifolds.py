import logging

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s %(levelname)-8s %(name)s: %(message)s",
    datefmt="%H:%M:%S",
)

import tanglepack
from tanglepack.examples import (
    HENON_P3,
    henon_map as _henon_map_factory,
    henon_map_inverse as _henon_map_inverse_factory,
    henon_jacobian as _henon_jacobian_factory,
)
import matplotlib.pyplot as plt


henon_map = _henon_map_factory(*HENON_P3)
henon_map_inverse = _henon_map_inverse_factory(*HENON_P3)
henon_jacobian = _henon_jacobian_factory(*HENON_P3)


henon = tanglepack.DynamicalSystem(henon_map, henon_map_inverse, henon_jacobian)
fp_solver = tanglepack.FixedPointSolver(henon)
man_maker = tanglepack.ManifoldInitializer(henon)
man_machine = tanglepack.ManifoldMachine(henon)

man_machine.area_cutoff = 1e-7

tangle = tanglepack.Tangle()

initial_guess = [[0, 1], [-1, 0], [-1, 1]]
initial_guess_zero = [4, -4]

fixed_point = fp_solver.construct_fixed_point(initial_guess)
fixed_point_zero = fp_solver.construct_fixed_point(initial_guess_zero)


print(f"fixed point: {fixed_point.coordinates}")

orbit_index = 0

approx_dir = {"unstable": [0, -1], "stable": [-1, -1]}

man_maker.orient_manifolds(fixed_point, approx_dir)

approx_dir = {"unstable": [-1, 0], "stable": [0, 1]}

man_maker.orient_manifolds(fixed_point_zero, approx_dir)

initial_unstable_segments = man_maker.construct_kevin_way(fixed_point, "unstable")
initial_stable_segments = man_maker.construct_kevin_way(fixed_point, "stable")

initial_zero_unstable = man_maker.construct_kevin_way(fixed_point_zero, "unstable")
initial_zero_stable = man_maker.construct_kevin_way(fixed_point_zero, "stable")

print(f"points! {initial_unstable_segments[(0, 0)].get_point_array()}")

# grow unstable manifold
num_iterations = 13
man_machine.grow_x_times(fixed_point, "unstable", num_iterations)
man_machine.grow_x_times(fixed_point, "stable", num_iterations)

man_machine.area_cutoff = 1e-4

# grow stable manifold
num_iterations = 7
man_machine.grow_x_times(fixed_point_zero, "stable", num_iterations)
man_machine.grow_x_times(fixed_point_zero, "unstable", num_iterations)

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


plt.xlim([-6, 6])
plt.ylim([-6, 6])

plt.title("k=2, b=1 Period 3 Nested Tangle")

plt.show()

import tanglepack
from tanglepack.examples import (
    HENON_K10,
    henon_map as _henon_map_factory,
    henon_map_inverse as _henon_map_inverse_factory,
)
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


henon_map = _henon_map_factory(*HENON_K10)
henon_map_inverse = _henon_map_inverse_factory(*HENON_K10)


# initialize the workbench
workbench = tanglepack.TangleWorkbench(henon_map, henon_map_inverse)

initial_guess = [4, -4]
# find the fixed point
fixed_point = workbench.construct_fixed_point(initial_guess)

print(f"The fixed point is: {fixed_point.coordinates[0]}")

approx_dirs = {"unstable": np.array([-1, 0]), "stable": np.array([0, 1])}
# orient the eigenvectors properly
workbench.orient_eigenvectors(fixed_point, approx_dirs)

# get the initial segments
unstable_segments, stable_segments = workbench.initialize_both_manifolds(fixed_point)

# grow the manifolds from the fixed point
# 8 growth iterations gives the 3+ bridges the rest of this script works with
# (the fsolve fixed-point fix made manifolds develop more slowly per iteration,
# so the old count of 6 no longer produced any crossings)
workbench.grow_n_times(fixed_point, "unstable", num_iterations=8)
workbench.grow_until_turnaround(fixed_point, "stable")

intersections = workbench.compute_intersections(fixed_point)
print(intersections)
print(len(intersections))

# for each bridge find the intersections that are the endpoints
# add a directed edge between the two intersections

# define bridge class with the endpoints as part of the class

# intersection point within structure should have which stable segment it is on
# should store which stable branch and which unstable branch it is on

# start thinking about how to write up description of algorithms that will be in the
# paper. Think about what style to write that in

# read Kevin Pip algorithm

bridges = workbench.create_bridges(fixed_point)
print(f"num bridges {len(bridges)}")
print(f"type bridges: {type(bridges)}")
print(bridges)

workbench.trim_stable_manifolds(fixed_point)

plt.figure()
workbench.plot_tangle(fixed_point, "stable", color="r")
workbench.plot_intersections(fixed_point)
workbench.plot_all_bridges()

plt.xlim([-15, 15])
plt.ylim([-15, 15])

new_bridges = workbench.iterate_bridge(bridges[2])

plt.figure()
workbench.plot_tangle(fixed_point, "stable", color="r")
workbench.plot_intersections(fixed_point)
workbench.plot_all_bridges()

plt.xlim([-15, 15])
plt.ylim([-15, 15])

workbench.iterate_bridge(new_bridges[0])
# plot!
plt.figure()
workbench.plot_tangle(fixed_point, "stable", color="r")
workbench.plot_intersections(fixed_point)
workbench.plot_all_bridges()

plt.xlim([-15, 15])
plt.ylim([-15, 15])

plt.show()

import tanglepack
import numpy as np
import networkx as nx
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
workbench.grow_n_times(fixed_point, "unstable", num_iterations=6)
# workbench.grow_until_arclength(fixed_point, "unstable", 60)
# workbench.grow_n_times(fixed_point, "stable", num_iterations=4)
workbench.grow_until_turnaround(fixed_point, "stable")
# workbench.grown_until_intersection(fixed_point, "unstable")

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

# print(f"These are the new bridges: {bridges}")
# new_bridges = workbench._man_machine.iterate_bridge(bridges[2])

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
# workbench.iterate_all_bridges()
# graph = workbench.build_intersection_graph(fixed_point)

# print(len(graph.nodes()))

# workbench.visualize_intersection_graph(graph)

# plot!
plt.figure()
# workbench.plot_tangle(fixed_point, "unstable", color="b")
workbench.plot_tangle(fixed_point, "stable", color="r")
workbench.plot_intersections(fixed_point)
workbench.plot_all_bridges()
# new_bridges.plot()
# bridges.plot()

plt.xlim([-15, 15])
plt.ylim([-15, 15])

plt.show()

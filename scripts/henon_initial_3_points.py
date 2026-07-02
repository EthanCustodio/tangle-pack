import tanglepack
import numpy as np
import matplotlib.pyplot as plt


def henon_map(point):
    """defines the henon map for binary horshoe parameters to test basic functionality"""

    k, b = (10, 1)

    x = point[0]
    y = point[1]

    return np.array([y - k + x ** 2, -b * x])


def henon_map_inverse(point):
    """defines the inverse henon map for"""

    k, b = (10, 1)

    x = point[0]
    y = point[1]

    return np.array([-y / b, x + k - (y ** 2) / (b ** 2)])


henon = tanglepack.DynamicalSystem(henon_map, henon_map_inverse)

initial_guess = [4, -4]

fp_solver = tanglepack.FixedPointSolver(henon)

fixed_point = fp_solver.construct_fixed_point(initial_guess, 1)

print(f'The fixed point is: {fixed_point.coordinates[0]}')
print(f'The fixed point is type: {type(fixed_point)}')

man_maker = tanglepack.ManifoldInitializer(henon)

initial_segment = man_maker.get_initial_fundamental_segment(fixed_point, 0, 0, 'unstable')
initial_points = initial_segment.get_point_array()
print(f'Initial 3 Points: {initial_points}')


initial_segment.plot(marker='o', ms=10)


# plt.figure()
# plt.scatter(*fixed_point.coordinates[0], c='k', s=1)
# plt.scatter(initial_points[:, 0], initial_points[:, 1], c='g', s=1)
# plt.show()


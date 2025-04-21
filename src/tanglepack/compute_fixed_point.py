import numpy as np
from scipy.optimize import newton as newton_method
from functools import partial
from scipy.differentiate import jacobian


def compute_fixed_point(initial_guess, transform):
    """Computes a fixed point based on the initial guess"""

    multipoint_difference_with_transform = partial(multipoint_shoot_flattened_difference, transform=transform)

    initial_guess_flattened = flatten_trajectory(initial_guess)

    fixed_point_flattened = newton_method(multipoint_difference_with_transform, initial_guess_flattened)

    fixed_point_full = unflatten_trajectory(fixed_point_flattened)

    return fixed_point_full


def multipoint_shoot_flattened_difference(trajectory, transform):
    """Takes the difference between an iterate and the current trajectory"""

    shoot = multipoint_shoot_flattened(trajectory, transform)

    difference = shoot - trajectory

    return difference


def multipoint_shoot_flattened(trajectory, transform):
    """takes a flattened trajectory and iterates it forward"""

    trajectory_full = unflatten_trajectory(trajectory)

    trajectory_full_transformed = multipoint_shoot(trajectory_full, transform)

    trajectory_flattened_transformed = flatten_trajectory(trajectory_full_transformed)

    return trajectory_flattened_transformed


def flatten_trajectory(trajectory):
    """takes a trajectory and makes it a 2n x 1 vector where n is the number of iterates"""

    number_iterates, _ = np.shape(trajectory)

    trajectory_reshaped = np.reshape(trajectory, (2 * number_iterates, 1))

    return trajectory_reshaped


def unflatten_trajectory(trajectory):
    """takes a flattened trajectory and reformats it back into a n x d matrix"""

    number_iterates_doubled, _ = np.shape(trajectory)

    number_iterates = number_iterates_doubled // 2

    trajectory_reshaped = np.reshape(trajectory, (number_iterates, 2))

    return trajectory_reshaped
    

def multipoint_shoot(trajectory, transform):
    """Takes a trajectory and transforms it forward once"""

    period, _ = np.shape(trajectory)

    transformed_trajectory = np.array([transform(trajectory[i, :]) for i in range(period)])

    return transformed_trajectory


def henon(point):
    """defines the henon map for classical parameters to test basic functionality"""

    a, b = (1.4, 0.3)

    x = point[0]
    y = point[1]

    return [1 - a * x ** 2 + y, b * x]

guess = np.array([[0.6, 0.2], [0.6, 0.2]])
# print(guess)
# print(multipoint_shoot(guess, henon))
# temp = flatten_trajectory(guess)
# print(temp)
# print(unflatten_trajectory(temp))
# print(multipoint_shoot_flattened(temp, henon))
print(compute_fixed_point(guess, henon))

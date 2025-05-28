import numpy as np
from .BranchPoint import BranchPoint


class FixedPoint:

    def __init__(self, period, num_branches):

        self.period = period
        self.num_branches = num_branches

        self.branch_points = [BranchPoint(num_branches) for _ in range(period)]
        self.coordinates = [np.empty((2, 1)) for _ in range(period)]

        self.unstable_eigenvectors = [np.empty((2, 1)) for _ in range(period)]
        self.unstable_eigenvalues = [0.0] * period

        self.stable_eigenvectors = [np.empty((2, 1)) for _ in range(period)]
        self.stable_eigenvalues = [0.0] * period

        self.accuracy = 0.0

        self.jacobians = [np.empty((2, 2)) for _ in range(period)]
        self.total_jacobian = np.empty((2, 2))

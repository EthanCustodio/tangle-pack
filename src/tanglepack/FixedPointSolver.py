import numpy as np
from scipy.optimize import newton as newton_method
from scipy.differentiate import jacobian as jacob


class FixedPointSolver():


    def __init__(self, dynamical_map, dynamical_map_inverse=None, jacobian_function=None):

        
        self.dynamical_map = dynamical_map
        self.dynamical_map_inverse = dynamical_map_inverse
        self.jacobian_function = jacobian_function

        
class DynamicalSystem:

    def __init__(
        self,
        dynamical_map,
        dynamical_map_inverse,
        jacobian_function=None,
        name="unnamed",
    ):

        self.map = dynamical_map
        self.map_inv = dynamical_map_inverse
        self.jacobian = jacobian_function
        self.name = name

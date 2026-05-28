from typing import Callable
import numpy.typing as npt
from typing_extensions import Annotated

from numpy import float64

# A 2D point: exactly shape (2,) and dtype float64
Point2D = Annotated[npt.NDArray[float64], (2,)]

# Function type: takes a 2D point, returns a 2D point
MapFunc = Callable[[Point2D], Point2D]

Matrix2D = Annotated[npt.NDArray[float64], (2, 2)]

# Function type for jacobian function: takes a 2D point, returns a 2D matrix
JacFunc = Callable[[Point2D], Matrix2D]


class DynamicalSystem:
    """
    Object which contains the map functions for a dynamical system. This object also
    stores and enforces the types that the mapping functions should handle.

    Attributes:
        map (MapFunc): The dynamical map of the system. Takes a point on the plane and
            returns a point on the plane.
        map_inv (MapFunc): The inverse dynamical map of the system. Takes a point on
            the plane and returns a point on the plane. map(map_inv()) is the identity.
        jacobian (JacFunc): Optional function which takes a point on the plane and
            returns the jacobian at that point as a 2x2 matrix. Useful to include if
            you have a fast method to compute the jacobian.
        name (str): Optional name of the dynamical system.
    """

    def __init__(
        self,
        dynamical_map: MapFunc,
        dynamical_map_inverse: MapFunc,
        jacobian_function: JacFunc | None = None,
        name: str = "unnamed",
    ):
        """
        Initalizes the system with the mapping functions.

        Args:
            dynamical_map (MapFunc): The dynamical map of the system.
            dynamical_map_inverse (MapFunc): The inverse dynamical map of the system.
            jacobian_function (JacFunc | None, optional): Optional function which takes
                a point on the plane and returns the jacobian at that point as a
                2x2 matrix. Defaults to None.
            name (str, optional): Optional name of the dynamical system.
                Defaults to "unnamed".
        """

        self.map = dynamical_map
        self.map_inv = dynamical_map_inverse
        self.jacobian = jacobian_function
        self.name = name

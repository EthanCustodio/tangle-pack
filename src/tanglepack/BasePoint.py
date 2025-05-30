from typing import Literal
import numpy as np


class BasePoint:

    def __init__(self, x=None, y=None, cdist=None, edist=None):

        self.x = x
        self.y = y

        self.cdist = cdist
        self.edist = edist

        self._set_coords()

    def get_point(self):

        return self._coords

    def get_cdist(self, stability: Literal["unstable", "stable"]) -> float:

        return self.cdist

    def set_x(self, x: float):

        self.x = x
        self._set_coords()

    def set_y(self, y: float):

        self.y = y
        self._set_coords()

    def _set_coords(self):

        self._coords = np.array([self.x, self.y])

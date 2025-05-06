import numpy as np

class BasePoint():

    def __init__(self, x=None, y=None, cdist=None, edist=None):

        self.x = x
        self.y = y

        self.cdist = cdist
        self.edist = edist


    def get_point(self):
        
        return np.array([self.x, self.y])

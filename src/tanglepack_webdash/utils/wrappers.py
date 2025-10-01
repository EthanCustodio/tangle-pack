# src/tanglepack_webdash/utils/wrappers.py
import numpy as np


def pointize(F):
    """
    Wrap F so it always sees a single 2-vector and returns a single 2-vector.
    Any batched/odd shapes get flattened to the first two numbers.
    """

    def G(p):
        a = np.asarray(p, float).reshape(-1)  # flatten
        a = a[:2]  # take first two
        return np.asarray(F(a), float).reshape(2)

    return G

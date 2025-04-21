import pytest
import numpy as np
from tanglepack.BranchPoint import BranchPoint

def test_point_creation():
    p = BranchPoint(2, 1.0, 2.0)

    assert p.x == 1.0
    assert p.y == 2.0
    assert p.forward == None
    assert p.backward == None


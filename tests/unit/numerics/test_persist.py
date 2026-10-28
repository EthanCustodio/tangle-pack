"""Saving and loading a workbench (N10a): the format, the file contents and the slot writer.

Load writes every slot through :meth:`FixedPoint.attach_chain`, the one writer
of the slot table, and its input is now a file: a misspelt stability raises
instead of filling a slot no reader looks up.
"""

from __future__ import annotations

import pytest

from helpers.fakes import bare_fixed_point
from tanglepack.numerics.PointStore import PointStore


def _chain():
    store = PointStore()
    return store.new_chain(store.append([(1.0, 0.0), (2.0, 0.0)], [1.0, 2.0]), 2.0)


def test_attach_chain_refuses_a_misspelt_stability():
    with pytest.raises(ValueError, match="'Unstable'"):
        bare_fixed_point().attach_chain(_chain(), "Unstable", 0, 0)


def test_attach_chain_refuses_a_third_branch():
    with pytest.raises(IndexError):
        bare_fixed_point().attach_chain(_chain(), "unstable", 0, 2)

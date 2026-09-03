"""Pins the Phase 8 split of ``TangleWorkbench`` into an orchestrator plus three
collaborator modules.

The split is a pure move: the workbench keeps every public method it had, with
the same signature, and simply forwards to the module that now owns the body.
These tests pin BOTH halves of that -- that the new modules exist and hold the
implementations, and that nothing a caller can see moved.
"""

from __future__ import annotations

import importlib
import inspect

import pytest

import tanglepack
from tanglepack.numerics import graphviz as graphviz_module
from tanglepack.numerics.TangleWorkbench import TangleWorkbench

# ``numerics/__init__`` re-exports the CLASS under each module's name, the way it
# does for every other one-class module here, so the modules themselves are
# reached through the import system rather than as package attributes.
bridge_iterator_module = importlib.import_module("tanglepack.numerics.BridgeIterator")
iterate_inference_module = importlib.import_module(
    "tanglepack.numerics.IterateInference"
)


def test_top_level_import_still_works() -> None:
    """The public entry point is unchanged by the split."""
    assert tanglepack.TangleWorkbench is TangleWorkbench
    assert tanglepack.numerics.TangleWorkbench is TangleWorkbench


def test_new_modules_import_cleanly() -> None:
    """The three collaborator modules exist and export their owners."""
    assert bridge_iterator_module.BridgeIterator is tanglepack.numerics.BridgeIterator
    assert (
        iterate_inference_module.IterateInference
        is tanglepack.numerics.IterateInference
    )
    assert callable(graphviz_module.build_intersection_graph)
    assert callable(graphviz_module.visualize_intersection_graph)


@pytest.mark.parametrize(
    "name, owner",
    [
        ("iterate_bridge", "BridgeIterator"),
        ("iterate_all_bridges", "BridgeIterator"),
        ("image_bridges", "BridgeIterator"),
        ("preimage_bridges", "BridgeIterator"),
        ("infer_iterates", "IterateInference"),
        ("infer_iterate_table", "IterateInference"),
        ("build_intersection_graph", "graphviz"),
        ("visualize_intersection_graph", "graphviz"),
    ],
)
def test_moved_bodies_live_in_their_new_module(name: str, owner: str) -> None:
    """Each moved method's body is defined outside ``TangleWorkbench.py``."""
    method = getattr(TangleWorkbench, name)
    source = inspect.getsource(method)
    assert "TangleWorkbench.py" in inspect.getsourcefile(method)
    # The workbench keeps only a thin forwarder: the real body is elsewhere.
    assert len(source.splitlines()) < 45, f"{name} still holds its body"
    module = {
        "BridgeIterator": bridge_iterator_module,
        "IterateInference": iterate_inference_module,
        "graphviz": graphviz_module,
    }[owner]
    assert any(
        name in vars(obj)
        for obj in (module, getattr(module, owner, None))
        if obj is not None
    ), f"{name} not found in {owner}"


@pytest.mark.parametrize(
    "name, params",
    [
        ("iterate_bridge", ["self", "bridge"]),
        ("iterate_all_bridges", ["self"]),
        ("image_bridges", ["self", "bridge_id", "n"]),
        ("preimage_bridges", ["self", "bridge_id", "n"]),
        ("infer_iterates", ["self", "cdist_rtol"]),
        ("infer_iterate_table", ["self", "cdist_rtol"]),
        ("build_intersection_graph", ["self"]),
    ],
)
def test_delegating_signatures_unchanged(name: str, params: list[str]) -> None:
    """A caller sees exactly the signature it saw before the split."""
    signature = inspect.signature(getattr(TangleWorkbench, name))
    assert list(signature.parameters) == params


def test_visualize_signature_unchanged() -> None:
    """The plotting entry point keeps all of its keyword options."""
    signature = inspect.signature(TangleWorkbench.visualize_intersection_graph)
    assert list(signature.parameters) == [
        "self",
        "G",
        "layout",
        "figsize",
        "display_mode",
        "compact_threshold",
        "node_size",
        "label_mode",
        "node_color_by",
        "show_iterate_edges",
        "save_path",
    ]


def test_collaborators_are_wired_to_their_workbench(workbench) -> None:
    """The workbench owns one collaborator of each kind, pointing back at it."""
    assert isinstance(workbench._bridge_iterator, bridge_iterator_module.BridgeIterator)
    assert isinstance(
        workbench._iterate_inference, iterate_inference_module.IterateInference
    )
    assert workbench._bridge_iterator.workbench is workbench
    assert workbench._iterate_inference.workbench is workbench


def test_iterate_bridge_still_patchable_on_the_workbench(
    henon_tangle_with_bridges, monkeypatch
) -> None:
    """``iterate_all_bridges`` goes through the workbench attribute, so the
    blast-recovery tests' ``monkeypatch.setattr(workbench, "iterate_bridge", ...)``
    keeps intercepting every iteration."""
    workbench, _fp = henon_tangle_with_bridges
    calls: list[object] = []

    def spy(bridge):
        calls.append(bridge)
        return []

    monkeypatch.setattr(workbench, "iterate_bridge", spy)
    assert workbench.iterate_all_bridges() == []
    assert len(calls) == len(workbench.bridges)

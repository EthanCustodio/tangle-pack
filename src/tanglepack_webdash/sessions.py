# src/tanglepack_webdash/session.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, Optional
import threading
import time

from tanglepack.TangleWorkbench import TangleWorkbench


@dataclass
class WBState:
    workbench: Optional[TangleWorkbench] = None
    fp: Optional[object] = None
    fig: Optional[object] = None

    # Bridges is keyed like (seg1.id, seg2.id)
    bridges: Optional[dict] = None
    created_at: float = field(default_factory=time.time)
    # extend here: e.g., last_manifolds: dict, ui_prefs: dict, etc.


# In-memory registry (single-process). Swap for Redis later if needed.
_REGISTRY: Dict[str, WBState] = {}
_LOCK = threading.Lock()


def get_state(session_id: str) -> WBState:
    """Get (or create) the state for a session id."""
    with _LOCK:
        st = _REGISTRY.get(session_id)
        if st is None:
            st = WBState()
            _REGISTRY[session_id] = st
        return st


def drop_state(session_id: str) -> None:
    with _LOCK:
        _REGISTRY.pop(session_id, None)


def sweep(max_items: int = 1000) -> None:
    """Optional: keep memory bounded in long dev sessions."""
    with _LOCK:
        if len(_REGISTRY) > max_items:
            # naive LRU-ish sweep by created_at
            oldest = sorted(_REGISTRY.items(), key=lambda kv: kv[1].created_at)[
                : len(_REGISTRY) - max_items
            ]
            for sid, _ in oldest:
                _REGISTRY.pop(sid, None)

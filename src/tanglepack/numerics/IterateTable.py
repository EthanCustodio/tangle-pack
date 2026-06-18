from __future__ import annotations

from typing import Optional
import numpy as np
from numpy.typing import NDArray


class IterateTable:
    """
    2D lookup: (intersection_id, n) → intersection_id.

    Supports:
        table[3, 2]  → ID of f^2(intersection 3)
        table[3, -1] → ID of f^{-1}(intersection 3)
        table[3, 0]  → 3  (identity)

    Setting one direction auto-records the reverse:
        table[3, 2] = 7  also records  table[7, -2] = 3

    Attributes:
        _forward: dict[int, dict[int, int]]  — _forward[id][n] = target_id  (n > 0)
        _backward: dict[int, dict[int, int]] — _backward[id][n] = target_id (n > 0)
    """

    def __init__(self):
        self._forward: dict[int, dict[int, int]] = {}
        self._backward: dict[int, dict[int, int]] = {}

    def __getitem__(self, key: tuple[int, int]) -> Optional[int]:
        source_id, n = key
        if n == 0:
            return source_id
        if n > 0:
            return self._forward.get(source_id, {}).get(n)
        else:
            return self._backward.get(source_id, {}).get(-n)

    def __setitem__(self, key: tuple[int, int], target_id: int):
        source_id, n = key
        if n == 0:
            return
        if n > 0:
            self._forward.setdefault(source_id, {})[n] = target_id
            self._backward.setdefault(target_id, {})[n] = source_id
        else:
            self._backward.setdefault(source_id, {})[-n] = target_id
            self._forward.setdefault(target_id, {})[-n] = source_id

    def __contains__(self, key: tuple[int, int]) -> bool:
        return self[key] is not None

    def forward_depth(self, source_id: int) -> int:
        """Maximum forward iterate depth recorded for this intersection."""
        return max(self._forward.get(source_id, {}).keys(), default=0)

    def backward_depth(self, source_id: int) -> int:
        """Maximum backward iterate depth recorded for this intersection."""
        return max(self._backward.get(source_id, {}).keys(), default=0)

    def forward_chain(self, source_id: int) -> list[int]:
        """Return [id, f(id), f^2(id), ...] for all recorded forward iterates."""
        chain = [source_id]
        n = 1
        while (nxt := self[source_id, n]) is not None:
            chain.append(nxt)
            n += 1
        return chain

    def backward_chain(self, source_id: int) -> list[int]:
        """Return [id, f^{-1}(id), f^{-2}(id), ...] for all recorded backward iterates."""
        chain = [source_id]
        n = 1
        while (prev := self[source_id, -n]) is not None:
            chain.append(prev)
            n += 1
        return chain

    def all_registered_ids(self) -> set[int]:
        """All intersection IDs that appear in any entry of this table."""
        ids = set(self._forward.keys()) | set(self._backward.keys())
        for sub in self._forward.values():
            ids.update(sub.values())
        for sub in self._backward.values():
            ids.update(sub.values())
        return ids

    def as_forward_array(self, ids: list[int], max_depth: int) -> NDArray[np.int64]:
        """
        Dense 2D array A where A[i, d-1] = ID of f^d(ids[i]), or -1 if unknown.

        Args:
            ids: Ordered list of intersection IDs (defines row order).
            max_depth: Number of forward iterate columns.

        Returns:
            Array of shape (len(ids), max_depth), dtype int64.
        """
        arr = np.full((len(ids), max_depth), fill_value=-1, dtype=np.int64)
        for i, iid in enumerate(ids):
            for d in range(1, max_depth + 1):
                target = self[iid, d]
                if target is not None:
                    arr[i, d - 1] = target
        return arr

    def as_backward_array(self, ids: list[int], max_depth: int) -> NDArray[np.int64]:
        """Dense array B[i, d-1] = ID of f^{-d}(ids[i]), or -1 if unknown."""
        arr = np.full((len(ids), max_depth), fill_value=-1, dtype=np.int64)
        for i, iid in enumerate(ids):
            for d in range(1, max_depth + 1):
                target = self[iid, -d]
                if target is not None:
                    arr[i, d - 1] = target
        return arr

    def register_iterate(self, source_id: int, n: int, target_id: int):
        """Explicit named method — delegates to __setitem__."""
        self[source_id, n] = target_id

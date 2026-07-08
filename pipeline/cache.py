"""On-disk + in-memory caching for KEGG REST lookups.

The notebook used a bare ``@lru_cache`` on compound lookups.  Here we keep the
in-process LRU but back it with a small JSON file on disk so repeat runs and
re-fetches don't re-hit KEGG.  Both compound *and* reaction entries live in the
same store, keyed by their KEGG id (``Cxxxxx`` / ``Rxxxxx`` / ``Mxxxxx``).
"""

from __future__ import annotations

import json
import os
import threading
from typing import Any, Optional

# Cache file lives next to the package by default; override with NETCONV_CACHE.
_DEFAULT_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), ".cache")
_CACHE_PATH = os.environ.get(
    "NETCONV_CACHE", os.path.join(_DEFAULT_DIR, "kegg_cache.json")
)

_lock = threading.Lock()
_store: Optional[dict] = None


def _load() -> dict:
    global _store
    if _store is None:
        try:
            with open(_CACHE_PATH, "r", encoding="utf-8") as fh:
                _store = json.load(fh)
        except (FileNotFoundError, json.JSONDecodeError):
            _store = {}
    return _store


def get(key: str) -> Optional[Any]:
    """Return a cached value for ``key`` or ``None`` if absent."""
    with _lock:
        return _load().get(key)


def put(key: str, value: Any) -> None:
    """Store ``value`` under ``key`` and flush to disk."""
    with _lock:
        store = _load()
        store[key] = value
        os.makedirs(os.path.dirname(_CACHE_PATH), exist_ok=True)
        tmp = _CACHE_PATH + ".tmp"
        with open(tmp, "w", encoding="utf-8") as fh:
            json.dump(store, fh)
        os.replace(tmp, _CACHE_PATH)

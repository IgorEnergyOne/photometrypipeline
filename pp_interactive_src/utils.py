# -*- coding: utf-8 -*-
"""
Shared utility / helper functions for pp_interactive_src.
"""
from __future__ import annotations

import os
import re
from typing import Callable

from .constants import DEBUG


def _combine(existing, new):
    """Accumulate matplotlib artist handles into a list.

    Used to merge multiple ghost series (e.g. the left-fold copy and the
    right-echo copy) into a single entry slot so the fast-path visibility
    toggle can handle them with one call to ``_set_visible``.
    """
    if new is None:
        return existing
    if existing is None:
        return new
    # Flatten into a list so _set_visible can iterate
    result = existing if isinstance(existing, list) else [existing]
    result = result + (new if isinstance(new, list) else [new])
    return result


def _safe(fn: Callable, *args, default=None, **kwargs):
    """Call *fn* with the given arguments, silently returning *default* on any
    exception.  Use this for fire-and-forget GUI operations (geometry tweaks,
    widget config updates, etc.) where a failure should never crash the app."""
    try:
        return fn(*args, **kwargs)
    except Exception:
        return default


def debug_print(*args, **kwargs):
    """Print debug messages to stdout only when the module-level DEBUG flag is True."""
    if DEBUG:
        print("[DEBUG]", *args, **kwargs)


def next_version(path: str) -> str:
    """Return a non-existing filename by appending or incrementing a ``_N`` suffix.

    Examples::

        next_version("out.csv")       # -> "out_1.csv"  (if out_1.csv is absent)
        next_version("out_1.csv")     # -> "out_2.csv"  (if out_2.csv is absent)
    """
    base, ext = os.path.splitext(path)
    # Check for existing regex pattern _(\d+)$ at the end of base
    match = re.search(r'_(\d+)$', base)
    if match:
        counter = int(match.group(1))
        prefix = base[:match.start()]
    else:
        counter = 1
        prefix = base

    while True:
        candidate = f"{prefix}_{counter}{ext}"
        if not os.path.exists(candidate):
            return candidate
        counter += 1


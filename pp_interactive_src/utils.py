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


_PROVISIONAL_DESIG = re.compile(r'^\d{4}\s+(?:[A-Z]{2}\d*|P-L|T-[123])$')


def format_asteroid_name(name: str) -> str:
    """Return *name* in the conventional form with the number in brackets.

    JPL target names put the provisional designation in brackets instead::

        format_asteroid_name("153201 (2000 WO107)")   # -> "(153201) 2000 WO107"
        format_asteroid_name("1024 Hale (A923 XC)")   # -> "(1024) Hale"
        format_asteroid_name("(2009 FF)")             # -> "2009 FF"
        format_asteroid_name("1627 Ivar")             # -> "(1627) Ivar"

    Anything not recognised as an asteroid name is returned unchanged.
    """
    s = str(name).strip()
    # numbered: "<number> [<name>] (<designation>)"
    m = re.match(r'^(\d+)\s+(.*?)\s*\(([^()]+)\)$', s)
    if m:
        number, proper, desig = m.groups()
        return f"({number}) {proper or desig}"
    # unnumbered: "(<designation>)"
    m = re.match(r'^\(([^()]*[^\d()][^()]*)\)$', s)
    if m:
        return m.group(1).strip()
    if s.isdigit():
        return f"({s})"
    # numbered without designation: "<number> <name>" (but not "2025 FA22")
    m = re.match(r'^(\d+)\s+(\S.*)$', s)
    if m and not _PROVISIONAL_DESIG.match(s):
        return f"({m.group(1)}) {m.group(2)}"
    return s


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


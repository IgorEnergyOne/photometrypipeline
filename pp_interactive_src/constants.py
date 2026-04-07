# -*- coding: utf-8 -*-
"""
Global constants and configuration for pp_interactive_src.
"""

# Debug output: set True to enable verbose console logging
DEBUG: bool = False

# Default rotation-period step (hours) used when incrementing/decrementing
# the period with the z/x hotkeys.
TIME_STEP: float = 0.02

# Default period (hours) shown in the rotation-phase controls on startup.
DEFAULT_PERIOD: float = 4.0

# Default maximum rotation phase shown on the X-axis (values > 1.0 create a
# wrapped "ghost" copy of the first part of the phase curve, drawn in a
# lighter shade so the user can distinguish repeated coverage).
DEFAULT_PHASE_MAX: float = 1.2

# Initial main-window dimensions (pixels).
WINDOW_WIDTH: int = 1400
WINDOW_HEIGHT: int = 720

# Sentinel key used in the per-alias offset dict when no per-band offset is
# set  -  the value stored under this key acts as a fallback for all bands.
OFFSET_ALL_KEY: str = "__ALL__"


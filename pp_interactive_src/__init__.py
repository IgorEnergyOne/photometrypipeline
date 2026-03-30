# -*- coding: utf-8 -*-
"""
pp_interactive_src package – interactive light-curve viewer & editor.
Usage
-----
As a script (via the thin launcher)::
    ./pp_interactive_src.py [--debug]
As a module::
    python -m pp_interactive_src [--debug]
Public re-exports
-----------------
The most commonly needed symbols are re-exported here so that other tools in
the photometrypipeline suite can from pp_interactive_src import LightCurveGUI
without knowing the internal module layout.
"""
from .constants import (  # noqa: F401
    DEBUG, TIME_STEP, DEFAULT_PERIOD, DEFAULT_PHASE_MAX,
    WINDOW_WIDTH, WINDOW_HEIGHT, OFFSET_ALL_KEY,
)
from .models import FitsContext, LightCurveData          # noqa: F401
from .plot_settings import PlotSettings                  # noqa: F401
from .plot import BlitManager, LightCurvePlot            # noqa: F401
from .image_viewer import AsteroidImageViewer            # noqa: F401
from .gui import LightCurveGUI                           # noqa: F401
from .main import main                                   # noqa: F401
__all__ = [
    "DEBUG", "TIME_STEP", "DEFAULT_PERIOD", "DEFAULT_PHASE_MAX",
    "WINDOW_WIDTH", "WINDOW_HEIGHT", "OFFSET_ALL_KEY",
    "FitsContext", "LightCurveData",
    "PlotSettings",
    "BlitManager", "LightCurvePlot",
    "AsteroidImageViewer",
    "LightCurveGUI",
    "main",
]

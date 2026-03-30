#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
pp_interactive_src – thin launcher.
All logic lives in the pp_interactive_src/ package.
Run directly::
    ./pp_interactive_src.py [--debug]
or as a module::
    python -m pp_interactive_src [--debug]
This file is kept in the project root so that the original invocation
pattern continues to work unchanged.
"""
from pp_interactive_src.main import main
if __name__ == "__main__":
    main()
# ── Legacy compatibility shim ──────────────────────────────────────────────
# Re-export the public API so any code that previously imported from here still works.
from pp_interactive_src import (  # noqa: F401
    DEBUG, TIME_STEP, DEFAULT_PERIOD, DEFAULT_PHASE_MAX,
    WINDOW_WIDTH, WINDOW_HEIGHT, OFFSET_ALL_KEY,
    FitsContext, LightCurveData,
    PlotSettings,
    BlitManager, LightCurvePlot,
    AsteroidImageViewer,
    LightCurveGUI,
)

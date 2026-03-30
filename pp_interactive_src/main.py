# -*- coding: utf-8 -*-
"""
Entry point for pp_interactive_src.
"""
from __future__ import annotations
import argparse
import ttkbootstrap as ttk
from .constants import DEBUG, WINDOW_WIDTH, WINDOW_HEIGHT
from .gui import LightCurveGUI
from . import constants as _consts
from . import utils as _utils
def main():
    parser = argparse.ArgumentParser(description="Interactive lightcurve viewer")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    args = parser.parse_args()
    # Propagate debug flag to the constants module so debug_print() works globally
    _consts.DEBUG = args.debug
    # Also patch the utils module which caches DEBUG at import time
    _utils.DEBUG = args.debug
    root = ttk.Window(themename="flatly")
    root.geometry(f"{WINDOW_WIDTH}x{WINDOW_HEIGHT}")
    app = LightCurveGUI(root)  # noqa: F841
    root.mainloop()
if __name__ == "__main__":
    main()

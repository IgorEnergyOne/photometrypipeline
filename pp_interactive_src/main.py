# -*- coding: utf-8 -*-
"""
Entry point for pp_interactive_src.
"""
from __future__ import annotations
import argparse
import ttkbootstrap as ttk
from .config import load_config, save_user_config
from .gui import LightCurveGUI
from . import constants as _consts
from . import utils as _utils


def main():
    parser = argparse.ArgumentParser(description="Interactive lightcurve viewer")
    parser.add_argument("--debug", action="store_true", help="Enable debug output")
    args = parser.parse_args()

    # Load merged config (default + user overrides)
    cfg = load_config()

    # CLI --debug overrides the config file value
    if args.debug:
        cfg.app.debug = True

    # Propagate debug flag to legacy modules that cache it as a module attribute
    _consts.DEBUG = cfg.app.debug
    _utils.DEBUG = cfg.app.debug

    # Also back-fill the legacy constant names so older import sites still work
    _consts.WINDOW_WIDTH   = cfg.window.width
    _consts.WINDOW_HEIGHT  = cfg.window.height
    _consts.DEFAULT_PERIOD = cfg.defaults.period
    _consts.TIME_STEP      = cfg.defaults.period_step
    _consts.DEFAULT_PHASE_MAX = cfg.defaults.phase_max

    root = ttk.Window(themename=cfg.app.theme)
    root.geometry(f"{cfg.window.width}x{cfg.window.height}")
    app = LightCurveGUI(root, cfg=cfg)  # noqa: F841
    root.mainloop()


if __name__ == "__main__":
    main()

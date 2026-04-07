# -*- coding: utf-8 -*-
"""
ConfigDialog - tabbed GUI for editing the application configuration.

Opened from  Plot Settings -> App Settings...

Four tabs
---------
* **App & Window** - theme, window size  (restart required)
* **Defaults**     - startup mode, time axis, legend, period, etc.
* **Plot Style**   - markers, error bars, ghost, grid, fonts, wrap line
* **Colors**       - flagged/rejected colours, default palette, band colours
"""
from __future__ import annotations

import copy
from dataclasses import asdict
from tkinter import colorchooser
from typing import Callable, Dict, List, Optional, Tuple

import tkinter as tk
import ttkbootstrap as ttk
from ttkbootstrap.constants import *

from .config import AppConfig, ColorsSection, save_user_config

# -- ttkbootstrap themes (light then dark) ------------------------------------
_THEMES = [
    "flatly", "cosmo", "journal", "litera", "lumen", "materia",
    "minty", "pulse", "sandstone", "simplex", "sketchy", "spacelab",
    "united", "yeti",
    "darkly", "cyborg", "slate", "solar", "superhero", "vapor",
]

_FONT_SIZES = [
    "xx-small", "x-small", "small", "medium",
    "large", "x-large", "xx-large",
]

_MODES        = ["target", "instrumental", "control", "relative"]
_TIME_MODES   = ["minutes", "julian_date", "mjd", "rotation_phase"]
_ERR_TYPES    = ["calibrated", "instrumental", "relative", "none"]
_LEGEND_LOCS  = ["upper left", "upper right", "lower left", "lower right", "center"]
_MARKERS      = [
    ("o", "circle"), ("s", "square"), ("p", "pentagon"), ("x", "x"),
    ("D", "diamond"), ("*", "star"), ("v", "triangle_down"), ("^", "triangle_up"),
    ("<", "triangle_left"), (">", "triangle_right"), ("+", "plus"),
    ("d", "thin_diamond"),
]
_MARKER_NAMES = [name for _, name in _MARKERS]
_MARKER_MAP   = {name: sym for sym, name in _MARKERS}
_MARKER_NAMES_MAP = {sym: name for sym, name in _MARKERS}


# ---------------------------------------------------------------------------
# Swatch helper widget
# ---------------------------------------------------------------------------

def _make_swatch(parent, color: str, size: Tuple[int, int] = (28, 18)) -> tk.Canvas:
    """Return a small Canvas widget painted with *color*."""
    w, h = size
    cv = tk.Canvas(parent, width=w, height=h, highlightthickness=1, bd=0,
                   cursor="hand2")
    cv.create_rectangle(0, 0, w, h, fill=color, outline="black", tags="rect")
    cv._color = color  # type: ignore[attr-defined]
    return cv


def _update_swatch(cv: tk.Canvas, color: str) -> None:
    """Repaint *cv* with the new *color*."""
    cv.itemconfig("rect", fill=color)
    cv._color = color  # type: ignore[attr-defined]


def _pick_color(parent, current: str, title: str = "Choose color") -> Optional[str]:
    """Open the system color-chooser and return the chosen hex, or None."""
    _, hx = colorchooser.askcolor(color=current, parent=parent, title=title)
    return hx or None


# ---------------------------------------------------------------------------
# Main dialog class
# ---------------------------------------------------------------------------

class ConfigDialog:
    """Tabbed *App Settings* dialog.

    Parameters
    ----------
    parent:
        Parent Tk window.
    cfg:
        Current :class:`~config.AppConfig` (will *not* be modified in-place;
        a copy is used internally).
    on_save:
        Callback invoked with the updated :class:`AppConfig` after the user
        clicks **OK** and the new config has been saved to disk.
    """

    def __init__(
        self,
        parent: tk.Tk,
        cfg: AppConfig,
        on_save: Callable[[AppConfig], None],
    ) -> None:
        self._parent = parent
        # Work on a deep copy so Cancel truly discards all changes
        self._cfg = copy.deepcopy(cfg)
        self._on_save = on_save

        self._top = ttk.Toplevel(parent)
        self._top.title("App Settings")
        self._top.transient(parent)
        self._top.grab_set()
        self._top.resizable(True, True)
        self._top.geometry(
            f"+{parent.winfo_x() + 60}+{parent.winfo_y() + 60}"
        )

        self._build()

    # ------------------------------------------------------------------
    # Top-level layout
    # ------------------------------------------------------------------

    def _build(self) -> None:
        outer = ttk.Frame(self._top, padding=10)
        outer.pack(fill=BOTH, expand=True)

        nb = ttk.Notebook(outer)
        nb.pack(fill=BOTH, expand=True)

        nb.add(self._build_app_tab(nb),      text="App & Window")
        nb.add(self._build_defaults_tab(nb), text="Defaults")
        nb.add(self._build_plot_tab(nb),     text="Plot Style")
        nb.add(self._build_colors_tab(nb),   text="Colors")

        sep = ttk.Separator(outer, orient=HORIZONTAL)
        sep.pack(fill=X, pady=(10, 5))

        btn_row = ttk.Frame(outer)
        btn_row.pack(fill=X)

        ttk.Button(btn_row, text="Reset to Defaults",
                   command=self._on_reset).pack(side=LEFT)
        ttk.Button(btn_row, text="Cancel",
                   command=self._top.destroy).pack(side=RIGHT, padx=(5, 0))
        ttk.Button(btn_row, text="OK", style="Accent.TButton",
                   command=self._on_ok).pack(side=RIGHT)

        self._top.bind("<Return>", lambda _e: self._on_ok())
        self._top.bind("<Escape>", lambda _e: self._top.destroy())

    # ------------------------------------------------------------------
    # Tab: App & Window
    # ------------------------------------------------------------------

    def _build_app_tab(self, parent) -> ttk.Frame:
        tab = ttk.Frame(parent, padding=12)

        # restart notice
        notice = ttk.Label(
            tab,
            text="Note:  Theme and window-size changes take effect on next launch.",
            bootstyle="info",
            wraplength=440,
        )
        notice.grid(row=0, column=0, columnspan=3, sticky=W, pady=(0, 10))

        # Theme
        ttk.Label(tab, text="Theme:").grid(row=1, column=0, sticky=E, padx=(0, 8), pady=5)
        self._v_theme = ttk.StringVar(value=self._cfg.app.theme)
        ttk.Combobox(tab, textvariable=self._v_theme, values=_THEMES,
                     state="readonly", width=20).grid(row=1, column=1, sticky=W)

        # Window width
        ttk.Label(tab, text="Window width (px):").grid(row=2, column=0, sticky=E, padx=(0, 8), pady=5)
        self._v_win_w = ttk.IntVar(value=self._cfg.window.width)
        ttk.Entry(tab, textvariable=self._v_win_w, width=8).grid(row=2, column=1, sticky=W)

        # Window height
        ttk.Label(tab, text="Window height (px):").grid(row=3, column=0, sticky=E, padx=(0, 8), pady=5)
        self._v_win_h = ttk.IntVar(value=self._cfg.window.height)
        ttk.Entry(tab, textvariable=self._v_win_h, width=8).grid(row=3, column=1, sticky=W)

        tab.columnconfigure(2, weight=1)
        return tab

    # ------------------------------------------------------------------
    # Tab: Defaults
    # ------------------------------------------------------------------

    def _build_defaults_tab(self, parent) -> ttk.Frame:
        tab = ttk.Frame(parent, padding=12)
        d = self._cfg.defaults
        r = 0

        def _row(label, widget_factory, *, pady=4):
            nonlocal r
            ttk.Label(tab, text=label).grid(row=r, column=0, sticky=E, padx=(0, 8), pady=pady)
            w = widget_factory(tab)
            w.grid(row=r, column=1, sticky=W, pady=pady)
            r += 1
            return w

        # Mode
        self._v_mode = ttk.StringVar(value=d.mode)
        _row("Photometry mode:",
             lambda p: ttk.Combobox(p, textvariable=self._v_mode,
                                    values=_MODES, state="readonly", width=18))

        # Time mode
        self._v_time_mode = ttk.StringVar(value=d.time_mode)
        _row("Time axis:",
             lambda p: ttk.Combobox(p, textvariable=self._v_time_mode,
                                    values=_TIME_MODES, state="readonly", width=18))

        # Error bar type
        self._v_errorbar = ttk.StringVar(value=d.errorbar_type)
        _row("Error bar type:",
             lambda p: ttk.Combobox(p, textvariable=self._v_errorbar,
                                    values=_ERR_TYPES, state="readonly", width=18))

        # Color legend location
        self._v_legend_loc = ttk.StringVar(value=d.color_legend_loc)
        _row("Color-info legend location:",
             lambda p: ttk.Combobox(p, textvariable=self._v_legend_loc,
                                    values=_LEGEND_LOCS, state="readonly", width=18))

        # Numeric entries ------------------------------------------------
        self._v_period      = ttk.DoubleVar(value=d.period)
        self._v_period_step = ttk.DoubleVar(value=d.period_step)
        self._v_phase_max   = ttk.DoubleVar(value=d.phase_max)
        self._v_offset_step = ttk.DoubleVar(value=d.offset_step)

        _row("Default period (h):",     lambda p: ttk.Entry(p, textvariable=self._v_period,      width=10))
        _row("Period step (h):",        lambda p: ttk.Entry(p, textvariable=self._v_period_step, width=10))
        _row("Phase max:",              lambda p: ttk.Entry(p, textvariable=self._v_phase_max,   width=10))
        _row("Offset step (mag):",      lambda p: ttk.Entry(p, textvariable=self._v_offset_step, width=10))

        # Boolean checkboxes ---------------------------------------------
        self._v_show_rej    = ttk.BooleanVar(value=d.show_rejected)
        self._v_show_legend = ttk.BooleanVar(value=d.show_legend)
        self._v_show_grid   = ttk.BooleanVar(value=d.show_grid)

        for label, var in [
            ("Show rejected points by default", self._v_show_rej),
            ("Show legend by default",          self._v_show_legend),
            ("Show grid by default",            self._v_show_grid),
        ]:
            ttk.Checkbutton(tab, text=label, variable=var).grid(
                row=r, column=0, columnspan=2, sticky=W, pady=3)
            r += 1

        tab.columnconfigure(2, weight=1)
        return tab

    # ------------------------------------------------------------------
    # Tab: Plot Style
    # ------------------------------------------------------------------

    def _build_plot_tab(self, parent) -> ttk.Frame:
        tab = ttk.Frame(parent, padding=12)
        p = self._cfg.plot
        r = 0

        def _label_row(text, *, pady=4):
            nonlocal r
            ttk.Label(tab, text=text, bootstyle="secondary").grid(
                row=r, column=0, columnspan=3, sticky=W, pady=(pady, 0))
            r += 1

        def _entry_row(label, var, *, pady=4):
            nonlocal r
            ttk.Label(tab, text=label).grid(row=r, column=0, sticky=E, padx=(0, 8), pady=pady)
            ttk.Entry(tab, textvariable=var, width=8).grid(row=r, column=1, sticky=W, pady=pady)
            r += 1

        def _scale_row(label, var, from_, to, *, pady=4):
            nonlocal r
            ttk.Label(tab, text=label).grid(row=r, column=0, sticky=E, padx=(0, 8), pady=pady)
            ttk.Scale(tab, from_=from_, to=to, variable=var,
                      orient=HORIZONTAL, length=130).grid(row=r, column=1, sticky=W, pady=pady)
            ttk.Entry(tab, textvariable=var, width=6).grid(row=r, column=2, sticky=W, pady=pady)
            r += 1

        # Figure size
        _label_row("Figure")
        self._v_fig_w = ttk.DoubleVar(value=p.figsize[0])
        self._v_fig_h = ttk.DoubleVar(value=p.figsize[1])
        _entry_row("Figure width (in):",  self._v_fig_w)
        _entry_row("Figure height (in):", self._v_fig_h)

        # Marker
        _label_row("Markers", pady=8)
        self._v_marker_style = ttk.StringVar(
            value=_MARKER_NAMES_MAP.get(p.marker_style, "circle"))
        ttk.Label(tab, text="Marker style:").grid(row=r, column=0, sticky=E, padx=(0, 8), pady=4)
        ttk.Combobox(tab, textvariable=self._v_marker_style, values=_MARKER_NAMES,
                     state="readonly", width=18).grid(row=r, column=1, sticky=W, pady=4)
        r += 1

        self._v_marker_size = ttk.DoubleVar(value=p.marker_size)
        _scale_row("Marker size (pt):", self._v_marker_size, 1, 20)

        # Error bars
        _label_row("Error bars", pady=8)
        self._v_cap_size   = ttk.DoubleVar(value=p.errorbar_capsize)
        self._v_cap_thick  = ttk.DoubleVar(value=p.errorbar_capthick)
        self._v_lw         = ttk.DoubleVar(value=p.errorbar_linewidth)
        _scale_row("Cap size (pt):",   self._v_cap_size,  0, 20)
        _scale_row("Cap thickness:",   self._v_cap_thick, 0, 10)
        _scale_row("Bar line width:",  self._v_lw,        0, 10)

        # Ghost / phase wrap
        _label_row("Phase-wrap ghost", pady=8)
        self._v_ghost_alpha = ttk.DoubleVar(value=p.ghost_alpha)
        self._v_ghost_ew    = ttk.DoubleVar(value=p.ghost_edge_width)
        _scale_row("Ghost alpha:",      self._v_ghost_alpha, 0.0, 1.0)
        _scale_row("Ghost edge width:", self._v_ghost_ew,    0.1, 5.0)

        # Phase wrap line
        _label_row("Phase-wrap boundary line", pady=8)
        self._v_wrap_lw    = ttk.DoubleVar(value=p.phase_wrap_linewidth)
        self._v_wrap_alpha = ttk.DoubleVar(value=p.phase_wrap_alpha)
        _scale_row("Line width:", self._v_wrap_lw,    0.1, 5.0)
        _scale_row("Alpha:",      self._v_wrap_alpha, 0.0, 1.0)

        # Phase wrap color swatch
        ttk.Label(tab, text="Line color:").grid(row=r, column=0, sticky=E, padx=(0, 8), pady=4)
        self._wrap_swatch = _make_swatch(tab, p.phase_wrap_color)
        self._wrap_swatch.grid(row=r, column=1, sticky=W)
        self._wrap_swatch.bind(
            "<Button-1>",
            lambda _e: self._pick_and_update(self._wrap_swatch, "Phase-wrap line color"),
        )
        r += 1

        # Grid
        _label_row("Grid", pady=8)
        self._v_grid_alpha = ttk.DoubleVar(value=p.grid_alpha)
        _scale_row("Grid alpha:", self._v_grid_alpha, 0.0, 1.0)

        # Fonts
        _label_row("Fonts", pady=8)
        self._v_legend_fs = ttk.StringVar(value=p.legend_fontsize)
        ttk.Label(tab, text="Legend font size:").grid(row=r, column=0, sticky=E, padx=(0, 8), pady=4)
        ttk.Combobox(tab, textvariable=self._v_legend_fs, values=_FONT_SIZES,
                     state="readonly", width=12).grid(row=r, column=1, sticky=W, pady=4)
        r += 1

        self._v_label_fs = ttk.IntVar(value=p.label_fontsize)
        _entry_row("Info label font size (pt):", self._v_label_fs)

        tab.columnconfigure(3, weight=1)
        return tab

    # ------------------------------------------------------------------
    # Tab: Colors
    # ------------------------------------------------------------------

    def _build_colors_tab(self, parent) -> ttk.Frame:
        tab = ttk.Frame(parent, padding=12)
        c = self._cfg.colors
        r = 0

        def _section(text):
            nonlocal r
            ttk.Label(tab, text=text, bootstyle="secondary").grid(
                row=r, column=0, columnspan=4, sticky=W, pady=(8, 2))
            ttk.Separator(tab, orient=HORIZONTAL).grid(
                row=r + 1, column=0, columnspan=4, sticky=EW, pady=(0, 4))
            r += 2

        def _color_row(label, color, swatch_store_key):
            nonlocal r
            ttk.Label(tab, text=label).grid(row=r, column=0, sticky=E, padx=(0, 8), pady=4)
            sw = _make_swatch(tab, color)
            sw.grid(row=r, column=1, sticky=W, padx=(0, 6), pady=4)
            btn = ttk.Button(
                tab, text="Pick...", width=6,
                command=lambda _sw=sw, _lbl=label: self._pick_and_update(_sw, _lbl),
            )
            btn.grid(row=r, column=2, sticky=W, pady=4)
            r += 1
            return sw

        # Flagged / rejected
        _section("Special colors")
        self._swatch_flagged  = _color_row("Flagged color:",  c.flagged,  "flagged")
        self._swatch_rejected = _color_row("Rejected color:", c.rejected, "rejected")

        self._v_flagged_use_filter = ttk.BooleanVar(value=c.flagged_use_filter)
        ttk.Checkbutton(
            tab,
            text="Flagged points use their filter band color",
            variable=self._v_flagged_use_filter,
        ).grid(row=r, column=0, columnspan=4, sticky=W, pady=4)
        r += 1

        # Default palette
        _section("Default color palette")
        palette = list(c.default_palette)
        # pad/trim to 10 entries
        while len(palette) < 10:
            palette.append("#888888")
        palette = palette[:10]

        self._palette_swatches: List[tk.Canvas] = []
        pal_frame = ttk.Frame(tab)
        pal_frame.grid(row=r, column=0, columnspan=4, sticky=W, pady=4)
        r += 1

        for idx, clr in enumerate(palette):
            sw = _make_swatch(pal_frame, clr, size=(34, 22))
            sw.grid(row=0, column=idx, padx=2)
            sw.bind(
                "<Button-1>",
                lambda _e, _i=idx, _sw=sw: self._pick_palette_entry(_i, _sw),
            )
            ttk.Label(pal_frame, text=str(idx + 1), font=("TkDefaultFont", 7)).grid(
                row=1, column=idx)
            self._palette_swatches.append(sw)

        # Band colors
        _section("Band colors")

        band_scroll = ttk.Frame(tab)
        band_scroll.grid(row=r, column=0, columnspan=4, sticky=NSEW, pady=4)
        tab.rowconfigure(r, weight=1)
        r += 1

        self._band_swatches: Dict[str, tk.Canvas] = {}
        for band_r, (band, color) in enumerate(sorted(c.band.items())):
            ttk.Label(band_scroll, text=band, width=4, anchor=E).grid(
                row=band_r, column=0, padx=(0, 6), pady=2)
            sw = _make_swatch(band_scroll, color)
            sw.grid(row=band_r, column=1, padx=(0, 6), pady=2)
            btn = ttk.Button(
                band_scroll, text="Pick...", width=6,
                command=lambda _b=band, _sw=sw: self._pick_and_update(
                    _sw, f"Color for band {_b}"),
            )
            btn.grid(row=band_r, column=2, pady=2)
            self._band_swatches[band] = sw

        tab.columnconfigure(3, weight=1)
        return tab

    # ------------------------------------------------------------------
    # Color-picking helpers
    # ------------------------------------------------------------------

    def _pick_and_update(self, swatch: tk.Canvas, title: str) -> None:
        new_color = _pick_color(self._top, swatch._color, title)  # type: ignore[attr-defined]
        if new_color:
            _update_swatch(swatch, new_color)

    def _pick_palette_entry(self, idx: int, swatch: tk.Canvas) -> None:
        new_color = _pick_color(self._top, swatch._color, f"Palette color #{idx + 1}")  # type: ignore[attr-defined]
        if new_color:
            _update_swatch(swatch, new_color)

    # ------------------------------------------------------------------
    # Button handlers
    # ------------------------------------------------------------------

    def _collect(self) -> AppConfig:
        """Read all widget vars back into a new AppConfig."""
        import copy
        cfg = copy.deepcopy(self._cfg)

        # App & Window
        cfg.app.theme   = self._v_theme.get()
        cfg.window.width  = int(self._v_win_w.get())
        cfg.window.height = int(self._v_win_h.get())

        # Defaults
        cfg.defaults.mode             = self._v_mode.get()
        cfg.defaults.time_mode        = self._v_time_mode.get()
        cfg.defaults.errorbar_type    = self._v_errorbar.get()
        cfg.defaults.color_legend_loc = self._v_legend_loc.get()
        cfg.defaults.period           = float(self._v_period.get())
        cfg.defaults.period_step      = float(self._v_period_step.get())
        cfg.defaults.phase_max        = float(self._v_phase_max.get())
        cfg.defaults.offset_step      = float(self._v_offset_step.get())
        cfg.defaults.show_rejected    = bool(self._v_show_rej.get())
        cfg.defaults.show_legend      = bool(self._v_show_legend.get())
        cfg.defaults.show_grid        = bool(self._v_show_grid.get())

        # Plot Style
        cfg.plot.figsize           = [float(self._v_fig_w.get()), float(self._v_fig_h.get())]
        cfg.plot.marker_style      = _MARKER_MAP.get(self._v_marker_style.get(), "o")
        cfg.plot.marker_size       = float(self._v_marker_size.get())
        cfg.plot.errorbar_capsize  = float(self._v_cap_size.get())
        cfg.plot.errorbar_capthick = float(self._v_cap_thick.get())
        cfg.plot.errorbar_linewidth = float(self._v_lw.get())
        cfg.plot.ghost_alpha       = float(self._v_ghost_alpha.get())
        cfg.plot.ghost_edge_width  = float(self._v_ghost_ew.get())
        cfg.plot.phase_wrap_color  = self._wrap_swatch._color  # type: ignore[attr-defined]
        cfg.plot.phase_wrap_linewidth = float(self._v_wrap_lw.get())
        cfg.plot.phase_wrap_alpha  = float(self._v_wrap_alpha.get())
        cfg.plot.grid_alpha        = float(self._v_grid_alpha.get())
        cfg.plot.legend_fontsize   = self._v_legend_fs.get()
        cfg.plot.label_fontsize    = int(self._v_label_fs.get())

        # Colors
        cfg.colors.flagged  = self._swatch_flagged._color   # type: ignore[attr-defined]
        cfg.colors.rejected = self._swatch_rejected._color  # type: ignore[attr-defined]
        cfg.colors.flagged_use_filter = bool(self._v_flagged_use_filter.get())
        cfg.colors.default_palette = [sw._color for sw in self._palette_swatches]  # type: ignore[attr-defined]
        for band, sw in self._band_swatches.items():
            cfg.colors.band[band] = sw._color  # type: ignore[attr-defined]

        return cfg

    def _on_ok(self) -> None:
        try:
            new_cfg = self._collect()
        except (ValueError, TypeError) as exc:
            import tkinter.messagebox as mb
            mb.showerror("Invalid value", str(exc), parent=self._top)
            return
        save_user_config(new_cfg)
        self._on_save(new_cfg)
        self._top.destroy()

    def _on_reset(self) -> None:
        """Reload defaults from the bundled TOML and repopulate all widgets."""
        from .config import load_config, _DEFAULT_TOML
        try:
            import tomli as _tomllib
        except ImportError:
            try:
                import tomllib as _tomllib  # type: ignore[no-redef]
            except ImportError:
                return
        with open(_DEFAULT_TOML, "rb") as fh:
            d = _tomllib.load(fh)
        from .config import _dict_to_config
        default_cfg = _dict_to_config(d)

        # Repopulate - close and re-open
        self._cfg = default_cfg
        self._top.destroy()
        ConfigDialog(self._parent, default_cfg, self._on_save)


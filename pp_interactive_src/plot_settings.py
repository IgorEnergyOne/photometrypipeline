# -*- coding: utf-8 -*-
"""
PlotSettings: manages plot style constants, per-band colour overrides, and the
Colors menu for pp_interactive_src.
"""
from __future__ import annotations

from typing import Callable, Dict, List, Optional

import tkinter as tk
import ttkbootstrap as ttk
from ttkbootstrap.constants import *
from tkinter import colorchooser

from .utils import _safe


class PlotSettings:
    """Manages plot settings, constants, and color/palette menus."""

    # Class-level fallbacks (used when no AppConfig is supplied)
    DEFAULT_COLORS = [
        "red", "orange", "olive", "green", "blue", "purple",
        "brown", "pink", "gray", "cyan"
    ]

    BAND_COLORS = {
        "U": "indigo", "B": "royalblue", "V": "limegreen",
        "R": '#7c3150', "I": "dimgray",
        "g": "#348034", "r": "#944d4d", "i": "#59327d", "z": "#753427",
        'G': 'green', 'BP': 'blue', 'RP': 'red',
    }

    SELECTION_MARKER = {'markersize': 5, 'zorder': 20, 'form': 'o', 'color': 'red'}

    MARKERS = [
        ('o', 'circle'), ('s', 'square'), ('p', 'pentagon'), ('x', 'x'), ('D', 'diamond'),
        ('*', 'star'), ('v', 'triangle_down'), ('^', 'triangle_up'), ('<', 'triangle_left'),
        ('>', 'triangle_right'), ('+', 'plus'), ('d', 'thin_diamond'),
    ]

    PALETTE_SWATCHES = DEFAULT_COLORS + list(dict.fromkeys([v for v in BAND_COLORS.values()])) + [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
    ]

    def __init__(
        self,
        root: tk.Tk,
        refresh_callback: Callable[[], None],
        get_bands_callback: Callable[[], List[str]],
        cfg=None,  # Optional[AppConfig]
    ):
        self.root = root
        self.refresh_callback = refresh_callback
        self.get_bands_callback = get_bands_callback

        # Apply config values as instance attributes (shadow class-level defaults)
        if cfg is not None:
            self.DEFAULT_COLORS = list(cfg.colors.default_palette)
            self.BAND_COLORS    = dict(cfg.colors.band)
            self.custom_colors: Dict[str, str] = {
                'flagged':  cfg.colors.flagged,
                'rejected': cfg.colors.rejected,
            }
            self.flagged_use_filter: bool = cfg.colors.flagged_use_filter
        else:
            self.custom_colors = {'flagged': 'orange', 'rejected': 'red'}
            self.flagged_use_filter = False

        self.PALETTE_SWATCHES = self._build_palette_swatches()

        self.color_menu: Optional[tk.Menu] = None
        self.color_menu_btn: Optional[ttk.Menubutton] = None

    #  -  -  -  palette helper  -  -  - 
    def _build_palette_swatches(self) -> List[str]:
        """Build the palette swatch list from current DEFAULT_COLORS and BAND_COLORS."""
        extra = [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
        ]
        return list(dict.fromkeys(
            self.DEFAULT_COLORS
            + list(self.BAND_COLORS.values())
            + extra
        ))

    def apply_config(self, cfg) -> None:
        """Re-apply settings from a new AppConfig (called after config is saved)."""
        self.DEFAULT_COLORS   = list(cfg.colors.default_palette)
        self.BAND_COLORS      = dict(cfg.colors.band)
        self.custom_colors['flagged']  = cfg.colors.flagged
        self.custom_colors['rejected'] = cfg.colors.rejected
        self.flagged_use_filter        = cfg.colors.flagged_use_filter
        self.PALETTE_SWATCHES = self._build_palette_swatches()

    #  -  -  -  colour resolution  -  -  - 
    def get_color(self, band: Optional[str], idx: int) -> str:
        """Resolve display colour for *band*, checking custom overrides first."""
        if band is None:
            return self.DEFAULT_COLORS[idx % len(self.DEFAULT_COLORS)]
        if str(band) in self.custom_colors:
            return self.custom_colors[str(band)]
        if band in self.BAND_COLORS:
            return self.BAND_COLORS[band]
        return self.DEFAULT_COLORS[idx % len(self.DEFAULT_COLORS)]

    #  -  -  -  menu construction  -  -  - 
    def build_menu(self, parent=None, pack_btn: bool = True):
        """Create the Colors menubutton (dynamically populated via rebuild_menu)."""
        master_for_menu = parent if parent is not None else self.root
        try:
            if hasattr(master_for_menu, 'add_command') and isinstance(master_for_menu, tk.Menu):
                menu_master = master_for_menu
            else:
                menu_master = self.root
        except Exception:
            menu_master = self.root

        self.color_menu = tk.Menu(menu_master, tearoff=0, postcommand=self.rebuild_menu)
        if pack_btn and parent is not None:
            self.color_menu_btn = ttk.Menubutton(parent, text="Colors")
            self.color_menu_btn["menu"] = self.color_menu
            self.color_menu_btn.pack(side=LEFT, padx=4)

        self.rebuild_menu()

    def rebuild_menu(self):
        """Populate the Colors menu with per-band entries and special controls."""
        if not hasattr(self, 'color_menu') or self.color_menu is None:
            return
        try:
            self.color_menu.delete(0, END)
        except Exception:
            return

        bands = self.get_bands_callback()
        if bands:
            self.color_menu.add_command(label="Apply BAND_COLORS defaults",
                                        command=lambda: self.apply_default_colors(True))
            self.color_menu.add_command(label="Apply simple palette defaults",
                                        command=lambda: self.apply_default_colors(False))
            self.color_menu.add_separator()

            for b in bands:
                sub = tk.Menu(self.color_menu, tearoff=0)
                sub.add_command(label="Pick color...", command=lambda _b=b: self.pick_color_for_band(_b))
                if b in self.BAND_COLORS:
                    sub.add_command(
                        label=f"Use canonical ({self.BAND_COLORS[b]})",
                        command=lambda _b=b, _c=self.BAND_COLORS[b]: self.apply_color_to_band(_b, _c),
                    )
                sub.add_command(label="Palette...", command=lambda _b=b: self.show_palette_picker(_b))
                self.color_menu.add_cascade(label=b, menu=sub)
            self.color_menu.add_separator()

        def _toggle_flagged_use():
            self.flagged_use_filter = not self.flagged_use_filter
            self.refresh_callback()
            self.rebuild_menu()

        chk_label = "Flagged use filter color"
        if self.flagged_use_filter:
            self.color_menu.add_command(label=chk_label + " [x]", command=_toggle_flagged_use)
        else:
            self.color_menu.add_command(label=chk_label, command=_toggle_flagged_use)

        self.color_menu.add_separator()
        self.color_menu.add_command(label="Set flagged color...",
                                    command=lambda: self.show_palette_picker('flagged'))
        self.color_menu.add_command(label="Set rejected color...",
                                    command=lambda: self.show_palette_picker('rejected'))

    #  -  -  -  colour application helpers  -  -  - 
    def apply_color_to_band(self, band: str, hexcolor: str):
        if not band:
            return
        try:
            self.custom_colors[str(band)] = str(hexcolor)
        except Exception:
            self.custom_colors[str(band)] = hexcolor
        self.refresh_callback()
        self.rebuild_menu()

    def apply_default_colors(self, use_band_colors: bool = True):
        bands = self.get_bands_callback()
        for i, b in enumerate(bands):
            if use_band_colors and b in self.BAND_COLORS:
                self.custom_colors[b] = self.BAND_COLORS[b]
            else:
                self.custom_colors[b] = self.DEFAULT_COLORS[i % len(self.DEFAULT_COLORS)]
        self.refresh_callback()
        self.rebuild_menu()

    def pick_color_for_band(self, band: str):
        try:
            cur = self.custom_colors.get(band, self.BAND_COLORS.get(band, '#000000'))
            rgb, hx = colorchooser.askcolor(color=cur, parent=self.root, title=f"Choose color for {band}")
            if hx:
                self.custom_colors[str(band)] = hx
                self.refresh_callback()
                self.rebuild_menu()
        except Exception:
            pass

    def show_palette_picker(self, band: Optional[str] = None, title: Optional[str] = None):
        if title is None:
            title = f"Choose color for {band}" if band else "Choose color"
        top = tk.Toplevel(self.root)
        top.transient(self.root)
        top.title(title)
        top.resizable(False, False)
        _safe(top.geometry, f"+{self.root.winfo_x() + 120}+{self.root.winfo_y() + 120}")

        frm = ttk.Frame(top, padding=6)
        frm.pack(fill='both', expand=True)

        cols = 8
        colors = list(dict.fromkeys(self.PALETTE_SWATCHES))
        sw_w, sw_h = 28, 18
        for i, c in enumerate(colors):
            r = i // cols
            col = i % cols
            sw = tk.Canvas(frm, width=sw_w, height=sw_h, highlightthickness=1, bd=0)
            sw.grid(row=r, column=col, padx=3, pady=3)
            rect = sw.create_rectangle(0, 0, sw_w, sw_h, fill=c, outline='black')
            sw.tag_bind(rect, "<Button-1>", lambda e, _c=c: (self.apply_color_to_band(band, _c), top.destroy()))
            sw.bind("<Button-1>", lambda e, _c=c: (self.apply_color_to_band(band, _c), top.destroy()))

        btn_row = (len(colors) + cols - 1) // cols
        ttk.Button(frm, text="More...",
                   command=lambda: (self.pick_color_for_band(band), top.destroy())).grid(
            row=btn_row, column=0, columnspan=cols, pady=(8, 0))

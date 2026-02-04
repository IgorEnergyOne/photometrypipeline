#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Interactive light‑curve viewer & editor (multi‑filter edition)

Key additions in this refactor:
1) Multi‑filter visibility menu with checkboxes (show any subset of filters).
2) Per‑filter vertical offsets per file (adjust Up/Down for the active scope).
3) Robust selection marker (always visible, consistent X transform) and
   fixes to accidental Y inversion on redraws.
4) Leaner structure and helper methods; safe, fast redraws via blitting.

Hotkeys:
  q – quit; r – toggle rejection; a – cancel selection
  S – show/hide image for selected point
  ←/→ – move selection point; ↑/↓ – nudge offset for active scope
  z/x – −/+ rotation period (when in rotation phase mode)

"""
from __future__ import annotations

import os
import re
import argparse
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Iterable, Set, Callable
import subprocess

import numpy as np
import pandas as pd
from astropy.time import Time

import matplotlib
import matplotlib.pyplot as plt

# Prevent the Matplotlib navigation toolbar from binding the 's' key to Save
# (the default rcParam 'keymap.save' contains ['s']). Clear it so plain 's' no longer opens the Save dialog.
matplotlib.rcParams['keymap.save'] = []
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import tkinter as tk
import ttkbootstrap as ttk
from ttkbootstrap.constants import *
from tkinter import filedialog, messagebox, simpledialog, colorchooser

# Import Pillow for image handling. Add a check in case it's not installed.
try:
    from PIL import Image, ImageTk
except ImportError:
    Image = None
    ImageTk = None
    print("Warning: Pillow library not found. Image display feature will not work.")
    print("Please install it using: pip install Pillow")

# ────────────────────────────── Globals & style ─────────────────────────────
# Global debug flag
DEBUG = False
TIME_STEP = 0.02  # hours
DEFAULT_PERIOD = 4.0  # hours
WINDOW_WIDTH = 1400  # px
WINDOW_HEIGHT = 720  # pix
OFFSET_ALL_KEY = "__ALL__"  # per‑file fallback offset when no per‑band offset is set


class PlotSettings:
    """
    Manages plot settings, constants, and color/palette menus.
    """
    DEFAULT_COLORS = [
        "red", "orange", "olive", "green", "blue", "purple",
        "brown", "pink", "gray", "cyan"
    ]
    BAND_COLORS = {"U": "indigo", "B": "royalblue", "V": "limegreen", "R": '#7c3150', "I": "dimgray",
                   "g": "#348034", "r": "#944d4d", "i": "#59327d", "z": "#753427"}

    SELECTION_MARKER = {'markersize': 5, 'zorder': 20, 'form': 'o', 'color': 'red'}

    MARKERS = [
        ('o', 'circle'), ('s', 'square'), ('p', 'pentagon'), ('x', 'x'), ('D', 'diamond'),
        ('*', 'star'), ('v', 'triangle_down'), ('^', 'triangle_up'), ('<', 'triangle_left'),
        ('>', 'triangle_right'), ('+', 'plus'), ('d', 'thin_diamond'),
    ]

    PALETTE_SWATCHES = DEFAULT_COLORS + list(dict.fromkeys([v for v in BAND_COLORS.values()])) + [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
    ]

    def __init__(self, root: tk.Tk, refresh_callback: Callable[[], None], get_bands_callback: Callable[[], List[str]]):
        self.root = root
        self.refresh_callback = refresh_callback
        self.get_bands_callback = get_bands_callback

        self.custom_colors: Dict[str, str] = {'flagged': 'orange', 'rejected': 'red'}
        self.flagged_use_filter: bool = False

        self.color_menu: Optional[tk.Menu] = None
        self.color_menu_btn: Optional[ttk.Menubutton] = None

    def get_color(self, band: Optional[str], idx: int) -> str:
        """Resolve color for a band, checking custom overrides first."""
        if band is None:
            return self.DEFAULT_COLORS[idx % len(self.DEFAULT_COLORS)]

        # Check custom overrides
        if str(band) in self.custom_colors:
            return self.custom_colors[str(band)]

        # Check canonical band colors
        if band in self.BAND_COLORS:
            return self.BAND_COLORS[band]

        return self.DEFAULT_COLORS[idx % len(self.DEFAULT_COLORS)]

    def build_menu(self, parent=None, pack_btn=True):
        """Create Colors menubutton (dynamically populated via rebuild_menu)."""
        master_for_menu = parent if parent is not None else self.root
        try:
            if hasattr(master_for_menu, 'add_command') and isinstance(master_for_menu, tk.Menu):
                menu_master = master_for_menu
            else:
                menu_master = self.root
        except Exception:
            menu_master = self.root

        self.color_menu = tk.Menu(menu_master, tearoff=0)
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
                    sub.add_command(label=f"Use canonical ({self.BAND_COLORS[b]})",
                                    command=lambda _b=b, _c=self.BAND_COLORS[b]: self.apply_color_to_band(_b, _c))
                sub.add_command(label="Palette...", command=lambda _b=b: self.show_palette_picker(_b))
                self.color_menu.add_cascade(label=b, menu=sub)
            self.color_menu.add_separator()

        def _toggle_flagged_use():
            self.flagged_use_filter = not self.flagged_use_filter
            self.refresh_callback()
            self.rebuild_menu()

        chk_label = "Flagged use filter color"
        if self.flagged_use_filter:
            self.color_menu.add_command(label=chk_label + " ✓", command=_toggle_flagged_use)
        else:
            self.color_menu.add_command(label=chk_label, command=_toggle_flagged_use)

        self.color_menu.add_separator()
        self.color_menu.add_command(label="Set flagged color...",
                                    command=lambda: self.show_palette_picker('flagged'))
        self.color_menu.add_command(label="Set rejected color...",
                                    command=lambda: self.show_palette_picker('rejected'))

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
        try:
            top = tk.Toplevel(self.root)
            top.transient(self.root)
            top.title(title)
            top.resizable(False, False)
            try:
                x = self.root.winfo_x() + 120
                y = self.root.winfo_y() + 120
                top.geometry(f"+{x}+{y}")
            except Exception:
                pass

            frm = ttk.Frame(top, padding=6)
            frm.pack(fill='both', expand=True)

            cols = 8
            colors = list(dict.fromkeys(self.PALETTE_SWATCHES))
            sw_w = 28
            sw_h = 18
            for i, c in enumerate(colors):
                r = i // cols;
                col = i % cols
                sw = tk.Canvas(frm, width=sw_w, height=sw_h, highlightthickness=1, bd=0)
                sw.grid(row=r, column=col, padx=3, pady=3)
                rect = sw.create_rectangle(0, 0, sw_w, sw_h, fill=c, outline='black')
                sw.tag_bind(rect, "<Button-1>", lambda e, _c=c: (self.apply_color_to_band(band, _c), top.destroy()))
                sw.bind("<Button-1>", lambda e, _c=c: (self.apply_color_to_band(band, _c), top.destroy()))

            btn_row = (len(colors) + cols - 1) // cols
            ttk.Button(frm, text="More...", command=lambda: (self.pick_color_for_band(band), top.destroy())).grid(
                row=btn_row, column=0, columnspan=cols, pady=(8, 0))
        except Exception:
            pass


def debug_print(*args, **kwargs):
    """Print debug messages only if DEBUG is True"""
    if DEBUG:
        print("[DEBUG]", *args, **kwargs)


def next_version(path: str) -> str:
    """Generate the next version of a file by incrementing the counter in the filename."""
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

# ────────────────────────────── Data layer ──────────────────────────────
class LightCurveData:
    """Container for one CSV light‑curve plus cached arrays for speed."""

    def __init__(self, filepath: Optional[str] = None):
        self.filename: Optional[str] = None
        self.df: Optional[pd.DataFrame] = None
        # cached numpy columns for fast plotting
        self.arr: Dict[str, np.ndarray] = {}
        if filepath:
            self.load(filepath)

    # ——— I/O ———
    def load(self, filepath: str) -> None:
        self.filename = filepath
        self.df = pd.read_csv(filepath)
        self._ensure_columns()
        self._fill_arrays_cache()

    def save(self) -> None:
        if self.df is not None and self.filename:
            self.df.to_csv(self.filename, index=False)

    # ——— internals ———
    def _ensure_columns(self):
        assert self.df is not None
        if 'rejected' not in self.df.columns:
            self.df['rejected'] = False
        if 'sextractor_flags' not in self.df.columns:
            self.df['sextractor_flags'] = 0
        if 'filename' not in self.df.columns:
            self.df['filename'] = os.path.basename(self.filename) if self.filename else "lightcurve.csv"

    def _fill_arrays_cache(self):
        """Fill the cached numpy arrays for fast plotting."""
        assert self.df is not None
        df = self.df
        self.arr['jd'] = df.get('julian_date', pd.Series(dtype=float)).to_numpy(dtype=float, copy=False)
        self.arr['mag'] = df.get('mag', pd.Series(dtype=float)).to_numpy(dtype=float, copy=False)
        self.arr['inst_mag'] = df.get('inst_mag', pd.Series(dtype=float)).to_numpy(dtype=float, copy=False)
        self.arr['mag_control'] = df.get('mag_control', pd.Series(dtype=float)).to_numpy(dtype=float, copy=False)
        self.arr['sig'] = df.get('sig', pd.Series(dtype=float)).to_numpy(dtype=float, copy=False)
        self.arr['inst_sig'] = df.get('inst_sig', pd.Series(dtype=float)).to_numpy(dtype=float, copy=False)
        self.arr['sig_control'] = df.get('sig_control', pd.Series(dtype=float)).to_numpy(dtype=float, copy=False)
        self.arr['rel_mag'] = df.get('rel_mag', pd.Series(dtype=float)).to_numpy(dtype=float, copy=False)
        self.arr['rel_sig'] = df.get('rel_sig', pd.Series(dtype=float)).to_numpy(dtype=float, copy=False)
        self.arr['band'] = df.get('band', pd.Series(dtype=object)).astype(str).to_numpy(copy=False)
        self.arr['flags'] = df.get('sextractor_flags', pd.Series(dtype=int)).to_numpy(dtype=int, copy=False)
        self.arr['rejected'] = df.get('rejected', pd.Series(dtype=bool)).to_numpy(dtype=bool, copy=False)

    def toggle_rejection(self, index: int) -> None:
        """Toggle the rejection flag for a given index."""
        if self.df is None:
            return
        self.df.loc[index, 'rejected'] = not bool(self.df.loc[index, 'rejected'])
        self.arr['rejected'] = self.df['rejected'].to_numpy(dtype=bool, copy=False)

    def get_bands(self) -> List[str]:
        """Get the unique bands (photometric filters) in the dataset."""
        if self.df is None or 'band' not in self.df.columns:
            return []
        vals = self.df['band'].dropna().astype(str).str.strip().unique().tolist()
        desired = ['U', 'B', 'V', 'R', 'I', 'g', 'r', 'i', 'z']
        return sorted(vals, key=lambda b: desired.index(b) if b in desired else 999)


# ──────────────────── Plotting layer ────────────────────────────
class BlitManager:
    """Minimal blitting helper for fast, flicker‑free updates."""

    def __init__(self, canvas: FigureCanvasTkAgg, ax: matplotlib.axes.Axes):
        self.canvas = canvas
        self.ax = ax
        self._bg = None
        self._artists: Set[matplotlib.artist.Artist] = set()
        self.cid_draw = canvas.mpl_connect("draw_event", self._on_draw)

    def add_artists(self, artists: Iterable[matplotlib.artist.Artist]):
        """Add artists to the blitting manager."""
        for a in artists:
            if a is None:
                continue
            try:
                a.set_animated(True)
            except Exception:
                pass
            self._artists.add(a)

    def clear(self):
        self._artists.clear()
        self._bg = None

    def _on_draw(self, _evt):
        self._bg = self.canvas.copy_from_bbox(self.ax.bbox)

    def draw(self):
        self.canvas.draw_idle()

    def quick_redraw(self):
        if self._bg is None:
            # First frame after rebuild: draw now so bg is captured immediately
            self.canvas.draw()
            try:
                self._bg = self.canvas.copy_from_bbox(self.ax.bbox)
            except Exception:
                return
        self.canvas.restore_region(self._bg)
        for a in list(self._artists):
            try:
                self.ax.draw_artist(a)
            except Exception:
                pass
        self.canvas.blit(self.ax.bbox)
        self.canvas.flush_events()


class LightCurvePlot:
    """Lightcurve plotter with zooming and blitting."""
    def __init__(self, fig: matplotlib.figure.Figure, ax: matplotlib.axes.Axes, canvas: FigureCanvasTkAgg):
        self.fig = fig
        self.ax = ax
        self.canvas = canvas

        self.xlabel = "Julian Date"
        self.ylabel = "Magnitude"
        self.title = "Lightcurve"
        self.auto_title = True
        self.y_limits: Optional[Tuple[float, float]] = None
        self.y_nticks: Optional[int] = None
        
        # User-defined zoom state (persistent across updates)
        self.user_xlim: Optional[Tuple[float, float]] = None
        self.user_ylim: Optional[Tuple[float, float]] = None

        # Asteroid image zoom level
        self.asteroid_zoom_level = 1.0
        self._current_asteroid_image_original = None
        self._current_asteroid_overlay_original = None

        self.marker_size = 4.0
        self.marker_style = 'o'
        self.errorbar_capsize = 2.0
        self.errorbar_capthick = 1.0
        self.errorbar_linewidth = 1.0

        self.available_markers = PlotSettings.MARKERS
        self.marker_dict = {name: marker for marker, name in self.available_markers}

        self.valid_color = "blue"
        self.rejected_color = "red"
        self.flagged_color = "orange"

        self.parent_gui: Optional["LightCurveGUI"] = None

        # entry registry: alias -> list of band entries
        self.plotted_handles: Dict[str, List[dict]] = {}
        self._layout_signature: Optional[str] = None

        # persistent selection marker
        (self.sel_artist,) = self.ax.plot([], [], PlotSettings.SELECTION_MARKER['form'], color=PlotSettings.SELECTION_MARKER['color'],
                                          markersize=PlotSettings.SELECTION_MARKER['markersize'],
                                          zorder=PlotSettings.SELECTION_MARKER['zorder'],
                                          visible=False,
                                          animated=True)
        # persistent selection text (filename) placed near the marker
        try:
            self.sel_text = self.ax.text(0, 0, '', color=PlotSettings.SELECTION_MARKER['color'], fontsize=9,
                                         zorder=PlotSettings.SELECTION_MARKER['zorder'] + 1, visible=False, animated=True)
        except Exception:
            self.sel_text = None

        self.blit = BlitManager(self.canvas, self.ax)
        self._register_blit_artists()

    def _register_blit_artists(self):
        artists: List[matplotlib.artist.Artist] = [self.sel_artist]
        if getattr(self, 'sel_text', None) is not None:
            artists.append(self.sel_text)
        for entries in self.plotted_handles.values():
            for e in entries:
                for key in ('valid', 'flagged', 'rejected'):
                    h = e.get(key)
                    if h is None:
                        continue
                    if hasattr(h, 'lines'):  # ErrorbarContainer
                        artists.extend(getattr(h, 'lines', []))
                        artists.extend(getattr(h, 'caplines', []))
                        artists.extend(getattr(h, 'barlinecols', []))
                    elif isinstance(h, (list, tuple)):
                        artists.extend(h)
                    else:
                        artists.append(h)
        self.blit.add_artists(artists)

    def _clear_axes(self):
        self.ax.cla()
        # respect plot setting for grid visibility
        grid_on = True
        try:
            grid_on = bool(self.parent_gui.show_grid)
        except Exception:
            grid_on = True
        # Avoid supplying line/grid properties when disabling the grid; passing
        # kwargs while the first arg is False triggers a matplotlib UserWarning.
        if grid_on:
            self.ax.grid(alpha=0.15)
        else:
            # Explicitly turn grid off
            self.ax.grid(False)
        self.plotted_handles.clear()
        self.blit.clear()
        # re-create the selection artist on the fresh Axes using SELECTION_MARKER defaults
        (self.sel_artist,) = self.ax.plot([], [], PlotSettings.SELECTION_MARKER.get('form', 'o'),
                                          color=PlotSettings.SELECTION_MARKER.get('color', 'red'),
                                          markersize=PlotSettings.SELECTION_MARKER.get('markersize', 5),
                                          zorder=PlotSettings.SELECTION_MARKER.get('zorder', 20), visible=False,
                                          animated=True)
        try:
            self.sel_text = self.ax.text(0, 0, '', color=PlotSettings.SELECTION_MARKER.get('color', 'red'), fontsize=9,
                                         zorder=PlotSettings.SELECTION_MARKER.get('zorder', 20) + 1, visible=False, animated=True)
        except Exception:
            self.sel_text = None
        # Clear any cached label text object references so they will be recreated
        self.xlabel_text = None
        self.ylabel_text = None
        self.title_text = None

    def _set_labels(self):
        # Always (re)apply the label/title strings and expose their Text artists
        try:
            # set_xlabel/set_ylabel/set_title return Text artists; keep references
            self.xlabel_text = self.ax.set_xlabel(self.xlabel)
            self.ylabel_text = self.ax.set_ylabel(self.ylabel)
            self.title_text = self.ax.set_title(self.title)
            # make them pickable so GUI can respond to clicks on them
            for t in (self.xlabel_text, self.ylabel_text, self.title_text):
                try:
                    # small tolerance for picking
                    t.set_picker(5)
                except Exception:
                    pass
        except Exception:
            # Fallback to previous behaviour if anything goes wrong
            if self.ax.get_xlabel() != self.xlabel:
                self.ax.set_xlabel(self.xlabel)
            if self.ax.get_ylabel() != self.ylabel:
                self.ax.set_ylabel(self.ylabel)
            if self.ax.get_title() != self.title:
                self.ax.set_title(self.title)
        # enforce magnitude axis direction exactly once
        lo, hi = self.ax.get_ylim()
        
        # Apply manual limits if set
        if self.y_limits is not None:
            req_lo, req_hi = self.y_limits
            self.ax.set_ylim(req_lo, req_hi)
            lo, hi = self.ax.get_ylim()

        # Apply nticks if set
        if self.y_nticks is not None and self.y_nticks > 1:
            from matplotlib.ticker import MaxNLocator
            self.ax.yaxis.set_major_locator(MaxNLocator(nbins=self.y_nticks))
        else:
            # Reset to auto locator if nticks is None or invalid
            from matplotlib.ticker import AutoLocator
            self.ax.yaxis.set_major_locator(AutoLocator())

        # Auto-scale logic correction for min extent
        if self.y_limits is None:
            # Check extent
            extent = abs(hi - lo)
            if extent < 0.2:
                mid = (hi + lo) / 2.0
                # Force at least 0.2 extent
                new_half = 0.1
                # We want to preserve direction
                if lo > hi: # currently inverted
                    self.ax.set_ylim(mid + new_half, mid - new_half)
                else:
                    self.ax.set_ylim(mid - new_half, mid + new_half)
                lo, hi = self.ax.get_ylim()

        if lo < hi:
            self.ax.set_ylim(hi, lo)

    def set_y_limits(self, ymin: Optional[float], ymax: Optional[float], nticks: Optional[int] = None):
        if ymin is None or ymax is None:
            self.y_limits = None
        else:
            self.y_limits = (ymin, ymax)
        self.y_nticks = nticks

    def _choose_color_for_band(self, b: Optional[str], idx: int) -> str:
        # This method is now an instance method to access self.parent_gui.plot_settings
        if self.parent_gui and hasattr(self.parent_gui, 'plot_settings'):
             return self.parent_gui.plot_settings.get_color(b, idx)
        # Fallback if no GUI or settings attached
        return PlotSettings.DEFAULT_COLORS[idx % len(PlotSettings.DEFAULT_COLORS)]

    def _compute_time_x(self, arr_jd: np.ndarray, time_mode: str, period_hours: Optional[float]) -> np.ndarray:
        if arr_jd is None or arr_jd.size == 0:
            return np.array([])
        if time_mode == 'julian_date':
            self.xlabel = "Julian Date"
            return arr_jd
        if time_mode == 'mjd':
            self.xlabel = "Modified Julian Date (MJD)"
            return arr_jd - 2400000.5
        if time_mode == 'minutes':
            # common origin across all LCs
            jmins = []
            if self.parent_gui is not None:
                for lc in self.parent_gui.lightcurves.values():
                    if lc.arr.get('jd') is not None and lc.arr['jd'].size:
                        jmins.append(np.nanmin(lc.arr['jd']))
            j0 = np.nanmin(jmins) if jmins else np.nanmin(arr_jd)
            x = (arr_jd - j0) * 24.0 * 60.0
            try:
                base = Time(j0, format='jd').to_value('iso', subfmt='date_hm')
                self.xlabel = f"Minutes from {base} UT"
            except Exception:
                self.xlabel = "Minutes"
            return x
        if time_mode == 'rotation_phase':
            self.xlabel = "Rotation Phase"
            if not period_hours or period_hours <= 0:
                return arr_jd
            jmins = []
            if self.parent_gui is not None:
                for lc in self.parent_gui.lightcurves.values():
                    if lc.arr.get('jd') is not None and lc.arr['jd'].size:
                        jmins.append(np.nanmin(lc.arr['jd']))
            jd0 = np.nanmin(jmins) if jmins else np.nanmin(arr_jd)
            phases = ((arr_jd - jd0) * 24.0 / period_hours)
            return np.mod(phases, 1.5)
        return arr_jd

    def _layout_sig(self, lc_dict, visible, mode, time_mode, errorbar_type, selected_band, show_rejected, offsets):
        """Signature that uniquely identifies the layout state (forces rebuild when it changes)."""
        visible_aliases = ','.join(sorted([a for a, v in visible.items() if v and a in lc_dict]))
        # debug_print(f"Offsets: {offsets}")
        offsets_data = ','.join(offset for offset in
                                (f"{alias}:" + ','.join(f"{band}={offsets[alias][band]:.4f}"
                                                        for band in sorted(offsets[alias].keys()))
                                 for alias in sorted(offsets.keys()) if alias in lc_dict))
        debug_print(f"Offset data: {offsets_data}")
        rp = ""
        try:
            if time_mode == 'rotation_phase' and self.parent_gui is not None and hasattr(self.parent_gui,
                                                                                         'rotation_period_var'):
                rp_val = float(self.parent_gui.rotation_period_var.get())
                rp = f"|P={rp_val:.6f}"
        except Exception:
            rp = "|P=?"
        # include show_rejected so toggling it can rebuild when needed
        # also include legend/grid toggles so toggling them forces a rebuild
        try:
            legend_flag = int(bool(getattr(self.parent_gui, 'show_legend', True)))
        except Exception:
            legend_flag = 1
        try:
            grid_flag = int(bool(getattr(self.parent_gui, 'show_grid', True)))
        except Exception:
            grid_flag = 1
        try:
            colors_flag = int(bool(getattr(self.parent_gui.show_colors_on_plot_var, 'get', lambda: False)()))
        except Exception:
            colors_flag = 0
        try:
            colors_loc = str(getattr(self.parent_gui.color_legend_loc_var, 'get', lambda: "upper left")())
        except Exception:
            colors_loc = "upper left"

        signature = (f"{visible_aliases}|{mode}|{time_mode}{rp}|{errorbar_type}"
                     f"|{selected_band}|rej={int(bool(show_rejected))}|legend={legend_flag}|grid={grid_flag}|colors={colors_flag}|loc={colors_loc}"
                     f"|{self.marker_style}|{self.marker_size}|{self.valid_color}|{offsets_data}")
        return signature

    def invalidate_layout(self) -> None:
        """Force a full rebuild on next update (e.g. when masks change)."""
        self._layout_signature = None

    # — main update ———
    def update(self,
               lc_dict: Dict[str, LightCurveData],
               offsets: Dict[str, Dict[str, float]],
               visible: Dict[str, bool],
               mode: str,
               time_mode: str,
               show_rejected: bool,
               errorbar_type: str,
               selected_bands: Set[str]):
        if not lc_dict:
            self._clear_axes()
            self._set_labels()
            self.canvas.draw()
            return

        new_sig = self._layout_sig(lc_dict, visible, mode, time_mode, errorbar_type, selected_bands, show_rejected,
                                   offsets)
        rebuild = (new_sig != self._layout_signature)

        if rebuild:
            debug_print(f"Requested full rebuild: {new_sig} != {self._layout_signature}")
            self._clear_axes()
            plotted_any = False
            for alias, lc in lc_dict.items():
                if not visible.get(alias, True) or lc.df is None or lc.df.empty:
                    continue
                arr = lc.arr
                jd = arr.get('jd')
                if jd is None or not jd.size:
                    continue

                # choose Y and Yerr arrays by mode
                # choose Y array by mode (which magnitude column to plot)
                if mode == 'target':
                    Yfull = arr.get('mag')
                    # check if not empty (give the message)
                    if Yfull is None or not Yfull.size:
                        messagebox.showwarning("No valid target magnitudes",
                                               f"No valid 'mag' values to plot.")
                        continue
                elif mode == 'instrumental':
                    Yfull = arr.get('inst_mag')
                    if Yfull is None or not Yfull.size:
                        messagebox.showwarning("No valid target magnitudes",
                                               f"No valid instrumental (inst_mag) values to plot.")
                        continue
                elif mode == 'relative':
                    Yfull = arr.get('rel_mag')
                    if Yfull is None or not Yfull.size:
                        messagebox.showwarning("No valid relative magnitudes",
                                               f"No valid relative (rel_mag) values to plot.")
                else:  # control
                    Yfull = arr.get('mag_control')
                    if Yfull is None or not Yfull.size:
                        messagebox.showwarning("No valid magnitudes for control star",
                                               f"No valid data (mag_control) for control star to plot.")

                # choose Yerr array based on the selected errorbar type regardless of mode
                # (previously Yerr was tied to the current mode which caused missing/incorrect
                # errorbars when users selected a different error source than the plotted mag)
                if errorbar_type == 'calibrated':
                    Yerr_full = arr.get('sig')
                elif errorbar_type == 'instrumental':
                    Yerr_full = arr.get('inst_sig')
                elif errorbar_type == 'relative':
                    Yerr_full = arr.get('rel_sig')
                else:  # 'none' or unknown
                    Yerr_full = None

                flags = arr.get('flags')
                rejected = arr.get('rejected')
                band_col = arr.get('band')
                if flags is None:
                    flags = np.zeros_like(jd, dtype=int)
                if rejected is None:
                    rejected = np.zeros_like(jd, dtype=bool)

                lc_bands = lc.get_bands() or [None]
                if selected_bands is None:
                    # None => show all available bands
                    bands_to_draw = lc_bands
                else:
                    # set (possibly empty) => show only those selected
                    bands_to_draw = [b for b in lc_bands
                                     if (b in selected_bands) or (b is None and 'None' in selected_bands)]

                # common X transform per band
                period = None
                if time_mode == 'rotation_phase' and self.parent_gui is not None:
                    try:
                        period = float(self.parent_gui.rotation_period_var.get())
                    except Exception:
                        period = None

                for bidx, b in enumerate(bands_to_draw):
                    band_mask = np.ones_like(jd, dtype=bool)
                    if b is not None and band_col is not None and len(band_col) == len(jd):
                        band_mask = (band_col.astype(str) == str(b))
                    if not np.any(band_mask):
                        continue

                    X = self._compute_time_x(jd[band_mask], time_mode, period)
                    Ybase = Yfull[band_mask]
                    Yerr = (Yerr_full[band_mask] if (Yerr_full is not None and len(Yerr_full) == len(jd)) else None)
                    flg = flags[band_mask]
                    rej = rejected[band_mask]

                    mask_valid = (~rej) & (flg == 0)
                    mask_flagged = (flg > 0) & (~rej)
                    mask_rejected = rej

                    # per‑band offset
                    off = offsets.get(alias, {}).get(str(b) if b is not None else OFFSET_ALL_KEY,
                                                     offsets.get(alias, {}).get(OFFSET_ALL_KEY, 0.0))

                    color = self._choose_color_for_band(b, bidx)

                    # determine flagged/rejected colors for this band (may follow filter color)
                    flagged_color_for_this = self.flagged_color
                    rejected_color_for_this = self.rejected_color

                    if self.parent_gui and hasattr(self.parent_gui, 'plot_settings'):
                        settings = self.parent_gui.plot_settings
                        if settings.flagged_use_filter:
                            flagged_color_for_this = color
                        else:
                            flagged_color_for_this = settings.get_color('flagged', 0)
                        rejected_color_for_this = settings.get_color('rejected', 0)

                    def make_series(m, c):
                        if not np.any(m):
                            return None
                        # Apply the offset to Y values
                        Y_plot = Ybase[m] + off

                        # show errorbars only if we actually have them and user didn't set 'none'
                        want_err = (Yerr is not None and
                                    getattr(self.parent_gui, 'errorbar_type', 'none') != 'none')

                        if want_err:
                            # Create error bars
                            container = self.ax.errorbar(X[m], Y_plot, yerr=Yerr[m], fmt='none', ecolor=c,
                                                         capsize=self.errorbar_capsize, capthick=self.errorbar_capthick,
                                                         elinewidth=self.errorbar_linewidth, zorder=5, alpha=1.0)
                            # Create separate marker line and ATTACH it to container for fast-path updates
                            line, = self.ax.plot(X[m], Y_plot, linestyle='None', marker=self.marker_style,
                                                 markersize=self.marker_size, color=c, zorder=10, picker=5)
                            container.marker_line = line  # <-- critical fix
                            return container
                        else:
                            # Marker-only series
                            line, = self.ax.plot(X[m], Y_plot, linestyle='None', marker=self.marker_style,
                                                 markersize=self.marker_size, color=c, zorder=10, picker=5)
                            return line

                    h_valid = make_series(mask_valid, color)
                    h_flagged = make_series(mask_flagged, flagged_color_for_this)
                    h_rejected = make_series(mask_rejected, rejected_color_for_this) if show_rejected else None

                    entry = {
                        'alias': alias, 'band': b, 'band_index': bidx,
                        'X': X, 'Ybase': Ybase, 'Yerr': Yerr,
                        'mask_valid': mask_valid, 'mask_flagged': mask_flagged, 'mask_rejected': mask_rejected,
                        # persist colors so fast-path updates can recolor without full rebuild
                        'color': color,
                        'flagged_color': flagged_color_for_this,
                        'rejected_color': rejected_color_for_this,
                    }
                    if h_valid is not None:
                        entry['valid'] = h_valid
                    if h_flagged is not None:
                        entry['flagged'] = h_flagged
                    if h_rejected is not None:
                        entry['rejected'] = h_rejected

                    self.plotted_handles.setdefault(alias, []).append(entry)

                    legend_label = alias if b is None else f"{alias} ({b})"
                    self.ax.plot([], [], marker=self.marker_style, linestyle='None', color=color, label=legend_label)
                    plotted_any = True

            if plotted_any:
                # only show legend when requested by the GUI
                try:
                    if getattr(self.parent_gui, 'show_legend', True):
                        self.ax.legend(title='Lightcurves', fontsize='small')
                    else:
                        # remove any existing legend
                        leg = self.ax.get_legend()
                        if leg is not None:
                            leg.remove()
                except Exception:
                    # fallback to showing legend if anything goes wrong
                    try:
                        self.ax.legend(title='Lightcurves', fontsize='small')
                    except Exception:
                        pass

            # ─── Draw Color Info on Plot ───
            if self.parent_gui and getattr(self.parent_gui.show_colors_on_plot_var, 'get', lambda: False)():
                info_lines = []
                # Only show for visible aliases
                vis_aliases = [a for a in lc_dict.keys() if visible.get(a, True)]
                
                # Check calculated colors
                if hasattr(self.parent_gui, 'calculated_colors'):
                    for alias in vis_aliases:
                        colors = self.parent_gui.calculated_colors.get(alias, {})
                        if colors:
                            for name, (val, err) in colors.items():
                                info_lines.append(f"  {name}: {val:.3f} ± {err:.3f}")
                
                if info_lines:
                    import matplotlib.offsetbox as offsetbox
                    loc_val = "upper left"
                    if self.parent_gui:
                        loc_val = self.parent_gui.color_legend_loc_var.get()
                    
                    text_str = "\n".join(info_lines)
                    at = offsetbox.AnchoredText(text_str, loc=loc_val, frameon=True, prop=dict(fontsize=9))
                    at.patch.set_boxstyle("round,pad=0.,rounding_size=0.2")
                    at.patch.set_alpha(0.7)
                    self.ax.add_artist(at)

            if self.auto_title:
                vis_aliases = [a for a in lc_dict.keys() if visible.get(a, True)]
                if len(vis_aliases) == 1:
                    only = vis_aliases[0]
                    lc = lc_dict[only]
                    if lc.df is not None and 'target' in lc.df.columns and not lc.df.empty:
                        try:
                            self.title = str(lc.df['target'].iloc[0])
                        except Exception:
                            pass

            self._set_labels()
            # Re-register blit artists for the new layout
            self._register_blit_artists()
            try:
                self.fig.tight_layout()
            except Exception:
                pass

            # ── Refresh blit state after rebuild ─────────────────────────────
            if hasattr(self, "blit") and self.blit is not None:
                # ensure a fresh background is captured for the new axes/limits
                self.blit._bg = None
                # do a guaranteed draw so draw_event can cache the background
                self.canvas.draw()
                # (optional) first fast repaint, so the user sees it immediately
                self.blit.quick_redraw()

            self._layout_signature = new_sig
            return

        debug_print("Only fast update needed")
        # fast path (update existing artists only)
        for alias, entries in self.plotted_handles.items():
            vis_alias = visible.get(alias, True)
            for e in entries:
                b = e['band']
                off = offsets.get(alias, {}).get(str(b) if b is not None else OFFSET_ALL_KEY,
                                                 offsets.get(alias, {}).get(OFFSET_ALL_KEY, 0.0))
                X = e['X']
                Y = e['Ybase'] + off

                def _set_visible(handle, visible):
                    def _vis(obj, vis):
                        if obj is None:
                            return
                        if hasattr(obj, 'set_visible'):
                            obj.set_visible(vis)
                            return
                        if isinstance(obj, (list, tuple)):
                            for o in obj:
                                _vis(o, vis)

                    if handle is None:
                        return
                    if hasattr(handle, 'lines') or hasattr(handle, 'caplines') or hasattr(handle, 'barlinecols'):
                        _vis(getattr(handle, 'lines', []), visible)
                        err_on = bool(visible and (self.parent_gui.errorbar_type != 'none'))
                        _vis(getattr(handle, 'caplines', []), err_on)
                        _vis(getattr(handle, 'barlinecols', []), err_on)
                        if hasattr(handle, 'marker_line'):
                            _vis(handle.marker_line, visible)
                        return
                    if isinstance(handle, (list, tuple)):
                        for o in handle:
                            _vis(o, visible)
                        return
                    _vis(handle, visible)

                def set_data(h, mask, color_override=None):
                    """Safely update X/Y for Line2D or ErrorbarContainer (masked)."""
                    if h is None:
                        return
                    Xm = X[mask]
                    Ym = Y[mask]

                    # ErrorbarContainer
                    if hasattr(h, 'lines'):
                        main = h.lines[0] if h.lines else None
                        if main is not None:
                            main.set_data(Xm, Ym)
                            main.set_marker(self.marker_style)
                            main.set_markersize(self.marker_size)
                            main.set_linestyle('None')
                            # set color for errorbar main/marker
                            try:
                                if color_override is not None:
                                    main.set_color(color_override)
                            except Exception:
                                pass

                        # rebuild vertical error segments only if we actually show errorbars now
                        show_err = (e['Yerr'] is not None) and (errorbar_type != 'none')
                        if show_err:
                            segs = []
                            idxs = np.nonzero(mask)[0]
                            for i, k in enumerate(idxs):
                                dy = e['Yerr'][k]
                                if not np.isfinite(dy):
                                    segs.append(np.array([[Xm[i], Ym[i]], [Xm[i], Ym[i]]]))
                                else:
                                    segs.append(np.array([[Xm[i], Ym[i] - dy], [Xm[i], Ym[i] + dy]]))
                            for lc in getattr(h, 'barlinecols', []):
                                lc.set_segments(segs)
                            for cap in getattr(h, 'caplines', []):
                                cap.set_visible(True)
                        else:
                            for lc in getattr(h, 'barlinecols', []):
                                lc.set_visible(False)
                            for cap in getattr(h, 'caplines', []):
                                cap.set_visible(False)

                        if hasattr(h, 'marker_line') and h.marker_line is not None:
                            ml = h.marker_line
                            ml.set_data(Xm, Ym)
                            ml.set_marker(self.marker_style)
                            ml.set_markersize(self.marker_size)
                            ml.set_linestyle('None')
                            try:
                                if color_override is not None:
                                    ml.set_color(color_override)
                            except Exception:
                                pass
                        return

                    # list/tuple of artists
                    if isinstance(h, (list, tuple)):
                        for hh in h:
                            set_data(hh, mask)
                        return

                    # Plain Line2D
                    h.set_data(Xm, Ym)
                    h.set_marker(self.marker_style)
                    h.set_markersize(self.marker_size)
                    h.set_linestyle('None')
                    try:
                        if color_override is not None:
                            h.set_color(color_override)
                    except Exception:
                        pass

                mask_v = e['mask_valid']
                mask_f = e['mask_flagged']
                mask_r = e['mask_rejected']

                # push updated coordinates into existing artists so offsets show immediately
                h = e.get('valid')
                if h is not None:
                    set_data(h, mask_v, e.get('color'))

                h = e.get('flagged')
                if h is not None:
                    set_data(h, mask_f, e.get('flagged_color'))

                h = e.get('rejected')
                if h is not None:
                    set_data(h, mask_r, e.get('rejected_color'))

                # visibility per entry
                want_r = vis_alias and show_rejected and np.any(e['mask_rejected'])
                want_v = vis_alias and np.any(e['mask_valid'])
                want_f = vis_alias and np.any(e['mask_flagged'])

                h = e.get('valid')
                if h is not None:
                    _set_visible(h, want_v)

                h = e.get('flagged')
                if h is not None:
                    _set_visible(h, want_f)

                h = e.get('rejected')
                if h is not None:
                    _set_visible(h, want_r)

        # selection marker: compute with same X transform + offsets
        sel = getattr(self, 'selection', {'alias': None, 'index': None})
        sel_alias = sel.get('alias')
        sel_idx = sel.get('index')
        if sel_alias in lc_dict and sel_idx is not None:
            lc = lc_dict[sel_alias]
            df = lc.df
            arr = lc.arr
            if df is not None and 0 <= sel_idx < len(df):
                jd_all = arr.get('jd')
                # band of selected point (for per‑band offset)
                band_val = None
                if 'band' in df.columns:
                    try:
                        band_val = str(df['band'].iloc[sel_idx])
                    except Exception:
                        band_val = None
                period = None
                if time_mode == 'rotation_phase' and self.parent_gui is not None:
                    try:
                        period = float(self.parent_gui.rotation_period_var.get())
                    except Exception:
                        period = None
                X_all = self._compute_time_x(jd_all, time_mode, period)
                sx = X_all[sel_idx]
                # choose Y for mode
                if mode == 'target' and 'mag' in df.columns:
                    sy = float(df['mag'].iloc[sel_idx])
                elif mode == 'instrumental' and 'inst_mag' in df.columns:
                    sy = float(df['inst_mag'].iloc[sel_idx])
                else:
                    sy = float(df.get('mag_control', pd.Series(np.zeros(len(df))))[sel_idx])
                off = offsets.get(sel_alias, {}).get(band_val, offsets.get(sel_alias, {}).get(OFFSET_ALL_KEY, 0.0))
                self.sel_artist.set_data([sx], [sy + off])
                self.sel_artist.set_visible(True)
                # show filename near the selection marker if available
                txt = ''
                try:
                    if 'filename' in df.columns:
                        txt = str(df['filename'].iloc[sel_idx])
                except Exception:
                    txt = ''
                if getattr(self, 'sel_text', None) is not None and txt:
                    # offset text slightly above marker in data coordinates
                    try:
                        lo, hi = self.ax.get_ylim()
                        # span may be inverted; use absolute span
                        yspan = abs(hi - lo) if hi is not None and lo is not None else 0.0
                        y_offset = 0.05 * (yspan if yspan > 0 else 1.0)

                        # Check if text would go outside the right edge
                        xlim = self.ax.get_xlim()
                        x_range = xlim[1] - xlim[0]
                        # Estimate text width as ~0.5 of plot width for typical filename lengths
                        text_width_estimate = 0.25 * x_range

                        # If text extends beyond right edge, place it to the left of marker
                        if sx + text_width_estimate > xlim[1]:
                            self.sel_text.set_position((sx, sy + off + y_offset))
                            self.sel_text.set_horizontalalignment('right')
                        else:
                            self.sel_text.set_position((sx, sy + off + y_offset))
                            self.sel_text.set_horizontalalignment('left')

                        # shorten long filenames
                        display = txt if len(txt) <= 40 else txt[:36] + '...'
                        self.sel_text.set_text(display)
                        self.sel_text.set_visible(True)
                    except Exception:
                        try:
                            self.sel_text.set_visible(False)
                        except Exception:
                            pass
                else:
                    self.sel_text.set_visible(False)
            else:
                self.sel_artist.set_visible(False)
                if getattr(self, 'sel_text', None) is not None:
                    try:
                        self.sel_text.set_visible(False)
                    except Exception:
                        pass
        else:
            self.sel_artist.set_visible(False)
            if getattr(self, 'sel_text', None) is not None:
                try:
                    self.sel_text.set_visible(False)
                except Exception:
                    pass

        self._set_labels()
        self.blit.quick_redraw()

    # ——— point picking ———
    def select_nearest_point(self, xdata: float, ydata: float,
                             lc_dict: Dict[str, LightCurveData],
                             offsets: Dict[str, Dict[str, float]],
                             mode: str, time_mode: str,
                             selected_bands: Set[str], show_rejected: bool) -> Tuple[
        Optional[str], Optional[int], Optional[str]]:
        pts = []  # (x_disp, y_disp, alias, df_idx, band)
        trans = self.ax.transData
        for alias, lc in lc_dict.items():
            if lc.df is None or lc.df.empty:
                continue
            df = lc.df;
            arr = lc.arr
            jd = arr.get('jd')
            if jd is None or not jd.size:
                continue
            bands = arr.get('band')
            flags = arr.get('flags') if arr.get('flags') is not None else np.zeros(len(df), dtype=int)
            rejected = arr.get('rejected') if arr.get('rejected') is not None else np.zeros(len(df), dtype=bool)

            band_mask = np.ones(len(df), dtype=bool)
            if (selected_bands is not None) and ('band' in df.columns):
                band_mask = df['band'].astype(str).isin(selected_bands).to_numpy()

            vis_mask = band_mask & (~np.isnan(jd))
            if not show_rejected:
                vis_mask &= (~rejected)

            period = None
            if time_mode == 'rotation_phase' and self.parent_gui is not None:
                try:
                    period = float(self.parent_gui.rotation_period_var.get())
                except Exception:
                    period = None
            X_all = self._compute_time_x(jd, time_mode, period)

            if mode == 'target':
                Y_all = arr.get('mag')
            elif mode == 'instrumental':
                Y_all = arr.get('inst_mag')
            else:
                Y_all = arr.get('mag_control')
            if Y_all is None:
                continue

            off_all = np.zeros(len(df))
            if bands is not None:
                # vectorized per‑row offset lookup
                for i in np.nonzero(vis_mask)[0]:
                    bval = str(bands[i]) if bands is not None else None
                    off_all[i] = offsets.get(alias, {}).get(bval, offsets.get(alias, {}).get(OFFSET_ALL_KEY, 0.0))
            else:
                off_all[:] = offsets.get(alias, {}).get(OFFSET_ALL_KEY, 0.0)

            Xv = X_all[vis_mask]
            Yv = (Y_all + off_all)[vis_mask]
            bands_v = (bands[vis_mask] if bands is not None else np.array([None] * len(Xv)))

            xy = trans.transform(np.column_stack((Xv, Yv)))
            vis_idx = np.nonzero(vis_mask)[0]
            for k, (xd, yd) in enumerate(xy):
                pts.append((xd, yd, alias, int(vis_idx[k]), str(bands_v[k]) if bands is not None else None))

        if not pts:
            return None, None, None
        click_disp = self.ax.transData.transform((xdata, ydata))
        d2 = [(xd - click_disp[0]) ** 2 + (yd - click_disp[1]) ** 2 for xd, yd, *_ in pts]
        i = int(np.argmin(d2))
        if d2[i] < (15.0 ** 2):
            _, _, alias, idx, b = pts[i]
            return alias, idx, b
        return None, None, None


# ───────────────────────────── Asteroid Viewer ────────────────────────────
class AsteroidImageViewer:
    """
    Manages the asteroid image window, including loading images,
    zooming, and persisting window state/zoom level.
    """
    def __init__(self, root: tk.Tk, data_provider):
        """
        :param root: The main application root window.
        :param data_provider: An object (LightCurveGUI) that provides access to
                              current_lc_alias, current_point_index, lightcurves, etc.
        """
        self.root = root
        self.data_provider = data_provider

        # State persistence
        self.zoom_level = 1.0
        self.window_geometry = None
        self.image_window = None

        # Image caching
        self._current_image_original = None
        self._current_overlay_original = None

        # Overlay toggle
        self.show_overlay_var = tk.BooleanVar(value=True)

    def toggle_visibility(self):
        """Toggles the visibility of the asteroid image window."""
        if self.image_window is not None and self.image_window.winfo_exists():
            try:
                if str(self.image_window.state()) == "withdrawn":
                    self._restore_window()
                else:
                    self._hide_window()
            except Exception:
                pass
        else:
            self._create_window()

    def _restore_window(self):
        if self.window_geometry:
            try:
                self.image_window.geometry(self.window_geometry)
            except Exception:
                pass
        self.image_window.deiconify()
        self.image_window.lift()
        # Keep focus on main window so keys work there (delayed to override window manager)
        self.root.after(100, lambda: self.root.focus_force())

    def _hide_window(self):
        try:
            self.window_geometry = self.image_window.geometry()
        except Exception:
            pass
        self.image_window.withdraw()

    def _create_window(self):
        try:
            top = tk.Toplevel(self.root)
            self.image_window = top
            top.title("Asteroid image")

            # Keep on top of main window
            top.transient(self.root)

            if self.window_geometry:
                try:
                    top.geometry(self.window_geometry)
                except Exception:
                    pass

            self._bind_events(top)
            self.update_image()

            # Overlay toggle checkbox
            chk = ttk.Checkbutton(top, text="Toggle Overlay (O)", variable=self.show_overlay_var,
                                  command=self.refresh_display)
            chk.pack(side=BOTTOM, pady=5)

            # 1. Immediate request
            self.root.focus_set()
            
            # 2. When the new window is actually mapped/shown, force focus back
            def return_focus(event=None):
                self.root.focus_force()
                # Unbind to avoid repeated calls
                top.unbind("<Map>")
            
            top.bind("<Map>", return_focus)
            
            # 3. Fail-safe delayed triggers
            self.root.after(50, lambda: self.root.focus_force())
            self.root.after(200, lambda: self.root.focus_force())

        except Exception:
            pass

    def _bind_events(self, win: tk.Toplevel):
        # Track geometry
        def _on_configure(evt=None):
            try:
                if evt is None or evt.widget is win:
                    self.window_geometry = win.geometry()
            except Exception:
                pass
        win.bind("<Configure>", _on_configure, add="+")

        # Zoom
        win.bind("<MouseWheel>", self._on_zoom)
        win.bind("<Button-4>", self._on_zoom)
        win.bind("<Button-5>", self._on_zoom)

        # Hotkey for overlay
        win.bind("<o>", self._toggle_overlay_hotkey)
        win.bind("<O>", self._toggle_overlay_hotkey)

        # Cleanup on destroy (if destroyed externally)
        # Note: We don't strictly need to bind destroy if we check winfo_exists,
        # but it's good practice to clear the ref.
        # However, binding <Destroy> can be tricky if it triggers on child widgets.
        # We'll rely on winfo_exists checks.

    def update_image(self):
        """Loads and displays the image for the current selection."""
        if not self.image_window or not self.image_window.winfo_exists():
            return

        # Clear previous content
        for widget in self.image_window.winfo_children():
            widget.destroy()

        alias = self.data_provider.current_lc_alias
        idx = self.data_provider.current_point_index

        if alias is None or idx is None:
            ttk.Label(self.image_window, text="No point selected.").pack(padx=20, pady=20)
            self.image_window.title("Asteroid Image")
            return

        try:
            lc = self.data_provider.lightcurves[alias]
            row_data = lc.df.iloc[idx]
        except (KeyError, IndexError):
            ttk.Label(self.image_window, text="Error retrieving data.").pack(padx=20, pady=20)
            return

        # Locate image
        cwd = os.getcwd()
        image_path = None
        overlay_path = None

        if 'filename' in row_data and pd.notna(row_data['filename']):
            try:
                if self.data_provider.mode == 'control':
                    image_name = Path(row_data['filename']).stem
                    image_path = [str(path) for path in Path(cwd).glob(f'**/Control_Star*{image_name}_thumb.png')][0]
                    overlay_path = [str(path) for path in Path(cwd).glob(f'**/Control_Star*{image_name}_thumb_overlay.png')][0]
                else:
                    image_name = Path(row_data['filename']).stem
                    image_path = [str(path) for path in Path(cwd).glob(f'**/*__{image_name}_thumb.png')][0]
                    overlay_path = [str(path) for path in Path(cwd).glob(f'**/*__{image_name}_thumb_overlay.png')][0]
            except IndexError:
                pass

        if not image_path or pd.isna(image_path):
            ttk.Label(self.image_window, text="No image path found.").pack(padx=20, pady=20)
            return

        # Resolve absolute path
        if not os.path.isabs(image_path) and lc.filename:
            csv_dir = os.path.dirname(lc.filename)
            full_path = os.path.join(csv_dir, image_path)
            full_path_overlay = os.path.join(csv_dir, overlay_path) if overlay_path else None
        else:
            full_path = image_path
            full_path_overlay = overlay_path

        if not os.path.exists(full_path):
            ttk.Label(self.image_window, text=f"File not found:\n{os.path.basename(full_path)}").pack(padx=20, pady=20)
            return

        # Load images
        try:
            img = Image.open(full_path)
            overlay = None
            if full_path_overlay and os.path.exists(full_path_overlay):
                try:
                    overlay = Image.open(full_path_overlay)
                except Exception:
                    pass

            self._current_image_original = img
            self._current_overlay_original = overlay

            # IMPORTANT: Do NOT reset zoom_level here to persist it across images.
            # self.zoom_level = 1.0

            self.refresh_display()
            self.image_window.title(f"Image: {os.path.basename(full_path)}")

        except Exception as e:
            ttk.Label(self.image_window, text=f"Error loading image:\n{e}").pack(padx=20, pady=20)

    def refresh_display(self):
        """Resizes and displays the stored image based on current zoom."""
        if not self.image_window or not self.image_window.winfo_exists():
            return
        if self._current_image_original is None:
            return

        # Find or create label
        img_label = None
        for widget in self.image_window.winfo_children():
            if isinstance(widget, ttk.Label) and hasattr(widget, 'image'):
                img_label = widget
                break

        if img_label is None:
            # Clear anything else (like error messages)
            for widget in self.image_window.winfo_children():
                widget.destroy()
            img_label = ttk.Label(self.image_window)
            img_label.pack(padx=10, pady=10)

        # Resize
        orig_w, orig_h = self._current_image_original.size
        new_w = int(orig_w * self.zoom_level)
        new_h = int(orig_h * self.zoom_level)

        # Use LANCZOS for quality
        resized_img = self._current_image_original.resize((new_w, new_h), Image.LANCZOS)

        if self._current_overlay_original and self.show_overlay_var.get():
            resized_overlay = self._current_overlay_original.resize((new_w, new_h), Image.LANCZOS)
            resized_img.paste(resized_overlay, (0, 0), resized_overlay)

        photo = ImageTk.PhotoImage(resized_img)
        img_label.configure(image=photo)
        img_label.image = photo

    def _on_zoom(self, event):
        if self._current_image_original is None:
            return

        scale_factor = 1.1
        # Windows/MacOS: event.delta, Linux: event.num
        if event.num == 4 or event.delta > 0:
            self.zoom_level *= scale_factor
        elif event.num == 5 or event.delta < 0:
            self.zoom_level /= scale_factor

        # Clamp
        if self.zoom_level < 0.1:
            self.zoom_level = 0.1
        elif self.zoom_level > 20.0:
            self.zoom_level = 20.0

        self.refresh_display()

    def _toggle_overlay_hotkey(self, event=None):
        self.show_overlay_var.set(not self.show_overlay_var.get())
        self.refresh_display()


# ───────────────────────────────── GUI layer ──────────────────────────────
class LightCurveGUI:
    def __init__(self, root: ttk.Window):
        self.root = root
        self.root.title("Lightcurve Viewer and Editor")

        # Data model: alias -> LightCurveData
        self.lightcurves: Dict[str, LightCurveData] = {}
        # Per‑file visibility / color / offsets
        self.lc_visible: Dict[str, bool] = {}
        # per‑file offsets: alias -> { band -> offset }
        self.lc_offsets: Dict[str, Dict[str, float]] = {}

        self._left_collapsed = False
        self._left_saved_width = 260  # last expanded width
        self.LEFT_DEFAULT_WIDTH = 260  # default width at startup
        self.LEFT_MINSIZE = 120  # minimum when expanded

        # selection
        self.current_lc_alias: Optional[str] = None
        self.current_point_index: Optional[int] = None
        self.current_point_band: Optional[str] = None

        # Asteroid Image Viewer
        self.asteroid_viewer = AsteroidImageViewer(self.root, self)

        # View state
        self.mode = 'target'  # 'target' | 'instrumental' | 'control' | 'relative'
        self.time_mode = 'minutes'  # 'minutes' | 'julian_date' | 'mjd' | 'rotation_phase'
        self.show_rejected = True
        self.errorbar_type = 'calibrated'  # 'instrumental' | 'calibrated' | 'none'

        # Multi‑filter state
        self.selected_bands: Set[str] = set()  # empty => all bands
        self.band_vars: Dict[str, tk.BooleanVar] = {}

        # Offset control scope (which band to nudge when ↑/↓)
        self.offset_scope_var = ttk.StringVar(value='auto')  # 'auto', 'ALL', or specific band label

        # Legend position
        self.color_legend_loc_var = ttk.StringVar(value='upper left')

        # Alignment controls
        self.align_reference_var = ttk.StringVar(value='')
        self.calculated_colors: Dict[str, Dict[str, Tuple[float, float]]] = {}  # alias -> { 'R-V': (val, err) }

        # Step used when nudging vertical offsets with ↑/↓ (mag)
        self.offset_step_var = ttk.DoubleVar(value=0.05)

        # Plot Settings
        self.plot_settings = PlotSettings(self.root, self._refresh_plot_callback, self._get_bands_callback)

        # Plot display toggles (exposed via Plot Settings)
        self.show_legend: bool = True
        self.show_grid: bool = True

        # Build UI
        self._build_ui()

        # connect plot backref
        self.plot.parent_gui = self

        # Bindings
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)
        self.root.bind("q", lambda e: self.confirm_exit())
        self.root.bind("r", lambda e: self.toggle_rejection())
        self.root.bind("a", lambda e: self.cancel_selection())
        self.root.bind("<Left>", lambda e: self.move_point(-1))
        self.root.bind("<Right>", lambda e: self.move_point(+1))
        self.root.bind("s", lambda e: self.adjust_lc_offset(-self._get_offset_step()))
        self.root.bind("d", lambda e: self.adjust_lc_offset(+self._get_offset_step()))
        self.root.bind("z", lambda e: self.adjust_rotation_period(-1))
        self.root.bind("x", lambda e: self.adjust_rotation_period(+1))
        # Bind Shift+S to show the image for the selected point
        self.root.bind("<S>", self.show_asteroid_image)

    def _refresh_plot_callback(self):
        if hasattr(self, 'plot') and getattr(self.plot, 'invalidate_layout', None):
            self.plot.invalidate_layout()
        self.request_plot_update()

    def _get_bands_callback(self) -> List[str]:
        return sorted({str(b) for lc in self.lightcurves.values() for b in lc.get_bands()})

    # ——— UI construction ———
    # Replace the top layout inside _build_ui()
    def _build_ui(self):
        master = ttk.Frame(self.root)
        master.pack(fill=BOTH, expand=True)
        self.master_frame = master

        # Create a persistent bottom control strip
        self.bottom_control = ttk.Frame(self.master_frame)
        self.bottom_control.pack(side=BOTTOM, fill=X)
        try:
            self.bottom_control.pack_propagate(False)
            self.bottom_control.configure(height=72)
        except Exception:
            pass

        # Main content area
        self.center_pane = ttk.Frame(master)
        self.center_pane.pack(side=TOP, fill=BOTH, expand=True)

        # Top controls
        top_frame = ttk.LabelFrame(self.center_pane, text='Controls', padding=0)
        top_frame.pack(side=TOP, fill=X, padx=6, pady=(6, 0))

        # FlowFrame: converts packed children into a wrapped grid layout on resize
        class FlowFrame(ttk.Frame):
            def __init__(self, master=None, spacing_x=6, spacing_y=4, **kwargs):
                super().__init__(master, **kwargs)
                self.spacing_x = spacing_x
                self.spacing_y = spacing_y
                self._in_layout = False
                self.bind('<Configure>', self._on_configure)

            def _on_configure(self, event=None):
                if self._in_layout:
                    return
                self._in_layout = True
                try:
                    w = (event.width if event is not None else self.winfo_width())
                    self._do_layout(w)
                finally:
                    self._in_layout = False

            def _do_layout(self, width):
                x = 0
                row = 0
                col = 0
                for child in self.winfo_children():
                    try:
                        child.update_idletasks()
                    except Exception:
                        pass
                    reqw = child.winfo_reqwidth()
                    try:
                        if child.winfo_manager() == 'pack':
                            child.pack_forget()
                    except Exception:
                        pass
                    if col > 0 and (x + reqw > max(10, width)):
                        row += 1
                        col = 0
                        x = 0
                    try:
                        child.grid(row=row, column=col, sticky='w', padx=(0, self.spacing_x), pady=(0, self.spacing_y))
                    except Exception:
                        pass
                    x += reqw + self.spacing_x
                    col += 1

        top = FlowFrame(top_frame, padding=6)
        top.pack(fill=X, expand=True)

        ttk.Button(top, text="Open CSV", command=self.open_csv).pack(side=LEFT, padx=4)
        self._build_save_menu(top)

        ttk.Label(top, text="Mode:").pack(side=LEFT, padx=(10, 0))
        self.mode_var = ttk.StringVar(value=self.mode)
        photo_mode_select = ttk.Combobox(top, textvariable=self.mode_var,
                                         values=['target', 'instrumental', 'relative', 'control'],
                                         state='readonly', width=14)
        photo_mode_select.pack(side=LEFT)
        self.mode_var.trace_add('write', lambda *_: self.set_mode())

        ttk.Label(top, text="Error Bars:").pack(side=LEFT, padx=(10, 0))
        self.errorbar_var = ttk.StringVar(value=self.errorbar_type)
        error_mode_select = ttk.Combobox(top, textvariable=self.errorbar_var,
                                         values=['calibrated', 'instrumental', 'relative', 'none'],
                                         state='readonly', width=14)
        error_mode_select.pack(side=LEFT)
        self.errorbar_var.trace_add('write', lambda *_: self.set_errorbar_type())

        ttk.Label(top, text="Time Axis:").pack(side=LEFT, padx=(10, 0))
        self.time_var = ttk.StringVar(value=self.time_mode)
        time_mode_select = ttk.Combobox(top, textvariable=self.time_var,
                                        values=['minutes', 'julian_date', 'mjd', 'rotation_phase'],
                                        state='readonly', width=18)
        time_mode_select.pack(side=LEFT)
        self.time_var.trace_add('write', lambda *_: self.set_time_mode())

        self.toggle_rejected_var = ttk.BooleanVar(value=self.show_rejected)

        self.filter_btn = ttk.Menubutton(top, text="Filters: All")
        self.filter_menu = tk.Menu(self.filter_btn, tearoff=0)
        self.filter_btn["menu"] = self.filter_menu
        self.filter_btn.pack(side=LEFT, padx=(10, 0))

        self.plot_settings_btn = ttk.Menubutton(top, text="Plot Settings")
        self.plot_settings_menu = tk.Menu(self.plot_settings_btn, tearoff=0)
        try:
            self.plot_settings.build_menu(parent=self.plot_settings_menu, pack_btn=False)
            self.plot_settings_menu.add_cascade(label="Colors", menu=self.plot_settings.color_menu)
        except Exception:
            pass
        self.plot_settings_btn["menu"] = self.plot_settings_menu
        self.plot_settings_btn.pack(side=LEFT, padx=(10, 0))

        self.plot_settings_menu.add_command(label="Marker Settings...", command=self.show_marker_settings)
        self.plot_settings_menu.add_separator()

        self.show_legend_var = ttk.BooleanVar(value=self.show_legend)
        self.show_grid_var = ttk.BooleanVar(value=self.show_grid)

        def _on_toggle_legend():
            self.show_legend = bool(self.show_legend_var.get())
            self.request_plot_update()

        def _on_toggle_grid():
            self.show_grid = bool(self.show_grid_var.get())
            self.request_plot_update()

        self.plot_settings_menu.add_checkbutton(label="Show Legend", variable=self.show_legend_var,
                                                command=_on_toggle_legend)
        self.plot_settings_menu.add_checkbutton(label="Show Grid", variable=self.show_grid_var,
                                                command=_on_toggle_grid)
        self.plot_settings_menu.add_checkbutton(label="Show Rejected", variable=self.toggle_rejected_var,
                                                command=self.set_show_rejected)

        ttk.Button(top, text="Rename Rejected", command=self.rename_rejected_files).pack(side=LEFT, padx=(10, 0))
        ttk.Button(top, text="Help (F1)", command=self.show_help).pack(side=LEFT, padx=(10, 0))
        self.root.bind("<F1>", lambda e: self.show_help())

        try:
            self.root.after(0, lambda: top._on_configure(None))
        except Exception:
            pass

        # Wrapper for side-by-side collapsible menus
        self.menus_wrapper = ttk.Frame(self.center_pane)
        self.menus_wrapper.pack(side=TOP, fill=X, padx=6, pady=(4, 0))

        # Collapsible Alignment / Colors Frame
        self.align_container = ttk.Frame(self.menus_wrapper)
        self.align_container.pack(side=LEFT, anchor='n', padx=(0, 10))
        
        # Header (Toggle Button + Label look)
        self.align_visible = tk.BooleanVar(value=False)
        
        header_frame = ttk.Frame(self.align_container)
        header_frame.pack(side=TOP, fill=X)
        
        self.toggle_btn = ttk.Button(header_frame, text="[ + ] Alignment & Colors", command=self.toggle_align_frame, style='Link.TButton')
        self.toggle_btn.pack(side=LEFT, anchor='w')
        
        # Inner Frame for Content
        self.align_inner = ttk.Frame(self.align_container, padding=6, borderwidth=1, relief="groove")
        # Initially HIDDEN
        
        # Use align_inner as parent for all prev widgets
        align_frame = self.align_inner
        
        # Merged Row: Alignment & Colors + Manual Offsets
        
        # 1. Alignment controls
        ttk.Label(align_frame, text="Ref Filter:").pack(side=LEFT, padx=(0, 5))
        self.align_ref_combo = ttk.Combobox(align_frame, textvariable=self.align_reference_var, 
                                            state='readonly', width=5)
        self.align_ref_combo.pack(side=LEFT, padx=(0, 5))
        
        ttk.Button(align_frame, text="Auto-Align", command=self.auto_align_bands).pack(side=LEFT, padx=5)
        ttk.Button(align_frame, text="Colors", command=self.show_colors_dialog).pack(side=LEFT, padx=5)
        
        self.show_colors_on_plot_var = tk.BooleanVar(value=False)
        def _toggle_colors_on_plot():
            self.request_plot_update()
            
        ttk.Checkbutton(align_frame, text="Show", variable=self.show_colors_on_plot_var, 
                        command=_toggle_colors_on_plot).pack(side=LEFT, padx=(5, 2))
                        
        loc_vals = ['upper left', 'upper right', 'lower left', 'lower right', 'center']
        ttk.Combobox(align_frame, textvariable=self.color_legend_loc_var, values=loc_vals, 
                     state='readonly', width=10).pack(side=LEFT, padx=(0, 10))
        def _on_loc_change(*_):
            self.request_plot_update()
        self.color_legend_loc_var.trace_add('write', _on_loc_change)

        # 2. Manual Offsets (Inline)
        ttk.Separator(align_frame, orient=VERTICAL).pack(side=LEFT, fill=Y, padx=5, pady=2)
        
        ttk.Label(align_frame, text="Offset:").pack(side=LEFT, padx=(5, 0))
        self.offset_var = ttk.DoubleVar(value=0.0)
        self.offset_entry_var = ttk.StringVar(value=f"{0.0:.3f}")
        
        # Compact slider
        self.offset_slider = ttk.Scale(align_frame, from_=-5.0, to=+5.0, variable=self.offset_var,
                                       orient=HORIZONTAL, length=80, command=lambda v: self.on_offset_change(v))
        self.offset_slider.pack(side=LEFT, padx=3)
        
        ttk.Entry(align_frame, textvariable=self.offset_entry_var, width=6).pack(side=LEFT)
        
        self.offset_scope_var = ttk.StringVar(value='auto')
        self.offset_scope_combo = ttk.Combobox(align_frame, textvariable=self.offset_scope_var, state='readonly', width=6)
        self.offset_scope_combo.pack(side=LEFT, padx=3)
        self.offset_scope_combo['values'] = ['auto', 'ALL']
        self.offset_scope_combo.set('auto')
        
        ttk.Label(align_frame, text="Step:").pack(side=LEFT, padx=(5,0))
        self.offset_step_var = ttk.DoubleVar(value=0.05)
        ttk.Entry(align_frame, textvariable=self.offset_step_var, width=5).pack(side=LEFT)

        # Collapsible Rotation Period Frame
        self.period_container = ttk.Frame(self.menus_wrapper)
        self.period_container.pack(side=LEFT, anchor='n')
        
        # Header (Toggle Button + Label look)
        self.period_visible = tk.BooleanVar(value=False)
        
        p_header_frame = ttk.Frame(self.period_container)
        p_header_frame.pack(side=TOP, fill=X)
        
        self.period_toggle_btn = ttk.Button(p_header_frame, text="[ + ] Rotation Period", command=self.toggle_period_frame, style='Link.TButton')
        self.period_toggle_btn.pack(side=LEFT, anchor='w')
        
        # Inner Frame for Content
        self.period_inner = ttk.Frame(self.period_container, padding=6, borderwidth=1, relief="groove")
        # Initially HIDDEN
        
        # Period controls
        rot = self.period_inner
        ttk.Label(rot, text="Period:").pack(side=LEFT, padx=(0, 0))
        self.rotation_period_var = ttk.StringVar(value=f"{DEFAULT_PERIOD:.3f}")
        ttk.Scale(rot, from_=0.1, to=50.0, variable=self.rotation_period_var,
                  orient=HORIZONTAL, length=120, command=lambda _v: self.update_rotation_period()).pack(side=LEFT, padx=3)
        
        e = ttk.Entry(rot, textvariable=self.rotation_period_var, width=7)
        e.pack(side=LEFT, padx=3)
        e.bind("<Return>", lambda _e: self.update_rotation_period())
        
        ttk.Label(rot, text="Step:").pack(side=LEFT, padx=(5, 0))
        self.rotation_step_var = ttk.DoubleVar(value=TIME_STEP)
        ttk.Entry(rot, textvariable=self.rotation_step_var, width=5).pack(side=LEFT)

        # Plot area
        plot_frame = ttk.Frame(self.center_pane)
        plot_frame.pack(side=TOP, fill=BOTH, expand=True, padx=6, pady=(4, 6))
        self.fig, self.ax = plt.subplots(figsize=(8.5, 5.2))
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)
        NavigationToolbar2Tk(self.canvas, plot_frame).update()
        self.plot = LightCurvePlot(self.fig, self.ax, self.canvas)
        self.plot.selection = {'alias': None, 'index': None}
        self.canvas.mpl_connect("button_press_event", self.on_click)
        try:
            self.canvas.mpl_connect("pick_event", self.on_pick_label)
        except Exception:
            pass
        self.plot.parent_gui = self

        # Bottom controls
        bottom = ttk.Frame(self.bottom_control)
        bottom.pack(side=BOTTOM, fill=X)
        try:
            bottom.pack_propagate(False)
            bottom.configure(height=60)
        except Exception:
            pass
        
        # Container for controls line
        ctrl = ttk.Frame(bottom)
        ctrl.pack(side=TOP, fill=X, pady=(2, 0))
        # Note: Period controls moved to top collapsible menu

        # Hotkeys Label
        self.hotkeys_label = ttk.Label(bottom,
                                       text="F1-help; q-quit; r-reject; a-cancel sel; S-show/hide image; <-/-> -move; z/x-period; s/d-offset")
        self.selected_point_label = ttk.Label(self.hotkeys_label, text="", foreground="black")
        self.selected_point_label.pack(side=RIGHT, padx=(15, 0))
        self.hotkeys_label.pack(side=BOTTOM, fill=X, padx=5, pady=(2, 2))
        self.master_frame.bind("<Configure>", lambda e: self.hotkeys_label.config(wraplength=e.width - 200))

    def open_y_scale_dialog(self):
        """Open a custom dialog to set Y-axis Min, Max, Nticks, and Label."""
        if not hasattr(self, 'plot'):
            return

        # Get current values
        try:
            lo, hi = self.plot.ax.get_ylim()
            # Inverted check: usually lo > hi for magnitudes
            cur_min = min(lo, hi)
            cur_max = max(lo, hi)
        except Exception:
            cur_min, cur_max = 0.0, 1.0
        
        cur_nticks = self.plot.y_nticks if self.plot.y_nticks else ""
        cur_label = getattr(self.plot, 'ylabel', 'Magnitude')

        # Create dialog
        dlg = tk.Toplevel(self.root)
        dlg.title("Y-Axis Settings")
        dlg.transient(self.root)
        dlg.resizable(False, False)
        
        # Center dialog
        try:
            x = self.root.winfo_x() + 150
            y = self.root.winfo_y() + 150
            dlg.geometry(f"+{x}+{y}")
        except Exception:
            pass

        frm = ttk.Frame(dlg, padding=10)
        frm.pack(fill=BOTH, expand=True)

        # Label
        ttk.Label(frm, text="Label:").grid(row=0, column=0, sticky='e', padx=5, pady=5)
        label_var = ttk.StringVar(value=cur_label)
        ttk.Entry(frm, textvariable=label_var, width=15).grid(row=0, column=1, padx=5, pady=5)

        # Min
        ttk.Label(frm, text="Min:").grid(row=1, column=0, sticky='e', padx=5, pady=5)
        min_var = ttk.StringVar(value=f"{cur_min:.3f}")
        ttk.Entry(frm, textvariable=min_var, width=10).grid(row=1, column=1, padx=5, pady=5)

        # Max
        ttk.Label(frm, text="Max:").grid(row=2, column=0, sticky='e', padx=5, pady=5)
        max_var = ttk.StringVar(value=f"{cur_max:.3f}")
        ttk.Entry(frm, textvariable=max_var, width=10).grid(row=2, column=1, padx=5, pady=5)

        # Nticks
        ttk.Label(frm, text="Ticks (optional):").grid(row=3, column=0, sticky='e', padx=5, pady=5)
        nticks_var = ttk.StringVar(value=str(cur_nticks))
        ttk.Entry(frm, textvariable=nticks_var, width=10).grid(row=3, column=1, padx=5, pady=5)

        def apply():
            try:
                v_min = float(min_var.get())
                v_max = float(max_var.get())
                
                nt_str = nticks_var.get().strip()
                nt = int(nt_str) if nt_str else None

                self.plot.ylabel = label_var.get()
                
                # Apply
                self.plot.set_y_limits(min(v_min, v_max), max(v_min, v_max), nt)
                self._refresh_plot_callback()
                dlg.destroy()
            except ValueError:
                messagebox.showerror("Invalid Input", "Please enter valid numbers.", parent=dlg)

        def auto():
            self.plot.ylabel = label_var.get()
            self.plot.set_y_limits(None, None, None)
            self._refresh_plot_callback()
            dlg.destroy()

        btn_frm = ttk.Frame(frm)
        btn_frm.grid(row=4, column=0, columnspan=2, pady=(10, 0))
        
        ttk.Button(btn_frm, text="Auto / Reset", command=auto).pack(side=LEFT, padx=5)
        ttk.Button(btn_frm, text="Apply", command=apply, style="Accent.TButton").pack(side=LEFT, padx=5)

        dlg.bind("<Return>", lambda e: apply())
        dlg.bind("<Escape>", lambda e: dlg.destroy())

    # Helpers for collapsing/expanding


    def _build_save_menu(self, parent):
        self.save_menu_btn = ttk.Menubutton(parent, text="Save")
        self.save_menu = tk.Menu(self.save_menu_btn, tearoff=0, font=("TkDefaultFont", 10))
        self.save_menu.add_command(label="Save CSV", command=self.save_csv)
        self.save_menu.add_command(label="Save Plot", command=self.save_plot)
        self.save_menu.add_command(label="Save Atlas", command=self.save_atlas)
        self.save_menu_btn["menu"] = self.save_menu
        self.save_menu_btn.pack(side=LEFT, padx=4)


    def show_help(self):
        msg = (
            "Hotkeys & Controls\n"
            "--------------------\n"
            "q              Quit application\n"
            "r              Toggle rejection for selected point\n"
            "a              Cancel/clear selection\n"
            "S              Show / Hide image for selected point\n"
            "o              Toggle overlay of the image\n"
            "<- / ->        Move selected point left/right\n"
            "s / d           Change vertical offset\n"
            "z / x          Decrease / Increase rotation period\n"
            "\n"
            "Mouse:\n"
            "Click on a point to select it.\n"
            "\n"
            "Options:\n"
            "- Mode: target / instrumental / control.\n"
            "- Error Bars: calibrated / instrumental / none.\n"
            "- Time Axis: minutes / JD / MJD / rotation phase.\n"
            "- Rotation Step (h): sets increment used by z/x.\n"
            "- Offset Step (mag): sets increment used by s/d.\n"

        )

        top = tk.Toplevel(self.root)
        top.title("Help")
        top.attributes("-topmost", True)
        top.resizable(True, True)

        txt = tk.Text(top, wrap="word", height=24, width=70)
        txt.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        txt.insert("1.0", msg)
        txt.configure(state="disabled")
        ttk.Button(top, text="Close", command=top.destroy).pack(pady=(0, 10))

    def _ensure_active_alias(self):
        """If nothing is selected but there is a single LC, select it."""
        if self.current_lc_alias:
            return self.current_lc_alias
        if len(self.lightcurves) == 1:
            alias = next(iter(self.lightcurves))
            self.current_lc_alias = alias
            return alias
        return None

    # ——— file ops ———
    def open_csv(self):
        files = filedialog.askopenfilenames(initialdir=os.getcwd(), filetypes=[("CSV Files", "*.csv")])
        if not files:
            return
        
        # Single-file mode: clear existing data
        self.lightcurves.clear()
        self.lc_visible.clear()
        self.lc_offsets.clear()
        
        # Load only the first selected file
        f = files[0]
        try:
            lc = LightCurveData(f)
        except Exception as e:
            messagebox.showerror("Load error", f"Failed to load {f}: {e}")
            return
            
        alias = self._unique_alias(os.path.basename(f))
        self.lightcurves[alias] = lc
        self.lc_visible[alias] = True
        
        # initialize per‑file offsets map
        bands = lc.get_bands()
        self.lc_offsets[alias] = {OFFSET_ALL_KEY: 0.0}
        for b in bands:
            self.lc_offsets[alias][b] = 0.0
            
        # Select this lightcurve
        self.current_lc_alias = alias
        
        self._rebuild_band_menu() if hasattr(self, "_rebuild_band_menu") else None
        
        # Sync offset UI
        self._sync_offset_var_to_scope()
        self._refresh_offset_scope_choices()

        # also refresh color menu to show newly available bands
        self._rebuild_color_menu() if hasattr(self, "_rebuild_color_menu") else None
        
        self.request_plot_update()

    def _unique_alias(self, base: str) -> str:
        alias = base
        k = 1
        while alias in self.lightcurves:
            stem, ext = os.path.splitext(base)
            alias = f"{stem}_{k}{ext}"
            k += 1
        return alias

    def save_csv(self):
        if not self.current_lc_alias:
            messagebox.showerror("Nothing to save", "Select a lightcurve entry to save.")
            return
        lc = self.lightcurves.get(self.current_lc_alias)
        if lc is None:
            messagebox.showerror("Nothing to save", "No lightcurve selected.")
            return
        file = filedialog.asksaveasfilename(
            initialdir=os.path.dirname(lc.filename) if lc.filename else os.getcwd(),
            initialfile=os.path.basename(lc.filename) if lc.filename else "lightcurve.csv",
            defaultextension=".csv",
            filetypes=[("CSV Files", "*.csv")],
        )
        if not file:
            file = next_version(lc.filename or "lightcurve.csv")
        lc.filename = file
        lc.save()
        messagebox.showinfo("Saved", f"CSV saved as {file}")

    def save_plot(self):
        if not self.lightcurves:
            messagebox.showerror("Nothing to save", "No data to save.")
            return
        # name the same as input csv file
        lc = self.lightcurves.get(self.current_lc_alias)
        default_png = Path(lc.filename).stem + ".png" if lc.filename else "lightcurve.png"
        file = filedialog.asksaveasfilename(initialfile=default_png, defaultextension=".png",
                                            filetypes=[("PNG Files", "*.png")])
        if not file:
            file = next_version(default_png)
        self.plot.fig.savefig(file, dpi=300, bbox_inches="tight")
        messagebox.showinfo("Saved", f"Plot saved as {file}")

    def save_atlas(self):
        if not self.current_lc_alias:
            messagebox.showerror("Nothing to save", "Select a lightcurve entry to save atlas for.")
            return
        lc = self.lightcurves.get(self.current_lc_alias)
        if lc is None or lc.df is None:
            messagebox.showerror("Nothing to save", "Selected lightcurve has no data.")
            return

        filetypes = [("FITS files", "*.fits *.fts"), ("All files", "*.*")]
        fits_filepath = filedialog.askopenfilename(title="Select Atlas FITS file", filetypes=filetypes)
        if not fits_filepath:
            messagebox.showerror("Nothing selected", "Select atlas FITS file")
            return
        atlas_name = f"{Path(lc.filename).stem}.ATL" if lc.filename else "output.ATL"
        file = filedialog.asksaveasfilename(initialfile=atlas_name, defaultextension=".ATL",
                                            filetypes=[("Atlas Files", "*.ATL")])
        if not file:
            messagebox.showerror("Nothing selected", "No atlas filename provided")
            return

        # check if multiple filters are present
        # Save all bands present in the file, not just the selected ones
        filters = lc.get_bands()
        tmps = []
        tmp_atls = []
        if len(filters) > 1:
            # create temp csvs for every filter and call pp_atlas for each
            for f in filters:
                df_filt = lc.df[lc.df['band'] == f]
                if df_filt.empty:
                    continue
                tmpf = f"tmp_atlas_{f}.csv"
                tmp_atl = f"{tmpf[:-4]}.ATL"
                tmp_atls.append(tmp_atl)
                df_filt.to_csv(tmpf, header=True, index=False)
                tmps.append(tmpf)
                atlas_cmd = f"pp_atlas -fname_header {fits_filepath} -fname_photo {tmpf} -fname_out {tmpf[:-4]}.ATL"
                subprocess.call(['/bin/sh', '-i', '-c', atlas_cmd])
            # use combine altas to get one altas file
            atlas_cmd = f"pp_atlas -combine " + " ".join([f"{tmp_file[:-4]}.ATL" for tmp_file in tmps]) + f" -fname_out {file}"
            try:
                subprocess.call(['/bin/sh', '-i', '-c', atlas_cmd])
            except Exception as e:
                messagebox.showerror("Atlas save failed", f"Error calling atlas command: {e}")
            # remove temp atlas files after combining
            finally:
                for tmp_atl, tmp_csv in zip(tmp_atls, tmps):
                    if os.path.exists(tmp_atl):
                        os.remove(tmp_atl)
                    if os.path.exists(tmp_csv):
                        os.remove(tmp_csv)
        else:
            tmp = 'tmp_atlas.csv'
            tmps.append(tmp)
            lc.df.to_csv(tmp, header=True, index=False)
            atlas_cmd = f"pp_atlas -fname_header {fits_filepath} -fname_photo {tmp} -fname_out {file}"
            try:
                subprocess.call(['/bin/sh', '-i', '-c', atlas_cmd])
            except Exception as e:
                messagebox.showerror("Atlas save failed", f"Error calling atlas command: {e}")
            finally:
                # cleanup temp files
                for tmp in tmps:
                    if os.path.exists(tmp):
                        os.remove(tmp)
            messagebox.showinfo("Saved", f"ATLAS file saved as {file}")

    def rename_rejected_files(self):
        """
        Renames files associated with rejected points to have a '._fits' extension,
        and ensures accepted points have a '.fits' extension.
        Searches for files in subdirectories relative to the CSV file location.
        """
        if not self.current_lc_alias:
            messagebox.showwarning("No Data", "No lightcurve loaded.")
            return

        lc = self.lightcurves[self.current_lc_alias]
        if lc.df is None or lc.df.empty:
            messagebox.showwarning("No Data", "Lightcurve has no data.")
            return

        if 'filename' not in lc.df.columns:
            messagebox.showwarning("Missing Data", "Dataframe is missing 'filename' column.")
            return

        # Confirm action
        if not messagebox.askyesno("Rename Files",
                                   "This operation will rename files on disk based on their rejection status:\n"
                                   "- Rejected points -> ._fits\n"
                                   "- Accepted points -> .fits\n\n"
                                   "The program will search for files matching the basenames in subdirectories (recursive) and extensions will be changed.\n"
                                   "Continue?"):
            return

        # Determine search root
        search_root = os.path.dirname(os.path.abspath(lc.filename)) if lc.filename else os.getcwd()

        # Only interested in extensions .fits, ._fits, .fts, ._fts
        relevant_exts = {'.fits', '._fits', '.fts', '._fts'}
        file_map = {}  # stem -> { ext -> [full_paths] }

        # Show busy cursor
        self.root.config(cursor="watch")
        self.root.update_idletasks()

        try:
            # 1. Index all files in subdirectories
            for root, dirs, files in os.walk(search_root):
                for name in files:
                    # check extension
                    p = Path(name)
                    ext = p.suffix
                    stem = p.stem
                    
                    if ext in relevant_exts:
                        if stem not in file_map:
                            file_map[stem] = {}
                        if ext not in file_map[stem]:
                            file_map[stem][ext] = []
                        file_map[stem][ext].append(os.path.join(root, name))

            renamed_count = 0
            errors = 0

            # 2. Iterate through DataFrame
            for idx, row in lc.df.iterrows():
                fname = str(row['filename'])
                is_rejected = bool(row['rejected'])

                stem = Path(fname).stem
                
                # If the file described in DF is not found in scan, skip
                if stem not in file_map:
                    continue

                variants = file_map[stem]
                
                # Check for existence of valid vs rejected variants
                # Prioritize .fits over .fts, ._fits over ._fts
                found_valid_list = variants.get('.fits', []) + variants.get('.fts', [])
                found_rejected_list = variants.get('._fits', []) + variants.get('._fts', [])

                if is_rejected:
                    # Ensure file is named ._fits
                    # If have valid files, rename them to ._fits
                    for old_path in found_valid_list:
                        dir_name = os.path.dirname(old_path)
                        name = os.path.basename(old_path)
                        
                        if name.endswith('.fits'):
                            new_name = name.replace('.fits', '._fits')
                        elif name.endswith('.fts'):
                             new_name = name.replace('.fts', '._fts')
                        else:
                             new_name = name + "._fits"
                        
                        new_path = os.path.join(dir_name, new_name)
                        try:
                            if os.path.exists(new_path):
                                print(f"Target {new_path} exists. Skipping rename of {old_path}")
                                continue
                            os.rename(old_path, new_path)
                            renamed_count += 1
                        except OSError as e:
                            print(f"Error renaming {old_path}: {e}")
                            errors += 1
                    
                else:
                    # Ensure file is named .fits 
                    # If there are rejected files, rename them to .fits
                    for old_path in found_rejected_list:
                        dir_name = os.path.dirname(old_path)
                        name = os.path.basename(old_path)
                        
                        if name.endswith('._fits'):
                            new_name = name.replace('._fits', '.fits')
                        elif name.endswith('._fts'):
                             new_name = name.replace('._fts', '.fts')
                        else:
                             base, ext = os.path.splitext(name)
                             if ext.startswith('._'):
                                 new_name = base + '.' + ext[2:]
                             else:
                                 new_name = base + '.fits'
                        
                        new_path = os.path.join(dir_name, new_name)
                        try:
                            if os.path.exists(new_path):
                                print(f"Target {new_path} exists. Skipping rename of {old_path}")
                                continue
                            os.rename(old_path, new_path)
                            renamed_count += 1
                        except OSError as e:
                            print(f"Error renaming {old_path}: {e}")
                            errors += 1

            messagebox.showinfo("Rename Complete", f"Renamed {renamed_count} files.\nErrors: {errors}")

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred during renaming: {e}")
        finally:
            self.root.config(cursor="")

    # ——— Alignment & Colors ———
    def _update_align_ref_choices(self):
        """Update the reference band combobox based on available bands."""
        if not self.lightcurves:
            self.align_ref_combo['values'] = []
            return
            
        # Collect bands from the current alias, or all aliases
        alias = self.current_lc_alias
        if alias and alias in self.lightcurves:
            bands = self.lightcurves[alias].get_bands()
        else:
            bands = sorted({str(b) for lc in self.lightcurves.values() for b in lc.get_bands()})
            
        self.align_ref_combo['values'] = bands
        if bands and self.align_reference_var.get() not in bands:
            self.align_reference_var.set(bands[0])

    def auto_align_bands(self):
        """Align other bands to the reference band using linear interpolation."""
        alias = self._ensure_active_alias()
        if not alias:
            messagebox.showwarning("No Data", "No lightcurve loaded/selected.")
            return

        ref_band = self.align_reference_var.get()
        if not ref_band:
            messagebox.showwarning("No Reference", "Please select a Reference Band.")
            return

        lc = self.lightcurves[alias]
        df = lc.df
        if df is None or df.empty:
            return

        # Prepare storage for colors
        # Only rewrite if we align properly
        if alias not in self.calculated_colors:
            self.calculated_colors[alias] = {}
        # clear this alias's colors because we are re-aligning
        self.calculated_colors[alias] = {}
        
        jd_all = lc.arr['jd']
        mag_all = lc.arr['mag']
        sig_all = lc.arr['sig']
        band_all = lc.arr['band'].astype(str)
        
        ref_mask = (band_all == ref_band)
        if not np.any(ref_mask):
            messagebox.showwarning("Error", f"Reference band {ref_band} has no data.")
            return

        # Sort reference data by time
        ref_idx = np.argsort(jd_all[ref_mask])
        ref_t = jd_all[ref_mask][ref_idx]
        ref_m = mag_all[ref_mask][ref_idx]
        ref_e = sig_all[ref_mask][ref_idx]

        # Which bands to align? All others.
        other_bands = [b for b in lc.get_bands() if b != ref_band]
        
        aligned_count = 0
        
        for b in other_bands:
            b_mask = (band_all == b)
            if not np.any(b_mask):
                continue
                
            b_t = jd_all[b_mask]
            b_m = mag_all[b_mask]
            b_e = sig_all[b_mask]
            
            
            # Use JD (Time) for alignment
            # If the user took data sequentially (R, V, R, V...), interpolation in time is correct.
            
            ref_m_interp = np.interp(b_t, ref_t, ref_m, left=np.nan, right=np.nan)
            
            # Calculate differences (residuals)
            # If B + Offset = Ref, then Offset = Ref - B.
            # Ref - B is also the color index (Ref-B).
            
            diff = ref_m_interp - b_m
            
            # Propagate errors: 
            # Let's interp variance for better estimate.
            ref_var = ref_e**2
            ref_var_interp = np.interp(b_t, ref_t, ref_var, left=np.nan, right=np.nan)
            
            var_diff = b_e**2 + ref_var_interp
            
            # Mask NaNs (points outside Ref range)
            valid = np.isfinite(diff) & np.isfinite(var_diff)
            
            if np.sum(valid) < 2:
                # Not enough overlapping data points
                continue
                
            # Weighted Mean
            d = diff[valid]
            v = var_diff[valid]
            weights = 1.0 / v
            
            w_mean = np.sum(d * weights) / np.sum(weights)
            w_err = np.sqrt(1.0 / np.sum(weights))
            
            # Also calculate StdDev of residuals to see if there is intrinsic variation
            # (e.g. color change during rotation)
            std_dev = np.std(d)
            
            # Use the larger of w_err or std_dev/sqrt(N)? 
            # Standard error of weighted mean is w_err.
            
            # Set offset
            # We want B + Offset = Ref => Offset as computed.
            # But wait, logic in plotting is: Y_plot = Ybase + off
            # So if we set off = w_mean, Y_plot = B + (Ref - B) = Ref. Correct.
            
            # Store offset
            self.lc_offsets[alias][b] = round(w_mean, 4)
            
            # Store Color Info
            # Color Indicies are usually M_band1 - M_band2.
            # Here we have Offset = Ref - B. 
            # So the color is (Ref - B).
            color_name = f"{ref_band}-{b}"
            self.calculated_colors[alias][color_name] = (w_mean, w_err)
            
            aligned_count += 1
            
        # Ensure Ref band offset is 0 
        self.lc_offsets[alias][ref_band] = 0.0
        self.lc_offsets[alias][OFFSET_ALL_KEY] = 0.0
        
        # Refresh
        self.request_plot_update()
        self._sync_offset_var_to_scope()
        
        msg = f"Aligned {aligned_count} filters to {ref_band}."
        if aligned_count > 0:
            msg += "\nColors have been calculated (View Colors)."
        messagebox.showinfo("Alignment", msg)

    def show_colors_dialog(self):
        alias = self.current_lc_alias
        if not alias or alias not in self.calculated_colors or not self.calculated_colors[alias]:
            messagebox.showinfo("Colors", "No calculated colors available.\nRun 'Auto-Align' first.")
            return
            
        colors = self.calculated_colors[alias]
        
        dlg = tk.Toplevel(self.root)
        dlg.title(f"Color Indices: {alias}")
        dlg.transient(self.root)
        
        try:
            dlg.geometry(f"+{self.root.winfo_x()+100}+{self.root.winfo_y()+100}")
        except:
            pass
            
        # Table
        cols = ("Color", "Value", "Error")
        tree = ttk.Treeview(dlg, columns=cols, show='headings', height=len(colors)+2)
        for c in cols:
            tree.heading(c, text=c)
            tree.column(c, width=100, anchor='center')
        tree.pack(fill=BOTH, expand=True, padx=10, pady=10)
        
        for name, (val, err) in colors.items():
            tree.insert("", END, values=(name, f"{val:.4f}", f"{err:.4f}"))
            
        # Buttons
        btn_frame = ttk.Frame(dlg)
        btn_frame.pack(fill=X, pady=10)
        
        def save_colors():
            fname = filedialog.asksaveasfilename(defaultextension=".txt", 
                                                 initialfile=f"{Path(self.lightcurves[alias].filename).stem}_colors.txt",
                                                 filetypes=[("Text Files", "*.txt")])
            if fname:
                try:
                    with open(fname, 'w') as f:
                        f.write(f"# Color indices for {alias}\n")
                        f.write(f"# Computed relative to reference filter \n")
                        f.write(f"Color | Value | Error\n")
                        for name, (val, err) in colors.items():
                            f.write(f"{name} | {val:.6f} | {err:.6f}\n")
                    messagebox.showinfo("Saved", f"Colors saved to {fname}")
                except Exception as e:
                    messagebox.showerror("Error", f"Failed to save: {e}")

        ttk.Button(btn_frame, text="Save to File", command=save_colors).pack(side=LEFT, padx=20)
        ttk.Button(btn_frame, text="Close", command=dlg.destroy).pack(side=RIGHT, padx=20)

    def toggle_align_frame(self):
        is_vis = self.align_visible.get()
        new_vis = not is_vis
        self.align_visible.set(new_vis)
        
        if new_vis:
            self.toggle_btn.config(text="[ - ] Alignment & Colors")
            self.align_inner.pack(side=TOP, fill=X)
        else:
            self.toggle_btn.config(text="[ + ] Alignment & Colors")
            self.align_inner.pack_forget()

    def toggle_period_frame(self):
        is_vis = self.period_visible.get()
        new_vis = not is_vis
        self.period_visible.set(new_vis)
        
        if new_vis:
            self.period_toggle_btn.config(text="[ - ] Rotation Period")
            self.period_inner.pack(side=TOP, fill=X)
        else:
            self.period_toggle_btn.config(text="[ + ] Rotation Period")
            self.period_inner.pack_forget()

    def on_offset_change(self, value):
        if not self._ensure_active_alias():
            return
        # value may come from scale (string) or entry (string)
        try:
            v = float(value)
        except Exception:
            v = 0.0
        # round to 3 decimals and store
        v = round(float(v), 3)
        # If you have per-band offsets:
        band_key = self._active_offset_band() if hasattr(self, "_active_offset_band") else None
        if band_key is not None and isinstance(self.lc_offsets.get(self.current_lc_alias), dict):
            self.lc_offsets[self.current_lc_alias][band_key] = v
        else:
            # ensure we store a mapping for this alias
            self.lc_offsets[self.current_lc_alias] = {OFFSET_ALL_KEY: float(v)}
        # reflect rounded value back into the slider and formatted entry
        try:
            self.offset_var.set(v)
        except Exception:
            pass
        try:
            self.offset_entry_var.set(f"{v:.3f}")
        except Exception:
            pass
        self.request_plot_update()
        
        # Sync with calculated colors if exists
        self._sync_manual_offset_to_colors(band_key, v)

        # draw immediately (fast blit)
        if hasattr(self.plot, "blit"):
            self.plot.blit.quick_redraw()

    def adjust_lc_offset(self, delta: float):
        if not self._ensure_active_alias():
            return
        # If you have per-band offsets:
        band_key = self._active_offset_band() if hasattr(self, "_active_offset_band") else None
        if band_key is not None and isinstance(self.lc_offsets.get(self.current_lc_alias), dict):
            cur = float(self.lc_offsets[self.current_lc_alias].get(band_key, 0.0))
            newv = round(cur + delta, 3)
            self.lc_offsets[self.current_lc_alias][band_key] = newv
        else:
            cur = float(self.lc_offsets.get(self.current_lc_alias, 0.0))
            newv = round(cur + delta, 3)
            # ensure alias has a band-mapping
            self.lc_offsets[self.current_lc_alias] = {OFFSET_ALL_KEY: newv}
        # update numeric slider var and formatted entry
        try:
            self.offset_var.set(newv)
        except Exception:
            pass
        try:
            self.offset_entry_var.set(f"{newv:.3f}")
        except Exception:
            pass
        self.request_plot_update()
        
        # Sync with calculated colors if exists
        self._sync_manual_offset_to_colors(band_key, newv)

        if hasattr(self.plot, "blit"):
            self.plot.blit.quick_redraw()


    def _sync_manual_offset_to_colors(self, band, offset_val):
        """Update the calculated color value for a band if the user manually changes its offset."""
        alias = self.current_lc_alias
        if not alias or alias not in self.calculated_colors:
            return
            
        # need to find which color corresponds to this band.
        # Exception: if band is the Reference itself? Usually Ref offset is 0.
        
        suffix = f"-{band}"
        
        # If we have a direct match (Ref-B), update it.
        # Note: offsets are stored as simple floats.
        for color_name in self.calculated_colors[alias]:
            if color_name.endswith(suffix):
                # This color is Ref-Band
                # The value should be the offset.
                old_val, err = self.calculated_colors[alias][color_name]
                self.calculated_colors[alias][color_name] = (offset_val, err)
                # Force plot info update (since color value changed txt)
                # But request_plot_update was already called. 
                # layout_sig includes loc/flags but might not include color VALUES if they are just text.
                # However, the plot logic rebuilds the text box every update() call.
                return


    # ——— Offsets ———
    def _refresh_offset_scope_choices(self):
        alias = self.current_lc_alias
        bands = []
        if alias and alias in self.lightcurves:
            bands = self.lightcurves[alias].get_bands()
        choices = ['ALL'] + bands
        self.offset_scope_combo['values'] = choices
        if self.offset_scope_var.get() not in choices:
            self.offset_scope_combo.set('ALL')

    def _active_offset_band(self) -> str:
        """Resolve user's current scope into a concrete band key or OFFSET_ALL_KEY."""
        scope = self.offset_scope_var.get()
        if scope == 'ALL':
            return OFFSET_ALL_KEY
        if scope != 'ALL':
            return scope  # explicit band chosen by user
        # auto: prefer selected point's band, else single visible band, else ALL
        if self.current_point_band:
            return self.current_point_band
        if len(self.selected_bands) == 1:
            return list(self.selected_bands)[0]
        return OFFSET_ALL_KEY

    def _get_offset_step(self) -> float:
        """Return the configured offset step (mag), validated and rounded."""
        try:
            step = float(self.offset_step_var.get())
            if not np.isfinite(step) or step <= 0:
                step = 0.05
        except Exception:
            step = 0.05
        return round(step, 4)

    def _sync_offset_var_to_scope(self):
        if not self.current_lc_alias:
            self.offset_var.set(0.0)
            self.offset_entry_var.set(f"{0.0:.3f}")
            return
        band_key = self._active_offset_band()
        off = self.lc_offsets.get(self.current_lc_alias, {}).get(band_key, 0.0)
        try:
            offf = float(off)
        except Exception:
            offf = 0.0
        # show rounded to 3 decimals in the input field
        v = round(offf, 3)
        try:
            self.offset_var.set(v)
        except Exception:
            pass
        try:
            self.offset_entry_var.set(f"{v:.4f}")
        except Exception:
            pass

    # ——— State setters ——
    def set_mode(self):
        self.mode = self.mode_var.get()
        if hasattr(self.plot, "invalidate_layout"):
            self.plot.invalidate_layout()
        self.request_plot_update()
        self.plot.blit.quick_redraw()

    def set_errorbar_type(self):
        self.errorbar_type = self.errorbar_var.get()
        if hasattr(self.plot, "invalidate_layout"):
            self.plot.invalidate_layout()
        self.request_plot_update()
        self.plot.blit.quick_redraw()

    def set_time_mode(self):
        self.time_mode = self.time_var.get()
        if hasattr(self.plot, "invalidate_layout"):
            self.plot.invalidate_layout()
        self.request_plot_update()
        self.plot.blit.quick_redraw()

    def set_show_rejected(self):
        self.show_rejected = bool(self.toggle_rejected_var.get())
        self.request_plot_update()

    # ——— Filter menu ———
    def _rebuild_band_menu(self):
        """Rebuild the Filters menu and keep GUI + plotting state in sync.
        Semantics:
          - self.selected_bands == set([...]) → only those bands are shown
          - self.selected_bands == set()      → show NONE (useful for 'Select None')
        """
        # 1 Collect all bands present across loaded LCs, as strings
        bands: list[str] = sorted({str(b) for lc in self.lightcurves.values() for b in lc.get_bands()})

        # 2 Default to ALL bands selected if nothing chosen yet
        if not getattr(self, "selected_bands", None):
            # if it's empty (or missing), pick 'all'
            self.selected_bands = set(bands)

        # 3 Drop any previously selected band that no longer exists
        self.selected_bands = {b for b in self.selected_bands if b in bands}

        # If we dropped everything (e.g., re-opened different files), fall back to ALL
        if not self.selected_bands and bands:
            self.selected_bands = set(bands)

        # 4 Rebuild the BooleanVars for each band
        self.band_vars = {b: tk.BooleanVar(value=(b in self.selected_bands)) for b in bands}

        # 5 Rebuild the menu UI
        self.filter_menu.delete(0, END)

        def set_all(val: bool):
            for b in bands:
                self.band_vars[b].set(val)
            self._on_band_menu_changed()

        self.filter_menu.add_command(label="Select All", command=lambda: set_all(True))
        self.filter_menu.add_command(label="Select None", command=lambda: set_all(False))
        self.filter_menu.add_separator()

        for b in bands:
            self.filter_menu.add_checkbutton(label=b, variable=self.band_vars[b],
                                             command=self._on_band_menu_changed)

        # 6 Update button text + redraw
        self._update_filter_btn_text()
        self._update_align_ref_choices()
        self.request_plot_update()

    def _update_filter_btn_text(self):
        all_bands = set(self.band_vars.keys())
        sel = getattr(self, "selected_bands", set())
        if not sel:
            self.filter_btn.configure(text="Filters: (none)")
        elif sel == all_bands and len(all_bands) > 0:
            self.filter_btn.configure(text="Filters: All")
        elif len(sel) > 4:
            self.filter_btn.configure(text=f"Filters: {len(sel)} selected")
        else:
            self.filter_btn.configure(text="Filters: " + ",".join(sorted(sel)))

    def _on_band_menu_changed(self):
        # recompute selection set from checkboxes
        sel = {b for b, var in self.band_vars.items() if var.get()}
        self.selected_bands = sel
        self._update_filter_btn_text()

        # changing the band set changes which artists exist -> force a rebuild
        if hasattr(self.plot, "invalidate_layout"):
            self.plot.invalidate_layout()

        # request the update then do an immediate paint
        self.request_plot_update()

        # IMPORTANT: after changing filters, make sure the first paint is a full draw
        # so the blitter recaptures a correct background
        if hasattr(self.plot, "blit"):
            self.plot.blit._bg = None
            self.plot.blit.quick_redraw()

    # ——— Plot interactions ———
    def request_plot_update(self):
        self.plot.update(self.lightcurves, self.lc_offsets, self.lc_visible, self.mode, self.time_mode,
                         self.show_rejected, self.errorbar_type, self.selected_bands)

    def on_click(self, event):
        if event.inaxes != self.ax or event.xdata is None or event.ydata is None:
            return

        alias, idx, b = self.plot.select_nearest_point(
            event.xdata, event.ydata,
            self.lightcurves,
            self.lc_offsets,
            self.mode,
            self.time_mode,
            getattr(self, "selected_bands", set()) if isinstance(getattr(self, "selected_bands", None), set) else set(),
            self.show_rejected,
        )

        if alias is None or idx is None:
            self.plot.selection = {'alias': None, 'index': None}
            self.current_point_index = None
            self.current_point_band = None
        else:
            self.current_lc_alias = alias
            self.current_point_index = idx
            self.current_point_band = b
            self._update_selected_point_label()
            if self.offset_scope_var.get() == 'auto':
                self.offset_scope_combo.set(b if b else 'ALL')
            self._sync_offset_var_to_scope()
            self.plot.selection = {'alias': alias, 'index': idx}

        self.request_plot_update()
        self.asteroid_viewer.update_image()  # Update image on any click

    def toggle_rejection(self):
        if self.current_lc_alias is None or self.current_point_index is None:
            return
        lc = self.lightcurves.get(self.current_lc_alias)
        if lc is None or lc.df is None:
            return
        if 0 <= self.current_point_index < len(lc.df):
            lc.toggle_rejection(self.current_point_index)
            self.plot.selection = {'alias': self.current_lc_alias, 'index': self.current_point_index}
            if hasattr(self.plot, "invalidate_layout"):
                self.plot.invalidate_layout()
            else:
                self.plot._layout_signature = None
            self.request_plot_update()

    def cancel_selection(self):
        self.plot.selection = {'alias': None, 'index': None}
        self.current_point_index = None
        self.current_point_band = None
        self._update_selected_point_label()
        self.request_plot_update()
        self.asteroid_viewer.update_image()  # Update (clear) image window

    def move_point(self, direction: int):
        """Move selection to the next/previous point using arrow keys."""
        if self.current_lc_alias is None:
            return
        lc = self.lightcurves.get(self.current_lc_alias)
        if lc is None or lc.df is None or lc.df.empty:
            return

        n = len(lc.df)
        if self.current_point_index is None:
            # If nothing is selected, select the first or last point
            self.current_point_index = 0 if direction > 0 else n - 1
        else:
            # Move index, clamping between 0 and n-1
            self.current_point_index = max(0, min(n - 1, self.current_point_index + direction))

        # Update plot selection and redraw
        self.plot.selection['alias'] = self.current_lc_alias
        self.plot.selection['index'] = self.current_point_index
        self._update_selected_point_label()
        self.request_plot_update()

        # Update the image if the window is open
        self.asteroid_viewer.update_image()

    def _update_selected_point_label(self):
        """Update the selected point footnote label."""
        if not hasattr(self, 'selected_point_label'):
            return

        if self.current_lc_alias is None or self.current_point_index is None:
            self.selected_point_label.config(text="")
            return

        lc = self.lightcurves.get(self.current_lc_alias)
        if lc is None or lc.df is None or self.current_point_index >= len(lc.df):
            self.selected_point_label.config(text="")
            return

        try:
            filename = str(lc.df['filename'].iloc[self.current_point_index])
            self.selected_point_label.config(text=f"Selected point: {filename}")
        except Exception:
            self.selected_point_label.config(text="")

    def show_asteroid_image(self, event=None):
        self.asteroid_viewer.toggle_visibility()

    # ——— Rotation period ———
    def adjust_rotation_period(self, delta: float):
        try:
            step = float(self.rotation_step_var.get())
            if not np.isfinite(step) or step <= 0:
                step = TIME_STEP
        except Exception:
            step = TIME_STEP

        try:
            current = float(self.rotation_period_var.get())
        except Exception:
            current = DEFAULT_PERIOD

        new_value = max(0.1, current + delta * step)
        self.rotation_period_var.set(f"{new_value:.4f}")
        self.update_rotation_period()

    def update_rotation_period(self, *_):
        try:
            val = float(self.rotation_period_var.get())
        except Exception:
            val = DEFAULT_PERIOD
        self.rotation_period_var.set(f"{val:.3f}")

        if self.time_mode != 'rotation_phase':
            return
        self.request_plot_update()

    def on_pick_label(self, event):
        try:
            artist = getattr(event, 'artist', None)
            if artist is None: return
            plot = getattr(self, 'plot', None)
            if plot is None: return

            if artist is getattr(plot, 'title_text', None):
                cur = getattr(plot, 'title', '')
                new = simpledialog.askstring("Edit Title", "Enter new title:", initialvalue=cur, parent=self.root)
                if new is None: return
                plot.title = new
                plot.auto_title = False
            elif artist is getattr(plot, 'xlabel_text', None):
                cur = getattr(plot, 'xlabel', '')
                new = simpledialog.askstring("Edit X label", "Enter new X axis label:", initialvalue=cur,
                                             parent=self.root)
                if new is None: return
                plot.xlabel = new
            elif artist is getattr(plot, 'ylabel_text', None):
                self.open_y_scale_dialog()
                return
            else:
                return

            if hasattr(plot, 'invalidate_layout'):
                plot.invalidate_layout()
            self.request_plot_update()
            if hasattr(plot, 'blit') and plot.blit is not None:
                plot.blit._bg = None
                plot.blit.quick_redraw()
        except Exception:
            pass

    # ---------- marker settings dialog ----------
    def show_marker_settings(self):
        marker_dialog = ttk.Toplevel()
        marker_dialog.title("Marker and Error Bar Settings")
        marker_dialog.transient(self.root)
        marker_dialog.grab_set()
        x = self.root.winfo_x() + 150
        y = self.root.winfo_y() + 150
        marker_dialog.geometry(f"+{x}+{y}")

        frame = ttk.Frame(marker_dialog, padding=10)
        frame.pack(fill=BOTH, expand=True)

        ttk.Label(frame, text="Marker Style:").grid(row=0, column=0, sticky=W, pady=2)
        marker_names = [name for _, name in self.plot.available_markers]
        current_marker_name = next((name for m, name in self.plot.available_markers if m == self.plot.marker_style),
                                   'circle')
        marker_var = ttk.StringVar(value=current_marker_name)
        ttk.Combobox(frame, textvariable=marker_var, values=marker_names, state='readonly').grid(row=0, column=1,
                                                                                                 sticky=EW, pady=2,
                                                                                                 padx=5)

        ttk.Label(frame, text="Marker Size:").grid(row=1, column=0, sticky=W, pady=2)
        size_var = ttk.DoubleVar(value=self.plot.marker_size)
        ttk.Scale(frame, from_=1, to=20, variable=size_var, orient=HORIZONTAL).grid(row=1, column=1, sticky=EW, pady=2,
                                                                                    padx=5)
        ttk.Entry(frame, textvariable=size_var, width=6).grid(row=1, column=2, sticky=W, pady=2, padx=5)

        ttk.Label(frame, text="Error Cap Size:").grid(row=2, column=0, sticky=W, pady=2)
        capsize_var = ttk.DoubleVar(value=self.plot.errorbar_capsize)
        ttk.Scale(frame, from_=0, to=20, variable=capsize_var, orient=HORIZONTAL).grid(row=2, column=1, sticky=EW,
                                                                                       pady=2, padx=5)
        ttk.Entry(frame, textvariable=capsize_var, width=6).grid(row=2, column=2, sticky=W, pady=2, padx=5)

        ttk.Label(frame, text="Cap Thickness:").grid(row=3, column=0, sticky=W, pady=2)
        capthick_var = ttk.DoubleVar(value=self.plot.errorbar_capthick)
        ttk.Scale(frame, from_=0.0, to=10, variable=capthick_var, orient=HORIZONTAL).grid(row=3, column=1, sticky=EW,
                                                                                          pady=2, padx=5)
        ttk.Entry(frame, textvariable=capthick_var, width=6).grid(row=3, column=2, sticky=W, pady=2, padx=5)

        ttk.Label(frame, text="Error Bar Width:").grid(row=4, column=0, sticky=W, pady=2)
        linewidth_var = ttk.DoubleVar(value=self.plot.errorbar_linewidth)
        ttk.Scale(frame, from_=0.0, to=10, variable=linewidth_var, orient=HORIZONTAL).grid(row=4, column=1, sticky=EW,
                                                                                           pady=2, padx=5)
        ttk.Entry(frame, textvariable=linewidth_var, width=6).grid(row=4, column=2, sticky=W, pady=2, padx=5)

        preview_frame = ttk.LabelFrame(frame, text="Preview", padding=5)
        preview_frame.grid(row=0, column=3, rowspan=5, padx=10, sticky=N + S)

        fig, ax = plt.subplots(figsize=(3, 2), dpi=80)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        x_demo, y_demo, yerr_demo = [0.2, 0.5, 0.8], [0.5, 0.5, 0.5], [0.15, 0.15, 0.15]
        ax.errorbar(x_demo, y_demo, yerr=yerr_demo, fmt='o', color='blue',
                    markersize=size_var.get(), capsize=capsize_var.get(), capthick=capthick_var.get(),
                    elinewidth=linewidth_var.get())
        canvas = FigureCanvasTkAgg(fig, master=preview_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=BOTH, expand=True)

        def update_preview(*_):
            try:
                marker_style = next((m for m, name in self.plot.available_markers if name == marker_var.get()), 'o')
                ax.clear()
                ax.set_xticks([])
                ax.set_yticks([])
                ax.set_xlim(0, 1)
                ax.set_ylim(0, 1)
                ax.errorbar(x_demo, y_demo, yerr=yerr_demo, fmt=marker_style, color='blue',
                            markersize=size_var.get(), capsize=capsize_var.get(), capthick=capthick_var.get(),
                            elinewidth=linewidth_var.get())
                canvas.draw_idle()
            except Exception as e:
                debug_print(f"Preview update error: {e}")

        marker_var.trace_add('write', update_preview)
        size_var.trace_add('write', update_preview)
        capsize_var.trace_add('write', update_preview)
        capthick_var.trace_add('write', update_preview)
        linewidth_var.trace_add('write', update_preview)

        btns = ttk.Frame(frame)
        btns.grid(row=5, column=0, columnspan=4, pady=10)

        def apply_settings():
            self.plot.marker_style = next((m for m, name in self.plot.available_markers if name == marker_var.get()),
                                          'o')
            self.plot.marker_size = size_var.get()
            self.plot.errorbar_capsize = capsize_var.get()
            self.plot.errorbar_capthick = capthick_var.get()
            self.plot.errorbar_linewidth = linewidth_var.get()
            self.plot._layout_signature = None
            self.request_plot_update()

        ttk.Button(btns, text="OK", command=apply_settings).pack(side=LEFT, padx=5)
        ttk.Button(btns, text="Cancel", command=marker_dialog.destroy).pack(side=LEFT, padx=5)
        marker_dialog.bind('<Return>', lambda e: apply_settings())
        marker_dialog.bind('<Escape>', lambda e: marker_dialog.destroy())
        update_preview()

    def confirm_exit(self):
        self.on_close()

    def on_close(self):
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            try:
                plt.close('all')
            except Exception:
                pass
            self.root.destroy()
            self.root.quit()


# ───────────────────────────── Entrypoint ────────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Interactive lightcurve viewer')
    parser.add_argument('--debug', action='store_true', help='Enable debug output')
    args = parser.parse_args()
    DEBUG = args.debug

    root = ttk.Window(themename="flatly")
    root.geometry(f"{WINDOW_WIDTH}x{WINDOW_HEIGHT}")
    app = LightCurveGUI(root)
    root.mainloop()

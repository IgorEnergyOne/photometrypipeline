# -*- coding: utf-8 -*-
"""
Main application window (LightCurveGUI) for pp_interactive_src.
"""
from __future__ import annotations
import os
import re
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple
import numpy as np
import pandas as pd
from astropy.io import fits
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['keymap.save'] = []
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import tkinter as tk
import ttkbootstrap as ttk
from ttkbootstrap.constants import *
from tkinter import filedialog, messagebox, simpledialog
import toolbox
from .constants import (DEFAULT_PERIOD, DEFAULT_PHASE_MAX, OFFSET_ALL_KEY,
                        TIME_STEP, WINDOW_WIDTH, WINDOW_HEIGHT)
from .jpl import iterative_lighttime_correction, calc_reduced_mag
from .models import FitsContext, LightCurveData
from .plot_settings import PlotSettings
from .plot import LightCurvePlot
from .image_viewer import AsteroidImageViewer
from .utils import _safe, next_version, debug_print

class LightCurveGUI:
    """Main application window.

    Owns the Tk root, all data (``lightcurves`` dict), view state, and wires
    together the plot, the asteroid-image viewer, and all menus/toolbars.

    Responsibilities
    ----------------
    * Opening / saving CSV files and ATLAS files.
    * Managing per-lightcurve visibility and vertical offsets.
    * Dispatching keyboard shortcuts and mouse clicks to the plot.
    * Fetching JPL ephemeris data on demand and storing results in the DataFrames.
    """
    def __init__(self, root: ttk.Window):
        self.root = root
        self.root.title("Lightcurve Viewer and Editor")

        # Data model: alias -> LightCurveData
        self.lightcurves: Dict[str, LightCurveData] = {}
        # Per-file visibility / color / offsets
        self.lc_visible: Dict[str, bool] = {}
        # per-file offsets: alias -> { band -> offset }
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

        # Multi-filter state
        self.selected_bands: Set[str] = set()  # empty => all bands
        self.band_vars: Dict[str, tk.BooleanVar] = {}

        # Offset control scope (which band to nudge when ↑/↓)
        self.offset_scope_var = ttk.StringVar(value='auto')  # 'auto', 'ALL', or specific band label

        # Legend position
        self.color_legend_loc_var = ttk.StringVar(value='upper left')

        # Alignment controls
        self.align_reference_var = ttk.StringVar(value='')
        self.calculated_colors: Dict[str, Dict[str, Tuple[float, float]]] = {}  # alias -> { 'R-V': (val, err) }

        # FITS contexts (alias -> FitsContext)
        self.fits_contexts: Dict[str, FitsContext] = {}

        # Step used when nudging vertical offsets with ↑/↓ (mag)
        self.offset_step_var = ttk.DoubleVar(value=0.05)

        # Plot Settings
        self.plot_settings = PlotSettings(self.root, self._refresh_plot_callback, self._get_bands_callback)

        # Plot display toggles (exposed via Plot Settings)
        self.show_legend: bool = True
        self.show_grid: bool = True

        # JPL-derived corrections toggles
        self.use_reduced_mag_var = tk.BooleanVar(value=False)
        self.use_lighttime_var = tk.BooleanVar(value=False)

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
        self.root.bind("<R>", lambda e: self.toggle_show_rejected())
        self.root.bind("<F>", lambda e: self.toggle_flagged_use_filter())
        # Bind Shift+S to show the image for the selected point
        self.root.bind("<S>", self.show_asteroid_image)

    def _refresh_plot_callback(self):
        """Invalidate the plot layout and schedule a redraw.  Passed as a
        callback to PlotSettings so color/style changes propagate immediately."""
        if hasattr(self, 'plot') and getattr(self.plot, 'invalidate_layout', None):
            self.plot.invalidate_layout()
        self.request_plot_update()

    def _get_bands_callback(self) -> List[str]:
        """Return a sorted list of all unique filter bands across loaded lightcurves."""
        return sorted({str(b) for lc in self.lightcurves.values() for b in lc.get_bands()})

    # ——— UI construction ———
    def _build_ui(self):
        """Build the entire widget hierarchy: toolbar, collapsible panels, plot canvas."""
        master = ttk.Frame(self.root)
        master.pack(fill=BOTH, expand=True)
        self.master_frame = master

        # Create a persistent bottom control strip
        self.bottom_control = ttk.Frame(self.master_frame)
        self.bottom_control.pack(side=BOTTOM, fill=X)
        _safe(self.bottom_control.pack_propagate, False)
        _safe(self.bottom_control.configure, height=72)

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
        _safe(self.plot_settings.build_menu, parent=self.plot_settings_menu, pack_btn=False)
        _safe(self.plot_settings_menu.add_cascade, label="Colors", menu=self.plot_settings.color_menu)
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
        self.plot_settings_menu.add_separator()
        self.plot_settings_menu.add_checkbutton(
            label="Use Reduced Magnitudes",
            variable=self.use_reduced_mag_var,
            command=self._on_toggle_reduced_mag,
        )
        self.plot_settings_menu.add_checkbutton(
            label="Use Lighttime Corrections",
            variable=self.use_lighttime_var,
            command=self._on_toggle_lighttime,
        )

        ttk.Button(top, text="Rename Rejected FITS", command=self.rename_rejected_files).pack(side=LEFT, padx=(10, 0))
        ttk.Button(top, text="Help (F1)", command=self.show_help).pack(side=LEFT, padx=(10, 0))
        self.root.bind("<F1>", lambda e: self.show_help())
        _safe(self.root.after, 0, lambda: top._on_configure(None))

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
        ttk.Button(align_frame, text="Reset Offsets", command=self.reset_all_offsets).pack(side=LEFT, padx=5)
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

        # Live chi2 display for current band offset
        ttk.Separator(align_frame, orient=VERTICAL).pack(side=LEFT, fill=Y, padx=5, pady=2)
        self.chi2_label_var = tk.StringVar(value="chi2: --")
        ttk.Label(align_frame, textvariable=self.chi2_label_var).pack(side=LEFT, padx=(5, 0))

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

        ttk.Separator(rot, orient=VERTICAL).pack(side=LEFT, fill=Y, padx=(8, 4), pady=2)
        ttk.Label(rot, text="Phase max:").pack(side=LEFT, padx=(0, 0))
        self.phase_max_var = ttk.StringVar(value=f"{DEFAULT_PHASE_MAX:.2f}")
        phase_max_entry = ttk.Entry(rot, textvariable=self.phase_max_var, width=5)
        phase_max_entry.pack(side=LEFT, padx=3)
        phase_max_entry.bind("<Return>", lambda _e: self._on_phase_max_changed())
        phase_max_entry.bind("<FocusOut>", lambda _e: self._on_phase_max_changed())

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
        _safe(self.canvas.mpl_connect, "pick_event", self.on_pick_label)
        self.plot.parent_gui = self

        # Bottom controls
        bottom = ttk.Frame(self.bottom_control)
        bottom.pack(side=BOTTOM, fill=X)
        _safe(bottom.pack_propagate, False)
        _safe(bottom.configure, height=60)

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

        # Ensure the plot knows the current mode before reading per-mode limits
        self.plot._current_mode = self.mode

        # Prefer the per-mode stored limits over the current axes limits
        stored_limits = self.plot._y_limits_per_mode.get(self.mode)
        if stored_limits is not None:
            cur_min, cur_max = stored_limits
        else:
            try:
                lo, hi = self.plot.ax.get_ylim()
                cur_min = min(lo, hi)
                cur_max = max(lo, hi)
            except Exception:
                cur_min, cur_max = 0.0, 1.0

        cur_nticks = self.plot._y_nticks_per_mode.get(self.mode) or ""
        # Show the user-set label if present, otherwise the current auto label
        cur_label = getattr(self.plot, '_user_ylabel', None) or getattr(self.plot, 'ylabel', 'Magnitude')

        # Create dialog
        dlg = tk.Toplevel(self.root)
        dlg.title("Y-Axis Settings")
        dlg.transient(self.root)
        dlg.resizable(False, False)
        _safe(dlg.geometry, f"+{self.root.winfo_x() + 150}+{self.root.winfo_y() + 150}")

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

                new_label = label_var.get().strip()
                # Store as user override so it persists across redraws in any mode
                self.plot._user_ylabel = new_label if new_label else None
                self.plot.ylabel = new_label or self.plot.ylabel

                # Apply per-mode limits
                self.plot.set_y_limits(min(v_min, v_max), max(v_min, v_max), nt)
                self.plot.invalidate_layout()
                self._refresh_plot_callback()
                dlg.destroy()
            except ValueError:
                messagebox.showerror("Invalid Input", "Please enter valid numbers.", parent=dlg)

        def auto():
            # Clear user label override so the auto-generated label is restored
            self.plot._user_ylabel = None
            # Reset per-mode Y limits for the current mode only
            self.plot.set_y_limits(None, None, None)
            self.plot.invalidate_layout()
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
        """Build the Save drop-down Menubutton and attach it to *parent*."""
        self.save_menu_btn = ttk.Menubutton(parent, text="Save")
        self.save_menu = tk.Menu(self.save_menu_btn, tearoff=0, font=("TkDefaultFont", 10))
        self.save_menu.add_command(label="Save CSV", command=self.save_csv)
        self.save_menu.add_command(label="Save Plot", command=self.save_plot)
        self.save_menu.add_command(label="Save Atlas", command=self.save_atlas)
        self.save_menu_btn["menu"] = self.save_menu
        self.save_menu_btn.pack(side=LEFT, padx=4)


    def show_help(self):
        """Open the keyboard-shortcut / help dialog."""
        msg = (
            "Hotkeys & Controls\n"
            "--------------------\n"
            "q              Quit application\n"
            "r              Toggle rejection for selected point\n"
            "R              Show / Hide rejected points\n"
            "F              Toggle flagged-point colour (filter colour ↔ flagged colour)\n"
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
        """Return the current alias, auto-selecting the sole lightcurve when
        nothing is explicitly selected.  Returns ``None`` if ambiguous."""
        if self.current_lc_alias:
            return self.current_lc_alias
        if len(self.lightcurves) == 1:
            alias = next(iter(self.lightcurves))
            self.current_lc_alias = alias
            return alias
        return None


    def open_csv(self):
        """Prompt for one or more CSV files and load them as lightcurves."""
        files = filedialog.askopenfilenames(initialdir=os.getcwd(), filetypes=[("CSV Files", "*.csv")])
        if not files:
            return
        
        # Single-file mode: clear existing data
        self.lightcurves.clear()
        self.lc_visible.clear()
        self.lc_offsets.clear()
        self.fits_contexts.clear()
        self.selected_bands.clear()

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
        
        # initialize per-file offsets map
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
        if hasattr(self, "_rebuild_color_menu"):
            self._rebuild_color_menu()
        elif hasattr(self, "plot_settings") and hasattr(self.plot_settings, "rebuild_menu"):
            self.plot_settings.rebuild_menu()

        self.request_plot_update()

    def _unique_alias(self, base: str) -> str:
        """Return a unique alias derived from *base*, appending ``_N`` if needed."""
        alias = base
        k = 1
        while alias in self.lightcurves:
            stem, ext = os.path.splitext(base)
            alias = f"{stem}_{k}{ext}"
            k += 1
        return alias

    def _fetch_jpl_data(self, lc, alias):
        """Fetch basic JPL ephemeris columns (r, delta, alpha, EclLon/Lat) for *lc*
        and append them to ``lc.df``.  Returns ``(True, [added_cols])`` or ``False``."""
        # 1. Get / ask for target + observatory
        ctx = self.fits_contexts.get(alias)
        if ctx:
            target   = ctx.target_object   or ""
            location = ctx.observatory_code or ""
        else:
            target   = ""
            location = ""

        res = self._ask_target_and_location(alias,
                                            prefill_target=target,
                                            prefill_location=location)
        if res is None:
            return False  # user cancelled
        target, location, fits_path = res

        # If the user loaded a FITS file, fits_contexts is already updated
        # inside _ask_target_and_location._load_from_fits.
        # For manual entry build a minimal context so the values persist.
        if alias not in self.fits_contexts:
            ctx = FitsContext.__new__(FitsContext)
            ctx.filepath = fits_path or ""
            ctx.header = None
            ctx.obsparam = None
            ctx.target_object = target
            ctx.observatory_code = location
            self.fits_contexts[alias] = ctx

        # 2. Epochs — the CSV uses 'julian_date'; fall back to 'jd' if needed
        if 'julian_date' in lc.df.columns:
            epochs = lc.df['julian_date'].values
        elif 'jd' in lc.df.columns:
            epochs = lc.df['jd'].values
        else:
            raise ValueError("CSV is missing a Julian Date column ('julian_date' or 'jd').")

        # 3. Query
        self.root.config(cursor="watch")
        self.root.update()
        try:
            eph_df = toolbox.jpl_query_eph(target, epochs, location)

            col_map = {
                'r': 'r',
                'delta': 'delta',
                'alpha': 'alpha_true',
                'EclLon': 'ObsEclLon',
                'EclLat': 'ObsEclLat'
            }

            added = []
            for src, dst in col_map.items():
                if src in eph_df.columns:
                    if len(eph_df) == len(lc.df):
                        lc.df[dst] = eph_df[src].values
                        added.append(dst)
                else:
                    print(f"JPL warning: Column {src} not returned. (Got: {eph_df.columns})")

            return True, added

        except Exception as e:
            raise e
        finally:
            self.root.config(cursor="")

    def save_csv(self):
        if not self.current_lc_alias:
            messagebox.showerror("Nothing to save", "Select a lightcurve entry to save.")
            return
        lc = self.lightcurves.get(self.current_lc_alias)
        if lc is None:
            messagebox.showerror("Nothing to save", "No lightcurve selected.")
            return

        JPL_COLS = ['r', 'delta', 'alpha_true', 'ObsEclLon', 'ObsEclLat']

        # Step 1 – determine what's already in the dataframe
        present = [c for c in JPL_COLS if c in lc.df.columns]
        missing = [c for c in JPL_COLS if c not in lc.df.columns]

        if not missing:
            # All JPL columns are already present — nothing to do, they will
            # be written to the file automatically.
            pass
        else:
            # Some (or all) columns are absent — ask the user whether to fetch.
            if present:
                detail = (f"Present : {', '.join(present)}\n"
                          f"Missing : {', '.join(missing)}\n\n")
            else:
                detail = f"Missing : {', '.join(missing)}\n\n"

            fetch_ans = messagebox.askyesno(
                "Add JPL Data to CSV?",
                f"The following JPL ephemeris columns are not in the current data:\n\n"
                f"{detail}"
                "Fetch missing columns from JPL Horizons now?\n"
                "(Requires internet connection and target/observatory info)\n\n"
                "Choose 'No' to save the CSV without JPL columns.",
            )

            if fetch_ans:
                try:
                    success, added_cols = self._fetch_jpl_data(lc, self.current_lc_alias)
                    if success and added_cols:
                        messagebox.showinfo(
                            "JPL Data Fetched",
                            f"Successfully added column(s):\n  {', '.join(added_cols)}",
                        )
                    elif success:
                        messagebox.showwarning("JPL Data", "Query succeeded but no new columns were added.")
                except Exception as e:
                    messagebox.showerror("JPL Error", f"Failed to fetch JPL data:\n{e}")
                    # user can still choose to save what we have

        # Step 2 – choose save path
        file = filedialog.asksaveasfilename(
            initialdir=os.path.dirname(lc.filename) if lc.filename else os.getcwd(),
            initialfile=os.path.basename(lc.filename) if lc.filename else "lightcurve.csv",
            defaultextension=".csv",
            filetypes=[("CSV Files", "*.csv")],
        )
        if not file:
            return  # user cancelled the save dialog

        # Step 3 – save; whatever columns are in lc.df (including any newly
        # fetched JPL ones) are written as-is.
        lc.filename = file
        lc.df.to_csv(file, index=False)
        messagebox.showinfo("Saved", f"CSV saved as {file}")


    def save_plot(self):
        """Save the current plot canvas as a PNG file."""
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

    def _show_atlas_options_dialog(self, lc) -> Optional[dict]:
        """
        Show a dialog asking which corrections to apply when exporting the
        ATLAS file.  Returns a dict with keys:
          'use_lighttime'   bool
          'use_reduced_mag' bool
        or None if the user cancelled.

        Both checkboxes are always freely toggleable.  No JPL data is fetched
        here - if the data is missing it will be fetched later in save_atlas
        once the FITS file path is known, as that path is required anyway.
        """
        has_lt = ('corrected_jd' in lc.df.columns)
        has_rm = ('reduced_mag'  in lc.df.columns)

        dlg = tk.Toplevel(self.root)
        dlg.title("Atlas Export Options")
        dlg.transient(self.root)
        dlg.resizable(False, False)
        dlg.grab_set()
        _safe(dlg.geometry, f"+{self.root.winfo_x()+200}+{self.root.winfo_y()+200}")

        frm = ttk.Frame(dlg, padding=14)
        frm.pack(fill=BOTH, expand=True)

        ttk.Label(
            frm,
            text="Select corrections to apply to the exported ATLAS file:",
            wraplength=380,
        ).pack(anchor='w', pady=(0, 8))

        lt_var = tk.BooleanVar(value=has_lt)
        rm_var = tk.BooleanVar(value=has_rm)

        ttk.Checkbutton(
            frm,
            text="Use Lighttime-corrected JD  (replaces julian_date with corrected_jd)",
            variable=lt_var,
        ).pack(anchor='w', pady=2)

        ttk.Checkbutton(
            frm,
            text="Use Reduced Magnitudes  (replaces mag with H(a) = mag - 5*log10(r*delta))",
            variable=rm_var,
        ).pack(anchor='w', pady=2)

        # Informational note about data availability
        notes = []
        if not has_lt:
            notes.append("  - Lighttime correction: JPL data not yet fetched (will be fetched if needed)")
        if not has_rm:
            notes.append("  - Reduced magnitudes: JPL data not yet fetched (will be fetched if needed)")
        if notes:
            ttk.Label(
                frm,
                text="Data status:\n" + "\n".join(notes),
                foreground='gray',
                wraplength=380,
            ).pack(anchor='w', pady=(8, 0))

        result: dict = {}

        def _ok():
            result['use_lighttime']   = bool(lt_var.get())
            result['use_reduced_mag'] = bool(rm_var.get())
            dlg.destroy()

        def _cancel():
            dlg.destroy()

        btn_frm = ttk.Frame(frm)
        btn_frm.pack(pady=(12, 0))
        ttk.Button(btn_frm, text="OK",     command=_ok,     style="Accent.TButton").pack(side=LEFT, padx=6)
        ttk.Button(btn_frm, text="Cancel", command=_cancel).pack(side=LEFT, padx=6)

        dlg.bind("<Return>", lambda e: _ok())
        dlg.bind("<Escape>", lambda e: _cancel())
        dlg.wait_window()

        return result if result else None

    def _build_atlas_df(self, lc, use_lighttime: bool, use_reduced_mag: bool) -> pd.DataFrame:
        """Return a copy of ``lc.df`` with column substitutions applied.

        * ``use_lighttime``   – replaces ``julian_date`` with ``corrected_jd``
        * ``use_reduced_mag`` – replaces ``mag`` with ``reduced_mag``

        The original DataFrame is never modified.
        """
        df = lc.df.copy()

        if use_lighttime and 'corrected_jd' in df.columns:
            jd_col = 'julian_date' if 'julian_date' in df.columns else 'jd'
            df[jd_col] = df['corrected_jd']

        if use_reduced_mag and 'reduced_mag' in df.columns:
            df['mag'] = df['reduced_mag']

        return df

    def save_atlas(self):
        if not self.current_lc_alias:
            messagebox.showerror("Nothing to save", "Select a lightcurve entry to save atlas for.")
            return
        lc = self.lightcurves.get(self.current_lc_alias)
        if lc is None or lc.df is None:
            messagebox.showerror("Nothing to save", "Selected lightcurve has no data.")
            return

        # Ask the user which corrections to apply
        opts = self._show_atlas_options_dialog(lc)
        if opts is None:          # user cancelled
            return
        use_lighttime   = opts['use_lighttime']
        use_reduced_mag = opts['use_reduced_mag']

        # Resolve FITS file (needed by pp_atlas and for JPL queries)
        alias = self.current_lc_alias
        ctx = self.fits_contexts.get(alias)
        prefill_target   = (ctx.target_object    or "") if ctx else ""
        prefill_location = (ctx.observatory_code or "") if ctx else ""
        fits_filepath    = ctx.filepath if ctx else ""

        # If we do not yet have a context (or the stored filepath is empty),
        # ask the user — they can either type the IDs manually or load a FITS file.
        if not ctx or not fits_filepath:
            res = self._ask_target_and_location(alias,
                                                prefill_target=prefill_target,
                                                prefill_location=prefill_location)
            if res is None:
                return   # cancelled
            target, location, loaded_fits = res
            fits_filepath = loaded_fits or ""

            if alias not in self.fits_contexts:
                # build minimal context from manual input
                new_ctx = FitsContext.__new__(FitsContext)
                new_ctx.filepath = fits_filepath
                new_ctx.header = None
                new_ctx.obsparam = None
                new_ctx.target_object = target
                new_ctx.observatory_code = location
                self.fits_contexts[alias] = new_ctx

        # pp_atlas requires a real FITS file for its header section.
        # If the user entered IDs manually (no FITS file), ask for the file now.
        if not fits_filepath or not os.path.isfile(fits_filepath):
            filetypes = [("FITS files", "*.fits *.fts"), ("All files", "*.*")]
            fits_filepath = filedialog.askopenfilename(
                title="Select the FITS file for the ATLAS header", filetypes=filetypes
            )
            if not fits_filepath:
                messagebox.showerror("Nothing selected",
                                     "A FITS file is required to build the ATLAS header.")
                return
            # update stored context with the real path
            ctx2 = self.fits_contexts.get(alias)
            if ctx2:
                ctx2.filepath = fits_filepath


        # Fetch JPL data if a correction was requested but data is missing
        needs_jpl = (use_lighttime and 'corrected_jd' not in lc.df.columns) or \
                    (use_reduced_mag and 'reduced_mag'  not in lc.df.columns)
        if needs_jpl:
            ok = self._fetch_jpl_data_for_corrections(alias)
            if not ok:
                # User aborted the fetch; ask whether to continue without corrections
                proceed = messagebox.askyesno(
                    "JPL Fetch Failed",
                    "JPL data could not be fetched.  The requested corrections will be\n"
                    "skipped.  Continue saving the ATLAS file without corrections?",
                )
                if not proceed:
                    return
                use_lighttime   = use_lighttime   and ('corrected_jd' in lc.df.columns)
                use_reduced_mag = use_reduced_mag and ('reduced_mag'  in lc.df.columns)

        # Output file path
        atlas_name = f"{Path(lc.filename).stem}.ATL" if lc.filename else "output.ATL"
        file = filedialog.asksaveasfilename(
            initialfile=atlas_name,
            defaultextension=".ATL",
            filetypes=[("Atlas Files", "*.ATL")],
        )
        if not file:
            messagebox.showerror("Nothing selected", "No atlas filename provided")
            return

        # ── 5. Build working DataFrame with requested substitutions ──────────
        work_df = self._build_atlas_df(lc, use_lighttime, use_reduced_mag)

        # ── 6. Write temp CSV(s) and call pp_atlas ───────────────────────────
        # Build the correction flags that will be forwarded to pp_atlas so that
        # the REDUCED MAG. and LT CORRECTED header fields are set to T/F correctly.
        correction_flags = ""
        if use_reduced_mag:
            correction_flags += " -use_reduced_mag"
        if use_lighttime:
            correction_flags += " -use_lt_corrected"

        filters = lc.get_bands()
        tmps = []
        tmp_atls = []
        if len(filters) > 1:
            for f in filters:
                df_filt = work_df[work_df['band'] == f]
                if df_filt.empty:
                    continue
                tmpf = f"tmp_atlas_{f}.csv"
                tmp_atl = f"{tmpf[:-4]}.ATL"
                tmp_atls.append(tmp_atl)
                df_filt.to_csv(tmpf, header=True, index=False)
                tmps.append(tmpf)
                atlas_cmd = (f"pp_atlas -fname_header {fits_filepath} "
                             f"-fname_photo {tmpf} -fname_out {tmpf[:-4]}.ATL"
                             f"{correction_flags}")
                subprocess.call(['/bin/sh', '-i', '-c', atlas_cmd])
            atlas_cmd = (f"pp_atlas -combine "
                         + " ".join([f"{tmp_file[:-4]}.ATL" for tmp_file in tmps])
                         + f" -fname_out {file}")
            try:
                subprocess.call(['/bin/sh', '-i', '-c', atlas_cmd])
            except Exception as e:
                messagebox.showerror("Atlas save failed", f"Error calling atlas command: {e}")
            finally:
                for tmp_atl, tmp_csv in zip(tmp_atls, tmps):
                    if os.path.exists(tmp_atl):
                        os.remove(tmp_atl)
                    if os.path.exists(tmp_csv):
                        os.remove(tmp_csv)
        else:
            tmp = 'tmp_atlas.csv'
            tmps.append(tmp)
            work_df.to_csv(tmp, header=True, index=False)
            atlas_cmd = (f"pp_atlas -fname_header {fits_filepath} "
                         f"-fname_photo {tmp} -fname_out {file}"
                         f"{correction_flags}")
            try:
                subprocess.call(['/bin/sh', '-i', '-c', atlas_cmd])
            except Exception as e:
                messagebox.showerror("Atlas save failed", f"Error calling atlas command: {e}")
            finally:
                for tmp in tmps:
                    if os.path.exists(tmp):
                        os.remove(tmp)

        # ── 7. Summary ───────────────────────────────────────────────────────
        notes = []
        if use_lighttime:
            notes.append("lighttime-corrected JD")
        if use_reduced_mag:
            notes.append("reduced magnitudes H(a)")
        extra = f"\nApplied: {', '.join(notes)}" if notes else ""
        messagebox.showinfo("Saved", f"ATLAS file saved as {file}{extra}")

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
        rejected_all = lc.arr.get('rejected')
        if rejected_all is None:
            rejected_all = np.zeros(len(jd_all), dtype=bool)

        # Exclude rejected points from reference band
        ref_mask = (band_all == ref_band) & (~rejected_all)
        if not np.any(ref_mask):
            messagebox.showwarning("Error", f"Reference band {ref_band} has no valid (non-rejected) data.")
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
            # Exclude rejected points from this band
            b_mask = (band_all == b) & (~rejected_all)
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
            n_pts = int(np.sum(valid))
            weights = 1.0 / v
            
            w_mean = np.sum(d * weights) / np.sum(weights)
            w_err = np.sqrt(1.0 / np.sum(weights))
            
            # Compute reduced chi2 as alignment quality metric
            residuals = d - w_mean
            chi2 = np.sum((residuals ** 2) / v)
            dof = max(n_pts - 1, 1)
            chi2_red = chi2 / dof

            # Set offset
            # We want B + Offset = Ref => Offset as computed.
            # But wait, logic in plotting is: Y_plot = Ybase + off
            # So if we set off = w_mean, Y_plot = B + (Ref - B) = Ref. Correct.
            
            # Store offset
            self.lc_offsets[alias][b] = round(w_mean, 4)
            
            # Store Color Info with metric
            # Tuple: (value, error, reduced_chi2, n_points)
            color_name = f"{ref_band}-{b}"
            self.calculated_colors[alias][color_name] = (w_mean, w_err, chi2_red, n_pts)

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

    def reset_all_offsets(self):
        """Reset every vertical offset for every loaded lightcurve back to zero.

        Also clears any computed colour indices, since those were derived from
        the (now-removed) alignment offsets.
        """
        if not self.lightcurves:
            return

        for alias, lc in self.lightcurves.items():
            bands = lc.get_bands()
            self.lc_offsets[alias] = {OFFSET_ALL_KEY: 0.0}
            for b in bands:
                self.lc_offsets[alias][b] = 0.0

        # Clear colour indices that depended on the old alignment
        self.calculated_colors.clear()

        # Reset chi2 display
        self.chi2_label_var.set("chi2: --")

        # Sync the offset slider/entry to the new zero state
        self._sync_offset_var_to_scope()
        self._refresh_offset_scope_choices()

        if hasattr(self, 'plot') and hasattr(self.plot, 'invalidate_layout'):
            self.plot.invalidate_layout()
        self.request_plot_update()

    def show_colors_dialog(self):
        """Open a Treeview dialog showing computed colour indices for the active alias."""
        alias = self.current_lc_alias
        if not alias or alias not in self.calculated_colors or not self.calculated_colors[alias]:
            messagebox.showinfo("Colors", "No calculated colors available.\nRun 'Auto-Align' first.")
            return
            
        colors = self.calculated_colors[alias]
        
        dlg = tk.Toplevel(self.root)
        dlg.title(f"Color Indices: {alias}")
        dlg.transient(self.root)
        _safe(dlg.geometry, f"+{self.root.winfo_x()+100}+{self.root.winfo_y()+100}")

        # Table
        cols = ("Color", "Value", "Error", "chi2_red", "N")
        tree = ttk.Treeview(dlg, columns=cols, show='headings', height=len(colors)+2)
        for c in cols:
            tree.heading(c, text=c)
            w = 80 if c in ("chi2_red", "N") else 100
            tree.column(c, width=w, anchor='center')
        tree.pack(fill=BOTH, expand=True, padx=10, pady=10)
        
        for name, cdata in colors.items():
            val, err = cdata[0], cdata[1]
            chi2_red = cdata[2] if len(cdata) > 2 else float('nan')
            n_pts = cdata[3] if len(cdata) > 3 else 0
            chi2_str = f"{chi2_red:.3f}" if np.isfinite(chi2_red) else "--"
            n_str = str(n_pts) if n_pts > 0 else "--"
            tree.insert("", END, values=(name, f"{val:.4f}", f"{err:.4f}", chi2_str, n_str))

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
                        f.write(f"Color | Value | Error | chi2_red | N\n")
                        for name, cdata in colors.items():
                            val, err = cdata[0], cdata[1]
                            chi2_red = cdata[2] if len(cdata) > 2 else float('nan')
                            n_pts = cdata[3] if len(cdata) > 3 else 0
                            chi2_str = f"{chi2_red:.6f}" if np.isfinite(chi2_red) else "N/A"
                            f.write(f"{name} | {val:.6f} | {err:.6f} | {chi2_str} | {n_pts}\n")
                    messagebox.showinfo("Saved", f"Colors saved to {fname}")
                except Exception as e:
                    messagebox.showerror("Error", f"Failed to save: {e}")

        ttk.Button(btn_frame, text="Save to File", command=save_colors).pack(side=LEFT, padx=20)
        ttk.Button(btn_frame, text="Close", command=dlg.destroy).pack(side=RIGHT, padx=20)

    def toggle_align_frame(self):
        """Toggle the Alignment & Colors collapsible panel open/closed."""
        is_vis = self.align_visible.get()
        new_vis = not is_vis
        self.align_visible.set(new_vis)
        
        if new_vis:
            self.toggle_btn.config(text="[ - ] Alignment & Colors")
            self.align_inner.pack(side=TOP, fill=X)
        else:
            self.toggle_btn.config(text="[ + ] Alignment & Colors")
            self.align_inner.pack_forget()

    # ——— FITS Context Helpers ———

    def _ask_target_and_location(
        self,
        alias: str,
        prefill_target: str = "",
        prefill_location: str = "",
    ):
        """
        Show a dialog with two text fields (Asteroid ID, Observatory Code) and a
        'Load from FITS file' button that auto-fills them.

        Returns a tuple (target, location, fits_filepath_or_None) on OK,
        or None if the user cancelled.

        fits_filepath_or_None is the path to the FITS file if the user loaded one,
        otherwise None (meaning the values were entered manually).
        """
        # Pre-populate from existing FitsContext if available and no explicit prefill
        ctx = self.fits_contexts.get(alias)
        if ctx:
            if not prefill_target:
                prefill_target   = ctx.target_object   or ""
            if not prefill_location:
                prefill_location = ctx.observatory_code or ""

        dlg = tk.Toplevel(self.root)
        dlg.title("Target & Observatory")
        dlg.transient(self.root)
        dlg.resizable(False, False)
        dlg.grab_set()
        _safe(dlg.geometry, f"+{self.root.winfo_x()+220}+{self.root.winfo_y()+180}")

        frm = ttk.Frame(dlg, padding=14)
        frm.pack(fill=tk.BOTH, expand=True)

        ttk.Label(frm, text="Provide target and observatory information:",
                  wraplength=360).grid(row=0, column=0, columnspan=3, sticky='w', pady=(0, 10))

        ttk.Label(frm, text="Asteroid ID:").grid(row=1, column=0, sticky='e', padx=(0, 6), pady=4)
        target_var = tk.StringVar(value=prefill_target)
        target_entry = ttk.Entry(frm, textvariable=target_var, width=28)
        target_entry.grid(row=1, column=1, columnspan=2, sticky='ew', pady=4)

        ttk.Label(frm, text="Observatory Code:").grid(row=2, column=0, sticky='e', padx=(0, 6), pady=4)
        loc_var = tk.StringVar(value=prefill_location)
        loc_entry = ttk.Entry(frm, textvariable=loc_var, width=28)
        loc_entry.grid(row=2, column=1, columnspan=2, sticky='ew', pady=4)

        ttk.Separator(frm, orient='horizontal').grid(row=3, column=0, columnspan=3, sticky='ew', pady=8)

        # --- status label shown after a FITS file is loaded ---
        fits_path_holder = [None]   # mutable cell
        fits_status_var = tk.StringVar(value="No FITS file loaded")
        ttk.Label(frm, textvariable=fits_status_var, foreground='gray',
                  wraplength=360).grid(row=4, column=0, columnspan=3, sticky='w')

        def _load_from_fits():
            lc = self.lightcurves.get(alias)
            init_dir = os.path.dirname(lc.filename) if lc and lc.filename else os.getcwd()
            ftypes = [("FITS files", "*.fits *.fts"), ("All files", "*.*")]
            fpath = filedialog.askopenfilename(
                parent=dlg,
                initialdir=init_dir,
                title="Select FITS header file",
                filetypes=ftypes,
            )
            if not fpath:
                return
            try:
                new_ctx = FitsContext(fpath)
            except Exception as e:
                messagebox.showerror("FITS Error", f"Could not read FITS file:\n{e}", parent=dlg)
                return
            # populate the fields
            if new_ctx.target_object and new_ctx.target_object != "Unknown":
                target_var.set(new_ctx.target_object)
            if new_ctx.observatory_code:
                loc_var.set(new_ctx.observatory_code)
            fits_path_holder[0] = fpath
            fits_status_var.set(f"Loaded: {os.path.basename(fpath)}")
            # store in fits_contexts so other operations can reuse it
            self.fits_contexts[alias] = new_ctx

        ttk.Button(frm, text="Load from FITS file", command=_load_from_fits).grid(
            row=5, column=0, columnspan=3, pady=(6, 2))

        result = [None]

        def _ok():
            t = target_var.get().strip()
            loc = loc_var.get().strip()
            if not t:
                messagebox.showwarning("Missing value", "Please enter the Asteroid ID.", parent=dlg)
                return
            if not loc:
                messagebox.showwarning("Missing value", "Please enter the Observatory Code.", parent=dlg)
                return
            result[0] = (t, loc, fits_path_holder[0])
            dlg.destroy()

        def _cancel():
            dlg.destroy()

        btn_frm = ttk.Frame(frm)
        btn_frm.grid(row=6, column=0, columnspan=3, pady=(10, 0))
        ttk.Button(btn_frm, text="OK",     command=_ok,     style="Accent.TButton").pack(side=tk.LEFT, padx=6)
        ttk.Button(btn_frm, text="Cancel", command=_cancel).pack(side=tk.LEFT, padx=6)

        frm.columnconfigure(1, weight=1)
        target_entry.focus_set()
        dlg.bind("<Return>", lambda e: _ok())
        dlg.bind("<Escape>", lambda e: _cancel())
        dlg.wait_window()

        return result[0]   # (target, location, fits_path) or None

    def ensure_fits_context(self) -> bool:
        """
        Checks if a FITS context is available for the current lightcurve.
        If not, shows the target/observatory dialog (with optional FITS load).
        Returns True if a context is available afterwards, False if cancelled.
        """
        alias = self.current_lc_alias
        if not alias:
            return False

        if alias in self.fits_contexts:
            return True

        res = self._ask_target_and_location(alias)
        if res is None:
            return False
        target, location, fits_path = res

        if alias not in self.fits_contexts:
            # Build a minimal FitsContext-like object from manual input
            ctx = FitsContext.__new__(FitsContext)
            ctx.filepath = fits_path or ""
            ctx.header = None
            ctx.obsparam = None
            ctx.target_object = target
            ctx.observatory_code = location
            self.fits_contexts[alias] = ctx

        return True


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
        """Handle slider/entry changes to the vertical offset for the current scope."""
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
        """Nudge the active band/alias offset by *delta* magnitudes and redraw."""
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
        """Update the calculated color value for a band if the user manually changes its offset.
        Also recomputes chi2_red from the actual data."""
        alias = self.current_lc_alias
        if not alias:
            self._update_chi2_label(None)
            return

        # Skip if band is ALL or reference band itself
        ref_band = self.align_reference_var.get() if hasattr(self, 'align_reference_var') else ''
        if not band or band == OFFSET_ALL_KEY or band == ref_band:
            self._update_chi2_label(None)
            return

        # Always compute chi2 from actual data, regardless of whether auto-align was run
        chi2_red, n_pts = self._compute_chi2_for_offset(alias, band, offset_val)

        # Update calculated_colors if entry exists
        if alias in self.calculated_colors:
            suffix = f"-{band}"
            for color_name in self.calculated_colors[alias]:
                if color_name.endswith(suffix):
                    cdata = self.calculated_colors[alias][color_name]
                    err = cdata[1] if len(cdata) > 1 else 0.0
                    self.calculated_colors[alias][color_name] = (offset_val, err, chi2_red, n_pts)
                    break

        # Update the live chi2 label
        self._update_chi2_label(chi2_red, n_pts)

    def _compute_chi2_for_offset(self, alias, band, offset_val):
        """Compute reduced chi2 for a given band at a given offset relative to the reference.
        Returns (chi2_red, n_points)."""
        lc = self.lightcurves.get(alias)
        if lc is None or lc.df is None or lc.df.empty:
            return float('nan'), 0

        ref_band = self.align_reference_var.get()
        if not ref_band:
            return float('nan'), 0

        jd_all = lc.arr['jd']
        mag_all = lc.arr['mag']
        sig_all = lc.arr['sig']
        band_all = lc.arr['band'].astype(str)
        rejected_all = lc.arr.get('rejected')
        if rejected_all is None:
            rejected_all = np.zeros(len(jd_all), dtype=bool)

        # Reference band data (excluding rejected)
        ref_mask = (band_all == ref_band) & (~rejected_all)
        if not np.any(ref_mask):
            return float('nan'), 0

        ref_idx = np.argsort(jd_all[ref_mask])
        ref_t = jd_all[ref_mask][ref_idx]
        ref_m = mag_all[ref_mask][ref_idx]
        ref_e = sig_all[ref_mask][ref_idx]

        # Target band data (excluding rejected)
        b_mask = (band_all == band) & (~rejected_all)
        if not np.any(b_mask):
            return float('nan'), 0

        b_t = jd_all[b_mask]
        b_m = mag_all[b_mask]
        b_e = sig_all[b_mask]

        # Interpolate reference to target times
        ref_m_interp = np.interp(b_t, ref_t, ref_m, left=np.nan, right=np.nan)
        ref_var = ref_e ** 2
        ref_var_interp = np.interp(b_t, ref_t, ref_var, left=np.nan, right=np.nan)

        # Residuals: (b_m + offset) should match ref_m_interp
        # offset_val = ref - b, so ref_m_interp - (b_m + 0) should equal offset
        # residual = (ref_m_interp - b_m) - offset_val
        diff = ref_m_interp - b_m
        var_diff = b_e ** 2 + ref_var_interp

        valid = np.isfinite(diff) & np.isfinite(var_diff)
        n_pts = int(np.sum(valid))
        if n_pts < 2:
            return float('nan'), n_pts

        d = diff[valid]
        v = var_diff[valid]

        residuals = d - offset_val
        chi2 = np.sum((residuals ** 2) / v)
        dof = max(n_pts - 1, 1)
        chi2_red = chi2 / dof

        return chi2_red, n_pts

    def _update_chi2_label(self, chi2_red=None, n_pts=None):
        """Refresh the chi2 status label in the alignment panel."""
        if not hasattr(self, 'chi2_label_var'):
            return
        if chi2_red is None or not np.isfinite(chi2_red):
            self.chi2_label_var.set(f"chi2: --")
        else:
            n_str = f" (N={n_pts})" if n_pts is not None else ""
            self.chi2_label_var.set(f"chi2: {chi2_red:.3f}{n_str}")

    def _refresh_chi2_for_current_band(self, band_key=None):
        """Recompute and display chi-squared for the current band and offset."""
        alias = self.current_lc_alias
        if not alias:
            self._update_chi2_label(None)
            return
        if band_key is None:
            band_key = self._active_offset_band()
        if band_key == OFFSET_ALL_KEY:
            self._update_chi2_label(None)
            return
        ref_band = self.align_reference_var.get() if hasattr(self, 'align_reference_var') else ''
        if not ref_band or band_key == ref_band:
            self._update_chi2_label(None)
            return
        # Try stored value first
        if alias in self.calculated_colors:
            suffix = f"-{band_key}"
            for color_name, cdata in self.calculated_colors[alias].items():
                if color_name.endswith(suffix) and len(cdata) > 2:
                    self._update_chi2_label(cdata[2], cdata[3] if len(cdata) > 3 else None)
                    return
        # Otherwise compute from current offset
        off = self.lc_offsets.get(alias, {}).get(band_key, 0.0)
        chi2_red, n_pts = self._compute_chi2_for_offset(alias, band_key, float(off))
        self._update_chi2_label(chi2_red, n_pts)


    # ——— Offsets ———
    def _refresh_offset_scope_choices(self):
        """Update the offset scope combobox values to match the loaded bands."""
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
        """Read and validate the offset step entry; returns 0.05 on error."""
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
            self._update_chi2_label(None)
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
        # Update chi2 label for the current band
        self._refresh_chi2_for_current_band(band_key)

    # ——— State setters ——
    def set_mode(self):
        """Apply the photometry mode selected in the Mode combobox."""
        self.mode = self.mode_var.get()
        if hasattr(self.plot, "invalidate_layout"):
            self.plot.invalidate_layout()
        self.request_plot_update()
        self.plot.blit.quick_redraw()

    def set_errorbar_type(self):
        """Apply the error-bar type selected in the Error Bars combobox."""
        self.errorbar_type = self.errorbar_var.get()
        if hasattr(self.plot, "invalidate_layout"):
            self.plot.invalidate_layout()
        self.request_plot_update()
        self.plot.blit.quick_redraw()

    def set_time_mode(self):
        """Apply the time-axis mode selected in the Time Axis combobox."""
        self.time_mode = self.time_var.get()
        # Reset any user-customised X label so the new mode's auto-label takes over
        if hasattr(self, 'plot') and hasattr(self.plot, '_user_xlabel'):
            self.plot._user_xlabel = None
        if hasattr(self.plot, "invalidate_layout"):
            self.plot.invalidate_layout()
        self.request_plot_update()
        self.plot.blit.quick_redraw()

    def set_show_rejected(self):
        """Toggle display of rejected points from the Plot Settings checkbutton."""
        self.show_rejected = bool(self.toggle_rejected_var.get())
        self.request_plot_update()

    def toggle_show_rejected(self):
        """Hotkey R: show/hide rejected points."""
        self.show_rejected = not self.show_rejected
        self.toggle_rejected_var.set(self.show_rejected)
        self.plot.invalidate_layout()
        self.request_plot_update()

    def toggle_flagged_use_filter(self):
        """Hotkey F: toggle whether flagged points use the filter colour or the
        dedicated flagged colour."""
        if hasattr(self, 'plot_settings'):
            self.plot_settings.flagged_use_filter = not self.plot_settings.flagged_use_filter
            self.plot_settings.rebuild_menu()
            self.plot.invalidate_layout()
            self.request_plot_update()

    def _on_toggle_reduced_mag(self):
        """Handle toggling of 'Use Reduced Magnitudes' in the Plot Settings menu."""
        if self.use_reduced_mag_var.get():
            if not self._ensure_jpl_data_for_corrections("Reduced Magnitudes"):
                self.use_reduced_mag_var.set(False)
                return
        if hasattr(self.plot, "invalidate_layout"):
            self.plot.invalidate_layout()
        self.request_plot_update()

    def _on_toggle_lighttime(self):
        """Handle toggling of 'Use Lighttime Corrections' in the Plot Settings menu."""
        if self.use_lighttime_var.get():
            if not self._ensure_jpl_data_for_corrections("Lighttime Corrections"):
                self.use_lighttime_var.set(False)
                return
        if hasattr(self.plot, "invalidate_layout"):
            self.plot.invalidate_layout()
        self.request_plot_update()

    def _ensure_jpl_data_for_corrections(self, feature_name: str) -> bool:
        """
        Check if JPL ephemeris data (r, delta, corrected_jd) is available for the
        current lightcurve.  If not, ask the user and fetch it.
        Returns True when data is available (now or already), False otherwise.
        """
        alias = self.current_lc_alias
        if not alias:
            messagebox.showwarning("No Data", "No lightcurve loaded/selected.")
            return False

        lc = self.lightcurves.get(alias)
        if lc is None or lc.df is None:
            return False

        has_data = ('r' in lc.df.columns and 'delta' in lc.df.columns)
        if has_data:
            return True

        reply = messagebox.askyesno(
            "JPL Data Missing",
            f"Calculating {feature_name} requires JPL Horizons data (r, delta, ...) which is not yet available.\n"
            "Fetch it now from JPL Horizons?  (Requires internet connection)",
        )
        if reply:
            return self._fetch_jpl_data_for_corrections(alias)
        return False

    def _fetch_jpl_data_for_corrections(self, alias: str) -> bool:
        """Fetch iterative-lighttime-corrected ephemeris data for *alias* and attach
        the result columns (``corrected_jd``, ``reduced_mag``, etc.) to ``lc.df``.
        Returns ``True`` on success."""
        lc = self.lightcurves.get(alias)
        if lc is None:
            return False

        # ── 1. Target + observatory ──────────────────────────────────────────
        ctx = self.fits_contexts.get(alias)
        prefill_target   = (ctx.target_object    or "") if ctx else ""
        prefill_location = (ctx.observatory_code or "") if ctx else ""

        res = self._ask_target_and_location(alias,
                                            prefill_target=prefill_target,
                                            prefill_location=prefill_location)
        if res is None:
            return False
        target, location, _fits_path = res

        # If user entered manually (no FITS loaded), build a minimal context
        if alias not in self.fits_contexts:
            ctx = FitsContext.__new__(FitsContext)
            ctx.filepath = _fits_path or ""
            ctx.header = None
            ctx.obsparam = None
            ctx.target_object = target
            ctx.observatory_code = location
            self.fits_contexts[alias] = ctx


        # ── 2. Epochs ───────────────────────────────────────────────────────
        jd_col = 'julian_date' if 'julian_date' in lc.df.columns else 'jd' if 'jd' in lc.df.columns else None
        if jd_col is None:
            messagebox.showerror("Missing Column", "The CSV is missing a 'julian_date' or 'jd' column.")
            return False
        epochs = lc.df[jd_col].dropna().tolist()
        if not epochs:
            messagebox.showerror("No Data", "No valid JD values found in the lightcurve.")
            return False

        # ── 3. Progress dialog ───────────────────────────────────────────────
        prog_dlg = tk.Toplevel(self.root)
        prog_dlg.title("Fetching JPL Data...")
        prog_dlg.transient(self.root)
        prog_dlg.resizable(False, False)
        _safe(prog_dlg.geometry, f"+{self.root.winfo_x()+200}+{self.root.winfo_y()+200}")
        ttk.Label(prog_dlg, text=f"Querying JPL Horizons for '{target}'...", padding=10).pack()
        prog_var = tk.DoubleVar(value=0)
        prog_bar = ttk.Progressbar(prog_dlg, variable=prog_var, maximum=100, length=320)
        prog_bar.pack(padx=20, pady=(0, 10))
        status_lbl = ttk.Label(prog_dlg, text="Starting...", padding=(10, 0))
        status_lbl.pack()
        prog_dlg.update()

        def progress_cb(current, total):
            try:
                pct = int(current / max(total, 1) * 100)
                prog_var.set(pct)
                status_lbl.config(text=f"Step {current}/{total}")
                prog_dlg.update()
            except Exception:
                pass

        # ── 4. Query ─────────────────────────────────────────────────────────
        self.root.config(cursor="watch")
        try:
            jpl_df = iterative_lighttime_correction(target, epochs, location, progress_callback=progress_cb)
        except Exception as e:
            prog_dlg.destroy()
            self.root.config(cursor="")
            messagebox.showerror("JPL Query Error", f"Failed to fetch JPL data:\n{e}")
            return False
        finally:
            self.root.config(cursor="")

        prog_dlg.destroy()

        # ── 5. Attach columns to lc.df ───────────────────────────────────────
        cols_to_add = [
            'r', 'delta', 'alpha_true', 'ObsEclLon', 'ObsEclLat',
            'corrected_jd', 'lighttime_days',
            'delta_uncorr', 'r_uncorr', 'alpha_true_uncorr',
            'ObsEclLon_uncorr', 'ObsEclLat_uncorr',
        ]
        for col in cols_to_add:
            if col in jpl_df.columns:
                lc.df[col] = jpl_df[col].values

        # Reduced magnitudes (corrected distances)
        if 'r' in lc.df.columns and 'delta' in lc.df.columns:
            lc.df['reduced_mag'] = calc_reduced_mag(lc.df['mag'], lc.df['r'], lc.df['delta'])

        # Reduced magnitudes (uncorrected distances – for reference)
        if 'r_uncorr' in lc.df.columns and 'delta_uncorr' in lc.df.columns:
            lc.df['reduced_mag_uncorr'] = calc_reduced_mag(
                lc.df['mag'], lc.df['r_uncorr'], lc.df['delta_uncorr']
            )

        lc._fill_arrays_cache()
        messagebox.showinfo("JPL Data", f"JPL data fetched successfully for '{target}'.")
        return True

    def _rebuild_band_menu(self):
        """Rebuild the Filters menu and keep GUI + plotting state in sync.
        Semantics:
          - self.selected_bands == set([...]) -> only those bands are shown
          - self.selected_bands == set()      -> show NONE (useful for 'Select None')
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
        """Update the Filters button label to summarise the active band selection."""
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
        """Handle a change in the Filters menu checkbuttons."""
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
        """Schedule a full plot update via ``self.plot.update()``."""
        self.plot.update(self.lightcurves, self.lc_offsets, self.lc_visible, self.mode, self.time_mode,
                         self.show_rejected, self.errorbar_type, self.selected_bands)

    def on_click(self, event):
        """Handle a matplotlib mouse-click: select the nearest data point."""
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
        """Toggle the rejection flag of the currently selected point and redraw."""
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
        """Deselect the current point and hide the selection marker."""
        self.plot.selection = {'alias': None, 'index': None}
        self.current_point_index = None
        self.current_point_band = None
        self._update_selected_point_label()
        self.request_plot_update()
        self.asteroid_viewer.update_image()  # Update (clear) image window

    def move_point(self, direction: int):
        """Move the selection cursor by *direction* steps (+1 right, -1 left)."""
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
        """Refresh the status-bar label showing info about the selected point."""
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
        """Toggle the asteroid image viewer for the selected point."""
        self.asteroid_viewer.toggle_visibility()

    # ——— Rotation period ———
    def adjust_rotation_period(self, delta: float):
        """Change the rotation period by *delta* steps and trigger a redraw."""
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
        """Read the period entry field and refresh the rotation-phase plot."""
        try:
            val = float(self.rotation_period_var.get())
        except Exception:
            val = DEFAULT_PERIOD
        self.rotation_period_var.set(f"{val:.3f}")

        if self.time_mode != 'rotation_phase':
            return
        self.plot.invalidate_layout()
        self.request_plot_update()

    def _on_phase_max_changed(self, *_):
        """Validate phase_max entry and trigger a plot rebuild."""
        try:
            val = float(self.phase_max_var.get())
            if not (1.0 < val <= 3.0):
                val = DEFAULT_PHASE_MAX
        except Exception:
            val = DEFAULT_PHASE_MAX
        self.phase_max_var.set(f"{val:.2f}")
        if self.time_mode != 'rotation_phase':
            return
        self.plot.invalidate_layout()
        self.request_plot_update()

    def on_pick_label(self, event):
        """Handle a Matplotlib pick event on axis labels/title to allow inline editing."""
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
                cur = getattr(plot, '_user_xlabel', None) or getattr(plot, 'xlabel', '')
                new = simpledialog.askstring("Edit X label", "Enter new X axis label:", initialvalue=cur,
                                             parent=self.root)
                if new is None: return
                # Empty string clears the override and reverts to auto-generated label
                plot._user_xlabel = new.strip() if new.strip() else None
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
        """Open the marker style/size settings dialog."""
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
        """Ask for confirmation and quit the application."""
        self.on_close()

    def on_close(self):
        """Handle the window close (WM_DELETE_WINDOW) protocol."""
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

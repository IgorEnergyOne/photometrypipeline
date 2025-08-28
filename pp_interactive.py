#!/usr/bin/env python3

import os
import subprocess
import numpy as np
import pandas as pd
from astropy.time import Time
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

import ttkbootstrap as ttk
from ttkbootstrap.constants import *
from tkinter import filedialog, messagebox, simpledialog

# ---------------------- Constants & Utilities ----------------------
DEFAULT_COLORS = [
    "blue", "red", "green", "orange", "purple",
    "brown", "pink", "gray", "olive", "cyan"
]
BAND_COLORS = {"B": "blue", "V": "green", "R": "red", "I": "indigo"}

# Marker map (marker symbol, human name)
MARKERS = [
    ('o', 'circle'), ('s', 'square'), ('p', 'pentagon'), ('x', 'x'), ('D', 'diamond'),
    ('*', 'star'), ('v', 'triangle_down'), ('^', 'triangle_up'), ('<', 'triangle_left'),
    ('>', 'triangle_right'), ('+', 'plus'), ('d', 'thin_diamond'),
]

TIME_STEP = 0.02 # hours
DEFAULT_PERIOD = 4.0 # hours


def next_version(path: str) -> str:
    """Return path like originalname_1.ext, originalname_2.ext … choosing the first free number."""
    base, ext = os.path.splitext(path)
    i = 1
    while os.path.exists(f"{base}_{i}{ext}"):
        i += 1
    return f"{base}_{i}{ext}"


# ---------------------- Data Layer ----------------------
class LightCurveData:
    """Single-responsibility: loading/saving and light data utilities."""

    def __init__(self):
        self.df: pd.DataFrame | None = None
        self.filename: str | None = None

    def load(self, filepath: str) -> None:
        self.filename = filepath
        self.df = pd.read_csv(filepath)
        # Add missing columns with defaults
        if 'rejected' not in self.df.columns:
            self.df['rejected'] = False
        if 'sextractor_flags' not in self.df.columns:
            self.df['sextractor_flags'] = 0

    def save(self) -> None:
        if self.df is not None and self.filename:
            self.df.to_csv(self.filename, index=False)

    def toggle_rejection(self, index: int) -> None:
        self.df.loc[index, 'rejected'] = not self.df.loc[index, 'rejected']

    # Band utilities
    def get_bands(self) -> list[str]:
        if self.df is None:
            return []
        if 'band' not in self.df.columns:
            return []
        vals = self.df['band'].dropna().astype(str).unique().tolist()
        # Stable ordering: common photometric order if present
        desired = ['U', 'B', 'V', 'R', 'I', 'g', 'r', 'i', 'z']
        # Preserve encountered order but bias by desired order
        present = {b: i for i, b in enumerate(vals)}
        ordered = sorted(vals, key=lambda b: (desired.index(b) if b in desired else 1_000 + present[b]))
        return ordered

    def subset_by_bands(self, bands: list[str]) -> dict[str, pd.DataFrame]:
        """Return a dict mapping band -> DataFrame (filtered) for provided bands present in df."""
        if self.df is None or 'band' not in self.df.columns:
            return {}
        out = {}
        for b in bands:
            sub = self.df[self.df['band'].astype(str) == str(b)]
            if not sub.empty:
                out[b] = sub.reset_index(drop=True)
        return out


# ---------------------- Plotting Layer ----------------------
class LightCurvePlot:
    """Owns matplotlib figure/axes, drawing, and interactive state."""

    def __init__(self, figure: matplotlib.figure.Figure, ax: matplotlib.axes.Axes, master_frame):
        self.fig = figure
        self.ax = ax
        self.master_frame = master_frame

        # Interaction state
        self.selected_index: int | None = None
        self.x_data: np.ndarray | None = None

        # Labels
        self.xlabel = "Julian Date"
        self.ylabel = "Magnitude"
        self.title = "Lightcurve"
        self.auto_title = True

        # Colors
        self.valid_color = "blue"      # used for single-band mode
        self.rejected_color = "red"
        self.flagged_color = "orange"

        # Marker & errorbar appearance
        self.marker_size = 2.0
        self.marker_style = 's'
        self.errorbar_capsize = 2.0
        self.errorbar_capthick = 1.0
        self.errorbar_linewidth = 1.0

        # Marker dictionary
        self.available_markers = MARKERS
        self.marker_dict = {name: marker for marker, name in self.available_markers}

        # Legend frame handle (for color pickers)
        self.legend_frame = None
        self.legend_patches = {}

        # Back-reference to GUI (set by GUI after creation)
        self.parent_gui = None

    # ---- basic helpers ----
    def clear(self):
        self.ax.clear()

    def set_labels(self):
        self.ax.set_xlabel(self.xlabel)
        self.ax.set_ylabel(self.ylabel)
        self.ax.set_title(self.title)
        # Allow clicking these text artists
        self.ax.title.set_picker(True)
        self.ax.xaxis.label.set_picker(True)
        self.ax.yaxis.label.set_picker(True)
        self.ax.invert_yaxis()

    def draw(self):
        self.fig.canvas.draw_idle()

    def target_name_parser(self, df: pd.DataFrame) -> str:
        try:
            return str(df['target'].iloc[0])
        except Exception:
            return 'Lightcurve'

    # ---- main update ----
    def update(self, df: pd.DataFrame, mode: str, time_mode: str, show_rejected: bool,
               errorbar_type: str = 'calibrated', selected_band: str = 'All', update_legend: bool = True) -> None:
        self.clear()
        if df is None or df.empty:
            self.draw()
            return

        # X axis handling
        try:
            jd = df['julian_date']
        except KeyError:
            messagebox.showwarning("Warning", "No date is available (wrong file?).")
            return

        # Build list of bands to plot
        if 'band' in df.columns:
            if selected_band == 'All':
                bands = df['band'].dropna().astype(str).unique().tolist()
            else:
                bands = [selected_band]
        else:
            # No band column: behave like single-band
            bands = [None]

        # Assign colors for bands (only used in All mode)
        band_colors = {b: DEFAULT_COLORS[i % len(DEFAULT_COLORS)] for i, b in enumerate(bands)}

        # For label management
        plotted_any = False

        for b in bands:
            sub = df if b is None else df[df['band'].astype(str) == str(b)]
            if sub.empty:
                continue

            if time_mode == 'julian_date':
                x = sub['julian_date']
                self.xlabel = "Julian Date"
            elif time_mode == 'mjd':
                x = sub['julian_date'] - 2400000.5
                self.xlabel = "Modified Julian Date (MJD)"
            elif time_mode == 'minutes':
                x = (sub['julian_date'] - sub['julian_date'].min()) * 24 * 60
                try:
                    self.xlabel = f"Minutes from {Time(sub['julian_date'].min(), format='jd').to_value('iso', subfmt='date_hm')} UT"
                except Exception:
                    self.xlabel = "Minutes"
            elif time_mode == 'rotation_phase':
                if hasattr(self.parent_gui, "rotation_period_var"):
                    period_hours = round(float(self.parent_gui.rotation_period_var.get()), 4)
                    jd0 = sub['julian_date'].min()
                    phases = ((sub['julian_date'] - jd0) * 24.0 / period_hours) * 2 * np.pi
                    phases = np.mod(phases, 2 * np.pi)
                    x = phases
                    self.xlabel = "Rotation Phase [rad]"
                else:
                    x = sub['julian_date']  # fallback
                    self.xlabel = "Julian Date"
            else:
                raise ValueError(f"Invalid time_mode: {time_mode}")

            # Y vector per mode
            if mode == 'target':
                y = sub['mag'].to_numpy()
            elif mode == 'instrumental':
                y = sub['inst_mag'].to_numpy()
            elif mode == 'control':
                if 'mag_control' not in sub.columns:
                    messagebox.showwarning("Warning", "No control star data available.")
                    return
                y = sub['mag_control'].to_numpy()
            else:
                raise ValueError(f"Invalid mode: {mode}")

            # Errorbars
            if errorbar_type == 'instrumental':
                yerr = sub['inst_sig'].to_numpy()
            elif errorbar_type == 'calibrated':
                col = 'sig' if mode != 'control' else 'sig_control'
                if col not in sub.columns:
                    yerr = None
                else:
                    yerr = sub[col].to_numpy()
            else:  # 'none'
                yerr = None

            # Masks
            mask_valid = ~sub['rejected'] & (sub['sextractor_flags'] == 0)
            mask_rejected = sub['rejected']
            mask_flagged = (sub['sextractor_flags'] > 0) & (~sub['rejected'])

            # Style
            if selected_band == 'All' and b is not None:
                valid_color = band_colors[b]
            else:
                valid_color = self.valid_color

            errorbar_kwargs = dict(
                fmt=self.marker_style,
                markersize=self.marker_size,
                capsize=self.errorbar_capsize,
                capthick=self.errorbar_capthick,
                elinewidth=self.errorbar_linewidth,
                picker=5,
            )
            plot_kwargs = dict(marker=self.marker_style, markersize=self.marker_size, linestyle='None', picker=5)

            X = x.to_numpy()
            if self.x_data is None or len(self.x_data) != len(df):
                # store total x for selection distance in single-band context; we'll refresh later
                self.x_data = df['julian_date'].to_numpy()

            # Plotting
            if yerr is not None:
                self.ax.errorbar(X[mask_valid], y[mask_valid], yerr=yerr[mask_valid], color=valid_color, **errorbar_kwargs)
                self.ax.errorbar(X[mask_flagged], y[mask_flagged], yerr=yerr[mask_flagged], color=self.flagged_color, **errorbar_kwargs)
                if show_rejected:
                    self.ax.errorbar(X[mask_rejected], y[mask_rejected], yerr=yerr[mask_rejected], color=self.rejected_color, **errorbar_kwargs)
            else:
                self.ax.plot(X[mask_valid], y[mask_valid], color=valid_color, **plot_kwargs)
                self.ax.plot(X[mask_flagged], y[mask_flagged], color=self.flagged_color, **plot_kwargs)
                if show_rejected:
                    self.ax.plot(X[mask_rejected], y[mask_rejected], color=self.rejected_color, **plot_kwargs)

            # Label once per band
            if selected_band == 'All' and b is not None:
                # add a tiny dummy point for legend label
                self.ax.plot([], [], marker=self.marker_style, linestyle='None', color=valid_color, label=f"{b}")

            plotted_any = True

        # Selected point marker/label (only meaningful if using a single series)
        if self.selected_index is not None and self.selected_index < len(df):
            # Choose x,y in current mode/time for the selected index
            try:
                x_click_jd = df['julian_date'].iloc[self.selected_index]
                if time_mode == 'julian_date':
                    sx = x_click_jd
                elif time_mode == 'mjd':
                    sx = x_click_jd - 2400000.5
                elif time_mode == 'minutes':
                    # minutes relative to min of its band; use df-wide min to keep consistent
                    sx = (x_click_jd - df['julian_date'].min()) * 24 * 60
                elif time_mode == 'rotation_phase':  # and hasattr(self.parent_gui, 'rotation_period_var'):
                    period_hours = float(self.parent_gui.rotation_period_var.get())
                    jd0 = df['julian_date'].min()
                    sx = np.mod((x_click_jd - jd0) * 24.0 / period_hours * 2 * np.pi, 2 * np.pi)
                if mode == 'target':
                    sy = df['mag'].iloc[self.selected_index]
                elif mode == 'instrumental':
                    sy = df['inst_mag'].iloc[self.selected_index]
                else:
                    sy = df['mag_control'].iloc[self.selected_index]

                self.ax.plot(sx, sy, 'o', color='orange', markersize=10)
                label = os.path.splitext(os.path.basename(df.loc[self.selected_index, 'filename']))[0] if 'filename' in df.columns else str(self.selected_index)
                self.ax.text(sx, sy + 0.1, label, color='orange')
            except Exception:
                pass

        # Legend
        if update_legend:
            self.update_legend(selected_band=selected_band)

        # Auto title
        if self.auto_title:
            self.title = self.target_name_parser(df)

        # Axes labels and draw
        self.set_labels()
        # If All bands, add a legend of bands
        if selected_band == 'All' and plotted_any:
            self.ax.legend(title='Band')
        self.draw()

    # ---- legend with clickable color pickers for valid/rejected/flagged ----
    def update_legend(self, selected_band: str = 'All') -> None:
        # Destroy pre-existing legend frame
        if self.legend_frame and self.legend_frame.winfo_exists():
            self.legend_frame.destroy()

        self.legend_frame = ttk.Frame(self.master_frame)
        self.legend_frame.pack(before=self.fig.canvas.get_tk_widget(), side="top", pady=5)

        legend_items_frame = ttk.Frame(self.legend_frame)
        legend_items_frame.pack()

        self.legend_handles = {
            'valid': self.valid_color,
            'rejected': self.rejected_color,
            'flagged': self.flagged_color,
        }

        legend_items = [
            ("Valid", 'valid'),
            ("Rejected", 'rejected'),
            ("Flagged", 'flagged'),
        ]

        # In All-bands mode, clarify that Valid color only applies in single-band mode
        if selected_band == 'All':
            note = ttk.Label(self.legend_frame, text="Note: Per-filter colors are auto-assigned in 'All' mode.")
            note.pack(pady=(0, 5))

        for label, color_key in legend_items:
            frame = ttk.Frame(legend_items_frame)
            frame.pack(side=LEFT, padx=10)

            color = self.legend_handles[color_key]
            patch_frame = ttk.Frame(frame, width=20, height=20)
            patch_frame.pack(side=LEFT, padx=2)
            patch_frame.pack_propagate(False)

            patch = ttk.Label(patch_frame, background=color, borderwidth=1, relief='solid')
            patch.pack(fill='both', expand=True)
            self.legend_patches[color_key] = patch

            ttk.Label(frame, text=label).pack(side=LEFT, padx=2)

            def on_patch_click(event, key=color_key):
                self.on_legend_click(key)

            patch.bind('<Button-1>', on_patch_click)
            patch_frame.bind('<Button-1>', on_patch_click)
            for widget in frame.winfo_children():
                if isinstance(widget, ttk.Label):
                    widget.bind('<Button-1>', on_patch_click)

    def on_legend_click(self, color_key: str) -> None:
        # Color picker dialog implemented via a combobox of DEFAULT_COLORS
        color_dialog = ttk.Toplevel()
        color_dialog.title(f"Select {color_key.capitalize()} Color")
        color_dialog.transient(self.master_frame)
        color_dialog.grab_set()

        x = self.master_frame.winfo_pointerx()
        y = self.master_frame.winfo_pointery()
        color_dialog.geometry(f"+{x}+{y}")

        frame = ttk.Frame(color_dialog, padding=10)
        frame.pack(fill=BOTH, expand=True)

        ttk.Label(frame, text=f"Select color for {color_key} points:").pack(pady=5)

        color_var = ttk.StringVar(value=self.legend_handles[color_key])
        color_combo = ttk.Combobox(frame, textvariable=color_var, values=DEFAULT_COLORS, state='readonly', width=15)
        color_combo.pack(pady=5)

        preview_frame = ttk.Frame(frame, height=30, width=100)
        preview_frame.pack_propagate(False)
        preview_frame.pack(pady=5)
        preview = ttk.Frame(preview_frame, style=f"{color_var.get().title()}.TFrame")
        preview.pack(fill=BOTH, expand=True)

        def update_preview(*_):
            preview.configure(style=f"{color_var.get().title()}.TFrame")

        color_var.trace_add('write', lambda *_: update_preview())

        for c in DEFAULT_COLORS:
            ttk.Style().configure(f"{c.title()}.TFrame", background=c)

        btns = ttk.Frame(frame)
        btns.pack(pady=5)

        def apply_color():
            c = color_var.get()
            if c in DEFAULT_COLORS:
                if color_key == 'valid':
                    self.valid_color = c
                elif color_key == 'rejected':
                    self.rejected_color = c
                elif color_key == 'flagged':
                    self.flagged_color = c

                if color_key in self.legend_patches:
                    self.legend_patches[color_key].configure(background=c)

                # Trigger a redraw with new colors
                if hasattr(self, 'parent_gui') and self.parent_gui and self.parent_gui.data.df is not None:
                    self.update(
                        self.parent_gui.data.df,
                        self.parent_gui.mode,
                        self.parent_gui.time_mode,
                        self.parent_gui.show_rejected,
                        self.parent_gui.errorbar_type,
                        self.parent_gui.selected_band,
                        update_legend=False,
                    )
            color_dialog.destroy()

        ttk.Button(btns, text="OK", command=apply_color).pack(side=LEFT, padx=5)
        ttk.Button(btns, text="Cancel", command=color_dialog.destroy).pack(side=LEFT, padx=5)

        color_dialog.bind('<Return>', lambda e: apply_color())
        color_dialog.bind('<Escape>', lambda e: color_dialog.destroy())
        update_preview()

    # ---- selection helpers ----
    def select_point(self, x_click: float, y_click: float, df: pd.DataFrame, mode: str, time_mode: str) -> None:
        # For distance calculation, convert x into JD-scale for stability
        if time_mode == 'julian_date':
            x_all = df['julian_date']
        elif time_mode == 'mjd':
            x_all = df['julian_date'] - 2400000.5
        elif time_mode == 'minutes':
            # minutes relative to min of its band; use df-wide min to keep consistent
            x_all = (df['julian_date'] - df['julian_date'].min()) * 24 * 60
        elif time_mode == 'rotation_phase': #and hasattr(self.parent_gui, 'rotation_period_var'):
            period_hours = float(self.parent_gui.rotation_period_var.get())
            jd0 = df['julian_date'].min()
            x_all = np.mod((df['julian_date'] - jd0) * 24.0 / period_hours * 2 * np.pi, 2 * np.pi)
        x_all = x_all.to_numpy()
        if mode == 'target':
            y_all = df['mag'].to_numpy()
        elif mode == 'instrumental':
            y_all = df['inst_mag'].to_numpy()
        else:
            y_all = df['mag_control'].to_numpy()

        # Normalize all distances to 1
        def normalize(x): return (x - np.min(x)) / (np.max(x) - np.min(x)), np.min(x), np.max(x)
        x_all, x_min, x_max = normalize(x_all)
        y_all, y_min, y_max = normalize(y_all)
        # normalize click coordinates
        x_click_norm = (x_click - x_min) / (x_max - x_min)
        y_click_norm = (y_click - y_min) / (y_max - y_min)
        # Note: This is a best-effort selection when multiple bands shown.
        distances = np.sqrt((x_all - x_click_norm)**2 + (y_all - y_click_norm)**2)
        idx = int(np.argmin(distances))
        self.selected_index = idx if distances[idx] < 0.05 else None

    def move_selection(self, direction: int) -> None:
        if self.parent_gui is None or self.parent_gui.data.df is None:
            return
        n = len(self.parent_gui.data.df)
        if n == 0:
            return
        if self.selected_index is None:
            self.selected_index = 0 if direction > 0 else n - 1
        else:
            self.selected_index = max(0, min(n - 1, self.selected_index + direction))


# ---------------------- GUI Layer ----------------------
class LightCurveGUI:
    def __init__(self, root: ttk.Window):
        self.root = root
        self.root.title("Lightcurve Viewer and Editor")

        # Model
        self.data = LightCurveData()

        # View state
        self.mode = 'target'              # 'target' | 'instrumental' | 'control'
        self.time_mode = 'minutes'         # 'minutes' | 'julian_date' | 'mjd'
        self.show_rejected = True
        self.errorbar_type = 'calibrated'  # 'instrumental' | 'calibrated' | 'none'
        self.selected_band = 'All'         # 'All' or concrete band value

        # Build UI
        self.create_widgets()

        # Window / key bindings
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)
        self.root.bind("q", self.confirm_exit)
        self.root.bind("r", self.toggle_rejection)
        self.root.bind("a", self.cancel_selection)
        self.root.bind("<Left>", self.move_left)
        self.root.bind("<Right>", self.move_right)

    # ---- UI construction ----
    def create_widgets(self) -> None:
        self.master_frame = ttk.Frame(self.root)
        self.master_frame.pack(fill=BOTH, expand=True)

        control_frame = ttk.Frame(self.master_frame, padding=10)
        control_frame.pack(side=TOP, fill=X)

        # File ops
        ttk.Button(control_frame, text="Open CSV", command=self.open_csv).pack(side=LEFT, padx=5)
        # Save menu button
        self.save_menu_btn = ttk.Menubutton(control_frame, text="Save")
        self.save_menu = ttk.Menu(self.save_menu_btn, tearoff=0)
        self.save_menu.add_command(label="Save CSV", command=self.save_csv)
        self.save_menu.add_command(label="Save Plot", command=self.save_plot)
        self.save_menu.add_command(label="Save Atlas", command=self.save_atlas)  # new placeholder
        self.save_menu_btn["menu"] = self.save_menu
        self.save_menu_btn.pack(side=ttk.LEFT)
        # ttk.Button(control_frame, text="Save CSV", command=self.save_csv).pack(side=LEFT, padx=5)
        # ttk.Button(control_frame, text="Save Plot", command=self.save_plot).pack(side=LEFT, padx=5)

        # Mode
        ttk.Label(control_frame, text="Mode:").pack(side=LEFT)
        self.mode_var = ttk.StringVar(value=self.mode)
        ttk.Combobox(control_frame, textvariable=self.mode_var, values=['target', 'instrumental', 'control'], state='readonly', width=13).pack(side=LEFT)
        self.mode_var.trace_add('write', lambda *_: self.set_mode())

        # Errorbar type
        ttk.Label(control_frame, text="Error Bars:").pack(side=LEFT, padx=(10, 0))
        self.errorbar_var = ttk.StringVar(value=self.errorbar_type)
        ttk.Combobox(control_frame, textvariable=self.errorbar_var, values=['calibrated', 'instrumental', 'none'], state='readonly', width=12).pack(side=LEFT)
        self.errorbar_var.trace_add('write', lambda *_: self.set_errorbar_type())

        # Time axis
        ttk.Label(control_frame, text="Time Axis:").pack(side=LEFT)
        self.time_var = ttk.StringVar(value=self.time_mode)
        ttk.Combobox(
            control_frame,
            textvariable=self.time_var,
            values=['minutes', 'julian_date', 'mjd', 'rotation_phase'],  # added
            state='readonly',
            width=13
        ).pack(side=LEFT)
        self.time_var.trace_add('write', lambda *_: self.set_time_mode())

        # Show rejected
        self.toggle_rejected_var = ttk.BooleanVar(value=self.show_rejected)
        ttk.Checkbutton(control_frame, text="Show Rejected", variable=self.toggle_rejected_var, command=self.set_show_rejected).pack(side=LEFT, padx=5)

        # Plot area
        plot_frame = ttk.Frame(self.master_frame)
        plot_frame.pack(fill=BOTH, expand=True, padx=10, pady=5)
        self.fig, self.ax = plt.subplots(figsize=(8, 4))
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=True)
        NavigationToolbar2Tk(self.canvas, plot_frame).update()

        self.plot = LightCurvePlot(self.fig, self.ax, self.master_frame)
        self.plot.parent_gui = self

        # Band selection
        ttk.Label(control_frame, text="Filter:").pack(side=LEFT, padx=(10, 0))
        self.band_var = ttk.StringVar(value=self.selected_band)
        self.band_combo = ttk.Combobox(control_frame, textvariable=self.band_var, state='readonly', width=10)
        self.band_combo.pack(side=LEFT)
        self.band_var.trace_add('write', lambda *_: self.set_band())
        # Initially only 'All'
        self.band_combo['values'] = ['All']
        self.band_var.set('All')

        # pack the matplotlib canvas
        self.canvas.get_tk_widget().pack(fill=BOTH, expand=True)

        # Rotation controls frame (under plot, hidden until needed)
        self.rotation_frame = ttk.Frame(self.master_frame)  # attach to master, not top bar
        #self.rotation_frame.pack(fill=X, pady=5)

        ttk.Label(self.rotation_frame, text="Rotation period (h):").pack(side=LEFT, padx=5)

        if self.data.df is not None and "julian_date" in self.data.df.columns:
            time_span = (self.data.df["julian_date"].max() - self.data.df["julian_date"].min()) / 24.0
        else:
            time_span = DEFAULT_PERIOD  # fallback default in hours
        self.rotation_period_var = ttk.StringVar(value=f"{time_span:.3f}")

        self.rotation_slider = ttk.Scale(
            self.rotation_frame, from_=0.1, to=100.0,
            variable=self.rotation_period_var, orient=HORIZONTAL,
            length=int(self.master_frame.winfo_screenwidth() / 2),  # ~half window width
            command=lambda v: self.update_rotation_period()
        )
        self.rotation_slider.pack(side=LEFT, padx=5, expand=True, fill=X)

        entry = ttk.Entry(self.rotation_frame, textvariable=self.rotation_period_var, width=8)
        entry.pack(side=LEFT, padx=5)
        entry.bind("<Return>", lambda e: self.update_rotation_period())

        self.root.bind("z", lambda e: self.adjust_rotation_period(-TIME_STEP))
        self.root.bind("x", lambda e: self.adjust_rotation_period(+TIME_STEP))
        entry.bind("<Return>", lambda e: self.update_rotation_period())

        # Marker settings
        ttk.Button(control_frame, text="Marker Settings", command=self.show_marker_settings).pack(side=LEFT, padx=8)


        # Matplotlib event connections
        self.canvas.mpl_connect("button_press_event", self.on_click)
        self.canvas.mpl_connect("pick_event", self.on_pick)

        # Help panel
        help_frame = ttk.Frame(self.master_frame, padding=10)
        help_frame.pack(side=BOTTOM, fill=X)
        ttk.Label(
            help_frame,
            text=(
                "Hotkeys: [r] toggle rejection | [a] cancel selection | [q] quit | "
                "[left]/[right] move selection | click = select/unselect or edit title/labels | "
                "Filter dropdown = select filter or 'All'"
            ),
        ).pack()

    # ---- UI actions ----
    def open_csv(self) -> None:
        file = filedialog.askopenfilename(initialdir=os.getcwd(), filetypes=[("CSV Files", "*.csv")])
        if not file:
            return
        self.data.load(file)
        # put filename in title bar
        self.root.title(f"Lightcurve Viewer – {os.path.basename(file)}")
        # Populate band combobox
        bands = self.data.get_bands()
        if not bands:
            self.band_combo['values'] = ['All']
            self.band_var.set('All')
        elif len(bands) == 1:
            # Only one band → no "All", auto-select that band
            self.band_combo['values'] = bands
            self.band_var.set(bands[0])
            self.selected_band = bands[0]
        else:
            # Multiple bands → include "All"
            self.band_combo['values'] = ['All'] + bands
            self.band_var.set('All')
            self.selected_band = 'All'
        # Initial plot
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, self.selected_band)

    def save_csv(self) -> None:
        if self.data.df is None or self.data.filename is None:
            messagebox.showerror("Nothing to save", "Load a CSV first.")
            return
        file = filedialog.asksaveasfilename(
            initialdir=os.path.dirname(self.data.filename),
            initialfile=os.path.basename(self.data.filename),
            defaultextension=".csv",
            filetypes=[("CSV Files", "*.csv")],
        )
        if not file:  # user cancelled -> write versioned next to current
            file = next_version(self.data.filename)
        self.data.filename = file
        self.data.save()
        messagebox.showinfo("Saved", f"CSV saved as {file}")

    def save_plot(self) -> None:
        if self.data.df is None:
            return
        default_png = os.path.splitext(self.data.filename or "plot")[0] + ".png"
        file = filedialog.asksaveasfilename(
            initialdir=os.path.dirname(default_png),
            initialfile=os.path.basename(default_png),
            defaultextension=".png",
            filetypes=[("PNG Files", "*.png")],
        )
        if not file:
            file = next_version(default_png)
        self.plot.fig.savefig(file, dpi=300, bbox_inches="tight")
        messagebox.showinfo("Saved", f"Plot saved as {file}")

    def save_atlas(self):
        if self.data.df is None:
            return
        # Ask for FITS/FTS file
        filetypes = [("FITS files", "*.fits *.fts"), ("All files", "*.*")]
        fits_filepath = filedialog.askopenfilename(
            title="Select Atlas FITS file",
            filetypes=filetypes
        )
        if fits_filepath is None:
            # give message that no files were selected
            messagebox.showerror("Nothing to save", "Select atlas FITS file")
            return

        atlas_name = f"{Path(self.data.filename).stem}.ATL"
        file = filedialog.asksaveasfilename(
            initialdir=os.path.dirname(atlas_name),
            initialfile=os.path.basename(atlas_name),
            defaultextension=".ATL",
            filetypes=[("Atlas Files", "*.ATL")],
        )
        # save data as a temporary file
        self.data.df.to_csv('tmp.csv', header=True)
        atlas_name = file
        atlas_cmd = f"pp_atlas -fname_header {fits_filepath} -fname_photo tmp.csv -fname_out {atlas_name}"
        subprocess.call(['/bin/sh', '-i', '-c', atlas_cmd])
        # remove tmp file
        os.remove('tmp.csv')
        # atlas = form_atlas(filename_header=fits_filepath, filename_photometry=self.data.df)
        # write_atlas(filename_atlas=self.data.filename, text_atlas=atlas)
        messagebox.showinfo("Saved", f"ATLAS file saved as {file}")

    def set_mode(self, _=None) -> None:
        self.mode = self.mode_var.get()
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, self.selected_band)

    def set_time_mode(self, _=None):
        self.time_mode = self.time_var.get()
        if self.time_mode == "rotation_phase":
            if not self.rotation_frame.winfo_ismapped():
                self.rotation_frame.pack(fill=X, pady=5)
        else:
            if self.rotation_frame.winfo_ismapped():
                self.rotation_frame.pack_forget()

        self.plot.update(
            self.data.df, self.mode, self.time_mode,
            self.show_rejected, self.errorbar_type, self.selected_band
        )

    def set_show_rejected(self) -> None:
        self.show_rejected = self.toggle_rejected_var.get()
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, self.selected_band)

    def set_errorbar_type(self, _=None) -> None:
        self.errorbar_type = self.errorbar_var.get()
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, self.selected_band)

    def set_band(self) -> None:
        self.selected_band = self.band_var.get()
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, self.selected_band)

    # ---- matplotlib events ----
    def on_click(self, event) -> None:
        if event.inaxes != self.ax or self.data.df is None:
            return
        if event.xdata is not None and event.ydata is not None:
            self.plot.select_point(event.xdata, event.ydata, self.data.df, self.mode, self.time_mode)
            self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, self.selected_band, update_legend=False)

    def on_pick(self, event) -> None:
        artist = event.artist
        if artist == self.ax.title:
            new = simpledialog.askstring("Edit Title", "New Title:", initialvalue=self.plot.title)
            if new:
                self.plot.title = new
                self.plot.auto_title = False
        elif artist == self.ax.xaxis.label:
            new = simpledialog.askstring("Edit X Label", "New X Label:", initialvalue=self.plot.xlabel)
            if new:
                self.plot.xlabel = new
        elif artist == self.ax.yaxis.label:
            new = simpledialog.askstring("Edit Y Label", "New Y Label:", initialvalue=self.plot.ylabel)
            if new:
                self.plot.ylabel = new
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, self.selected_band, update_legend=False)

    # ---- keyboard conveniences ----
    def toggle_rejection(self, _=None) -> None:
        if self.plot.selected_index is not None and self.data.df is not None:
            self.data.toggle_rejection(self.plot.selected_index)
            self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, self.selected_band, update_legend=False)

    def move_left(self, _=None) -> None:
        self.plot.move_selection(-1)
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, self.selected_band, update_legend=False)

    def move_right(self, _=None) -> None:
        self.plot.move_selection(1)
        self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, self.selected_band, update_legend=False)

    def cancel_selection(self, _=None) -> None:
        if self.plot.selected_index is not None:
            self.plot.selected_index = None
            self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, self.selected_band, update_legend=False)

    def adjust_rotation_period(self, delta: float):
        current = float(self.rotation_period_var.get())
        new_value = max(0.1, current + delta)
        self.rotation_period_var.set(f"{new_value:.4f}")
        self.update_rotation_period()

    def update_rotation_period(self, *_):
        if self.time_mode != "rotation_phase" or self.data.df is None:
            return
        period = round(float(self.rotation_period_var.get()), 4)

        # recompute phase column
        self.data.df["phase"] = ((self.data.df["julian_date"] % period) / period)

        # update scatter directly if exists
        if hasattr(self.plot, "sc"):
            self.plot.sc.set_offsets(
                np.c_[self.data.df["phase"], self.data.df["mag"]]
            )
            self.canvas.draw_idle()
        else:
            # fallback full redraw
            self.plot.update(
                self.data.df, self.mode, self.time_mode,
                self.show_rejected, self.errorbar_type, self.selected_band
            )

    def confirm_exit(self, _=None) -> None:
        self.on_close()

    # ---- marker settings dialog ----
    def show_marker_settings(self) -> None:
        dlg = ttk.Toplevel()
        dlg.title("Marker and Error Bar Settings")
        dlg.transient(self.root)
        dlg.grab_set()

        x = self.root.winfo_x() + 50
        y = self.root.winfo_y() + 50
        dlg.geometry(f"+{x}+{y}")

        frame = ttk.Frame(dlg, padding=10)
        frame.pack(fill=BOTH, expand=True)

        # Marker style
        ttk.Label(frame, text="Marker Style:").grid(row=0, column=0, sticky=W, pady=2)
        marker_names = [name for _, name in self.plot.available_markers]
        current_marker_name = next((name for m, name in self.plot.available_markers if m == self.plot.marker_style), 'circle')
        marker_var = ttk.StringVar(value=current_marker_name)
        ttk.Combobox(frame, textvariable=marker_var, values=marker_names, state='readonly').grid(row=0, column=1, sticky=EW, pady=2, padx=5)

        # Marker size
        ttk.Label(frame, text="Marker Size:").grid(row=1, column=0, sticky=W, pady=2)
        size_var = ttk.DoubleVar(value=self.plot.marker_size)
        ttk.Scale(frame, from_=1, to=30, variable=size_var, orient=HORIZONTAL).grid(row=1, column=1, sticky=EW, pady=2, padx=5)
        ttk.Entry(frame, textvariable=size_var, width=6).grid(row=1, column=2, sticky=W, pady=2, padx=5)

        # Error bar cap size
        ttk.Label(frame, text="Error Cap Size:").grid(row=2, column=0, sticky=W, pady=2)
        capsize_var = ttk.DoubleVar(value=self.plot.errorbar_capsize)
        ttk.Scale(frame, from_=0, to=20, variable=capsize_var, orient=HORIZONTAL).grid(row=2, column=1, sticky=EW, pady=2, padx=5)
        ttk.Entry(frame, textvariable=capsize_var, width=6).grid(row=2, column=2, sticky=W, pady=2, padx=5)

        # Error bar cap thickness
        ttk.Label(frame, text="Cap Thickness:").grid(row=3, column=0, sticky=W, pady=2)
        capthick_var = ttk.DoubleVar(value=self.plot.errorbar_capthick)
        ttk.Scale(frame, from_=0.0, to=10, variable=capthick_var, orient=HORIZONTAL).grid(row=3, column=1, sticky=EW, pady=2, padx=5)
        ttk.Entry(frame, textvariable=capthick_var, width=6).grid(row=3, column=2, sticky=W, pady=2, padx=5)

        # Error bar line width
        ttk.Label(frame, text="Error Bar Width:").grid(row=4, column=0, sticky=W, pady=2)
        linewidth_var = ttk.DoubleVar(value=self.plot.errorbar_linewidth)
        ttk.Scale(frame, from_=0.0, to=10, variable=linewidth_var, orient=HORIZONTAL).grid(row=4, column=1, sticky=EW, pady=2, padx=5)
        ttk.Entry(frame, textvariable=linewidth_var, width=6).grid(row=4, column=2, sticky=W, pady=2, padx=5)

        # Preview
        preview_frame = ttk.LabelFrame(frame, text="Preview", padding=5)
        preview_frame.grid(row=0, column=3, rowspan=5, padx=10, sticky=N+S)

        fig, ax = plt.subplots(figsize=(3, 2), dpi=80)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
        x_demo = [0.2, 0.5, 0.8]; y_demo = [0.5, 0.5, 0.5]; yerr_demo = [0.2, 0.2, 0.2]
        ax.errorbar(x_demo, y_demo, yerr=yerr_demo, fmt='o', color='blue',
                    markersize=size_var.get(), capsize=capsize_var.get(), capthick=capthick_var.get(),
                    elinewidth=linewidth_var.get())
        canvas = FigureCanvasTkAgg(fig, master=preview_frame)
        canvas.draw(); canvas.get_tk_widget().pack(fill=BOTH, expand=True)

        def update_preview(*_):
            try:
                marker_style = next((m for m, name in self.plot.available_markers if name == marker_var.get()), 'o')
                ax.clear(); ax.set_xticks([]); ax.set_yticks([]); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
                ax.errorbar(x_demo, y_demo, yerr=yerr_demo, fmt=marker_style, color='blue',
                            markersize=size_var.get(), capsize=capsize_var.get(), capthick=capthick_var.get(),
                            elinewidth=linewidth_var.get())
                canvas.draw_idle()
            except Exception as e:
                print(f"Preview update error: {e}")

        marker_var.trace_add('write', update_preview)
        size_var.trace_add('write', update_preview)
        capsize_var.trace_add('write', update_preview)
        capthick_var.trace_add('write', update_preview)
        linewidth_var.trace_add('write', update_preview)

        btns = ttk.Frame(frame); btns.grid(row=5, column=0, columnspan=4, pady=10)

        def apply_settings():
            self.plot.marker_style = next((m for m, name in self.plot.available_markers if name == marker_var.get()), 'o')
            self.plot.marker_size = size_var.get()
            self.plot.errorbar_capsize = capsize_var.get()
            self.plot.errorbar_capthick = capthick_var.get()
            self.plot.errorbar_linewidth = linewidth_var.get()
            if self.data.df is not None:
                self.plot.update(self.data.df, self.mode, self.time_mode, self.show_rejected, self.errorbar_type, self.selected_band)

        ttk.Button(btns, text="Apply", command=apply_settings).pack(side=LEFT, padx=5)
        ttk.Button(btns, text="Close", command=dlg.destroy).pack(side=LEFT, padx=5)

        dlg.resizable(True, False)
        dlg.focus_set()
        dlg.wait_window()

    # ---- lifecycle ----
    def on_close(self) -> None:
        if messagebox.askokcancel("Quit", "Do you want to quit?"):
            plt.close('all')
            self.root.destroy()
            self.root.quit()


# ---------------------- Entrypoint ----------------------
if __name__ == '__main__':
    root = ttk.Window(themename="flatly")
    app = LightCurveGUI(root)
    root.mainloop()

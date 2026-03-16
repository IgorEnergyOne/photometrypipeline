#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys

if "PHOTPIPEDIR" not in os.environ:
    os.environ["PHOTPIPEDIR"] = os.path.dirname(os.path.abspath(__file__))

import argparse
import glob
from pathlib import Path
import numpy as np
import pandas as pd
import astropy.io.fits as fits
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from astropy.visualization import ZScaleInterval, ImageNormalize

import tkinter as tk
from tkinter import filedialog, messagebox
import ttkbootstrap as ttk
from ttkbootstrap.constants import *

from catalog import catalog


# ─────────────────────────────────── constants ──────────────────────────────

FLAG_HELP = """\
SExtractor FLAG descriptions

    1   Aperture photometry is likely to be biased by neighbouring sources
    or by more than 10% of bad pixels in any aperture.
    2   The object has been deblended.
    4   At least one object pixel is saturated.
    8   The isophotal footprint of the detected object is truncated
        (too close to an image boundary).
    16  At least one photometric aperture is incomplete or corrupted
        (hitting buffer or memory limits).
    32  The isophotal footprint is incomplete or corrupted
        (hitting buffer or memory limits).
    64  A memory overflow occurred during deblending.
    128 A memory overflow occurred during extraction.

Flags are bitwise-combined, so a value of 3 means flags 1 AND 2 are set.
"""

# Flag → display colour
_FLAG_COLORS = {
    0:  'cyan',
    1:  'green',
    2:  'yellow',
    3:  'orange',          # user-specified
    4:  'purple',
    8:  'blue',
    16: '#c0392b',       # firebrick-red, user-specified
}


def get_flag_color(flag):
    return _FLAG_COLORS.get(int(flag), 'blue')


# ─────────────────────────────────── helpers ────────────────────────────────

def read_ldac(ldac_path):
    cat = catalog(ldac_path)
    if str(ldac_path).endswith('.db'):
        cat.read_database(filename=ldac_path)
    else:
        cat.read_ldac(filename=ldac_path)
    df = cat.data.to_pandas()
    if df is not None and 'ident' in df.columns:
        df = df[df['ident'].notna()].copy()
        df['ident'] = df['ident'].astype('Int64')
    return df


def show_flag_help(parent=None):
    """Open a simple Toplevel with the flag description text."""
    top = ttk.Toplevel(parent)
    top.title("Flag Descriptions (F1)")
    top.resizable(False, False)
    txt = tk.Text(top, wrap="word", width=62, height=22,
                  font=("Courier", 10), relief="flat", padx=10, pady=8)
    txt.insert("1.0", FLAG_HELP)
    txt.config(state="disabled")
    txt.pack(fill=BOTH, expand=YES, padx=6, pady=6)
    ttk.Button(top, text="Close", bootstyle=SECONDARY,
               command=top.destroy).pack(pady=(0, 8))


def _ask_save_file(parent, default_name):
    """Open a Save-As dialog rooted at *parent* (appears in front).

    The user sees only the filename field pre-filled with *default_name*;
    the dialog starts in the current working directory.  Returns the full
    chosen path as a string, or None if the user cancels.
    """
    path = filedialog.asksaveasfilename(
        parent=parent,
        title="Save CSV as…",
        initialdir=str(Path.cwd()),
        initialfile=default_name,
        defaultextension=".csv",
        filetypes=[("CSV files", "*.csv"), ("All files", "*")],
    )
    return path if path else None


# ──────────────────────────────── mode: single FITS ─────────────────────────

def mode_fits(fits_file, ldac_file=None):
    if ldac_file is None:
        ldac_file = fits_file.replace('.fits', '.ldac.db')
        if not os.path.exists(ldac_file):
            ldac_file = fits_file.replace('.fits', '.ldac')
            if not os.path.exists(ldac_file):
                print("Could not find matching LDAC file automatically.")
                return

    df = read_ldac(ldac_file)
    fits_data = fits.open(fits_file)
    data = fits_data[0].data

    # ── ttkbootstrap window ───────────────────────────────────────────────
    root = ttk.Window(themename="flatly")
    root.title(f"Stars – {Path(ldac_file).name}")
    root.geometry("1050x820")

    # Top toolbar
    toolbar_frame = ttk.Frame(root, padding=(6, 4))
    toolbar_frame.pack(side=TOP, fill=X)

    def export_csv():
        default = f"stars_mode_fits_{Path(fits_file).stem}.csv"
        out_path = _ask_save_file(root, default)
        if out_path is None:
            return
        df.to_csv(out_path, index=False)
        print(f"Exported data to {out_path}")

    ttk.Button(toolbar_frame, text="Export CSV", bootstyle=PRIMARY,
               command=export_csv).pack(side=LEFT, padx=4)

    ttk.Button(toolbar_frame, text="Help (F1)", bootstyle=SECONDARY,
               command=lambda: show_flag_help(root)).pack(side=LEFT, padx=4)

    root.bind("<F1>", lambda _e: show_flag_help(root))

    # ── matplotlib figure embedded in tk window ───────────────────────────
    fig, ax = plt.subplots(figsize=(9, 8))
    fig.tight_layout()

    try:
        norm = ImageNormalize(data, interval=ZScaleInterval())
    except Exception:
        norm = None

    if norm is not None:
        ax.imshow(data, cmap='gray', origin='lower', norm=norm)
    else:
        ax.imshow(data, cmap='gray', origin='lower', vmin=0, vmax=6500)

    flags = df['FLAGS'].unique()
    for flag in flags:
        flag_data = df[df['FLAGS'] == flag]
        s = flag_data['ISOAREA_IMAGE'] ** 0.5 * np.pi if 'ISOAREA_IMAGE' in flag_data else 15
        color = get_flag_color(flag)
        ax.scatter(flag_data['XWIN_IMAGE'], flag_data['YWIN_IMAGE'], s=s,
                   label=f'FLAG={int(flag)}', alpha=1.0,
                   edgecolors=color, facecolors='none',
                   picker=True, pickradius=5)

    ax.legend()
    ax.set_title(f"Stars in {Path(ldac_file).name}")

    annot = ax.annotate("", xy=(0, 0), xytext=(20, 20),
                        textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)

    def update_annot(ind, points, flag_data):
        pos = points.get_offsets()[ind[0]]
        annot.xy = pos
        if hasattr(fig, 'highlight_circle'):
            fig.highlight_circle.remove()
        fig.highlight_circle = ax.scatter(
            pos[0], pos[1], s=300, facecolors='none',
            edgecolors='red', linewidth=2, zorder=5)
        star = flag_data.iloc[ind[0]]
        mag = star.get('MAG_AUTO', np.nan)
        flux = star.get('FLUX_AUTO', np.nan)
        ident = star.get('ident', 0)
        ident_str = str(int(ident)) if pd.notnull(ident) else "unknown"
        annot.set_text(f"ID: {ident_str}\nMag: {mag:.2f}\nFlux: {flux:.2f}")
        annot.get_bbox_patch().set_alpha(0.8)

    def on_pick(event):
        if isinstance(event.artist, plt.matplotlib.collections.PathCollection):
            points = event.artist
            if hasattr(fig, 'highlight_circle') and points == fig.highlight_circle:
                return
            label = points.get_label()
            flag = float(label.split('=')[1])
            flag_data = df[df['FLAGS'] == flag]
            update_annot(event.ind, points, flag_data)
            annot.set_visible(True)
            fig.canvas.draw_idle()

    fig.canvas.mpl_connect('pick_event', on_pick)

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=BOTH, expand=YES)

    nav_bar = ttk.Frame(root)
    nav_bar.pack(side=BOTTOM, fill=X)
    NavigationToolbar2Tk(canvas, nav_bar)

    def _on_close_fits():
        plt.close('all')
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", _on_close_fits)
    root.mainloop()


# ──────────────────────────────── mode: directory ───────────────────────────

# Mag-type option labels and the column pairs they draw
#   key  → (inst_col, inst_err_col, label_text, axis_label)
_MAG_OPTIONS = {
    "Both (AUTO)":         ("MAG_AUTO",  "MAGERR_AUTO",  "Instrumental (MAG_AUTO)",  "Instrumental Magnitude"),
    "Both (APER)":         ("MAG_APER",  "MAGERR_APER",  "Instrumental (MAG_APER)",  "Instrumental Magnitude (APER)"),
    "Instrumental AUTO":   ("MAG_AUTO",  "MAGERR_AUTO",  "Instrumental (MAG_AUTO)",  "Instrumental Magnitude"),
    "Instrumental APER":   ("MAG_APER",  "MAGERR_APER",  "Instrumental (MAG_APER)",  "Instrumental Magnitude (APER)"),
    "Calibrated":          (None,        None,            None,                       None),
}

_MAG_COMBO_VALUES = [
    "Both (AUTO)",
    "Both (APER)",
    "Instrumental AUTO",
    "Instrumental APER",
    "Calibrated",
]

# All recognised calibrated-magnitude columns in order of preference.
# Each entry: (display_label, mag_col, err_col)
_CAL_FILTERS = [
    # Johnson-Cousins
    ("U",  "_Umag",  "_e_Umag"),
    ("B",  "_Bmag",  "_e_Bmag"),
    ("V",  "_Vmag",  "_e_Vmag"),
    ("R",  "_Rmag",  "_e_Rmag"),
    ("I",  "_Imag",  "_e_Imag"),
    # SDSS
    ("u",  "_umag",  "_e_umag"),
    ("g",  "_gmag",  "_e_gmag"),
    ("r",  "_rmag",  "_e_rmag"),
    ("i",  "_imag",  "_e_imag"),
    ("z",  "_zmag",  "_e_zmag"),
    # Gaia
    ("G",   "_Gmag",  "_e_Gmag"),
    ("G_BP","_BPmag", "_e_BPmag"),
    ("G_RP","_RPmag", "_e_RPmag"),
]


class StarDetailWindow:
    """ttkbootstrap window showing per-star photometry across all frames."""

    def __init__(self, parent_root, star_id, dfs, star_rows, ldac_files):
        self.dfs        = dfs
        # Store bare filenames for the CSV 'Filename' column
        self._filenames = [Path(f).name for f in ldac_files]

        star_id_str = str(int(star_id)) if pd.notnull(star_id) else str(star_id)

        # ── window ────────────────────────────────────────────────────────
        self.win = ttk.Toplevel(parent_root)
        self.win.title(f"Star ID: {star_id_str}")
        self.win.geometry("960x680")

        # ── controls toolbar ──────────────────────────────────────────────
        ctrl = ttk.Frame(self.win, padding=(8, 6))
        ctrl.pack(side=TOP, fill=X)

        # Mag Type drop-down
        ttk.Label(ctrl, text="Mag Type:").pack(side=LEFT, padx=(0, 4))
        self.mag_var = ttk.StringVar(value="Both (AUTO)")
        mag_cb = ttk.Combobox(ctrl, textvariable=self.mag_var,
                              values=_MAG_COMBO_VALUES,
                              state="readonly", width=18)
        mag_cb.pack(side=LEFT, padx=(0, 12))
        mag_cb.bind("<<ComboboxSelected>>", lambda _e: self._redraw())

        # Cal Filter drop-down (populated with detected filters)
        ttk.Label(ctrl, text="Cal. Filter:").pack(side=LEFT, padx=(0, 4))
        self.cal_filter_var = ttk.StringVar()
        self.cal_cb = ttk.Combobox(ctrl, textvariable=self.cal_filter_var,
                              state="readonly", width=7)
        self.cal_cb.pack(side=LEFT, padx=(0, 12))
        self.cal_cb.bind("<<ComboboxSelected>>", lambda _e: self._redraw())

        # Show Error drop-down
        ttk.Label(ctrl, text="Show Error:").pack(side=LEFT, padx=(0, 4))
        self.err_var = ttk.StringVar(value="Yes")
        err_cb = ttk.Combobox(ctrl, textvariable=self.err_var,
                              values=["Yes", "No"],
                              state="readonly", width=6)
        err_cb.pack(side=LEFT, padx=(0, 16))
        err_cb.bind("<<ComboboxSelected>>", lambda _e: self._redraw())

        # Export button
        ttk.Button(ctrl, text="Export CSV", bootstyle=PRIMARY,
                   command=self._export).pack(side=LEFT, padx=4)

        # Help button
        ttk.Button(ctrl, text="Help (F1)", bootstyle=SECONDARY,
                   command=lambda: show_flag_help(self.win)).pack(side=LEFT, padx=4)
        self.win.bind("<F1>", lambda _e: show_flag_help(self.win))

        # ── matplotlib figure ─────────────────────────────────────────────
        self.fig, self.ax_left = plt.subplots(figsize=(8, 5))
        self.ax_right = self.ax_left.twinx()

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.win)
        self.canvas.get_tk_widget().pack(fill=BOTH, expand=YES)

        nav_bar = ttk.Frame(self.win)
        nav_bar.pack(side=BOTTOM, fill=X)
        NavigationToolbar2Tk(self.canvas, nav_bar)

        self.update_data(star_id, star_rows)

    def update_data(self, star_id, star_rows):
        self.star_id   = star_id
        self.star_rows = star_rows

        star_id_str = str(int(star_id)) if pd.notnull(star_id) else str(star_id)
        self._star_id_str = star_id_str
        self.win.title(f"Star ID: {star_id_str}")

        # Pre-extract arrays (NaN-safe)
        def _col(rows, col):
            return [float(r.get(col, np.nan)) if col in r.index else np.nan
                    for r in rows]

        self.inst_mag_auto  = _col(star_rows, 'MAG_AUTO')
        self.inst_err_auto  = _col(star_rows, 'MAGERR_AUTO')
        self.inst_mag_aper  = _col(star_rows, 'MAG_APER')
        self.inst_err_aper  = _col(star_rows, 'MAGERR_APER')
        self.fluxes         = _col(star_rows, 'FLUX_AUTO')

        # Detect which calibrated filters are actually present and non-empty.
        self._cal_data: dict = {}
        for label, mag_col, err_col in _CAL_FILTERS:
            mags = _col(star_rows, mag_col)
            errs = _col(star_rows, err_col)
            if any(not np.isnan(m) for m in mags):
                self._cal_data[label] = (mags, errs)

        self._cal_labels = list(self._cal_data.keys()) or ["(none)"]
        default_cal = "R" if "R" in self._cal_data else self._cal_labels[0]

        # Update the Cal Filter combobox
        self.cal_cb.configure(values=self._cal_labels)
        if self._cal_labels == ["(none)"]:
            self.cal_cb.configure(state="disabled")
            self.cal_filter_var.set("(none)")
        else:
            self.cal_cb.configure(state="readonly")
            if self.cal_filter_var.get() not in self._cal_labels:
                self.cal_filter_var.set(default_cal)

        self._redraw()

    # ── drawing ───────────────────────────────────────────────────────────

    def _redraw(self):
        ax_l = self.ax_left
        ax_r = self.ax_right

        ax_l.cla()
        ax_r.cla()
        ax_r.set_visible(False)

        mag_type = self.mag_var.get()
        show_err = self.err_var.get() == "Yes"
        x = np.arange(len(self.dfs))

        show_both  = mag_type.startswith("Both")
        show_inst  = show_both or mag_type.startswith("Instrumental")
        show_cal   = show_both or mag_type == "Calibrated"
        use_aper   = "APER" in mag_type

        inst_mags = self.inst_mag_aper if use_aper else self.inst_mag_auto
        inst_errs = self.inst_err_aper if use_aper else self.inst_err_auto
        inst_label = "Instrumental Mag (MAG_APER)" if use_aper else "Instrumental Mag (MAG_AUTO)"
        inst_ylabel = "Instrumental Magnitude (APER)" if use_aper else "Instrumental Magnitude"

        # When showing both series simultaneously, shift them ±1/5 of one
        # index step so points don't sit on top of each other.
        JITTER = 0.1  # fraction of one index unit
        x_inst = x - JITTER if show_both else x
        x_cal  = x + JITTER if show_both else x

        # ── instrumental (left axis) ──────────────────────────────────────
        if show_inst and any(not np.isnan(m) for m in inst_mags):
            if show_err:
                ax_l.errorbar(x_inst, inst_mags, yerr=inst_errs,
                              fmt='o', linestyle='none', color='tab:blue',
                              label=inst_label, capsize=3)
            else:
                ax_l.scatter(x_inst, inst_mags, color='tab:blue',
                             label=inst_label, zorder=3)
            ax_l.set_ylabel(inst_ylabel, color='tab:blue')
            ax_l.tick_params(axis='y', labelcolor='tab:blue')

        # ── calibrated (right axis when "Both", else left) ────────────────
        # Resolve data for the currently selected calibrated filter
        sel_cal = self.cal_filter_var.get()
        if sel_cal in self._cal_data:
            cal_mags, cal_magerrs = self._cal_data[sel_cal]
        else:
            cal_mags, cal_magerrs = [], []
        cal_label = f"Calibrated Mag ({sel_cal})"

        if show_cal and any(not np.isnan(m) for m in cal_mags):
            ax_to_use = ax_r if show_both else ax_l
            if show_both:
                ax_r.set_visible(True)
                ax_r.yaxis.set_label_position("right")
                ax_r.yaxis.tick_right()
            if show_err:
                ax_to_use.errorbar(x_cal, cal_mags, yerr=cal_magerrs,
                                   fmt='o', linestyle='none', color='tab:red',
                                   label=cal_label, capsize=3)
            else:
                ax_to_use.scatter(x_cal, cal_mags, color='tab:red',
                                  label=cal_label, zorder=3)
            ax_to_use.set_ylabel(f"Calibrated Magnitude ({sel_cal})", color='tab:red')
            ax_to_use.tick_params(axis='y', labelcolor='tab:red')

        ax_l.set_xlabel("Image Index")
        ax_l.set_title(f"Star ID: {self._star_id_str}")

        # Restore integer ticks on the x-axis regardless of jitter
        ax_l.set_xticks(x)

        h1, l1 = ax_l.get_legend_handles_labels()
        h2, l2 = ax_r.get_legend_handles_labels() if show_both else ([], [])
        if h1 or h2:
            ax_l.legend(h1 + h2, l1 + l2, loc='upper left')

        self.fig.tight_layout()
        self.canvas.draw_idle()

    # ── export ────────────────────────────────────────────────────────────

    def _export(self):
        default = f"star_{self._star_id_str}_data.csv"
        out_path = _ask_save_file(self.win, default)
        if out_path is None:
            return

        # Build base columns
        data = {
            'Image_Index':           list(np.arange(len(self.dfs))),
            'Filename':              self._filenames,
            'Instrumental_MAG_AUTO': self.inst_mag_auto,
            'Err_MAG_AUTO':          self.inst_err_auto,
            'Instrumental_MAG_APER': self.inst_mag_aper,
            'Err_MAG_APER':          self.inst_err_aper,
            'Flux_AUTO':             self.fluxes,
        }
        # Append only the currently selected calibrated filter
        sel_cal = self.cal_filter_var.get()
        if sel_cal in self._cal_data:
            mags, errs = self._cal_data[sel_cal]
            data[f'Cal_{sel_cal}_mag']     = mags
            data[f'Cal_{sel_cal}_mag_err'] = errs

        export_df = pd.DataFrame(data)
        # Round all float columns to 4 decimal places
        float_cols = export_df.select_dtypes(include='float').columns
        export_df[float_cols] = export_df[float_cols].round(4)
        export_df.to_csv(out_path, index=False)
        print(f"Exported data to {out_path}")


# ────────────────────────────────── mode_dir ────────────────────────────────

def mode_dir(directory):
    directory = Path(directory)
    fits_files = sorted(glob.glob(str(directory / '*.fits')))
    ldac_files = sorted(glob.glob(str(directory / '*.ldac.db')))
    if not ldac_files:
        ldac_files = sorted(glob.glob(str(directory / '*.ldac')))

    if not fits_files or not ldac_files:
        print("No matches for fits or ldac files in the specified directory.")
        return

    dfs = []
    print(f"Reading {len(ldac_files)} LDAC files...")
    for lf in ldac_files:
        dfs.append(read_ldac(lf))

    common_ids = set.intersection(*(set(df["ident"]) for df in dfs))
    common_ids = sorted(list(common_ids))
    print(f"Found {len(common_ids)} stars present in all LDACs.")
    if len(common_ids) == 0:
        return

    df_first = dfs[0]
    df_first = df_first[df_first['ident'].isin(common_ids)]
    available_cols = [c for c in [
        'ident', 'FLUX_APER', 'FLUXERR_APER', 'MAG_APER', 'MAGERR_APER',
        'FLUX_AUTO', 'FLUXERR_AUTO', 'MAG_AUTO', 'MAGERR_AUTO', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE',
        'XWIN_IMAGE', 'YWIN_IMAGE', 'ra_deg', 'dec_deg', 'FLAGS',
        'mag', '_Bmag', '_e_Bmag', '_Vmag', '_e_Vmag',
        '_Imag', '_e_Imag', '_Rmag', '_e_Rmag',
    ] if c in df_first.columns]
    df_first = df_first[available_cols]
    df_first['ident'] = df_first['ident'].astype(int)

    fits_data = fits.open(fits_files[0])
    data = fits_data[0].data

    # ── ttkbootstrap main window ──────────────────────────────────────────
    root = ttk.Window(themename="flatly")
    root.title(f"Common Stars – {directory.name}")
    root.geometry("1050x900")

    # Top toolbar
    toolbar_frame = ttk.Frame(root, padding=(6, 4))
    toolbar_frame.pack(side=TOP, fill=X)

    def export_csv_all():
        out_path = _ask_save_file(root, "Common_stars_photometry_data.csv")
        if out_path is None:
            return
        df_first.to_csv(out_path, index=False)
        print(f"Exported common stars data to {out_path}")

    ttk.Button(toolbar_frame, text="Export All to CSV", bootstyle=PRIMARY,
               command=export_csv_all).pack(side=LEFT, padx=4)

    ttk.Button(toolbar_frame, text="Help (F1)", bootstyle=SECONDARY,
               command=lambda: show_flag_help(root)).pack(side=LEFT, padx=4)

    root.bind("<F1>", lambda _e: show_flag_help(root))

    status_var = ttk.StringVar(
        value=f"{len(common_ids)} common stars  |  Click a star to view its light-curve")
    ttk.Label(toolbar_frame, textvariable=status_var,
              bootstyle=SECONDARY).pack(side=LEFT, padx=12)

    # ── matplotlib overview figure ────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 8))
    fig.tight_layout()

    try:
        norm = ImageNormalize(data, interval=ZScaleInterval())
    except Exception:
        norm = None

    if norm is not None:
        ax.imshow(data, cmap='gray', origin='lower', norm=norm)
    else:
        ax.imshow(data, cmap='gray', origin='lower', vmin=0, vmax=6500)

    flags = df_first['FLAGS'].unique()
    for flag in flags:
        flag_data = df_first[df_first['FLAGS'] == flag]
        if 'ISOAREA_IMAGE' in flag_data.columns:
            s = flag_data['ISOAREA_IMAGE'] ** 0.5 * np.pi
        else:
            s = 15
        color = get_flag_color(flag)
        ax.scatter(flag_data['XWIN_IMAGE'], flag_data['YWIN_IMAGE'], s=s,
                   label=f'FLAG={int(flag)}', alpha=1.0,
                   edgecolors=color, facecolors='none',
                   picker=True, pickradius=5)

    ax.legend()
    ax.set_title(
        f"Common stars from {len(ldac_files)} frames  |  {Path(fits_files[0]).name}")

    # keep child windows alive
    _child_windows = []

    def on_pick(event):
        if not isinstance(event.artist,
                          plt.matplotlib.collections.PathCollection):
            return
        points = event.artist
        if hasattr(fig, 'highlight_circle') and points == fig.highlight_circle:
            return

        label = points.get_label()
        flag = float(label.split('=')[1])
        flag_data_pick = df_first[df_first['FLAGS'] == flag]
        ind = event.ind[0]
        star_id = flag_data_pick.iloc[ind]['ident']

        pos = points.get_offsets()[ind]
        if hasattr(fig, 'highlight_circle'):
            fig.highlight_circle.remove()
        fig.highlight_circle = ax.scatter(
            pos[0], pos[1], s=300, facecolors='none',
            edgecolors='red', linewidth=2, zorder=5)
        fig.canvas.draw_idle()

        star_rows = [d[d['ident'] == star_id].iloc[0] for d in dfs]

        # Update existing window if there is one and it hasn't been closed
        if _child_windows and _child_windows[0].win.winfo_exists():
            win = _child_windows[0]
            win.update_data(star_id, star_rows)
            # Bring existing window to the front
            win.win.lift()
            win.win.focus_force()
        else:
            _child_windows.clear()
            win = StarDetailWindow(root, star_id, dfs, star_rows, ldac_files)
            _child_windows.append(win)

        star_id_str = str(int(star_id)) if pd.notnull(star_id) else str(star_id)
        status_var.set(f"Selected Star ID: {star_id_str}")

    fig.canvas.mpl_connect('pick_event', on_pick)

    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=BOTH, expand=YES)

    nav_bar = ttk.Frame(root)
    nav_bar.pack(side=BOTTOM, fill=X)
    NavigationToolbar2Tk(canvas, nav_bar)

    def _on_close_dir():
        plt.close('all')
        root.destroy()

    root.protocol("WM_DELETE_WINDOW", _on_close_dir)
    root.mainloop()


# ─────────────────────────────────── entry point ────────────────────────────

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="View stars and their photometry from FITS and LDAC files.")
    parser.add_argument(
        'path', nargs='?', default='.',
        help="File or directory path. If .fits file → single FITS mode. "
             "If directory → multi-frame mode.")
    args = parser.parse_args()
    path_obj = Path(args.path)

    if path_obj.is_file() and str(path_obj).endswith('.fits'):
        mode_fits(str(path_obj))
    elif path_obj.is_dir() or str(args.path) == '.':
        mode_dir(str(path_obj))
    else:
        print("Invalid path provided. Please provide a valid .fits file or a directory.")

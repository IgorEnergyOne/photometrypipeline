# -*- coding: utf-8 -*-
"""
star_picker.py - integrated star-selector for pp_interactive_src.

Provides:
  * read_ldac_df            - read any LDAC / .ldac.db file -> DataFrame
  * collect_fits_ldac_pairs - find (fits, ldac) pairs in a directory
  * compute_common_stars    - stars present in every LDAC file (intersection)
  * apply_new_control_star  - replace control-star columns in a LightCurveData
  * reset_control_star      - undo the replacement, restoring original data
  * StarPickerDialog        - Toplevel that shows a FITS image + star overlay;
                              only stars present in ALL frames are selectable.

Logic adapted from pp_stars.py (left unchanged by this module).
"""
from __future__ import annotations

import glob
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from astropy.io import fits

import tkinter as tk
import ttkbootstrap as ttk
from ttkbootstrap.constants import *

try:
    from astropy.visualization import ZScaleInterval, ImageNormalize
    _HAS_ZSCALE = True
except ImportError:
    _HAS_ZSCALE = False


def _read_cat(ldac_path: str) -> pd.DataFrame:
    from catalog import catalog as _Catalog
    cat = _Catalog(ldac_path)
    if str(ldac_path).endswith('.db'):
        cat.read_database(filename=ldac_path)
    else:
        cat.read_ldac(filename=ldac_path)
    return cat.data.to_pandas()


# -- flag colours (identical to pp_stars.py) ----------------------------------
_FLAG_COLORS: Dict[int, str] = {
    0: 'cyan', 1: 'green', 2: 'yellow', 3: 'orange',
    4: 'purple', 8: 'blue', 16: '#c0392b',
}

def _flag_color(flag: int) -> str:
    return _FLAG_COLORS.get(int(flag), 'blue')


# -- calibrated filter columns (same list as pp_stars.py) --------------------
_CAL_FILTERS = [
    ("U",   "_Umag",  "_e_Umag"),
    ("B",   "_Bmag",  "_e_Bmag"),
    ("V",   "_Vmag",  "_e_Vmag"),
    ("R",   "_Rmag",  "_e_Rmag"),
    ("I",   "_Imag",  "_e_Imag"),
    ("u",   "_umag",  "_e_umag"),
    ("g",   "_gmag",  "_e_gmag"),
    ("r",   "_rmag",  "_e_rmag"),
    ("i",   "_imag",  "_e_imag"),
    ("z",   "_zmag",  "_e_zmag"),
    ("G",   "_Gmag",  "_e_Gmag"),
    ("G_BP","_BPmag", "_e_BPmag"),
    ("G_RP","_RPmag", "_e_RPmag"),
]

_MAG_COMBO_VALUES = [
    "Both (AUTO)",
    "Both (APER)",
    "Instrumental AUTO",
    "Instrumental APER",
    "Calibrated",
]


# -- LDAC I/O -----------------------------------------------------------------

def read_ldac_df(ldac_path: str) -> pd.DataFrame:
    """Read an LDAC or .ldac.db file and return a cleaned DataFrame."""
    df = _read_cat(ldac_path)
    if df is None:
        return pd.DataFrame()
    if 'ident' in df.columns:
        df = df[df['ident'].notna()].copy()
        df['ident'] = df['ident'].astype('Int64')
    return df.reset_index(drop=True)


def collect_fits_ldac_pairs(directory: str) -> List[Tuple[str, str]]:
    """Return sorted (fits_path, ldac_path) pairs found in *directory*."""
    pairs: List[Tuple[str, str]] = []
    for ff in sorted(glob.glob(str(Path(directory) / '*.fits'))):
        stem = Path(ff).stem
        db   = Path(directory) / (stem + '.ldac.db')
        ldac = Path(directory) / (stem + '.ldac')
        if db.exists():
            pairs.append((ff, str(db)))
        elif ldac.exists():
            pairs.append((ff, str(ldac)))
    return pairs


# -- star-matching helpers ----------------------------------------------------

def find_nearest_px(df: pd.DataFrame, x: float, y: float,
                    max_dist: float = 30.0) -> Optional[int]:
    """Return row-index of the nearest star within *max_dist* pixels."""
    if df is None or df.empty or 'XWIN_IMAGE' not in df.columns:
        return None
    dist = np.hypot(df['XWIN_IMAGE'].values - x, df['YWIN_IMAGE'].values - y)
    idx = int(np.argmin(dist))
    return idx if dist[idx] <= max_dist else None


def find_nearest_radec(df: pd.DataFrame, ra: float, dec: float,
                       tol_arcsec: float = 5.0) -> Optional[int]:
    """Return row-index of the nearest star within *tol_arcsec*."""
    if df is None or df.empty or 'ra_deg' not in df.columns:
        return None
    tol = tol_arcsec / 3600.0
    cos_dec = np.cos(np.radians(dec))
    dra  = (df['ra_deg'].values  - ra) * cos_dec
    ddec =  df['dec_deg'].values - dec
    dist = np.hypot(dra, ddec)
    idx  = int(np.argmin(dist))
    return idx if dist[idx] <= tol else None


def compute_common_stars(
    pairs: List[Tuple[str, str]],
    tol_arcsec: float = 5.0,
) -> pd.DataFrame:
    """Return a DataFrame of stars that appear in EVERY LDAC file.

    Uses the first LDAC as the reference catalogue and cross-matches every
    subsequent LDAC against it.  Only stars with a match within *tol_arcsec*
    in all other frames are kept.

    The returned DataFrame has the columns of the first LDAC plus a column
    ``n_frames`` = number of frames the star was found in (always == len(pairs)
    for every row in the result).
    """
    if not pairs:
        return pd.DataFrame()

    # Load all LDACs
    dfs: List[pd.DataFrame] = []
    for _, lp in pairs:
        try:
            dfs.append(read_ldac_df(lp))
        except Exception:
            dfs.append(pd.DataFrame())

    ref = dfs[0]
    if ref is None or ref.empty or 'ra_deg' not in ref.columns:
        return pd.DataFrame()

    # Indices in the reference frame that survive all cross-matches
    surviving = list(range(len(ref)))

    for other in dfs[1:]:
        if other is None or other.empty or 'ra_deg' not in other.columns:
            continue
        still_alive = []
        for i in surviving:
            star = ref.iloc[i]
            match = find_nearest_radec(
                other,
                float(star['ra_deg']),
                float(star['dec_deg']),
                tol_arcsec=tol_arcsec,
            )
            if match is not None:
                still_alive.append(i)
        surviving = still_alive
        if not surviving:
            break

    if not surviving:
        return pd.DataFrame()

    result = ref.iloc[surviving].copy().reset_index(drop=True)
    result['n_frames'] = len(pairs)
    return result


# -- control-star application -------------------------------------------------

_CONTROL_COLS = ('mag_control', 'sig_control',
                 'control_mag_inst', 'control_sig_inst',
                 'rel_mag', 'rel_sig')


def apply_new_control_star(
    lc,
    directory: str,
    ref_ra: float,
    ref_dec: float,
    photo_col: str = 'MAG_APER',
    photo_err_col: str = 'MAGERR_APER',
    tol_arcsec: float = 5.0,
) -> Tuple[int, int]:
    """Replace control-star columns using the star at (ref_ra, ref_dec).

    Returns (n_matched, n_total) - frames successfully cross-matched.
    """
    df = lc.df
    if df is None:
        return 0, 0

    # Back up original columns on first call only
    if not hasattr(lc, '_original_control_backup'):
        lc._original_control_backup: Dict[str, pd.Series] = {}
        for col in _CONTROL_COLS:
            lc._original_control_backup[col] = (
                df[col].copy() if col in df.columns
                else pd.Series(np.nan, index=df.index, dtype=float)
            )

    # Build ldac-filename -> ldac-path look-up
    pairs = collect_fits_ldac_pairs(directory)
    ldac_map: Dict[str, str] = {}
    for _ff, lf in pairs:
        p = Path(lf)
        stem = p.stem if not p.name.endswith('.ldac.db') else Path(p.stem).stem
        ldac_map[stem + '.ldac'] = lf
        ldac_map[p.name] = lf

    if 'photo_method' in df.columns:
        pm = str(df['photo_method'].iloc[0]).upper()
        if pm == 'AUTO':
            photo_col, photo_err_col = 'MAG_AUTO', 'MAGERR_AUTO'

    n_total   = len(df)
    new_inst  = np.full(n_total, np.nan)
    new_isig  = np.full(n_total, np.nan)
    n_matched = 0
    ldac_cache: Dict[str, pd.DataFrame] = {}

    for i, (row_i, row) in enumerate(df.iterrows()):
        csv_fname = str(row.get('filename', ''))
        ldac_path = (
            ldac_map.get(csv_fname)
            or ldac_map.get(Path(csv_fname).name)
        )
        if ldac_path is None:
            stem = Path(csv_fname).stem
            ldac_path = next((v for k, v in ldac_map.items() if stem in k), None)
        if ldac_path is None:
            continue

        if ldac_path not in ldac_cache:
            try:
                ldac_cache[ldac_path] = read_ldac_df(ldac_path)
            except Exception:
                ldac_cache[ldac_path] = pd.DataFrame()
        fdf = ldac_cache[ldac_path]
        if fdf is None or fdf.empty:
            continue

        idx = find_nearest_radec(fdf, ref_ra, ref_dec, tol_arcsec=tol_arcsec)
        if idx is None:
            continue

        star = fdf.iloc[idx]
        pc  = photo_col     if photo_col     in fdf.columns else 'MAG_AUTO'
        pce = photo_err_col if photo_err_col in fdf.columns else 'MAGERR_AUTO'

        mag = float(star.get(pc, np.nan))
        err = float(star.get(pce, 0.0))
        if np.isnan(mag):
            continue

        new_inst[i] = mag
        new_isig[i] = 0.0 if np.isnan(err) else err
        n_matched += 1

    target_inst = df.get('inst_mag',      pd.Series(np.nan, index=df.index)).to_numpy(float)
    target_isig = df.get('inst_sig',      pd.Series(0.0,   index=df.index)).to_numpy(float)
    zp          = df.get('zeropoint',     pd.Series(0.0,   index=df.index)).to_numpy(float)
    zps         = df.get('zeropoint_sig', pd.Series(0.0,   index=df.index)).to_numpy(float)

    df['control_mag_inst'] = new_inst
    df['control_sig_inst'] = new_isig
    df['mag_control']      = new_inst + zp
    df['sig_control']      = np.sqrt(new_isig**2 + zps**2)
    df['rel_mag']          = target_inst - new_inst
    df['rel_sig']          = np.sqrt(target_isig**2 + new_isig**2)

    lc._fill_arrays_cache()
    return n_matched, n_total


def reset_control_star(lc) -> bool:
    """Restore the original control-star columns. Returns True on success."""
    if not hasattr(lc, '_original_control_backup'):
        return False
    df = lc.df
    if df is None:
        return False
    for col, series in lc._original_control_backup.items():
        df[col] = series.values if len(series) == len(df) else np.nan
    del lc._original_control_backup
    lc._fill_arrays_cache()
    return True


# -- Star Photometry Window ---------------------------------------------------

class StarPhotometryWindow:
    """Popup showing per-star photometry across all frames.

    Adapted from pp_stars.StarDetailWindow.  Opens automatically when the
    user selects a star in StarPickerDialog; updates in-place if already open.
    """

    def __init__(
        self,
        parent: tk.Tk,
        star_id,
        dfs: List[pd.DataFrame],
        star_rows: List[pd.Series],
        ldac_files: List[str],
    ) -> None:
        self.dfs = dfs
        self._filenames = [Path(f).name for f in ldac_files]

        self.win = ttk.Toplevel(parent)
        self.win.resizable(True, True)

        # ── controls toolbar ──────────────────────────────────────────────
        ctrl = ttk.Frame(self.win, padding=(8, 6))
        ctrl.pack(side=TOP, fill=X)

        ttk.Label(ctrl, text="Mag Type:").pack(side=LEFT, padx=(0, 4))
        self.mag_var = ttk.StringVar(value="Both (AUTO)")
        mag_cb = ttk.Combobox(ctrl, textvariable=self.mag_var,
                              values=_MAG_COMBO_VALUES,
                              state="readonly", width=18)
        mag_cb.pack(side=LEFT, padx=(0, 12))
        mag_cb.bind("<<ComboboxSelected>>", lambda _e: self._redraw())

        ttk.Label(ctrl, text="Cal. Filter:").pack(side=LEFT, padx=(0, 4))
        self.cal_filter_var = ttk.StringVar()
        self.cal_cb = ttk.Combobox(ctrl, textvariable=self.cal_filter_var,
                                   state="readonly", width=7)
        self.cal_cb.pack(side=LEFT, padx=(0, 12))
        self.cal_cb.bind("<<ComboboxSelected>>", lambda _e: self._redraw())

        ttk.Label(ctrl, text="Show Error:").pack(side=LEFT, padx=(0, 4))
        self.err_var = ttk.StringVar(value="Yes")
        err_cb = ttk.Combobox(ctrl, textvariable=self.err_var,
                              values=["Yes", "No"],
                              state="readonly", width=6)
        err_cb.pack(side=LEFT, padx=(0, 16))
        err_cb.bind("<<ComboboxSelected>>", lambda _e: self._redraw())

        # ── matplotlib figure ─────────────────────────────────────────────
        self.fig, self.ax_left = plt.subplots(figsize=(8, 4.5))
        self.ax_right = self.ax_left.twinx()

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.win)
        self.canvas.get_tk_widget().pack(fill=BOTH, expand=YES)

        nav_bar = ttk.Frame(self.win)
        nav_bar.pack(side=BOTTOM, fill=X)
        NavigationToolbar2Tk(self.canvas, nav_bar)

        self.update_data(star_id, star_rows)

    # ── data update ───────────────────────────────────────────────────────

    def update_data(self, star_id, star_rows: List[pd.Series]) -> None:
        self.star_id   = star_id
        self.star_rows = star_rows

        star_id_str = (str(int(star_id))
                       if star_id is not None and pd.notnull(star_id)
                       else str(star_id))
        self._star_id_str = star_id_str
        self.win.title(f"Star Photometry – ID {star_id_str}")
        self.win.geometry("900x560")

        def _col(rows: List[pd.Series], col: str) -> List[float]:
            out = []
            for r in rows:
                if r is None or not hasattr(r, 'index') or len(r.index) == 0:
                    out.append(np.nan)
                elif col in r.index:
                    v = r[col]
                    out.append(float(v) if pd.notnull(v) else np.nan)
                else:
                    out.append(np.nan)
            return out

        self.inst_mag_auto = _col(star_rows, 'MAG_AUTO')
        self.inst_err_auto = _col(star_rows, 'MAGERR_AUTO')
        self.inst_mag_aper = _col(star_rows, 'MAG_APER')
        self.inst_err_aper = _col(star_rows, 'MAGERR_APER')

        # Detect which calibrated filters are actually present
        self._cal_data: dict = {}
        for label, mag_col, err_col in _CAL_FILTERS:
            mags = _col(star_rows, mag_col)
            errs = _col(star_rows, err_col)
            if any(not np.isnan(m) for m in mags):
                self._cal_data[label] = (mags, errs)

        self._cal_labels = list(self._cal_data.keys()) or ["(none)"]
        default_cal = "R" if "R" in self._cal_data else self._cal_labels[0]

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

    def _redraw(self) -> None:
        ax_l = self.ax_left
        ax_r = self.ax_right
        ax_l.cla()
        ax_r.cla()
        ax_r.set_visible(False)

        mag_type  = self.mag_var.get()
        show_err  = self.err_var.get() == "Yes"
        x         = np.arange(len(self.dfs))

        show_both = mag_type.startswith("Both")
        show_inst = show_both or mag_type.startswith("Instrumental")
        show_cal  = show_both or mag_type == "Calibrated"
        use_aper  = "APER" in mag_type

        inst_mags  = self.inst_mag_aper if use_aper else self.inst_mag_auto
        inst_errs  = self.inst_err_aper if use_aper else self.inst_err_auto
        inst_label = ("Instrumental Mag (MAG_APER)"
                      if use_aper else "Instrumental Mag (MAG_AUTO)")

        JITTER = 0.1
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
            ax_l.set_ylabel(inst_label, color='tab:blue')
            ax_l.tick_params(axis='y', labelcolor='tab:blue')

        # ── calibrated (right axis when "Both", else left) ────────────────
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
            ax_to_use.set_ylabel(f"Calibrated Magnitude ({sel_cal})",
                                 color='tab:red')
            ax_to_use.tick_params(axis='y', labelcolor='tab:red')

        ax_l.set_xlabel("Image Index")
        ax_l.set_title(f"Star ID: {self._star_id_str}")
        ax_l.set_xticks(x)

        h1, l1 = ax_l.get_legend_handles_labels()
        h2, l2 = (ax_r.get_legend_handles_labels()
                  if show_both else ([], []))
        if h1 or h2:
            ax_l.legend(h1 + h2, l1 + l2, loc='upper left')

        self.fig.tight_layout()
        self.canvas.draw_idle()


# -- Star Picker Dialog -------------------------------------------------------

class StarPickerDialog:
    """Toplevel showing a FITS image with SExtractor stars overlaid.

    Only stars that are present in EVERY LDAC file in the directory are
    shown as selectable (coloured circles).  Stars absent from one or more
    frames are drawn as small grey crosses for context only.

    Confirm the selection with the "Use as Control Star" button or Enter.
    """

    def __init__(
        self,
        parent: tk.Tk,
        directory: str,
        callback: Callable[[float, float, str, str], None],
        tol_arcsec: float = 5.0,
    ) -> None:
        self._parent    = parent
        self._directory = directory
        self._callback  = callback

        self._pairs = collect_fits_ldac_pairs(directory)
        if not self._pairs:
            from tkinter import messagebox
            messagebox.showerror(
                "No data",
                f"No FITS/LDAC pairs found in:\n{directory}",
                parent=parent,
            )
            return

        self._frame_idx  = 0
        self._selected: Optional[Dict] = None
        self._highlight  = None
        self._ldac_cache: Dict[str, pd.DataFrame] = {}
        self._current_df: Optional[pd.DataFrame]  = None  # full frame stars
        self._common_df: Optional[pd.DataFrame]   = None  # intersection
        self._phot_window: Optional[StarPhotometryWindow] = None  # photometry popup

        self._build()

        # Compute the cross-frame star intersection (may take a moment)
        self._compute_common(tol_arcsec)

        self._load_frame(0)

    # -- UI construction --------------------------------------------------

    def _build(self) -> None:
        self._top = ttk.Toplevel(self._parent)
        self._top.title("Select Control Star")
        self._top.transient(self._parent)
        self._top.resizable(True, True)
        self._top.geometry(
            f"960x760+{self._parent.winfo_x()+40}+{self._parent.winfo_y()+40}")

        # ---- TOP toolbar ------------------------------------------------
        tb = ttk.Frame(self._top, padding=(6, 4))
        tb.pack(side=TOP, fill=X)

        ttk.Button(tb, text="<< Prev", width=7,
                   command=lambda: self._nav(-1)).pack(side=LEFT)
        self._frame_var = ttk.StringVar(value="Loading...")
        ttk.Label(tb, textvariable=self._frame_var, width=38).pack(side=LEFT, padx=4)
        ttk.Button(tb, text="Next >>", width=7,
                   command=lambda: self._nav(+1)).pack(side=LEFT)

        ttk.Separator(tb, orient=VERTICAL).pack(side=LEFT, fill=Y, padx=8)
        ttk.Label(tb, text="Phot:").pack(side=LEFT)
        self._photo_var = ttk.StringVar(value="APER")
        ttk.Combobox(tb, textvariable=self._photo_var,
                     values=["APER", "AUTO"], state="readonly",
                     width=6).pack(side=LEFT, padx=4)

        ttk.Separator(tb, orient=VERTICAL).pack(side=LEFT, fill=Y, padx=8)
        ttk.Label(tb, text="Match radius (arcsec):").pack(side=LEFT)
        self._tol_var = ttk.StringVar(value="5.0")
        ttk.Entry(tb, textvariable=self._tol_var, width=5).pack(side=LEFT, padx=4)

        # ---- BOTTOM status bar (packed BEFORE canvas so it is never hidden)
        status_frm = ttk.Frame(self._top, padding=(6, 4))
        status_frm.pack(side=BOTTOM, fill=X)

        self._use_btn = ttk.Button(
            status_frm,
            text="Use as Control Star",
            style="Accent.TButton",
            state=DISABLED,
            command=self._on_use,
        )
        self._use_btn.pack(side=RIGHT, padx=4)
        ttk.Button(status_frm, text="Close",
                   command=self._top.destroy).pack(side=RIGHT, padx=4)
        ttk.Button(status_frm, text="Cancel",
                   command=self._top.destroy).pack(side=RIGHT, padx=4)

        self._info_var = ttk.StringVar(
            value="Computing star intersection across all frames...")
        ttk.Label(status_frm, textvariable=self._info_var,
                  wraplength=600, anchor=W).pack(side=LEFT, fill=X, expand=True)

        # ---- Matplotlib canvas (packed LAST so it fills all remaining space)
        plot_frm = ttk.Frame(self._top)
        plot_frm.pack(side=TOP, fill=BOTH, expand=True)

        self._fig, self._ax = plt.subplots(figsize=(9, 6.5))
        self._fig.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01)

        self._canvas = FigureCanvasTkAgg(self._fig, master=plot_frm)
        self._canvas.get_tk_widget().pack(fill=BOTH, expand=True)
        NavigationToolbar2Tk(self._canvas, plot_frm).update()

        self._canvas.mpl_connect("button_press_event", self._on_click)

        # Enter key confirms selection
        self._top.bind("<Return>", lambda _e: self._on_use())
        self._top.bind("<Escape>", lambda _e: self._top.destroy())

    # -- common-star intersection -----------------------------------------

    def _compute_common(self, tol_arcsec: float) -> None:
        """Load all LDACs and find the intersection of stars."""
        self._top.config(cursor="watch")
        self._top.update_idletasks()
        try:
            self._common_df = compute_common_stars(self._pairs, tol_arcsec=tol_arcsec)
            n = len(self._common_df) if self._common_df is not None else 0
            self._info_var.set(
                f"{n} stars present in all {len(self._pairs)} frame(s).  "
                f"Click a coloured circle to select."
            )
        except Exception as exc:
            self._info_var.set(f"Error computing common stars: {exc}")
            self._common_df = pd.DataFrame()
        finally:
            self._top.config(cursor="")

    # -- frame navigation -------------------------------------------------

    def _nav(self, delta: int) -> None:
        self._load_frame((self._frame_idx + delta) % len(self._pairs))

    def _load_frame(self, idx: int) -> None:
        self._frame_idx = idx
        fits_path, ldac_path = self._pairs[idx]
        self._frame_var.set(
            f"[{idx+1}/{len(self._pairs)}]  {Path(ldac_path).name}")

        if ldac_path not in self._ldac_cache:
            try:
                self._ldac_cache[ldac_path] = read_ldac_df(ldac_path)
            except Exception as exc:
                self._info_var.set(f"Error reading LDAC: {exc}")
                return
        self._current_df = self._ldac_cache[ldac_path]

        try:
            with fits.open(fits_path) as hdu:
                self._current_image = hdu[0].data
        except Exception as exc:
            self._info_var.set(f"Error reading FITS: {exc}")
            return

        self._redraw()

    # -- drawing ----------------------------------------------------------

    def _redraw(self) -> None:
        self._ax.cla()
        data = self._current_image
        all_df   = self._current_df
        common   = self._common_df

        # FITS image
        try:
            if _HAS_ZSCALE:
                norm = ImageNormalize(data, interval=ZScaleInterval())
                self._ax.imshow(data, cmap='gray', origin='lower', norm=norm)
            else:
                self._ax.imshow(data, cmap='gray', origin='lower')
        except Exception:
            self._ax.imshow(data, cmap='gray', origin='lower')

        # All detected stars: small grey crosses (context only, not selectable)
        if all_df is not None and not all_df.empty and 'XWIN_IMAGE' in all_df.columns:
            self._ax.scatter(
                all_df['XWIN_IMAGE'], all_df['YWIN_IMAGE'],
                s=6, marker='+', c='gray', linewidths=0.5,
                alpha=0.5, zorder=2, label='Not in all frames',
            )

        # Common stars: coloured by flag, selectable
        if common is not None and not common.empty and 'XWIN_IMAGE' in common.columns:
            flag_groups = (common.groupby('FLAGS')
                           if 'FLAGS' in common.columns else [(0, common)])
            for flag, gdf in flag_groups:
                sz = (np.sqrt(gdf['ISOAREA_IMAGE'].values) * np.pi
                      if 'ISOAREA_IMAGE' in gdf.columns
                      else np.full(len(gdf), 18))
                sz = np.clip(sz, 12, 60)
                self._ax.scatter(
                    gdf['XWIN_IMAGE'], gdf['YWIN_IMAGE'],
                    s=sz, edgecolors=_flag_color(flag),
                    facecolors='none', linewidths=1.2,
                    zorder=4, label=f'All frames, FLAG={int(flag)}',
                )
            self._ax.legend(loc='lower right', fontsize='x-small', framealpha=0.6)

        # Re-draw selection highlight
        if self._selected is not None:
            self._highlight = self._ax.scatter(
                [self._selected['x']], [self._selected['y']],
                s=400, facecolors='none',
                edgecolors='red', linewidths=2.5, zorder=10,
            )

        self._ax.set_axis_off()
        self._canvas.draw_idle()

    # -- click handler ----------------------------------------------------

    def _on_click(self, event) -> None:
        if event.inaxes != self._ax or event.button != 1:
            return
        try:
            if self._canvas.toolbar.mode != '':
                return   # zoom / pan active
        except Exception:
            pass

        # Search only within the common-star set
        common = self._common_df
        if common is None or common.empty:
            self._info_var.set(
                "No common stars found. Check that the LDAC files match the "
                "CSV filenames and that the match radius is appropriate.")
            return

        idx = find_nearest_px(common, event.xdata, event.ydata, max_dist=35)
        if idx is None:
            self._info_var.set(
                "No common star near the clicked position.  "
                "Only coloured circles are selectable (grey crosses are "
                "stars not present in all frames).")
            return

        star = common.iloc[idx]
        pm  = self._photo_var.get()
        pc  = 'MAG_APER'    if pm == 'APER' else 'MAG_AUTO'
        pce = 'MAGERR_APER' if pm == 'APER' else 'MAGERR_AUTO'
        pc  = pc  if pc  in common.columns else 'MAG_AUTO'
        pce = pce if pce in common.columns else 'MAGERR_AUTO'

        mag   = float(star.get(pc,  np.nan))
        err   = float(star.get(pce, np.nan))
        flags = int(star.get('FLAGS', 0)) if 'FLAGS' in star.index else 0
        x     = float(star.get('XWIN_IMAGE', event.xdata))
        y     = float(star.get('YWIN_IMAGE', event.ydata))
        ra    = float(star.get('ra_deg',  0.0))
        dec   = float(star.get('dec_deg', 0.0))

        if ra == 0.0 and dec == 0.0:
            self._info_var.set(
                "Selected star has no RA/Dec in this LDAC - cannot cross-match.")
            return

        self._selected = dict(ra=ra, dec=dec, x=x, y=y,
                              mag=mag, err=err, flags=flags, pc=pc, pce=pce)

        if self._highlight is not None:
            try:
                self._highlight.remove()
            except Exception:
                pass
        self._highlight = self._ax.scatter(
            [x], [y], s=400, facecolors='none',
            edgecolors='red', linewidths=2.5, zorder=10,
        )
        self._canvas.draw_idle()

        ident = star.get('ident', None)
        id_str = (str(int(ident))
                  if ident is not None and pd.notnull(ident) else 'N/A')
        self._info_var.set(
            f"[OK] Selected  ident={id_str}  "
            f"X={x:.1f}  Y={y:.1f}  "
            f"RA={ra:.5f}  Dec={dec:.5f}  "
            f"{pc}={mag:.3f} +/- {err:.3f}  FLAGS={flags}  "
            f"-- Press Enter or click 'Use as Control Star' to confirm."
        )
        self._use_btn.configure(state=NORMAL)
        self._use_btn.focus_set()   # focus so Enter key works immediately

        # Open / update the photometry window automatically
        self._show_photometry(star)

    # -- photometry window ------------------------------------------------

    def _show_photometry(self, star: pd.Series) -> None:
        """Open or update StarPhotometryWindow for the given common-star row."""
        ra  = float(star.get('ra_deg',  0.0))
        dec = float(star.get('dec_deg', 0.0))
        if ra == 0.0 and dec == 0.0:
            return

        try:
            tol = float(self._tol_var.get())
        except ValueError:
            tol = 5.0

        # Ensure all LDACs are loaded (may not all be cached yet)
        dfs: List[pd.DataFrame] = []
        star_rows: List[pd.Series] = []
        ldac_files: List[str] = []

        for _, ldac_path in self._pairs:
            if ldac_path not in self._ldac_cache:
                try:
                    self._ldac_cache[ldac_path] = read_ldac_df(ldac_path)
                except Exception:
                    self._ldac_cache[ldac_path] = pd.DataFrame()
            df = self._ldac_cache[ldac_path]
            dfs.append(df)
            ldac_files.append(ldac_path)

            if df is not None and not df.empty:
                idx = find_nearest_radec(df, ra, dec, tol_arcsec=tol)
                star_rows.append(df.iloc[idx] if idx is not None
                                 else pd.Series(dtype=float))
            else:
                star_rows.append(pd.Series(dtype=float))

        ident = star.get('ident', None)

        # Reuse existing window if still open, otherwise create a new one
        if self._phot_window is not None:
            try:
                if self._phot_window.win.winfo_exists():
                    self._phot_window.dfs = dfs
                    self._phot_window._filenames = [Path(f).name
                                                    for f in ldac_files]
                    self._phot_window.update_data(ident, star_rows)
                    self._phot_window.win.lift()
                    return
            except Exception:
                pass
        self._phot_window = StarPhotometryWindow(
            self._top, ident, dfs, star_rows, ldac_files)

    # -- confirm ----------------------------------------------------------

    def _on_use(self) -> None:
        if self._selected is None:
            return
        s = self._selected
        try:
            tol = float(self._tol_var.get())
        except ValueError:
            tol = 5.0
        self._callback(s['ra'], s['dec'], s['pc'], s['pce'])
        # Update status bar – do NOT close the window so the user can
        # keep trying different control stars.
        self._info_var.set(
            f"✓ Applied control star  RA={s['ra']:.5f}  Dec={s['dec']:.5f}  "
            f"{s['pc']}={s['mag']:.3f} ± {s['err']:.3f}  "
            f"-- Select another star or click Close."
        )


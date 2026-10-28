# -*- coding: utf-8 -*-
"""
Plotting layer for pp_interactive_src: BlitManager and LightCurvePlot.
"""
from __future__ import annotations

from typing import TYPE_CHECKING, Dict, Iterable, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from astropy.time import Time
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import messagebox

from .constants import DEFAULT_PHASE_MAX, OFFSET_ALL_KEY
from .models import LightCurveData
from .plot_settings import PlotSettings
from .utils import _combine, _safe, debug_print, format_asteroid_name

if TYPE_CHECKING:
    from .gui import LightCurveGUI


class BlitManager:
    """Minimal blitting helper for fast, flicker-free Matplotlib updates.

    Captures an axes background snapshot after every full draw and restores it
    before re-drawing only the animated artists, avoiding a full canvas redraw
    on every selection-marker or offset change.
    """

    def __init__(self, canvas: FigureCanvasTkAgg, ax: matplotlib.axes.Axes):
        self.canvas = canvas
        self.ax = ax
        self._bg = None
        self._artists: Set[matplotlib.artist.Artist] = set()
        self.cid_draw = canvas.mpl_connect("draw_event", self._on_draw)

    def add_artists(self, artists: Iterable[matplotlib.artist.Artist]):
        """Register *artists* as animated so they are included in quick redraws."""
        for a in artists:
            if a is None:
                continue
            try:
                a.set_animated(True)
            except Exception:
                pass
            self._artists.add(a)

    def clear(self):
        """Remove all tracked artists and invalidate the cached background."""
        self._artists.clear()
        self._bg = None

    def _on_draw(self, _evt):
        """Cache the axes background immediately after a full canvas draw."""
        self._bg = self.canvas.copy_from_bbox(self.ax.bbox)

    def draw(self):
        """Request a full canvas redraw (queued, non-blocking)."""
        self.canvas.draw_idle()

    def quick_redraw(self):
        """Restore the cached background and re-draw only the animated artists.

        Falls back to a full ``canvas.draw()`` on the first call (before any
        background has been captured).
        """
        if self._bg is None:
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
    """Renders one or more LightCurveData objects onto a Matplotlib axes.

    Owns the axes, the BlitManager, and all plotted artist handles.  Keeps a
    *layout signature* string; when the signature changes (different visible
    aliases, mode, time axis, offsets, ...) the axes are fully rebuilt.
    For changes that only move existing artists (offset nudge, selection move)
    a fast blit-only path is used instead.
    """

    def __init__(self, fig: matplotlib.figure.Figure, ax: matplotlib.axes.Axes, canvas: FigureCanvasTkAgg,
                 cfg=None):  # cfg: Optional[AppConfig]
        self.fig = fig
        self.ax = ax
        self.canvas = canvas
        self._cfg = cfg  # kept for style lookups inside update() / _clear_axes()

        self.xlabel = "Julian Date"
        self.ylabel = "Magnitude"
        self.title = "Lightcurve"
        self.auto_title = True

        self._y_limits_per_mode: Dict[str, Optional[Tuple[float, float]]] = {}
        # (brighter, fainter) offsets from the curve centre, alternative to fixed limits
        self._y_offsets_per_mode: Dict[str, Tuple[float, float]] = {}
        self._y_nticks_per_mode: Dict[str, Optional[int]] = {}
        self._current_mode: str = 'target'
        self._user_xlabel: Optional[str] = None
        self._user_ylabel: Optional[str] = None
        self.user_xlim: Optional[Tuple[float, float]] = None
        self.user_ylim: Optional[Tuple[float, float]] = None
        self.asteroid_zoom_level = 1.0
        self._current_asteroid_image_original = None
        self._current_asteroid_overlay_original = None

        # Style values - initialised from config, may be overridden at runtime
        _p = cfg.plot if cfg is not None else None
        self.marker_size          = _p.marker_size          if _p else 4.0
        self.marker_style         = _p.marker_style         if _p else 'o'
        self.errorbar_capsize     = _p.errorbar_capsize     if _p else 2.0
        self.errorbar_capthick    = _p.errorbar_capthick    if _p else 1.0
        self.errorbar_linewidth   = _p.errorbar_linewidth   if _p else 1.0

        self.available_markers = PlotSettings.MARKERS
        self.marker_dict = {name: marker for marker, name in self.available_markers}

        self.valid_color = "blue"
        self.rejected_color = "red"
        self.flagged_color = "orange"

        self.parent_gui: Optional["LightCurveGUI"] = None

        self.plotted_handles: Dict[str, List[dict]] = {}
        self._layout_signature: Optional[str] = None

        _label_fs = _p.label_fontsize if _p else 9

        (self.sel_artist,) = self.ax.plot([], [], PlotSettings.SELECTION_MARKER['form'],
                                          color=PlotSettings.SELECTION_MARKER['color'],
                                          markersize=PlotSettings.SELECTION_MARKER['markersize'],
                                          zorder=PlotSettings.SELECTION_MARKER['zorder'],
                                          visible=False, animated=True)
        (self.multi_sel_artist,) = self.ax.plot([], [], 's',
                                                color='#00d4ff',
                                                markersize=PlotSettings.SELECTION_MARKER['markersize'] + 3,
                                                markeredgewidth=2.0,
                                                markerfacecolor='none',
                                                zorder=PlotSettings.SELECTION_MARKER['zorder'] - 1,
                                                visible=False, animated=True)
        try:
            self.sel_text = self.ax.text(0, 0, '', color=PlotSettings.SELECTION_MARKER['color'],
                                         fontsize=_label_fs,
                                         zorder=PlotSettings.SELECTION_MARKER['zorder'] + 1,
                                         visible=False, animated=True)
        except Exception:
            self.sel_text = None

        self.blit = BlitManager(self.canvas, self.ax)
        self._register_blit_artists()

    #  -  -  -  internal helpers  -  -  - 

    def _register_blit_artists(self):
        """Collect all current plot artists into the BlitManager for quick redraws."""
        artists: List[matplotlib.artist.Artist] = [self.sel_artist, self.multi_sel_artist]
        if getattr(self, 'sel_text', None) is not None:
            artists.append(self.sel_text)

        def _flatten(obj):
            if obj is None:
                return []
            if isinstance(obj, (list, tuple)):
                flat = []
                for item in obj:
                    flat.extend(_flatten(item))
                return flat
            if hasattr(obj, 'lines'):
                flat = []
                for item in getattr(obj, 'lines', ()):
                    flat.extend(_flatten(item))
                if hasattr(obj, 'marker_line') and obj.marker_line is not None:
                    flat.append(obj.marker_line)
                return flat
            return [obj]

        for entries in self.plotted_handles.values():
            for e in entries:
                for key in ('valid', 'flagged', 'rejected', 'ghost_valid', 'ghost_flagged', 'ghost_rejected'):
                    h = e.get(key)
                    if h is not None:
                        artists.extend(_flatten(h))
        self.blit.add_artists(artists)

    def _clear_axes(self):
        """Clear the axes and all tracked handles, ready for a full rebuild."""
        self.ax.cla()
        grid_on = True
        try:
            grid_on = bool(self.parent_gui.show_grid)
        except Exception:
            grid_on = True
        _p = self._cfg.plot if self._cfg is not None else None
        _grid_alpha  = _p.grid_alpha   if _p else 0.15
        _label_fs    = _p.label_fontsize if _p else 9
        if grid_on:
            self.ax.grid(alpha=_grid_alpha)
        else:
            self.ax.grid(False)
        self.plotted_handles.clear()
        self.blit.clear()
        (self.sel_artist,) = self.ax.plot(
            [], [], PlotSettings.SELECTION_MARKER.get('form', 'o'),
            color=PlotSettings.SELECTION_MARKER.get('color', 'red'),
            markersize=PlotSettings.SELECTION_MARKER.get('markersize', 5),
            zorder=PlotSettings.SELECTION_MARKER.get('zorder', 20),
            visible=False, animated=True,
        )
        (self.multi_sel_artist,) = self.ax.plot(
            [], [], 's',
            color='#00d4ff',
            markersize=PlotSettings.SELECTION_MARKER.get('markersize', 5) + 3,
            markeredgewidth=2.0,
            markerfacecolor='none',
            zorder=PlotSettings.SELECTION_MARKER.get('zorder', 20) - 1,
            visible=False, animated=True,
        )
        try:
            self.sel_text = self.ax.text(
                0, 0, '', color=PlotSettings.SELECTION_MARKER.get('color', 'red'),
                fontsize=_label_fs,
                zorder=PlotSettings.SELECTION_MARKER.get('zorder', 20) + 1,
                visible=False, animated=True,
            )
        except Exception:
            self.sel_text = None
        self.xlabel_text = None
        self.ylabel_text = None
        self.title_text = None

    def curve_center(self, offsets: Dict[str, Dict[str, float]]) -> Optional[float]:
        """Centre of the plotted non-rejected points (offsets included), or
        ``None`` if nothing is plotted.

        Taken as the midpoint of the 5th-95th percentile range so that single
        outliers do not pull it away from the lightcurve itself."""
        values = []
        for alias, entries in self.plotted_handles.items():
            for e in entries:
                b = e['band']
                off = offsets.get(alias, {}).get(str(b) if b is not None else OFFSET_ALL_KEY,
                                                 offsets.get(alias, {}).get(OFFSET_ALL_KEY, 0.0))
                y = np.asarray(e['Ybase'], dtype=float)[~e['mask_rejected']] + off
                values.append(y[np.isfinite(y)])
        y = np.concatenate(values) if values else np.empty(0)
        if not y.size:
            return None
        lo, hi = np.percentile(y, [5, 95])
        return float(lo + hi) / 2.0

    def _set_labels(self, offsets: Optional[Dict[str, Dict[str, float]]] = None):
        """Apply axis labels/title, enforce magnitude-axis inversion, and honour
        the user-defined Y limits (fixed, or as offsets from the curve centre)
        and tick count."""
        # Prefer explicit user overrides; fall back to auto-generated values
        effective_xlabel = self._user_xlabel if self._user_xlabel is not None else self.xlabel
        effective_ylabel = self._user_ylabel if self._user_ylabel is not None else self.ylabel
        try:
            self.xlabel_text = self.ax.set_xlabel(effective_xlabel)
            self.ylabel_text = self.ax.set_ylabel(effective_ylabel)
            self.title_text = self.ax.set_title(self.title)
            for t in (self.xlabel_text, self.ylabel_text, self.title_text):
                try:
                    t.set_picker(5)
                except Exception:
                    pass
        except Exception:
            if self.ax.get_xlabel() != effective_xlabel:
                self.ax.set_xlabel(effective_xlabel)
            if self.ax.get_ylabel() != effective_ylabel:
                self.ax.set_ylabel(effective_ylabel)
            if self.ax.get_title() != self.title:
                self.ax.set_title(self.title)

        lo, hi = self.ax.get_ylim()

        requested = self.y_limits
        if requested is None and self.y_offsets is not None and offsets is not None:
            center = self.curve_center(offsets)
            if center is not None:
                brighter, fainter = self.y_offsets
                requested = (center - brighter, center + fainter)

        if requested is not None:
            req_lo, req_hi = requested
            self.ax.set_ylim(req_lo, req_hi)
            lo, hi = self.ax.get_ylim()

        if self.y_nticks is not None and self.y_nticks > 1:
            from matplotlib.ticker import MaxNLocator
            self.ax.yaxis.set_major_locator(MaxNLocator(nbins=self.y_nticks))
        else:
            from matplotlib.ticker import AutoLocator
            self.ax.yaxis.set_major_locator(AutoLocator())

        if requested is None:
            # Auto-scaling: never show less than ``min_y_range`` magnitudes
            min_range = self._cfg.plot.min_y_range if self._cfg is not None else 0.2
            extent = abs(hi - lo)
            if extent < min_range:
                mid = (hi + lo) / 2.0
                new_half = min_range / 2.0
                if lo > hi:
                    self.ax.set_ylim(mid + new_half, mid - new_half)
                else:
                    self.ax.set_ylim(mid - new_half, mid + new_half)
                lo, hi = self.ax.get_ylim()

        if lo < hi:
            self.ax.set_ylim(hi, lo)

    #  -  -  -  per-mode Y-limit properties  -  -  - 

    @property
    def y_limits(self) -> Optional[Tuple[float, float]]:
        """Y-axis limits for the currently active mode, or ``None`` for auto-scaling."""
        return self._y_limits_per_mode.get(self._current_mode)

    @property
    def y_offsets(self) -> Optional[Tuple[float, float]]:
        """(brighter, fainter) Y offsets from the curve centre for the current
        mode, or ``None`` when not in offset mode."""
        return self._y_offsets_per_mode.get(self._current_mode)

    @property
    def y_nticks(self) -> Optional[int]:
        """Y-axis tick count for the currently active mode, or ``None`` for auto."""
        return self._y_nticks_per_mode.get(self._current_mode)

    def set_y_limits(self, ymin: Optional[float], ymax: Optional[float], nticks: Optional[int] = None):
        """Set manual Y-axis bounds and optional tick count for the **current mode**.
        Pass ``None`` for both bounds to revert to automatic scaling for this mode."""
        mode = self._current_mode
        self._y_offsets_per_mode.pop(mode, None)
        if ymin is None or ymax is None:
            self._y_limits_per_mode.pop(mode, None)
        else:
            self._y_limits_per_mode[mode] = (ymin, ymax)
        self._set_y_nticks(nticks)

    def set_y_offsets(self, brighter: float, fainter: float, nticks: Optional[int] = None):
        """Show the **current mode** as the curve centre minus *brighter* to the
        centre plus *fainter* magnitudes; the centre follows the data on redraw."""
        mode = self._current_mode
        self._y_limits_per_mode.pop(mode, None)
        self._y_offsets_per_mode[mode] = (abs(brighter), abs(fainter))
        self._set_y_nticks(nticks)

    def reset_y_ranges(self):
        """Drop the manual Y ranges (limits and centre offsets) of **all** modes,
        e.g. when a new document is opened."""
        self._y_limits_per_mode.clear()
        self._y_offsets_per_mode.clear()

    def _set_y_nticks(self, nticks: Optional[int]):
        mode = self._current_mode
        if nticks is None:
            self._y_nticks_per_mode.pop(mode, None)
        else:
            self._y_nticks_per_mode[mode] = nticks

    def _choose_color_for_band(self, b: Optional[str], idx: int) -> str:
        """Return the display colour for *b*, delegating to PlotSettings when
        a parent GUI is available."""
        if self.parent_gui and hasattr(self.parent_gui, 'plot_settings'):
            return self.parent_gui.plot_settings.get_color(b, idx)
        return PlotSettings.DEFAULT_COLORS[idx % len(PlotSettings.DEFAULT_COLORS)]

    def _choose_marker_for_batch(self, batch_val, batches: List) -> str:
        """Return the marker shape for a specific batch identifier."""
        available = [m[0] for m in PlotSettings.MARKERS]
        try:
            idx = batches.index(batch_val)
        except Exception:
            try:
                idx = int(batch_val) - 1
            except Exception:
                idx = 0
        return available[idx % len(available)]

    def _compute_time_x(self, arr_jd: np.ndarray, time_mode: str, period_hours: Optional[float]) -> np.ndarray:
        """Convert raw JD values to the X-axis representation selected by *time_mode*.

        Supported modes: ``'julian_date'``, ``'mjd'``, ``'minutes'``,
        ``'rotation_phase'``.  Also updates ``self.xlabel`` as a side effect.
        """
        if arr_jd is None or arr_jd.size == 0:
            return np.array([])
        if time_mode == 'julian_date':
            self.xlabel = "Julian Date"
            return arr_jd
        if time_mode == 'mjd':
            self.xlabel = "Modified Julian Date (MJD)"
            return arr_jd - 2400000.5
        if time_mode == 'minutes':
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
            return ((arr_jd - jd0) * 24.0 / period_hours)
        return arr_jd

    @staticmethod
    def _ghost_color(color: str, amount: float = 0.45) -> str:
        """Return a lightened version of *color* for ghost phase copies."""
        try:
            import matplotlib.colors as mcolors
            r, g, b, *_ = mcolors.to_rgba(color)
        except Exception:
            return color
        r2 = r + (1.0 - r) * amount
        g2 = g + (1.0 - g) * amount
        b2 = b + (1.0 - b) * amount
        return '#{:02x}{:02x}{:02x}'.format(int(r2 * 255), int(g2 * 255), int(b2 * 255))

    # Keep old name as alias so any other call sites still work
    _lighten_color = _ghost_color

    def _layout_sig(self, lc_dict, visible, mode, time_mode, errorbar_type, selected_band, show_rejected, offsets):
        """Signature that uniquely identifies the layout state (forces rebuild when it changes)."""
        visible_aliases = ','.join(sorted([a for a, v in visible.items() if v and a in lc_dict]))
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
                try:
                    pm_val = float(self.parent_gui.phase_max_var.get())
                except Exception:
                    pm_val = DEFAULT_PHASE_MAX
                rp = f"|P={rp_val:.6f}|PM={pm_val:.4f}"
        except Exception:
            rp = "|P=?"
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
        try:
            use_rm_flag = int(bool(getattr(self.parent_gui.use_reduced_mag_var, 'get', lambda: False)()))
        except Exception:
            use_rm_flag = 0
        try:
            use_lt_flag = int(bool(getattr(self.parent_gui.use_lighttime_var, 'get', lambda: False)()))
        except Exception:
            use_lt_flag = 0
        try:
            edit_mode_flag = int(bool(getattr(self.parent_gui, 'edit_mode', False)))
        except Exception:
            edit_mode_flag = 0
        try:
            batch_markers_flag = int(bool(getattr(self.parent_gui, 'show_batch_markers', False)))
        except Exception:
            batch_markers_flag = 0

        signature = (f"{visible_aliases}|{mode}|{time_mode}{rp}|{errorbar_type}"
                     f"|{selected_band}|rej={int(bool(show_rejected))}|legend={legend_flag}|grid={grid_flag}"
                     f"|colors={colors_flag}|loc={colors_loc}"
                     f"|rm={use_rm_flag}|lt={use_lt_flag}|edit={edit_mode_flag}|bm={batch_markers_flag}"
                     f"|{self.marker_style}|{self.marker_size}|{self.valid_color}|{offsets_data}")
        return signature

    def invalidate_layout(self) -> None:
        """Force a full rebuild on next update (e.g. when masks change)."""
        self._layout_signature = None

    #  -  -  -  main update  -  -  - 

    def update(
        self,
        lc_dict: Dict[str, LightCurveData],
        offsets: Dict[str, Dict[str, float]],
        visible: Dict[str, bool],
        mode: str,
        time_mode: str,
        show_rejected: bool,
        errorbar_type: str,
        selected_bands: Set[str],
    ):
        # Track the active mode so y_limits / y_nticks properties resolve correctly
        self._current_mode = mode

        if not lc_dict:
            self._clear_axes()
            self._set_labels()
            self.canvas.draw()
            return

        new_sig = self._layout_sig(lc_dict, visible, mode, time_mode, errorbar_type,
                                   selected_bands, show_rejected, offsets)
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

                # -- Lighttime-corrected JD (if enabled) ---------------------
                use_lt = False
                use_rm = False
                if self.parent_gui is not None:
                    try:
                        use_lt = bool(self.parent_gui.use_lighttime_var.get())
                    except Exception:
                        pass
                    try:
                        use_rm = bool(self.parent_gui.use_reduced_mag_var.get())
                    except Exception:
                        pass

                if use_lt and lc.df is not None and 'corrected_jd' in lc.df.columns:
                    jd = lc.df['corrected_jd'].to_numpy(dtype=float, copy=False)
                else:
                    jd = arr.get('jd')

                # -- Y array selection by mode --------------------------------
                if mode == 'target':
                    if use_rm and lc.df is not None and 'reduced_mag' in lc.df.columns:
                        Yfull = lc.df['reduced_mag'].to_numpy(dtype=float, copy=False)
                    else:
                        Yfull = arr.get('mag')
                    if Yfull is None or not Yfull.size:
                        messagebox.showwarning("No valid target magnitudes",
                                               "No valid 'mag' values to plot.")
                        continue
                elif mode == 'instrumental':
                    Yfull = arr.get('inst_mag')
                    if Yfull is None or not Yfull.size:
                        messagebox.showwarning("No valid target magnitudes",
                                               "No valid instrumental (inst_mag) values to plot.")
                        continue
                elif mode == 'relative':
                    Yfull = arr.get('rel_mag')
                    if Yfull is None or not Yfull.size:
                        messagebox.showwarning("No valid relative magnitudes",
                                               "No valid relative (rel_mag) values to plot.")
                else:  # control
                    Yfull = arr.get('mag_control')
                    if Yfull is None or not Yfull.size:
                        messagebox.showwarning("No valid magnitudes for control star",
                                               "No valid data (mag_control) for control star to plot.")

                # -- Y error array by errorbar type ---------------------------
                if errorbar_type == 'calibrated':
                    Yerr_full = arr.get('sig')
                elif errorbar_type == 'instrumental':
                    Yerr_full = arr.get('inst_sig')
                elif errorbar_type == 'relative':
                    Yerr_full = arr.get('rel_sig')
                else:
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
                    bands_to_draw = lc_bands
                else:
                    bands_to_draw = [b for b in lc_bands
                                     if (b in selected_bands) or (b is None and 'None' in selected_bands)]

                period = None
                phase_max = DEFAULT_PHASE_MAX
                if time_mode == 'rotation_phase' and self.parent_gui is not None:
                    try:
                        period = float(self.parent_gui.rotation_period_var.get())
                    except Exception:
                        period = None
                    try:
                        phase_max = float(self.parent_gui.phase_max_var.get())
                        if not (1.0 < phase_max <= 3.0):
                            phase_max = DEFAULT_PHASE_MAX
                    except Exception:
                        phase_max = DEFAULT_PHASE_MAX

                show_batch_markers = bool(getattr(self.parent_gui, 'show_batch_markers', False))
                use_batch_markers = show_batch_markers and lc.has_batches()
                batch_arr_full = lc.arr.get('batch_number')
                batches_list = lc.get_batches() if use_batch_markers else []

                for bidx, b in enumerate(bands_to_draw):
                    band_mask = np.ones_like(jd, dtype=bool)
                    if b is not None and band_col is not None and len(band_col) == len(jd):
                        band_mask = (band_col.astype(str) == str(b))
                    if not np.any(band_mask):
                        continue

                    batch_col = batch_arr_full[band_mask] if (use_batch_markers and batch_arr_full is not None) else None

                    X_raw = self._compute_time_x(jd[band_mask], time_mode, period)

                    if time_mode == 'rotation_phase':
                        X = np.mod(X_raw, 1.0)
                        is_cycle1 = X_raw < 1.0
                    else:
                        X = X_raw
                        is_cycle1 = np.ones(len(X_raw), dtype=bool)

                    Ybase = Yfull[band_mask]
                    Yerr = (Yerr_full[band_mask]
                            if (Yerr_full is not None and len(Yerr_full) == len(jd))
                            else None)
                    flg = flags[band_mask]
                    rej = rejected[band_mask]

                    mask_valid = (~rej) & (flg == 0)
                    mask_flagged = (flg > 0) & (~rej)
                    mask_rejected = rej

                    off = offsets.get(alias, {}).get(
                        str(b) if b is not None else OFFSET_ALL_KEY,
                        offsets.get(alias, {}).get(OFFSET_ALL_KEY, 0.0),
                    )

                    color = self._choose_color_for_band(b, bidx)

                    flagged_color_for_this = self.flagged_color
                    rejected_color_for_this = self.rejected_color

                    if self.parent_gui and hasattr(self.parent_gui, 'plot_settings'):
                        settings = self.parent_gui.plot_settings
                        if settings.flagged_use_filter:
                            flagged_color_for_this = color
                        else:
                            flagged_color_for_this = settings.get_color('flagged', 0)
                        rejected_color_for_this = settings.get_color('rejected', 0)

                    ghost_face = self._ghost_color(color, amount=0.60)
                    ghost_flagged_face = self._ghost_color(flagged_color_for_this, amount=0.60)
                    ghost_rejected_face = self._ghost_color(rejected_color_for_this, amount=0.60)
                    _p = self._cfg.plot if self._cfg is not None else None
                    GHOST_ALPHA = _p.ghost_alpha      if _p else 0.45
                    GHOST_ZORDER = 7
                    GHOST_EW    = _p.ghost_edge_width if _p else 1.2

                    def make_series(m, c, X_override=None, face_color=None, edge_color=None,
                                    edge_width=None, alpha=1.0, zorder_marker=10):
                        if not np.any(m):
                            return None

                        def _plot_sub_series(sub_m, marker_form):
                            if not np.any(sub_m):
                                return None
                            Xplot = X_override if X_override is not None else X
                            Y_plot = Ybase[sub_m] + off
                            mfc = face_color if face_color is not None else c
                            mec = edge_color if edge_color is not None else c
                            mew = edge_width if edge_width is not None else 0.5
                            want_err = (Yerr is not None and
                                        getattr(self.parent_gui, 'errorbar_type', 'none') != 'none')
                            if want_err:
                                container = self.ax.errorbar(
                                    Xplot[sub_m], Y_plot, yerr=Yerr[sub_m], fmt='none', ecolor=mec,
                                    capsize=self.errorbar_capsize, capthick=self.errorbar_capthick,
                                    elinewidth=self.errorbar_linewidth,
                                    zorder=zorder_marker - 5, alpha=alpha,
                                )
                                line, = self.ax.plot(
                                    Xplot[sub_m], Y_plot, linestyle='None', marker=marker_form,
                                    markersize=self.marker_size, color=mfc,
                                    markerfacecolor=mfc, markeredgecolor=mec,
                                    markeredgewidth=mew, zorder=zorder_marker, picker=5, alpha=alpha,
                                )
                                container.marker_line = line
                                return container
                            else:
                                line, = self.ax.plot(
                                    Xplot[sub_m], Y_plot, linestyle='None', marker=marker_form,
                                    markersize=self.marker_size, color=mfc,
                                    markerfacecolor=mfc, markeredgecolor=mec,
                                    markeredgewidth=mew, zorder=zorder_marker, picker=5, alpha=alpha,
                                )
                                return line

                        if use_batch_markers and batch_col is not None and batches_list:
                            res_handles = None
                            for b_val in batches_list:
                                sub_m = m & (batch_col == b_val)
                                if np.any(sub_m):
                                    b_marker = self._choose_marker_for_batch(b_val, batches_list)
                                    h_sub = _plot_sub_series(sub_m, b_marker)
                                    res_handles = _combine(res_handles, h_sub)
                            sub_nan = m & pd.isna(batch_col)
                            if np.any(sub_nan):
                                h_sub = _plot_sub_series(sub_nan, self.marker_style)
                                res_handles = _combine(res_handles, h_sub)
                            return res_handles
                        else:
                            return _plot_sub_series(m, self.marker_style)

                    if time_mode == 'rotation_phase' and phase_max > 1.0:
                        c1 = is_cycle1
                        lc_mask = ~is_cycle1
                        lc_in_range = lc_mask & (X_raw <= phase_max + 1e-9)

                        # LEFT SIDE
                        h_valid = make_series(mask_valid & c1, color)
                        h_flagged = make_series(mask_flagged & c1, flagged_color_for_this,
                                                face_color=color,
                                                edge_color=flagged_color_for_this, edge_width=0.5)
                        h_rejected = (make_series(mask_rejected & c1, rejected_color_for_this,
                                                  face_color=color,
                                                  edge_color=rejected_color_for_this, edge_width=0.5)
                                      if show_rejected else None)

                        h_ghost_valid = make_series(mask_valid & lc_mask, ghost_face,
                                                    face_color=ghost_face, edge_color=color,
                                                    edge_width=GHOST_EW,
                                                    alpha=GHOST_ALPHA, zorder_marker=GHOST_ZORDER)
                        h_ghost_flagged = make_series(mask_flagged & lc_mask, ghost_flagged_face,
                                                      face_color=ghost_flagged_face,
                                                      edge_color=flagged_color_for_this,
                                                      edge_width=GHOST_EW,
                                                      alpha=GHOST_ALPHA, zorder_marker=GHOST_ZORDER)
                        h_ghost_rejected = (make_series(mask_rejected & lc_mask, ghost_rejected_face,
                                                        face_color=ghost_rejected_face,
                                                        edge_color=rejected_color_for_this,
                                                        edge_width=GHOST_EW,
                                                        alpha=GHOST_ALPHA, zorder_marker=GHOST_ZORDER)
                                            if show_rejected else None)

                        # RIGHT SIDE
                        h_valid = _combine(h_valid,
                                           make_series(mask_valid & lc_in_range, color,
                                                       X_override=X_raw))
                        h_flagged = _combine(h_flagged,
                                             make_series(mask_flagged & lc_in_range, flagged_color_for_this,
                                                         X_override=X_raw,
                                                         face_color=color,
                                                         edge_color=flagged_color_for_this, edge_width=0.5))
                        if show_rejected:
                            h_rejected = _combine(h_rejected,
                                                  make_series(mask_rejected & lc_in_range, rejected_color_for_this,
                                                              X_override=X_raw,
                                                              face_color=color,
                                                              edge_color=rejected_color_for_this, edge_width=0.5))

                        # Cycle-1 echo -> ghost
                        Xg = X + 1.0
                        in_ext = (Xg <= phase_max + 1e-9) & c1
                        h_ghost_valid = _combine(h_ghost_valid,
                                                 make_series(mask_valid & in_ext, ghost_face,
                                                             X_override=Xg,
                                                             face_color=ghost_face, edge_color=color,
                                                             edge_width=GHOST_EW,
                                                             alpha=GHOST_ALPHA, zorder_marker=GHOST_ZORDER))
                        h_ghost_flagged = _combine(h_ghost_flagged,
                                                   make_series(mask_flagged & in_ext, ghost_flagged_face,
                                                               X_override=Xg,
                                                               face_color=ghost_flagged_face,
                                                               edge_color=flagged_color_for_this,
                                                               edge_width=GHOST_EW,
                                                               alpha=GHOST_ALPHA, zorder_marker=GHOST_ZORDER))
                        if show_rejected:
                            h_ghost_rejected = _combine(h_ghost_rejected,
                                                        make_series(mask_rejected & in_ext, ghost_rejected_face,
                                                                    X_override=Xg,
                                                                    face_color=ghost_rejected_face,
                                                                    edge_color=rejected_color_for_this,
                                                                    edge_width=GHOST_EW,
                                                                    alpha=GHOST_ALPHA,
                                                                    zorder_marker=GHOST_ZORDER))
                    else:
                        # Normal (non-phase) mode
                        h_valid = make_series(mask_valid, color)
                        h_flagged = make_series(mask_flagged, flagged_color_for_this,
                                                face_color=color,
                                                edge_color=flagged_color_for_this, edge_width=0.5)
                        h_rejected = (make_series(mask_rejected, rejected_color_for_this,
                                                  face_color=color,
                                                  edge_color=rejected_color_for_this, edge_width=0.5)
                                      if show_rejected else None)
                        h_ghost_valid = h_ghost_flagged = h_ghost_rejected = None

                    entry = {
                        'alias': alias, 'band': b, 'band_index': bidx,
                        'X': X, 'Ybase': Ybase, 'Yerr': Yerr,
                        'mask_valid': mask_valid, 'mask_flagged': mask_flagged, 'mask_rejected': mask_rejected,
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
                    if h_ghost_valid is not None:
                        entry['ghost_valid'] = h_ghost_valid
                    if h_ghost_flagged is not None:
                        entry['ghost_flagged'] = h_ghost_flagged
                    if h_ghost_rejected is not None:
                        entry['ghost_rejected'] = h_ghost_rejected

                    self.plotted_handles.setdefault(alias, []).append(entry)

                    legend_label = alias if b is None else f"{alias} ({b})"
                    self.ax.plot([], [], marker=self.marker_style, linestyle='None',
                                 color=color, label=legend_label)
                    plotted_any = True

            if plotted_any:
                try:
                    _legend_fs = self._cfg.plot.legend_fontsize if self._cfg else 'small'
                    if getattr(self.parent_gui, 'show_legend', True):
                        self.ax.legend(title='Lightcurves', fontsize=_legend_fs)
                    else:
                        leg = self.ax.get_legend()
                        if leg is not None:
                            leg.remove()
                except Exception:
                    try:
                        self.ax.legend(title='Lightcurves', fontsize='small')
                    except Exception:
                        pass

            # -- Draw Color Info on Plot --------------------------------------
            if self.parent_gui and getattr(self.parent_gui.show_colors_on_plot_var, 'get', lambda: False)():
                info_lines = []
                vis_aliases = [a for a in lc_dict.keys() if visible.get(a, True)]
                if hasattr(self.parent_gui, 'calculated_colors'):
                    for alias in vis_aliases:
                        colors = self.parent_gui.calculated_colors.get(alias, {})
                        if colors:
                            for name, cdata in colors.items():
                                val, err = cdata[0], cdata[1]
                                info_lines.append(f"  {name}: {val:.3f} +/- {err:.3f}")
                if info_lines:
                    import matplotlib.offsetbox as offsetbox
                    loc_val = "upper left"
                    if self.parent_gui:
                        loc_val = self.parent_gui.color_legend_loc_var.get()
                    _info_fs = self._cfg.plot.label_fontsize if self._cfg else 9
                    text_str = "\n".join(info_lines)
                    at = offsetbox.AnchoredText(text_str, loc=loc_val, frameon=True,
                                                prop=dict(fontsize=_info_fs))
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
                            self.title = format_asteroid_name(lc.df['target'].iloc[0])
                        except Exception:
                            pass

            # Y label based on active mode/corrections (only update the auto-label;
            # user override is respected in _set_labels via self._user_ylabel)
            if mode == 'target':
                use_rm = False
                if self.parent_gui is not None:
                    try:
                        use_rm = bool(self.parent_gui.use_reduced_mag_var.get())
                    except Exception:
                        pass
                self.ylabel = "Reduced Magnitude H(a)" if use_rm else "Calibrated Magnitude"
            elif mode == 'instrumental':
                self.ylabel = "Instrumental Magnitude"
            elif mode == 'relative':
                self.ylabel = "Relative Magnitude"
            elif mode == 'control':
                self.ylabel = "Calibrated Magnitude"
            else:
                self.ylabel = "Magnitude"

            # Lighttime note in xlabel
            use_lt = False
            if self.parent_gui is not None:
                try:
                    use_lt = bool(self.parent_gui.use_lighttime_var.get())
                except Exception:
                    pass
            if use_lt and "(LT-corr)" not in self.xlabel:
                self.xlabel = self.xlabel + " (LT-corr)"

            # -- Draw batch indicators in edit mode ---------------------------
            if getattr(self.parent_gui, 'edit_mode', False):
                for alias, lc in lc_dict.items():
                    if not visible.get(alias, True) or lc.df is None or not lc.has_batches():
                        continue
                    batches = lc.get_batches()
                    if len(batches) <= 1:
                        continue
                    cmap = plt.get_cmap('tab10')
                    for b_i, b_val in enumerate(batches):
                        b_idx = lc.get_batch_indices(b_val)
                        if not b_idx:
                            continue
                        jd_b = lc.arr['jd'][b_idx]
                        if jd_b.size == 0 or np.all(np.isnan(jd_b)):
                            continue
                        X_b = self._compute_time_x(jd_b, time_mode, period)
                        if time_mode == 'rotation_phase':
                            X_b = np.mod(X_b, 1.0)
                        valid_xb = X_b[np.isfinite(X_b)]
                        if valid_xb.size == 0:
                            continue
                        xmin_b, xmax_b = float(np.min(valid_xb)), float(np.max(valid_xb))
                        b_color = cmap(b_i % 10)
                        self.ax.axvspan(xmin_b, xmax_b, color=b_color, alpha=0.08, zorder=0)

            # rotation_phase x limits and wrap boundary
            if time_mode == 'rotation_phase':
                try:
                    _pm = DEFAULT_PHASE_MAX
                    if self.parent_gui is not None:
                        try:
                            _pm = float(self.parent_gui.phase_max_var.get())
                            if not (1.0 < _pm <= 3.0):
                                _pm = DEFAULT_PHASE_MAX
                        except Exception:
                            pass
                    self.ax.set_xlim(0.0, _pm)
                    if _pm > 1.0:
                        _p = self._cfg.plot if self._cfg is not None else None
                        self.ax.axvline(
                            x=1.0,
                            color=_p.phase_wrap_color      if _p else 'gray',
                            linestyle='--',
                            linewidth=_p.phase_wrap_linewidth if _p else 0.8,
                            alpha=_p.phase_wrap_alpha        if _p else 0.6,
                            zorder=1,
                        )
                except Exception:
                    pass

            self._set_labels(offsets)
            self._update_selection_artists(lc_dict, offsets, mode, time_mode)
            self._register_blit_artists()
            _safe(self.fig.tight_layout)

            if hasattr(self, "blit") and self.blit is not None:
                self.blit._bg = None
                self.canvas.draw()
                self.blit.quick_redraw()

            self._layout_signature = new_sig
            return

        # -- fast path (update existing artists only) -------------------------
        debug_print("Only fast update needed")
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
                    if hasattr(handle, 'lines'):
                        err_on = bool(visible and (getattr(self.parent_gui, 'errorbar_type', 'none') != 'none'))
                        if len(handle.lines) > 0 and handle.lines[0] is not None:
                            _vis(handle.lines[0], visible)
                        if len(handle.lines) > 1:
                            _vis(handle.lines[1], err_on)
                        if len(handle.lines) > 2:
                            _vis(handle.lines[2], err_on)
                        if hasattr(handle, 'marker_line'):
                            _vis(handle.marker_line, visible)
                        return
                    if isinstance(handle, (list, tuple)):
                        for o in handle:
                            _vis(o, visible)
                        return
                    _vis(handle, visible)

                def set_data(h, mask, color_override=None, face_color=None,
                             edge_color=None, edge_width=None):
                    if h is None:
                        return
                    Xm = X[mask]
                    Ym = Y[mask]

                    def _apply_marker_colors(artist):
                        if artist is None:
                            return
                        try:
                            if face_color is not None:
                                artist.set_markerfacecolor(face_color)
                                artist.set_color(face_color)
                            elif color_override is not None:
                                artist.set_color(color_override)
                            if edge_color is not None:
                                artist.set_markeredgecolor(edge_color)
                            if edge_width is not None:
                                artist.set_markeredgewidth(edge_width)
                        except Exception:
                            pass

                    if hasattr(h, 'lines'):
                        main = h.lines[0] if h.lines else None
                        if main is not None:
                            main.set_data(Xm, Ym)
                            main.set_marker(self.marker_style)
                            main.set_markersize(self.marker_size)
                            main.set_linestyle('None')
                            _apply_marker_colors(main)

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
                            bar_cols = h.lines[2] if len(h.lines) > 2 else ()
                            for lc_col in bar_cols:
                                lc_col.set_segments(segs)
                                lc_col.set_visible(True)
                            caps = h.lines[1] if len(h.lines) > 1 else ()
                            for cap in caps:
                                cap.set_visible(True)
                        else:
                            bar_cols = h.lines[2] if len(h.lines) > 2 else ()
                            for lc_col in bar_cols:
                                lc_col.set_visible(False)
                            caps = h.lines[1] if len(h.lines) > 1 else ()
                            for cap in caps:
                                cap.set_visible(False)

                        if hasattr(h, 'marker_line') and h.marker_line is not None:
                            ml = h.marker_line
                            ml.set_data(Xm, Ym)
                            ml.set_marker(self.marker_style)
                            ml.set_markersize(self.marker_size)
                            ml.set_linestyle('None')
                            _apply_marker_colors(ml)
                        return

                    if isinstance(h, (list, tuple)):
                        for hh in h:
                            set_data(hh, mask)
                        return

                    h.set_data(Xm, Ym)
                    h.set_marker(self.marker_style)
                    h.set_markersize(self.marker_size)
                    h.set_linestyle('None')
                    _apply_marker_colors(h)

                mask_v = e['mask_valid']
                mask_f = e['mask_flagged']
                mask_r = e['mask_rejected']

                h = e.get('valid')
                if h is not None:
                    set_data(h, mask_v, e.get('color'))

                h = e.get('flagged')
                if h is not None:
                    set_data(h, mask_f, face_color=e.get('color'),
                             edge_color=e.get('flagged_color'), edge_width=2.0)

                h = e.get('rejected')
                if h is not None:
                    set_data(h, mask_r, face_color=e.get('color'),
                             edge_color=e.get('rejected_color'), edge_width=2.0)

                want_r = vis_alias and show_rejected and np.any(e['mask_rejected'])
                want_v = vis_alias and np.any(e['mask_valid'])
                want_f = vis_alias and np.any(e['mask_flagged'])

                for key, want in (('valid', want_v), ('flagged', want_f), ('rejected', want_r),
                                   ('ghost_valid', want_v), ('ghost_flagged', want_f),
                                   ('ghost_rejected', want_r)):
                    h = e.get(key)
                    if h is not None:
                        _set_visible(h, want)

        # -- selection marker update ------------------------------------------
        self._update_selection_artists(lc_dict, offsets, mode, time_mode)
        self._set_labels(offsets)
        self.blit.quick_redraw()

    def _update_selection_artists(self, lc_dict, offsets, mode, time_mode):
        """Update positions and visibility of single selection and multi-selection markers."""
        sel = getattr(self, 'selection', {'alias': None, 'index': None})
        sel_alias = sel.get('alias')
        sel_idx = sel.get('index')

        # 1. Multi-selection points
        sel_indices = set()
        if self.parent_gui is not None and hasattr(self.parent_gui, 'selected_point_indices'):
            sel_indices = set(self.parent_gui.selected_point_indices)

        if sel_alias in lc_dict and sel_indices:
            lc = lc_dict[sel_alias]
            df = lc.df
            arr = lc.arr
            if df is not None and len(df) > 0:
                valid_sel_idx = [i for i in sorted(sel_indices) if 0 <= i < len(df)]
                if valid_sel_idx:
                    jd_all = arr.get('jd')
                    period = None
                    if time_mode == 'rotation_phase' and self.parent_gui is not None:
                        try:
                            period = float(self.parent_gui.rotation_period_var.get())
                        except Exception:
                            period = None
                    X_all = self._compute_time_x(jd_all, time_mode, period)
                    if time_mode == 'rotation_phase':
                        X_all = np.mod(X_all, 1.0)

                    if mode == 'target' and 'mag' in df.columns:
                        Y_all = arr.get('mag')
                    elif mode == 'instrumental' and 'inst_mag' in df.columns:
                        Y_all = arr.get('inst_mag')
                    elif mode == 'relative' and 'rel_mag' in df.columns:
                        Y_all = arr.get('rel_mag')
                    else:
                        Y_all = arr.get('mag_control', pd.Series(np.zeros(len(df))))

                    bands = arr.get('band')
                    multi_x = []
                    multi_y = []
                    for idx_val in valid_sel_idx:
                        b_val = str(bands[idx_val]) if bands is not None and idx_val < len(bands) else None
                        off = offsets.get(sel_alias, {}).get(
                            b_val, offsets.get(sel_alias, {}).get(OFFSET_ALL_KEY, 0.0))
                        multi_x.append(X_all[idx_val])
                        multi_y.append(Y_all[idx_val] + off)
                    self.multi_sel_artist.set_data(multi_x, multi_y)
                    self.multi_sel_artist.set_visible(True)
                else:
                    self.multi_sel_artist.set_visible(False)
            else:
                self.multi_sel_artist.set_visible(False)
        else:
            self.multi_sel_artist.set_visible(False)

        # 2. Single focused selection point and text label
        if sel_alias in lc_dict and sel_idx is not None:
            lc = lc_dict[sel_alias]
            df = lc.df
            arr = lc.arr
            if df is not None and 0 <= sel_idx < len(df):
                jd_all = arr.get('jd')
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
                if time_mode == 'rotation_phase':
                    X_all = np.mod(X_all, 1.0)
                sx = X_all[sel_idx]
                if mode == 'target' and 'mag' in df.columns:
                    sy = float(df['mag'].iloc[sel_idx])
                elif mode == 'instrumental' and 'inst_mag' in df.columns:
                    sy = float(df['inst_mag'].iloc[sel_idx])
                elif mode == 'relative' and 'rel_mag' in df.columns:
                    sy = float(df['rel_mag'].iloc[sel_idx])
                else:
                    sy = float(df.get('mag_control', pd.Series(np.zeros(len(df))))[sel_idx])
                off = offsets.get(sel_alias, {}).get(
                    band_val, offsets.get(sel_alias, {}).get(OFFSET_ALL_KEY, 0.0))
                self.sel_artist.set_data([sx], [sy + off])
                self.sel_artist.set_visible(True)
                txt = ''
                try:
                    if 'batch_number' in df.columns and not pd.isna(df['batch_number'].iloc[sel_idx]):
                        bnum = df['batch_number'].iloc[sel_idx]
                        txt += f"[B:{bnum}] "
                    elif 'batch' in df.columns and not pd.isna(df['batch'].iloc[sel_idx]):
                        bnum = df['batch'].iloc[sel_idx]
                        txt += f"[B:{bnum}] "
                    if 'filename' in df.columns:
                        txt += str(df['filename'].iloc[sel_idx])
                except Exception:
                    txt = ''
                if getattr(self, 'sel_text', None) is not None and txt:
                    try:
                        lo, hi = self.ax.get_ylim()
                        yspan = abs(hi - lo) if hi is not None and lo is not None else 0.0
                        y_offset = 0.05 * (yspan if yspan > 0 else 1.0)
                        xlim = self.ax.get_xlim()
                        x_range = xlim[1] - xlim[0]
                        text_width_estimate = 0.25 * x_range
                        if sx + text_width_estimate > xlim[1]:
                            self.sel_text.set_position((sx, sy + off + y_offset))
                            self.sel_text.set_horizontalalignment('right')
                        else:
                            self.sel_text.set_position((sx, sy + off + y_offset))
                            self.sel_text.set_horizontalalignment('left')
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

    #  -  -  -  point picking  -  -  - 

    def select_points_in_region(
        self,
        xmin: float,
        xmax: float,
        ymin: float,
        ymax: float,
        lc_dict: Dict[str, LightCurveData],
        offsets: Dict[str, Dict[str, float]],
        mode: str,
        time_mode: str,
        selected_bands: Set[str],
        show_rejected: bool,
    ) -> List[Tuple[str, int, Optional[str]]]:
        """Find all data points whose (X, Y) fall within [xmin, xmax] and [ymin, ymax].

        Returns a list of (alias, dataframe_index, band) tuples.
        """
        x_lo, x_hi = min(xmin, xmax), max(xmin, xmax)
        y_lo, y_hi = min(ymin, ymax), max(ymin, ymax)

        matched = []
        for alias, lc in lc_dict.items():
            if lc.df is None or lc.df.empty:
                continue
            df = lc.df
            arr = lc.arr
            jd = arr.get('jd')
            if jd is None or not jd.size:
                continue
            bands = arr.get('band')
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
            if time_mode == 'rotation_phase':
                X_all = np.mod(X_all, 1.0)

            if mode == 'target':
                Y_all = arr.get('mag')
            elif mode == 'instrumental':
                Y_all = arr.get('inst_mag')
            elif mode == 'relative':
                Y_all = arr.get('rel_mag')
            else:
                Y_all = arr.get('mag_control')
            if Y_all is None:
                continue

            off_all = np.zeros(len(df))
            if bands is not None:
                for i in np.nonzero(vis_mask)[0]:
                    bval = str(bands[i]) if bands is not None else None
                    off_all[i] = offsets.get(alias, {}).get(
                        bval, offsets.get(alias, {}).get(OFFSET_ALL_KEY, 0.0))
            else:
                off_all[:] = offsets.get(alias, {}).get(OFFSET_ALL_KEY, 0.0)

            Y_plotted = Y_all + off_all

            in_box = vis_mask & (X_all >= x_lo) & (X_all <= x_hi) & (Y_plotted >= y_lo) & (Y_plotted <= y_hi)
            matched_indices = np.nonzero(in_box)[0]
            for idx in matched_indices:
                b_val = str(bands[idx]) if bands is not None else None
                matched.append((alias, int(idx), b_val))

        return matched 

    def select_nearest_point(
        self,
        xdata: float,
        ydata: float,
        lc_dict: Dict[str, LightCurveData],
        offsets: Dict[str, Dict[str, float]],
        mode: str,
        time_mode: str,
        selected_bands: Set[str],
        show_rejected: bool,
    ) -> Tuple[Optional[str], Optional[int], Optional[str]]:
        """Find the data point closest (in display pixels) to a mouse click.

        Returns ``(alias, dataframe_index, band)`` for the nearest visible
        point within 15 px, or ``(None, None, None)`` when none is close enough.
        """
        pts = []
        trans = self.ax.transData
        for alias, lc in lc_dict.items():
            if lc.df is None or lc.df.empty:
                continue
            df = lc.df
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
            if time_mode == 'rotation_phase':
                X_all = np.mod(X_all, 1.0)

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
                for i in np.nonzero(vis_mask)[0]:
                    bval = str(bands[i]) if bands is not None else None
                    off_all[i] = offsets.get(alias, {}).get(
                        bval, offsets.get(alias, {}).get(OFFSET_ALL_KEY, 0.0))
            else:
                off_all[:] = offsets.get(alias, {}).get(OFFSET_ALL_KEY, 0.0)

            Xv = X_all[vis_mask]
            Yv = (Y_all + off_all)[vis_mask]
            bands_v = (bands[vis_mask] if bands is not None else np.array([None] * len(Xv)))

            xy = trans.transform(np.column_stack((Xv, Yv)))
            vis_idx = np.nonzero(vis_mask)[0]
            for k, (xd, yd) in enumerate(xy):
                pts.append((xd, yd, alias, int(vis_idx[k]),
                            str(bands_v[k]) if bands is not None else None))

        if not pts:
            return None, None, None
        click_disp = self.ax.transData.transform((xdata, ydata))
        d2 = [(xd - click_disp[0]) ** 2 + (yd - click_disp[1]) ** 2 for xd, yd, *_ in pts]
        i = int(np.argmin(d2))
        if d2[i] < (15.0 ** 2):
            _, _, alias, idx, b = pts[i]
            return alias, idx, b
        return None, None, None


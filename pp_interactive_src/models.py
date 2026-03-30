# -*- coding: utf-8 -*-
"""
Data-layer classes for pp_interactive_src: FitsContext and LightCurveData.
"""
from __future__ import annotations

import os
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from astropy.io import fits

import toolbox

from .constants import OFFSET_ALL_KEY


class FitsContext:
    """FITS-header context for a single observation night.

    Parses an on-disk FITS file once and caches the subset of header values
    needed by the rest of the application: target name and MPC observatory
    code.  Can also be constructed without a file (via ``__new__``) when the
    user provides the values manually.
    """

    def __init__(self, filepath: str):
        self.filepath = filepath
        self.header = None
        self.target_object: Optional[str] = None
        self.observatory_code: Optional[str] = None
        self.obsparam: Optional[dict] = None
        self._load()

    def _load(self):
        """Read FITS header and extract target / observatory information."""
        try:
            self.header = toolbox.get_fits_header(self.filepath)
            self.obsparam = toolbox.get_obsparam(self.header)
            self.target_object = self.header.get('OBJECT', 'Unknown')
            if self.obsparam and 'mpc_code' in self.obsparam:
                self.observatory_code = str(self.obsparam['mpc_code'])
            else:
                self.observatory_code = str(
                    self.header.get('MPCCODE', self.header.get('OBSERVAT', '500'))
                )
        except Exception as e:
            print(f"Error loading FITS header from {self.filepath}: {e}")
            self.header = fits.Header()


class LightCurveData:
    """Container for one CSV light-curve plus cached numpy arrays for fast plotting.

    Attributes:
        filename:  Path to the source CSV file.
        df:        Full DataFrame (the authoritative data store).
        arr:       Dict of pre-cast numpy arrays mirroring key DataFrame columns;
                   kept in sync after every mutation so the plot layer never has
                   to re-allocate on each redraw.
    """

    def __init__(self, filepath: Optional[str] = None):
        self.filename: Optional[str] = None
        self.df: Optional[pd.DataFrame] = None
        self.arr: Dict[str, np.ndarray] = {}
        if filepath:
            self.load(filepath)

    # ——— I/O ———
    def load(self, filepath: str) -> None:
        """Load a CSV file, ensure required columns exist, and prime the array cache."""
        self.filename = filepath
        self.df = pd.read_csv(filepath)
        self._ensure_columns()
        self._fill_arrays_cache()

    def save(self) -> None:
        """Persist the current DataFrame back to the original CSV path."""
        if self.df is not None and self.filename:
            self.df.to_csv(self.filename, index=False)

    # ——— internals ———
    def _ensure_columns(self):
        """Add mandatory columns with sensible defaults if they are absent."""
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
        # JPL-derived optional columns
        self.arr['reduced_mag'] = df.get('reduced_mag', pd.Series(dtype=float)).to_numpy(dtype=float, copy=False)
        self.arr['corrected_jd'] = df.get('corrected_jd', pd.Series(dtype=float)).to_numpy(dtype=float, copy=False)

    def toggle_rejection(self, index: int) -> None:
        """Toggle the rejection flag for a given row index."""
        if self.df is None:
            return
        self.df.loc[index, 'rejected'] = not bool(self.df.loc[index, 'rejected'])
        self.arr['rejected'] = self.df['rejected'].to_numpy(dtype=bool, copy=False)

    def get_bands(self) -> List[str]:
        """Return the unique photometric filter bands, in the canonical order."""
        if self.df is None or 'band' not in self.df.columns:
            return []
        vals = self.df['band'].dropna().astype(str).str.strip().unique().tolist()
        desired = ['U', 'B', 'V', 'R', 'I', 'g', 'r', 'i', 'z']
        return sorted(vals, key=lambda b: desired.index(b) if b in desired else 999)


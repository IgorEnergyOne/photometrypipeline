#!/usr/bin/env python3

"""
Interactive cone-search viewer PAN-STARRS DR1 catalog
"""

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astropy.table import Table, Column
from astropy.coordinates import Angle
from scipy import spatial
import re

MAX_SOURCES = 50_000

DISPLAY_COLS = [
    ("ID GAIA",        lambda r: int(r.get("id_gaia", 0))),
    ("ID_PS1",         lambda r: int(r.get("id_panstarrs_PS1", 0))),
    ("RA [deg]",       lambda r: f"{r['ra_deg']:.6f}"),
    ("DEC [deg]",      lambda r: f"{r['dec_deg']:.6f}"),
    ("offset [arcsec]",lambda r: f"{r['dist']:.3f}"),
    ("B [mag]",        lambda r: f"{r['Bmag']:.3f}"),
    ("V [mag]",        lambda r: f"{r['Vmag']:.3f}"),
    ("R [mag]",        lambda r: f"{r['Rmag']:.3f}"),
    ("I [mag]",        lambda r: f"{r['Imag']:.3f}"),
    ("B-V",            lambda r: f"{r['B-V']:.3f}"),
    ("V-R",            lambda r: f"{r['V-R']:.3f}"),
    ("R-I",            lambda r: f"{r['R-I']:.3f}"),
    ("B(GAIA)-B(PS1)", lambda r: f"{r['B_GAIA-B_PS1']:.3f}"),
    ("V(GAIA)-V(PS1)", lambda r: f"{r['V_GAIA-V_PS1']:.3f}"),
    ("R(GAIA)-R(PS1)", lambda r: f"{r['R_GAIA-R_PS1']:.3f}"),
    ("I(GAIA)-I(PS1)", lambda r: f"{r['I_GAIA-I_PS1']:.3f}"),
]

# Map label → lambda for fast lookup
_GETVAL = dict(DISPLAY_COLS)


def _as_plain(values):
    """Return a plain ndarray, replacing masked entries with NaN."""
    if hasattr(values, "filled"):      # MaskedArray or MaskedColumn
        return values.filled(np.nan)
    return values

_DMS_RE = re.compile(
    r"""^\s*           # optional leading spaces
        (?P<sgn>[+-])? # optional sign
        (?P<d>\d+)(?:[d°\s:]|\s+|$)     # degrees
        (?:(?P<m>\d+)(?:[m'\s:]|\s+|$))?# minutes
        (?:(?P<s>\d+(?:\.\d*)?)
            (?:[s"]|\s*)                # seconds
        )?                              # seconds optional
        \s*$""",
    re.VERBOSE)

def _parse_angle(text: str, type='ra') -> float:
    """
    Accepts any of:
      • decimal degrees   →  187.12345
      • DMS               →  187d07m24.4s   or  187 07 24.4
      • HMS (RA)          →   12h27m33.1s   or   12 27 33.1
    Returns the value in **degrees**.
    """
    text = text.strip()
    # plain decimal?  (fast path)
    try:
        return float(text)
    except ValueError:
        pass

    # give astropy a shot – it understands most astronomical notations
    try:
        return Angle(text).degree
    except Exception:
        # try “12 27 33.1” → “12d27m33.1s”   or “…h…”
        parts = text.replace(':', ' ').split()
        if len(parts) == 3 and all(p.replace('.', '', 1).isdigit()
                                   or (p[0] in '+-' and p[1:].replace('.', '', 1).isdigit())
                                   for p in parts):
            # assume dms unless explicitly > 24 → treat as degrees
            first = float(parts[0])
            if type == 'ra':
                if abs(first) <= 24:                       # could be RA in hours
                    guess = f"{parts[0]}h{parts[1]}m{parts[2]}s"
                else:                                      # definitely degrees
                    guess = f"{parts[0]}d{parts[1]}m{parts[2]}s"
                try:
                    return Angle(guess).degree
                except Exception:
                    pass
            elif type == 'dec':
                guess = f"{parts[0]}d{parts[1]}m{parts[2]}s"
                try:
                    return Angle(guess).degree
                except Exception:
                    pass
        raise ValueError(f"Bad angle format: “{text}”")


def _parse_center(text: str):
    """
    Split a free-form “RA Dec” string into two parts and convert each to degrees.
    Works for:
      • decimal  → 187.12 +2.34
      • DMS      → 187 07 12.0  +02 20 24
      • HMS/DMS  → 12 27 33.1  -03:22:10.1
      • comma-separated, etc.
    """
    tokens = text.replace(',', ' ').split()
    if len(tokens) < 2:
        raise ValueError("Enter both RA and Dec (separated by space or comma).")

    # try every possible split: tokens[:i] = RA , tokens[i:] = Dec
    for i in range(1, len(tokens)):
        ra_str  = ' '.join(tokens[:i])
        dec_str = ' '.join(tokens[i:])
        try:
            ra  = _parse_angle(ra_str, type='ra')
            dec = _parse_angle(dec_str, type='dec')
            return ra, dec
        except ValueError:
            continue

    raise ValueError(f"Cannot interpret RA/Dec in “{text}”.")


def _is_number(txt: str):
    try:
        float(txt)
        return True
    except ValueError:
        return False


def sources_with(data, condition):
    """
    reject sources based on condition
    input: condition
    return: number of sources left
    """

    data_select = data[condition]
    return data_select

# ------------------------------------------------------------
def fetch_catalog_ps1(ra_deg: float, dec_deg: float, radius: float,
                  mag_max: 21, min_mag: 0, band: 'Rmag') -> pd.DataFrame:
    """Generic Vizier cone search returning a pandas.DataFrame."""
    field  = SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg),
                       frame='icrs')
    radius_deg = radius / 3600.0

    vquery = Vizier(columns=['objID', 'RAJ2000', 'DEJ2000',
                             'e_RAJ2000', 'e_DEJ2000',
                             'gmag', 'e_gmag',
                             'rmag', 'e_rmag',
                             'imag', 'e_imag',
                             'zmag', 'e_zmag',
                             'ymag', 'e_ymag'],
                    column_filters={"rmag":
                                        (">={:f} AND <={:f}".format(-5, 30))
                                    },
                    row_limit=MAX_SOURCES,
                    timeout=300)
    data = vquery.query_region(field,
                               radius=radius_deg * u.deg,
                               catalog="II/349/ps1",
                               cache=False)[0]

    # calculate distance from query point (in arcseconds)
    data['dist'] = np.sqrt((data['RAJ2000'] - ra_deg)**2
                           + (data['DEJ2000'] - dec_deg)**2) * 3600
    # sort by distance fro the query poin
    data.sort('dist')

    # rename column names using PP conventions
    data.rename_column('objID', 'id_panstarrs')
    data.rename_column('RAJ2000', 'ra_deg')
    data.rename_column('DEJ2000', 'dec_deg')
    data.rename_column('e_RAJ2000', 'e_ra_deg')
    data['e_ra_deg'].convert_unit_to(u.deg)
    data.rename_column('e_DEJ2000', 'e_dec_deg')
    data['e_dec_deg'].convert_unit_to(u.deg)
    data.rename_column('gmag', 'gp1mag')
    data.rename_column('e_gmag', 'e_gp1mag')
    data.rename_column('rmag', 'rp1mag')
    data.rename_column('e_rmag', 'e_rp1mag')
    data.rename_column('imag', 'ip1mag')
    data.rename_column('e_imag', 'e_ip1mag')
    data.rename_column('zmag', 'zp1mag')
    data.rename_column('e_zmag', 'e_zp1mag')
    data.rename_column('ymag', 'yp1mag')
    data.rename_column('e_ymag', 'e_yp1mag')
    data['mag'] = data['rp1mag']  # use rmag for astrometry

    g = data['gp1mag'].data
    e_g = data['e_gp1mag'].data
    r = data['rp1mag'].data
    e_r = data['e_rp1mag'].data
    i = data['ip1mag'].data
    e_i = data['e_ip1mag'].data
    z = data['zp1mag'].data
    e_z = data['e_zp1mag'].data

    # transform to Johnson Cousins BVRI
    B = (g + 0.212 + 0.556 * (g - r) + 0.034 * (g - r) ** 2)
    Berr = np.sqrt(e_g ** 2 + 0.032 ** 2)
    V = (g + 0.005 - 0.536 * (g - r) + 0.011 * (g - r) ** 2)
    Verr = np.sqrt(e_g ** 2 + 0.012 ** 2)
    R = (r - 0.137 - 0.108 * (g - r) - 0.029 * (g - r) ** 2)
    Rerr = np.sqrt(e_r ** 2 + 0.015 ** 2)
    I = (i - 0.366 - 0.136 * (g - r) - 0.018 * (g - r) ** 2)
    Ierr = np.sqrt(e_i ** 2 + 0.017 ** 2)

    data.add_column(Column(data=B, name='Bmag', unit=u.mag))
    data.add_column(Column(data=Berr, name='e_Bmag',
                           unit=u.mag))
    data.add_column(Column(data=V, name='Vmag', unit=u.mag))
    data.add_column(Column(data=Verr, name='e_Vmag',
                           unit=u.mag))
    data.add_column(Column(data=R, name='Rmag', unit=u.mag))
    data.add_column(Column(data=Rerr, name='e_Rmag',
                           unit=u.mag))
    data.add_column(Column(data=I, name='Imag', unit=u.mag))
    data.add_column(Column(data=Ierr, name='e_Imag',
                           unit=u.mag))


    # transform to SDSS to select solar-like stars
    g_sdss = (g + 0.013 + 0.145 * (g - r) + 0.019 * (g - r) ** 2)
    gerr_sdss = np.sqrt(e_g ** 2 + 0.008 ** 2)
    r_sdss = (r - 0.001 + 0.004 * (g - r) + 0.007 * (g - r) ** 2)
    rerr_sdss = np.sqrt(e_r ** 2 + 0.004 ** 2)
    i_sdss = (i - 0.005 + 0.011 * (g - r) + 0.010 * (g - r) ** 2)
    ierr_sdss = np.sqrt(e_i ** 2 + 0.004 ** 2)
    z_sdss = (z + 0.013 - 0.039 * (g - r) - 0.012 * (g - r) ** 2)
    zerr_sdss = np.sqrt(e_z ** 2 + 0.01 ** 2)

    data.add_column(Column(data=_as_plain(g_sdss), name='_gmag',
                           unit=u.mag))
    data.add_column(Column(data=_as_plain(gerr_sdss), name='_e_gmag',
                           unit=u.mag))
    data.add_column(Column(data=_as_plain(r_sdss), name='_rmag',
                           unit=u.mag))
    data.add_column(Column(data=_as_plain(rerr_sdss), name='_e_rmag',
                           unit=u.mag))
    data.add_column(Column(data=_as_plain(i_sdss), name='_imag',
                           unit=u.mag))
    data.add_column(Column(data=_as_plain(ierr_sdss), name='_e_imag',
                           unit=u.mag))
    data.add_column(Column(data=_as_plain(z_sdss), name='_zmag',
                           unit=u.mag))
    data.add_column(Column(data=_as_plain(zerr_sdss), name='_e_zmag',
                           unit=u.mag))

    data = data.to_pandas()
    # filter catalogue to get stars of magnitude in given range
    data = data[data[f'{band}'].between(min_mag, mag_max)]
    return data

def fetch_catalog_gaia(ra_deg: float, dec_deg: float, radius: float,
                  mag_max: 21, min_mag: 0, band: 'Rmag') -> pd.DataFrame:
    """Generic Vizier cone search returning a pandas.DataFrame."""
    field  = SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg),
                       frame='icrs')
    radius_deg = radius / 3600.0

    vquery = Vizier(columns=['Source', 'RA_ICRS', 'DE_ICRS',
                             'e_RA_ICRS', 'e_DE_ICRS', 'pmRA',
                             'pmDE', 'Epoch',
                             'Gmag', 'e_Gmag',
                             'BPmag', 'e_BPmag',
                             'RPmag', 'eRPmag', 'VarFlag'],
                    column_filters={"phot_g_mean_mag":
                                        (">={:f} AND <={:f}".format(-5, 30)),
                                    "VarFlag": "!=VARIABLE"},
                    row_limit=MAX_SOURCES,
                    timeout=300)
    data = vquery.query_region(field,
                               radius=radius_deg * u.deg,
                               catalog="I/355/gaiadr3",
                               cache=False)[0]

    # calculate distance from query point (in arcseconds)
    data['dist'] = np.sqrt((data['RA_ICRS'] - ra_deg)**2
                           + (data['DE_ICRS'] - dec_deg)**2) * 3600
    # sort by distance for the query point
    data.sort('dist')

    # rename column names using PP conventions
    data.rename_column('Source', 'id_gaia')
    data.rename_column('RA_ICRS', 'ra_deg')
    data.rename_column('DE_ICRS', 'dec_deg')
    # add angles in hms and dms
    data.add_column(Column(data=Angle(data["ra_deg"]).to_string(
        unit=u.hourangle, sep=':', pad=True, precision=2),
                           name="RA_hms"))
    data.add_column(Column(data=Angle(data["dec_deg"]).to_string(
        unit=u.deg, sep=':', pad=True, precision=2, alwayssign=True),
                           name="Dec_dms"))
    g = data['Gmag']
    e_g = data['e_Gmag']
    bp = data['BPmag']
    rp = data['RPmag']

    # transform to Johnson Cousins BVRI
    data['Bmag'] = g - (0.01448 - 0.6874 * (bp - rp) - 0.3604 * (bp - rp) ** 2
                        + 0.06718 * (bp - rp) ** 3 - 0.006061 * (bp - rp) ** 4)
    data['e_Bmag'] = np.sqrt(e_g ** 2 + 0.0633 ** 2)
    data['Vmag'] = g - (-0.02704 + 0.01424 * (bp - rp) - 0.2156 * (bp - rp) ** 2 + 0.01426 * (bp - rp) ** 3)
    data['e_Vmag'] = np.sqrt(e_g ** 2 + 0.03017 ** 2)
    data['Rmag'] = g - (-0.02275 + 0.3961 * (bp - rp) - 0.1243 * (bp - rp) ** 2
                        - 0.01396 * (bp - rp) ** 3 + 0.003775 * (bp - rp) ** 4)
    data['e_Rmag'] = np.sqrt(e_g ** 2 + 0.03167 ** 2)
    data['Imag'] = g - (0.01753 + 0.76 * (bp - rp) - 0.0991 * (bp - rp) ** 2)
    data['e_Imag'] = np.sqrt(e_g ** 2 + 0.03765 ** 2)
    data['B-V'] = data['Bmag'] - data['Vmag']
    data['V-R'] = data['Vmag'] - data['Rmag']
    data['R-I'] = data['Rmag'] - data['Imag']

    # transform to SDSS to select solar-like stars
    g_sdss = g - (0.2199 - 0.6365 * (bp - rp)
                  - 0.1548 * (bp - rp) ** 2 + 0.0064 * (bp - rp) ** 3)
    e_g_sdss = np.sqrt(e_g ** 2 + 0.0745 ** 2)
    r_sdss = g - (-0.09837 + 0.08592 * (bp - rp) + 0.1907 * (bp - rp) ** 2
                  - 0.1701 * (bp - rp) ** 3 + 0.02263 * (bp - rp) ** 4)
    e_r_sdss = np.sqrt(e_g ** 2 + 0.03776 ** 2)
    i_sdss = g - (-0.293 + 0.6404 * (bp - rp) - 0.09609 * (bp - rp) ** 2 - 0.002104 * (bp - rp) ** 3)
    e_i_sdss = np.sqrt(e_g ** 2 + 0.04092 ** 2)
    z_sdss = g - (-0.4619 + 0.8992 * (bp - rp) - 0.08271 * (bp - rp) ** 2 + 0.005029 * (bp - rp) ** 3)
    e_z_sdss = np.sqrt(e_g ** 2 + 0.041161 ** 2)

    data.add_column(Column(data=_as_plain(g_sdss), name='_gmag', unit=u.mag))
    data.add_column(Column(data=_as_plain(e_g_sdss), name='_e_gmag',
                           unit=u.mag))
    data.add_column(Column(data=_as_plain(r_sdss), name='_rmag',
                           unit=u.mag))
    data.add_column(Column(data=_as_plain(e_r_sdss), name='_e_rmag',
                           unit=u.mag))
    data.add_column(Column(data=_as_plain(i_sdss), name='_imag',
                           unit=u.mag))
    data.add_column(Column(data=_as_plain(e_i_sdss), name='_e_imag',
                           unit=u.mag))
    data.add_column(Column(data=_as_plain(z_sdss), name='_zmag',
                           unit=u.mag))
    data.add_column(Column(data=_as_plain(e_z_sdss), name='_e_zmag',
                           unit=u.mag))
    data = data.to_pandas()
    # filter catalogue to get stars of magnitude in given range
    data = data[data[f'{band}'].between(min_mag, mag_max)]
    return data

def solar_mask(df: pd.DataFrame, solar_coef: float = 0.2) -> pd.DataFrame:
    """select solar-like stars based on their colors in sdss and likeness coefficient"""
    sol_gr, sol_ri = 0.44, 0.11
    mask = ((df['_gmag'] - df['_rmag']).between(sol_gr - solar_coef, sol_gr + solar_coef) &
            (df['_rmag'] - df['_imag']).between(sol_ri - solar_coef, sol_ri + solar_coef))
    return df[mask]


def skycoord_match(gaia_cat: pd.DataFrame, ps1_cat: pd.DataFrame) -> pd.DataFrame:
    gaia_coord = SkyCoord(gaia_cat['ra_deg'].values * u.deg, gaia_cat['dec_deg'].values * u.deg)
    ps1_coord = SkyCoord(ps1_cat['ra_deg'].values * u.deg, ps1_cat['dec_deg'].values * u.deg)

    idx, sep, _ = gaia_coord.match_to_catalog_sky(ps1_coord)
    mask = sep < 1.0 * u.arcsec  # keep only good pairs

    gaia_matched = gaia_cat[mask].reset_index(drop=True)
    ps1_matched = ps1_cat.iloc[idx[mask]].reset_index(drop=True)

    # Merge the two DataFrames side-by-side
    both = pd.concat([gaia_matched, ps1_matched.add_suffix('_PS1')], axis=1)
    return both


def calculate_cat_diffs(cat_matched: pd.DataFrame) -> pd.DataFrame:
    """calculated differences between color filters and between catalogs"""
    cat_matched['B_GAIA-B_PS1'] = cat_matched['Bmag'] - cat_matched['Bmag_PS1']
    cat_matched['V_GAIA-V_PS1'] = cat_matched['Vmag'] - cat_matched['Vmag_PS1']
    cat_matched['R_GAIA-R_PS1'] = cat_matched['Rmag'] - cat_matched['Rmag_PS1']
    cat_matched['I_GAIA-I_PS1'] = cat_matched['Imag'] - cat_matched['Imag_PS1']
    cat_matched['B-V'] = cat_matched['Bmag'] - cat_matched['Vmag']
    cat_matched['V-R'] = cat_matched['Vmag'] - cat_matched['Rmag']
    cat_matched['R-I'] = cat_matched['Rmag'] - cat_matched['Imag']
    return cat_matched



class ConeGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Stars cone-search viewer")
        self.geometry("1200x650")

        self.current_df   = None     # full DataFrame behind the grid
        self.filtered_df = None  # view after filter
        self._sort_state  = {}       # column → bool   (True = descending)
        self._build_widgets()

    # ------------------------------------------------------------
    def _build_widgets(self):
        # top-row  –  query controls
        frm = ttk.Frame(self);
        frm.pack(anchor='w', pady=4)

        self.ent_center = self._labeled_entry(frm, "RA, Dec")
        self.ent_center.config(width=25)  # ← after it’s created
        self.ent_rad = self._labeled_entry(frm, 'Radius"')


        # magnitude filter selector  (default = R)
        self.ent_magmin = self._labeled_entry(frm, "mag min")
        self.ent_magmax = self._labeled_entry(frm, "mag max")
        ttk.Label(frm, text="band").pack(side='left', padx=(4, 0))
        self.cmb_magband = ttk.Combobox(frm, state='readonly', width=3,
                                        values = ('B', 'V', 'R', 'I'))
        self.cmb_magband.current(2)  # index 2 → 'R'
        self.cmb_magband.pack(side='left', padx=(0, 4))
        self.ent_solar = self._labeled_entry(frm, "Solar coeff.")
        self.ent_solar.insert(0, "0.2")
        self.ent_magmin.insert(0, "0")
        self.ent_magmax.insert(0, "21")
        self.ent_rad.insert(0, "60.0")

        # ── universal shortcuts for all Entry widgets ─────────────
        self.bind_class("Entry", "<Control-a>",
                        lambda e: (e.widget.select_range(0, "end"),
                                   e.widget.icursor("end"), "break"))
        self.bind_class("Entry", "<Control-A>",
                        lambda e: (e.widget.select_range(0, "end"),
                                   e.widget.icursor("end"), "break"))
        self.bind_class("Entry", "<Command-a>",  # macOS
                        lambda e: (e.widget.select_range(0, "end"),
                                   e.widget.icursor("end"), "break"))

        self.bind_class("Entry", "<Control-z>",
                        lambda e: e.widget.event_generate("<<Undo>>"))
        self.bind_class("Entry", "<Control-y>",
                        lambda e: e.widget.event_generate("<<Redo>>"))
        self.bind_class("Entry", "<Command-z>",
                        lambda e: e.widget.event_generate("<<Undo>>"))
        self.bind_class("Entry", "<Command-y>",
                        lambda e: e.widget.event_generate("<<Redo>>"))

        self.var_solar = tk.BooleanVar()
        ttk.Checkbutton(frm, text="Solar-like", variable=self.var_solar) \
            .pack(side='left', padx=6)

        ttk.Button(frm, text="Query", command=self.query) \
            .pack(side='left', padx=10)
        ttk.Button(frm, text="Save", command=self._save_csv) \
            .pack(side='left')

        # second row – free-text filter
        f2 = ttk.Frame(self);
        f2.pack(anchor='w', pady=(0, 4))
        ttk.Label(f2, text="Filter ").pack(side='left', padx=(4, 0))
        self.ent_filter = ttk.Entry(f2, width=25)
        self.ent_filter.pack(side='left', padx=2)
        self.ent_filter.bind("<Return>", self._apply_filter)
        ttk.Button(f2, text="Apply", command=self._apply_filter) \
            .pack(side='left', padx=10)

        # ─── NEW: RA/Dec format switch ─────────────────────────────
        self.var_hms = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            f2, text="RA°->h / Dec°->dms", variable=self.var_hms,
            command=self._refresh_format
        ).pack(side='left', padx=8)

        # third row – check-boxes for columns
        cb_frame = ttk.Frame(self);
        cb_frame.pack(anchor='w', padx=4, pady=(0, 4))
        self.col_vars = {}
        for name, _ in DISPLAY_COLS:
            var = tk.BooleanVar(value=True)
            self.col_vars[name] = var
            ttk.Checkbutton(cb_frame, text=name, variable=var,
                            command=self._update_columns) \
                .pack(side='left', padx=2)

        # table itself
        self._build_treeview()

    # ------------------------------------------------------------
    # helpers + UI utilities
    # ------------------------------------------------------------

    def _labeled_entry(self, parent, text) -> ttk.Entry:
        fr = ttk.Frame(parent);
        fr.pack(side='left', padx=4)
        ttk.Label(fr, text=text).pack(side='left')
        e = tk.Entry(fr, width=8)
        e.pack(side='left')
        return e
    # ------------------------------------------------------------------
    # Treeview rebuild helpers
    # ------------------------------------------------------------------
    def _selected_columns(self):
        return [n for n, v in self.col_vars.items() if v.get()]

    def _build_treeview(self):
        if hasattr(self, 'tree'):
            self.tree.destroy()

        cols = self._selected_columns()
        self.tree = ttk.Treeview(self, columns=cols, show='headings',
                                 selectmode='extended')
        for c in cols:
            self.tree.heading(c, text=c,
                              command=lambda _c=c: self._sort_by(_c))
            self.tree.column(c, width=110 if "ID" in c else 90,
                             anchor='center')
        self.tree.pack(fill='both', expand=True, padx=4, pady=4)

        vsb = ttk.Scrollbar(self, orient='vertical', command=self.tree.yview)
        self.tree.configure(yscroll=vsb.set)
        vsb.place(relx=1.0, rely=0, relheight=1.0, anchor='ne')

        # context menu (unchanged)
        self.menu = tk.Menu(self, tearoff=0)
        self.menu.add_command(label="Copy Cell", command=self._copy_cell)
        self.menu.add_command(label="Copy Row", command=self._copy_row)
        self.tree.bind('<Button-3>', self._show_context_menu)
        self.tree.bind('<Control-c>', self._copy_selection)
        self.tree.bind('<Control-C>', self._copy_selection)

    def _sort_by(self, col):
        """Sort displayed rows by the clicked column."""
        data = []
        for iid in self.tree.get_children(''):
            val = self.tree.set(iid, col)
            key = float(val) if _is_number(val) else val
            data.append((key, iid))

        descending = self._sort_state.get(col, False)
        data.sort(reverse=not descending)
        for pos, (_, iid) in enumerate(data):
            self.tree.move(iid, '', pos)
        self._sort_state[col] = not descending

    def _update_columns(self):
        self._build_treeview()
        if self.filtered_df is not None:
            self._populate(self.filtered_df)
        elif self.current_df is not None:
            self._populate(self.current_df)

    # ──────────────────────────────────────────────────────────────────────
    #  FILTERING
    # ──────────────────────────────────────────────────────────────────────
    def _apply_filter(self, _event=None):
        """
        Update self.filtered_df and redraw the table.

        • If the filter box is empty, self.filtered_df == self.current_df
        • Otherwise keep only the rows that contain the pattern (case-insensitive)
          anywhere in their string representation.
        """
        if self.current_df is None:
            return

        pattern = self.ent_filter.get().strip()
        if not pattern:
            self.filtered_df = self.current_df
        else:
            pat = re.escape(pattern)
            self.filtered_df = self.current_df[
                self.current_df.apply(
                    lambda row: row.astype(str)
                    .str.contains(pat, case=False, na=False)
                    .any(),
                    axis=1)
            ]
        self._populate(self.filtered_df)

    # ------------------------------------------------------------------
    # context-menu / clipboard helpers (unchanged)
    # ------------------------------------------------------------------
    def _show_context_menu(self, event):
        iid = self.tree.identify_row(event.y)
        col = self.tree.identify_column(event.x)
        if iid:
            self._menu_iid = iid
            self._menu_col = col
            self.menu.tk_popup(event.x_root, event.y_root)

    def _copy_cell(self):
        iid = getattr(self, '_menu_iid', None)
        col = getattr(self, '_menu_col', None)
        if not iid or not col:
            return
        idx = int(col.replace('#', '')) - 1
        vals = self.tree.item(iid, 'values')
        self.clipboard_clear()
        self.clipboard_append(vals[idx])

    def _copy_row(self):
        iid = getattr(self, '_menu_iid', None)
        if not iid:
            return
        vals = self.tree.item(iid, 'values')
        self.clipboard_clear()
        self.clipboard_append('\t'.join(vals))

    def _copy_selection(self, _event=None):
        items = self.tree.selection()
        if not items:
            return
        lines = ['\t'.join(self.tree.item(iid, 'values')) for iid in items]
        self.clipboard_clear()
        self.clipboard_append('\n'.join(lines))

    # ------------------------------------------------------------------
    # CSV saving
    # ------------------------------------------------------------------
    def _save_csv(self):
        if self.current_df is None or self.current_df.empty:
            messagebox.showinfo("Save CSV", "No data to save.")
            return

        filename = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files","*.csv"), ("All files","*.*")]
        )
        if filename:
            # save *every* column; change to self._selected_columns() if you
            # prefer only the visible subset
            self.current_df.to_csv(filename, index=False)
            messagebox.showinfo("Save CSV", f"Table saved to {filename}")

    # ------------------------------------------------------------------
    # query + populate
    # ------------------------------------------------------------------
    def _read_inputs(self):
        try:
            ra, dec = _parse_center(self.ent_center.get())  # CHANGED
            rad = float(self.ent_rad.get())
            mmin = float(self.ent_magmin.get())
            mmax = float(self.ent_magmax.get())
            solar_coef = float(self.ent_solar.get())
        except ValueError as exc:
            raise ValueError(f"Invalid entry: {exc}")
        if mmin > mmax:
            mmin, mmax = mmax, mmin
        return dict(ra=ra, dec=dec, radius=rad,
                    magmin=mmin, magmax=mmax,
                    solar=self.var_solar.get(),
                    solar_coef=solar_coef)

    def query(self):
        try:
            p = self._read_inputs()
        except ValueError as e:
            messagebox.showerror("Input error", str(e))
            return

        self.config(cursor='watch')
        self.update()
        # get the color band for star filtering
        band_col = {'B': 'Bmag', 'V': 'Vmag', 'R': 'Rmag', 'I': 'Imag'} \
            [self.cmb_magband.get()]
        try:
            df_gaia = fetch_catalog_gaia(p['ra'], p['dec'],
                                         p['radius'], p['magmax'], p['magmin'], band=band_col)
            df_ps1 = fetch_catalog_ps1(p['ra'], p['dec'],
                                       p['radius'], p['magmax'], p['magmin'], band=band_col)
            df = calculate_cat_diffs(
                skycoord_match(df_gaia, df_ps1))
        except Exception as exc:
            messagebox.showerror("Query failed", str(exc))
            self.config(cursor='')
            return

        if p['solar']:
            df = solar_mask(df, solar_coef=p['solar_coef'])

        # after the optional solar-masking
        self.current_df = df  # full table
        self.filtered_df = df  # currently shown view
        self.config(cursor='')
        self._populate(df)


    # ------------------------------------------------------------
    # ------------------------------------------------------------------
    def _populate(self, df: pd.DataFrame):
        cols = self._selected_columns()
        self.tree.delete(*self.tree.get_children())

        show_hms = self.var_hms.get()  # ← switch state

        for _, r in df.iterrows():
            row = []
            for c in cols:
                if c.startswith("RA"):
                    row.append(
                        Angle(r['ra_deg'], unit=u.deg)
                        .to_string(unit=u.hourangle if show_hms else u.deg,
                                   sep=':', pad=True, precision=2)
                        if show_hms else f"{r['ra_deg']:.6f}"
                    )
                elif c.startswith("DEC"):
                    row.append(
                        Angle(r['dec_deg'], unit=u.deg)
                        .to_string(unit=u.deg, sep=':', pad=True,
                                   precision=2, alwayssign=True)
                        if show_hms else f"{r['dec_deg']:.6f}"
                    )
                else:
                    row.append(_GETVAL[c](r))
            self.tree.insert("", "end", values=row)

    def _refresh_format(self):
        """Called when the RA/Dec format checkbox is toggled."""
        # update headings
        if 'RA [deg]' in self.tree["columns"]:
            self.tree.heading('RA [deg]',
                              text='RA [hms]' if self.var_hms.get() else 'RA [deg]')
        if 'DEC [deg]' in self.tree["columns"]:
            self.tree.heading('DEC [deg]',
                              text='DEC [dms]' if self.var_hms.get() else 'DEC [deg]')

        # repaint the rows that are *currently visible*
        # → simply re-apply the active filter (if any)
        self._apply_filter()

# ----------------------------------------------------------------------
if __name__ == "__main__":
    ConeGUI().mainloop()

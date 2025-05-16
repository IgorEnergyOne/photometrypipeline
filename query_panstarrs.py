#!/usr/bin/env python3

"""
Interactive cone-search viewer PAN-STARRS DR1 catalog
"""

import tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.vizier import Vizier
from astropy.table import Table, Column

MAX_SOURCES = 50_000

def sources_with(data, condition):
    """
    reject sources based on condition
    input: condition
    return: number of sources left
    """

    data_select = data[condition]
    return data_select

# ------------------------------------------------------------
def fetch_catalog(ra_deg: float, dec_deg: float, radius: float,
                  mag_max: 21, min_mag: 0) -> pd.DataFrame:
    """Generic Vizier cone search returning a pandas.DataFrame."""
    field  = SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg),
                       frame='icrs')
    radius_deg = radius / 60.0

    vquery = Vizier(columns=['objID', 'RAJ2000', 'DEJ2000',
                             'e_RAJ2000', 'e_DEJ2000',
                             'gmag', 'e_gmag',
                             'rmag', 'e_rmag',
                             'imag', 'e_imag',
                             'zmag', 'e_zmag',
                             'ymag', 'e_ymag'],
                    column_filters={"rmag":
                                        ("<{:f}".format(mag_max)),
                                    "rmag":
                                        (">{:f}".format(min_mag))},
                    row_limit=MAX_SOURCES,
                    timeout=300)
    data = vquery.query_region(field,
                               radius=radius_deg * u.deg,
                               catalog="II/349/ps1",
                               cache=False)[0]

    # rename column names using PP conventions
    data.rename_column('objID', 'ident')
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

    g_sdss = (g + 0.013 + 0.145 * (g - r) + 0.019 * (g - r) ** 2)
    gerr_sdss = np.sqrt(e_g ** 2 + 0.008 ** 2)
    r_sdss = (r - 0.001 + 0.004 * (g - r) + 0.007 * (g - r) ** 2)
    rerr_sdss = np.sqrt(e_r ** 2 + 0.004 ** 2)
    i_sdss = (i - 0.005 + 0.011 * (g - r) + 0.010 * (g - r) ** 2)
    ierr_sdss = np.sqrt(e_i ** 2 + 0.004 ** 2)
    z_sdss = (z + 0.013 - 0.039 * (g - r) - 0.012 * (g - r) ** 2)
    zerr_sdss = np.sqrt(e_z ** 2 + 0.01 ** 2)

    data.add_column(Column(data=g_sdss, name='_gmag',
                           unit=u.mag))
    data.add_column(Column(data=gerr_sdss, name='_e_gmag',
                           unit=u.mag))
    data.add_column(Column(data=r_sdss, name='_rmag',
                           unit=u.mag))
    data.add_column(Column(data=rerr_sdss, name='_e_rmag',
                           unit=u.mag))
    data.add_column(Column(data=i_sdss, name='_imag',
                           unit=u.mag))
    data.add_column(Column(data=ierr_sdss, name='_e_imag',
                           unit=u.mag))
    data.add_column(Column(data=z_sdss, name='_zmag',
                           unit=u.mag))
    data.add_column(Column(data=zerr_sdss, name='_e_zmag',
                           unit=u.mag))
    return data.to_pandas()

def solar_mask(df: pd.DataFrame, solar_coef: float = 0.2) -> pd.Series:
    """Solar-analogue colour selection (for cat columns g/r/i)."""
    sol_gr = 0.44  # g-r
    sol_ri = 0.11  # r-i
    stars_solar = sources_with(df, ((df['_gmag'] - df['_rmag']) > sol_gr - solar_coef)
                               & (df['_gmag'] - df['_rmag'] < sol_gr + solar_coef)
                               & (df['_rmag'] - df['_imag'] > sol_ri - solar_coef)
                               & (df['_rmag'] - df['_imag'] < sol_ri + solar_coef))
    stars = stars_solar
    return stars

# ------------------------------------------------------------
class ConeGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("PAN-STARRS Cone-search viewer")
        self.geometry("900x600")
        self._build_widgets()

        # cache (per catalogue)
        self.cache = {}   # key = cat_key -> (centre SkyCoord, radius′, DataFrame)

    # --------------------------------------------------------
    def _build_widgets(self):
        frm = ttk.Frame(self); frm.pack(anchor="w", pady=4)

        self.cat = tk.StringVar(value="PANSTARRS")
        ttk.Label(frm, text="Catalog").pack(side="left", padx=(0,4))
        # ttk.Combobox(frm, textvariable=self.cat,
        #              values=list(VIZIER_CATS.keys()), width=10,
        #              state="readonly").pack(side="left", padx=4)

        self.ent_ra     = self._labeled_entry(frm, "RA °")
        self.ent_dec    = self._labeled_entry(frm, "Dec °")
        self.ent_rad    = self._labeled_entry(frm, "Radius'")
        self.ent_magmin = self._labeled_entry(frm, "mag min")
        self.ent_magmax = self._labeled_entry(frm, "mag max")

        self.var_solar = tk.BooleanVar()
        ttk.Checkbutton(frm, text="Solar-like", variable=self.var_solar
                        ).pack(side="left", padx=6)

        self.ent_solar_coef = self._labeled_entry(frm, "Solar coef")
        self.ent_solar_coef.insert(0, "0.2")  # default value
        self.ent_magmin.insert(0, "0")  # default value
        self.ent_magmax.insert(0, "21")  # default value
        self.ent_rad.insert(0, "1.0")  # default value
        ttk.Button(frm, text="Query", command=self.query
                   ).pack(side="left", padx=10)

        # table
        cols = ("ident","ra_deg","dec_deg","B_mag", "V_mag", "R_mag", "I_mag")
        self.tree = ttk.Treeview(self, columns=cols, show="headings", height=26)
        for c in cols:
            self.tree.heading(c, text=c)
            self.tree.column(c, width=110 if c=="ident" else 90, anchor="center")
        self.tree.pack(fill="both", expand=True, padx=4, pady=4)

        vsb = ttk.Scrollbar(self, orient="vertical",
                            command=self.tree.yview)
        self.tree.configure(yscroll=vsb.set)
        vsb.place(relx=1.0, rely=0, relheight=1.0, anchor="ne")

    # --------------------------------------------------------
    def _labeled_entry(self, parent, text) -> ttk.Entry:
        fr = ttk.Frame(parent); fr.pack(side="left", padx=4)
        ttk.Label(fr, text=text).pack(side="left")
        e = ttk.Entry(fr, width=7); e.pack(side="left")
        return e

    # --------------------------------------------------------
    def _read_inputs(self):
        try:
            ra   = float(self.ent_ra.get())
            dec  = float(self.ent_dec.get())
            rad  = float(self.ent_rad.get())
            mmin = float(self.ent_magmin.get())
            mmax = float(self.ent_magmax.get())
            solar_coef = float(self.ent_solar_coef.get())
        except ValueError:
            raise ValueError("Invalid numeric entry.")
        if mmin > mmax:
            mmin, mmax = mmax, mmin
        return dict(cat_key=self.cat.get(), ra=ra, dec=dec,
                    radius=rad, magmin=mmin, magmax=mmax,
                    solar=self.var_solar.get(), solar_coef=solar_coef)

    # --------------------------------------------------------
    def query(self):
        try:
            p = self._read_inputs()
        except ValueError as e:
            messagebox.showerror("Input error", str(e)); return

        key = p["cat_key"]
        centre = SkyCoord(p["ra"], p["dec"], unit="deg")

        # # decide if we can reuse cache
        # if key in self.cache:
        #     c_solar_coef = float(self.ent_solar_coef.get())
        #     c_centre, c_radius, c_df = self.cache[key]
        #     sep = centre.separation(c_centre).arcminute
        #     if (sep + p["radius"] <= c_radius) \
        #         and (p["magmax"] <= c_df["mag"].max()) \
        #         and (p['solar_coef'] <= c_solar_coef):
        #         # filter cached df
        #         df = c_df[
        #             (centre.separation(SkyCoord(c_df["ra_deg"],
        #                                         c_df["dec_deg"],
        #                                         unit="deg")).arcminute <= p["radius"])
        #             & (c_df["mag"].between(p["magmin"], p["magmax"]))
        #         ]
        #         if p["solar"]:
        #             df = solar_mask(df, solar_coef=c_solar_coef)
        #         self._populate(df); return

        # fresh query
        self.config(cursor="watch"); self.update()
        try:
            df = fetch_catalog(p["ra"], p["dec"],
                               p["radius"], p["magmax"], p["magmin"])
        except Exception as exc:
            messagebox.showerror("Query failed", str(exc))
            self.config(cursor=""); return

        #df = df[df["mag"] >= p["magmin"]]
        if p["solar"]:
            df = solar_mask(df, solar_coef= float(self.ent_solar_coef.get()))

        self.cache[key] = (centre, p["radius"], df)  # save complete set
        self.config(cursor="")
        self._populate(df)

    # --------------------------------------------------------
    def _populate(self, df: pd.DataFrame):
        self.tree.delete(*self.tree.get_children())
        for _, r in df.iterrows():
            ident = r.get("Source") or r.get("objID") or r.get("ident") or ""
            self.tree.insert("", "end", values=(
                int(ident),
                f"{r['ra_deg']:.6f}", f"{r['dec_deg']:.6f}",
                f"{r['Bmag']:.3f}", f"{r['Vmag']:.3f}",
                f"{r['Rmag']:.3f}", f"{r['Imag']:.3f}"
            ))
        df.to_csv(f"selection_{self.cat.get()}.csv", index=False)

# ------------------------------------------------------------
if __name__ == "__main__":
    ConeGUI().mainloop()

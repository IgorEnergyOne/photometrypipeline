# -*- coding: utf-8 -*-
"""
JPL Horizons query helpers for pp_interactive_src.
"""
from __future__ import annotations

import re

import numpy as np
import pandas as pd
from astropy.table import vstack
from astroquery.jplhorizons import Horizons
from astroquery.jplsbdb import SBDB  # noqa: F401 – imported for completeness


def check_object_name(name: str) -> str:
    """Check body name for unwanted symbols and normalise spacing."""
    has_whitespace = bool(re.search(r'\s+', name))
    only_letters = name.isalpha()  # noqa: F841
    has_numbers = any(c.isdigit() for c in name)
    has_letters = re.search(r"[a-zA-Z]", name)
    # check if it is provisional designation with no whitespace
    if has_letters and has_numbers and not has_whitespace:
        name = name[:4] + ' ' + name[4:]
    # check if there is more than one whitespace
    elif has_whitespace:
        name = re.sub(r'\s+', ' ', name)
    return name


def jpl_query_eph(body: str, epochs, location: str, progress_callback=None):
    """Query JPL Horizon system for ephemeris data.

    The query is split into chunks of 50 elements to stay within the API limits.
    """
    step = 50
    end = len(epochs)
    body = check_object_name(body)
    full_ephemerides = []

    if progress_callback:
        progress_callback(0, end)

    for i in range(0, end, step):
        obj = Horizons(id="{}".format(body), location=location, epochs=epochs[i:i + step])
        chunk_ephemerides = obj.ephemerides()
        full_ephemerides = vstack([full_ephemerides, chunk_ephemerides])

        if progress_callback:
            progress_callback(min(i + step, end), end)

    full_ephemerides = full_ephemerides.to_pandas().drop(columns="col0")
    return full_ephemerides


def calc_reduced_mag(app_mag, r, delta):
    """
    Calculate reduced magnitude H(alpha) = V - 5*log10(r*delta).
    *r* and *delta* must be in AU.
    """
    return app_mag - 5 * np.log10(r * delta)


def iterative_lighttime_correction(body: str, epochs, location: str, progress_callback=None):
    """
    Perform iterative lighttime correction:
    1. Query JPL at JD -> get r1, d1.
    2. Calc LT1 = 499.0047838361 * d1 (sec).
    3. JD1_corr = JD - LT1.
    4. Query JPL at JD1_corr -> get r2, d2.
    5. Calc average distance d_avg = (d1 + d2) / 2.
    6. Calc LT2 = 499.0047838361 * d_avg (sec).
    7. JD2_corr = JD - LT2.
    8. Query JPL at JD2_corr -> get final r, d, PAB parameters.

    Returns a DataFrame with corrected time, final distances, and PAB parameters.
    """
    C_AU_S = 499.0047838361  # 1 AU in seconds (1/c)

    jd_0 = np.array(epochs)

    def step_progress(current, total, offset_percent, scale_percent):
        if progress_callback:
            percent = offset_percent + (current / total) * scale_percent
            progress_callback(int(percent), 100)

    # Step 1 – initial query for d1
    df1 = jpl_query_eph(body, jd_0, location,
                        progress_callback=lambda c, t: step_progress(c, t, 0, 33))

    if 'delta' in df1.columns:
        d1 = df1['delta'].values
    else:
        raise ValueError("JPL query did not return 'delta'")

    lt1_days = d1 * C_AU_S / 86400.0
    jd_1 = jd_0 - lt1_days

    # Step 2 – second query for d2
    df2 = jpl_query_eph(body, jd_1, location,
                        progress_callback=lambda c, t: step_progress(c, t, 33, 33))

    if 'delta' in df2.columns:
        d2 = df2['delta'].values
    else:
        raise ValueError("Second JPL query failed to return distance data.")

    d_avg = (d1 + d2) / 2
    lt2_days = d_avg * C_AU_S / 86400.0
    jd_final = jd_0 - lt2_days

    # Step 3 – final query at light-time-corrected epochs
    df_final = jpl_query_eph(body, jd_final, location,
                             progress_callback=lambda c, t: step_progress(c, t, 66, 34))

    df_final['corrected_jd'] = jd_final
    df_final['lighttime_days'] = lt2_days
    df_final['delta_uncorr'] = df1['delta']
    df_final['r_uncorr'] = df1['r']
    df_final['alpha_true_uncorr'] = df1['alpha_true']
    df_final['ObsEclLon_uncorr'] = df1['ObsEclLon']
    df_final['ObsEclLat_uncorr'] = df1['ObsEclLat']

    return df_final


def get_lighttime(jpl_query_data):
    """Get lighttime from JPL query data and convert it to days."""
    lighttime_jd = jpl_query_data['lighttime']
    return lighttime_jd


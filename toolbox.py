"""
Toolbox for the Photometry Pipeline
2016-03-09, mommermiscience@gmail.com
"""
from __future__ import print_function
from __future__ import division

# Photometry Pipeline
# Copyright (C) 2016-2018  Michael Mommert, mommermiscience@gmail.com

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see
# <http://www.gnu.org/licenses/>.

import os
import re
import sys
from pathlib import Path
from urllib.parse import urlparse

import pandas as pd
from astropy.table import vstack
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from astroquery.jplsbdb import SBDB
from astroquery.vizier import Vizier
from astroquery.exceptions import RemoteServiceError

try:
    from past.utils import old_div
except ImportError:
    print('Module future not found. Please install with: pip install future')
    sys.exit()

import math
import time
import requests
import json
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from photutils.aperture import CircularAperture, EllipticalAperture, aperture_photometry


# only import if Python3 is used
if sys.version_info > (3, 0):
    from future import standard_library
    standard_library.install_aliases()
    from builtins import range

import logging




# WORKING WITH VIZIER MIRRORS

def _mirror_origin(mirror):
    """Return the scheme+host of *mirror*, discarding any portal sub-path.
    Examples
    --------
    >>> _mirror_origin("http://vizier.nao.ac.jp/vizier/")
    'http://vizier.nao.ac.jp'
    >>> _mirror_origin("https://vizier.cfa.harvard.edu/vizier/")
    'https://vizier.cfa.harvard.edu'
    """
    parsed = urlparse(mirror)
    return f"{parsed.scheme}://{parsed.netloc}"


def test_mirror(mirror, timeout=6, logger=None):
    """Query a small record from one catalog on a given mirror.

    The test URL is always built from the server *root* (scheme + host),
    because the VizieR CGI lives at ``/viz-bin/VizieR`` on every mirror
    regardless of what portal sub-path the mirror URL may contain.
    """
    cat_id = "I/355/gaiadr3"  # Gaia DR3
    origin = _mirror_origin(mirror)
    test_url = f"{origin}/viz-bin/VizieR?-source={cat_id}&-out.max=1"
    log = logger.info if logger else print
    try:
        start = time.perf_counter()
        r = requests.get(test_url, timeout=timeout)
        elapsed = time.perf_counter() - start
        if r.status_code == 200 and "VizieR" in r.text:
            log(f"Mirror {mirror} responded in {elapsed:.2f}s")
            return elapsed
        else:
            log(f"Mirror {mirror} returned invalid response (code {r.status_code})")
    except requests.RequestException as e:
        log(f"Mirror {mirror} failed: {e}")
    return None


def find_fastest_vizier_mirror(mirrors, logger=None):
    """Return fastest available mirror; returns None if none accessible."""
    results = {}
    log = logger.info if logger else print
    log("Testing Vizier mirrors...\n")
    for m in mirrors:
        log(f"Testing {m} ...")
        latency = test_mirror(m, logger=logger)
        if latency is not None:
            results[m] = latency
            log(f"OK ({latency:.2f}s)")
        else:
            log("failed")

    if not results:
        log("Warning: No Vizier mirrors are currently accessible!")
        return None

    fastest = min(results, key=results.get)
    log(f"Fastest mirror: {fastest} ({results[fastest]:.2f}s)")
    return fastest


def save_cache(mirror, cache_file, logger=None):
    """Save fastest mirror to cache file with timestamp."""
    log = logger.info if logger else print
    try:
        cache_file.write_text(
            json.dumps({"fastest_mirror": mirror, "timestamp": time.time()}, indent=2)
        )
        log(f"Saved mirror {mirror} to cache: {cache_file}")
    except Exception as e:
        if logger:
            logger.warning(f"Failed to save cache file {cache_file}: {e}")
        else:
            print(f"Warning: failed to save cache file {cache_file}: {e}")


def load_cache(path_cache, max_age_hours=12, logger=None):
    """Load cached fastest mirror if not older than max_age_hours."""
    cache_file = Path(path_cache)
    log = logger.info if logger else print
    if not cache_file.exists():
        log(f"No cache file found at {path_cache}")
        return None

    try:
        data = json.loads(cache_file.read_text())
        age_h = (time.time() - data.get("timestamp", 0)) / 3600
        if age_h < max_age_hours:
            log(f"Loaded cached mirror {data.get('fastest_mirror')} (age {age_h:.2f}h)")
            return data.get("fastest_mirror")
        else:
            log(f"Cache expired ({age_h:.1f}h old)")
    except Exception as e:
        if logger:
            logger.warning(f"Failed to load cache file {path_cache}: {e}")
        else:
            print(f"Warning: failed to load cache file {path_cache}: {e}")

    return None


def load_vizier_mirrors(path):
    """Read mirror URLs from a plain-text file (one URL per line, '#' comments).

    Parameters
    ----------
    path : str or Path
        Path to the mirrors file (e.g. ``$PHOTPIPEDIR/setup/vizier_mirrors.dat``).

    Returns
    -------
    list[str]
        List of mirror URL strings.
    """
    mirrors = []
    try:
        for raw in Path(path).read_text().splitlines():
            line = raw.split('#')[0].strip()   # strip inline comments
            if line:
                mirrors.append(line)
    except Exception as e:
        print(f"Warning: could not load VizieR mirrors from {path}: {e}")
    return mirrors


class ResilientVizier:
    """Transparent wrapper around :class:`astroquery.vizier.Vizier` that
    automatically switches to an alternative mirror on connection failure.

    Class-level state (shared across all instances within a process):

    ``_mirrors`` : list[str]
        Mirror URLs loaded from ``setup/vizier_mirrors.dat``.
        Set once in ``_pp_conf.py`` via
        ``ResilientVizier._mirrors = load_vizier_mirrors(...)``.

    ``_cache_path`` : Path
        Path to the JSON mirror-selection cache file.
        Set once in ``_pp_conf.py`` via
        ``ResilientVizier._cache_path = Path(rootpath) / '.vizier_mirror_cache.json'``.

    ``_active_mirror`` : str or None
        Currently active mirror URL, or ``None`` to use astroquery's default.
        Pre-seeded from disk cache in ``_pp_conf.py``.

    ``max_retries`` : int
        Number of times to retry after a mirror switch (default ``1``).
        Override in ``_pp_conf.py`` if needed (e.g. ``ResilientVizier.max_retries = 2``).
    """

    # ── class-level state, populated by _pp_conf.py ──────────────────────
    _mirrors: list = []
    _cache_path: Path = None
    _active_mirror: str = None  # None → astroquery uses its built-in default
    max_retries: int = 1
    # ─────────────────────────────────────────────────────────────────────

    # Exceptions that indicate a network / server problem worth retrying
    _NETWORK_ERRORS = (
        requests.ConnectionError,
        requests.Timeout,
        requests.exceptions.ReadTimeout,
    )

    def __init__(self, **kwargs):
        """Accept the same keyword arguments as :class:`~astroquery.vizier.Vizier`."""
        self._kwargs = kwargs

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _build_vizier(self):
        """Return a real ``Vizier`` instance pointed at the active mirror.

        ``astroquery.Vizier(server=...)`` expects a bare hostname (e.g.
        ``'vizier.nao.ac.jp'``).  We extract that from the full mirror URL so
        that mirrors stored with a portal sub-path (e.g. ``/vizier/``) are
        handled correctly.
        """
        kwargs = dict(self._kwargs)
        if self.__class__._active_mirror:
            # Extract hostname only – astroquery builds the full URL internally.
            # The correct constructor kwarg is 'vizier_server' (not 'server').
            kwargs['vizier_server'] = urlparse(self.__class__._active_mirror).hostname
        return Vizier(**kwargs)

    @classmethod
    def _switch_mirror(cls, logger=None):
        """Reset active mirror, re-probe all mirrors, save winner to cache."""
        log = logger.info if logger else print
        log("ResilientVizier: current mirror failed – probing all mirrors …")
        # Invalidate stale state first so a second concurrent failure doesn't
        # re-use the dead mirror while probing is in progress.
        cls._active_mirror = None

        if not cls._mirrors:
            log("ResilientVizier: no mirror list configured – using astroquery default")
            return

        winner = find_fastest_vizier_mirror(cls._mirrors, logger=logger)
        cls._active_mirror = winner  # may be None if all mirrors are down

        if winner and cls._cache_path:
            save_cache(winner, Path(cls._cache_path), logger=logger)
        elif not winner:
            log("ResilientVizier: WARNING – no VizieR mirror is reachable!")

    # ------------------------------------------------------------------
    # Public query methods (mirror the Vizier API)
    # ------------------------------------------------------------------

    def _run_with_failover(self, method_name, *args, **kwargs):
        """Execute ``method_name`` on a real Vizier, retrying after mirror switch."""
        last_exc = None
        for attempt in range(self.__class__.max_retries + 1):
            try:
                v = self._build_vizier()
                return getattr(v, method_name)(*args, **kwargs)
            except self.__class__._NETWORK_ERRORS as exc:
                last_exc = exc
                logging.warning(
                    f"ResilientVizier: network error on attempt {attempt + 1} "
                    f"({type(exc).__name__}: {exc}) – switching mirror"
                )
                self.__class__._switch_mirror()
            except RemoteServiceError as exc:
                last_exc = exc
                logging.warning(
                    f"ResilientVizier: remote service error on attempt {attempt + 1} "
                    f"({exc}) – switching mirror"
                )
                self.__class__._switch_mirror()
        raise last_exc  # re-raise if all retries exhausted

    def query_region(self, *args, **kwargs):
        return self._run_with_failover('query_region', *args, **kwargs)

    def query_object(self, *args, **kwargs):
        return self._run_with_failover('query_object', *args, **kwargs)

    def query_catalog(self, *args, **kwargs):
        return self._run_with_failover('query_catalog', *args, **kwargs)

# TIME AND DATE

def jd_to_gregorian(jd, is_mjd=False):
    """ convert a julian date into a gregorian data """
    if is_mjd:
        mjd = jd
    else:
        mjd = jd - 2400000.5

    MJD0 = 2400000.5  # 1858 November 17, 00:00:00 hours

    modf = math.modf
    a = int(mjd+MJD0+0.5)
    b = int(old_div((a-1867216.25), 36524.25))
    c = a + b - int(modf(old_div(b, 4))[1]) + 1525

    d = int(old_div((c-122.1), 365.25))
    e = 365*d + int(modf(old_div(d, 4))[1])
    f = int(old_div((c-e), 30.6001))

    day = int(c - e - int(30.6001*f))
    month = int(f - 1 - 12*int(modf(old_div(f, 14))[1]))
    year = int(d - 4715 - int(modf(old_div((7+month), 10))[1]))
    fracofday = mjd - math.floor(mjd)
    hour = int(math.floor(fracofday * 24.0))
    minute = int(math.floor(((fracofday*24.0)-hour)*60.))
    second = int(math.floor(((((fracofday*24.0)-hour)*60.)-minute)*60.))

    return (year, month, day, hour, minute, second)


def dateobs_to_jd(date):
    """convert a string of the format YYYY-MM-DDTHH:MM:SS into a julian
        date; 'T' is used as a separator between date and time
    """
    if 'T' in date:
        date = date.split('T')
    if ' ' in date:
        date = date.split(' ')
    time = date[1].split(':')
    date = date[0].split('-')

    # check if date is yyyy-mm-dd or dd-mm-yyyy
    if len(date[2]) == 4 and len(date[0]) < 3:
        date = date[::-1]

    a = (14 - float(date[1]))//12
    y = float(date[0]) + 4800 - a
    m = float(date[1]) + 12*a - 3
    return float(date[2]) + ((153*m + 2)//5) + 365*y + y//4 - y//100 \
        + y//400 - 32045.5 + old_div(float(time[0]), 24.) + old_div(float(time[1]), 1440.) \
        + old_div(float(time[2]), 86400.)


def jd_to_fractionalyear(jd, is_mjd=False):
    """ convert a julian date into a fractional year, e.g., 2000.123456 """
    if is_mjd:
        jd += 2400000.5
    date = jd_to_gregorian(jd)
    year = date[0]+old_div(date[1], 12.)+old_div(date[2], 365.) + \
        old_div(date[3], 8760.)+old_div(date[4], 525600.)
    return year


def fractionalyear_to_jd(date):
    """ convert a fractional year into a julian date """
    jd_jan1 = dateobs_to_jd('%4d-01-01T00:00:00' % math.floor(date))
    return jd_jan1 + 365*(date-math.floor(date))


# ASTROMATIC tools

def read_scamp_output(xml_filename='scamp_output.xml'):
    """ routine to read in the 'scamp.xml' file """
    raw = open(xml_filename, 'r').readlines()
    headers, hdr_idx, data, data_idx = {}, 0, [], 0
    read_this, idx = False, 0
    while idx < len(raw):
        # read header
        if read_this and raw[idx].find('<FIELD name=') > -1:
            headers[raw[idx][raw[idx].find('<FIELD name')+13:
                             raw[idx].find('" datatype')]] = hdr_idx
            hdr_idx += 1
        # read data
        # new data line
        if read_this and raw[idx].find('<TR>') > -1:
            this_data = []
        # flush data line
        if read_this and raw[idx].find('</TR>') > -1:
            data.append(np.hstack(this_data))
        # actually read data line
        if read_this and raw[idx].find('<TD>') > -1:
            line = raw[idx].replace('</TD>', '<TD>').split('<TD>')
            for item in line:
                if len(item.strip()) > 0 and item.find('\n') == -1:
                    this_data.append(item)
        # control reading
        # activate reading
        if not read_this and \
           raw[idx].find('<TABLE ID="Fields" name="Fields">') > -1:
            read_this = True
        # deactivate reading
        if read_this and raw[idx].find('</TABLEDATA></DATA>') > -1:
            read_this = False
        idx += 1

    # check if data rows have same length as header
    for i in range(len(data)):
        if len(headers) != len(data[i]):
            raise (RuntimeError,
                   ('data and header lists from SCAMP output file have '
                    'different lengths for image %s; do the FITS files have the '
                    'OBJECT keyword populated?') % data[i][headers['Catalog_Name']])
    return (headers, data)


# PP tools

def get_binning(header, obsparam):
    """ derive binning from image header
        use obsparam['binning'] keywords, unless both keywords are set to 1
        return: tuple (binning_x, binning_y)"""

    if (isinstance(obsparam['binning'][0], int) and
            isinstance(obsparam['binning'][1], int)):
        binning_x = obsparam['binning'][0]
        binning_y = obsparam['binning'][1]
    elif '#' in obsparam['binning'][0]:
        if '#blank' in obsparam['binning'][0]:
            binning_x = float(header[obsparam['binning'][0].
                                     split('#')[0]].split()[0])
            binning_y = float(header[obsparam['binning'][1].
                                     split('#')[0]].split()[1])
        elif '#x' in obsparam['binning'][0]:
            binning_x = float(header[obsparam['binning'][0].
                                     split('#')[0]].split('x')[0])
            binning_y = float(header[obsparam['binning'][1].
                                     split('#')[0]].split('x')[1])
        elif '#_' in obsparam['binning'][0]:
            binning_x = float(header[obsparam['binning'][0].
                                     split('#')[0]].split('_')[0])
            binning_y = float(header[obsparam['binning'][1].
                                     split('#')[0]].split('_')[1])
        elif '#CH#' in obsparam['binning'][0]:
            # only for RATIR
            channel = header['INSTRUME'].strip()[1]
            binning_x = float(header[obsparam['binning'][0].
                                     replace('#CH#', channel)])
            binning_y = float(header[obsparam['binning'][1].
                                     replace('#CH#', channel)])
    else:
        binning_x = header[obsparam['binning'][0]]
        binning_y = header[obsparam['binning'][1]]
    return (binning_x, binning_y)


def skycenter(catalogs, ra_key='ra_deg', dec_key='dec_deg', display=True):
    """derive center position and radius from catalogs"""
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    # using percentiles instead of min/max to get better handle
    # on outliers
    percent = 5 # outliers percentile value
    min_ra = min([np.percentile(cat[ra_key], percent)
                  for cat in catalogs])
    max_ra = max([np.percentile(cat[ra_key], 100-percent)
                  for cat in catalogs])
    min_dec = min([np.percentile(cat[dec_key], percent)
                   for cat in catalogs])
    max_dec = max([np.percentile(cat[dec_key], 100-percent)
                   for cat in catalogs])

    ra, dec = (np.rad2deg(np.angle(np.exp(1j*np.deg2rad(min_ra)) +
                                   np.exp(1j*np.deg2rad(max_ra)))),
               np.rad2deg(np.angle(np.exp(1j*np.deg2rad(min_dec)) +
                                   np.exp(1j*np.deg2rad(max_dec)))))
    lower_left = SkyCoord(ra=min_ra, dec=min_dec, frame='icrs', unit='deg')
    upper_right = SkyCoord(ra=max_ra, dec=max_dec, frame='icrs', unit='deg')

    if display:
        # print data to the console
        print('\n#####################################################')
        print("Derive center position and radius from catalogs(toolbox.skycenter):")
        print('Outliers percentile value: {} %'.format(percent))
        print('Derived lower left and upper right coordinates:')
        print(f'lower_left = {lower_left},\nupper_right = {upper_right}')
        print('#####################################################\n')

    rad = lower_left.separation(upper_right).deg / 2

    return ra, dec, rad


# miscellaneous tools

def if_val_in_dict(target_val, dic):
    """check if a value appears in a nested dict structure"""
    result = False
    for key, val in dic.items():
        if type(val) is dict:
            if if_val_in_dict(target_val, val):
                result = True
        elif type(val) is list:
            if target_val in val:
                result = True
        else:
            if target_val == val:
                result = True
    return result


def load_image(fitsfile):
    hdu = fits.open(fitsfile)[0]
    data = hdu.data.astype(float)
    hdr  = hdu.header
    wcs  = WCS(hdr)
    exptime = hdr.get('EXPTIME', hdr.get('EXPOSURE', 1.0))
    dateobs = hdr.get('DATE-OBS')
    return data, hdr, wcs, exptime, dateobs


def find_asteroid_ldac(ldac, ast_coords):
    src = SkyCoord(ldac['ra_deg'], ldac['dec_deg'], unit=u.deg)
    sep = src.separation(ast_coords)
    i_min = np.argmin(sep)
    return i_min, sep[i_min]


def _parse_aperture_param(token: str):
    """Parse a single aperture parameter token.

    Returns
    -------
    float
        Fixed value (e.g. ``"5.0"`` or ``"-18"``).
    (float, float)
        Search range (e.g. ``"2:8"``).
    None
        Optimise automatically (``"auto"``, ``"_"``, ``"?"``, or ``""``).
    """
    s = token.strip().lower()
    if s in ('auto', '_', '?', ''):
        return None
    if ':' in s:
        parts = s.split(':')
        if len(parts) != 2:
            raise ValueError(f"Range token must be 'lo:hi', got '{s}'")
        lo, hi = float(parts[0]), float(parts[1])
        if lo >= hi:
            raise ValueError(f"Range lo ({lo}) must be strictly less than hi ({hi})")
        return (lo, hi)
    return float(s)


def parse_aperture_string(aperture_str):
    """Parse an aperture specification string into a parameter dictionary.

    Each numeric parameter can be:

    * a **fixed value** → stored as ``float``  (e.g. ``"5.0"``)
    * a **search range** → stored as ``(lo, hi)`` tuple  (e.g. ``"2:8"``)
    * ``"auto"`` / ``"_"`` / **omitted** → stored as ``None`` (optimise)

    When **all** parameters are fixed floats the returned ``'mode'`` key is
    ``'manual'`` and the aperture is applied directly.  If any parameter is
    ``None`` or a range the mode is ``'optimal'``, which triggers the
    SNR-grid search in :func:`~optimal_aperture.find_optimal_aperture`.

    An optional ``fwhm=<value>`` token may appear anywhere in the string
    (typically at the end).  When present the parsed FWHM (in pixels) is
    stored under the ``'fwhm'`` key; otherwise ``'fwhm'`` is ``None``.

    Supported shapes
    ----------------
    ``c`` – circular::

        "c 5"            → fixed radius = 5 px  (manual)
        "c"              → optimise radius  (optimal)
        "c auto"         → same as above
        "c 2:8"          → search radius ∈ [2, 8] px  (optimal)
        "c 2:8 fwhm=3.5" → same, FWHM seed = 3.5 px

    ``e`` – elliptical::

        "e 5 2 -18"          → a=5, b=2, θ=-18°  (manual)
        "e"                  → all parameters optimised  (optimal)
        "e auto 2 -18"       → optimise a; b=2, θ=-18° fixed
        "e 2:8 1:4 auto"     → a ∈ [2,8], b ∈ [1,4], θ auto
        "e fwhm=4.0"         → all optimised with FWHM seed = 4.0 px

    ``p`` – pill::

        "p 3 2 -15"              → w=3, h=2, θ=-15°  (manual)
        "p"                      → all parameters optimised  (optimal)
        "p auto auto -15"        → optimise w and h; θ fixed
        "p auto auto -15 fwhm=3" → same with FWHM seed = 3 px

    Parameters
    ----------
    aperture_str : str

    Returns
    -------
    dict
        Always has ``'type'`` (``'circular'``, ``'elliptical'``, or ``'pill'``),
        ``'mode'`` (``'manual'`` or ``'optimal'``), ``'fwhm'`` (``float`` or
        ``None``), plus shape-specific keys whose values are ``float``,
        ``(float, float)``, or ``None``.

    Raises
    ------
    ValueError
        On invalid format or out-of-range fixed values.
    """
    tokens = aperture_str.strip().lower().split()
    if not tokens:
        raise ValueError("Empty aperture string.")

    # ── Extract optional fwhm=<value> token (may appear anywhere) ──────────
    fwhm_val = None
    shape_tokens = []
    for t in tokens:
        m = re.match(r'^fwhm=(\S+)$', t)
        if m:
            try:
                fwhm_val = float(m.group(1))
            except ValueError:
                raise ValueError(f"Invalid fwhm value in aperture string: '{t}'")
            if fwhm_val <= 0:
                raise ValueError(f"fwhm must be > 0, got {fwhm_val}")
        else:
            shape_tokens.append(t)
    tokens = shape_tokens
    if not tokens:
        raise ValueError("Aperture string contains only 'fwhm=...' – shape type is missing.")

    kind = tokens[0]
    rest = tokens[1:]          # may be shorter than expected → missing ≡ None

    def _get(idx):
        return _parse_aperture_param(rest[idx]) if idx < len(rest) else None

    def _check_positive(v, name):
        if isinstance(v, float) and v <= 0:
            raise ValueError(f"{name} must be > 0, got {v}")

    def _check_theta(v, name='theta'):
        if isinstance(v, float) and not (-180.0 <= v <= 180.0):
            raise ValueError(f"{name} must be in [-180, 180], got {v}")

    try:
        if kind == 'c':
            radius = _get(0)
            _check_positive(radius, 'radius')
            result = {'type': 'circular', 'radius': radius}

        elif kind == 'e':
            a, b, theta = _get(0), _get(1), _get(2)
            _check_positive(a, 'a')
            _check_positive(b, 'b')
            _check_theta(theta)
            result = {'type': 'elliptical', 'a': a, 'b': b, 'theta': theta}

        elif kind == 'p':
            width, height, theta = _get(0), _get(1), _get(2)
            _check_positive(width, 'width')
            _check_positive(height, 'height')
            _check_theta(theta)
            result = {'type': 'pill', 'width': width, 'height': height, 'theta': theta}

        else:
            raise ValueError(f"Unknown aperture type '{kind}'. Use 'c', 'e', or 'p'.")

    except ValueError as exc:
        raise ValueError(f"Invalid aperture string '{aperture_str}': {exc}") from exc

    # Determine mode: 'manual' only if every parameter is a plain float
    param_values = [v for k, v in result.items() if k != 'type']
    result['mode'] = 'manual' if all(isinstance(v, float) for v in param_values) else 'optimal'

    # Attach the optional FWHM seed extracted earlier
    result['fwhm'] = fwhm_val

    return result



def lister(path: Path, name_pattern: str, return_type='name', object_type="file") -> list:
    """
    Lists all objects that match the naming pattern by the given path
    :param path: path to the directory
    :param name_pattern: name and format of the files: (e.g. "sim_*.dat")
    :param return_type: "name" - returns only the names of the files in the directory
                        "path" - returns full paths to the files in the directory
    :return: list of paths (or names) of the files that satisfy given conditions
    """
    # get paths to every object with specified name pattern
    objects = sorted(list(Path(path).glob('{}'.format(name_pattern))))
    # check if the object is a file or a directory
    if object_type == "file":
        objects = [obj_path for obj_path in objects if os.path.isfile(obj_path)]
    elif object_type == "dir":
        objects = [obj_path for obj_path in objects if os.path.isdir(obj_path)]
    # get only the names of the objects
    if return_type == 'name':
        objects = sorted([obj.name for obj in objects])
    return objects


def get_full_name(asteroid_id) -> str:
    """gets the funll name of the asteroid from JPL SBDB"""
    jpl_query = SBDB.query("{}".format(asteroid_id), phys=False)
    # check if shortname exists (exists for asteroids with names)
    shortname = jpl_query['object'].get('shortname')
    if shortname:
        return shortname
    else:
        name = jpl_query['object'].get('fullname')
        return name


def init_obs_dict(dict_path: str = os.environ.get('PHOTPIPEDIR') + '/user_scripts/observatories.dat') -> dict:
    """read the data with observatories locations and their codes"""
    obs_dict = {}
    with open(dict_path, 'r') as file:
        obs_file = file.readlines()[1:]
    for obs_site in obs_file:
        code, site = obs_site.strip('\n').split(maxsplit=1)
        obs_dict.update({code: site})
    return obs_dict


def init_mpc_obs_dict(dict_path: str = os.environ.get('PHOTPIPEDIR') + '/user_scripts/observatories_mpc.dat') -> dict:
    """read the data with mpc observatories locations and their codes"""
    obs_dict = {}
    obs_file = pd.read_fwf(dict_path, colspecs=[(0, 4), (4, 14), (14, 23), (23, 33), (33, 200)])
    obs_file['Lat'] = obs_file.apply(lambda x: math.degrees(math.atan2(float(x['sin']), float(x['cos']))), axis=1)
    for obs_site in obs_file.iloc:
        code, long, lat, site_name = obs_site[['Code', 'Long.', 'Lat', 'Name']]
        long_letter = 'E' if float(long) > 0 else 'W'
        lat_sign = "+" if lat > 0 else ""
        full_name = f"{long_letter} {long:.2f} {lat_sign}{lat:.2f} {site_name}"
        obs_dict.update({code: full_name})
    return obs_dict


def detect_phot_system(filter: str) -> str:
    """detects which photometric system is used for the photometry"""
    if filter in ['U', 'B', 'V', 'R', 'I', 'C', 'Clear']:
        return 'Johnson-Cousins'
    elif filter in ['u', 'g', 'r', 'i', 'z']:
        return 'Sloan'
    else:
        return 'Unknown'


def get_fits_header(filename: str) -> dict:
    """gets the header of the fits file"""
    # name of the fits image file
    # open image file
    hdulist = fits.open(filename, mode='update', verify='silentfix',
                        ignore_missing_end=True)
    header = hdulist[0].header
    return header


def get_obsparam(header: dict) -> dict:
    from setup.telescopes import telescope_parameters, instrument_identifiers
    """gets the correct telescope parameters from the pipeline database"""
    instrument_keys = ['TELESCOP', 'INSTRUME', 'PPINSTRU', 'LCAMMOD', 'FPA', 'CAM_NAME',
                       ]
    instruments = []
    for key in instrument_keys:
        if key in header:
            # check the header entry is not empty
            if header[key].strip():
                instruments.append(header[key])
                break
    telescope = instrument_identifiers[instruments[0]]
    obsparam = telescope_parameters[telescope]
    return obsparam


def julian_to_ymd(julian_date):
    """Formats a Julian date as a string in the format "YYYY MON DD.D"""
    month_dict = {'January': 'JAN', 'February': 'FEB', 'March': 'MAR',
                  'April': 'APR', 'May': 'MAY', 'June': 'JUN',
                  'July': 'JUL', 'August': 'AUG', 'September': 'SEP',
                  'October': 'OCT', 'November': 'NOV', 'December': 'DEC'}
    # Create an astropy Time object from the Julian date
    t = Time(julian_date, format='jd', scale='utc')

    # Extract the decimal day
    decimal_day = t.datetime.day + t.datetime.hour / 24 + t.datetime.minute / 1440
    # Format the month name
    month_name = month_dict[t.datetime.strftime('%B')]

    # Format the string with the year, month name, and decimal day
    formatted_date = f"{t.datetime.year} {month_name} {decimal_day:.1f}"

    return formatted_date


def check_object_name(name):
    """check body name for unwanted symbols"""
    # check name
    has_whitespace = bool(re.search(r'\s+', name))
    only_letters = name.isalpha()
    has_numbers = any(c.isdigit() for c in name)
    has_letters = re.search(r"[a-zA-Z]", name)
    # check if it is provisional designation with no whitespace
    if has_letters and has_numbers and not has_whitespace:
        name = name[:4] + ' ' + name[4:]
    # check if there is more than one whitespace
    elif has_whitespace:
        name = re.sub(r'\s+', ' ', name)
    return name


def jpl_query_eph(body, epochs, location):
    """query JPL Horizon system for the data"""
    # query is split into chunks of 50 elements
    step = 50
    # ===============================================
    end = len(epochs)
    body = check_object_name(body)
    full_ephemerides = []
    for i in range(0, end, step):
        obj = Horizons(id="{}".format(body), location=location, epochs=epochs[i:i + step])
        chunk_ephemerides = obj.ephemerides()
        full_ephemerides = vstack([full_ephemerides, chunk_ephemerides])

    full_ephemerides = full_ephemerides.to_pandas().drop(columns="col0")
    return full_ephemerides


def lighttime_to_au(lighttime):
    """
    Converts light-time in minutes to distance in Astronomical Units (AU).
    """
    # Speed of light
    C = 299792458 # m/s
    AU_METERS = 149597870700 # meters to 1 AU
    return (lighttime * 60 * C) / AU_METERS


def calc_reduced_mag(app_mag, r, delta):
    """calculates reduced magnitude for the object (as seen from 1 AU from the Sun and 1 AU from the Earth)"""
    # $$m_red = m - 5 \log_{10}(r \Delta)$$
    red_mag = app_mag - 5 * np.log10(r * delta)
    return red_mag


def get_lighttime(jpl_query_data):
    """gets lighttime from the JPL query data and converts it to days"""
    lighttime_jd = jpl_query_data['lighttime']
    return lighttime_jd

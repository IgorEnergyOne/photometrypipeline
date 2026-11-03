"""
CATALOG_PROVIDERS - query backends for reference catalogs with fallback

Every reference catalog is queried through an ordered list of providers
(_pp_conf.catalog_providers; default: VizieR only). The first provider
that answers wins; a provider that fails is skipped and the next one is
tried. Each non-VizieR provider returns its table with the same raw
column names and units as the corresponding VizieR query, so the
catalog-specific post-processing in catalog._download_catalog() is
identical regardless of where the data came from.
"""

import http.client
import logging
import time

import numpy as np
import requests
import astropy.units as u
from astropy.io import ascii
from astropy.table import MaskedColumn
from astroquery.exceptions import RemoteServiceError
from astroquery.exceptions import TimeoutError as AstroqueryTimeoutError

import _pp_conf
from toolbox import ResilientVizier as Vizier, CatalogUnavailableError

try:
    from pyvo.dal import DALServiceError, DALQueryError
    _PYVO_ERRORS = (DALServiceError, DALQueryError)
except ImportError:
    _PYVO_ERRORS = ()

# exceptions meaning "this provider could not answer" -> try the next one
UNAVAILABLE_ERRORS = (CatalogUnavailableError, RemoteServiceError,
                      AstroqueryTimeoutError, http.client.HTTPException,
                      OSError) + _PYVO_ERRORS  # OSError covers requests

# ---------------------------------------------------------------------------
# VizieR query definitions (column names are original VizieR names)

VIZIER_SPECS = {
    'PANSTARRS': dict(catalog='II/349/ps1',
                      columns=['objID', 'RAJ2000', 'DEJ2000',
                               'e_RAJ2000', 'e_DEJ2000',
                               'gmag', 'e_gmag', 'rmag', 'e_rmag',
                               'imag', 'e_imag', 'zmag', 'e_zmag',
                               'ymag', 'e_ymag'],
                      mag_column='rmag', timeout=300),
    'SkyMapper': dict(catalog='II/358/smss', columns=['all'],
                      mag_column='rPSF', timeout=300),
    'GAIA': dict(catalog='I/345/gaia2',
                 columns=['Source', 'RA_ICRS', 'DE_ICRS',
                          'e_RA_ICRS', 'e_DE_ICRS', 'pmRA', 'pmDE', 'Epoch',
                          'Gmag', 'e_Gmag', 'BPmag', 'e_BPmag',
                          'RPmag', 'e_RPmag', 'Var'],
                 mag_column='Gmag',
                 filters={'Var': '!=VARIABLE'},  # remove variable stars
                 timeout=300),
    'GAIA3': dict(catalog='I/355/gaiadr3',
                  columns=['Source', 'RA_ICRS', 'DE_ICRS',
                           'e_RA_ICRS', 'e_DE_ICRS', 'pmRA', 'pmDE', 'Epoch',
                           'Gmag', 'e_Gmag', 'BPmag', 'e_BPmag',
                           'RPmag', 'e_RPmag', 'VarFlag'],
                  mag_column='Gmag',
                  filters={'VarFlag': '!=VARIABLE'},  # remove variable stars
                  timeout=300),
    'USNO-B1': dict(catalog='I/284/out',
                    columns=['USNO-B1.0', 'RAJ2000', 'DEJ2000',
                             'e_RAJ2000', 'e_DEJ2000', 'R2mag'],
                    mag_column='R2mag', timeout=300),
    'TGAS': dict(catalog='I/337/tgas',
                 columns=['Source', 'RA_ICRS', 'DE_ICRS',
                          'e_RA_ICRS', 'e_DE_ICRS', 'pmRA', 'pmDE',
                          'phot_g_mean_mag'],
                 mag_column='phot_g_mean_mag', timeout=300),
    '2MASS': dict(catalog='II/246/out',
                  columns=['2MASS', 'RAJ2000', 'DEJ2000', 'errMaj',
                           'errMin', 'errPA', 'Jmag', 'e_Jmag',
                           'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag',
                           'Qflg', 'Rflg'],
                  mag_column='Jmag'),
    'URAT-1': dict(catalog='I/329/urat1',
                   columns=['URAT1', 'RAJ2000', 'DEJ2000',
                            'sigm', 'f.mag', 'e_f.mag'],
                   mag_column='f.mag'),
    'APASS9': dict(catalog='II/336/apass9',
                   columns=['recno', 'RAJ2000', 'DEJ2000',
                            'e_RAJ2000', 'e_DEJ2000', 'Vmag', 'e_Vmag',
                            'Bmag', 'e_Bmag', "g'mag", "e_g'mag",
                            "r'mag", "e_r'mag", "i'mag", "e_i'mag"],
                   mag_column='Vmag'),
    'SDSS-R9': dict(catalog='V/139/sdss9',
                    columns=['SDSS9', 'RA_ICRS', 'DE_ICRS',
                             'e_RA_ICRS', 'e_DE_ICRS', 'umag', 'e_umag',
                             'gmag', 'e_gmag', 'rmag', 'e_rmag',
                             'imag', 'e_imag', 'zmag', 'e_zmag'],
                    mag_column='gmag',
                    filters={'mode': '1', 'q_mode': '+'}),
    'SDSS-R16': dict(catalog='V/154/sdss16',
                     columns=['objID', 'RA_ICRS', 'DE_ICRS',
                              'e_RA_ICRS', 'e_DE_ICRS', 'umag', 'e_umag',
                              'gmag', 'e_gmag', 'rmag', 'e_rmag',
                              'imag', 'e_imag', 'zmag', 'e_zmag'],
                     mag_column='gmag'),
}
VIZIER_SPECS['PANSTARRS_KB'] = VIZIER_SPECS['PANSTARRS']


# ---------------------------------------------------------------------------
# helpers

def _set_units(tab, units):
    """attach astropy units {column: unit} in place"""
    for col, unit in units.items():
        if col in tab.colnames:
            tab[col].unit = unit


def _mask_sentinel(tab, columns, sentinel):
    """mask values <= sentinel (archives' 'no data' markers) in place"""
    for col in columns:
        if col not in tab.colnames:
            continue
        data = np.ma.masked_invalid(np.ma.asarray(tab[col], dtype=float))
        data = np.ma.masked_where(data <= sentinel, data)
        tab[col] = MaskedColumn(data, name=col, unit=tab[col].unit)


def _strings(tab, columns):
    """convert object columns (as returned by TAP services) to str"""
    for col in columns:
        if col in tab.colnames:
            tab[col] = tab[col].astype(str)


# ---------------------------------------------------------------------------
# providers: fn(catalogname, field, rad_deg, max_sources, max_mag)
#            -> (Table or None if the service found no sources, origin)

def _vizier(catalogname, field, rad_deg, max_sources, max_mag):
    """query VizieR (with mirror failover) using VIZIER_SPECS"""
    spec = VIZIER_SPECS[catalogname]
    column_filters = {spec['mag_column']: '<{:f}'.format(max_mag)}
    column_filters.update(spec.get('filters', {}))
    kwargs = dict(columns=spec['columns'], column_filters=column_filters,
                  row_limit=max_sources)
    if 'timeout' in spec:
        kwargs['timeout'] = spec['timeout']

    result = Vizier(**kwargs).query_region(field, radius=rad_deg * u.deg,
                                           catalog=spec['catalog'],
                                           cache=False)
    origin = 'VizieR {:s}'.format(spec['catalog'])
    return (result[0] if len(result) > 0 else None), origin


# MAST PS1 cone searches are limited to this radius
MAST_MAX_RADIUS_DEG = 0.5


def _mast_ps1(catalogname, field, rad_deg, max_sources, max_mag):
    """Pan-STARRS1 mean-object table from MAST (raw VizieR II/349 names)"""
    from astroquery.mast import Catalogs

    if rad_deg > MAST_MAX_RADIUS_DEG:
        logging.warning(
            'MAST PS1 radius {:.3f} deg exceeds {:.1f} deg limit; '
            'clamping query radius to {:.1f} deg'.format(
                rad_deg, MAST_MAX_RADIUS_DEG, MAST_MAX_RADIUS_DEG))
        rad_deg = MAST_MAX_RADIUS_DEG

    release = getattr(_pp_conf, 'mast_ps1_release', 'dr2')
    bands = 'grizy'
    columns = ['objID', 'raMean', 'decMean', 'raMeanErr', 'decMeanErr']
    for b in bands:
        columns += ['{:s}MeanPSFMag'.format(b), '{:s}MeanPSFMagErr'.format(b)]

    tab = Catalogs.query_region(field, radius=rad_deg * u.deg,
                                catalog='Panstarrs', data_release=release,
                                table='mean', columns=columns,
                                nDetections=[('gte', 2)],
                                rMeanPSFMag=[('lte', max_mag)],
                                pagesize=max_sources)
    origin = 'MAST PS1 {:s} mean'.format(release)
    if len(tab) == 0:
        return None, origin

    tab = tab[columns]
    tab.rename_columns(
        columns,
        ['objID', 'RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000'] +
        [n for b in bands for n in ('{:s}mag'.format(b),
                                    'e_{:s}mag'.format(b))])
    mags = [n for b in bands for n in ('{:s}mag'.format(b),
                                       'e_{:s}mag'.format(b))]
    _mask_sentinel(tab, mags + ['e_RAJ2000', 'e_DEJ2000'], -999)
    _set_units(tab, {'RAJ2000': u.deg, 'DEJ2000': u.deg,
                     'e_RAJ2000': u.arcsec, 'e_DEJ2000': u.arcsec})
    _set_units(tab, {m: u.mag for m in mags})
    return tab, origin


# zero-point uncertainties VizieR adds to Gaia magnitude errors
_GAIA_SIGMA0 = {'g': 0.0027553202, 'bp': 0.0027901700, 'rp': 0.0037793818}


def _gaia_tap(catalogname, field, rad_deg, max_sources, max_mag):
    """ESA Gaia archive TAP (raw VizieR I/345 or I/355 names)"""
    from astroquery.gaia import Gaia

    table, epoch, varcol = {'GAIA': ('gaiadr2.gaia_source', 2015.5, 'Var'),
                            'GAIA3': ('gaiadr3.gaia_source', 2016.0,
                                      'VarFlag')}[catalogname]
    query = ("SELECT TOP {:d} source_id, ra, dec, ra_error, dec_error, "
             "pmra, pmdec, "
             "phot_g_mean_mag, phot_g_mean_flux_over_error, "
             "phot_bp_mean_mag, phot_bp_mean_flux_over_error, "
             "phot_rp_mean_mag, phot_rp_mean_flux_over_error, "
             "phot_variable_flag FROM {:s} "
             "WHERE 1=CONTAINS(POINT('ICRS', ra, dec), "
             "CIRCLE('ICRS', {:.8f}, {:.8f}, {:.8f})) "
             "AND phot_g_mean_mag < {:f} "
             "AND phot_variable_flag != 'VARIABLE'").format(
                 int(max_sources), table, field.ra.deg, field.dec.deg,
                 rad_deg, max_mag)
    tab = Gaia.launch_job_async(query, verbose=False).get_results()
    origin = 'ESA Gaia archive {:s}'.format(table)
    if len(tab) == 0:
        return None, origin

    for band, name in (('g', 'G'), ('bp', 'BP'), ('rp', 'RP')):
        foe = np.ma.asarray(tab['phot_{:s}_mean_flux_over_error'.format(band)],
                            dtype=float)
        tab['e_{:s}mag'.format(name)] = MaskedColumn(
            np.ma.sqrt((2.5 / np.log(10) / foe) ** 2
                       + _GAIA_SIGMA0[band] ** 2), unit=u.mag)
        tab.rename_column('phot_{:s}_mean_mag'.format(band),
                          '{:s}mag'.format(name))
        tab.remove_column('phot_{:s}_mean_flux_over_error'.format(band))
    tab.rename_columns(['source_id', 'ra', 'dec', 'ra_error', 'dec_error',
                        'pmra', 'pmdec', 'phot_variable_flag'],
                       ['Source', 'RA_ICRS', 'DE_ICRS', 'e_RA_ICRS',
                        'e_DE_ICRS', 'pmRA', 'pmDE', varcol])
    _strings(tab, [varcol])
    tab['Epoch'] = epoch
    tab['Epoch'].unit = u.yr
    _set_units(tab, {'RA_ICRS': u.deg, 'DE_ICRS': u.deg,
                     'e_RA_ICRS': u.mas, 'e_DE_ICRS': u.mas,
                     'Gmag': u.mag, 'BPmag': u.mag, 'RPmag': u.mag})
    return tab, origin


def _gaia_apass9(catalogname, field, rad_deg, max_sources, max_mag):
    """APASS DR9 copy in the ESA Gaia archive (raw VizieR II/336 names)"""
    from astroquery.gaia import Gaia

    query = ("SELECT TOP {:d} recno, raj2000, dej2000, e_ra, e_dec, "
             "vmag, e_vmag, bmag, e_bmag, g_mag, e_g_mag, "
             "r_mag, e_r_mag, i_mag, e_i_mag FROM external.apassdr9 "
             "WHERE 1=CONTAINS(POINT('ICRS', raj2000, dej2000), "
             "CIRCLE('ICRS', {:.8f}, {:.8f}, {:.8f})) "
             "AND vmag < {:f}").format(int(max_sources), field.ra.deg,
                                       field.dec.deg, rad_deg, max_mag)
    tab = Gaia.launch_job_async(query, verbose=False).get_results()
    origin = 'ESA Gaia archive external.apassdr9'
    if len(tab) == 0:
        return None, origin

    tab.rename_columns(['raj2000', 'dej2000', 'e_ra', 'e_dec',
                        'vmag', 'e_vmag', 'bmag', 'e_bmag',
                        'g_mag', 'e_g_mag', 'r_mag', 'e_r_mag',
                        'i_mag', 'e_i_mag'],
                       ['RAJ2000', 'DEJ2000', 'e_RAJ2000', 'e_DEJ2000',
                        'Vmag', 'e_Vmag', 'Bmag', 'e_Bmag',
                        "g'mag", "e_g'mag", "r'mag", "e_r'mag",
                        "i'mag", "e_i'mag"])
    mags = ['Vmag', 'e_Vmag', 'Bmag', 'e_Bmag', "g'mag", "e_g'mag",
            "r'mag", "e_r'mag", "i'mag", "e_i'mag"]
    _set_units(tab, {'RAJ2000': u.deg, 'DEJ2000': u.deg,
                     'e_RAJ2000': u.arcsec, 'e_DEJ2000': u.arcsec})
    _set_units(tab, {m: u.mag for m in mags})
    return tab, origin


def _irsa_2mass(catalogname, field, rad_deg, max_sources, max_mag):
    """2MASS point source catalog from IRSA (raw VizieR II/246 names)"""
    from astroquery.ipac.irsa import Irsa

    query = ("SELECT TOP {:d} designation, ra, dec, err_maj, err_min, "
             "err_ang, j_m, j_msigcom, h_m, h_msigcom, k_m, k_msigcom, "
             "ph_qual, rd_flg FROM fp_psc "
             "WHERE CONTAINS(POINT('ICRS', ra, dec), "
             "CIRCLE('ICRS', {:.8f}, {:.8f}, {:.8f})) = 1 "
             "AND j_m < {:f}").format(int(max_sources), field.ra.deg,
                                      field.dec.deg, rad_deg, max_mag)
    tab = Irsa.query_tap(query).to_table()
    origin = 'IRSA fp_psc'
    if len(tab) == 0:
        return None, origin

    tab.rename_columns(['designation', 'ra', 'dec', 'err_maj', 'err_min',
                        'err_ang', 'j_m', 'j_msigcom', 'h_m', 'h_msigcom',
                        'k_m', 'k_msigcom', 'ph_qual', 'rd_flg'],
                       ['2MASS', 'RAJ2000', 'DEJ2000', 'errMaj', 'errMin',
                        'errPA', 'Jmag', 'e_Jmag', 'Hmag', 'e_Hmag',
                        'Kmag', 'e_Kmag', 'Qflg', 'Rflg'])
    _strings(tab, ['2MASS', 'Qflg', 'Rflg'])
    _set_units(tab, {'RAJ2000': u.deg, 'DEJ2000': u.deg,
                     'errMaj': u.arcsec, 'errMin': u.arcsec, 'errPA': u.deg,
                     'Jmag': u.mag, 'e_Jmag': u.mag, 'Hmag': u.mag,
                     'e_Hmag': u.mag, 'Kmag': u.mag, 'e_Kmag': u.mag})
    return tab, origin


SKYSERVER_URL = ('https://skyserver.sdss.org/dr{:d}/SkyServerWS/'
                 'SearchTools/SqlSearch')


def _skyserver(catalogname, field, rad_deg, max_sources, max_mag):
    """SDSS PhotoPrimary model magnitudes from SkyServer
    (raw VizieR V/139 or V/154 names)"""
    # SkyServer has no DR9 web service; DR9-DR12 share the same imaging
    # and calibration, so DR12 stands in for SDSS-R9
    release, idcol, clean = {'SDSS-R9': (12, 'SDSS9', 'AND p.clean = 1 '),
                             'SDSS-R16': (16, 'objID', '')}[catalogname]
    query = ("SELECT TOP {:d} p.objID, p.ra, p.dec, p.raErr, p.decErr, "
             "p.u, p.err_u, p.g, p.err_g, p.r, p.err_r, "
             "p.i, p.err_i, p.z, p.err_z "
             "FROM dbo.fGetNearbyObjEq({:.8f}, {:.8f}, {:.6f}) AS n "
             "JOIN PhotoPrimary AS p ON p.objID = n.objID "
             "WHERE p.g < {:f} {:s}").format(
                 int(max_sources), field.ra.deg, field.dec.deg,
                 rad_deg * 60, max_mag, clean)
    response = requests.get(SKYSERVER_URL.format(release),
                            params={'cmd': query, 'format': 'csv'},
                            timeout=300)
    response.raise_for_status()
    origin = 'SDSS SkyServer DR{:d} PhotoPrimary'.format(release)

    # a valid answer is '#Table1' followed by CSV; errors come back as text
    lines = response.text.splitlines()
    if not lines or lines[0].strip() != '#Table1':
        raise RemoteServiceError('SkyServer DR{:d}: {:s}'.format(
            release, response.text[:200].strip()))
    if len(lines) < 3:
        return None, origin
    tab = ascii.read(lines[1:], format='csv')

    mags = [n for b in 'ugriz' for n in (b, 'err_' + b)]
    tab.rename_columns(
        ['objID', 'ra', 'dec', 'raErr', 'decErr'] + mags,
        [idcol, 'RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_DE_ICRS'] +
        [n for b in 'ugriz' for n in (b + 'mag', 'e_' + b + 'mag')])
    if idcol == 'SDSS9':
        _strings(tab, [idcol])
    mags = [n for b in 'ugriz' for n in (b + 'mag', 'e_' + b + 'mag')]
    _mask_sentinel(tab, mags, -9999)
    _set_units(tab, {'RA_ICRS': u.deg, 'DE_ICRS': u.deg,
                     'e_RA_ICRS': u.arcsec, 'e_DE_ICRS': u.arcsec})
    _set_units(tab, {m: u.mag for m in mags})
    return tab, origin


PROVIDERS = {'vizier': (_vizier, 'VizieR'),
             'mast': (_mast_ps1, 'MAST'),
             'gaia_tap': (_gaia_tap, 'Gaia archive'),
             'gaia_apass9': (_gaia_apass9, 'Gaia archive'),
             'irsa': (_irsa_2mass, 'IRSA'),
             'skyserver': (_skyserver, 'SkyServer')}


# ---------------------------------------------------------------------------
# provider selection with fallback

# (provider, catalogname) -> time until which the provider is skipped
_down_until = {}


def providers_for(catalogname):
    """ordered provider names configured for catalogname"""
    return list(getattr(_pp_conf, 'catalog_providers', {}).get(
        catalogname, ['vizier']))


def fetch_raw(catalogname, field, rad_deg, max_sources, max_mag,
              display=False):
    """
    query catalogname through its configured providers, in order
    input: catalogname, field (SkyCoord), rad_deg, max_sources, max_mag
    return: (astropy Table with raw VizieR column names, or None if the
             first provider to answer found no sources; origin string)
    raises CatalogUnavailableError if no provider could answer
    """
    cooldown = getattr(_pp_conf, 'provider_cooldown_s', 900)
    names = providers_for(catalogname)
    errors = []

    for idx, name in enumerate(names):
        fn, label = PROVIDERS[name]
        key = (name, catalogname)

        if _down_until.get(key, 0) > time.time():
            errors.append('{:s}: skipped (failed recently)'.format(label))
            logging.info('{:s} skipped for {:s}: failed within the last '
                         '{:d}s'.format(label, catalogname, int(cooldown)))
            continue

        try:
            tab, origin = fn(catalogname, field, rad_deg, max_sources,
                             max_mag)
        except UNAVAILABLE_ERRORS as exc:
            _down_until[key] = time.time() + cooldown
            errors.append('{:s}: {:s}: {}'.format(label,
                                                  type(exc).__name__, exc))
            logging.warning('{:s} unavailable for {:s}: {}'.format(
                label, catalogname, exc))
            if display:
                nxt = (' trying {:s}'.format(PROVIDERS[names[idx + 1]][1])
                       if idx + 1 < len(names) else '')
                print('{:s} unavailable,{:s}'.format(label, nxt),
                      end=' ', flush=True)
            continue

        _down_until.pop(key, None)
        logging.info('{:s} served by {:s}'.format(catalogname, origin))
        return tab, origin

    raise CatalogUnavailableError('; '.join(errors))

#!/usr/bin/env python3

import os
import argparse
from pathlib import Path

from astropy.io import fits
from astropy.table import vstack
from astropy.time import Time
from astroquery.jplhorizons import Horizons
import pandas as pd

# row order in the atlas header
row_order = ['object', 'reference', 'info_observer', 'info_reducer', 'observing_site', 'telescope',
             'detector', 'info_aspect', 'aspect_data', 'columns', 'phot_system', 'relative_phot',
             'reduced_mag', 'lt_corrected', 'info_correction', 'obs_time',
             'zero_time', 'zero_mag', 'time_unit']

# atlas dictionary (transforms keywords to atlas header entries)
atlas_dict = {
    'object': 'OBJECT',
    'reference': 'REFERENCE',
    'info': 'INFORMATION',
    'info_observer': 'INFORMATION',
    'info_reducer': 'INFORMATION',
    'info_aspect': 'INFORMATION',
    'info_correction': 'INFORMATION',
    'observing_site': 'OBSERVING SITE',
    'telescope': 'TELESCOPE',
    'detector': 'DETECTOR',
    'columns': 'COLUMNS',
    'phot_system': 'PHOT. SYSTEM',
    'relative_phot': 'RELATIVE PHOT.',
    'aspect_data': 'ASPECT DATA',
    'reduced_mag': 'REDUCED MAG.',
    'lt_corrected': 'LT CORRECTED',
    'obs_time': 'OBSERVING TIME',
    'zero_time': 'ZERO TIME',
    'zero_mag': 'ZERO MAG',
    'time_unit': 'UNIT OF TIME'
}


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


def init_obs_dict(dict_path: str = os.environ.get('PHOTPIPEDIR') + '/user_scripts/observatories.dat') -> dict:
    """read the data with observatories locations and their codes"""
    obs_dict = {}
    with open(dict_path, 'r') as file:
        obs_file = file.readlines()[1:]
    for obs_site in obs_file:
        code, site = obs_site.strip('\n').split(maxsplit=1)
        obs_dict.update({code: site})
    return obs_dict


def form_atlas_entry(entry: str, atlas_entry_len=16) -> str:
    """forms header entries for atlas files by adding ...: at the end"""
    lendiff = atlas_entry_len - len(entry)
    formed_entry = entry + "." * (lendiff - 1) + ":"
    return formed_entry


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
    """gets the correct telescope parameters from the pipeline database"""
    instrument_keys = ['PPINSTRU', 'LCAMMOD', 'FPA', 'CAM_NAME', 'INSTRUME',
                   'TELESCOP']
    instruments = []
    for key in instrument_keys:
        if key in header:
            instruments.append(header[key])
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


def jpl_query_eph(body, epochs, location):
    """query JPL Horizon system for the data"""
    # query is split into chunks of 50 elements
    step = 50
    # ===============================================
    end = len(epochs)
    full_ephemerides = []
    for i in range(0, end, step):
        obj = Horizons(id="{}".format(body), location=location, epochs=epochs[i:i + step])
        chunk_ephemerides = obj.ephemerides()
        full_ephemerides = vstack([full_ephemerides, chunk_ephemerides])

    full_ephemerides = full_ephemerides.to_pandas().drop(columns="col0")
    return full_ephemerides


def midtime_aspect_data(date: str, target: str, obs_code: str):
    """calculates midtime and aspect data for that moment"""
    # calculate the midpoint of the observation
    # query for the aspect data
    # columns = ['r', 'delta', 'alpha_true', 'PABLon', 'PABLat']
    columns = ['r', 'delta', 'alpha_true', 'ObsEclLon', 'ObsEclLat']
    query_data = jpl_query_eph(body=target,
                               location=obs_code,
                                  epochs=[date])

    query_data = query_data[columns]
    asp = query_data.values.tolist()[0]
    formatted_aspect = f'{asp[0]:.4f} {asp[1]:.4f} {asp[2]:.2f} {asp[3]:.2f} {asp[4]:.2f}'
    return formatted_aspect

def write_atlas(filename_atlas: str, text_atlas: str):
    """writes atlas file"""
    with open(filename_atlas, 'w') as file:
        file.write(text_atlas)

end_atlas = """\n===============------------------------========================
END OF OBJECT   """

def form_atlas(filename_header, filename_photometry):
    """forms atlas file from the resulting pipeline data and fits header"""
    header = get_fits_header(filename_header)
    obsparam = get_obsparam(header)
    obs_dict = init_obs_dict()
    # get photometry data
    photometry_data = pd.read_csv(filename_photometry)
    # zero time of observations (int of julian date)
    zero_time = int(photometry_data['julian_date'].values[0])
    # (observing time) - mean time of observation
    observing_time = (photometry_data['julian_date'].iloc[0]
                      + photometry_data['julian_date'].iloc[-1]) / 2
    # reducing time
    photometry_data['reduc_time'] = photometry_data['julian_date'].values - zero_time
    data_formatted = photometry_data[['reduc_time', 'mag', 'sig']].to_string(header=False, index=False, formatters={
        'reduc_time': '   {:.4f}'.format, 'mag': '{:.4f}'.format, 'sig': '{:.4f}'.format})

    fits_dict = {
        "object": header.get(obsparam.get('object')),
        "observer": 'Observer(s): ' + header.get(obsparam.get('observer')),
        "reference": obsparam.get('reference', 'Krugly et al. in prep.'),
        "info_observer": 'Observer(s): ' + header.get(obsparam['observer']),
        "info_reducer": 'Reducer(s): Yu. Krugly, pipeline',
        "info_aspect": f"aspect data on observing midtime {julian_to_ymd(observing_time)}",
        'aspect_data': midtime_aspect_data(observing_time,
                                           header.get(obsparam.get('object')),
                                           obsparam.get('observatory_code')),
        "info_correction": 'Every point was l.t. corrected, mag reduced, sol. ph. angle corrected to midtime',
        "info_reduc": 'reduced to midtime of the night',
        'observing_site': obsparam.get('observatory_code') + ' ' + obs_dict[obsparam.get('observatory_code')],
        "telescope": obsparam.get('telescope_keyword'),
        "detector": 'CCD',  # header.get(obsparam['detector']),
        "columns": f"#{header.get(obsparam.get('filter'))}.",
        "exptime": obsparam.get('exptime'),
        "airmass": obsparam.get('airmass'),
        "filter": header.get(obsparam.get('filter')),
        "phot_system": detect_phot_system(header.get(obsparam.get('filter'))),
        'relative_phot': 'F',
        'reduced_mag': 'F',
        'lt_corrected': 'F',
        'obs_time': f'{observing_time:.1f} ({julian_to_ymd(observing_time)})',
        'zero_time': f'{zero_time:.1f} ({julian_to_ymd(zero_time)})',
        'zero_mag': 0.0,
        'time_unit': '1 day'
    }

    formatted_atlas = ''
    for idx, row in enumerate(row_order):
        formatted_atlas += f"{form_atlas_entry(atlas_dict.get(row))} {fits_dict.get(row)}\n"
    formatted_atlas += 'DATA:\n'

    formatted_atlas += data_formatted
    formatted_atlas += end_atlas

    return formatted_atlas


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='automated ATLAS file creation')
    parser.add_argument('-fname_header', help='which fits file to use for header',
                        default=None)
    parser.add_argument('-fname_photo', help='which csv photometry file to use for data',
                        default=None)
    parser.add_argument('-out_fname', help='name for the resulting atlas file',
                        default=None)

    args = parser.parse_args()
    filename_header = args.fname_header
    if filename_header is None:
        filename_header = lister(os.getcwd(), '*.fit*', 'path', 'file')[0]
    filename_photo = args.fname_photo
    if filename_photo is None:
        filename_photo = lister(os.getcwd(), '*_.csv', 'path', 'file')[0]
    filename_atlas = args.out_fname
    if filename_atlas is None:
        filename_atlas = str(filename_photo).replace('.csv', '.ATL')

    rootpath = os.environ.get('PHOTPIPEDIR')  # os.environ.get('PHOTPIPEDIR')
    exec(open(rootpath + '/setup/telescopes.py').read())
    text_atlas = form_atlas(filename_header, filename_photo)
    write_atlas(filename_atlas, text_atlas)
    print(f'ATLAS file created: {filename_atlas}')
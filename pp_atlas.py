#!/usr/bin/env python3

import warnings
warnings.filterwarnings("ignore")

import os
import re
import argparse
from pathlib import Path

from astropy.io import fits
from astropy.table import vstack
from astropy.time import Time
from astroquery.jplhorizons import Horizons
from astroquery.jplsbdb import SBDB
import pandas as pd
import numpy as np
import math



# row order in the atlas header
row_order = ['object', 'reference', 'info_observer', 'info_reducer', 'info_add', 'observing_site', 'telescope',
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
    'info_add': 'INFORMATION',
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

def get_full_name(asteroid_id) -> str:
    jpl_query = SBDB.query("{}".format(asteroid_id), phys=False)
    # check if shortname exists (exists for asteroids with names)
    shortname = jpl_query['object'].get('shortname')
    if shortname:
        return shortname
    else:
        name = jpl_query['object'].get('fullname')
        return name


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


def form_atlas_entry(entry: str, atlas_entry_len=15) -> str:
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
    instrument_keys = ['TELESCOP', 'INSTRUME', 'PPINSTRU', 'LCAMMOD', 'FPA', 'CAM_NAME' ,
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
    obsparam = get_obsparam(header) # inst_sigma (reduced sigma) * 2**0.5- cal_sigma - zeropoint_sigma
    obs_dict = init_obs_dict()
    obs_dict_mpc = init_mpc_obs_dict()
    # get photometry data
    photometry_data = pd.read_csv(filename_photometry)
    # zero time of observations (int of julian date - 0.5)
    zero_time = int(photometry_data['julian_date'].values[0]) - 0.5
    # (observing time) - mean time of observation
    observing_time = (photometry_data['julian_date'].iloc[0]
                      + photometry_data['julian_date'].iloc[-1]) / 2
    # get rejected flag
    photometry_data.replace(to_replace=[True, False], value=['!', ''], inplace=True, regex=True)
    # reducing time
    photometry_data['reduc_time'] = photometry_data['julian_date'].values - zero_time
    data_formatted = photometry_data[['rejected', 'reduc_time', 'mag', 'inst_sig', 'sig', 'sextractor_flags']].to_string(header=False,
                                                                                                     index=False,
                                                                                                     formatters={
                                    'rejected': '{:s}'.format, 'reduc_time': '  {:.7f}'.format,
                                    'mag': '{:.4f}'.format, 'inst_sig': '{:.4f}'.format, 'sig': '{:.4f}'.format
                                                                                                     })
    # calculate the median value of sig (error for the object's magnitude and percentiles)
    sig_median = np.median(photometry_data['sig'])
    sig_percentiles = np.percentile(photometry_data['sig'], [16, 84]) - sig_median
    # get the method of photometry analysis that was conducted
    photometry_method = photometry_data['photo_method'].iloc[0]
    # get the photometric catalog used, drop the '_transformed' from the catalog name
    catalog = photometry_data['catalog'].iloc[0].replace('_transformed', '')

    orig_filter = header.get(obsparam.get('filter'))
    # get the photo filter in which the images were processed (not the one in fits header)
    reduc_filter = photometry_data['band'].iloc[0]

    # try to get the name of the observatory by the observatory code
    try:
        observatory = obs_dict[obsparam.get('observatory_code')]
    except KeyError:
        try:
            observatory = obs_dict_mpc[obsparam.get('observatory_code')]
        except KeyError:
            print(f"Observatory code {obsparam.get('observatory_code')} not found in the database")
            observatory = obsparam.get('observatory_code')




    fits_dict = {
        "object": get_full_name(header.get(obsparam.get('object'))),
        "observer": 'Observer(s): ' + header.get(obsparam.get('observer', 'observer'), 'no data'),
        "reference": obsparam.get('reference', 'Krugly et al. in prep.'),
        "info_observer": 'Observer(s): ' + header.get(obsparam.get('observer', 'observer'), 'no data'),
        "info_reducer": 'Reducer(s): Yu. Krugly, pipeline',
        "info_aspect": f"aspect data on observing midtime {julian_to_ymd(observing_time)}",
        "info_add": f"filter: {orig_filter}, band: {reduc_filter}, method: {photometry_method}",
        'aspect_data': midtime_aspect_data(observing_time,
                                           header.get(obsparam.get('object')),
                                           obsparam.get('observatory_code')),
        "info_correction": 'Corrected to midtime',
        "info_reduc": 'reduced to midtime of the night',
        'observing_site': observatory + f", code {obsparam.get('observatory_code')}",
        "telescope": obsparam.get('telescope_keyword') + f", {header.get(obsparam.get('telescope_diameter', 'diameter'), '')}",
        "detector": 'CCD',  # header.get(obsparam['detector']),
        "columns": f"#{reduc_filter}-.f", # new #R-.
        "exptime": obsparam.get('exptime'),
        "airmass": obsparam.get('airmass'),
        "filter": reduc_filter,
        "phot_system": detect_phot_system(header.get(obsparam.get('filter'))),
        'relative_phot': 'F',
        'reduced_mag': 'F',
        'lt_corrected': 'F',
        'obs_time': f'{observing_time:.1f} ({julian_to_ymd(observing_time)})',
        'zero_time': f'{zero_time:.1f} ({julian_to_ymd(zero_time)})',
        'zero_mag': f'{0.0} sigma = {sig_median:.4f} {sig_percentiles[0]:.4f} +{sig_percentiles[1]:.4f}, catalog: {catalog}',
        'time_unit': '1 day'
    }

    formatted_atlas = ''
    for idx, row in enumerate(row_order):
        formatted_atlas += f"{form_atlas_entry(atlas_dict.get(row))} {fits_dict.get(row)}\n"
    formatted_atlas += 'DATA:\n'

    formatted_atlas += data_formatted
    formatted_atlas += end_atlas

    return formatted_atlas

def combine_atlas(fname_batch: str, fname_out: str):
    core_path = os.getcwd()
    # read paths from the file
    with open(fname_batch, 'r') as file:
        paths = file.readlines()
    all_atlas = []
    # check every path if correct
    for path in paths:
        path = path.replace('\n', '')
        path = Path(path)
        # check if the path is absolute or not
        if path.is_absolute():
            atlas_path = path
        else:
            atlas_path = Path(core_path) / path

        # check if it is a directory
        if atlas_path.is_dir():
            # look for the atlas file in the directory
            atlas_files = list(atlas_path.glob('*.ATL'))
            if len(atlas_files) == 0:
                raise FileNotFoundError('no atlas files found in the directory: %s' % atlas_path)
            elif len(atlas_files) > 1:
                raise ValueError('multiple atlas files found in the directory: %s' % atlas_path)
            atlas_path = atlas_files[0]
        # check if the file exists
        if not atlas_path.is_file():
            raise FileNotFoundError('file does not exist: %s' % atlas_path)

        with open(atlas_path, 'r') as file:
            atlas = file.readlines()
            # remove "END OF OBJECT" character
            atlas = atlas[:-1]
            atlas = ''.join(atlas)
            all_atlas.append(atlas)

    atlas_whole = ''.join(all_atlas)
    atlas_whole += 'END OF OBJECT'
    with open(fname_out, 'w') as file:
        file.write(atlas_whole)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='automated ATLAS file creation')
    parser.add_argument('-fname_header', help='which fits file to use for header',
                        default=None)
    parser.add_argument('-fname_photo', help='which csv photometry file to use for data',
                        default=None)
    parser.add_argument('-fname_out', help='name for the resulting atlas file',
                        default=None)
    parser.add_argument('-combine', help='combine multiple ATLAS files into one by the paths from file',
                        default=None)

    args = parser.parse_args()


    filename_header = args.fname_header
    filename_atlas = args.fname_out

    if args.combine is None:
        rootpath = os.environ.get('PHOTPIPEDIR')
        if filename_header is None:
            filename_header = lister(os.getcwd(), '*.fit*', 'path', 'file')[0]
            # check if path is a directory
        elif os.path.isdir(filename_header):
            filename_header = lister(filename_header, '*.fit*', 'path', 'file')[0]
        filename_photo = args.fname_photo
        if filename_photo is None:
            filename_photo = lister(os.getcwd(), '*_.csv', 'path', 'file')[0]
        elif os.path.isdir(filename_photo):
            filename_photo = lister(filename_photo, '*_.csv', 'path', 'file')[0]
        if filename_atlas is None:
            filename_atlas = str(os.path.basename(filename_photo)).replace('.csv', '.ATL')
        exec(open(rootpath + '/setup/telescopes.py').read())
        text_atlas = form_atlas(filename_header, filename_photo)
        # transform text from UNIX to DOS format
        text_atlas = text_atlas.replace('\n', '\r\n')
        print('\n#-----------------------\nresulting ATLAS file:\n\n' + text_atlas)
        write_atlas(filename_atlas, text_atlas)
        print(f'ATLAS file created: {filename_atlas}')
    else:
        if filename_atlas is None:
            filename_atlas = "combined_atlas.ATL"
        combine_atlas(args.combine, filename_atlas)
        print(f'combine ATLAS file created: {filename_atlas}')

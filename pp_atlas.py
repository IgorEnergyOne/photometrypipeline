#!/usr/bin/env python3

import warnings

warnings.filterwarnings("ignore")

import os
import argparse
from pathlib import Path

import pandas as pd
import numpy as np
import toolbox

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


def form_atlas_entry(entry: str, atlas_entry_len=15) -> str:
    """forms header entries for atlas files by adding ...: at the end"""
    lendiff = atlas_entry_len - len(entry)
    formed_entry = entry + "." * (lendiff - 1) + ":"
    return formed_entry


def midtime_aspect_data(date: str, target: str, obs_code: str):
    """calculates midtime and aspect data for that moment"""
    # calculate the midpoint of the observation
    # query for the aspect data
    # columns = ['r', 'delta', 'alpha_true', 'PABLon', 'PABLat']
    columns = ['r', 'delta', 'alpha_true', 'ObsEclLon', 'ObsEclLat']
    query_data = toolbox.jpl_query_eph(body=target,
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


def form_atlas(filename_header, filename_photometry, use_reduced_mag=False, use_lt_corrected=False):
    """forms atlas file from the resulting pipeline data and fits header"""
    header = toolbox.get_fits_header(filename_header)
    obsparam = toolbox.get_obsparam(header)  # inst_sigma (reduced sigma) * 2**0.5- cal_sigma - zeropoint_sigma
    obs_dict = toolbox.init_obs_dict()
    obs_dict_mpc = toolbox.init_mpc_obs_dict()
    # get photometry data
    if type(filename_photometry) == pd.DataFrame:
        photometry_data = filename_photometry
    else:
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
    data_formatted = photometry_data[
        ['rejected', 'reduc_time', 'mag', 'inst_sig', 'sig', 'sextractor_flags']].to_string(header=False,
                                                                                            index=False,
                                                                                            formatters={
                                                                                                'rejected': '{:s}'.format,
                                                                                                'reduc_time': '  {:.7f}'.format,
                                                                                                'mag': '{:.4f}'.format,
                                                                                                'inst_sig': '{:.4f}'.format,
                                                                                                'sig': '{:.4f}'.format
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
        "object": toolbox.get_full_name(header.get(obsparam.get('object'))),
        "observer": 'Observer(s): ' + header.get(obsparam.get('observer', 'observer'), 'no data'),
        "reference": obsparam.get('reference', 'Krugly et al. in prep.'),
        "info_observer": 'Observer(s): ' + header.get(obsparam.get('observer', 'observer'), 'no data'),
        "info_reducer": 'Reducer(s): Yu. Krugly, pipeline',
        "info_aspect": f"aspect data on observing midtime {toolbox.julian_to_ymd(observing_time)}",
        "info_add": f"filter: {orig_filter}, band: {reduc_filter}, method: {photometry_method}",
        'aspect_data': midtime_aspect_data(observing_time,
                                           header.get(obsparam.get('object')),
                                           obsparam.get('observatory_code')),
        "info_correction": 'Corrected to midtime',
        "info_reduc": 'reduced to midtime of the night',
        'observing_site': observatory + f", code {obsparam.get('observatory_code')}",
        "telescope": obsparam.get(
            'telescope_keyword') + f", {header.get(obsparam.get('telescope_diameter', 'diameter'), '')}",
        "detector": 'CCD',  # header.get(obsparam['detector']),
        "columns": f"#{reduc_filter}.-f",  # new #R-.
        "exptime": obsparam.get('exptime'),
        "airmass": obsparam.get('airmass'),
        "filter": reduc_filter,
        "phot_system": toolbox.detect_phot_system(header.get(obsparam.get('filter'))),
        'relative_phot': 'F',
        'reduced_mag': 'T' if use_reduced_mag else 'F',
        'lt_corrected': 'T' if use_lt_corrected else 'F',
        'obs_time': f'{observing_time:.1f} ({toolbox.julian_to_ymd(observing_time)})',
        'zero_time': f'{zero_time:.1f} ({toolbox.julian_to_ymd(zero_time)})',
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


def combine_atlas(input_paths, fname_out: str):
    """
    Combines multiple ATLAS files into one.

    Args:
        input_paths: Either a path to a file containing a list of paths, or a list of directory paths
        fname_out: Output file path for the combined ATLAS file
    """
    core_path = os.getcwd()
    all_atlas = []

    # Determine if input_paths is a file (original behavior) or a list of directories (new behavior)
    if isinstance(input_paths, str):
        # Original behavior: read paths from a file
        with open(input_paths, 'r') as file:
            paths = file.readlines()
    else:
        # New behavior: input_paths is already a list of directory paths
        paths = input_paths
    # Process each path
    for path in paths:
        path = path.replace('\n', '')
        path = Path(path)
        # check if the path is absolute or not
        if not path.is_absolute():
            path = Path(core_path) / path
        # check if it is a directory
        if path.is_dir():
            # look for the atlas file in the directory
            atlas_files = list(path.glob('*.ATL'))
            if not atlas_files:
                raise FileNotFoundError('no atlas files found in the directory: %s' % path)
            if len(atlas_files) > 1:
                raise ValueError('multiple atlas files found in the directory: %s' % path)
            path = atlas_files[0]

        # check if the file exists
        if not path.is_file():
            raise FileNotFoundError('file does not exist: %s' % path)

        with open(path, 'r') as file:
            atlas = file.readlines()
            # remove "END OF OBJECT" character
            atlas = atlas[:-1]
            atlas = ''.join(atlas)
            # transform text from UNIX to DOS format
            atlas = atlas.replace('\n', '\r\n')
            all_atlas.append(atlas)

        atlas_whole = ''.join(all_atlas)
        atlas_whole += 'END OF OBJECT'
        with open(fname_out, 'w') as file:
            file.write(atlas_whole)

    print(f'Successfully combined {len(all_atlas)} ATLAS files into {fname_out}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='automated ATLAS file creation')
    parser.add_argument('-fname_header', help='which fits file to use for header',
                        default=None)
    parser.add_argument('-fname_photo', help='which csv photometry file to use for data',
                        default=None)
    parser.add_argument('-fname_out', help='name for the resulting atlas file',
                        default=None)
    parser.add_argument('-use_reduced_mag', action='store_true', default=False,
                        help='mark the ATLAS file as containing reduced magnitudes (sets REDUCED MAG.: T)')
    parser.add_argument('-use_lt_corrected', action='store_true', default=False,
                        help='mark the ATLAS file as lighttime-corrected (sets LT CORRECTED: T)')
    parser.add_argument('-combine', nargs='+',
                        help='combine multiple ATLAS files. Can be a single file containing a list of paths, '
                             'or a list of directories containing .ATL files')

    args = parser.parse_args()

    filename_header = args.fname_header
    filename_atlas = args.fname_out

    if args.combine is None:
        rootpath = os.environ.get('PHOTPIPEDIR')
        if filename_header is None:
            filename_header = toolbox.lister(os.getcwd(), '*.fit*', 'path', 'file')[0]
            # check if path is a directory
        elif os.path.isdir(filename_header):
            filename_header = toolbox.lister(filename_header, '*.fit*', 'path', 'file')[0]
        filename_photo = args.fname_photo
        if filename_photo is None:
            filename_photo = toolbox.lister(os.getcwd(), '*_.csv', 'path', 'file')[0]
        elif os.path.isdir(filename_photo):
            filename_photo = toolbox.lister(filename_photo, '*_.csv', 'path', 'file')[0]
        if filename_atlas is None:
            filename_atlas = str(os.path.basename(filename_photo)).replace('.csv', '.ATL')
        exec(open(rootpath + '/setup/telescopes.py').read())
        text_atlas = form_atlas(filename_header, filename_photo,
                                use_reduced_mag=args.use_reduced_mag,
                                use_lt_corrected=args.use_lt_corrected)
        # transform text from UNIX to DOS format
        text_atlas = text_atlas.replace('\n', '\r\n')
        print('\n#-----------------------\nresulting ATLAS file:\n\n' + text_atlas)
        write_atlas(filename_atlas, text_atlas)
        print(f'ATLAS file created: {filename_atlas}')
    else:
        if filename_atlas is None:
            filename_atlas = "combined_atlas.ATL"

        # If only one path is provided and it's a file, use original behavior
        if len(args.combine) == 1 and os.path.isfile(args.combine[0]):
            combine_atlas(args.combine[0], filename_atlas)
        else:
            # Treat as list of directories
            combine_atlas(args.combine, filename_atlas)

        print(f'Combined ATLAS file created: {filename_atlas}')

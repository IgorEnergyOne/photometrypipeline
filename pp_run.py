#!/usr/bin/env python3

""" PP_RUN - wrapper for automated data analysis
    v1.0: 2016-02-10, mommermiscience@gmail.com
"""
from __future__ import print_function

# Photometry Pipeline
# Copyright (C) 2016-2018 Michael Mommert, mommermiscience@gmail.com

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

import re
import os
import gc
import sys
import yaml
try:
    import numpy as np
except ImportError:
    print('Module numpy not found. Please install with: pip install numpy')
    sys.exit()
import logging
import argparse
try:
    from astropy.io import fits
except ImportError:
    print('Module astropy not found. Please install with: pip install astropy')
    sys.exit()

# only import if Python3 is used
if sys.version_info > (3, 0):
    from builtins import str
    from builtins import range

# pipeline-specific modules
import _pp_conf
from catalog import *
import pp_prepare
import pp_extract
import pp_register
import pp_photometry
import pp_calibrate
import pp_distill
from diagnostics import registration as diag

# setup logging
logging.basicConfig(filename=_pp_conf.log_filename,
                    level=_pp_conf.log_level,
                    format=_pp_conf.log_formatline,
                    datefmt=_pp_conf.log_datefmt)


def run_the_pipeline(filenames, man_targetname, man_filtername,
                     fixed_aprad, source_tolerance, solar,
                     rerun_registration, asteroids, keep_wcs):
    """
    wrapper to run the photometry pipeline
    """
    # increment pp process idx
    _pp_conf.pp_process_idx += 1

    # # reset diagnostics for this data set
    # _pp_conf.dataroot, _pp_conf.diagroot, \
    #     _pp_conf.index_filename, _pp_conf.reg_filename, _pp_conf.cal_filename, \
    #     _pp_conf.res_filename = _pp_conf.setup_diagnostics()

    # setup logging again (might be a different directory)
    logging.basicConfig(filename='LOG',
                        level=_pp_conf.log_level,
                        format=_pp_conf.log_formatline,
                        datefmt=_pp_conf.log_datefmt)

    # read telescope information from fits headers
    # check that they are the same for all images
    logging.info('##### new pipeline process in {:s} #####'.format(
        os.getcwd()))
    logging.info(('check for same telescope/instrument for %d ' +
                  'frames') % len(filenames))
    instruments = []
    for idx, filename in enumerate(filenames):
        try:
            hdulist = fits.open(filename, ignore_missing_end=True)
        except IOError:
            logging.error('cannot open file %s' % filename)
            print('ERROR: cannot open file %s' % filename)
            filenames.pop(idx)
            continue

        header = hdulist[0].header
        for key in _pp_conf.instrument_keys:
            if key in header:
                instruments.append(header[key])
                break

    if len(filenames) == 0:
        raise IOError('cannot find any data...')

    if len(instruments) == 0:
        raise KeyError('cannot identify telescope/instrument; please update' +
                       '_pp_conf.instrument_keys accordingly')

    # check if there is only one unique instrument
    if len(set(instruments)) > 1:
        print('ERROR: multiple instruments used in dataset: %s' %
              str(set(instruments)))
        logging.error('multiple instruments used in dataset: %s' %
                      str(set(instruments)))
        for i in range(len(filenames)):
            logging.error('%s %s' % (filenames[i], instruments[i]))
        sys.exit()

    telescope = _pp_conf.instrument_identifiers[instruments[0]]
    obsparam = _pp_conf.telescope_parameters[telescope]
    logging.info('%d %s frames identified' % (len(filenames), telescope))

    # read filter information from fits headers
    # check that they are the same for all images
    logging.info(('check for same filter for %d ' +
                  'frames') % len(filenames))
    filters = []
    filenames_filter = []

    # if to select specific filter from series with multiple filters
    if man_filtername:
        for idx, filename in enumerate(filenames):
            try:
                hdulist = fits.open(filename, ignore_missing_end=True)
            except IOError:
                logging.error('cannot open file %s' % filename)
                print('ERROR: cannot open file %s' % filename)
                filenames.pop(idx)
                continue

            header = hdulist[0].header
            if man_filtername == header[obsparam['filter']]:
                filters.append(header[obsparam['filter']])
                filenames_filter.append(filename)

    # default regime (processing all available data)
    else:
        for idx, filename in enumerate(filenames):
            try:
                hdulist = fits.open(filename, ignore_missing_end=True)
            except IOError:
                logging.error('cannot open file %s' % filename)
                print('ERROR: cannot open file %s' % filename)
                filenames.pop(idx)
                continue

            header = hdulist[0].header
            filters.append(header[obsparam['filter']])

    if len(filters) == 0:
        raise KeyError('cannot identify filter; please update' +
                       'setup/telescopes.py accordingly')

    if len(set(filters)) > 1:
        print('ERROR: multiple filters used in dataset: %s' % str(set(filters)))
        logging.error('multiple filters used in dataset: %s' %
                      str(set(filters)))
        for i in range(len(filenames)):
            logging.error('%s %s' % (filenames[i], filters[i]))
        sys.exit()

    if man_filtername is None:
        try:
            filtername = obsparam['filter_translations'][filters[0]]
        except KeyError:
            print(('Cannot translate filter name (%s); please adjust ' +
                   'keyword "filter_translations" for %s in ' +
                   'setup/telescopes.py') % (filters[0], telescope))
            logging.error(('Cannot translate filter name (%s); please adjust ' +
                           'keyword "filter_translations" for %s in ' +
                           'setup/telescopes.py') % (filters[0], telescope))
            return None
    else:
        filtername = man_filtername
        filenames = filenames_filter
    logging.info('%d %s frames identified' % (len(filenames), filtername))

    print('run photometry pipeline on %d %s %s frames' %
          (len(filenames), telescope, filtername))

    change_header = {}
    if man_targetname is not None:
        change_header['OBJECT'] = man_targetname

    print('\n-----preparing images (pp_prepare.prepare)')
    # prepare fits files for photometry pipeline
    preparation = pp_prepare.prepare(filenames, obsparam,
                                     change_header,
                                     diagnostics=True, display=True,
                                     keep_wcs=keep_wcs)

    # run wcs registration
    if not keep_wcs:
        # default sextractor/scamp parameters
        if auto:
            snr, source_minarea = obsparam['source_snr'], obsparam['source_minarea']
            aprad = obsparam['aprad_default']
            src_tol = obsparam['source_tolerance']
            mancat = None
            max_rad = 5.0
            source_maxarea = 0
        # parameters from config file
        else:
            snr, source_minarea, source_maxarea = reg_snr, reg_minarea, reg_maxarea
            aprad = reg_aprad
            src_tol = source_tolerance
            mancat = reg_cat
            max_rad = reg_max_rad
        registration_run_number = 0
        while True:
            print('\n----- run image registration (pp_register.register)')
            registration = pp_register.register(filenames=filenames, telescope=telescope,
                                                sex_snr=snr,
                                                source_minarea=source_minarea,
                                                source_maxarea=source_maxarea,
                                                aprad=aprad,
                                                mancat=mancat,
                                                obsparam=obsparam,
                                                source_tolerance=src_tol,
                                                nodeblending=True,
                                                max_rad=max_rad,
                                                display=True,
                                                diagnostics=True)
            if len(registration['badfits']) == len(filenames):
                summary_message = "<FONT COLOR=\"red\">registration failed</FONT>"
            elif len(registration['goodfits']) == len(filenames):
                summary_message = "<FONT COLOR=\"green\">all images registered" + \
                    "</FONT>; "
                break
            else:
                summary_message = "<FONT COLOR=\"orange\">registration failed for " + \
                    ("%d/%d images</FONT>; " %
                     (len(registration['badfits']),
                      len(filenames)))
            # break from loop if maximum number of iterations (2) achieved
            registration_run_number += 1
            if registration_run_number == 2:
                break

        # add information to summary website, if requested
        if _pp_conf.use_diagnostics_summary:
            diag.insert_into_summary(summary_message)

        # in case not all image were registered successfully
        filenames = registration['goodfits']

    # stop here if registration failed for all images
    if len(filenames) == 0:
        logging.info('Nothing else to do for this image set')
        print('Nothing else to do for this image set')
        diag.abort('pp_registration')
        return None

    # run photometry (curve-of-growth analysis)
    if auto:
        snr, source_minarea = 1.5, obsparam['source_minarea']
        source_maxarea = obsparam.get('source_maxarea', 0)
        background_only = False
        target_only = False
    else:
        snr, source_minarea, source_maxarea = photo_snr, photo_minarea, photo_maxarea
        background_only = photo_background
        target_only = photo_target
    if fixed_aprad == 0:
        aprad = None  # force curve-of-growth analysis
    else:
        aprad = fixed_aprad  # skip curve_of_growth analysis

    print('\n----- derive optimum photometry aperture (pp_photometry.photometry)\n')
    phot = pp_photometry.photometry(filenames, snr, source_minarea, source_maxarea, aprad,
                                    man_targetname, background_only,
                                    target_only,
                                    telescope, obsparam, display=True,
                                    diagnostics=True)

    # data went through curve-of-growth analysis
    if phot is not None:
        summary_message = ("<FONT COLOR=\"green\">aprad = %5.1f px, " +
                           "</FONT>") % phot['optimum_aprad']
        if phot['n_target'] > 0:
            summary_message += "<FONT COLOR=\"green\">based on target and " + \
                               "background</FONT>; "
        else:
            summary_message += "<FONT COLOR=\"orange\">based on background " + \
                               "only </FONT>; "
    # a fixed aperture radius has been used
    else:
        if _pp_conf.photmode == 'APER':
            summary_message += "using a fixed aperture radius of %.1f px;" % aprad

    # add information to summary website, if requested
    if _pp_conf.use_diagnostics_summary:
        diag.insert_into_summary(summary_message)

    # run photometric calibration
    if auto:
        minstars = _pp_conf.minstars
        maxstars = 300
        manualcatalog = None
        maxflag = 3
    # taking pp_calibrate parameters from configfile
    else:
        minstars = cal_minstars
        maxstars = cal_maxstars
        manualcatalog = cal_catalog
        maxflag = cal_maxflag

    print('\n----- run photometric calibration (pp_calibrate.calibrate)\n')

    while True:
        calibration = pp_calibrate.calibrate(filenames=filenames,
                                             minstars=minstars,
                                             maxstars=maxstars,
                                             manfilter=filtername,
                                             manualcatalog=manualcatalog,
                                             obsparam=obsparam,
                                             maxflag=maxflag,
                                             solar=solar,
                                             display=True,
                                             diagnostics=True)

        try:
            zps = [frame['zp'] for frame in calibration['zeropoints']]
            zp_errs = [frame['zp_sig']
                       for frame in calibration['zeropoints']]

            # rerun calibration
            if solar and any(np.isnan(zps)):
                logging.warning(('Photometric calibration failed for one '
                                 'or more frames; re-try without the '
                                 'solar option'))
                print(('Warning: Photometric calibration '
                       'failed for one or more frames; '
                       're-try without the -solar option'))
                solar = False
                continue

            if calibration['ref_cat'] is not None:
                refcatname = calibration['ref_cat'].catalogname
            else:
                refcatname = 'instrumental magnitudes'
                summary_message = "<FONT COLOR=\"green\">average zeropoint = " + \
                    ("%5.2f+-%5.2f using %s</FONT>; " %
                     (np.average(zps),
                      np.average(zp_errs),
                      refcatname))
        except TypeError:
            summary_message = "<FONT COLOR=\"red\">no phot. calibration</FONT>; "
        break

    # add information to summary website, if requested
    if _pp_conf.use_diagnostics_summary:
        diag.insert_into_summary(summary_message)

    # distill photometry results
    print('\n----- distill photometry results (pp_distill.distill)\n')
    distillate = pp_distill.distill(catalogs=calibration['catalogs'],
                                    man_targetname=man_targetname,
                                    offset=[0, 0],
                                    fixed_targets_file=None,
                                    posfile=None,
                                    rejectionfilter=rejectionfilter,
                                    asteroids=asteroids,
                                    display=True, diagnostics=True)

    targets = np.array(list(distillate['targetnames'].keys()))
    try:
        target = targets[targets != 'control_star'][0]
        mags = [frame[7] for frame in distillate[target]]
        summary_message = ("average target brightness and std: " +
                           "%5.2f+-%5.2f\n" % (np.average(mags),
                                               np.std(mags)))
    except IndexError:
        summary_message = "no primary target extracted"

    # add information to summary website, if requested
    if _pp_conf.use_diagnostics_summary:
        diag.insert_into_summary(summary_message)

    print('\nDone!\n')
    logging.info('----- successfully done with this process ----')

    gc.collect()  # collect garbage; just in case, you never know...

def read_yml_config(path_config: str) -> dict:
    """
    reads YAML config into a dictionary
    """
    with open('{}'.format(path_config)) as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
        print("\n----- Configfile parameters (pp_run.read_yml_config):")
        for key, value in config.items():
            print("{}:\n{}".format(key, str(value)))
    return config



if __name__ == '__main__':
    import os
    # command line arguments
    parser = argparse.ArgumentParser(description='automated WCS registration')
    parser.add_argument('-auto',
                        help='use automatic (default) pipeline procedure instead\
                         of the config file, works with console arguments,\
                         else work with the default config file',
                        action='store_true', default=False)
    parser.add_argument('-config',
                        help='use a custom configfile provided by the user\
                         instead of a default one', default=False)
    parser.add_argument('-prefix', help='data prefix',
                        default=None)
    parser.add_argument('-target', help='primary targetname override',
                        default=None)
    parser.add_argument('-filter', help='process images only with specified filter',
                        default=None)
    parser.add_argument('-fixed_aprad', help='fixed aperture radius (px)',
                        default=0)
    parser.add_argument('-source_tolerance',
                        help='tolerance on source properties for registration',
                        choices=['none', 'low', 'medium', 'high'],
                        default='high')
    parser.add_argument('-solar',
                        help='restrict to solar-color stars',
                        action="store_true", default=False)
    parser.add_argument('-rerun_registration',
                        help=('rerun registration step if not '
                              'successful for all images'),
                        action="store_true", default=False)
    parser.add_argument('-asteroids',                           # NOT WORKING CORRECTLY
                        help='extract all known asteroids',
                        action="store_true", default=False)
    parser.add_argument('-reject',
                        help='schemas for target rejection',
                        nargs=1, default='pos')
    parser.add_argument('-keep_wcs',
                        help='keep wcs information and skip registration',
                        action="store_true", default=False)
    parser.add_argument('images', help='images to process or \'all\'',
                        nargs='+')

    core_path = os.environ.get('PHOTPIPEDIR')
    args = parser.parse_args()
    auto = args.auto
    config = args.config
    # if auto - then reading console arguments
    if auto:
        prefix = args.prefix
        man_targetname = args.target
        man_filtername = args.filter
        fixed_aprad = float(args.fixed_aprad)
        source_tolerance = args.source_tolerance
        solar = args.solar
        rerun_registration = args.rerun_registration
        asteroids = args.asteroids
        rejectionfilter = args.reject
        keep_wcs = args.keep_wcs
        filenames = sorted(args.images)
    # if path to config is provided - use it instead of default config
    else:
        if config:
            configpath = args.config
        else:
            # path to the default configfile
            configpath = core_path + '/user_scripts/pp_config.yml'
            print('Using default config at {}'.format(configpath))
        # ======user modification ================
        # read configfile with pp_run arguments
        # path to photometrypipeline directory
        core_path = os.environ.get('PHOTPIPEDIR')
        # print('using configfile parameters, file {}'.format(configpath))
        config = read_yml_config(path_config=configpath)
        # =========== pp_run arguments ==========================
        prefix = config['pp_run'].get('prefix')
        man_targetname = str(config['pp_run'].get('target'))
        man_filtername = config['pp_run'].get('filter')
        fixed_aprad = config['pp_run'].get('fixed_aprad')
        solar = config['pp_run'].get('solar')
        rerun_registration = config['pp_run'].get('rerun_registration')
        asteroids = config['pp_run'].get('asteroids')
        rejectionfilter = config['pp_run'].get('reject')
        keep_wcs = config['pp_run'].get('keep_wcs')
        filenames = sorted(args.images)
        # ========== pp_prepare arguments =======================
        ra = config['pp_prepare'].get('ra')
        dec = config['pp_prepare'].get('dec')
        flipx = config['pp_prepare'].get('flipx')
        flipy = config['pp_prepare'].get('flipy')
        rotate = config['pp_prepare'].get('rotate')
        # ========= pp_register arguments =====================
        reg_max_rad = config['pp_register'].get('max_rad')
        reg_snr = config['pp_register'].get('snr')
        reg_minarea = config['pp_register'].get('minarea')
        reg_maxarea = config['pp_register'].get('maxarea')
        reg_aprad = config['pp_register'].get('aprad')
        source_tolerance = config['pp_register'].get('src_tol')
        reg_cat = config['pp_register'].get('reg_cat')
        # ========= pp_photometry argumetns ====================
        photo_snr = config['pp_photometry'].get('snr')
        photo_minarea = config['pp_photometry'].get('minarea')
        photo_maxarea = config['pp_photometry'].get('maxarea')
        photo_background = config['pp_photometry'].get('background_only')
        photo_target = config['pp_photometry'].get('target_only')
        photomode = config['pp_photometry'].get('photmode')
        # ========= pp_calibrate ===============================
        cal_minstars = config['pp_calibrate'].get('minstars')
        cal_maxstars = config['pp_calibrate'].get('maxstars')
        cal_catalog = config['pp_calibrate'].get('catalog')
        cal_maxflag = config['pp_calibrate'].get('maxflag')
    # if filenames = ['all'], walk through directories and run pipeline
    # each dataset
    _masterroot_directory = os.getcwd()

    if len(filenames) == 1 and filenames[0] == 'all':

        # dump data set information into summary file
        _pp_conf.use_diagnostics_summary = True
        diag.create_summary()

        # turn prefix and fits suffixes into regular expression
        if prefix is None:
            prefix = ''
        regex = re.compile('^'+prefix+'.*[fits|FITS|fit|FIT|Fits|fts|FTS]$')

        # walk through directories underneath
        for root, dirs, files in os.walk(_masterroot_directory):

            # ignore .diagnostics directories
            if '.diagnostics' in root:
                continue

            # identify data frames
            filenames = sorted([s for s in files if re.match(regex, s)])

            # call run_the_pipeline for each directory separately
            if len(filenames) > 0:
                print('\n RUN PIPELINE IN %s' % root)
                os.chdir(root)

                run_the_pipeline(filenames, man_targetname, man_filtername,
                                 fixed_aprad, source_tolerance, solar,
                                 rerun_registration, asteroids, keep_wcs)
                os.chdir(_masterroot_directory)
            else:
                print('\n NOTHING TO DO IN %s' % root)

    else:
        # call run_the_pipeline only on filenames
        run_the_pipeline(filenames, man_targetname, man_filtername,
                         fixed_aprad, source_tolerance, solar,
                         rerun_registration, asteroids, keep_wcs)
        pass

"""
Personal Photometry Pipeline Configuation File
2016-11-01, mommermiscience@gmail.com
"""

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

# telescope/instrument configurations

# MYTELESCOPE setup parameters
mytelescope_param = {
    'telescope_instrument': '70-cm reflector AZT-8/FLI ML47-10',  # telescope/instrument name
    'telescope_keyword': 'AZT-8',  # telescope/instrument keyword
    'observatory_code': '121',  # MPC observatory code
    'secpix': (0.94, 0.94),  # pixel size (arcsec) before binning

    # image orientation preferences
    'flipx': True,
    'flipy': True,
    'rotate': 0,

    # instrument-specific FITS header keywords
    'binning': ('XBINNING', 'YBINNING'),  # binning in x/y
    'extent': ('NAXIS1', 'NAXIS2'),  # N_pixels in x/y
    'ra': 'OBJCTRA',  # telescope pointing, RA
    'dec': 'OBJCTDEC',  # telescope pointin, Dec
    'radec_separator': ' ',  # RA/Dec hms separator, use 'XXX'
    # if already in degrees
    'date_keyword': 'DATE-OBS',  # obs date/time
    # keyword; use
    # 'date|time' if
    # separate
    'obsmidtime_jd': 'MIDTIMJD',  # obs midtime jd keyword
    # (usually provided by
    # pp_prepare
    'object': 'OBJECT',  # object name keyword
    'filter': 'FILTER',  # filter keyword
    'filter_translations': {'V': 'V', 'R': 'R',
                            'I': 'I', 'B': 'B',
                            'None': 'Clear'},
    # filtername translation dictionary
    'exptime': 'EXPTIME',  # exposure time keyword (s)
    #'airmass': 'AIRMASS',  # airmass keyword

    # source extractor settings
    'source_minarea': 12,  # default sextractor source minimum N_pixels
    'source_maxarea': 30,  # default sextractor source minimum N_pixels
    'source_snr': 7,  # default sextractor source snr for registration
    'aprad_default': 5,  # default aperture radius in px
    'aprad_range': [2, 10],  # [minimum, maximum] aperture radius (px)
    'sex-config-file': rootpath + '/setup/AZT-8.sex',
    'mask_file': {},
    #                        mask files as a function of x,y binning

    # scamp settings
    'scamp-config-file': rootpath + '/setup/AZT-8.scamp',
    'reg_max_mag': 18.5,
    'reg_search_radius': 1.0,  # deg
    'source_tolerance': 'high',

 # swarp settings
    #'copy_keywords': ('TELESCOP,INSTRUME,FILTER,EXPTIME,OBJECT,' +
     #                 'DATE-OBS,TIME-OBS,RA,DEC,SECPIX,AIRMASS,' +
      #                'TEL_KEYW,CCDBIN1,CCDBIN2,MIDTIMJD'),
    #                         keywords to be copied in image
    #                         combination using swarp
#    'swarp-config-file': rootpath+'/setup/vatt4k.swarp',

    # default catalog settings
    'astrometry_catalogs': ['GAIA'],
    'photometry_catalogs': ['PANSTARRS','SDSS-R9', 'APASS9', '2MASS']
}

# add telescope configurations to 'official' telescopes.py

implemented_telescopes.append('AZT-8')

# translate INSTRUME (or others, see _pp_conf.py) header keyword into
#   PP telescope keyword
# example: INSTRUME keyword in header is 'mytel'
instrument_identifiers['instrume_identifier'] = 'AZT-8'

# translate telescope keyword into parameter set defined here
telescope_parameters['AZT-8'] = AZT-8_param

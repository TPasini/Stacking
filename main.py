# -*- coding: utf-8 -*-
#
# Copyright (C) 2021 - Thomas Pasini
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# Credits for all libraries go to F. de Gasperin and H. W. Edler

import os, sys, shutil
from astropy.table import Table
from lib import lib_fits as libfits
import glob
import lib.lib_util as lib_util
import lib.lib_log as lib_log


def print_title():
    print("""

          _____   _                    _      _
         / ____| | |                  | |    (_)
        | (___   | |_    __ _    ___  | | __  _   _ __     __ _
         \___ \  | __|  / _` |  / __| | |/ / | | | '_ \   / _` |
         ____) | | |_  | (_| | | (__  |   <  | | | | | | | (_| |
        |_____/   \__|  \__,_|  \___| |_|\_\ |_| |_| |_|  \__, |
                                                           __/ |
                                                          |___/

      """)

    return

print_title()

logger_obj = lib_log.Logger('pipeline-stacking.logger')
logger = lib_log.logger
s = lib_util.Scheduler(log_dir=logger_obj.log_dir, dry = False)
w = lib_util.Walker('pipeline-stacking.walker')

#INPUT: DATASET OF GALAXY CLUSTERS/GROUPS
#import subtraction#: SUBTRACT COMPACT SOURCES AND CREATE NEW COLUMN
#import injection#: INJECT MOCK HALOS IN GIVEN POSITION

with w.if_todo('FITSgen'):
    logger.info('Generation of table fits file...')
    print('')
    import run_scripts.createfits #: CREATE A TABLE WITH THEIR PROPERTIES AND IMAGES NAMES

file_cat = 'LISTstacking.fits'

cat = Table.read(file_cat)
image = cat['Imagename']

with w.if_todo('Stacking'):
    logger.info('Smoothing images to common beam and stacking...')
    objectlist = libfits.AllImages(image)
    objectlist.convolve_to(circbeam=True)
    objectlist.write(suffix='smooth', inflate=True)
    import run_scripts.stacking #: STACK IMAGES OF MOCK HALOS AND PROVIDE PLOTS AND IMAGES OF RESULTS

if os.path.exists('Smoothed'):
    shutil.rmtree('Smoothed')

os.mkdir('Smoothed')
filelist = glob.glob('*-smooth.fits')
for el in filelist:
    shutil.move(el, 'Smoothed')

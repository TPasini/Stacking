import os, sys, shutil
from astropy.table import Table
from lib import lib_fits as libfits
import glob

#INPUT: DATASET OF GALAXY CLUSTERS/GROUPS
#import subtraction#: SUBTRACT COMPACT SOURCES AND CREATE NEW COLUMN
#import injection#: INJECT MOCK HALOS IN GIVEN POSITION

print('')
print('Creating fits table..')
import run_scripts.createfits #: CREATE A TABLE WITH THEIR PROPERTIES AND IMAGES NAMES

file_cat = 'LISTstacking.fits'

cat = Table.read(file_cat)
image = cat['Imagename']

print('Smoothing images to common beam..')
objectlist = libfits.AllImages(image)
objectlist.convolve_to(circbeam=True)
objectlist.write(suffix='smooth', inflate=True)

print('Stacking images..')
print('')


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

import run_scripts.stacking #: STACK IMAGES OF MOCK HALOS AND PROVIDE PLOTS AND IMAGES OF RESULTS

if os.path.exists('Smoothed'):
    shutil.rmtree('Smoothed')

os.mkdir('Smoothed')
filelist = glob.glob('*-smooth.fits')
for el in filelist:
    shutil.move(el, 'Smoothed')

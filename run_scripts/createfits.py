import numpy as np
import os, sys
import astropy.io.fits as fits

if os.path.exists('LISTstacking.fits'):
    os.remove('LISTstacking.fits')
    
#res = np.genfromtxt('LISTstacking.txt', delimiter=' ', dtype=None, names=('Name', 'z', 'M500', 'RA', 'DEC', 'Imagename', 'Inj_flux'))

res = np.loadtxt('LISTstacking.txt', delimiter=' ', usecols=[1,2,3,4])
names = np.loadtxt('LISTstacking.txt', delimiter=' ', usecols=[0,5,6], dtype=np.str)


a1 = np.array(names[:,0])
a2 = np.array(res[:,0])
a3 = np.array(res[:,1])
a4 = np.array(res[:,2])
a5 = np.array(res[:,3])
a6 = np.array(names[:,1])
a7 = np.array(names[:,2])
col1 = fits.Column(name='name', array=a1, format='A17')
col2 = fits.Column(name='z', array=a2, format='F', unit = '')
col3 = fits.Column(name='M500', array=a3, format='F', unit = '10^14 Msun')
col4 = fits.Column(name='RAinj', array=a4, format='F', unit = 'deg')
col5 = fits.Column(name='DECinj', array=a5, format='F', unit = 'deg')
col6 = fits.Column(name='Imagename', array=a6, format='A40')
col7 = fits.Column(name='Imagename-smooth', array=a7, format='A40')
cols = fits.ColDefs([col1, col2, col3, col4, col5, col6, col7])
hdu = fits.BinTableHDU.from_columns(cols)
hdu.writeto('LISTstacking.fits')

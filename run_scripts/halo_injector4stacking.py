#!/usr/bin/env python
# create radio halo lorenzian model, inject it into the uv data, make new image 

#usage: halo_injector.py center_RA(deg) center_DEC(deg) injected_flux(Jy) re(kpc) do_0inj ms_file(s)
#center_RA(deg) center_DEC(deg): position to inject mock halo
#injected_flux(Jy): flux density of mock halo
#re(kpc): e-folding radius of mock halo
#do_0inj(True/False): if True, image ONLY low resolution map without injection
#ms_file(s): name of ms file(s)


from astropy.io import fits
from astropy import units as u
import casacore.tables as tables
import glob
import os
import sys
import numpy as np
from pylab import meshgrid, imshow
from astropy.wcs import WCS
from astropy.wcs import WCS as pywcs
from astropy.cosmology import FlatLambdaCDM
from astroquery.simbad import Simbad
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
import csv



RA = float(sys.argv[1])
DEC = float(sys.argv[2])
flux = float(sys.argv[3])
#re = float(sys.argv[4])
do_0inj = str(sys.argv[4])
ms_files = sys.argv[5:]
#do_0inj = str(sys.argv[3])
#ms_files = sys.argv[4:]

"""
CLUSTER NAME & HEADER
"""

root_name=str(glob.glob('*.calibrated'))
root_name=str(root_name[10:27])

data_fits, header_fits = fits.getdata(root_name+'_image_9-MFS-image.fits', header=True)

imsize=int(header_fits['NAXIS1'])
ref_freq=float(header_fits['CRVAL3'])
pixscale = float(1.5)


"""
IMPORT REDSHIFT and MASS500
"""
def get_data(fname,colname):
    data=fits.open(fname)
    data=data[1].data
    return data[colname]

catalogue = '../clusterDR2.fits'
source_name = get_data(catalogue,'Name')
for i in source_name:
    source_name = source_name.replace(' ', '')
z = get_data(catalogue,'z')
M500 = get_data(catalogue,'MSZ')


for i, element in enumerate(source_name):
    if root_name==element:
        source_name=element
        z=z[i]
        M500=M500[i]*1.e14



re = 167.0

 

"""
ANGULAR/LINEAR CONVERSION FACTOR
"""

#arcsec/kpc factor
asecperkpc = cosmo.arcsec_per_kpc_proper(z).value
#kpc/arcsec factor
kpc_per_arcsec=1./asecperkpc
#luminosity distance
DL=cosmo.luminosity_distance(z).value


#e-folding radius in arcsec
re_as = re*asecperkpc

#arcsec corresponding to 50 kpc to be used as taper
taper50kpc=50/kpc_per_arcsec
taper50kpcSTRING=str(taper50kpc)+'asec'

"""
IMAGING PARAMETERS
"""

pixscale_tapered = int(taper50kpc/5.)
if int(pixscale_tapered) == 1:
    pixscale_tapered = 2
imsize_tapered = int(imsize*pixscale/pixscale_tapered)

#baseline averaging to be used in wsclean
baselineav = 150.e6*2.*np.pi*np.float(pixscale_tapered)/(24.*3600*np.float(imsize_tapered))


# limit baseline averaging to 10
if baselineav > 10.0:
    baselineav = 10.


niter=1000000

print('1 arcsec = ', round(kpc_per_arcsec,1), 'kpc' )
print('10 arcmin = ', round(kpc_per_arcsec*60*10,1), 'kpc' )



spidx = -1.3
DL = DL*(3.e22) #from Mpc to m


"""
FLATTEN: convert a fits file into a 2D image (with new header and data)
"""

def flatten(filename, channel=0, freqaxis=0):
    f = fits.open(filename)
    naxis=f[0].header['NAXIS']
    if naxis<2:
        raise RadioError('Cannot make map from this')
    if naxis==2:
        pass
        #return f[0].header,f[0].data

    w = pywcs(f[0].header)
    wn = pywcs(naxis=2)

    wn.wcs.crpix[0]=w.wcs.crpix[0]
    wn.wcs.crpix[1]=w.wcs.crpix[1]
    wn.wcs.cdelt=w.wcs.cdelt[0:2]
    wn.wcs.crval=w.wcs.crval[0:2]
    wn.wcs.ctype[0]=w.wcs.ctype[0]
    wn.wcs.ctype[1]=w.wcs.ctype[1]

    header = wn.to_header()
    header['NAXIS']=2
    header['NAXIS1']=f[0].header['NAXIS1']
    header['NAXIS2']=f[0].header['NAXIS2']
    copy=('EQUINOX','EPOCH')
    for k in copy:
        r=f[0].header.get(k)
        if r:
            header[k]=r

    dataslice=[]
    for i in range(naxis,0,-1):
        if i<=2:
            dataslice.append(np.s_[:],)
        elif i==freqaxis:
            dataslice.append(channel)
        else:
            dataslice.append(0)
# add freq
    header["FREQ"] = f[0].header['CRVAL3']

    # add beam if present
    try:
        header["BMAJ"]=f[0].header['BMAJ']
        header["BMIN"]=f[0].header['BMIN']
        header["BPA"]=f[0].header['BPA']
    except:
        pass
    return header, f[0].data[tuple(dataslice)]



"""
SUMMARY OF INPUT PARAMETERS
"""

root_img=root_name




print('object =', root_name)
print('redshift =' ,z)
#print('imsize =' ,imsize)
#print('pixel size =' ,pixscale, 'arcsec')
print('center of mock halo =' ,RA, DEC, 'deg')
print('flux density of mock halo =' ,flux*1e3, 'mJy')
#print('reference frequency =' ,round(ref_freq/1e6,0), 'MHz')
print('1 arcsec = ', round(kpc_per_arcsec,1), 'kpc' )
print('50 kpc = ', round(taper50kpc,1), 'arcsec' )
print('luminosity distance =', round(DL/(3.e22),1), 'Mpc')
print('re =' ,re, 'kpc')
#print('2.6Re = ' , round(2.6*re_as,0), ' arcsec')
print('3Re = ' , round(3.*re_as,0), ' arcsec')
print('D = 6Re = ' , round(6.*re_as,0), ' arcsec')
#print('5Re = ' , round(5.*re_as,0), ' arcsec')
#print('D = 10Re = ' , round(10.*re_as,0), ' arcsec')
#print('spectral index of mock halo =', spidx)
print('M500 = ', round(M500,2), 'e14.9 Msun')
#print('expected Pradio_Lband = ', round(Pradio_Lband,1), ' W/Hz')
print('exptected R_e = ', round(re,1), ' kpc')
#print('flux_Lband = ', round(flux_Lband*1000,1), ' mJy')
#print('flux_144 = ', round(flux*1000,1), ' mJy')

#sys.exit()

ms_files_wsclean = ' '.join(ms_files)



"""
LOW RESOLUTION IMAGE WITHOUT INJECTION
"""


if do_0inj=='True':
    os.system('wsclean -no-update-model-required -minuv-l 80.0 -size '+str(int(imsize_tapered))+' '+str(int(imsize_tapered))+' -reorder -weight briggs -0.5 -multiscale -taper-gaussian '+str(taper50kpcSTRING)+' -weighting-rank-filter 3 -clean-border 1 -mgain 0.8 -fit-beam -data-column DIFFUSE_SUB -join-channels -channels-out 6 -padding 1.4 -auto-mask 2.5 -auto-threshold 1.0 -fit-spectral-pol 3 -pol i -baseline-averaging '+str(baselineav)+' -name '+root_img+'_ORIGINAL_T50 -scale '+str(pixscale_tapered)+'arcsec -niter '+str(niter)+' '+ms_files_wsclean+' >wsclean0.log')
    
    #remove useless files (need only MFS image and residual)
    
    os.system('rm -rf *ORIGINAL_T50*000*.fits') #channel images
    os.system('rm -rf *ORIGINAL_T50*-dirty.fits') #dirty images
    os.system('rm -rf *ORIGINAL_T50*-psf.fits') #psf images
    sys.exit()


"""
CREATE MODELS: run wsclean with niter 1
"""


if not os.path.exists(root_img+'_mod'+str(flux*1000)+'mJy-MFS-model.fits'):
    os.system('wsclean -no-update-model-required -size '+str(imsize)+' '+str(imsize)+' -reorder -weight briggs 0.0 -mgain 0. -data-column DATA -join-channels -channels-out 6 -pol i -baseline-averaging '+str(baselineav)+' -name '+root_img+'_mod'+str(flux*1000)+'mJy -scale '+str(pixscale)+'arcsec -niter 1 '+ms_files_wsclean+' >wsclean1.log')
    
    #remove useless files (need only models of channel maps, e.g. A2472_mod0.01Jy-0000-model.fits)


models = glob.glob(root_img+'_mod'+str(flux*1000)+'mJy-*-model.fits')




"""
UPDATE MODELS: add the mock halo
"""

head, data = flatten(models[0])
wcs=WCS(head)
center_ra_px, center_dec_px = wcs.wcs_world2pix(RA*u.deg, DEC*u.deg, 1, ra_dec_order=True)

x = np.linspace(0, imsize, imsize)
X,Y = meshgrid(x, x) 

#model with a Lorentz function
def lorentz_2D(flux, spidx, freq, re_as, center_ra, center_dec, x, y):
    re_px = re_as/pixscale
    flux_freq = flux*(freq/ref_freq)**spidx
    I_0 = flux_freq/(2.*np.pi*re_as**2) #Jy/arcsec^2
    r=np.sqrt((x-center_ra)**2+(y-center_dec)**2)
    lorentz = I_0 * np.exp(- r/re_px) * pixscale**2 #Jy/pixel

    
    if abs(ref_freq-freq) < 3.0:
        print('I0 = ' ,  round(I_0*1000000,1), ' muJy/arcsec2 at ', round(freq/(1.e6),0), ' MHz')
        print('e-folding radius = ' , re_as, ' arcsec')

    return lorentz



for model in models:
    hdr = fits.getheader(model)
    freq = hdr['CRVAL3']
    mock_halo = lorentz_2D(flux, spidx, freq, re_as, center_ra_px, center_dec_px, X, Y)
    model_update = 'inject_'+model
    os.system('cp '+model+' '+model_update)

    fits.update(model_update, mock_halo, hdr)
    
        

# FT the model to obtain mock visibilities
os.system('wsclean -predict -name inject_'+root_img+'_mod'+str(flux*1000)+'mJy -channels-out 6 '+ms_files_wsclean+' >wsclean_predict.log')



# add model visibilities to the data
for ms_file in ms_files:
    tables.taql("update $ms_file set MODEL_DATA = DIFFUSE_SUB + MODEL_DATA")
#    tables.taql("update $ms_file set MODEL_DATA =  MODEL_DATA")




"""
IMAGE OF THE MOCK HALO
"""


os.system('wsclean -no-update-model-required -minuv-l 80.0 -size '+str(int(imsize_tapered))+' '+str(int(imsize_tapered))+' -reorder -weight briggs -0.5 -multiscale -taper-gaussian '+str(taper50kpcSTRING)+' -weighting-rank-filter 3 -clean-border 1 -mgain 0.8 -fit-beam -data-column MODEL_DATA -join-channels -channels-out 6 -padding 1.4 -auto-mask 2.5 -auto-threshold 1.0 -fit-spectral-pol 3 -pol i -baseline-averaging '+str(baselineav)+' -name '+root_img+'_MOCK_'+str(flux*1000)+'mJy_T50 -scale '+str(pixscale_tapered)+'arcsec -niter '+str(niter)+' '+ms_files_wsclean+' >wsclean2.log')




# remove useless files

os.system('rm -rf '+str(root_img)+'_mod'+str(flux*1000)+'*000*') #niter 1 channel images
os.system('rm -rf '+str(root_img)+'_mod'+str(flux*1000)+'*MFS*') #niter 1 MFS images
os.system('rm -rf inject_'+str(root_img)+'_mod'+str(flux*1000)+'*000*') #injected model images
os.system('rm -rf '+str(root_img)+'_MOCK_'+str(flux*1000)+'*T50*000*') #mock image channel maps
os.system('rm -rf '+str(root_img)+'_MOCK_'+str(flux*1000)+'*T50*psf.fits') #mock image psf map
os.system('rm -rf '+str(root_img)+'_MOCK_'+str(flux*1000)+'*T50*dirty.fits') #mock image dirty map


os.system('cp '+str(root_img)+'_MOCK_'+str(flux*1000)+'mJy_T50-MFS-image.fits '+str(root_img)+'_MOCK.fits ')



#if os.path.exists('../LISTstacking.txt'):
#    os.remove('../LISTstacking.txt')
with open('../LISTstacking.txt', 'a') as f:
    writer = csv.writer(f, delimiter=" ", quoting=csv.QUOTE_NONE, escapechar=' ')
#    writer.writerow(['#name', 'z', 'M500', 'RAinj', 'DECinj', 'imagename', 'imagenamesmooth', 'Sinj'])
    writer.writerow([str(root_img), z, M500, RA, DEC, str(root_img)+'_MOCK.fits', str(root_img)+'_MOCK-smooth.fits', flux*1000])

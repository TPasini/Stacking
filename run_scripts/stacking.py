import os, sys, glob, itertools, shutil
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from scipy.ndimage.filters import minimum_filter, maximum_filter
import lib.lib_coordinates_mode as cm
from lib.lib_fits import flatten
import lib.lib_fits
from astropy import wcs as pywcs
from astropy.cosmology import FlatLambdaCDM
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.nddata import Cutout2D
from astropy import constants as const
import astropy.io.fits as aif
import lib.lib_radio as radiolib
import lib.scipy_modified as scimod
import astropy.wcs.utils as astrowcs
from pathlib import Path

cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

do_plot = True # True to save a plot for each source

if os.path.exists('/Stacking_plots'):
    os.remove('/Stacking_plots')
    
if os.path.exists('Cut_images'):
    shutil.rmtree('Cut_images')
os.mkdir('Cut_images')

if os.path.exists('Regridded'):
    shutil.rmtree('Regridded')
os.mkdir('Regridded')

file_cat = 'LISTstacking.fits'

cat = Table.read(file_cat)
image = cat['Imagename-smooth']
cat['flux1mpc'] = np.zeros(len(cat))
cat['std'] = np.zeros(len(cat))
cat['lum1mpc'] = np.zeros(len(cat))
cat['lum_dist'] = np.zeros(len(cat))
print('The catalog has', len(cat), 'initial entries.')

# get min max z
min_z = min(cat['z'])
max_z = max(cat['z'])
mean_z = np.mean(cat['z'])
print('The minimum redshift is', min_z)
print('The maximum redshift is', max_z)
print('The mean redshift is', "%.2f" % mean_z)

class Source(object):

    def __init__(self, c):
        self.c = c
        self.file_surv = None
        self.data_regrid = None
        self.data_open = None
    
    def reproject(self, file_surv, size, beamspacing):
        """
        create data_regrid from data_surv
        """
        header_surv, data_surv = flatten(file_surv)
        w = pywcs.WCS(header_surv)
        coord_source = SkyCoord((self.c['RAinj']*u.deg).astype(np.float64), (self.c['DECinj']*u.deg).astype(np.float64), frame='fk5')
        
        data_cut = Cutout2D(data_surv, coord_source, size=[size,size], wcs=w, mode='partial', fill_value=np.nan)

        header_surv['CRPIX1'] = size/2
        header_surv['CRPIX2'] = size/2
        
        aif.writeto('Cut_images/'+imagetouse+'-cut.fits', data_cut.data, header = header_surv)

        del data_surv
        print("Regrid", data_cut.shape[0], "-> ", size_in_pixel, 'pixel')

        prova = scimod.imresize(data_cut.data, (size_in_pixel,size_in_pixel), mode='F', interp='bilinear')
        header_surv['CRPIX1'] = size_in_pixel / 2
        header_surv['CRPIX2'] = size_in_pixel / 2
        header_surv['CDELT2'] = beamspacing.value * (size / size_in_pixel)
        header_surv['CDELT1'] = - beamspacing.value * (size / size_in_pixel)
        aif.writeto('Regridded/' + imagetouse + '-regrid.fits', prova, header=header_surv)

        return scimod.imresize(data_cut.data, (size_in_pixel,size_in_pixel), mode='F', interp='bilinear')

    def calc_flux(self, data):
        """
        Get flux
        """
 
        pixel_per_arcmin = ((80.*u.kpc)/cosmo.kpc_proper_per_arcmin(self.c['z']))/(beamspacing*u.arcmin/u.pixel) #How many pixel per arcmin
        beam_area_arcmin = (radiores*u.arcmin/60.)*(radiores*u.arcmin/60.)*np.pi/(4*np.log(2.)) # beam area in arcmin**2
        beam_area_pixel = beam_area_arcmin*pixel_per_arcmin**2 # beam area in pixel**2 for 50 kpc radius
        flux = np.nansum(data[mask_flux1mpc])*u.Jy/beam_area_pixel

        print('Flux density inside the mask:', "%.2f" % flux.value, 'Jy')
        return flux

    def calc_lumdist(self, data):
        """
        Get luminosity distance
        """
        lumdist = cosmo.luminosity_distance(self.c['z'])
        return lumdist


    def calc_noise(self, data, niter=100, eps=1e-6):
        """
        Return the rms of the mask_rms region
        """
        
        rms = 1.; oldrms = 1.
        for i in range(niter):
            rms = np.nanstd(data[mask_rms & (np.abs(data)<3*rms)])
            if np.abs(oldrms-rms)/rms < eps:
                print("Std:", "%.2e" % (rms*1e+6), 'uJy/beam')
                return rms

            oldrms = rms

        raise Exception('Noise estimation failed to converge.')
        
#    def make_mask(image_name, mask_name=None, threshpix=5, atrous_do=False, rmsbox=(100,10), adaptive_thresh=50, write_srl=False, write_gaul=False, write_ds9=False, mask_combine=None):
#
#        # wavelets are required to fit gaussians
#        if atrous_do or write_srl or write_ds9 or write_gaul: stop_at = None
#        else: stop_at = 'isl'
#
#        # DO THE SOURCE DETECTION
#        img = bdsf.process_image(image_name, rms_box=rmsbox, frequency=54e6,
#            thresh_isl=float(threshpix*4/5), thresh_pix=float(threshpix), rms_map=True, mean_map='zero', atrous_do=atrous_do, atrous_jmax=4,
#            adaptive_rms_box=True, adaptive_thresh=adaptive_thresh, rms_box_bright=(30,5),
#            flagging_opts=True, flag_maxsize_fwhm=0.5, stop_at=stop_at, quiet=True, debug=False)
#
#        # WRITE THE MASK FITS
#        if mask_name == None: mask_name = image_name+'.newmask'
#        if os.path.exists(mask_name): os.system('rm -r ' + mask_name)
#        img.export_image(img_type='island_mask', img_format='fits', outfile=mask_name, clobber=True)
        
class Source_group(object):

    def __init__(self, name=None):
        self.sources = []
        self.sources_selection = []
        self.name = name
        self.zrange = [0,999]

    def add_source(self, source):
        self.sources.append(source)

    def select_sources(self, zrange=None, valname='', valstdname=''):
        if zrange is None: zrange=[0,999]
        #print("Selection: - zrange=",zrange,"on",valname)
        self.zrange = zrange
        self.sources_selection = []

        sources = []; fluxes = []
        # first cut in z and get fluxes
        for source in self.sources:
            if source.c['z'] <= zrange[0]: continue
            if source.c['z'] > zrange[1]: continue
            sources.append(source)
            fluxes.append(source.c[valname])
        
        #ONLY FOR NON-RESIDUAL IMAGES
        for source in sources:
#            if source.c[valname] < 5e-3 or source.c[valname] > 10:
#            if source.c[valname] > 1.:
#                print("Remove source: too high contribution - flux is %.2f" % (source.c[valname]))
#            elif source.c[valstdname] > 1.5*0.45e-3:
#                print("Remove source: too high std - std is %.3f mJy/b" % (source.c[valstdname]*1.e3))
#            else:
            self.sources_selection.append(source)

        #print("Selecting: %i/%i" % (len(self.sources_selection), len(self.sources)))

    def get_val(self, order=None, valname=''):
        """
        Return fluxes/lums ordered by sources' redshift
        """

        vals = []; redshifts = []
        for source in self.sources_selection:
            vals.append(source.c[valname])
            redshifts.append(source.c['z'])
        if order == 'redshift': vals = [val for (redshift, val) in sorted(zip(redshifts,vals))]
        return vals
        
        
    def stack(self, mode=''):
        data_stack_all = np.zeros((size_in_pixel,size_in_pixel))
        data_stack_all_rms = np.zeros((size_in_pixel,size_in_pixel))
        for source in self.sources_selection:
            if mode == "open":
                source.data_open[ np.isnan(source.data_open) ] = 0.
                data_stack_all += source.data_open*((source.c['std']**2)+ (source.c['lum_dist']**2))
                data_stack_all_rms += (source.c['std']**2)+(source.c['lum_dist']**2)
            elif mode == "orig":
                source.data_regrid[ np.isnan(source.data_regrid) ] = 0.
                data_stack_all += source.data_regrid*((1/(source.c['std']**2))*(1/((1+(source.c['z']))**4)))
                data_stack_all_rms += (1/(source.c['std']**2))*(1/((1+(source.c['z']))**4))
            else:
                print("ERROR: wrong mode")
                sys.exit(1)

        return data_stack_all/data_stack_all_rms
        

    def plot(self, filename):
        fig = plt.figure(figsize=(18, 10))
        fig2 = plt.figure(figsize=(15, 15))
        fig3 = plt.figure(figsize=(15, 15))
        fig4 = plt.figure(figsize=(13, 13))
        fig5 = plt.figure(figsize=(15, 15))

        self.select_sources(self.zrange, 'flux1mpc', 'std')
        img_stack_orig = self.stack(mode='orig')

        def blockPrint():
            sys.stdout = open(os.devnull, 'w')

        def enablePrint():
            sys.stdout = sys.__stdout__

        blockPrint()
        stackednoise = Source.calc_noise(self, img_stack_orig)

        enablePrint()
        
        print('')
        print('#####################')
        print('### STACKED IMAGE ###')
        print('#####################')
        print('')
        
        print('The rms noise of the stacked image is', "%2.3e" % (stackednoise * 1e+6), 'uJy/beam')
        #print('This is', "%.2f" % float((np.mean(cat['std'])) / stackednoise), 'times the starting noise')
        #print('Theoretical improvement:', "%.2f" % float(np.sqrt(len(cat))))

        detection_threshold = 2 * stackednoise

        i = 0
        j = 0
        for n, cell in enumerate(img_stack_orig[int(size_in_pixel/2)]):
            j = j + 1
            if cell > detection_threshold:
                i = i + 1

        detection = 0
        if i > j/10:
            detection = 1


        if detection == 1:
            print('')
            print('#########################')
            print('DIFFUSE EMISSION DETECTED')
            print('#########################')
        else:
            print('')
            print('#############################')
            print('DIFFUSE EMISSION NOT DETECTED')
            print('#############################')

        ##STACKED IMAGE
        ax1 = fig.add_subplot(111)
        # ax1.set_title(r'Stack of %i images' % len(self.sources_selection))
        #ax1.set_xlabel(r'Ra [Mpc]', fontsize=18)
        #ax1.set_ylabel(r'Dec [Mpc]', fontsize=18)
        ax1.set_xlabel(r'$\Delta$ Ra [pixel]', fontsize=18)
        ax1.set_ylabel(r'$\Delta$ Dec [pixel]', fontsize=18)
        # circle1 = plt.Circle((0, 0), 0.5, color='r', fill=False)
        # circle2 = plt.Circle((0, 0), 1.0, color='b', fill=False)
        # ax1.add_artist(circle1)
        # ax1.add_artist(circle2)
        ax1.tick_params(labelsize=17)
        norm = colors.SymLogNorm(linthresh=0.00003, vmin=-0.00016, vmax=0.0025)
        im = ax1.imshow(img_stack_orig, extent=(-size_in_pixel, size_in_pixel, -size_in_pixel, size_in_pixel), norm=norm, origin='lower', cmap='viridis')
        cbar = fig.colorbar(im)
        cbar.set_label('Jy/beam', size=16)
        cbar.ax.tick_params(labelsize=14)

        if os.path.exists('Stacking_plots/stack.fits'):
            os.remove('Stacking_plots/stack.fits')

        header_orig = hdul[0].header
        
        header_orig['CRPIX1'] = size_in_pixel/2
        header_orig['CRPIX2'] = size_in_pixel/2
        header_orig['CDELT2'] = beamspacing.value/60 * (size / size_in_pixel)
        header_orig['CDELT1'] = - beamspacing.value/60 * (size / size_in_pixel)

        aif.writeto('Stacking_plots/stack.fits', img_stack_orig, header=header_orig)

        ##FLUXES
        ax2 = fig2.add_subplot(111)
        ax2.set_xlabel(r'source #', fontsize=20)
        ax2.set_ylabel(r'Flux density [Jy]', fontsize=20)
        fluxes = self.get_val(order='redshift', valname='flux1mpc')
        ax2.plot(fluxes, 'bo', label='Ordered in redshift')
        ax2.plot(fluxes, 'b-', linewidth=2.5)
        ax2.plot(sorted(fluxes), 'ro', label='Ordered in flux')
        ax2.plot(sorted(fluxes), 'r-', linewidth=2.5)
        ax2.legend(loc=2, prop={'size': 20})
        ax2.tick_params(labelsize=19)
        # ax2.set_title('Fluxes - avg: %.1f mJy' % (np.mean(fluxes)*1e3))

        ###RADIAL PROFILE
        ax3 = fig3.add_subplot(111)
        # ax3.set_title('Radial profile')
        ax3.set_xlabel(r'Dist (pixel)', fontsize=20)
        ax3.set_ylabel(r'Flux density', fontsize=20)
        # ax3.axvline(25, color='r')
        # ax3.axvline(50, color='b')
        ax3.set_xlim(0, size_in_pixel / 2)
        # ax.plot(self.radial_profile(img_stack_open), 'k-', label='Open image')
        ax3.plot(self.radial_profile(img_stack_orig), 'k:')
        ax3.tick_params(labelsize=19)
        # ax3.legend(loc=1)
        ax3.set_yscale('log')

        ##REDSHIFT DISTRIBUTION
        ax4 = fig4.add_subplot(111)
        # ax4.set_title('Redshift distribution')
        ax4.set_xlabel(r'Redshift', fontsize=24)
        ax4.set_ylabel(r'Number of groups', fontsize=24)
        ax4.hist([source.c['z'] for source in self.sources_selection], color='red', stacked=True, edgecolor='k', alpha=1,fill=True, linewidth=2, bins=8)
        ax4.tick_params(labelsize=23)

        # CUMULATIVE LUMINOSITY DISTRIBUTION
#        ax5 = fig5.add_subplot(111)
#        ax5.set_xlabel(r'source #', fontsize=20)
#        ax5.set_ylabel(r'Cumulative Radio Luminosity [10$^{21}$ W/Hz]', fontsize=20)
#        lumin = self.get_val(valname='lum1mpc')
#        weight = self.get_val(valname='std')
#        ax5.plot(sorted(np.cumsum(lumin)*weight / 1e+21), 'ko')
#        ax5.plot(sorted(np.cumsum(lumin)*weight / 1e+21), 'k-', linewidth=2.5)
#        ax5.tick_params(labelsize=19)
        # ax5.set_title('Luminosity cumulative distribution')

        fig.savefig('Stacking_plots/Stack.pdf', bbox_inches='tight')
        fig2.savefig('Stacking_plots/Fluxes.png', bbox_inches='tight')
        fig3.savefig('Stacking_plots/Radial profile.png', bbox_inches='tight')
        fig4.savefig('Stacking_plots/z_distr.pdf', bbox_inches='tight')
        #fig5.savefig('Stacking_plots/cumdistr.png', bbox_inches='tight')
        plt.close()



    def radial_profile(self, data):
        """
        Return radial profile of a 2d array,
        center is a len=2 list
        """
        
        center=[size_in_pixel/2,size_in_pixel/2]
        
        y, x = np.indices((data.shape))
        r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
        r = r.astype(int)
        
        tbin = np.bincount(r.ravel(), data.ravel())
        nr = np.bincount(r.ravel())
        radialprofile = tbin / nr
        return radialprofile[0:int(size_in_pixel/2)] # cut at size/2 to avoid corners

##################
#####MAIN CODE####
##################

sources = Source_group()

for i, c in enumerate(cat):

    imagetouse = c['Imagename-smooth']
    z = c['z']
    ra = c['RAinj']
    dec = c['DECinj']
    
    source = Source(c)
    
    print('')
    print('###################')
    print('##### IMAGE', i+1, '#####')
    print('###################')
    print('')

    source.file_surv = os.path.basename(imagetouse)
    hdul = aif.open(imagetouse)
    print('Image:', imagetouse)

    beamspacing = ((hdul[0].header['CDELT2'])*60)*u.arcmin
    radiores = ((hdul[0].header['BMAJ'])*3600)*u.arcsec
    resmin = ((hdul[0].header['BMIN'])*3600)*u.arcsec
    print('The beam is', "%.2f" %  radiores.value, 'x', "%.2f" % resmin.value, 'arcsec')
    #print('The beam spacing is', "%.2f" % (beamspacing.value*60), 'arcsec, leading to: 1 pixel =', "%.2f" % (beamspacing.value*60), 'arcsec')
    print('The pixel spacing is: 1 pixel =', "%.2f" % (beamspacing.value*60), 'arcsec')

    sizemean = 2000*u.kpc/cosmo.kpc_proper_per_arcmin(mean_z)
    size_in_pixel_withunit = (sizemean/beamspacing)*u.pixel #pixel for 2 Mpc
    size_in_pixel = int(size_in_pixel_withunit.value)

    if not os.path.exists('Stacking_plots'):
        os.makedirs('Stacking_plots')

#    #Flux mask
    mask_flux1mpc = np.zeros((size_in_pixel,size_in_pixel), dtype=bool)
    y,x = np.ogrid[-size_in_pixel/2:size_in_pixel/2, -size_in_pixel/2:size_in_pixel/2]
    #mask_flux1mpc[ x**2+y**2 <= ((size_in_pixel/2)/10)**2 ] = True
    mask_flux1mpc[x ** 2 + y ** 2 <= ((size_in_pixel / 4) / 1) ** 2] = True
    
    # rms mask
    mask_rms = np.zeros((size_in_pixel,size_in_pixel), dtype=bool)
    y,x = np.ogrid[-size_in_pixel/2:size_in_pixel/2, -size_in_pixel/2:size_in_pixel/2]
    mask_rms[ x**2+y**2 <= (size_in_pixel/1.5)**2 ] = True
    mask_rms[ x**2+y**2 <= ((size_in_pixel/4)/1)**2 ] = False
    
    fig = plt.figure(figsize=(18,6))
    
    size = 800

    source.data_regrid = source.reproject(source.file_surv, size, beamspacing/60)
    #header_orig = hdul[0].header
    
    if source.data_regrid is None:
        continue

    source.c['lum_dist'] = source.calc_lumdist(source.data_regrid).value
    source.c['flux1mpc'] =  source.calc_flux(source.data_regrid).value
    
    source.c['std'] = source.calc_noise(source.data_regrid)

    lum1mpc = source.c['flux1mpc']*u.Jy * 4*np.pi * (cosmo.luminosity_distance(z))**2
    #print('Lum within flux mask:', lum1mpc.to('W/Hz'))
    #source.c['lum1mpc'] = lum1mpc.to('W/Hz').value

    sources.add_source(source)
    
    if do_plot:
        z = source.c['z']
        fig = plt.figure(figsize=(12,6))
        ax = fig.add_subplot(121)
        ax.set_title('Original - rms noise='+str("%.3f" % (source.c['std']*1e3))+' mJy/b')
        ax.set_xlabel(r'Ra [Mpc]')
        ax.set_ylabel(r'Dec [Mpc]')
        #circle1 = plt.Circle((0, 0), 0.2, color='r', fill=False)
        #ax.add_artist(circle1)
        #norm = colors.SymLogNorm(linthresh=450e-6, vmin=-450e-6, vmax=0.005)
        norm = colors.SymLogNorm(linthresh=0.00016, vmin=-0.00016, vmax=0.0025)
        ax.imshow(source.data_regrid, extent=(-1,1,-1,1), norm=norm, origin='lower')
        print('Saving: Stacking_plots/'+str(i)+'.png')
        fig.savefig('Stacking_plots/'+str(i)+'.png', bbox_inches='tight')
        plt.close()

sources.plot(filename='Stacking_plots/stack.png')

print('')
print('Images and plots are in the running directory.')
print('')




#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import os, sys
from astroquery.simbad import Simbad
import numpy as np
import glob
from astropy.io import fits
import pyrap.tables as pt
import os.path
import bdsf
import pyregion
import argparse
import pickle
import aplpy
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def cleanup():
    os.system('rm -rf *000*.fits') # remove channel maps
    os.system('rm -rf *_subROBUST*.fits') # remove non-masked images
    os.system('rm -rf *_ROBUST*fits') # remove non-masked images
    os.system('rm -rf *-dirty.fits') # remove all dirty images
    return


def getimsize(boxfile, cellsize=1.5):
   """
   find imsize need to image a DS9 boxfile region
   """
   r = pyregion.open(boxfile)
   
   xs = np.ceil((r[0].coord_list[2])*1.6*3600./cellsize)
   ys = np.ceil((r[0].coord_list[3])*1.6*3600./cellsize)

   imsize = np.ceil(xs) # // Round up decimals to an integer
   if(imsize % 2 == 1): 
       imsize = imsize + 1
   return np.int(imsize)

def compute_uvmin(redshift,sourceLLS=1.0):
    '''
    sourceLLS in units of Mpc# 
    taper output for WSClean in arcsec
    '''
    oneradinmpc = cosmo.angular_diameter_distance(redshift)/(360./(2.*np.pi))
    scalebarlengthdeg    = sourceLLS/oneradinmpc.value
    
    return 1./(scalebarlengthdeg*np.pi/180.)

def compute_taper(redshift,taperscale):
    '''
    taperscale in units of kpc# 
    '''
    oneradinmpc = cosmo.angular_diameter_distance(redshift)/(360./(2.*np.pi))
    taper    = 1e-3*taperscale/(oneradinmpc.value)
    
    return taper*3600

def adjustniter_for_taper(taper, niter):
  if taper < 5:
    return np.int(niter)
  if taper >= 5 and taper < 15:
    return np.int(niter/2)  
  if taper >= 15:
    return np.int(niter/4)


def makeimage(mslist, imageout, pixsize, imsize, channelsout=6, niter=15000, robust=-0.5, minuv=80, uvtaper=None, multiscale=False, predict=True,fitsmask=None, deepmultiscale=False, cluster_redshift=None, column=None):

    # some setup
#    username = os.getlogin()
#    if username == 'rvweerenold':
#       wsclean = '/net/lofar1/data1/rvweeren/software/wsclean-code-2.6oct12/wsclean/build/wsclean'
#    else:
    wsclean = 'wsclean'


    if uvtaper != None:
       if uvtaper < 0:
           print('Not supported uvtaper', uvtaper)
       else:
           print(imsize, pixsize, uvtaper)
           imsizein  = np.int(np.float(imsize)*(pixsize/(uvtaper/5.)))
           pixsizein = np.int(uvtaper/5.)
           if np.float(pixsizein) < pixsize: # to deal with rounding issues which cause a 0arcsec pixelsize
              pixsizein = pixsize
              imsizein  = imsize
              
    else:
       imsizein  = imsize 
       pixsizein = pixsize

    if int(imsizein) < 511: # otherwise images get too small for multiscales
        imsizein = 512
    
    baselineav = 2.5e3*60000.*2.*np.pi *np.float(pixsizein)/(24.*60.*60*np.float(imsizein)) 

    # limit baseline averaging to 10, fixes prob
    if baselineav > 10.0:
        baselineav = 10.
    
    baselineav = str (baselineav)
 
    # few simple checks to make sure we have useful data
    msliststring = ' '.join(map(str, mslist))
    os.system('rm -f ' + imageout + '-*.fits')
    imcol = 'CORRECTED_DATA'
    t = pt.table(mslist[0],readonly=True) # just test for first ms in mslist
    colnames =t.colnames()
    if 'CORRECTED_DATA' not in colnames: # check which column to image
      imcol = 'DATA' 
    t.close()
    
    if column != None:
      imcol = column


    
    # build wsclean command      
    cmd = wsclean + ' '
    cmd += '-no-update-model-required -minuv-l ' + str(minuv) + ' '
    cmd += '-size ' + str(imsizein) + ' ' + str(imsizein) + ' -reorder '
    cmd += '-weight briggs ' + str(robust) + ' -weighting-rank-filter 3 -clean-border 1 '
    cmd += '-mgain 0.8 -fit-beam -data-column ' + imcol +' -join-channels -channels-out '
    cmd += str(channelsout) + ' -padding 1.4 '
    #cmd += '-parallel-deconvolution ' + str(np.int(imsizein)/2) + ' ' 
    if multiscale:
       if predict:
         cmd += '-multiscale '+' -multiscale-scales 0,2,4,8,16 '  
       else:
         cmd += '-multiscale '+' -multiscale-scales 0,4,8,16,32,64 '

    if fitsmask != None:
      if os.path.isfile(fitsmask): 
        cmd += '-fits-mask '+ fitsmask + ' '
      else:
        print('fitsmask: ', fitsmask, 'does not exist')
        sys.exit()
    else:
        cmd += '-auto-mask 2.5 -auto-threshold 1.0 '
    
    if uvtaper != None:
       cmd += '-taper-gaussian ' +  str(uvtaper) + 'arcsec '
    
    cmd += '-fit-spectral-pol 3 '# -beam-shape 6arcsec 6arcsec 0deg '   
    cmd += '-pol i '
    cmd += '-baseline-averaging ' + baselineav + ' '
      
    
    cmd += '-name ' + imageout + ' -scale ' + str(pixsizein) + 'arcsec ' 

    print('WSCLEAN: ', cmd + '-niter ' + str(niter) + ' ' + msliststring)
    #logging.info(cmd + '-niter ' + str(niter) + ' ' + msliststring)
    os.system(cmd + '-niter ' + str(niter) + ' ' + msliststring)

    if deepmultiscale:
        
      # predict first to fill MODEL_DATA so we can continue with clean
      cmdp = wsclean + ' -size ' 
      cmdp += str(imsizein) + ' ' + str(imsizein) + ' -channels-out ' + str(channelsout) + ' -padding 1.4 -predict ' 

      
      cmdp += '-name ' + imageout + ' -scale ' + str(pixsizein) + 'arcsec ' + msliststring
      print('PREDICT STEP for continue: ', cmdp)
      os.system(cmdp)
       
      # NOW continue cleaning  
      cmd += '-niter ' + str(niter/15) + ' -multiscale -continue ' + msliststring
      print('WSCLEAN continue: ', cmd)
      os.system(cmd)

    # REMOVE nagetive model components, these are artifacts (only for Stokes I)
    #if idg:
    #  removenegativefrommodel(sorted(glob.glob(imageout +'-????-I-model*.fits')))  # only Stokes I
    #else:    
    #  removenegativefrommodel(sorted(glob.glob(imageout + '-????-model.fits')))

    if predict:
      cmd = wsclean + ' -size ' 
      cmd += str(imsizein) + ' ' + str(imsizein) + ' -channels-out ' + str(channelsout) + ' -padding 1.4 -predict ' 

      
      cmd += '-name ' + imageout + ' -scale ' + str(pixsizein) + 'arcsec ' + msliststring
      print('PREDICT STEP: ', cmd)
      os.system(cmd)


def subtractcompact(mslist, imageout, pixsize, imsize, minuv, channelsout=6, niter=15000, robust=-0.5, outcolumn='DIFFUSE_SUB'):

   # some setup
#   username = os.getlogin()
#   if username == 'rvweerenold':
 #    makemask = '/net/para10/data1/shimwell/software/killmsddf/new-install/DDFacet/SkyModel/MakeMask.py'
#   else:
    makemask = 'MakeMask.py'
     
    makeimage(mslist, imageout +'_compact', pixsize, imsize, channelsout=channelsout, niter=niter, robust=robust, minuv=minuv, predict=False)

   # make a mask
    imagename  = imageout +'_compact' + '-MFS-image.fits'
    cmdm  = makemask + ' --Th=3.0 --RestoredIm=' + imagename
    print(cmdm)
    os.system(cmdm)
    fitsmask = imagename + '.mask.fits'

   # re-image with mask
    makeimage(mslist, imageout +'_compactmask', pixsize, imsize, channelsout=channelsout, niter=niter, robust=-0.5, minuv=minuv, multiscale=True, predict=True, fitsmask=fitsmask, deepmultiscale=False)
    
   # now subtract the columns 
    
    for ms in mslist:
        ts  = pt.table(ms, readonly=False)
        colnames = ts.colnames()
        if outcolumn not in colnames:
            desc = ts.getcoldesc('DATA')
            desc['name']=outcolumn
            ts.addcols(desc)
            ts.close() # to write results
        else:
            print(outcolumn, ' already exists')
            ts.close()

    for ms in mslist:
        ts  = pt.table(ms, readonly=False)
        colnames = ts.colnames()
        if 'CORRECTED_DATA' in colnames:
            data = ts.getcol('CORRECTED_DATA')
        else:
            data = ts.getcol('DATA')
            model = ts.getcol('MODEL_DATA')
            ts.putcol(outcolumn,data-model)
            ts.close()
        return



# some setup
#username = os.getlogin()
#if username == 'rvweerenold':
#  makemask = '/net/para10/data1/shimwell/software/killmsddf/new-install/DDFacet/SkyModel/MakeMask.py'
#else:
makemask = 'MakeMask.py'


parser = argparse.ArgumentParser(description='Make images from extraction run. Requires working version of the DR2-pipeline software and WSClean (Oct 2018 or newer')
parser.add_argument('-b','--boxfile', help='optional boxfile to set imsize automatically', type=str)
#parser.add_argument('--fitsmask', help='fitsmask for deconvolution, if not provided use automasking', type=str)
parser.add_argument('--imsize', help='image size, you can take it from selfcal.log', type=int)
parser.add_argument('-n', '--niter', help='niter, default=25000', default=25000, type=int)
#parser.add_argument('--robust', help='Briggs robust paramter, default=-0.5', default=-0.5, type=float)
parser.add_argument('--channelsout', help='channelsout, default=6', default=6, type=int)
parser.add_argument('--minuv', help='inner uv-cut for image in lambda, default=80', default=80., type=float)
parser.add_argument('--nodosub', help='in addition to normal imaging, also subtract compact sources', action='store_true')
parser.add_argument('--onlydosub', help='only produce images subtracting compact sources', action='store_true')
parser.add_argument('--pixelscale', help='pixels size in arcsec, deafult=1.5', default=1.5, type=float)
parser.add_argument('--sourceLLS', help='size in Mpc of diffuse emission for uvcut, default=0.4', default=0.4, type=float)
parser.add_argument('--z', help='redshift of cluster, not required if --nodosub is used', default=-1.0, type=float)
parser.add_argument('-i','--imagename', help='imagename, default=image', required=True, type=str)
parser.add_argument('--maskthreshold', help='threshold for MakeMask.py, default=3.0', default=3.0, type=int)
parser.add_argument('ms', nargs='*', help='msfile(s)')

args = vars(parser.parse_args())

minuv    = args['minuv']  
pixsize  = args['pixelscale']  
niter    = args['niter']  
mslist   = sorted(args['ms'])
imageout = args['imagename']
#imsize   = args['imsize']

if args['boxfile'] == None and args['imsize'] == None:
  print('Incomplete input detected, either boxfile or imsize is required')
  sys.exit()
if args['boxfile'] != None and args['imsize'] != None:
  print('Wrong input detected, both boxfile and imsize are set')
  sys.exit()

if args['boxfile'] != None:
  imsize   = str(getimsize(args['boxfile'], args['pixelscale']))
if args['imsize'] != None:
  imsize = str(args['imsize']) 


#if args['imsize'] == None:
#    print 'Error: imsize not provided'
#    sys.exit()


if args['z'] < 0: # if no redshift provided try to find it automatically
  customSimbad = Simbad()
  customSimbad.add_votable_fields('z_value')
  #obj=customSimbad.query_object( args['imagename'] )
  #print obj
  #sys.exit()
  try:
    obj=customSimbad.query_object( args['imagename'] )
  #except: # in this case the cluster name is recognized by Simbad as no error was raised above

    if obj['Z_VALUE'][0] > 0.0:
        print('Found redshift',  obj['Z_VALUE'][0])
        args['z'] = obj['Z_VALUE'][0]
        os.system('echo ' + str(args['z']) + ' > redshift-used.log')
    else:
        print('Warning: Cluster known but redshift not found, assuming z=0.2')
        args['z'] = 0.2
        os.system('echo ' + str(args['z']) + ' > redshift-used.log')
  except:
    print('Cluster name is not known by Simbad' )

if not args['nodosub']:
  
  if args['z'] < 0:
    print('You need provide a redshift, none was given')
    sys.exit()
    
  minuv_forsub = compute_uvmin(args['z'], sourceLLS=args['sourceLLS'])

  subtractcompact(mslist, imageout, pixsize, imsize, minuv_forsub, channelsout=args['channelsout'],\
                  niter=np.int(niter/1.25), robust=-0.5, outcolumn='DIFFUSE_SUB')





if not args['nodosub']:  

    #  -----------------------------------------------------------------
    #  --- make the taper 50 kpc image, compact source subtracted ----
    #  -----------------------------------------------------------------

    makeimage(mslist, imageout +'_subROBUST-0.5TAPER50kpc', pixsize, imsize, channelsout=args['channelsout'], niter=adjustniter_for_taper(compute_taper(args['z'],50.), niter), robust=-0.5, minuv=minuv, predict=False, column='DIFFUSE_SUB', uvtaper=compute_taper(args['z'],50.))

    # make a mask
    imagename  = imageout +'_subROBUST-0.5TAPER50kpc' + '-MFS-image.fits'
    cmdm  = makemask + ' --Th='+ str(args['maskthreshold']) + ' --RestoredIm=' + imagename
    print(cmdm)
    os.system(cmdm)
    fitsmask = imagename + '.mask.fits'

    # re-image with mask
    makeimage(mslist, imageout +'_masksubROBUST-0.5TAPER50kpc', pixsize, imsize, channelsout=args['channelsout'], niter=adjustniter_for_taper(compute_taper(args['z'],50.), niter), robust=-0.5, minuv=minuv, multiscale=True, predict=False, fitsmask=fitsmask, deepmultiscale=False, column='DIFFUSE_SUB', uvtaper=compute_taper(args['z'],50.))    


if (not args['onlydosub'] == True) and (args['z'] >0.): # otherwise cannot computer compute_taper(z)


    #  --------------------------------------------------
    #  --- make the taper 50 kpc image ----
    #  --------------------------------------------------
    makeimage(mslist, imageout + '_ROBUST-0.5TAPER50kpc', pixsize, imsize, channelsout=args['channelsout'], niter=adjustniter_for_taper(compute_taper(args['z'],50.), niter), robust=-0.5, minuv=minuv, predict=False, uvtaper=compute_taper(args['z'],50.))

    # make a mask
    imagename  = imageout +  '_ROBUST-0.5TAPER50kpc' +'-MFS-image.fits'
    cmdm  = makemask + ' --Th='+ str(args['maskthreshold']) + ' --RestoredIm=' + imagename
    print(cmdm)
    os.system(cmdm)
    fitsmask = imagename + '.mask.fits'

    # re-image with mask
    makeimage(mslist, imageout +'_maskROBUST-0.5TAPER50kpc', pixsize, imsize, channelsout=args['channelsout'], niter=adjustniter_for_taper(compute_taper(args['z'],50.), niter), robust=-0.5, minuv=minuv, multiscale=True, predict=False, fitsmask=fitsmask, deepmultiscale=False,uvtaper=compute_taper(args['z'],50.))



cleanup()

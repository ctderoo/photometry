'''
Purpose: python module facilitating photometry using ds9 regions and astropy/photutils.

Written by: Caroline Roberts & Casey DeRoo
'''

#setup
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style
plt.style.use(astropy_mpl_style)
from astropy.io import fits
from matplotlib.colors import LogNorm

import pdb

import photutils

#from photutils import CircularAperture
#from photutils import CircularAnnulus
#from photutils import aperture_photometry

def create_aperture_from_ds9_region(region_file):
    '''
    Loads a ds9 region ('reg') file saved in physical format and returns an astropy aperture instance.
    '''
    
    # 
    f = open(region_file,'rb')
    for i in range(3):
        f.readline()
    region_text = f.read()

    pdb.set_trace()
    regions = region_text.split('\n')

    src_apertures,bkg_apertures = [],[] 

    for region in regions[:-1]:
        # This is a poor example of string handling, but it works.
        value_text = region.split('(')[1].split(')')[0]
        pos_x,pos_y,r_in,r_out = np.array(value_text.split(','),dtype = float)
        
        src_apertures.append(photutils.CircularAperture((pos_x,pos_y), r=r_in))
        bkg_aperture.append(photutils.CircularAnnulus((pos_x,pos_y), r_in=r_in, r_out=r_out))

    return src_apertures,bkg_apertures

def perform_aperture_photometry(fits_file, ann_reg_fn, exp_time = None):
    
    # Opening the fits file itself and getting the image data.
    hdu = fits.open(fits_file)
    image = hdu[0].data
    # This is poor approximation to what's happening in terms of errors -- if each 
    # detector value represents some number of photons, then the error will be proportional
    # to the square root of that number of photons. But this shouldn't be interpreted as 
    # rigorously correct.
    error = np.sqrt(image)

    # Get the exposure time from the .fits file.    
    if exp_time is None:
        try:
            exp_time = hdu[0].header['EXPTIME']
        except:
            print('Exposure time is not defined in the .fits file header. Please specify during function call.')

    src_aps,bkg_aps = create_aperture_from_ds9_region(ann_reg_fn)

    src_table = photutils.aperture_photometry(img, src_aps, error=error)
    bkg_table = photutils.aperture_photometry(img, bkg_aps, error=error)

    # Calculating the mean background and associated error.
    bkg_mean = bkg_table['aperture_sum'] / bkg_aps.area
    bkg_mean_err = bkg_table['aperture_sum_err'] / bkg_aps.area

    # Calculating the source aperture brightness, minus the expected background values.
    src_sum = src_table['aperture_sum'] - bkg_mean*src_aps.area
    src_error = np.sqrt(src_table['aperture_sum_err']**2 + (bkg_aps.area*bkg_mean_err)**2)

    return src_sum,src_error


from _future__ import division, print_function
import numpy as np 
import pandas as pd

import tic_code

## constants
# original values from email by Josh Pepper from 06/03/2016
tess_scale = 20.25 / 3600.0        # arcsec / pixel --> deg / pixel
tess_fwhm  = 1.88  * tess_scale    # new fwhm from camera image (sigma = 0.80)
tess_aper  = 4.0  * tess_scale     # 4 pixel aperture
tess_srad  = 10.0 * tess_scale     # 10 pixel search radius 
tess_sigma = tess_fwhm / (2.0 * sqrt(2.0 * log(2))) # by definition

def get_something(starParams):

    if not np.isfinite(starParams.tmage):
        starParams.tmage = 0.25

    aper, apere = tic_code.aperture(starParams.tessmag, starParams.tmage)
    aper  *= tess_scale
    apere *= tess_scale

    tflux, tfluxe = tic_code.tmag2flux(starParams.tessmag, starParams.tmage)

    # do something to get nearby stars




def get_nearby(srad, TIC):
    """
    find targets in the TIC that are within a given distance
    srad is in degrees
    TIC is a pandas dataframe
    returns the indices of the matching rows
    """

    # only works for if search radius is less than 1 degreess
    assert srad < 1

    #do something coarse first
    rows = np.where(TIC.ra_deg -  )
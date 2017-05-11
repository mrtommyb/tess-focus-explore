from __future__ import division, print_function
import numpy as np 
import pandas as pd

from tqdm import tqdm

from scipy.special import erf

# ## constants
# # original values from email by Josh Pepper from 06/03/2016

# tess_scale = 20.25 / 3600.0  # arcsec / pixel --> deg / pixel
# # tess_fwhm  = 2.0  * tess_scale  # 2 pixel
# tess_fwhm = 1.88 * tess_scale  # new fwhm from camera image (sigma = 0.80)
# tess_aper = 4.0 * tess_scale  # 4 pixel aperture
# tess_srad = 10.0 * tess_scale  # 10 pixel search radius
# tess_sigma = tess_fwhm / (2.0 * np.sqrt(2.0 * np.log(2)))  # by definition




class CalculateContamination(object):

    def __init__(self):

        self.pixScale = 20.25 / 3600.0
        self.tessmagError = 0.2  # give every star the same uncertainty

    def findContamSingle(self, starParams, TIC, **kwargs):
        self.starParams = starParams.copy()
        self.tessid = self.starParams.loc[:, 'TICID']

        if 'nearbyRad' in kwargs:
            self.nearbyRad = self.pixScale * kwargs['nearbyRad']
        else:
            self.nearbyRad = self.pixScale * 15

        if 'psfFwhm' in kwargs:
            self.psfFwhm = kwargs['psfFwhm'] * self.pixScale
        else:
            self.psfFwhm = 1.88 * self.pixScale
        self.psfSigma = self.psfFwhm / (2.0 * np.sqrt(2.0 * np.log(2)))

        self.find_nearby_stars(TIC)
        self.nearbyCat = TIC.loc[self.nearbyStars, :].copy()
        self.nearbyCat.loc[:, 'dist'] = self.dist

        self.calc_contam()

    def find_nearby_stars(self, TIC):
        """
        find targets in the TIC that are within a given distance
        TIC is a pandas dataframe
        returns the indices of the matching rows
        """

        dist = angSepVincenty(self.starParams.loc[:, 'RA_DEG'].values,
                              self.starParams.loc[:, 'DEC_DEG'].values,
                              TIC.loc[:, 'RA_DEG'],
                              TIC.loc[:, 'DEC_DEG'])

        self.nearbyStars = dist < self.nearbyRad

        # remove the star itself
        # search 0.05 arcsec
        self.nearbyStars[np.abs(dist) < (0.05 / 3600.)] = False
        self.dist = dist[dist < self.nearbyRad]

    def calc_tflux(self):
        aper = aperture(self.starParams.loc[:, 'TESSMAG'].values)
        aper *= self.pixScale
        self.pixAper = aper

        tflux = tmag2flux(self.starParams.loc[:, 'TESSMAG'].values)
        assert not np.any(np.isnan(tflux))
        self.starParams.loc[:, 'tflux'] = tflux

        tflux_nearby = tmag2flux(
            self.nearbyCat.loc[:, 'TESSMAG'].values)
        self.nearbyCat.loc[:, 'tflux'] = tflux_nearby

    def calc_contam(self):
        self.calc_tflux()

        # i'm rewriting the tic code here
        xb = yb = self.pixAper / 2.
        x0 = angSepVincenty(self.starParams.loc[:, 'RA_DEG'].values,
                            self.starParams.loc[:, 'DEC_DEG'].values,
                            self.nearbyCat.loc[:, 'RA_DEG'],
                            self.starParams.loc[:, 'DEC_DEG'].values).values
        y0 = angSepVincenty(self.starParams.loc[:, 'RA_DEG'].values,
                            self.starParams.loc[:, 'DEC_DEG'].values,
                            self.starParams.loc[:, 'RA_DEG'].values,
                            self.nearbyCat.loc[:, 'DEC_DEG']).values
        sq2 = np.sqrt(2)
        s = self.psfSigma

        contx = erf((xb + x0) / (sq2 * s)) + erf((xb - x0) / (sq2 * s))
        conty = erf((yb + y0) / (sq2 * s)) + erf((yb - y0) / (sq2 * s))

        cont = 0.25 * contx * conty

        cflx = cont * self.nearbyCat.loc[:,'tflux']

        self.totalContamFlux = np.sum(cflx)
        self.fluxRatio = self.totalContamFlux / self.starParams.loc[:, 'tflux']


def angSepVincenty(ra1, dec1, ra2, dec2):
    """
    Vincenty formula for distances on a sphere
    """
    ra1_rad = np.radians(ra1)
    dec1_rad = np.radians(dec1)
    ra2_rad = np.radians(ra2)
    dec2_rad = np.radians(dec2)

    sin_dec1, cos_dec1 = np.sin(dec1_rad), np.cos(dec1_rad)
    sin_dec2, cos_dec2 = np.sin(dec2_rad), np.cos(dec2_rad)
    delta_ra = ra2_rad - ra1_rad
    cos_delta_ra, sin_delta_ra = np.cos(delta_ra), np.sin(delta_ra)

    diffpos = np.arctan2(np.sqrt((cos_dec2 * sin_delta_ra) ** 2 +
                         (cos_dec1 * sin_dec2 -
                         sin_dec1 * cos_dec2 * cos_delta_ra) ** 2),
                         sin_dec1 * sin_dec2 + cos_dec1 * cos_dec2 *
                         cos_delta_ra)

    return np.degrees(diffpos)


def aperture(tmag):
    tmag[tmag < 4.0] = 4.0
    npix = np.zeros_like(tmag)

    npix = (274.2898 - 77.7918 * tmag +
            7.7410 * tmag**2 - 0.2592 * tmag**3)

    npix[npix < 1] = 1
    aper1 = np.sqrt(npix)
    aper2 = 2 * np.sqrt(npix / np.pi)

    return (aper1 + aper2) / 2


def tmag2flux(tmag):
    """
    convert tmag to flux
    the tic code uses imag
    """
    flux = 2635.0 * np.power(10.0, -tmag / 2.5)
    return flux


if __name__ == '__main__':
    # get rid of those annyoing SettingWithCopyWarning
    import warnings

    with warnings.catch_warnings():
        warnings.simplefilter('error')


        # ticid = 88416102
        ticloc = 'tic_at_50_minus30.h5'
        TIC = pd.read_hdf(ticloc)
        # starParams = TIC.loc[TIC.loc[:, 'TICID'] == ticid]

        # C = CalculateContamination()

        # C.findContamSingle(starParams, TIC)#, psfFwhm=5.0)
        # print(C.fluxRatio)

        #let's trying calculating for high priority stars
        pr = TIC.PRIORITY.argsort()[::-1]

        width_factors = [1.0, 1.1, 1.2, 1.3, 1.5, 2.0, 3.0, 5.0] 

        dfout = pd.DataFrame(np.zeros([len(width_factors),4000]))

        psf0 = 1.88
        for i, w in enumerate(width_factors):
            for j, ticid in tqdm(
                    enumerate(TIC.loc[:, 'TICID'].iloc[pr][0:4000])):
                starParams = TIC.loc[TIC.loc[:, 'TICID'] == ticid]
                C = CalculateContamination()
                C.findContamSingle(starParams, TIC, psfFwhm=psf0 * w)
                fr1 = C.fluxRatio.copy()
                dfout.iloc[i,j] = fr1.values



                # C = CalculateContamination()
                # C.findContamSingle(starParams, TIC, psfFwhm=3.6)
                # fr2 = C.fluxRatio.copy()
                # # # print([fr1.values, fr2.values])
                # # if ((1-fr1)/(1-fr2)).values > 1.01:
                # #     print(((1-fr1)/(1-fr2)).values)
            #     # #     i +=1

            # print()
            # print(i)






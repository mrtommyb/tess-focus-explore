import numpy as np
import pandas as pd

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


if __name__ == '__main__':

    # TIC is stored at
    tic_dir = '/Users/tom/Dropbox/TIC4/CTL/'

    # we are going to make a version of the TIC
    # star is a 6 degrees radius circle

    # we are going to center the coordinates on
    # ra = 50 degrees
    # dec = -30 degrees
    # that means we use the files:
    tic_file = '02-04.csv'
    usecols = [0, 1, 2, 3, 4, 5, 6, 7, 20]

    header_file = 'header.txt'
    h = pd.read_csv(tic_dir + header_file, usecols=usecols)

    # i don't like some column names
    h.rename(columns={'RA': 'RA_DEG',
                      'DEC': 'DEC_DEG'}, inplace=True)

    TIC = pd.read_csv(tic_dir + tic_file, header=0, usecols=usecols)

    TIC.columns = h.columns

    deg_from_center = angSepVincenty(TIC.RA_DEG, TIC.DEC_DEG, 50, -30)





    savetic = TIC[deg_from_center < 12.].copy()
    print('saving file with {} lines'.format(savetic.size))

    savetic.to_hdf('tic_at_50_minus30.h5', key='data', mode='w')

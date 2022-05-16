import os
import sys
import re
import logging
import numpy as np
import pandas as pd

from datetime import datetime

import magnetopost.util as mputil

#from datetick import datetick
def datetick(arg):
    pass

OVERWRITE_CACHE=False
rootdir = '/home/gary/media_sunspot'
ned = ('north','east','down')

if '-O' in sys.argv:
    OVERWRITE_CACHE=True
for arg in sys.argv:
    if '--rootdir=' == arg[:len('--rootdir=')]:
        rootdir = arg[len('--rootdir='):]

def read_ccmc_printout(filename):
    '''
    Reads files with text like ::
# Data printout from CCMC-simulation: version 1D-1.3
# Data type:  CalcDeltaB  PP
# Run name:   Gary_Quaresima_20210809_PP_1
# Missing data:  -1.09951e+12
# year month day hour minute second lon lat alt B_north B_east B_down B_north_mag B_east_mag B_down_mag B_north_fac B_east_fac B_down_fac B_north_iono,SigP B_east_iono,SigP B_down_iono,SigP B_north_iono,SigH B_east_iono,SigH B_down_iono,SigH 
# year month day hr min s deg deg km nT nT nT nT nT nT nT nT nT nT nT nT nT nT nT 
  2.01900E+03  9.00000E+00  2.00000E+00  4.00000E+00  1.00000E+01  0.00000E+00  7.28150E+01  1.89070E+01  1.00000E+00  6.62460E+00 -2.46310E+00  2.40840E+00  6.11990E+00 -3.47180E+00  2.15050E+00 -9.67200E-01  2.64940E+00  1.41000E-01  4.60600E-01  1.05700E-01 -2.07000E-02 -1.93000E-01  3.18400E-01 -8.70000E-01
  ...
    '''
    with open(filename) as f:
        head = [f.readline() for _ in range(5)]
        arr = np.loadtxt(f)
    headers = tuple(head[-1].split(' ')[1:-1])
    return headers, arr


def extract_from_swmf_magnetometer_files(info, surface_location, n_steps=None):
    """
    This reads the values for the different contributions to the
    magnetic field that SWMF can calculate natively.
    These are stored in the text files with names of the form mag_grid____.out

    NOTE:

    lats = np.unique(data[:,1])
      -> array([-87.5       , -86.49425287, -85.48850575, -84.48275862,
                ...
                85.48850575,  86.49425287,  87.5       ])
    lats[1]-lats[0]
      -> 1.0057471260000028
    lats[6]-lats[5]
      -> 1.0057471270000065
    lats.shape
      -> (175,)
    87.5-(-87.5)
      -> 175.0
    175./(175-1)
      -> 1.0057471264367817

    lons = np.unique(data[:,0])
      -> array([  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   8.,   9.,  10.,
                ...
                352., 353., 354., 355., 356., 357., 358., 359.])
    lons.shape
      -> (360,)
    (359.-0.)/(360-1)
      -> 1.0

    so while lons in steps of 1. , lats are in steps of 175./174 
    (to some decimal aprox that is not completely consistent between elements)
    so someone probably screwed up the linspace and meant to use (176,)
    """

    logging.info("Extracting point dB information from mag_grid files in {}".format(info["dir_run"]))

    cachepath = f'{info["dir_run"]}/derived/extract_from_swmf_magnetometer_files'
    if not OVERWRITE_CACHE:
        try:
            dBMhd = pd.read_pickle(f'{cachepath}/dBMhd-{surface_location}.pkl')
            dBFac = pd.read_pickle(f'{cachepath}/dBFac-{surface_location}.pkl')
            dBHal = pd.read_pickle(f'{cachepath}/dBHal-{surface_location}.pkl')
            dBPed = pd.read_pickle(f'{cachepath}/dBPed-{surface_location}.pkl')
            logging.info('Using cached calculations found in {}'.format(cachepath))
            return dBMhd, dBFac, dBHal, dBPed
        except FileNotFoundError:
            pass

    dBMhd = pd.DataFrame()
    dBFac = pd.DataFrame()
    dBHal = pd.DataFrame()
    dBPed = pd.DataFrame()
    t = 1
    for time, cdfname in info['files']['magnetosphere'].items():

        if n_steps is not None and t > n_steps:
            break
        t = t + 1
        
        if not '3d__var_' in cdfname:
            raise RuntimeError('FILENAME DOES NOT FOLLOW 3d__var_ convention')
        filename = cdfname.replace('.cdf','')
        filename = re.sub('3d__var_._','mag_grid_', f'{filename[:-8]}.out')

        logging.info("Reading {}".format(filename))
        with open(filename, 'r') as f:
            first = f.readline()
            second = f.readline()
            third = f.readline()
            headerline = f.readline()
        arr = np.genfromtxt(filename, skip_header=4)

        headers = tuple(headerline[:-1].split(' '))
        csyst = first[19:22]

        assert(first[:19] == 'Magnetometer grid (')
        assert(headers[0]=='Lon' and headers[1]=='Lat')
        assert(headers[5:] == ( 'dBnMhd', 'dBeMhd', 'dBdMhd',
                                'dBnFac', 'dBeFac', 'dBdFac',
                                'dBnHal', 'dBeHal', 'dBdHal',
                                'dBnPed', 'dBePed', 'dBdPed') )

        _r, LAT, LON = mputil.GetMagnetometerCoordinates(surface_location, time, csyst, 'sph')
        if not 0.99 < _r < 1.01: raise ValueError

        if LON < 0:
            LON = LON + 360.
        if csyst != 'GEO':
            print(f'csyst={csyst}')
        if abs(LAT - 18.907) > 1e-8 or LON != 72.815:
            print(f'WARNING: LAT,LON = {LAT},{LON}')

        Tr = np.all([LON-0.5 <= arr[:, 0], 
                     arr[:, 0] <= LON+0.5, LAT-0.5 <= arr[:, 1],
                     arr[:, 1] <= LAT+0.5], axis=0)
        k = np.where(Tr==True)[0][0]

        dBMhd = dBMhd.append(pd.Series(data=arr[k,5 :8 ], index=ned, name=datetime(*time)))
        dBFac = dBFac.append(pd.Series(data=arr[k,8 :11], index=ned, name=datetime(*time)))
        dBHal = dBHal.append(pd.Series(data=arr[k,11:14], index=ned, name=datetime(*time)))
        dBPed = dBPed.append(pd.Series(data=arr[k,14:17], index=ned, name=datetime(*time)))

    os.makedirs(cachepath, exist_ok=True)

    logging.info("Saving parsed content of dB files to {cachepath}")
    dBMhd.to_pickle(f'{cachepath}/dBMhd-{surface_location}.pkl')
    dBFac.to_pickle(f'{cachepath}/dBFac-{surface_location}.pkl')
    dBHal.to_pickle(f'{cachepath}/dBHal-{surface_location}.pkl')
    dBPed.to_pickle(f'{cachepath}/dBPed-{surface_location}.pkl')

    return dBMhd, dBFac, dBHal, dBPed


def extract_from_swmf_ccmc_printout_file(info, surface_location, n_steps=None):

    year = list(info['files']['magnetosphere'].keys())[0][0]
    filename = f'{info["dir_run"]}/derived/{year}_{surface_location}_pointdata.txt'

    logging.info("Reading point dB information from " + filename)

    headers, arr = read_ccmc_printout(filename)
    assert( headers[10:] ==('sumBn', 'sumBe', 'sumBd',
                             'dBn', 'dBe', 'dBd',
                            'facdBn', 'facdBe', 'facdBd',
                            'JhdBn', 'JhdBe', 'JhdBd',
                            'JpBn', 'JpBe', 'JpBd') )

    if n_steps is not None:
        arr = arr[0:n_steps,:]

    times = np.array(arr[:,0:6], dtype=int)
    dtimes = [datetime(*time) for time in times]

    dBMhd = pd.DataFrame(data=arr[:,13:16], columns=ned, index=dtimes)
    dBFac = pd.DataFrame(data=arr[:,16:19], columns=ned, index=dtimes)
    dBHal = pd.DataFrame(data=arr[:,19:22], columns=ned, index=dtimes)
    dBPed = pd.DataFrame(data=arr[:,22:25], columns=ned, index=dtimes)

    return dBMhd, dBFac, dBHal, dBPed


def extract_from_CalcDeltaB_file(filename):

    logging.info("Reading dB information from " + filename)

    headers, arr = read_ccmc_printout(filename)
    assert( headers[9:] ==('B_north', 'B_east', 'B_down',
                           'B_north_mag', 'B_east_mag', 'B_down_mag',
                           'B_north_fac', 'B_east_fac', 'B_down_fac',
                           'B_north_iono,SigP', 'B_east_iono,SigP', 'B_down_iono,SigP',
                           'B_north_iono,SigH', 'B_east_iono,SigH', 'B_down_iono,SigH') )
    times = np.array(arr[:,0:6], dtype=int)

    dtimes = [datetime(*time) for time in times]
    B_mag = pd.DataFrame(data=arr[:,9:12], columns=ned, index=dtimes)
    B_fac = pd.DataFrame(data=arr[:,12:15], columns=ned, index=dtimes)
    B_ionoSigP = pd.DataFrame(data=arr[:,15:18], columns=ned, index=dtimes)
    B_ionoSigH = pd.DataFrame(data=arr[:,18:21], columns=ned, index=dtimes)
    return B_mag, B_fac, B_ionoSigH, B_ionoSigP # flipped order H an P


def extract_from_magnetopost_files(info, surface_location, n_steps=None):

    msph_times = info['files']['magnetosphere'].keys()
    iono_times = info['files']['ionosphere'].keys()

    msph_dtimes = [datetime(*time) for time in msph_times]
    iono_dtimes = [datetime(*time) for time in iono_times]

    if n_steps is not None:
        msph_dtimes = msph_dtimes[0:n_steps]
        iono_dtimes = iono_dtimes[0:n_steps]

    def get(ftag, dtimes):
        df = pd.DataFrame()
        infile = f'{info["dir_run"]}/derived/timeseries/{ftag}-{surface_location}.npy'
        logging.info(f"Reading {infile}")
        dB = np.load(infile)

        df['north'] = pd.Series(data=dB[:,0], index=dtimes)
        df['east']  = pd.Series(data=dB[:,1], index=dtimes)
        df['down']  = pd.Series(data=dB[:,2], index=dtimes)
        return df

    bs_msph     = get('bs_msph', msph_dtimes)
    bs_fac      = get('bs_fac', msph_dtimes)

    bs_hall     = get('bs_hall', iono_dtimes)
    bs_pedersen = get('bs_pedersen', iono_dtimes)
    bs_hall     = get('bs_hall', iono_dtimes)

    cl_msph     = get('cl_msph', msph_dtimes)
    helm_outer  = get('helm_outer', msph_dtimes)
    helm_rCurrents_gapSM  = get('helm_rCurrents_gapSM', msph_dtimes)

    probe = get('probe', msph_dtimes)

    return bs_msph, bs_fac, bs_hall, bs_pedersen, cl_msph,helm_outer, helm_rCurrents_gapSM, probe


def extract_all():
    rundir = '{rootdir}/DIPTSUR2/'
    surface_location = 'colaba'

    dBMhd, dBFac, dBHal, dBPed = extract_from_swmf_magnetometer_files(rundir, surface_location)
    B_mag, B_fac, B_ionoSigH, B_ionoSigP = extract_from_CalcDeltaB_file('pointdata_751815008468.txt')
    bs_msph, bs_fac, bs_hall, bs_pedersen = extract_from_magnetopost_files(rundir, surface_location)


#def main():
    #surface_point('DIPTSUR2','colaba')
    #msph_point('DIPTSUR2','GMpoint1')
    #surface_point('SWPC_SWMF_052811_2','YKC')
    #msph_point('SWPC_SWMF_052811_2','GMpoint6')


#if __name__ == '__main__':
#    main()

import numpy as np
from numba import njit
#import pandas as pd
import datetime

import os
import sys

from magnetopost import util
from magnetopost.units_and_constants import phys

@njit # Biot Savart integrand (4pi taken outside)
def _integrand_bs(x, y, z, obs_point, jx, jy, jz):
    integrand = np.empty((3,),dtype=np.float32); integrand[:]=np.nan

    r_x = obs_point[0] - x
    r_y = obs_point[1] - y
    r_z = obs_point[2] - z
    r = np.sqrt(r_x**2 + r_y**2 + r_z**2)
    if r > 1e-5:
        integrand[2] = jx*r_y - jy*r_x
        integrand[0] = jy*r_z - jz*r_y
        integrand[1] = jz*r_x - jx*r_z
        integrand[:] =  integrand[:]/r**3
    else:
        integrand[:] = 0.

    return integrand

@njit
def _integral_bs(X, Y, Z, obs_point, JX, JY, JZ, Measure):
    ret = np.zeros((3,), dtype=np.float32)
    for i in range(X.size):

        # 'measure' refers to integration measure ( dx, dV, dA, whatever apropriate)
        ret[:] = ret[:] + _integrand_bs(X[i], Y[i], Z[i], obs_point, JX[i], JY[i], JZ[i]) * Measure[i]

    return ret/(4.*np.pi)


def slice_bs_pedersen(run, time, ie_slice, obs_point):
    obs_point_str = obs_point
    if obs_point == "origin":
        obs_point = np.zeros(3)
    else:
        obs_point = util.GetMagnetometerCoordinates(obs_point_str, time, 'SM', 'car')

    data_arr, varidx, units = ie_slice

    X  = data_arr[varidx['X'], :]
    Y  = data_arr[varidx['Y'], :]
    Z  = data_arr[varidx['Z'], :]

    # note, surface current density (not vol current density)
    KX = data_arr[varidx['Ex'],:] * data_arr[varidx['SigmaP'],:]
    KY = data_arr[varidx['Ey'],:] * data_arr[varidx['SigmaP'],:]
    KZ = data_arr[varidx['Ez'],:] * data_arr[varidx['SigmaP'],:]

    Measure = data_arr[varidx['measure'],:]

    print(units)
    #assert( units['SigmaP'] == "" and units['Ex'] == "" and units['X'] == "R" )
    scalefact = phys['mu0']*phys['Siemens']*phys['mV']/phys['m']

    integral = scalefact*_integral_bs(X,Y,Z,obs_point,KX,KY,KZ,Measure)

    outname = '/home/gquaresi/timeseries/slices/' \
        + f'BS_pedersen_{util.Tstr(time)}_obs_point={obs_point_str}.npy'

    np.save(outname, integral)


def slice_integral_bs_bulkiono(ie_slice, obs_point):
    data_arr, varidx, units = ie_slice

    X  = data_arr[varidx['X'] , :]
    Y  = data_arr[varidx['Y'] , :]
    Z  = data_arr[varidx['Z'] , :]
    Jx = data_arr[varidx['Jx'], :]
    Jy = data_arr[varidx['Jy'], :]
    Jz = data_arr[varidx['Jz'], :]

    Measure = data_arr[varidx['measure'],:]

    assert( units['Jx'] == "`mA/m^2" and units['X'] == "R" )
    scalefact = phys['mu0']*phys['muA']/(phys['m']**2)

    integral = scalefact*_integral_bs(X,Y,Z,obs_point,Jx,Jy,Jz,Measure)
    print(integral)



if __name__ == '__main__':
    slice_integral_bs_pedersen('DIPTSUR2', (2019,9,2,4,11,0), (0,0,0))


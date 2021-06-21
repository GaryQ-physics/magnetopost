import numpy as np
from numba import njit
import datetime
import pandas as pd

from magnetopost import util
from magnetopost.units_and_constants import phys

@njit
def _jit_B_biotsavart(ms_slice, x0, rcut):
    integral = np.zeros((3,),dtype=np.float32)

    for ind in range(ms_slice.data_arr.shape[0]):
        measure = ms_slice.data_arr[ind, ms_slice.varidx['measure']]
        x       = ms_slice.data_arr[ind, ms_slice.varidx['x']]
        y       = ms_slice.data_arr[ind, ms_slice.varidx['y']]
        z       = ms_slice.data_arr[ind, ms_slice.varidx['z']]

        distanceSquared = x**2 + y**2 + z**2
        if distanceSquared < rcut**2:
            continue

        partials_bx = ms_slice.get_native_partial_derivatives(ind, 'bx')
        partials_by = ms_slice.get_native_partial_derivatives(ind, 'by')
        partials_bz = ms_slice.get_native_partial_derivatives(ind, 'bz')

        curl_B1_z = partial_by[0] - partial_bx[1]
        curl_B1_x = partial_bz[1] - partial_by[2]
        curl_B1_y = partial_bx[2] - partial_bz[0]

        r_x = x0[0] - x
        r_y = x0[1] - y
        r_z = x0[2] - z
        r = np.sqrt(r_x**2 + r_y**2 + r_z**2)
        if r < 1e-5:
            continue

        integrand = np.empty((3,),dtype=np.float32)
        integrand[2] = curl_B1_x*r_y - curl_B1_y*r_x
        integrand[0] = curl_B1_y*r_z - curl_B1_z*r_y
        integrand[1] = curl_B1_z*r_x - curl_B1_x*r_z
        integrand[:] =  integrand[:]/r**3

        integral[:] = integral[:] + measure*integrand
        #ret[i_eps,:] = ret[i_eps,:] + integrand

    #for i_eps in range(n_eps):
    #    ret[i_eps,:] = ( (unique_epsilons[i_eps]**3)/(4*np.pi) ) * ret[i_eps,:]
    return integral


def slice_bs_msph(run, time, ms_slice, obs_point):
    funcnameStr = 'bs_msph'

    if obs_point == "origin":
        x0 = np.zeros(3)
    else:
        x0 = util.GetMagnetometerCoordinates(obs_point, time, 'GSM', 'car')

    integral = _jit_B_biotsavart(ms_slice, x0, run['rCurrents'])
    
    outname = f'{run["rundir"]}/derived/timeseries/slices/' \
        + f'{funcnameStr}-{obs_point}-{util.Tstr(time)}.npy'
    np.save(outname, integral)


def stitch_bs_msph(run, times, obs_point):
    funcnameStr = 'bs_msph'

    integrals = []
    for time in times:
        outname = f'{run["rundir"]}/derived/timeseries/slices/' \
            + f'{funcnameStr}-{obs_point}-{util.Tstr(time)}.npy'

        integrals.append(np.load(outname))

    arr_name = f'{run["rundir"]}/derived/timeseries/' \
            + f'{funcnameStr}-{obs_point}.npy'
    arr = np.array(integrals)
    np.save(arr_name, arr)

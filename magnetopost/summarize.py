import numpy as np
from numba import njit
import datetime

import os
import sys

from magnetopost import util
from magnetopost.units_and_constants import phys
from magnetopost.gap_integrals import get_dipole_field

unitmu0 = np.float32( phys['mu0']*(phys['muA']/phys['m']**2) )#!!!!
_count = 0
_mean  = 1
_sndMo = 2
_std   = 3
_min   = 4
_max   = 5

@njit(error_model='numpy')
def _jit_stats_summary(ms_slice, rcut):
    #unique_epsilons = np.array([  0.0625,
    #                              0.1250,
    #                              0.2500,
    #                              0.5000,
    #                              1.0000,
    #                              2.0000,
    #                              4.0000,
    #                              8.0000 ])
    #n_eps = unique_epsilons.size
    #summary_arr = np.empty((nVarTot,n_eps,6), dtype=np.float32)
    #nBlock, nI, nJ, nK = DataArray.shape[1:]

    darr = ms_slice.data_arr
    varidx = ms_slice.varidx.copy()
    #ms_slice.total_measure

    nVar = len(varidx)
    varidx['div_b1'] = len(varidx)
    varidx['norm_curl_b1'] = len(varidx)
    varidx['div_b1_over_norm_curl_b1'] = len(varidx)
    varidx['norm_jR'] = len(varidx)
    varidx['jR_error'] = len(varidx)
    varidx['jR_fractional_error'] = len(varidx)
    nVarTot = len(varidx)
    #print(varidx)

    summary_arr = np.empty((nVarTot,6), dtype=np.float32)
    summary_arr[:,_count] = 0.
    summary_arr[:,_mean ] = 0.
    summary_arr[:,_sndMo] = 0.
    summary_arr[:,_std  ] = 0.
    summary_arr[:,_min  ] = np.inf
    summary_arr[:,_max  ] = -np.inf

    store = np.empty((nVarTot,), dtype=np.float32)
    for ind in range(darr.shape[0]):
        store[:]=np.nan

        distanceSquared = ( darr[ind, varidx['x']]**2 \
                          + darr[ind, varidx['y']]**2 \
                          + darr[ind, varidx['z']]**2 )
        if distanceSquared <= rcut**2:
            continue

        for iVar in range(nVar):
            store[iVar] = darr[ind, iVar]

        partials_b1x = ms_slice.get_native_partial_derivatives(ind, 'b1x')
        partials_b1y = ms_slice.get_native_partial_derivatives(ind, 'b1y')
        partials_b1z = ms_slice.get_native_partial_derivatives(ind, 'b1z')
        div_B1 = partials_b1x[0] + partials_b1y[1] + partials_b1z[2]
        curl_B1_z = partials_b1y[0] - partials_b1x[1]
        curl_B1_x = partials_b1z[1] - partials_b1y[2]
        curl_B1_y = partials_b1x[2] - partials_b1z[0]
        norm_curl_B1 = np.sqrt(curl_B1_x**2 + curl_B1_y**2 + curl_B1_z**2) 

        store[varidx['div_b1']] = div_B1
        store[varidx['norm_curl_b1']] = norm_curl_B1
        store[varidx['div_b1_over_norm_curl_b1']] = div_B1/norm_curl_B1

        jRx = (1./unitmu0)*curl_B1_x
        jRy = (1./unitmu0)*curl_B1_y
        jRz = (1./unitmu0)*curl_B1_z
        norm_jR = np.sqrt(jRx**2 + jRy**2 + jRz**2)
        jR_error = np.sqrt( (store[varidx['jx']]-jRx)**2 \
                          + (store[varidx['jy']]-jRy)**2 \
                          + (store[varidx['jz']]-jRz)**2 )

        store[varidx['norm_jR']] = norm_jR
        store[varidx['jR_error']] = jR_error
        store[varidx['jR_fractional_error']] = jR_error/norm_jR

        for iVar in range(nVarTot):
            if not np.isnan(store[iVar]):
                summary_arr[iVar,_count] = summary_arr[iVar,_count] + 1
                summary_arr[iVar,_mean ] = summary_arr[iVar,_mean ] + store[iVar]
                summary_arr[iVar,_sndMo] = summary_arr[iVar,_sndMo] + store[iVar]**2
                if store[iVar] < summary_arr[iVar,_min]:
                    summary_arr[iVar,_min] = store[iVar]
                if store[iVar] > summary_arr[iVar,_max]:
                    summary_arr[iVar,_max] = store[iVar]

    summary_arr[:,_mean ] = summary_arr[:,_mean ]/summary_arr[:,_count]
    summary_arr[:,_sndMo] = summary_arr[:,_sndMo]/summary_arr[:,_count]
    summary_arr[:,_std  ] = np.sqrt(summary_arr[:,_sndMo] - summary_arr[:,_mean]**2)#!! could suffer catastrophic cancelation
    return summary_arr


def slice_summary(run, time, ms_slice):
    funcnameStr = 'summary'

    summary = _jit_stats_summary(ms_slice, run['rCurrents'])

    outname = f'{run["rundir"]}/derived/timeseries/timesteps/' \
        + f'{funcnameStr}-{util.Tstr(time)}.npy'
    np.save(outname, summary)

def stitch_summary(run, times):
    funcnameStr = 'summary'

    summarys = []
    for time in times:
        outname = f'{run["rundir"]}/derived/timeseries/timesteps/' \
            + f'{funcnameStr}-{util.Tstr(time)}.npy'

        summarys.append(np.load(outname))

    arr_name = f'{run["rundir"]}/derived/timeseries/' \
            + f'{funcnameStr}.npy'
    arr = np.array(summarys)
    np.save(arr_name, arr)
    #with open(f'{arr}.head', 'w') as f:
    #    f.write(varidx)


if __name__ == '__main__':
    from magnetopost.model_patches import SWMF2
    #fn = '/home/gary/temp/3d__var_3_e20031120-070000-000.out'
    fn = '/home/gary/media_sunspot/SWPC_SWMF_052811_2/GM_CDF/3d__var_1_t00001001_n0002710.out.cdf'
    sl = SWMF2.get_ms_slice_class(fn)
    slice_summary({'rCurrents':1.}, (2019,9,2,4,11,0), sl)

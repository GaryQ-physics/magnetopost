import numpy as np
from numba import njit
import datetime

from hxform import hxform as hx
from magnetopost import util
from magnetopost.units_and_constants import phys

@njit
def _jit_B_biotsavart(ms_slice, x0, rcut, include):
    integral = np.zeros((3,),dtype=np.float32)

    for ind in range(ms_slice.data_arr.shape[0]):
        measure = ms_slice.data_arr[ind, ms_slice.varidx['measure']]
        x       = ms_slice.data_arr[ind, ms_slice.varidx['x']]
        y       = ms_slice.data_arr[ind, ms_slice.varidx['y']]
        z       = ms_slice.data_arr[ind, ms_slice.varidx['z']]

        distanceSquared = x**2 + y**2 + z**2
        if distanceSquared < rcut**2:
            continue

        if include is not None and not include[ind]: 
            continue

        partials_b1x = ms_slice.get_native_partial_derivatives(ind, 'b1x')
        partials_b1y = ms_slice.get_native_partial_derivatives(ind, 'b1y')
        partials_b1z = ms_slice.get_native_partial_derivatives(ind, 'b1z')

        curl_B1_z = partials_b1y[0] - partials_b1x[1]
        curl_B1_x = partials_b1z[1] - partials_b1y[2]
        curl_B1_y = partials_b1x[2] - partials_b1z[0]

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
    return ( 1./(4*np.pi) )* integral


@njit
def _jit_B_coulomb(ms_slice, x0, rcut, include):
    integral = np.zeros((3,),dtype=np.float32)

    for ind in range(ms_slice.data_arr.shape[0]):
        measure = ms_slice.data_arr[ind, ms_slice.varidx['measure']]
        x       = ms_slice.data_arr[ind, ms_slice.varidx['x']]
        y       = ms_slice.data_arr[ind, ms_slice.varidx['y']]
        z       = ms_slice.data_arr[ind, ms_slice.varidx['z']]

        distanceSquared = x**2 + y**2 + z**2
        if distanceSquared < rcut**2:
            continue

        if include is not None and not include[ind]: 
            continue

        partials_b1x = ms_slice.get_native_partial_derivatives(ind, 'b1x')
        partials_b1y = ms_slice.get_native_partial_derivatives(ind, 'b1y')
        partials_b1z = ms_slice.get_native_partial_derivatives(ind, 'b1z')

        div_B1 = partials_b1x[0] + partials_b1y[1] + partials_b1z[2]

        r_x = x0[0] - x
        r_y = x0[1] - y
        r_z = x0[2] - z
        r = np.sqrt(r_x**2 + r_y**2 + r_z**2)
        if r < 1e-5:
            continue

        integrand = np.empty((3,),dtype=np.float32)
        integrand[0] = div_B1*r_x
        integrand[1] = div_B1*r_y
        integrand[2] = div_B1*r_z
        integrand[:] =  integrand[:]/r**3

        integral[:] = integral[:] + measure*integrand
        #ret[i_eps,:] = ret[i_eps,:] + integrand

    #for i_eps in range(n_eps):
    #    ret[i_eps,:] = ( (unique_epsilons[i_eps]**3)/(4*np.pi) ) * ret[i_eps,:]
    return ( 1./(4*np.pi) )* integral


def slice_bs_msph(run, time, ms_slice, obs_point, insubset=None, subsetStr=None):
    '''
    insubset is None (default) or a python function which acts on (N,3) float arrays returning (N,) boolean arrays
    '''
    funcnameStr = 'bs_msph'

    x0 = util.GetMagnetometerCoordinates(obs_point, time, 'GSM', 'car')

    if insubset is None:
        include = None
        subsetStr = ''
    else:
        include = insubset(ms_slice.data_arr[:, [ms_slice.varidx['x'],ms_slice.varidx['y'],ms_slice.varidx['z']]])
        if subsetStr is None:
            subsetStr = ''

    integral = _jit_B_biotsavart(ms_slice, x0, run['rCurrents'], include)
    integral = hx.GSMtoSM(integral, time, ctype_in='car', ctype_out='car')
    x0 = hx.GSMtoSM(x0, time, ctype_in='car', ctype_out='car')
    integral = hx.get_NED_vector_components(integral.reshape(1,3), x0.reshape(1,3)).ravel()

    outname = f'{run["rundir"]}/derived/timeseries/slices/' \
        + f'{funcnameStr}{subsetStr}-{obs_point}-{util.Tstr(time)}.npy'
    np.save(outname, integral)

def slice_cl_msph(run, time, ms_slice, obs_point):
    funcnameStr = 'cl_msph'

    x0 = util.GetMagnetometerCoordinates(obs_point, time, 'GSM', 'car')

    #includeStr = '_x_gt_0'
    #include = ms_slice.data_arr[:, ms_slice.varidx['x']] > 0.
    includeStr = ''
    include = None

    integral = _jit_B_coulomb(ms_slice, x0, run['rCurrents'], include)
    integral = hx.GSMtoSM(integral, time, ctype_in='car', ctype_out='car')
    x0 = hx.GSMtoSM(x0, time, ctype_in='car', ctype_out='car')
    integral = hx.get_NED_vector_components(integral.reshape(1,3), x0.reshape(1,3)).ravel()

    outname = f'{run["rundir"]}/derived/timeseries/slices/' \
        + f'{funcnameStr}{includeStr}-{obs_point}-{util.Tstr(time)}.npy'
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


def stitch_cl_msph(run, times, obs_point):
    funcnameStr = 'cl_msph'

    integrals = []
    for time in times:
        outname = f'{run["rundir"]}/derived/timeseries/slices/' \
            + f'{funcnameStr}-{obs_point}-{util.Tstr(time)}.npy'

        integrals.append(np.load(outname))

    arr_name = f'{run["rundir"]}/derived/timeseries/' \
            + f'{funcnameStr}-{obs_point}.npy'
    arr = np.array(integrals)
    np.save(arr_name, arr)

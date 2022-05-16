import logging
import numpy as np
from numba import njit

from hxform import hxform as hx
from magnetopost import util

@njit
def get_integrand(x, normal, x0, ms_slice):
    B = np.empty(3)
    B[0] = ms_slice.interpolate(x, 'b1x')
    B[1] = ms_slice.interpolate(x, 'b1y')
    B[2] = ms_slice.interpolate(x, 'b1z')

    r = x0 - x
    normr =  np.linalg.norm(r)
    return (1./normr**3) * (r*np.dot(B,normal) + np.cross(r, np.cross(B,normal)) )

@njit
def _jit_helm_outer(ms_slice, x0):
    xmin = ms_slice.xGlobalMin
    ymin = ms_slice.yGlobalMin
    zmin = ms_slice.zGlobalMin
    xmax = ms_slice.xGlobalMax
    ymax = ms_slice.yGlobalMax
    zmax = ms_slice.zGlobalMax

    if xmax-xmin != 256 or ymax-ymin != 256 or zmax-zmin != 256:
        print('WARNING: unexpected bounds')#!!!

    N = 512 #!!! hard coded
    ax = 0.5 + np.arange(N)

    dx = (xmax-xmin)/N 
    dy = (ymax-ymin)/N 
    dz = (zmax-zmin)/N 

    x_ax = dx*ax + xmin
    y_ax = dy*ax + ymin
    z_ax = dz*ax + zmin

    integral = np.zeros(3)
    for i in range(N):
        for j in range(N):

            point  = np.array([x_ax[i], y_ax[j], zmin])
            normal = np.array([0.     , 0.     , -1. ])
            integral[:] += get_integrand(point, normal, x0, ms_slice)*dx*dy

            point  = np.array([x_ax[i], y_ax[j], zmax])
            normal = np.array([0.     , 0.     , +1. ])
            integral[:] += get_integrand(point, normal, x0, ms_slice)*dx*dy

            point  = np.array([x_ax[i], ymin, z_ax[j]])
            normal = np.array([0.     , -1. , 0.     ])
            integral[:] += get_integrand(point, normal, x0, ms_slice)*dx*dz

            point  = np.array([x_ax[i], ymax, z_ax[j]])
            normal = np.array([0.     , +1. , 0.     ])
            integral[:] += get_integrand(point, normal, x0, ms_slice)*dx*dz

            point  = np.array([xmin, y_ax[i], z_ax[j]])
            normal = np.array([-1. , 0.     , 0.     ])
            integral[:] += get_integrand(point, normal, x0, ms_slice)*dz*dy

            point  = np.array([xmax, y_ax[i], z_ax[j]])
            normal = np.array([+1. , 0.     , 0.     ])
            integral[:] += get_integrand(point, normal, x0, ms_slice)*dz*dy

    integral[:] = (-1./(4.*np.pi)) * integral
    return integral


def helm_outer(info, time, ms_slice, obs_point):

    funcnameStr = 'helm_outer'

    x0 = util.GetMagnetometerCoordinates(obs_point, time, 'GSM', 'car')

    integral = _jit_helm_outer(ms_slice, x0)
    integral = hx.GSMtoSM(integral, time, ctype_in='car', ctype_out='car')
    x0 = hx.GSMtoSM(x0, time, ctype_in='car', ctype_out='car')
    integral = hx.get_NED_vector_components(integral.reshape(1,3), x0.reshape(1,3)).ravel()

    outname = f'{info["dir_run"]}/derived/timeseries/timesteps/' \
        + f'{funcnameStr}-{obs_point}-{util.Tstr(time)}.npy'
    np.save(outname, integral)
    logging.info(f"Writing {outname}")


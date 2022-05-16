import numpy as np
import logging
from numba import njit

from hxform import hxform as hx
from magnetopost import util
from magnetopost.units_and_constants import phys
from magnetopost.gap_integrals import get_dipole_field

@njit()
def get_dipole_field_V(xyz):
    # Xyz_D and returned b_D in SMG (SM) coordinates
    b = np.empty(xyz.shape, dtype=xyz.dtype)
    for i in range(xyz.shape[0]):
        b[i,:] = get_dipole_field(xyz[i,:])
    return b

@njit # Biot Savart integrand (4pi taken outside)
def _integrand_bs(x, y, z, x0, jx, jy, jz):
    integrand = np.empty((3,),dtype=np.float32); integrand[:]=np.nan

    r_x = x0[0] - x
    r_y = x0[1] - y
    r_z = x0[2] - z
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
def _integral_bs(X, Y, Z, x0, JX, JY, JZ, Measure):
    ret = np.zeros((3,), dtype=np.float32)
    for i in range(X.size):

        # 'measure' refers to integration measure ( dx, dV, dA, whatever apropriate)
        ret[:] = ret[:] + _integrand_bs(X[i], Y[i], Z[i], x0, JX[i], JY[i], JZ[i]) * Measure[i]

    return ret/(4.*np.pi)

# in ModIonoMagPerturb.f90
#          call get_planet_field(Time_simulation, XyzIono_D, 'SMG', bIono_D)
#          bIono_D = bIono_D/sqrt(sum(bIono_D**2))
#          ! get the Hall and Perdersen currents in xyz coords
#          Jh_IID(i,j,:) = cross_product(bIono_D, eIono_IID(i,j,:))*SigmaH(i,j)
#          Jp_IID(i,j,:) = eIono_IID(i,j,:) * SigmaP(i,j)

def bs_pedersen(info, time, ie_slice, obs_point):
    funcnameStr = 'bs_pedersen'

    x0 = util.GetMagnetometerCoordinates(obs_point, time, 'SM', 'car')

    data_arr, varidx, units = ie_slice

    X  = data_arr[varidx['X'], :]
    Y  = data_arr[varidx['Y'], :]
    Z  = data_arr[varidx['Z'], :]

    # note, surface current density (not vol current density)
    KX = data_arr[varidx['Ex'],:] * data_arr[varidx['SigmaP'],:]
    KY = data_arr[varidx['Ey'],:] * data_arr[varidx['SigmaP'],:]
    KZ = data_arr[varidx['Ez'],:] * data_arr[varidx['SigmaP'],:]

    Measure = data_arr[varidx['measure'],:]

    assert( units['SigmaP'] == 'S' and units['Ex'] == 'mV/m' and units['X'] == 'R' )
    scalefact = phys['mu0']*phys['Siemens']*phys['mV']/phys['m']

    integral = scalefact*_integral_bs(X,Y,Z,x0,KX,KY,KZ,Measure)
    integral = hx.get_NED_vector_components(integral.reshape(1,3), x0.reshape(1,3)).ravel()

    outname = f'{info["dir_run"]}/derived/timeseries/timesteps/' \
        + f'{funcnameStr}-{obs_point}-{util.Tstr(time)}.npy'
    np.save(outname, integral)
    logging.info(f"Wrote {outname}")


def bs_hall(info, time, ie_slice, obs_point):
    funcnameStr = 'bs_hall'

    x0 = util.GetMagnetometerCoordinates(obs_point, time, 'SM', 'car')

    data_arr, varidx, units = ie_slice

    XYZ = data_arr[[varidx['X'],varidx['Y'],varidx['Z']], :].transpose()
    E_iono = data_arr[[varidx['Ex'],varidx['Ey'],varidx['Ez']], :].transpose()
    # note, surface current density (not vol current density)
    unit_b_dipole = get_dipole_field_V(XYZ)
    unit_b_dipole = unit_b_dipole/np.linalg.norm(unit_b_dipole, axis=1)[:,None]
    K = np.cross(unit_b_dipole, E_iono) * data_arr[varidx['SigmaH'],:][:,None]
    Measure = data_arr[varidx['measure'],:]

    assert( units['SigmaH'] == 'S' and units['Ex'] == 'mV/m' and units['X'] == 'R' )
    scalefact = phys['mu0']*phys['Siemens']*phys['mV']/phys['m']

    integral = scalefact*_integral_bs(XYZ[:,0],XYZ[:,1],XYZ[:,2],x0,K[:,0],K[:,1],K[:,2],Measure)
    integral = hx.get_NED_vector_components(integral.reshape(1,3), x0.reshape(1,3)).ravel()

    outname = f'{info["dir_run"]}/derived/timeseries/timesteps/' \
        + f'{funcnameStr}-{obs_point}-{util.Tstr(time)}.npy'
    np.save(outname, integral)
    logging.info(f"Wrote {outname}")


def slice_integral_bs_bulkiono(ie_slice, obs_point):
    # Not used
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

    integral = scalefact*_integral_bs(X,Y,Z,x0,Jx,Jy,Jz,Measure)
    print(integral)


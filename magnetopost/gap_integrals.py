import numpy as np
from numba import njit
import datetime
import pandas as pd

from hxform import hxform as hx
from magnetopost import util
from magnetopost.units_and_constants import phys

@njit
def matvectprod(A, v):
    ret = np.zeros(A.shape[0])
    for i in range(A.shape[0]):
        for j in range(A.shape[0]):
            ret[i] += A[i,j]*v[j]
    return ret

@njit
def sph_to_xyz(r, theta, phi):
    Xyz_D = np.empty(3)
    Xyz_D[0] = r*np.cos(phi)*np.sin(theta)
    Xyz_D[1] = r*np.sin(phi)*np.sin(theta)
    Xyz_D[2] = r*np.cos(theta)
    return Xyz_D

@njit
def map_along_dipole_lines(Xyz_D, rMap):
    # Xyz_D and returned XyzMap_D in SMG (SM) coordinates
    # Solution of the vector potential equation (proportional to (x^2+y^2)/r^3)
    # so sqrt(xMap^2+yMap^2)/sqrt(x^2+y^2) = sqrt(rMap^3/r^3)
    iHemisphere = int(np.sign(Xyz_D[2]))
    XyzMap_D = np.empty(3, dtype=np.float32)

    r = np.linalg.norm(Xyz_D)
    XyRatio = np.sqrt(rMap/r)**3 # ratio of input and mapped X and Y components
    XyzMap_D[0:2] = XyRatio*Xyz_D[0:2]
    XyMap2 = XyzMap_D[0]**2 + XyzMap_D[1]**2

    if rMap**2 < XyMap2:
       # The point does not map to the given radius
       iHemisphere = 0
       # Put mapped point to the magnetic equator
       XyzMap_D[0:2] = (rMap/np.sqrt(Xyz_D[0]**2 + Xyz_D[1]**2))*Xyz_D[0:2]
       XyzMap_D[2] = 0
    else:
       XyzMap_D[2] = iHemisphere*np.sqrt(rMap**2 - XyMap2)

    return XyzMap_D, iHemisphere

@njit
def get_dipole_field(xyz):
    # Xyz_D and returned b_D in SMG (SM) coordinates
    b = np.empty(3)
    r = np.linalg.norm(xyz)
    DipoleStrength = 3.12e+4 #"dipole moment"(not really) in  nT * R_e**3  # https://en.wikipedia.org/wiki/Dipole_model_of_the_Earth%27s_magnetic_field
    Term1      = DipoleStrength*xyz[2]*3/r**2
    b[0:2] = Term1*xyz[0:2]/r**3
    b[2]    = (Term1*xyz[2]-DipoleStrength)/r**3
    return b


@njit
def _jit_mhd_SurfaceIntegral(ms_slice, GM_2_gap, x0, nTheta,nPhi, rCurrents):
    # x0 and returned dB_mhd_SurfaceIntegral are in cartesian gap_csys coordinates (default SMG)
    dB_mhd_SurfaceIntegral = np.zeros(3)

    dTheta = np.pi    / (nTheta-1)
    dPhi   = 2.*np.pi / nPhi

    for iTheta in range(nTheta):
        Theta = iTheta * dTheta
        # the area of the triangle formed by pole and the latitude segment at Theta=dTheta/2
        # is approximately dTheta/4*dTheta/2, so sin(theta) replaced with dTheta/8.
        SinTheta = max(np.sin(Theta), dTheta/8.)
        dSurface = rCurrents**2*SinTheta*dTheta*dPhi

        for iPhi in range(nPhi):
            Phi = iPhi * dPhi

            xyz_Currents = sph_to_xyz(rCurrents, Theta, Phi)

            #b1_Currents = np.empty(3,dtype='f8')
            b1_Currents = np.empty(3)
            xyz_inGM = matvectprod(GM_2_gap.transpose(), xyz_Currents)
            # use GM interpolator, which is in GM_csys coordinates, to get b1 and convert to gap_csys coordinates
            b1_Currents[0] = ms_slice.interpolate(xyz_inGM, 'b1x')
            b1_Currents[1] = ms_slice.interpolate(xyz_inGM, 'b1y')
            b1_Currents[2] = ms_slice.interpolate(xyz_inGM, 'b1z')
            b1_Currents[:] = matvectprod(GM_2_gap, b1_Currents)#?? doesnt seem to work without [:] for some reason

            Unit_xyz_Currents = xyz_Currents / rCurrents
            InvDist2_D = dSurface*(xyz_Currents - x0)/(4*np.pi*np.sqrt(np.sum((xyz_Currents - x0)**2))**3)

            B1r   = np.dot(Unit_xyz_Currents, b1_Currents)
            B1t = np.cross(Unit_xyz_Currents, b1_Currents)
            dB_mhd_SurfaceIntegral[:] = dB_mhd_SurfaceIntegral[:] + B1r*InvDist2_D + np.cross(B1t, InvDist2_D)

            #Br == n.b
            #Bt == n&b
            #InvDist2_D == (dA/4pi) * (x - x0)/(|x - x0|**3)
            #dB == (dA/4pi)*[ Br*((x - x0)/(|x - x0|**3)) + Bt&((x - x0)/(|x - x0|**3)) ]
            #   == (dA/4pi)*[ (n.b)*(x - x0) + (n&b)&(x - x0) ]/(|x - x0|**3)

    return dB_mhd_SurfaceIntegral


@njit
def _jit_fac_integral(ms_slice, GM_2_gap, x0, nTheta,nPhi,nR, rCurrents):
    # x0 and returned dB_fac are in cartesian gap_csys coordinates (default SMG)
    rIonosphere = 1.01725 # rEarth + iono_height #!!! hard coded

    dB_fac                 = np.zeros(3)

    dTheta = np.pi    / (nTheta-1)
    dPhi   = 2.*np.pi / nPhi
    dR     = (rCurrents - rIonosphere) / nR

    for iTheta in range(nTheta):
        Theta = iTheta * dTheta
        # the area of the triangle formed by pole and the latitude segment at Theta=dTheta/2
        # is approximately dTheta/4*dTheta/2, so sin(theta) replaced with dTheta/8.
        SinTheta = max(np.sin(Theta), dTheta/8.)
        dSurface = rCurrents**2*SinTheta*dTheta*dPhi

        for iPhi in range(nPhi):
            Phi = iPhi * dPhi

            xyz_Currents = sph_to_xyz(rCurrents, Theta, Phi)

            #b_Currents = np.empty(3,dtype='f8');
            b_Currents = np.empty(3);
            #j_Currents = np.empty(3,dtype='f8')
            j_Currents = np.empty(3)
            xyz_inGM = matvectprod(GM_2_gap.transpose(), xyz_Currents)
            # use GM interpolator, which is in GM_csys coordinates, to get b and j and convert to gap_csys coordinates
            b_Currents[0] = ms_slice.interpolate(xyz_inGM, 'bx')
            b_Currents[1] = ms_slice.interpolate(xyz_inGM, 'by')
            b_Currents[2] = ms_slice.interpolate(xyz_inGM, 'bz')
            j_Currents[0] = ms_slice.interpolate(xyz_inGM, 'jx')
            j_Currents[1] = ms_slice.interpolate(xyz_inGM, 'jy')
            j_Currents[2] = ms_slice.interpolate(xyz_inGM, 'jz')
            b_Currents[:] = matvectprod(GM_2_gap, b_Currents)#?? doesnt seem to work without [:] for some reason
            j_Currents[:] = matvectprod(GM_2_gap, j_Currents)

            Unit_xyz_Currents = xyz_Currents / rCurrents
            Unit_b_Currents = b_Currents / np.linalg.norm(b_Currents)
            _Fac_rCurrents = np.dot(Unit_b_Currents,j_Currents)#!!!!!!!!!!!
            Fac_term = _Fac_rCurrents * np.dot(Unit_b_Currents, Unit_xyz_Currents)

            for k in range(nR):
                R = rCurrents - dR*(k+0.5)

                xyz_Map, iHemisphere = map_along_dipole_lines(xyz_Currents, R)
                b0_Map = get_dipole_field(xyz_Map)

                Unit_xyz_Map = xyz_Map / R
                Unit_b0_Map = b0_Map / np.linalg.norm(b0_Map)

                # The volume element is proportional to 1/Br. The sign
                # should be preserved (not yet!!!),
                # because the sign is also there in the radial
                # component of the field aligned current: Br/B*FAC.
                # In the end j_D = b_D/Br*[(Br/B)*(j.B)]_rcurr  !!!!!!!!!!!! NOT DIMENSIONALLY CONSISTENT

                dVol_FACcoords = dSurface * (dR/np.dot(Unit_xyz_Map, Unit_b0_Map))
                J_fac = Fac_term * Unit_b0_Map

                dB_fac[:] = dB_fac[:] + dVol_FACcoords* \
                  np.cross(J_fac, x0-xyz_Map)/(4*np.pi*(np.linalg.norm(xyz_Map-x0))**3)

    return dB_fac


def slice_bs_fac(run, time, ms_slice, obs_point, nTheta=181,nPhi=180,nR=30, gap_csys='SM'):
    funcnameStr = 'bs_fac'

    assert(gap_csys=='SM')
    x0 = util.GetMagnetometerCoordinates(obs_point, time, 'SM', 'car') #!!!! gap_csys?

    GM_csys = 'GSM'

    GM_2_gap = hx.get_transform_matrix(time, GM_csys, gap_csys)

    dB_fac = _jit_fac_integral(ms_slice, GM_2_gap, x0, nTheta,nPhi,nR, run['rCurrents'])
    dB_fac = (phys['mu0']*phys['muA']/phys['m']**2) * dB_fac
    dB_fac = hx.get_NED_vector_components(dB_fac.reshape(1,3), x0.reshape(1,3)).ravel()

    outname = f'{run["rundir"]}/derived/timeseries/slices/' \
        + f'{funcnameStr}-{obs_point}-{util.Tstr(time)}.npy'

    np.save(outname, dB_fac)


def slice_helm_rCurrents(run, time, ms_slice, obs_point, nTheta=181,nPhi=180, gap_csys='SM'):
    funcnameStr = 'helm_rCurrents'

    x0 = util.GetMagnetometerCoordinates(obs_point, time, gap_csys, 'car')

    GM_csys = 'GSM'

    GM_2_gap = hx.get_transform_matrix(time, GM_csys, gap_csys)

    dB_mhd_SurfaceIntegral = _jit_mhd_SurfaceIntegral(ms_slice, GM_2_gap, x0, nTheta,nPhi, run['rCurrents'])
    dB_mhd_SurfaceIntegral = hx.transform(dB_mhd_SurfaceIntegral, time, gap_csys, 'SM')
    x0 = hx.transform(x0, time, gap_csys, 'SM')
    dB_mhd_SurfaceIntegral = hx.get_NED_vector_components(dB_mhd_SurfaceIntegral.reshape(1,3), x0.reshape(1,3)).ravel()

    outname = f'{run["rundir"]}/derived/timeseries/slices/' \
        + f'{funcnameStr}_gap{gap_csys}-{obs_point}-{util.Tstr(time)}.npy'

    np.save(outname, dB_mhd_SurfaceIntegral)


def stitch_bs_fac(run, times, obs_point):
    funcnameStr = 'bs_fac'

    integrals = []
    for time in times:
        outname = f'{run["rundir"]}/derived/timeseries/slices/' \
            + f'{funcnameStr}-{obs_point}-{util.Tstr(time)}.npy'

        integrals.append(np.load(outname))

    arr_name = f'{run["rundir"]}/derived/timeseries/' \
            + f'{funcnameStr}-{obs_point}.npy'
    arr = np.array(integrals)
    np.save(arr_name, arr)


def stitch_helm_rCurrents(run, times, obs_point, gap_csys='SM'):
    funcnameStr = 'helm_rCurrents'

    integrals = []
    for time in times:
        outname = f'{run["rundir"]}/derived/timeseries/slices/' \
            + f'{funcnameStr}_gap{gap_csys}-{obs_point}-{util.Tstr(time)}.npy'

        integrals.append(np.load(outname))

    arr_name = f'{run["rundir"]}/derived/timeseries/' \
            + f'{funcnameStr}_gap{gap_csys}-{obs_point}.npy'
    arr = np.array(integrals)
    np.save(arr_name, arr)


if __name__ == '__main__':
    from magnetopost.model_patches import SWMF
    sl = SWMF.get_ms_slice_class('/home/gary/temp/3d__var_3_e20031120-070000-000.out')
    slice_bs_fac(None, (2019,9,2,4,11,0), sl, "origin")

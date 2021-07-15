import numpy as np
from numba import njit
import datetime
import pandas as pd

from hxform import hxform as hx
from magnetopost import util
from magnetopost.units_and_constants import phys
import magnetopost.tmp_cxtransform as cx

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
def _jit_gap_region_integrals(ms_slice, GM_2_gap, x0, nTheta,nPhi,nR, rCurrents):
    # x0 and returned dB_fac, dB_mhd_SurfaceIntegral are in cartesian gap_csys coordinates (default SMG)
    rIonosphere = 1.01725 # rEarth + iono_height #!!! hard coded

    dB_fac                 = np.zeros(3)
    dB_mhd_SurfaceIntegral = np.zeros(3)

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
            b0_Currents = get_dipole_field(xyz_Currents)

            b_Currents = np.empty(3,dtype='f8'); j_Currents = np.empty(3,dtype='f8')
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

            #####################
            Br   = np.dot(Unit_xyz_Currents, b_Currents)
            Bt = np.cross(Unit_xyz_Currents, b_Currents)
            InvDist2_D = dSurface*(xyz_Currents - x0)/(4*np.pi*np.sqrt(np.sum((xyz_Currents - x0)**2))**3)
            dB_mhd_SurfaceIntegral[:] = dB_mhd_SurfaceIntegral[:] + Br*InvDist2_D + np.cross(Bt, InvDist2_D)
            #####################

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

    return dB_fac, dB_mhd_SurfaceIntegral

def slice_bs_fac(run, time, ms_slice, obs_point, nTheta=181,nPhi=180,nR=30, gap_csys='SM'):
    funcnameStr = 'bs_fac'

    if obs_point == "origin":
        x0 = np.zeros(3)
    else:
        x0 = util.GetMagnetometerCoordinates(obs_point, time, 'SM', 'car')

    GM_csys = 'GSM'
    assert(gap_csys=='SM')

    GM_2_gap = np.empty((3,3))
    _ = cx.transform([1.,0.,0.], (2019,9,2,6,30,0), 'GSM', 'MAG')#!!!!!!!! DOESNT WORK UNLESS UNCOMMENTED
    GM_2_gap[:, 0] = cx.transform([1.,0.,0.], time, GM_csys, gap_csys)
    GM_2_gap[:, 1] = cx.transform([0.,1.,0.], time, GM_csys, gap_csys)
    GM_2_gap[:, 2] = cx.transform([0.,0.,1.], time, GM_csys, gap_csys)

    dB_fac, dB_mhd_SurfaceIntegral = _jit_gap_region_integrals(ms_slice, GM_2_gap, x0, nTheta,nPhi,nR, run['rCurrents'])
    dB_fac = (phys['mu0']*phys['muA']/phys['m']**2) * dB_fac
    dB_fac = hx.get_NED_vector_components(dB_fac.reshape(1,3), x0.reshape(1,3)).ravel()
    dB_mhd_SurfaceIntegral = hx.get_NED_vector_components(dB_mhd_SurfaceIntegral.reshape(1,3), x0.reshape(1,3)).ravel()

    outname = f'{run["rundir"]}/derived/timeseries/slices/' \
        + f'{funcnameStr}-{obs_point}-{util.Tstr(time)}.npy'

    outname_SURF = f'{run["rundir"]}/derived/timeseries/slices/' \
        + f'helm_rCurrent-{obs_point}-{util.Tstr(time)}.npy'

    np.save(outname, dB_fac)
    np.save(outname_SURF, dB_mhd_SurfaceIntegral)


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

    integrals_SURF = []
    for time in times:
        outname_SURF = f'{run["rundir"]}/derived/timeseries/slices/' \
            + f'helm_rCurrent-{obs_point}-{util.Tstr(time)}.npy'
        integrals_SURF.append(np.load(outname_SURF))

    arr_name_SURF = f'{run["rundir"]}/derived/timeseries/' \
            + f'helm_rCurrent-{obs_point}.npy'
    arr_SURF = np.array(integrals_SURF)
    np.save(arr_name_SURF, arr_SURF)

if __name__ == '__main__':
    from magnetopost.model_patches import SWMF
    sl = SWMF.get_ms_slice_class('/home/gary/temp/3d__var_3_e20031120-070000-000.out')
    slice_bs_fac(None, (2019,9,2,4,11,0), sl, "origin")

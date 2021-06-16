import sys
import os
import numpy as np

import spacepy.coordinates as sc
from spacepy.time import Ticktock

def tpad(time, length=6):
    assert(length==6)
    return tuple(time)


def transform(v, time, csys_in, csys_out, ctype_in=None, ctype_out=None):
    """Transfrom between coordinates systems using SpacePy library.

    Parameters
    ----------
    v : array-like
        list of three floats
        list containing lists of three floats
        np.array of three floats
        (Nv, 3) float np.array, where Nv=1 or Nv=Nt
        
    time : array-like
        list of 3+ ints
        list containing lists of 3+ ints
        np.array of 3+ ints
        (Nt, 3) float np.array, where Nt = 1 or Nt = Nv

        The 3+ ints are [year, month, day, [hours, [minutes, [seconds]]]]

    csys_in : str

    csys_out : str
    
    ctype_in : str
               'car' (default) or 'sph'

    ctype_out : str
               'car' (default) or 'sph'
               
    Returns
    -------
    array-like with structure matching either `time` (if `Nt` != 1) or
    `v` (if `Nv` =! 1). If `Nv` and `Nt` != 1, structure matches `v`. 

    Examples
    --------
    t1 = [2000, 1, 1, 0, 0, 0]
    t2 = [2000, 1, 1, 2, 0, 0]
    v1 = [0, 0, 1]
    v2 = [0, 1, 0]

    # All equivalent and return a list with three floats    
    transform([0, 1, 1], time1)
    transform([0, 1, 1], time1, ctype_in='car')
    transform([0, 1, 1], time1, ctype_in='car', ctype_out='car')

    # The following 3 calls return a list with two lists of 3 elements
    # 1. Transform two vectors at same time t1
    transform([v1, v2], t1) 

    # 2. Transform two vectors at different times
    transform([v1, v2], [t1, t2])

    # 3. Transform one vector at different times
    transform(v1, [t1, t2])

    # The following returns a (2, 3) np.array
    # Transform one vector at three times
    transform(np.array([v1, v2]), t1) 

    """
    
    import numpy.matlib
    
    if ctype_in is None:
        ctype_in = 'car'
    if ctype_out is None:
        ctype_out = 'car'

    in_type = type(v)
    v = np.array(v)
    t = np.array(time)

    if len(t.shape) > 1 and len(v.shape) > 1:
        if t.shape[0] != v.shape[0]:
            raise ValueError("t and v cannot be different lengths")
    if len(v.shape) == 1 and len(t.shape) > 1:
        v = numpy.matlib.repmat(v, t.shape[0], 1)
    if len(t.shape) == 1 and len(v.shape) > 1:
        t = numpy.matlib.repmat(t, v.shape[0], 1)

    cvals = sc.Coords(v, csys_in, ctype_in)
    if len(t.shape)==1:
        t_str = '%04d-%02d-%02dT%02d:%02d:%02d' % tpad(t, length=6)
    else:
        t_str = []
        for i in range(t.shape[0]):
            t_str.append('%04d-%02d-%02dT%02d:%02d:%02d' % tpad(t[i,:], length=6))
        t_str = np.array(t_str)
            

    cvals.ticks = Ticktock(t_str, 'ISO')
    newcoord = cvals.convert(csys_out, ctype_out)

    ret = newcoord.data
    if len(t.shape) == 1 and len(v.shape) == 1:
        ret = newcoord.data[0, :]

    if in_type == np.ndarray:
        return ret
    else:
        return ret.tolist()


def MAGtoGSM(v_MAG, time, ctype_in, ctype_out):
    return transform(v_MAG, time, 'MAG', 'GSM', ctype_in=ctype_in, ctype_out=ctype_out)

def GSMtoMAG(v_GSM, time, ctype_in, ctype_out):
    return transform(v_GSM, time, 'GSM', 'MAG', ctype_in=ctype_in, ctype_out=ctype_out)

def GEOtoGSM(v_GEO, time, ctype_in, ctype_out):
    return transform(v_GEO, time, 'GEO', 'GSM', ctype_in=ctype_in, ctype_out=ctype_out)

def GEOtoMAG(v_GEO, time, ctype_in, ctype_out):
    """
    note: no interpolation is done on the dipole position between years,
        so this function is actually only dependent on time[0]
    """
    return transform(v_GEO, time, 'GEO', 'MAG', ctype_in=ctype_in, ctype_out=ctype_out)

def MAGtoGEO(v_MAG, time, ctype_in, ctype_out):
    """
    note: no interpolation is done on the dipole position between years,
        so this function is actually only dependent on time[0]
    """
    return transform(v_MAG, time, 'MAG', 'GEO', ctype_in=ctype_in, ctype_out=ctype_out)

def MAGtoSM(v_MAG, time, ctype_in, ctype_out):
    return transform(v_MAG, time, 'MAG', 'SM', ctype_in=ctype_in, ctype_out=ctype_out)

def MAGtoGEI(v_MAG, time, ctype_in, ctype_out):
    return transform(v_MAG, time, 'MAG', 'GEI', ctype_in=ctype_in, ctype_out=ctype_out)

def GEOtoGEI(v_GEO, time, ctype_in, ctype_out):
    return transform(v_GEO, time, 'GEO', 'GEI', ctype_in=ctype_in, ctype_out=ctype_out)

def GEOtoSM(v_GEO, time, ctype_in, ctype_out, fixedYear=False):
    if fixedYear:
        time = [2000] + list(time)[1:]
    return transform(v_GEO, time, 'GEO', 'SM', ctype_in=ctype_in, ctype_out=ctype_out)

def dipole(time):
    return MAGtoGEO(np.array([0.,0.,1.]), time, 'car', 'sph')

def StoC(r, theta, phi):
    """Convert from spherical to cartesian coordinates

    r, theta, phi: array-like

    """
    x = r*np.cos(phi)*np.sin(theta)
    y = r*np.sin(phi)*np.sin(theta)
    z = r*np.cos(theta)
    return x, y, z


def UTtoHMS(UT, **kwargs):
    """Convert universal time in fractional hours into integer hour, minutes, seconds
    
    from cxtransform import UTtoHMS

    print(UTtoHMS(12))              # [12, 0, 0]

    print(UTtoHMS(24))              # [0, 0, 0]

    print(UTtoHMS(24, keep24=True)) # [24, 0, 0]
    """

    keep24 = False
    if 'keep24' in kwargs:
        keep24 = kwargs['keep24']
        
    if UT > 24 or UT < 0:
        raise ValueError('Required: 0 <= UT <= 24.')

    hours = int(UT)
    minutes = int((UT-hours)*60.)
    seconds = int(round((UT-hours-minutes/60.)*3600.))
    if seconds == 60:
        seconds = 0
        minutes = minutes + 1
    if minutes == 60:
        minutes = 0
        hours = hours + 1

    if hours == 24 and keep24 == False:
        return [0, 0, 0]

    return [hours, minutes, seconds]


def MAGtoMLT(pos, time, csys='sph', debug=False):
    """Compute magnetic local time given a UT and MAG position or longitude

    Uses equation 93 in https://arxiv.org/abs/1611.10321

    Usage: 
    ------
    import cxtransform as cx
    mlt = cx.MAGtoMLT(MAGlong, time)
    mlt = cx.MAGtoMLT([MAGlong1, Mlong2, ...], time)

    mlt = cx.MAGtoMLT([MAGx, MAGy, MAGz], time, csys='car')
    mlt = cx.MAGtoMLT([[MAGx1, MAGy1, MAGz1],...], time, csys='car')

    Returns:
    --------
    mlt: float or array-like

    Examples:
    --------
    import cxtransform as cx

    mlt = cx.MAGtoMLT(0., [2000, 1, 1, 0, 0, 0])
    print(mlt) # 18.869936573301775

    mlt = cx.MAGtoMLT([0., 0.], [2000, 1, 1, 0, 0, 0])
    print(mlt) # [18.86993657 18.86993657]

    mlt = cx.MAGtoMLT([-1., 0., 0.], [2000, 1, 1, 0, 0, 0], csys='car')
    print(mlt) # 6.869936573301775

    mlt = cx.MAGtoMLT([[-1., 0., 0.],[-1., 0., 0.]], [2000, 1, 1, 0, 0, 0], csys='car')
    print(mlt) # [6.86993657 6.86993657]
    """

    assert(csys == 'car' or csys == 'sph')
    pos = np.array(pos)
    time = np.array(time)
    if not isinstance(pos, float):
        pos = np.array(pos)

    if csys == 'sph':
        phi = pos*np.pi/180.
    else:
        if pos.shape == (3, ):
            phi = np.arctan2(pos[1], pos[0])
        else:
            phi = np.arctan2(pos[:, 1], pos[:, 0])

    if debug:
        print('phi =' + str(phi))

    subsol_pt = transform(np.array([1, 0, 0]), time, 'GSM', 'MAG')

    #import pdb;pdb.set_trace()
    if len(subsol_pt.shape) == 1:
        phi_cds = np.arctan2(subsol_pt[1], subsol_pt[0])
    else:
        phi_cds = np.arctan2(subsol_pt[:, 1], subsol_pt[:, 0])

    if debug:
        print('phi_cds =' + str(phi_cds))

    delta = phi - phi_cds # note np.array([a1, a2, ...])+b == np.array([a1+b, a2+b, ...])

    if debug:
        print('delta =' + str(delta))

    if isinstance(delta,float):
        delta = np.array([delta])

    idx = np.where(delta > np.pi)
    delta[idx] = delta[idx] - 2.*np.pi
    idx = np.where(delta <= -np.pi)
    delta[idx] = delta[idx] + 2.*np.pi

    if delta.size == 1:
        delta = delta[0]

    MLT = 12. + delta*24./(2.*np.pi)
    return MLT



def GSMtoMAGLocalComponents(time, mlat, mlon, dB):
    if dB.shape == (3,):
        return toMAGLocalComponents(time, mlat, mlon, np.array([dB]))[0,:]

    time = np.array(time, dtype=int)

    if len(time.shape) == 1:
        time = np.array([time])

    N = time.shape[0]
    M = dB.shape[0]
    assert(dB.shape == (M,3))

    station_pos = MAGtoGSM(np.array([1., mlat, mlon]), time, 'sph', 'car')

    Pole = MAGtoGSM(np.array([0., 0., 1.]), time, 'car', 'car')

    U3 = station_pos
    U3norm = np.sqrt(U3[:,0]**2 + U3[:,1]**2 + U3[:,2]**2) #( not really needed since norm is 1 in these units where R_e=1, but just to be consistent in all units)
    divU3norm = 1./U3norm
    U3 = U3*divU3norm[:,np.newaxis]
    U1 = np.cross(Pole, U3)
    U1norm = np.sqrt(U1[:,0]**2 + U1[:,1]**2 + U1[:,2]**2)
    divU1norm = 1./U1norm
    #https://stackoverflow.com/questions/5795700/multiply-numpy-array-of-scalars-by-array-of-vectors
    U1 = U1*divU1norm[:,np.newaxis]
    U2 = np.cross(U3, U1)
    assert(U1.shape == (N, 3)) #check

    #print(U1)
    #print(U2)
    #print(U3)

    R = np.empty((N, 3, 3))
    R[:, :, 0] = U2
    R[:, :, 1] = U1
    R[:, :, 2] = -U3
    
    R = np.linalg.inv(R)
    assert(R.shape == (N,3,3))

    if N == 1:
        R = np.repeat(R, dB.shape[0], axis=0)
    elif N == M:
        pass
    else:
        raise ValueError('dimensions of time and dB dont match')

    
    #for i in range(N):
        #R[i,:,0] = U2[i,:]
        #R[i,:,1] = U1[i,:]
        #R[i,:,2] = -U3[i,:]
    #    M = np.column_stack([U2[i,:], U1[i,:], -U3[i,:]])
    #    R[i, :, :] = np.linalg.inv(M)

    dB_rot = np.einsum('ijk,ik->ij', R, dB)
    return dB_rot

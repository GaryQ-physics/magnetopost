import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../')
from config import conf

from scipy.integrate import odeint
import cxtransform as cx
from scipy.interpolate import RegularGridInterpolator

import util
from probe import probe
from fieldlines import fieldlines
import biot_savart_kameleon_interpolated_grid as bsk


xlims = [-100., 15.]
ylims = [-10., 10.] #!!!!!!!!!!!!!!!!
zlims = [-15., 15.]
dx = 0.3
dy = 0.3
dz = 0.3


#def fetch(filename, X, var, mlat=0., mlon=0.):

#    if not 'dB' in var:

#        return probe(filename, X, var=var, library='kameleon')
#    else:
#        bsk.integrate(run, time_fname, mlat, mlon, para=True,
#        xlims=(-48., 16.), ylims=(-32., 32.), zlims=(-32., 32.), d=0.125, returnAll=False)



def data_in_U(run, time, variable, u, v, U, mlat=0., mlon=0., Vchar=0.125**3):
    """Data in U coordinates"""

    U1 = U[0]
    U2 = U[1]
    U3 = U[2]
    if isinstance(u, float):
        X = u*U1 + v*U2
    else:
        X = np.einsum('i,j->ij', u, U1) + np.einsum('i,j->ij', v, U2)
    
    if 'bu' in variable:
        filename = util.time2CDFfilename(run, time)
        B = fetch(filename, X, ['bx', 'by', 'bz'])
        if variable == 'bu1':
            return np.dot(B, U1)
        if variable == 'bu2':
            return np.dot(B, U2)
        if variable == 'bu3':
            return np.dot(B, U3)
    elif 'dB' in variable:
        print('X.shape =' + str(X.shape))
        ret = bsk.dB_dV_slice(run, time, mlat, mlon, u, v, U, returnAll=True)
        if variable == 'dB_Magnitude':
            print('dB.shape =' + str(ret[0].shape))
            return np.sqrt(np.einsum('ij,ij->i',Vchar*ret[0], Vchar*ret[0]))
        else:
            dB_loc = bsk.toMAGLocalComponents(time, mlat, mlon, Vchar*ret[0])
            if variable == 'dB_north':
                return dB_loc[:, 0]
            if variable == 'dB_east':
                return dB_loc[:, 1]
            if variable == 'dB_down':
                return dB_loc[:, 2]
    else:
# Since _CCMC used in file filemeta, in order to run plotting, cannot
# pass library='kameleonV' since then both kameleon and kameleonV
# library='kameleon' in the same shell.
        ret = probe(util.time2CDFfilename(run,time), X, var=variable, library='kameleon')
        #print('HELLO THERE\n\n\n\n\n\n\n\n\n\n')
        #print(np.min(ret))
        #print(np.max(ret))
        return ret


def data2d(run, time, parameter, X, Y, U, debug=False, mlat=0., mlon=0.):

    # grid of the corresponding values of variable. To be color plotted
    Z = data_in_U(run, time, parameter,
                        X.flatten(), Y.flatten(), U, mlat=mlat, mlon=mlon)

    #print('HELLO again\n\n\n\n\n\n\n\n\n\n')
    #print(np.min(Z))
    #print(np.max(Z))

    return Z.reshape(X.shape)


def unitvector(run, time, mag, debug=False):

    # Points where field line is returned on are points a distance
    # 0, 0.5, 1.0 along the field line from the starting point.
    s_grid = np.array([0., 0.5, 1.])
    # Compute (x, y, z) of points at s_grid values.
    sols = fieldlines(run, time, mag, s_grid=s_grid, debug=debug, max_iterations=1)
    ret = []
    for sol in sols:
        if debug: print(sol.shape)
        # initialize vectors for defining field line cut plane
        v1 = (np.nan)*np.empty((3, ))
        v2 = (np.nan)*np.empty((3, ))
        v3 = (np.nan)*np.empty((3, ))
        U1 = (np.nan)*np.empty((3, ))
        U2 = (np.nan)*np.empty((3, ))
        U3 = (np.nan)*np.empty((3, ))

        # Three vectors in from origin to point on field line.
        v1 = sol[0, :]
        v2 = sol[2, :]
        v3 = sol[1, :]

        # Define cut plane coordinates based on field line 
        # (U3 is normal to the plane)
        U2 = (v1 - v2)/np.linalg.norm(v1-v2)
        U3 = np.cross(v3 - v1, U2)

        if np.linalg.norm(U3) < 1e-3:
            print("WARNING: close to straight line")
        U3 = U3/np.linalg.norm(U3)
        U1 = np.cross(U2, U3)
        ret.append([U1, U2, U3])
    return ret


def writedata(time, mlat, mlon, debug=False):
    """Write output of unitvector() to file
    
    Calling unitvector() from ParaView does not work, so write output to
    txt file.
    """

    # Compute centered dipole unit vector in GSM at given time
    Mdipole = cx.MAGtoGSM([0., 0., 1.], time[0:6], 'car', 'car')
    U1, U2, U3 = unitvector(time, np.array([1., mlat, mlon]))[0]

    subdir = '%04d%02d%02dT%02d%02d/' % tuple(time[0:5])
    if not os.path.exists(conf["run_path_derived"] + subdir):
        os.mkdir(conf["run_path_derived"] + subdir)

    out_fname = conf["run_path_derived"] + subdir + \
                'cut_plane_info_%.2f_%.2f' % (mlat, mlon) + util.maketag(time) + '.txt'
    f = open(out_fname, 'w')
    
    print('Writing ' + out_fname)
    f.write('%.7e %.7e %.7e\n' % (Mdipole[0], Mdipole[1], Mdipole[2]))
    f.write('%.7e %.7e %.7e\n' % (U1[0], U1[1], U1[2]))
    f.write('%.7e %.7e %.7e\n' % (U2[0], U2[1], U2[2]))
    f.write('%.7e %.7e %.7e\n' % (U3[0], U3[1], U3[2]))
    f.close()

    if debug:
        print('Wrote ' + out_fname)
        print(time)
        print('Mdipole = ', Mdipole)
        print('U1 = ', U1)
        print('U2 = ', U2)
        print('U3 = ', U3)

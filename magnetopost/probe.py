import logging
import numpy as np
from magnetopost import util


def slice_probe(info, time, ms_slice, obs_point):
    funcnameStr = 'probe'

    x0 = util.GetMagnetometerCoordinates(obs_point, time, 'GSM', 'car')

    b1x = ms_slice.interpolate(x0, 'b1x')
    b1y = ms_slice.interpolate(x0, 'b1y')
    b1z = ms_slice.interpolate(x0, 'b1z')
    bx = ms_slice.interpolate(x0, 'bx')
    by = ms_slice.interpolate(x0, 'by')
    bz = ms_slice.interpolate(x0, 'bz')
    jx = ms_slice.interpolate(x0, 'jx')
    jy = ms_slice.interpolate(x0, 'jy')
    jz = ms_slice.interpolate(x0, 'jz')
    ux = ms_slice.interpolate(x0, 'ux')
    uy = ms_slice.interpolate(x0, 'uy')
    uz = ms_slice.interpolate(x0, 'uz')
    p = ms_slice.interpolate(x0, 'p')
    rho = ms_slice.interpolate(x0, 'rho')

    arr = np.array([b1x, b1y, b1z, bx, by, bz, jx, jy, jz, ux, uy, uz, p, rho])

    outname = f'{info["dir_derived"]}/timeseries/timesteps/' \
        + f'{funcnameStr}-{obs_point}-{util.Tstr(time)}.npy'
    np.save(outname, arr)
    logging.info(f"Wrote {outname}")


def stitch_probe(info, times, obs_point):
    funcnameStr = 'probe'

    arrs = []
    for time in times:
        outname = f'{info["dir_derived"]}/timeseries/timesteps/' \
            + f'{funcnameStr}-{obs_point}-{util.Tstr(time)}.npy'

        arrs.append(np.load(outname))

    outname = f'{info["dir_derived"]}/timeseries/' \
            + f'{funcnameStr}-{obs_point}.npy'
    arrs = np.array(arrs)
    np.save(outname, arrs)
    logging.info(f"Wrote {outname}")

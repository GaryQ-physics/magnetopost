import os
import ast
import logging
import numpy as np

from magnetopost.config import defined_magnetometers
from hxform import hxform as hx

def Tstr(time, length=6):
    return '%.4d%.2d%.2dT%.2d%.2d%.2d'%(time[:6])


def setup(dir_run, info):
    assert os.path.exists(dir_run), "dir_run = " + dir_run + " not found"
    
    dir_derived = os.path.join(dir_run, "derived")
    dir_slices = os.path.join(dir_derived, "timeseries", "slices")
    dir_figures = os.path.join(dir_derived, "figures")
    
    if not os.path.exists(dir_derived):
        os.mkdir(os.path.join(dir_derived))
        logging.info("Created " + dir_derived)
    
    if not os.path.exists(dir_slices):
        os.makedirs(os.path.join(dir_slices))
        print("Created " + dir_slices)
    
    if not os.path.exists(dir_figures):
        os.makedirs(os.path.join(dir_figures))
        logging.info("Created " + dir_figures)
    
    from magnetopost.model_patches import SWMF
    SWMF.generate_filelist_txts(dir_run, info)


def prep_run(dir_run):

    with open(os.path.join(dir_run, 'derived', 'run.info.py'), 'r') as f:
        run = ast.literal_eval(f.read())

    run['rundir'] = dir_run

    run['magnetosphere_files'] = {}
    with open(os.path.join(run['rundir'], 'derived', 'magnetosphere_files.txt'), 'r') as f:
        for line in f.readlines():
            items = line.split(' ')
            time = tuple([int(ti) for ti in items[:6]])
            run['magnetosphere_files'][time] = os.path.join(run['rundir'], items[-1][:-1])

    run['ionosphere_files'] = {}
    with open(os.path.join(run['rundir'], 'derived', 'ionosphere_files.txt'), 'r') as f:
        for line in f.readlines():
            items = line.split(' ')
            time = tuple([int(ti) for ti in items[:6]])
            run['ionosphere_files'][time] = os.path.join(run['rundir'], items[-1][:-1])

    return run


def GetMagnetometerCoordinates(magnetometer, time, csys, ctype, hxway=True):
    if isinstance(magnetometer, str):
        magnetometer = defined_magnetometers[magnetometer]

    if hxway:
        return hx.transform(np.array(magnetometer.coords), time, magnetometer.csys, csys, ctype_in=magnetometer.ctype, ctype_out=ctype)
    else: 
        import spacepy.coordinates as sc
        from spacepy.time import Ticktock
        cvals = sc.Coords(magnetometer.coords, magnetometer.csys, magnetometer.ctype)
        t_str = '%04d-%02d-%02dT%02d:%02d:%02d'%(time[:6])
        cvals.ticks = Ticktock(t_str, 'ISO')
        _ = cvals.convert('MAG', 'car')
        newcoord = cvals.convert(csys, ctype)
        return newcoord.data[0, :]

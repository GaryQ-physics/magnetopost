import os
import ast
from magnetopost.config import defined_magnetometers

def Tstr(time, length=6):
    return '%.4d%.2d%.2dT%.2d%.2d%.2d'%(time[:6])


def prep_run(rundir):
    if not rundir[-1]=='/': rundir = f'{rundir}/'

    with open(f'{rundir}derived/run.info.py', 'r') as f:
        run = ast.literal_eval(f.read())

    run['rundir'] = rundir

    run['magnetosphere_files'] = {}
    with open(rundir+'derived/magnetosphere_files.txt', 'r') as f:
        for line in f.readlines():
            items = line.split(' ')
            time = tuple([int(ti) for ti in items[:6]])
            filename = f'{rundir}{items[-1][:-1]}'
            run['magnetosphere_files'][time] = filename

    run['ionosphere_files'] = {}
    with open(rundir+'derived/ionosphere_files.txt', 'r') as f:
        for line in f.readlines():
            items = line.split(' ')
            time = tuple([int(ti) for ti in items[:6]])
            filename = f'{rundir}{items[-1][:-1]}'
            run['ionosphere_files'][time] = filename

    return run


import spacepy.coordinates as sc
from spacepy.time import Ticktock
def GetMagnetometerCoordinates(magnetometer, time, csys, ctype):
    if isinstance(magnetometer, str):
        magnetometer = defined_magnetometers[magnetometer]

    cvals = sc.Coords(magnetometer.coords, magnetometer.csys, magnetometer.ctype)
    t_str = '%04d-%02d-%02dT%02d:%02d:%02d'%(time[:6])

    cvals.ticks = Ticktock(t_str, 'ISO')
    newcoord = cvals.convert(csys, ctype)

    return newcoord.data[0, :]

import os
import ast
from magnetosphere.config import defined_magnetometers

def Tstr(time, length=6):
	return '%.4d%.2d%.2dT%.2d%.2d%.2d'%(time[:6])


def prep_run(confpath):
    if not confpath[-1]=='/': confpath = f'{confpath}/'

    with open(f'{confpath}derived/run.info.py', 'r') as f:
        run = ast.literal_eval(f.read())

    run['magnetosphere_files'] = {}
    with open(confpath+'derived/magnetosphere_files.txt', 'r') as f:
        for line in f.readlines():
            items = line.split(' ')
            time = tuple(items[:6])
            filename = f'{confpath}{items[-1][:-1]}'
            run['magnetosphere_files'][time] = filename

    run['ionosphere_files'] = {}
    with open(confpath+'derived/ionosphere_files.txt', 'r') as f:
        for line in f.readlines():
            items = line.split(' ')
            time = tuple(items[:6])
            filename = f'{confpath}{items[-1][:-1]}'
            run['ionosphere_files'][time] = filename

    return run


#def get_available_times(run):
    #times = [md.filename2time(run, fname) for fname in get_available_ms_files()]
    #return run['magnetosphere_files'].keys()

def GetMagnetometerCoordinates(magnetometer, time, csys, ctype):
	if isinstance(magnetometer, str):
		magnetometer = defined_magnetometers[magnetometer]

	magnetometer.csys 
	assert(False)

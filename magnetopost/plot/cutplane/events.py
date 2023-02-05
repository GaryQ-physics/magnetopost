import os
import sys
import numpy as np

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

import cxtransform as cx
from util import time2filename

def eventlist():
    """List of events
    
    Returns
    -------
    np.ndarray of [year, month, day, hour, minute, mag_lat, mag_lon, mlt]
    """

    infile = conf["run_path_derived"] + 'LOCALIZED.txt'

    data = np.genfromtxt(infile, skip_header=1)    

    mlt = cx.MAGtoMLT(data[:, 5], data[:, 0:5])

    # Swap mlat and mlon colums so in expected order (lat then long)
    data[:, [6,5]] = data[:, [5,6]]
    
    data = np.hstack((data, np.reshape(mlt, (mlt.shape[0], 1))))
    
    return data

def events(time):
    """Returns events that ocurred at a given time
    
    Example
    -------
    
    events([2003,11,20])
    events([2003, 11, 20, 17, 46])
    """

    event_list = eventlist()
    idx = np.all(time == event_list[:, 0:len(time)], axis=1)
    return event_list[idx,:]

def count(time):
    """Returns # of events at a given time"""
    
    return len(events(time))

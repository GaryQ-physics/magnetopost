def _logger(level=None):
    import sys
    import logging
    logger = logging.getLogger('magnetopost')

    for handler in logger.handlers[:]:
        logger.removeHandler(handler)

    handler = logging.StreamHandler(stream=sys.stdout)
    formatter = logging.Formatter('%(asctime)s.%(msecs)03d:%(filename)s:%(funcName)s(): %(message)s', datefmt='%H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if level is None:
        logger.setLevel(logging.INFO)

    return logger

logger = _logger()

from magnetopost import util
from magnetopost.postproc import job_ms
from magnetopost.postproc import job_ie
from magnetopost.plot.timeseries import plot
from magnetopost.extract_magnetometer_data import extract_from_swmf_ccmc_printout_file
from magnetopost.extract_magnetometer_data import extract_from_swmf_magnetometer_files
from magnetopost.extract_magnetometer_data import extract_from_magnetopost_files


import os
import resource

if os.path.exists('/home/gary/'):
    #https://stackoverflow.com/questions/16779497/how-to-set-memory-limit-for-thread-or-process-in-python
    soft, hard = int(13*2**30), int(13*2**30)
    resource.setrlimit(resource.RLIMIT_AS,(soft, hard))
elif os.path.exists('/home/gquaresi/'):
    soft, hard = 90*2**30, 90*2**30
    resource.setrlimit(resource.RLIMIT_AS,(soft, hard))

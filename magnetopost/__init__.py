if False:
    def logger():
        import sys
        import logging
        logger = logging.getLogger('magnetospost')
        handler = logging.StreamHandler(stream=sys.stdout)
        formatter = logging.Formatter('%(asctime)s.%(msecs)03d:%(filename)s:%(funcName)s(): %(message)s', datefmt='%H:%M:%S')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    
        for handler in logger.handlers[:]:
            logger.removeHandler(handler)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)
    
        return logger
    
    logger = logger()

import logging
logging.basicConfig(
    format='%(filename)s:%(funcName)s(): %(message)s',
    level=logging.INFO,
    datefmt='%S')

from magnetopost import util
from magnetopost.postproc import job_ms
from magnetopost.postproc import job_ie
from magnetopost.plot.timeseries import plot
from magnetopost.extract_magnetometer_data import extract_from_swmf_ccmc_printout_file
from magnetopost.extract_magnetometer_data import extract_from_swmf_magnetometer_files
from magnetopost.extract_magnetometer_data import extract_from_magnetopost_files


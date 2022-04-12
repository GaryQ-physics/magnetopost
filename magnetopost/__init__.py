import logging
logging.basicConfig(
    format='%(filename)s:%(funcName)s(): %(message)s',
    level=logging.INFO,
    datefmt='%S')

from magnetopost import util
from magnetopost.postproc import job_ms
from magnetopost.postproc import job_ie
from magnetopost.plot import msph_point
from magnetopost.plot import surf_point
from magnetopost.extract_magnetometer_data import extract_from_swmf_ccmc_printout_file
from magnetopost.extract_magnetometer_data import extract_from_swmf_magnetometer_files
from magnetopost.extract_magnetometer_data import extract_from_magnetopost_files

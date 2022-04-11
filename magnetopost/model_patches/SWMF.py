import logging

from swmf_file_reader import batsrus_class as bats
from swmf_file_reader.read_ie_files import read_iono_tec


def get_ms_slice_class(filename):
    assert(filename[-4:] == '.out')
    return bats.get_class_from_native(filename[:-4])


def get_iono_slice(filename):
    return read_iono_tec(filename)


def generate_filelist_txts(dir_run, info):

    import os
    import re
    import json

    fn = os.path.join(dir_run, 'derived/run.info.py')
    with open(fn, 'w') as outfile:
        outfile.write(json.dumps(info))

    logging.info("Wrote {}".format(fn))

    magnetosphere_outs = sorted(os.listdir(os.path.join(dir_run, 'GM/IO2')))
    fn = os.path.join(dir_run, 'derived/magnetosphere_files.txt')
    k = 0
    with open(fn,'w') as fl:
        regex = r"3d__var_.*\.out$"
        for fname in magnetosphere_outs:
            if re.search(regex, fname):
                k = k + 1
                assert(fname[:8] == '3d__var_')
                assert(fname[9] == '_')
                if fname[10] == 'e':
                    Y = int(fname[11:15])
                    M = int(fname[15:17])
                    D = int(fname[17:19])
                    assert(fname[19] == '-')
                    h = int(fname[20:22])
                    m = int(fname[22:24])
                    s = int(fname[24:26])
                    assert(fname[26] == '-')
                    mil = int(fname[27:30])
                    assert(fname[30:] == '.out')
                    fl.write(f'{Y} {M} {D} {h} {m} {s} {mil} GM/IO2/{fname}\n')

    logging.info("Wrote {} file names to {}".format(k, fn))

    ionosphere_outs = sorted(os.listdir(os.path.join(dir_run, 'IE/ionosphere')))
    fn = os.path.join(dir_run, 'derived/ionosphere_files.txt')
    k = 0
    with open(fn,'w') as fl:
        regex = r"i_.*\.tec$"

        for fname in ionosphere_outs:
            if re.search(regex, fname):
                k = k + 1
                assert(fname[:2] == 'i_')
                if fname[2] == 'e':
                    Y = int(fname[3:7])
                    M = int(fname[7:9])
                    D = int(fname[9:11])
                    assert(fname[11] == '-')
                    h = int(fname[12:14])
                    m = int(fname[14:16])
                    s = int(fname[16:18])
                    assert(fname[18] == '-')
                    mil = int(fname[19:22])
                    assert(fname[22:] == '.tec')
                    fl.write(f'{Y} {M} {D} {h} {m} {s} {mil} IE/ionosphere/{fname}\n')

    logging.info("Wrote {} file names to {}".format(k, fn))

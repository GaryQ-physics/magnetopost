import subprocess
from swmf_file_reader import read_swmf_files as rswmf
#from swmf_file_reader.read_ie_files import read_iono_tec
from magnetopost.model_patches.tmp_read_ie_files import read_iono_tec

def get_ms_slice_class(filename):
    assert(filename[-4:] == '.out')
    return rswmf.return_class(filename[:-4])


def get_iono_slice(filename):
    return read_iono_tec(filename)


##########################################################
#########################################################

def generate_filelist_txts():
    import os
    import re

    magnetosphere_outs = sorted(os.listdir('./GM/IO2'))
    with open('/tmp/magnetosphere_files.txt','w') as fl:
        regex = r"3d__var_.*\.out$"

        for fname in magnetosphere_outs:
            if re.search(regex, fname):
                assert(fname[:8] == '3d__var_')
                PlotNumber = int(fname[8]) # hopefully this is only in the single didgets or this will break
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
                    fl.write(f'{Y} {M} {D} {h} {m} {s} {mil} ./GM/IO2/{fname}\n')


    ionosphere_outs = sorted(os.listdir('./IE/ionosphere'))
    with open('/tmp/ionosphere_files.txt','w') as fl:
        regex = r"i_.*\.tec$"

        for fname in ionosphere_outs:
            if re.search(regex, fname):
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
                    fl.write(f'{Y} {M} {D} {h} {m} {s} {mil} ./IE/ionosphere/{fname}\n')


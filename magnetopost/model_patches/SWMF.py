import subprocess
from swmf_file_reader import read_swmf_files as rswmf
#from swmf_file_reader.read_ie_files import read_iono_tec
from magnetopost.model_patches.tmp_read_ie_files import read_iono_tec

def get_ms_slice_class(filename):
    assert(filename[-4:] == '.out')
    return rswmf.return_class(filename[:-4])


def get_iono_slice(filename):
    return read_iono_tec(filename)


def generate_filelist_txts(): # TODO
    assert(fname[:11] == '3d__var_3_e')
    Y = int(fname[11:15])
    M = int(fname[15:17])
    D = int(fname[17:19])
    assert(fname[19] == '-')
    h = int(fname[20:22])
    m = int(fname[22:24])
    s = int(fname[24:26])
    assert(fname[26] == '-')
    f = int(fname[27:30])
    assert(fname[30:] == '.out')
    return (Y,M,D,h,m,s,f)

#def get_magnetosphere_filetag(run, time):
#	if len(time) == 7:
#		return '%s3d__var_2_e%.4d%.2d%.2d-.2d%.2d%.2d-%.3d.out'%(run["magnetosphere_files"], *time)
#	elif len(time) == 6:
#		fname = '%s3d__var_*_e%.4d%.2d%.2d-.2d%.2d%.2d-*.out'%(run["magnetosphere_files"], *time)
#		stdout = subprocess.check_output(f'ls -1 {fname}', shell=True, text=True)
#		return stdout.split('\n')[0]
#
#def get_ionosphere_filename(run, time):
#    if len(time) == 7:
#        return '%si_e%.4d%.2d%.2d-.2d%.2d%.2d-%.3d.tec'%(run["ionosphere_files"], *time)
#    elif len(time) == 6:
#        fname = '%si_e%.4d%.2d%.2d-.2d%.2d%.2d-*.tec'%(run["ionosphere_files"], *time)
#        stdout = subprocess.check_output(f'ls -1 {fname}', shell=True, text=True)
#        return stdout.split('\n')[0]
#
#def get_available_ms_files(run)
#	stdout = subprocess.check_output(f'ls -1 {run["magnetosphere_files"]}3d__var_3_e*.out', shell=True, text=True)
#	return stdout.split('\n')

from swmf_file_reader import batsrus_class as bats
from swmf_file_reader.read_ie_files import read_iono_cdf

def get_ms_slice_class(filename):
    assert(filename[-4:] == '.cdf')
    return bats.get_class_from_cdf(filename)

def get_iono_slice(filename):
    return read_iono_cdf(filename)

def generate_filelist_txts(info):
    import subprocess
    script = '''
#!/bin/sh
awk -F '[ /:]' 'NR>1 {print $5" "$6" "$7" "$10" "$11" "$12" 0 ./GM_CDF/"$1}' GM_CDF/SWPC_SWMF_052811_2_GM_cdf_list > derived/magnetosphere_files.txt
awk -F '[ /:]' 'NR>1 {print $5" "$6" "$7" "$10" "$11" "$12" 0 ./IONO-2D_CDF/SWPC_SWMF_052811_2.swmf."$1}' IONO-2D_CDF/SWPC_SWMF_052811_2_IE_CDF_list > derived/ionosphere_files.txt
'''
    return subprocess.call(script, shell=True)

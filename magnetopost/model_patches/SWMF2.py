from swmf_file_reader import batsrus_class as bats
from swmf_file_reader.read_ie_files import read_iono_cdf

def get_ms_slice_class(filename):
    assert(filename[-4:] == '.cdf')
    return bats.return_class(filename)

def get_iono_slice(filename):
    return read_iono_cdf(filename)





##########################################################k
# #!/bin/bash
# cd GM_CDF/
# awk -F '[ /:]' '{print $5" "$6" "$7" "$10" "$11" "$12" 0 ./GM_CDF/"$1}' SWPC_SWMF_052811_2_GM_cdf_list > /tmp/magnetosphere_files.txt
# cd ../IONO-2D_CDF/
# awk -F '[ /:]' '{print $5" "$6" "$7" "$10" "$11" "$12" 0 ./IONO-2D_CDF/SWPC_SWMF_052811_2.swmf."$1}' SWPC_SWMF_052811_2_IE_CDF_list > /tmp/ionosphere_files.txt
# 
# then remove first line manually

def generate_filelist_txts():
    import os
    import re
    with open('/media/sunspot/git-data/sblake/SWPC_SWMF_052811_2/GM_CDF/SWPC_SWMF_052811_2_GM_cdf_list','r') as f:
        pass
    with open('/tmp/magnetosphere_files.txt','w') as fl:
        pass

try:
    from magnetopost.model_patches import SWMF
except:
    pass
try:
    from magnetopost.model_patches import SWMF2
except:
    pass
try:
    from magnetopost.model_patches import LFM
except:
    pass


def get_ms_slice_class(run, time):
    fname = run['magnetosphere_files'][time]

    if   run['model'] == 'SWMF':
        return SWMF.get_ms_slice_class(fname)
    elif run['model'] == 'SWMF2':
        return SWMF2.get_ms_slice_class(fname)
    elif run['model'] == 'LFM':
        return LFM.get_ms_slice_class(fname)

def get_iono_slice(run, time):
    fname = run['ionosphere_files'][time] 

    if   run['model'] == 'SWMF':
        return SWMF.get_iono_slice(fname)
    elif run['model'] == 'SWMF2':
        return SWMF2.get_iono_slice(fname)
    elif run['model'] == 'LFM':
        return LFM.get_iono_slice(fname)

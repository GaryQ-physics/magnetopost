def get_ms_slice_class(info, run, time):

    fname = run['magnetosphere_files'][time]

    if info['model'] == 'SWMF':
        if info['file_type'] == "out":
            from magnetopost.model_patches import SWMF_out
            return SWMF_out.get_ms_slice_class(fname)
        else:
            from magnetopost.model_patches import SWMF_cdf
            return SWMF_cdf.get_ms_slice_class(fname)


def get_iono_slice(info, run, time):
    fname = run['ionosphere_files'][time] 

    if info['model'] == 'SWMF':
        if info['file_type'] == "out":
            from magnetopost.model_patches import SWMF_out
            return SWMF_out.get_iono_slice(fname)
        else:
            from magnetopost.model_patches import SWMF_cdf
            return SWMF_cdf.get_iono_slice(fname)

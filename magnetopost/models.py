def get_ms_slice_class(info, time):

    fname = info['files']['magnetosphere'][time]

    if info['model'] == 'SWMF':
        import swmfio
        return swmfio.read_batsrus(fname) # Call with .out removed


def get_iono_slice(info, time):

    fname = info['files']['ionosphere'][time]

    if info['model'] == 'SWMF':
        import swmfio
        return swmfio.read_rim(fname)

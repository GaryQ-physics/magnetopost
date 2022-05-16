def get_ms_slice_class(info, time):

    fname = info['files']['magnetosphere'][time]

    if info['model'] == 'SWMF':
        if info['file_type'] == "out":
            from swmf_file_reader import batsrus_class as bats
            return bats.get_class_from_native(fname[:-4])
        else:
            from swmf_file_reader import batsrus_class as bats
            return bats.get_class_from_cdf(fname)


def get_iono_slice(info, time):

    fname = info['files']['ionosphere'][time]

    if info['model'] == 'SWMF':
        if info['file_type'] == "out":
            from swmf_file_reader.read_ie_files import read_iono_tec
            return read_iono_tec(fname)
        else:
            from swmf_file_reader.read_ie_files import read_iono_cdf
            return read_iono_cdf(fname)

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


#def filename2time(fname):
#	if   run['model'] == 'SWMF':
#		return SWMF.filename2time(fname)
#	elif run['model'] == 'SWMF2':
#		return SWMF2.filename2time(fname)
#	elif run['model'] == 'LFM':
#		return LFM.filename2time(fname)
#
#def get_magnetosphere_filetag(run, time):
#	if   run['model'] == 'SWMF':
#		return f"{run['magnetosphere_files']}{SWMF.get_magnetosphere_filetag(fname)}"
#	elif run['model'] == 'SWMF2':
#		return f"{run['magnetosphere_files']}{SWMF2.get_magnetosphere_filetag(fname)}"
#	elif run['model'] == 'LFM':
#		return f"{run['magnetosphere_files']}{LFM.get_magnetosphere_filetag(fname)}"
#
#def get_ionosphere_filename(run, time):
#	if   run['model'] == 'SWMF':
#		return f"{run['ionosphere_files']}{SWMF.get_ionosphere_filename(time)}"
#	elif run['model'] == 'SWMF2':
#		return f"{run['ionosphere_files']}{SWMF2.get_ionosphere_filename(time)}"
#	elif run['model'] == 'LFM':
#		return f"{run['ionosphere_files']}{LFM.get_ionosphere_filename(time)}"
#
#def get_available_ms_files(run)
#	if   run['model'] == 'SWMF':
#		return SWMF.get_available_ms_files(fname)
#	elif run['model'] == 'SWMF2':
#		return SWMF2.get_available_ms_files(fname)
#	elif run['model'] == 'LFM':
#		return LFM.get_available_ms_files(fname)

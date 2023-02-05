import os
import sys
import numpy as np
import pickle

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

from niceticks import niceticks
import cutplane_plot as cp
import util
import magnetometers as mg

def configure_dict(var):

    var['str_id'] = var['varname']
    if var['mlat'] is not None:
        var['str_id'] = var['str_id'] + '_mlat_%.3f_mlon_%.3f'%(var['mlat'],var['mlon'])

    #var['var_path'] = conf(var['run']+'_derived') + 'cutplane/' + var['str_id'] + '/'

    var['min'] = 0.
    var['max'] = 0.
    #var['fileid'] = 0



def process_one(filename, var, fileid=0, pkl='', debug=True):

    filename_path = conf[var['run']+"_derived"] + "cutplanes/" + var['str_id'] + '/'
    if not os.path.exists(filename_path):
        os.makedirs(filename_path)

    filename_png = ('{0:s}-{1:s}-type_{2:d}_delta_{3:.4f}_'
                    'xlims_{4:.5f},{5:.5f}_ylims_{6:.5f},{7:.5f}.png'
                    .format(filename, var['str_id'], 1, var['delta'],
                    var['xlims'][0], var['xlims'][1],
                    var['ylims'][0], var['ylims'][1]) )
    filename_png = filename_path + filename_png

    zticks = var['zticks']
    if os.path.exists(pkl):
        with open(pkl, 'rb') as handle:
            if debug:
                print("Reading " + pkl)
            minmax = pickle.load(handle)
            if debug:
                print("Read " + pkl)

        if zticks is None:
            var_min = np.nanmin(minmax['min'])
            var_max = np.nanmax(minmax['max'])
            if not np.isnan(var_min) and not np.isnan(var_max):
                zticks = niceticks(var_min, var_max, 8, debug=debug)

    if not var['regen'] and os.path.exists(filename_png):
        # Image already exists
        print('Image exists for {0:d}th file {1:s} for {2:s} at resoluton {3:.2f}. Will not regenerate.'\
              .format(k, filename, var, var['delta']))
        return None

    util.dlfile(conf[var['run'] + '_cdf'] + filename, debug=debug)

    info = cp.plot(var['run'], util.CDFfilename2time(var['run'], filename), var['varname'], var['plane'],
                         dx=var['delta'], dy=var['delta'],
                         xlims=var['xlims'], ylims=var['ylims'],
                         dpi=var['dpi'], showplot=var['showplot'],
                         zticks=zticks, logz=True,
                         png=True, pngfile=filename_png, debug=debug,
                         mlat_dB=var['mlat'], mlon_dB=var['mlon'])

    pkl_fid = pkl + '-fid_%d.pkl'%(fileid)
    if not os.path.exists(pkl_fid):
        minmax_single = {'min' : info['min'], 'max' : info['max']}

        if debug:
            print("Writing " + pkl_fid)
        with open(pkl_fid, 'wb') as handle:
            pickle.dump(minmax_single, handle)
        if debug:
            print("Wrote " + pkl_fid)

    print('Finished {0:d}th file {1:s} for {2:s} at resoluton {3:.5f}'\
           .format(fileid, filename, var, var['delta']))


def process(files, variables, para=False, process_type=1, debug=True):
    # type(vars) == list
    for var in variables:
        # type(var) == dict
        configure_dict(var)
        times = []

        pkl_path = conf[var['run']+"_derived"] + "cutplanes/minmax/"
        if not os.path.exists(pkl_path):
            os.makedirs(pkl_path)
        pkl = pkl_path + var['str_id'] + '.pkl'

        if process_type == 1:
            var['delta'] = var['delta1']
            var['dpi'] = var['dpi1']
            var['xlims'] = var['xlimsO']
            var['ylims'] = var['ylimsO']

        if process_type == 2:
            var['delta'] = var['delta2']
            var['dpi'] = var['dpi2']
            var['xlims'] = var['xlimsO']
            var['ylims'] = var['ylimsO']

            minmax = {'min': np.nan*np.empty(len(files)),
                      'max': np.nan*np.empty(len(files)),
                      'probe': {
                                  'ux': np.nan*np.empty(len(files)),
                                  'bz': np.nan*np.empty(len(files)),
                                  'n_events': np.nan*np.empty(len(files))
                            }
                      }
            for k in range(len(files)):
                pkl_fid = pkl + '-fid_%d.pkl'%(k)
                if os.path.exists(pkl_fid):
                    with open(pkl_fid, 'rb') as handle_fid:

                        if debug:
                            print("Reading " + pkl_fid)
                        minmax_single = pickle.load(handle_fid)
                        if debug:
                            print("Read " + pkl_fid)

                    minmax['min'][k] = minmax_single['min']
                    minmax['max'][k] = minmax_single['max']

            with open(pkl, 'wb') as handle:
                pickle.dump(minmax, handle)

        if process_type == 3:
            var['delta'] = var['delta3']
            var['dpi'] = var['dpi2']
            var['xlims'] = var['xlimsI']
            var['ylims'] = var['ylimsI']

        def wrap(fname, fileid):
            process_one(fname, var, fileid=fileid, pkl=pkl)

        if para:
            from joblib import Parallel, delayed
            import multiprocessing
            num_cores = multiprocessing.cpu_count()
            if num_cores is not None and num_cores > len(files):
                num_cores = len(files)
            print('Parallel processing {0:d} file(s) using {1:d} cores'\
                  .format(len(files), num_cores))
            Parallel(n_jobs=num_cores)(\
                    delayed(wrap)(files[i], i) for i in range(len(files)))

        else:
            for i in range(len(files)):
                wrap(files[i], i)


def main(run):

    debug         = True  # Show extra logging information
    para          = False  # Process in parallel using all CPUs
    showplot      = False  # Show plot on screen. Set to false for long runs.
                           # Does not always work in Spyder/IPython, especially when
                           # para=True. Starting a new console sometimes fixes.
                           # TODO: Read PNG and display using PIL instead.

    # Testing options
    #low_only    = False  # Do only low-res first processing
    #high_only   = False  # Execute only high-res second processing
    process_types = [1, 2, 3]

    test_serial   = False   # Process few files in serial
    test_parallel = True  # Process few files in parallel

    opts = {
            'run' : run,
            'regen': True,       # Regenerate image even if found
            'plane': 'xz',       # Cut plane to plot
            'xlimsO': [-30, 15], # zoomed out Horizontal axis limits in R_E
            'ylimsO': [-20, 20], # zoomed out Vertical axis limits in R_E
            'xlimsI': [-7.53125, 7.53125], # zoomed in Horizontal axis limits in R_E
            'ylimsI': [-7.53125, 7.53125], # zoomed in Vertical axis limits in R_E
            'zticks': None,      # z-axis ticks for variable
            'dpi1': 96,          # Make multiple of 16 (for mimwrite animation)
            'dpi2': 96*3,        # Make multiple of 16 (for mimwrite animation)
            'delta1': 0.5,       # Low-res cut plane resolution in R_E
            'delta2': 0.125,     # High-res cut plane resolution in R_E
            'delta3': 0.0625,    # Zoomed cut plane resolution in R_E
            'showplot': showplot,# Show the plot on screen
            'mlat' : None,
            'mlon' : None
            }
    nf = None

    built_in_vars = ['bx','by','bz','ux','uy','uz','jx','jy','jz','rho','p','e']
    dB_vars = ['dB_Magnitude', 'dB_north', 'dB_east', 'dB_down']
    MAG_locations = [(0.,0.)] # list of (mlat, mlot) tuples

    ###
    if run in ['CARR_IMPULSE', 'DIPTSUR2', 'TESTANALYTIC', 'IMP10_RUN_SAMPLE']:
        fixed_time = (2019,9,2,6,30,0) # only the year 2019 matters
        pos = mg.GetMagnetometerLocation('colaba', fixed_time, 'MAG', 'sph')
        MAG_locations.append((pos[1], pos[2]))

    if test_serial:
        para = False
        #built_in_vars = ['jy']
        built_in_vars = ['bx','by','bz','ux','uy','uz','jx','jz','rho','e']
        nf = 4
        opts["showplot"] = False

    if test_parallel:
        para = True
        #built_in_vars = ['bx','by','bz','ux','uy','uz','jx','jz','rho','e']
        #dB_vars = []
        nf = None
        opts["showplot"] = False
        process_types = [3]
    ###

    files = list(util.get_available_slices(run)[0])
    if nf is not None:
        files = files[:nf]


    variables = []
    for variable in built_in_vars:
        vardict = opts.copy()
        vardict['varname'] = variable
        vardict['zticks'] = None
        variables.append(vardict.copy())

    for variable in dB_vars:
        for loc in MAG_locations:
            if len(loc) != 2:
                raise ValueError('MAG locations have two components, mlat, mlon (radius is earth surface)')
            vardict = opts.copy()
            vardict['varname'] = variable
            vardict['zticks'] = None
            vardict['mlat'] = loc[0]
            vardict['mlon'] = loc[1]
            variables.append(vardict.copy())

    for process_type in process_types:
        print('%d processing'%(process_type))
        process(files, variables, para=para, process_type=process_type)


if __name__ == '__main__':
    #util.generate_TESTANALYTIC_cdflist()
    main('DIPTSUR2')


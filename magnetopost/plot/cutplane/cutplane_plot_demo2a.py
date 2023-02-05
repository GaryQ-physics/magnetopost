"""
Generate cut plane for each file and for each variable associate with BATSRUS
run in conf['run_path'] specified in config.py. Files are downloaded as-needed.

Use cutplane_plot_demo2_animate.py to create .mp4 from images.

First set of cut planes at low resolution is used to compute min/max of
each variable. Second set is at high resolution with fixed axes limits
with limits determined using low resolution min/max values.

For long runs, execute from command line - Spyder often freezes after 50+ files
are processed.

When running in parallel wih Spyder/IPython, sometimes there is a broken pipe
error message. Program always runs after second attempt, however. In addition,
the disabiling of output buffering will not work (see line that starts with
sys.stdout = ...). 

TODO: If number of cores > number of vars, parallelize by files instead of
variables?

TODO: Generalize to use OpenGGCM and LFM. Main modification needed is the
generalization of CDF file names.
"""

import os
import sys
import numpy as np
import pickle

sys.path.append(os.path.dirname(os.path.abspath(__file__)) + '/../' )
from config import conf

from niceticks import niceticks
import cutplane_plot as cp
import util

run = "SWPC"
#run = "SCARR5"

debug         = True  # Show extra logging information
para          = False  # Process in parallel using all CPUs
showplot      = False  # Show plot on screen. Set to false for long runs.
                       # Does not always work in Spyder/IPython, especially when
                       # para=True. Starting a new console sometimes fixes.
                       # TODO: Read PNG and display using PIL instead.
plot_type     = 2

# Testing options
first_only    = False  # Do only low-res first processing
second_only   = False  # Execute only high-res second processing
test_serial   = True   # Process few files in serial
test_parallel = False  # Process few files in parallel

variables = ['bx','by','bz','ux','uy','uz','jx','jy','jz','rho','p','e']



opts = {
        'regen': True,       # Regenerate image even if found
        'plane': 'xz',       # Cut plane to plot
        'xlims': [-30, 15],  # Horizontal axis limits in R_E
        'ylims': [-20, 20],  # Vertical axis limits in R_E
        'zticks': None,      # z-axis ticks for variable
        'nf': None,          # Number of files to process. None => all files
        'dpi1': 96,          # Make multiple of 16 (for mimwrite animation)
        'dpi2': 96*3,        # Make multiple of 16 (for mimwrite animation)
        'delta1': 0.5,       # Low-res cut plane resolution in R_E
        'delta2': 0.125,     # High-res cut plane resolution in R_E
        'showplot': showplot # Show the plot on screen
        }

if first_only and second_only:
    raise ValueError("Both first_only and second_only are True. Only one may be True.")

#(u, idx, ne) = np.unique(event_list[:, 0:5], axis=0, return_index=True, return_counts=True)
#event_list = event_list[idx, :]
#event_list = np.hstack((event_list, np.reshape(ne, (ne.size, 1))))

if test_serial:
    para = False

    #variables = ['p', 'dB_Magnitude']
    variables = ['jx','jy','jz','p', 'dB_Magnitude', 'dB_north', 'dB_east', 'dB_down']

    opts["nf"] = 2
    opts["showplot"] = True

if test_parallel:
    para = True
    variables = ['p', 'jy']
    opts["nf"] = 3
    opts["showplot"] = True
    #opts["nf"] = None
    #opts["showplot"] = False
    

def process_var(var, opts):

    if 'dB' in var:
        var_name = var + '_mlat_%.3f_mlon_%.3f'%(opts['mlat_dB'], opts['mlon_dB'])
    else:
        var_name = var

    filename_path = conf[run + "_derived"] + "cutplanes/" + var_name
    if not os.path.exists(filename_path):
        os.makedirs(filename_path)

    k = 1
    times = []
    minmax = {'min': np.nan*np.empty(len(files)),
              'max': np.nan*np.empty(len(files)),
              'probe': {
                          'ux': np.nan*np.empty(len(files)),
                          'bz': np.nan*np.empty(len(files)),
                          'n_events': np.nan*np.empty(len(files))
                    }
              }

    if plot_type == 2: # todo: optimize
        dto = util.time2datetime(util.CDFfilename2time(run, files[0]))
        dtf = util.time2datetime(util.CDFfilename2time(run, files[-1]))
        dts = []
        print('1st dts = ' + str(dts))
        data = {'ux': [], 'bz': [], 'n_events': []}

    for filename in files:
        print('k = ' + str(k))
        if k > opts['nf']:
            break


        filename_png = conf[run + "_derived"] + "cutplanes/" + var_name + "/" \
                        + '{0:s}-{1:s}-type_{2:d}_delta_{3:.3f}.png' \
                        .format(filename, var_name, plot_type, opts['delta'])

        pkl = conf[run + "_derived"] + "cutplanes/minmax/" + var_name + '.pkl'

        zticks = opts['zticks']
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

        if not opts['regen'] and os.path.exists(filename_png):
            # Image already exists
            print('Image exists for file {0:d}/{1:d} for {2:s} at resoluton {3:.2f}. Will not regenerate.'\
                  .format(k, opts['nf'], var, opts['delta']))
            k = k + 1
            continue

        times.append(util.CDFfilename2time(run, filename))

        util.dlfile(conf[run + '_cdf'] + filename, debug=debug)

        if plot_type == 1:
            info = cp.plot(run, times[-1], var, opts['plane'],
                             dx=opts['delta'], dy=opts['delta'],
                             xlims=opts['xlims'], ylims=opts['ylims'],
                             dpi=opts['dpi'], showplot=opts['showplot'],
                             zticks=zticks, logz=True,
                             png=True, pngfile=filename_png, debug=debug,
                             mlat_dB=opts['mlat_dB'], mlon_dB=opts['mlon_dB'])
        elif plot_type == 2:

            def set_ylims(p, axes=None):

                if type(p) != list:
                    ps = [p]
                else:
                    ps = p

                import warnings

                var_mins = []
                var_maxes = []
                for p in ps:
                    with warnings.catch_warnings():
                        warnings.filterwarnings('ignore', r'All-NaN slice encountered')
                        var_mins.append(np.nanmin(minmax['probe'][p]))
                        var_maxes.append(np.nanmax(minmax['probe'][p]))
    
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore', r'All-NaN axis encountered')
                    var_min = np.nanmin(var_mins)
                    var_max = np.nanmax(var_maxes)

                if var_min != var_max and not np.isnan(var_min):
                    yticks = niceticks(var_min, var_max, 8, debug=debug)
                    axes.set_yticks(yticks)
                    axes.set_ylim([var_min, var_max])

            # Not sure why next two code lines needed. Is it due to pandas import
            # in another hapiclient function? Issue this addresses described
            # at https://www.gitmemory.com/issue/facebook/prophet/999/500035319
            #from pandas.plotting import register_matplotlib_converters
            #register_matplotlib_converters()
            
            print('2nd dts = ' + str(dts))
            print('len(dts) = ' + str(len(dts)))
            if len(dts) < 1: ##!!! just continues everything
                k = k+1
                continue
            print('3rd dts = ' + str(dts))

            from matplotlib import pyplot as plt
            from hapiclient.plot.datetick import datetick
            from cutplane import probe
        
            fig = plt.figure(figsize=(7, 9), dpi=opts['dpi'])
            axes = fig.add_axes([0.07, 0.55, 0.9, 0.4])
            info = cp.plot(times[-1], var, opts['plane'],
                             axes=axes, 
                             dx=opts['delta'], dy=opts['delta'],
                             xlims=opts['xlims'], ylims=opts['ylims'],
                             logz=True,
                             showplot=opts['showplot'],
                             zticks=zticks,
                             debug=debug, 
                             mlat_dB=opts['mlat_dB'], mlon_dB=opts['mlon_dB'])
            
            time = filename2time(filename)
            dts.append(util.time2datetime(time))
            
            # Sample upstream solar wind
            # TODO: Generalize code to allow plotting arbitrary number
            # of variables.
            d = probe(time, (31.5, 0, 0), ['ux', 'bz'], dictionary=True)

            #import pdb;pdb.set_trace()
            # Index of event with matching yr, mo, day, hr
            idx = np.all(times[-1][0:5] == event_list[:, 0:5], axis=1)
            if np.any(idx):
                # At least one event occured at this time step.
                d['n_events'] = event_list[idx==True, -1]
            else:
                d['n_events'] = 0;

            data['ux'].append(d['ux'])
            data['bz'].append(d['bz'])
            data['n_events'].append(d['n_events'])

            meta = util.filemeta(filename)
        
            axes = fig.add_axes([0.07, 0.1, 0.9, 0.18])
            axes.plot(dts, data['bz'], 'k')
            axes.set_xlim([dto, dtf])
            axes.grid()
            name = meta["parameters"]['bz']['plot_name']
            units = meta["parameters"]['bz']['vis_unit']
            datetick('x', axes=axes, debug=False)
            set_ylims(['bz', 'n_events'], axes=axes)
            
            for i in range(len(dts)):
                axes.plot([dts[i], dts[i]], [0, data['n_events'][i]], 'b-')
            axes.legend([name + ' [' + units + ']', '# Events'])

                
            axes = fig.add_axes([0.07, 0.28, 0.9, 0.18])
            axes.plot(dts, data['ux'], 'k')
            axes.grid()
            name = meta["parameters"]['ux']['plot_name']
            units = meta["parameters"]['ux']['vis_unit']
            axes.legend([name + ' [' + units + ']'])
            axes.set_xlim([dto, dtf])
            datetick('x', axes=axes)
            axes.axes.xaxis.set_ticklabels([])
            set_ylims('ux', axes=axes)

                            
            plt.savefig(filename_png, dpi=opts['dpi']) 
            
            if opts['showplot']:
                plt.show()
            else:
                plt.close()

        else:
            raise ValueError("plot_type = {0:d} is not valid".format(plot_type))

        # This information will be the same for all vars but will
        # be stored in each variable's .pkl file.

        if plot_type == 2:
            minmax['probe']['ux'][k-1] = d['ux']
            minmax['probe']['bz'][k-1] = d['bz']
            minmax['probe']['n_events'][k-1] = d['n_events']

        minmax['min'][k-1] = info['min']
        minmax['max'][k-1] = info['max']

        if debug:
            print("Writing " + pkl)
        with open(pkl, 'wb') as handle:
            pickle.dump(minmax, handle)
        if debug:
            print("Wrote " + pkl)

        print('Finished file {0:d}/{1:d} for {2:s} at resoluton {3:.5f}'\
               .format(k, opts['nf'], var, opts['delta']))
        
        k = k + 1

    return minmax

def process_all(opts):
    if para:
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        if num_cores is not None and num_cores > len(variables):
            num_cores = len(variables)
        print('Parallel processing {0:d} variable(s) using {1:d} cores'\
              .format(len(opts.keys()), num_cores))
        results = Parallel(n_jobs=num_cores)(\
                    delayed(process_var)(var, opts[var]) for var in variables)
    else:
        print('Serial processing {0:d} variable(s).'.format(len(opts.keys())))
        i = 0
        results = []
        for var in variables:
            results.append(process_var(var, opts[var]))
            i = i + 1
    return results

# Disable output buffering (sometimes causes problems with parallel runs)
# https://stackoverflow.com/questions/107705/disable-output-buffering
try:
    # Throws error in IPython, so put in try/catch.
    sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
except:
    pass


# Create output directory if needed
filename_path = conf[run + "_derived"] + "cutplanes/"
if not os.path.exists(filename_path):
    os.makedirs(filename_path)

# Create directory for variable min/max pkl files if needed
filename_path = conf[run + "_derived"] + "cutplanes/minmax/"
if not os.path.exists(filename_path):
    os.makedirs(filename_path)

print('Image files will be written to subdirectory of\n' + filename_path)
# Create directory for each variable
#for var in variables:
#    filename_path = conf[run + "_derived"] + "cutplanes/" + var
#    if not os.path.exists(filename_path):
#        os.makedirs(filename_path)


files = list(util.get_available_slices(run)[0])

print('{0:d} files found.'.format(len(files)))

if opts['nf'] is None:
    # Process all files
    opts['nf'] = len(files)
else:
    print('Processing only first {0:d} files.'.format(opts['nf']))


# Create options dict for each variable.
Opts = {}
for var in variables:
    # .copy to ensure thread safety.
    # Currently not needed b/c opts is not modified.
    Opts[var] = opts.copy()
    if 'dB' in var:
        import magnetometers as mg
        if run == 'CARR_IMPULSE':
            fixed_time = (2019,9,2,6,30,0)
            pos = mg.GetMagnetometerLocation('colaba', fixed_time, 'MAG', 'sph')
        Opts[var]['mlat_dB'] = 0. #pos[1]  # HERES WHERE MLAT MLON ARE SET
        Opts[var]['mlon_dB'] = 0. #pos[2]
    else:
        Opts[var]['mlat_dB'] = None
        Opts[var]['mlon_dB'] = None


if not second_only:
    print('First processing.')
    for var in variables:
        Opts[var]['delta'] = opts['delta1']
        Opts[var]['dpi'] = opts['dpi1']

    print(Opts)
    process_all(Opts)

if not first_only:
    print('Second processing.')        
    for var in variables:
        Opts[var]['delta'] = opts['delta2']
        Opts[var]['dpi'] = opts['dpi2']

    process_all(Opts)

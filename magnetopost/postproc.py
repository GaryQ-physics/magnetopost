import os
from magnetopost import util, models
#from magnetopost.summarize import slice_summary, stitch_summary
from magnetopost.gap_integrals import slice_bs_fac, stitch_bs_fac
from magnetopost.ionosphere_integrals import slice_bs_pedersen, stitch_bs_pedersen, slice_bs_hall, stitch_bs_hall


def job_msph(points, do_summary=False, cutplanes=None):
    run = util.parseConfig(f'{os.path.abspath(".")}/') #returns dictionary
    times = util.get_available_times(run)

    def wrap(time):
        ms_slice = models.get_ms_slice_class(run, time)
        #if do_summary:
            #slice_summary(run, time, clss =ms_slice)

        for point in points:
            slice_bs_fac(run, time, ms_slice, point)

    for time in times:
        wrap(time)

    #if do_summary: stitch_summary(run, times)

    for point in points:
        stitch_bs_fac(run, times, point)


def job_iono(points):
    run = util.prep_run(f'{os.path.abspath(".")}/') #returns dictionary
    times = run['ionosphere_files'].keys()

    def wrap(time):
        ie_slice = models.get_iono_slice(run, time)
        for point in points:
            slice_bs_pedersen(run,time, ie_slice, point)
            slice_bs_hall(run,time, ie_slice, point)

    for time in times:
        wrap(time)

    for point in points:
        stitch_bs_pedersen(run, times, point)
        stitch_bs_hall(run, times, point)


def job(points, do_summary=False, cutplanes=None, stitch_only=False):
    run = util.prep_run(f'{os.path.abspath(".")}/') #returns dictionary
    times = run['ionosphere_files'].keys()

    def wrap(time):
        try:
            ms_slice = models.get_ms_slice_class(run, time)
            ie_slice = models.get_iono_slice(run, time)
        except:
            print('FAILED')

        for point in points:
            slice_bs_pedersen(run,time, ie_slice, point)
            slice_bs_hall(run, time, ie_slice, point)
            slice_bs_fac(run, time, ms_slice, point)

        #if do_summary:
        #    slice_summary(run, time, ms_slice)

    if not stitch_only:
        if False:
            for time in times:
                wrap(time)
        else:
            from joblib import Parallel, delayed
            import multiprocessing
            num_cores = multiprocessing.cpu_count()
            if num_cores is None:
                assert(False)
            num_cores = min(num_cores, len(times), 20)
            print(f'Parallel processing {len(times)} time slices using {num_cores} cores')
            Parallel(n_jobs=num_cores)(\
                      delayed(wrap)(time) for time in list(times))

    for point in points:
        stitch_bs_pedersen(run, times, point)
        stitch_bs_hall(run, times, point)
        stitch_bs_fac(run, times, point)


    #if do_summary:
    #    stitch_summary(run, times)

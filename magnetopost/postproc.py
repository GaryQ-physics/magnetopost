import os
from magnetopost import util, models
from magnetopost.summarize import slice_summary, stitch_summary
from magnetopost.gap_integrals import slice_bs_fac, stitch_bs_fac, slice_helm_rCurrents, stitch_helm_rCurrents
from magnetopost.ionosphere_integrals import slice_bs_pedersen, stitch_bs_pedersen, slice_bs_hall, stitch_bs_hall
from magnetopost.magnetosphere_integrals import slice_bs_msph, stitch_bs_msph, slice_cl_msph, stitch_cl_msph
from magnetopost.boundary_integrals import slice_helm_outer, stitch_helm_outer
from magnetopost.probe import slice_probe, stitch_probe

def job_ie(points, stitch_only=False):
    run = util.prep_run(f'{os.path.abspath(".")}/') #returns dictionary
    times = run['ionosphere_files'].keys()

    def wrap(time):
        ie_slice = models.get_iono_slice(run, time)
        #print("Working on ionosphere at time = {}".format(time))

        for point in points:
            slice_bs_pedersen(run,time, ie_slice, point)
            slice_bs_hall(run, time, ie_slice, point)

    if not stitch_only:
        if False:
            for time in times:
                wrap(time)
        else:
            from joblib import Parallel, delayed
            import multiprocessing
            from tqdm import tqdm
            # https://stackoverflow.com/a/49950707
            num_cores = multiprocessing.cpu_count()
            num_cores = min(num_cores, len(times), 20)
            print(f'Parallel processing {len(times)} ionosphere slices using {num_cores} cores')
            Parallel(n_jobs=num_cores)(\
                      delayed(wrap)(time) for time in tqdm(times))

    for point in points:
        stitch_bs_hall(run, times, point)
        stitch_bs_pedersen(run, times, point)


def job_ms(points, do_summary=False, cutplanes=None, stitch_only=False):

    run = util.prep_run(f'{os.path.abspath(".")}/') #returns dictionary
    times = run['magnetosphere_files'].keys()

    def wrap(time):
        ms_slice = models.get_ms_slice_class(run, time)
        
        #print("Working on magnetosphere at time = {}".format(time))
        for point in points:
            slice_bs_fac(run,time, ms_slice, point)
            slice_bs_msph(run,time, ms_slice, point)
            slice_cl_msph(run,time, ms_slice, point)
            slice_helm_outer(run, time, ms_slice, point)
            slice_helm_rCurrents(run,time, ms_slice, point)
            slice_helm_rCurrents(run,time, ms_slice, point, gap_csys='GSM')
            slice_probe(run, time, ms_slice, point)

        if do_summary:
            slice_summary(run, time, ms_slice)

    if not stitch_only:
        if False:
            for time in times:
                wrap(time)
        else:
            from joblib import Parallel, delayed
            import multiprocessing
            from tqdm import tqdm
            num_cores = multiprocessing.cpu_count()
            num_cores = min(num_cores, len(times), 20)
            print(f'Parallel processing {len(times)} magnetosphere slices using {num_cores} cores')
            Parallel(n_jobs=num_cores)(\
                      delayed(wrap)(time) for time in tqdm(times))

    for point in points:
        stitch_bs_fac(run, times, point)
        stitch_bs_msph(run, times, point)
        stitch_cl_msph(run, times, point)
        stitch_helm_outer(run, times, point)
        stitch_helm_rCurrents(run, times, point)
        stitch_helm_rCurrents(run, times, point, gap_csys='GSM')
        stitch_probe(run, times, point)

    if do_summary:
        stitch_summary(run, times)

def job(points):
    ''' example points:
     points = ("YKC","YKC_N","YKC_S","OTT","FRD")
     points = ("GMpoint4","GMpoint5")
    '''
    job_ie(points)
    job_ms(points)
    print('DONE JOB')

import os
from magnetopost import util, models
#from magnetopost.summarize import slice_summary, stitch_summary
from magnetopost.gap_integrals import slice_bs_fac, stitch_bs_fac
from magnetopost.ionosphere_integrals import slice_bs_pedersen, stitch_bs_pedersen, slice_bs_hall, stitch_bs_hall
from magnetopost.magnetosphere_integrals import slice_bs_msph, stitch_bs_msph

def job_ie(points, stitch_only=False):
    run = util.prep_run(f'{os.path.abspath(".")}/') #returns dictionary
    times = run['ionosphere_files'].keys()

    def wrap(time):
        ie_slice = models.get_iono_slice(run, time)

        for point in points:
            slice_bs_pedersen(run,time, ie_slice, point)
            slice_bs_hall(run, time, ie_slice, point)

    if not stitch_only:
        if False:
            for time in times:
                wrap_ie(time)
        else:
            from joblib import Parallel, delayed
            import multiprocessing
            num_cores = multiprocessing.cpu_count()
            num_cores = min(num_cores, len(times), 20)
            print(f'Parallel processing {len(times)} ionosphere slices using {num_cores} cores')
            Parallel(n_jobs=num_cores)(\
                      delayed(wrap)(time) for time in times)

    for point in points:
        stitch_bs_hall(run, times, point)
        stitch_bs_pedersen(run, times, point)


def job_ms(points, do_summary=False, cutplanes=None, stitch_only=False):
    run = util.prep_run(f'{os.path.abspath(".")}/') #returns dictionary
    times = run['magnetosphere_files'].keys()

    def wrap(time):
        ms_slice = models.get_ms_slice_class(run, time)

        for point in points:
            slice_bs_fac(run,time, ms_slice, point)
            slice_bs_msph(run,time, ms_slice, point)

    if not stitch_only:
        if False:
            for time in times:
                wrap_ie(time)
        else:
            from joblib import Parallel, delayed
            import multiprocessing
            num_cores = multiprocessing.cpu_count()
            num_cores = min(num_cores, len(times), 20)
            print(f'Parallel processing {len(times)} ionosphere slices using {num_cores} cores')
            Parallel(n_jobs=num_cores)(\
                      delayed(wrap)(time) for time in times)

    for point in points:
        stitch_bs_fac(run, times, point)
        stitch_bs_msph(run, times, point)


def main():
    job_ie(("colaba",))
    job_ms(("colaba",))

if __name__ == '__main__':
    main()
    print('DONE')

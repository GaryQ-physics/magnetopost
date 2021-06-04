import os
from magnetopost import util, models
#from magnetopost.summarize import slice_summary, stitch_summary
from magnetopost.ionosphere_integrals import slice_bs_pedersen


def job(points, do_summary=False, cutplanes=None):
    run = util.parseConfig(f'{os.path.abspath(".")}/') #returns dictionary
    times = util.get_available_times(run)

    def wrap(time):
        ms_slice = models.get_ms_slice_class(run, time)
        if do_summary:
            slice_summary(run, time, clss = ms_slice)

        for point in points:
            slice_biotsavart(point)

    for time in times:
        wrap(time)

    if do_summary: stitch_summary(run, times)


def job_iono(points):
    run = util.prep_run(f'{os.path.abspath(".")}/') #returns dictionary
    times = run['ionosphere_files'].keys()

    def wrap(time):
        ie_slice = models.get_iono_slice(run, time)
        for point in points:
            slice_bs_pedersen(run,time, ie_slice, point)

    for time in times:
        wrap(time)

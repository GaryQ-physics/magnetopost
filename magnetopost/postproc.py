import logging
from magnetopost import util, models
from magnetopost.summarize import slice_summary, stitch_summary
from magnetopost.gap_integrals import slice_bs_fac, stitch_bs_fac, slice_helm_rCurrents, stitch_helm_rCurrents
from magnetopost.ionosphere_integrals import slice_bs_pedersen, stitch_bs_pedersen, slice_bs_hall, stitch_bs_hall
from magnetopost.magnetosphere_integrals import slice_bs_msph, stitch_bs_msph, slice_cl_msph, stitch_cl_msph
from magnetopost.boundary_integrals import slice_helm_outer, stitch_helm_outer
from magnetopost.probe import slice_probe, stitch_probe


def job_ie(info, points, n_steps=None, stitch_only=False):

    run = util.prep_run(info["dir_run"]) # returns dictionary
    times = list(run['ionosphere_files'].keys())

    if n_steps is not None:
        times = times[0:n_steps]

    def wrap(time):
        logging.info("Working on ionosphere at time = {}".format(time))
        logging.info("Reading data for = {}".format(time))

        ie_slice = models.get_iono_slice(info, run, time)
        for point in points:
            slice_bs_pedersen(run, time, ie_slice, point)
            slice_bs_hall(run, time, ie_slice, point)

    if not stitch_only:
        if False:
            for time in times:
                wrap(time)
        else:
            from joblib import Parallel, delayed
            import multiprocessing
            num_cores = multiprocessing.cpu_count()
            num_cores = min(num_cores, len(times), 20)
            logging.info(f'Parallel processing {len(times)} ionosphere slices using {num_cores} cores')
            Parallel(n_jobs=num_cores)(delayed(wrap)(time) for time in times)

    for point in points:
        stitch_bs_pedersen(run, times, point)
        stitch_bs_hall(run, times, point)

    logging.info("Finished ionosphere processing")


def job_ms(info, points, n_steps=None, stitch_only=False, do_summary=False):


    run = util.prep_run(info["dir_run"]) # returns dict
    times = list(run['magnetosphere_files'].keys())

    if n_steps is not None:
        times = times[0:n_steps]

    def wrap(time):
        logging.info("Working on magnetosphere at time = {}".format(time))
        logging.info("Reading data for = {}".format(time))
        ms_slice = models.get_ms_slice_class(info, run, time)

        logging.info("Computing integrals for = {}".format(time))
        for point in points:
            slice_bs_fac(run, time, ms_slice, point)
            slice_bs_msph(run, time, ms_slice, point)
            slice_cl_msph(run, time, ms_slice, point)
            slice_helm_outer(run, time, ms_slice, point)
            slice_helm_rCurrents(run, time, ms_slice, point)
            slice_helm_rCurrents(run, time, ms_slice, point, gap_csys='GSM')
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
            num_cores = multiprocessing.cpu_count()
            num_cores = min(num_cores, len(times), 20)
            logging.info(f'Parallel processing {len(times)} magnetosphere slices using {num_cores} cores')
            Parallel(n_jobs=num_cores)(delayed(wrap)(time) for time in times)

    logging.info("Finished magnetosphere processing")


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


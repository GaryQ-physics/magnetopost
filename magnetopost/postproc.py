import logging
from magnetopost import models
from magnetopost.summarize import slice_summary, stitch_summary
from magnetopost.gap_integrals import slice_bs_fac, slice_helm_rCurrents, stitch_helm_rCurrents
from magnetopost.magnetosphere_integrals import slice_bs_msph, stitch_bs_msph, slice_cl_msph, stitch_cl_msph
from magnetopost.probe import slice_probe, stitch_probe


def cat(info, times, points, function_str):

    import numpy as np
    from magnetopost import util

    for point in points:

        integrals = []
        for time in times:
            outname = f'{info["dir_run"]}/derived/timeseries/timesteps/' \
                + f'{function_str}-{point}-{util.Tstr(time)}.npy'

            integrals.append(np.load(outname))

        outfile = f'{info["dir_run"]}/derived/timeseries/' \
                + f'{function_str}-{point}.npy'
        print(f"Writing {outfile}")
        np.save(outfile, np.array(integrals))


def job_ie(info, points, n_steps=None, stitch_only=False, parallel=True):

    from magnetopost.ionosphere_integrals import bs_pedersen, bs_hall

    times = list(info['files']['ionosphere'].keys())

    if n_steps is not None:
        times = times[0:n_steps]

    def wrap(time):

        logging.info("Working on ionosphere at time = {}".format(time))

        ie_slice = models.get_iono_slice(info, time)
        for point in points:
            bs_pedersen(info, time, ie_slice, point)
            bs_hall(info, time, ie_slice, point)

    if parallel:
        from joblib import Parallel, delayed
        import multiprocessing
        num_cores = multiprocessing.cpu_count()
        num_cores = min(num_cores, len(times), 20)
        logging.info(f'Parallel processing {len(times)} ionosphere timesteps using {num_cores} cores')
        Parallel(n_jobs=num_cores)(delayed(wrap)(time) for time in times)
    else:
        logging.info(f'Serial processing {len(times)} ionosphere timesteps')
        for time in times:
            wrap(time)

    logging.info("Finished ionosphere processing")

    cat(info, times, points, 'bs_pedersen')
    cat(info, times, points, 'bs_hall')



def job_ms(info, points, n_steps=None, stitch_only=False, do_summary=False):

    from magnetopost.boundary_integrals import helm_outer

    times = list(info['files']['magnetosphere'].keys())

    if n_steps is not None:
        times = times[0:n_steps]

    def wrap(time):
        logging.info("Working on magnetosphere at time = {}".format(time))
        logging.info("Reading data for = {}".format(time))
        ms_slice = models.get_ms_slice_class(info, time)

        logging.info("Computing integrals for = {}".format(time))
        for point in points:
            slice_bs_fac(info, time, ms_slice, point)
            slice_bs_msph(info, time, ms_slice, point)
            slice_cl_msph(info, time, ms_slice, point)
            helm_outer(info, time, ms_slice, point)
            slice_helm_rCurrents(info, time, ms_slice, point)
            slice_helm_rCurrents(info, time, ms_slice, point, gap_csys='GSM')
            slice_probe(info, time, ms_slice, point)

        if do_summary:
            slice_summary(info, time, ms_slice)

    if not stitch_only:
        if False:
            for time in times:
                wrap(time)
        else:
            from joblib import Parallel, delayed
            import multiprocessing
            num_cores = multiprocessing.cpu_count()
            num_cores = min(num_cores, len(times), 20)
            logging.info(f'Parallel processing {len(times)} magnetosphere timesteps using {num_cores} cores')
            Parallel(n_jobs=num_cores)(delayed(wrap)(time) for time in times)

    logging.info("Finished magnetosphere processing")

    cat(info, times, points, 'helm_outer')
    cat(info, times, points, 'bs_fac')

    for point in points:
        #stitch_bs_fac(info, times, point)
        stitch_bs_msph(info, times, point)
        stitch_cl_msph(info, times, point)
        #stitch_helm_outer(info, times, point)
        stitch_helm_rCurrents(info, times, point)
        stitch_helm_rCurrents(info, times, point, gap_csys='GSM')
        stitch_probe(info, times, point)

    if do_summary:
        stitch_summary(info, times)


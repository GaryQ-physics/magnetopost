import logging
import numpy as np
import matplotlib.pyplot as plt

import magnetopost as mp

def write_plot(fig, outfile):
    logging.info("Writing {}.[svg,pdf,png]".format(outfile))
    fig.savefig(outfile + '.svg')
    fig.savefig(outfile + '.pdf')
    fig.savefig(outfile + '.png')


def norm(df):
    return np.sqrt(df['north']**2+df['east']**2+df['down']**2)


def gen_rmse(diff):
    return np.sqrt( np.sum(diff**2)/diff.size ) #!!! NANS


def surf_point(info, surface_location, n_steps=None):

    #logging.info("rundir = {}".format(dir_run))

    if info['file_type'] == "cdf":
        dBMhd, dBFac, dBHal, dBPed = mp.extract_magnetometer_data.extract_from_swmf_ccmc_printout_file(info, surface_location, n_steps=n_steps)
    if info['file_type'] == "out":
        dBMhd, dBFac, dBHal, dBPed = mp.extract_magnetometer_data.extract_from_swmf_magnetometer_files(info, surface_location, n_steps=n_steps)

    bs_msph, bs_fac, bs_hall, bs_pedersen, cl_msph, helm_outer, helm_rCurrents_gapSM, probe = mp.extract_magnetometer_data.extract_from_magnetopost_files(info, surface_location, n_steps=n_steps)

    B_G  = bs_fac + bs_hall + bs_pedersen
    #B_G2 = dBFac + dBHal+ dBPed

    fig, axs = plt.subplots(nrows=4, ncols=2, sharex=True, figsize=(12,12), dpi=100)

    def foo(i, swmf, ours, title):
        norm(swmf).plot(ax=axs[i,0], label='SWMF magnetometer files', color='Orange')
        norm(ours).plot(ax=axs[i,0], label='Our calculation', color='Blue')
        diff = norm(ours) - norm(swmf)
        rmse = gen_rmse(diff)
        diff.plot(ax=axs[i,1], label=f'Difference (RMSE={rmse:.1f})')

        axs[i,0].set_title(title)
        axs[i,0].legend(title_fontsize=1)
        axs[i,1].legend(title_fontsize=1)
        axs[i,0].set_ylabel('nT')

    foo(0, dBMhd, bs_msph    , 'dBMhd')
    foo(1, dBFac, bs_fac     , 'dBFac')
    foo(2, dBHal, bs_hall    , 'dBHal')
    foo(3, dBPed, bs_pedersen, 'dBPed')

    #fig.suptitle(dir_run)
    #datetick('x')

    outfile = f'{info["dir_run"]}/derived/figures/{info["run_name"]}-compare_with_swmf'
    write_plot(fig, outfile)


    fig.clf(); del fig

    fig, axs = plt.subplots(nrows=5, ncols=1, sharex=True, figsize=(12,12), dpi=100)

    method_3 = B_G + bs_msph
    method_2 = B_G + bs_msph + cl_msph + helm_outer
    method_1 = B_G + helm_rCurrents_gapSM

    norm(method_1).plot(ax=axs[0], label='Method 1.', color='Blue')
    norm(method_2).plot(ax=axs[0], label='Method 2.', color='Orange')

    diff = norm(method_1)-norm(method_2)
    rmse = gen_rmse(diff)
    diff.plot(ax=axs[1], label=f'(Method 1.) - (Method 2.) (RMSE={rmse:.1f})')

    norm(method_1).plot(ax=axs[2], label='Method 1.', color='Blue')
    norm(method_3).plot(ax=axs[2], label='Method 3.', color='Orange')

    diff = norm(method_1)-norm(method_3)
    rmse = gen_rmse(diff)
    diff.plot(ax=axs[3], label=f'(Method 1.) - (Method 3.) (RMSE={rmse:.1f})')

    norm(cl_msph).plot(ax=axs[4], label=r'$\int_{\mathcal{M}}$Coulomb', color='Orange')
    norm(helm_outer).plot(ax=axs[4], label=r'$\oint_{\mathcal{O}}$', color='Green')
    norm(bs_msph).plot(ax=axs[4], label=r'$\int_{\mathcal{M}}$Biot_Savart',  color='Blue')

    [ax.legend() for ax in axs]
    [ax.set_ylabel('nT') for ax in axs]
    #fig.suptitle(runname)
    #datetick('x')

    outfile = f'{info["dir_run"]}/derived/figures/{info["run_name"]}-compare_methods_123'
    write_plot(fig, outfile)


def msph_point(info, surface_location, n_steps=None):

    bs_msph, bs_fac, bs_hall, bs_pedersen,  cl_msph, helm_outer, helm_rCurrents_gapSM, probe = mp.extract_magnetometer_data.extract_from_magnetopost_files(info, surface_location, n_steps=n_steps)
    B_G  = bs_fac + bs_hall + bs_pedersen
    #B_G2 = dBFac + dBHal+ dBPed

    fig, axs = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(12,12), dpi=100)

    method_A  = bs_msph + cl_msph + helm_outer - helm_rCurrents_gapSM
    method_B  = bs_msph + cl_msph + helm_outer + B_G
    #method_B2 = bs_msph + cl_msph + helm_outer + B_G2

    norm(method_A).plot(ax=axs[0],
                                label='Method A.',  color='Blue')
    norm(probe).plot(ax=axs[0],
                                label=r'$\mathbf{B}$',  color='Orange')

    diff = norm(method_A)-norm(probe)
    rmse = gen_rmse(diff)
    diff.plot(ax=axs[1],
                                label=r'(Method A.) - $\mathbf{B}$'+f' (RMSE={rmse:.1f})')

    norm(method_B).plot(ax=axs[2],
                                label='Method B.',  color='Blue')
    norm(probe).plot(ax=axs[2],
                                label=r'$\mathbf{B}$',  color='Orange')

    diff = norm(method_B)-norm(probe)
    rmse = gen_rmse(diff)
    diff.plot(ax=axs[3],
                                label=r'(Method B.) - $\mathbf{B}$'+f' (RMSE={rmse:.1f})')

    [ax.legend() for ax in axs]
    [ax.set_ylabel('nT') for ax in axs]
    fig.suptitle(info["run_name"])
    #datetick('x')

    outfile = f'{info["dir_run"]}/derived/figures/{info["run_name"]}-compare_AB'
    logging.info(f"Writing {outfile}-compareAB.[pdf,png]")
    write_plot(fig, outfile)


    fig, axs = plt.subplots(nrows=2, ncols=1, sharex=True, figsize=(12,12), dpi=100)
    norm(helm_rCurrents_gapSM).plot(ax=axs[0],
                                label=r'$\oint_{\mathcal{I}}$',  color='Blue')
    norm(B_G).plot(ax=axs[0],
                                label=r'$\mathbf{B}_{\mathcal{G}}$',  color='Orange')
    diff = norm(helm_rCurrents_gapSM) - norm(B_G)
    rmse = gen_rmse(diff)
    diff.plot(ax=axs[1],
                                label='difference')
    #datetick('x')

    outfile = f'{info["dir_run"]}/derived/figures/{info["run_name"]}-consistency'
    logging.info(f"Writing {outfile}-compareAB.[pdf,png]")
    write_plot(fig, outfile)

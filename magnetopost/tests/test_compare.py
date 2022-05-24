#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 2 12:42:32 2022

@author: Dean Thomas

Compare results from separate runs of magnetopost, allows regression 
analysis to ensure changes to code do not lead to unexpected results.
"""

import logging
import numpy as np
import matplotlib.pyplot as plt

import magnetopost as mp
import magnetopost.util as util
    
def _write_plot(fig, outfile):
    logging.info("Writing {}.[svg,pdf,png]".format(outfile))
    fig.savefig(outfile + '.svg')
    fig.savefig(outfile + '.pdf')
    fig.savefig(outfile + '.png')

def _norm(df):
    return np.sqrt(df['north']**2+df['east']**2+df['down']**2)

def _gen_rms_diff(diff):
    return np.sqrt( np.sum(diff**2)/diff.size )

def term_comparison(old_info, new_info, point, n_steps=None):
    """Generates plots to compare numeric integrations from different runs
    
    In this routine, we compare the individual terms such as 
    bs_fac (Biot-Savart for FAC), bs_hall (Biot-Savart for Hall), ...

    Arguments:
        
    dir_old - directory containing derived directory with old run
    dir_new - directory containing derived directory with new run
    point = location used in runs
    n_steps = how many time slices to read, None = all slices
    """
    
    logging.info('Terms for ionosphere comparison')
    logging.info("old directory ={}".format(old_info["dir_run"]))
    logging.info("new directory ={}".format(new_info["dir_run"]))

    # Read data from old run for comparison to new run
    bs_msph_old, bs_fac_old, bs_hall_old, bs_pedersen_old, cl_msph_old, \
        helm_outer_old, helm_rCurrents_gapSM_old, probe_old \
            = mp.extract_magnetometer_data.extract_from_magnetopost_files(old_info, point, n_steps=n_steps)

    # Read data from new run
    bs_msph_new, bs_fac_new, bs_hall_new, bs_pedersen_new, cl_msph_new, \
        helm_outer_new, helm_rCurrents_gapSM_new, probe_new \
            = mp.extract_magnetometer_data.extract_from_magnetopost_files(new_info, point, n_steps=n_steps)

    # Create figure for plots
    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(12,12), dpi=300)

    # Two routines for repetitive plot creation
    def add_plot(i, old, new, title, point):
        _norm(old).plot(ax=axs[i,0], label='Original', color='Orange')
        _norm(new).plot(ax=axs[i,0], label='New', color='Blue', linestyle='dotted')
        diff = _norm(new) - _norm(old)
        rms_diff = _gen_rms_diff(diff)
        diff.plot(ax=axs[i,1], label=f'Difference (RMS={rms_diff:.1f})')

        axs[i,0].set_title(title + " for " + point)
        axs[i,0].legend(title_fontsize=1)
        axs[i,1].legend(title_fontsize=1)
        axs[i,0].set_ylabel('nT')

    def save_plot(i, point):
        outfile = f'{new_info["dir_run"]}/derived/figures/terms_compare_{point}_{i}'
        _write_plot(fig, outfile)
    
    # Create plots comparing new to old and showing RMS differences
    # Compare all eight variables
    add_plot(0, bs_msph_old, bs_msph_new, 'bs_msph', point)
    add_plot(1, bs_fac_old, bs_fac_new, 'bs_fac', point)
    save_plot(1, point)

    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(12,12), dpi=300)
    add_plot(0, bs_hall_old, bs_hall_new, 'bs_hall', point)
    add_plot(1, bs_pedersen_old, bs_pedersen_new, 'bs_pedersen', point)
    save_plot(2, point)
    
    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(12,12), dpi=300)
    add_plot(0, cl_msph_old, cl_msph_new, 'cl_msph', point)
    add_plot(1, helm_outer_old, helm_outer_old, 'helm_outer', point)
    save_plot(3, point)
   
    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(12,12), dpi=300)
    add_plot(0, helm_rCurrents_gapSM_old, helm_rCurrents_gapSM_new, 'helm_rCurrents', point)
    add_plot(1, probe_old, probe_new, 'probe', point)
    save_plot(4, point)
    
    

def method_comparison(old_info, new_info, point, n_steps=None):
    """Generates plots to compare methods from different runs
    
    In this routine, we compare the methods described in the paper
    Method A, Method B, Method 1, ...
    
    Arguments:
        
    dir_old - directory containing derived directory with old run
    dir_new - directory containing derived directory with new run
    point = location used in runs
    n_steps = how many time slices to read, None = all slices
    """

    logging.info('Methods for ionosphere comparison')
    logging.info("old directory ={}".format(old_info["dir_run"]))
    logging.info("new directory ={}".format(new_info["dir_run"]))

    # Read data from old run for comparison to new run
    bs_msph_old, bs_fac_old, bs_hall_old, bs_pedersen_old, cl_msph_old, \
        helm_outer_old, helm_rCurrents_gapSM_old, probe_old \
            = mp.extract_magnetometer_data.extract_from_magnetopost_files(old_info, point, n_steps=n_steps)

    # Read data from new run
    bs_msph_new, bs_fac_new, bs_hall_new, bs_pedersen_new, cl_msph_new, \
        helm_outer_new, helm_rCurrents_gapSM_new, probe_new \
            = mp.extract_magnetometer_data.extract_from_magnetopost_files(new_info, point, n_steps=n_steps)

    # Create figure for plots
    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(12,12), dpi=300)

    # Two routines for repetitive plot creation
    def add_plot(i, old, new, title, point):
        _norm(old).plot(ax=axs[i,0], label='Original', color='Orange')
        _norm(new).plot(ax=axs[i,0], label='New', color='Blue', linestyle='dotted')
        diff = _norm(new) - _norm(old)
        rms_diff = _gen_rms_diff(diff)
        diff.plot(ax=axs[i,1], label=f'Difference (RMS={rms_diff:.1f})')

        axs[i,0].set_title(title + ' for ' + point)
        axs[i,0].legend(title_fontsize=1)
        axs[i,1].legend(title_fontsize=1)
        axs[i,0].set_ylabel('nT')

    def save_plot(i, point):
        outfile = f'{new_info["dir_run"]}/derived/figures/method_compare_{point}_{i}'
        _write_plot(fig, outfile)
        
    B_G_old  = bs_fac_old + bs_hall_old + bs_pedersen_old
    B_G_new  = bs_fac_new + bs_hall_new + bs_pedersen_new

    method_A_old  = bs_msph_old + cl_msph_old + helm_outer_old - helm_rCurrents_gapSM_old
    method_B_old  = bs_msph_old + cl_msph_old + helm_outer_old + B_G_old

    method_A_new  = bs_msph_new + cl_msph_new + helm_outer_new - helm_rCurrents_gapSM_new
    method_B_new  = bs_msph_new + cl_msph_new + helm_outer_new + B_G_new

    method_1_old = B_G_old + helm_rCurrents_gapSM_old
    method_2_old = B_G_old + bs_msph_old + cl_msph_old + helm_outer_old
    method_3_old = B_G_old + bs_msph_old

    method_1_new = B_G_new + helm_rCurrents_gapSM_new
    method_2_new = B_G_new + bs_msph_new + cl_msph_new + helm_outer_new
    method_3_new = B_G_new + bs_msph_new

    # Create plots comparing new to old and showing RMS differences
    # Compare all eight variables
    add_plot(0, B_G_old, B_G_new, 'B_G', point)
    add_plot(1, method_A_old, method_A_new, 'Method A', point)
    save_plot(1, point)
    
    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(12,12), dpi=300)
    add_plot(0, method_B_old, method_B_new, 'Method B', point)
    add_plot(1, method_1_old, method_1_new, 'Method 1', point)
    save_plot(2, point)
    
    fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(12,12), dpi=300)
    add_plot(0, method_2_old, method_2_new, 'Method 2', point)
    add_plot(1, method_3_old, method_3_new, 'Method 3', point)
    save_plot(3, point)
        
    
if __name__ == "__main__":
    points = ["colaba", "GMpoint1"]

    #DT Need to think of a better way to initialize the info blocks
    #DT This is a hack to do the comparison
    old_info = {
            "model": "SWMF",
            "run_name": "DIPTSUR2",
            "rCurrents": 1.8,
            "file_type": "out",
            "dir_run": "/Volumes/Physics HD/DIPTSUR2_original",
            "dir_plots": "/Volumes/Physics HD/DIPTSUR2_original/DIPTSUR2.plots"
    }

    new_info = {
            "model": "SWMF",
            "run_name": "DIPTSUR2",
            "rCurrents": 1.8,
            "file_type": "out",
            "dir_run": "/Volumes/Physics HD/DIPTSUR2_float64",
            "dir_plots": "/Volumes/Physics HD/DIPTSUR2_float64/DIPTSUR2.plots"
    }
    
    util.setup(old_info)
    util.setup(new_info)

    #DT Need to think about what is being compared.  As is, this compares
    #DT every term (e.g., bs_fac) and every method (e.g., Method A) for
    #DT all points.  This may not make sense for the methods.
    for point in points:
        term_comparison(old_info, new_info, point)
        method_comparison(old_info, new_info, point)


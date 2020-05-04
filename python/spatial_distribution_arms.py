#!/usr/bin/env python3
"""
Spatial distribution of energetic particles for ARMS MHD runs
"""
from __future__ import print_function

import argparse
import itertools
import json
import math
import multiprocessing

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from matplotlib.colors import LogNorm

from sde_util import load_mhd_config, mkdir_p

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

def plot_dist_2d(plot_config, mhd_config, show_plot=True):
    """Plot 2D spatial distribution

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    tframe = plot_config["tframe"]
    nreduce = plot_config["nreduce"]
    nr_mhd, = mhd_config["nx"]  # r
    nt_mhd, = mhd_config["ny"]  # theta
    np_mhd, = mhd_config["nz"]  # phi
    nrr = nr_mhd // nreduce
    ntr = nt_mhd // nreduce
    npr = np_mhd // nreduce

    mhd_run_dir = plot_config["mhd_run_dir"]
    fname = mhd_run_dir + "bin_data/xpos.dat"
    rpos = np.fromfile(fname)
    fname = mhd_run_dir + "bin_data/ypos.dat"
    tpos = np.fromfile(fname)
    fname = mhd_run_dir + "bin_data/zpos.dat"
    ppos = np.fromfile(fname)

    ng = 2  # number of ghost cells
    rpos_r = rpos[ng:-ng:nreduce]
    tpos_r = tpos[ng:-ng:nreduce] + math.pi * 0.5
    ppos_r = ppos[ng:-ng:nreduce]

    run_name = plot_config["run_name"]
    tframe_str = str(tframe).zfill(4)
    fname = '../data/' + run_name + '/fxy-' + tframe_str + '_sum.dat'
    fdata = np.fromfile(fname).reshape([-1, npr, ntr, nrr])
    p3d, t3d, r3d = np.meshgrid(ppos_r, tpos_r, rpos_r, indexing='ij')
    fdata /= r3d**2 * np.sin(t3d)
    nbands, _, _, _ = fdata.shape
    fdata_phi = np.sum(fdata, axis=1)
    fdata_theta = np.sum(fdata, axis=2)

    fname = '../data/' + run_name + '/fp-' + str(tframe).zfill(4) + '_sum.dat'
    data = np.fromfile(fname, dtype=np.float64)
    dsz, = data.shape
    nbins = dsz // 2
    pinit = 0.1
    pmom = data[0:nbins]
    pmin, pmax = pmom[0], pmom[-1]
    emin = pmin**2 / pinit**2
    emax = pmax**2 / pinit**2
    ebins_edges = np.logspace(math.log10(emin), math.log10(emax), nbands+1)
    nmin, nmax = 10, 2E3

    rr, tt = np.meshgrid(rpos_r, tpos_r)
    xx = rr * np.cos(tt)
    yy = rr * np.sin(tt)

    fnorms = [1, 4, 16, 64]

    fig = plt.figure(figsize=[15, 10])
    rect0 = [0.07, 0.55, 0.40, 0.40]
    hgap, vgap = 0.05, 0.05
    bands, bande = plot_config["high_bands"]
    for iband, eband in enumerate(range(bands, bande+1)):
        row = iband // 2
        col = iband % 2
        rect = np.copy(rect0)
        rect[0] = rect0[0] + col * (rect0[2] + hgap)
        rect[1] = rect0[1] - row * (rect0[3] + vgap)
        ax = fig.add_axes(rect)
        fdata_sum = fdata_phi[eband] * fnorms[iband]
        fdata1 = fdata_sum + 1E-5
        im = ax.pcolormesh(xx, yy, fdata1,
                           norm=LogNorm(vmin=nmin, vmax=nmax))
        ax.tick_params(labelsize=12)
        if row == 1:
            ax.set_xlabel(r'$x/R_e$', fontsize=16)
        else:
            ax.tick_params(axis='x', labelleft=False)
        if col == 0:
            ax.set_ylabel(r'$y/R_e$', fontsize=16)
        else:
            ax.tick_params(axis='y', labelleft=False)
        norm = fnorms[iband]
        e0 = plot_config["e0"]
        enel = "{%0.0f}" % (ebins_edges[eband] * e0)
        eneh = "{%0.0f}" % (ebins_edges[eband+1] * e0)
        ftext = (r'$' + enel + r'$keV' + r'$' +
                 r'< \varepsilon <' +
                 eneh + r'$keV')
        if norm > 1:
            fig_text = r'$' + str(fnorms[iband]) + 'n($' + ftext + r'$)$'
        else:
            fig_text = r'$n($' + ftext + r'$)$'
        ax.text(0.05, 0.85, fig_text, color='white', fontsize=12,
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
    rect_cbar = np.copy(rect0)
    rect_cbar[0] += rect0[2] * 2 + hgap + 0.02
    rect_cbar[2] = 0.01
    rect_cbar[1] = rect0[1] - vgap - rect0[3]*0.5
    rect_cbar[3] = rect0[3] + vgap
    cbar_ax = fig.add_axes(rect_cbar)
    cbar1 = fig.colorbar(im, cax=cbar_ax, extend='both')
    cbar1.ax.tick_params(labelsize=16)

    L0 = 6.9634E8  # m
    v0 = 1E3       # m/s
    t0 = L0 / v0
    tva = mhd_config["dt_out"][0] * t0
    title = 'Frame: %d' % plot_config["tframe"]
    fig.suptitle(title, fontsize=24)

    fdir = '../img/nrho/' + run_name + '/'
    mkdir_p(fdir)
    fname = fdir + 'nrho_bands_phi_' + tframe_str + '.jpg'
    fig.savefig(fname)

    rr, pp = np.meshgrid(rpos_r, ppos_r)
    xx = rr * np.sin(pp)
    yy = rr * np.cos(pp)

    fnorms = [1, 4, 16, 64]

    fig = plt.figure(figsize=[15, 10])
    rect0 = [0.07, 0.55, 0.40, 0.40]
    hgap, vgap = 0.05, 0.05
    bands, bande = plot_config["high_bands"]
    for iband, eband in enumerate(range(bands, bande+1)):
        row = iband // 2
        col = iband % 2
        rect = np.copy(rect0)
        rect[0] = rect0[0] + col * (rect0[2] + hgap)
        rect[1] = rect0[1] - row * (rect0[3] + vgap)
        ax = fig.add_axes(rect)
        fdata_sum = fdata_theta[eband] * fnorms[iband]
        fdata1 = fdata_sum + 1E-5
        im = ax.pcolormesh(xx, yy, fdata1,
                           norm=LogNorm(vmin=nmin, vmax=nmax))
        ax.tick_params(labelsize=12)
        if row == 1:
            ax.set_xlabel(r'$x/R_e$', fontsize=16)
        else:
            ax.tick_params(axis='x', labelleft=False)
        if col == 0:
            ax.set_ylabel(r'$z/R_e$', fontsize=16)
        else:
            ax.tick_params(axis='y', labelleft=False)
        norm = fnorms[iband]
        e0 = plot_config["e0"]
        enel = "{%0.0f}" % (ebins_edges[eband] * e0)
        eneh = "{%0.0f}" % (ebins_edges[eband+1] * e0)
        ftext = (r'$' + enel + r'$keV' + r'$' +
                 r'< \varepsilon <' +
                 eneh + r'$keV')
        if norm > 1:
            fig_text = r'$' + str(fnorms[iband]) + 'n($' + ftext + r'$)$'
        else:
            fig_text = r'$n($' + ftext + r'$)$'
        ax.text(0.05, 0.7, fig_text, color='white', fontsize=12,
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
    rect_cbar = np.copy(rect0)
    rect_cbar[0] += rect0[2] * 2 + hgap + 0.02
    rect_cbar[2] = 0.01
    rect_cbar[1] = rect0[1] - vgap - rect0[3]*0.5
    rect_cbar[3] = rect0[3] + vgap
    cbar_ax = fig.add_axes(rect_cbar)
    cbar1 = fig.colorbar(im, cax=cbar_ax, extend='both')
    cbar1.ax.tick_params(labelsize=16)

    L0 = 6.9634E8  # m
    v0 = 1E3       # m/s
    t0 = L0 / v0
    tva = mhd_config["dt_out"][0] * t0
    title = 'Frame: %d' % plot_config["tframe"]
    fig.suptitle(title, fontsize=24)

    fdir = '../img/nrho/' + run_name + '/'
    mkdir_p(fdir)
    fname = fdir + 'nrho_bands_theta_' + tframe_str + '.jpg'
    fig.savefig(fname)

    if show_plot:
        plt.show()
    else:
        plt.close("all")


def get_cmd_args():
    """Get command line arguments
    """
    default_mhd_run = 'arms_test'
    default_mhd_run_dir = ('/net/scratch3/xiaocan/mhd/' + default_mhd_run + '/')
    default_sde_run = 'test_run'
    parser = argparse.ArgumentParser(description='Spatial distribution')
    parser.add_argument('--mhd_run', action="store",
                        default=default_mhd_run, help='MHD run name')
    parser.add_argument('--mhd_run_dir', action="store",
                        default=default_mhd_run_dir, help='MHD run directory')
    parser.add_argument('--sde_run', action="store",
                        default=default_sde_run, help='SDE run name')
    parser.add_argument('--tframe', action="store", default='100', type=int,
                        help='Time frame')
    parser.add_argument('--ene_band', action="store", default='2', type=int,
                        help='energy band')
    parser.add_argument('--multi_frames', action="store_true", default=False,
                        help='whether to analyze multiple frames')
    parser.add_argument('--time_loop', action="store_true", default=False,
                        help='whether to use a time loop to analyze multiple frames')
    parser.add_argument('--multi_bands', action="store_true", default=False,
                        help='whether to plot multiple energy bands')
    parser.add_argument('--tstart', action="store", default='0', type=int,
                        help='starting time frame')
    parser.add_argument('--tend', action="store", default='10', type=int,
                        help='ending time frame')
    parser.add_argument('--dist_2d', action="store_true", default=False,
                        help='whether to plot 2D spatial distributions')
    return parser.parse_args()


def analysis_single_frames(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    if args.dist_2d:
        plot_dist_2d(plot_config, mhd_config, show_plot=True)


def process_input(plot_config, mhd_config, args, tframe):
    """process one time frame"""
    plot_config["tframe"] = tframe
    pass


def analysis_multi_frames(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    tframes = range(plot_config["tmin"], plot_config["tmax"] + 1)
    if args.time_loop:
        for tframe in tframes:
            print("Time frame: %d" % tframe)
            plot_config["tframe"] = tframe
            if args.dist_2d:
                plot_dist_2d(plot_config, mhd_config, show_plot=False)
    else:
        ncores = multiprocessing.cpu_count()
        ncores = 16
        Parallel(n_jobs=ncores)(delayed(process_input)(plot_config, mhd_config,
                                                       args, tframe)
                                for tframe in tframes)


def main():
    """business logic for when running this module as the primary one!"""
    args = get_cmd_args()
    with open('config/spectrum_config_10ta.json', 'r') as file_handler:
        config = json.load(file_handler)
    mhd_config = load_mhd_config(args.mhd_run_dir)
    plot_config = {}
    run_name = args.mhd_run + "/" + args.sde_run
    sde_run_config = config[run_name]
    nreduce = sde_run_config["nreduce"]
    plot_config["nx"] = mhd_config["nx"] // nreduce
    plot_config["ny"] = mhd_config["ny"] // nreduce
    plot_config["nz"] = mhd_config["nz"] // nreduce
    plot_config["nreduce"] = nreduce
    plot_config["tmax"] = sde_run_config["tmax"]
    plot_config["tmin"] = sde_run_config["tmin"]
    if "e0" in sde_run_config:
        plot_config["e0"] = sde_run_config["e0"]
    else:
        plot_config["e0"] = 1.0
    if "high_bands" in sde_run_config:
        plot_config["high_bands"] = sde_run_config["high_bands"]
    else:
        plot_config["high_bands"] = [2, 5]
    plot_config["eband"] = args.ene_band
    plot_config["run_name"] = run_name
    plot_config["tframe"] = args.tframe
    plot_config["mhd_run_dir"] = args.mhd_run_dir
    if args.multi_frames:
        analysis_multi_frames(plot_config, mhd_config, args)
    else:
        analysis_single_frames(plot_config, mhd_config, args)


if __name__ == "__main__":
    main()

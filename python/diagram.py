#!/usr/bin/env python3
"""
module for diagrams for different presentations
"""
from __future__ import print_function

import argparse
import itertools
import json
import multiprocessing

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import palettable
from joblib import Parallel, delayed
from matplotlib.colors import LogNorm

from sde_util import load_mhd_config, mkdir_p

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = \
[r"\usepackage{amsmath, bm}",
 r"\DeclareMathAlphabet{\mathsfit}{\encodingdefault}{\sfdefault}{m}{sl}",
 r"\SetMathAlphabet{\mathsfit}{bold}{\encodingdefault}{\sfdefault}{bx}{sl}",
 r"\newcommand{\tensorsym}[1]{\bm{\mathsfit{#1}}}"]


def figure_text(eband):
    """get figure text for different energy band

    Args:
        eband: energy band
    """
    if eband == 0:
        eneh = "{%0.0f}" % (2**eband * 0.75)**2
        fig_text = r'$\varepsilon <' + eneh + r'$keV'
    elif eband == 4:
        enel = "{%0.0f}" % (2**(eband - 2) * 1.5)**2
        fig_text = r'$\varepsilon >' + enel + r'$keV'
    else:
        enel = "{%0.0f}" % (2**(eband - 1) * 0.75)**2
        eneh = "{%0.0f}" % (2**(eband - 1) * 1.5)**2
        fig_text = (r'$' + enel + r'$keV' + r'$'
                    r'< \varepsilon <' +
                    eneh + r'$keV')

    return fig_text


def plot_box(center, lenx, leny, ax, color):
    """Plot a box in figure
    """
    xl = center[0] - lenx / 2
    xr = center[0] + lenx / 2
    yb = center[1] - leny / 2
    yt = center[1] + leny / 2
    xbox = [xl, xr, xr, xl, xl]
    ybox = [yb, yb, yt, yt, yb]
    ax.plot(xbox, ybox, color=color, linewidth=2)


def rhessi18(plot_config, mhd_config, show_plot=True):
    """Plot spatial distributions for multiple energy band

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    # mpl.rc('text', usetex=True)
    fig = plt.figure(figsize=[15, 15])
    rect0 = [0.07, 0.53, 0.41, 0.42]
    rect = np.copy(rect0)
    hgap, vgap = 0.04, 0.05
    run_name = plot_config["run_name"]
    tframe_str = str(plot_config["tframe"]).zfill(4)
    nx, = plot_config["nx"]
    ny, = plot_config["ny"]
    fname = '../data/' + run_name + '/fxy-' + tframe_str + '_sum.dat'
    fdata = np.fromfile(fname)
    nx, = plot_config["nx"]
    ny, = plot_config["ny"]
    print("data size: %d %d" % (nx, ny))
    dists = fdata.reshape([fdata.shape[0] // (nx * ny), -1])
    print("Total number of particles: %d" % np.sum(dists))

    axs = [None] * 4
    ims = [None] * 4
    vmin = plot_config['nmin']
    vmax = plot_config['nmax']
    sizes = [mhd_config["xmin"][0], mhd_config["xmax"][0],
             mhd_config["ymin"][0], mhd_config["ymax"][0]]
    lnorm = 75.0  # in Mm
    sizes = np.asarray(sizes) * lnorm

    # High-energy particle distributions
    for i in range(2):
        rect[0] = rect0[0] + (hgap + rect0[2]) * i
        eband = i + 3
        fband = dists[eband, :]
        fnorm = 4**(eband - 3)
        fband = np.reshape(fband, (ny, nx)) * fnorm
        fband += vmin * 0.01
        print("min, max and mean of the data: %f %f %f %f" %
              (np.min(fband), np.max(fband),
               np.mean(fband), np.std(fband)))
        axs[i] = fig.add_axes(rect)
        ims[i] = axs[i].imshow(fband, cmap=plt.cm.plasma,
                               aspect='auto', origin='lower', extent=sizes,
                               norm=LogNorm(vmin=vmin, vmax=vmax),
                               interpolation='bicubic')
        axs[i].tick_params(labelsize=16)
        if fnorm > 1:
            fig_text = r'$' + str(fnorm) + r'$' + r'$n($' + figure_text(eband) + r'$)$'
        else:
            fig_text = r'$n($' + figure_text(eband) + r'$)$'
        color = 'k' if eband == 0 else 'w'
        axs[i].text(0.02, 0.92, fig_text, color=color, fontsize=20,
                    bbox=dict(facecolor='none', alpha=1.0,
                              edgecolor='none', pad=10.0),
                    horizontalalignment='left', verticalalignment='bottom',
                    transform=axs[i].transAxes)

    rect_cbar = np.copy(rect0)
    rect_cbar[0] = rect0[0] + rect0[2] * 2 + 0.05
    rect_cbar[2] = 0.015
    cbar_ax = fig.add_axes(rect_cbar)
    cbar1 = fig.colorbar(ims[0], cax=cbar_ax, extend='both')
    cbar1.ax.tick_params(labelsize=16)

    # Current density
    mhd_run_dir = plot_config["mhd_run_dir"]
    fpath = mhd_run_dir + 'bin_data/'
    fname = fpath + 'mhd_data_' + tframe_str
    nx = mhd_config["nx"][0] + 4
    ny = mhd_config["ny"][0] + 4
    dxm = mhd_config["dx"][0]
    dym = mhd_config["dy"][0]
    mhd_fields = np.fromfile(fname, dtype=np.float32)
    mhd_fields = mhd_fields.reshape((ny, nx, 8))
    bx = mhd_fields[:, :, 4]
    by = mhd_fields[:, :, 5]
    jz = np.gradient(by, dxm, axis=1) - np.gradient(bx, dym, axis=0)
    rect[1] = rect0[1] - rect0[3] - vgap
    rect[0] = rect0[0]
    axs[2] = fig.add_axes(rect)
    ims[2] = axs[2].imshow(jz, cmap=plt.cm.seismic,
                           aspect='auto', origin='lower', extent=sizes,
                           vmin=-500, vmax=500,
                           interpolation='none')
    axs[2].tick_params(labelsize=16)
    axs[2].text(0.02, 0.92, r'$j_z/j_0$', color='k', fontsize=20,
                bbox=dict(facecolor='none', alpha=1.0,
                          edgecolor='none', pad=10.0),
                horizontalalignment='left', verticalalignment='bottom',
                transform=axs[2].transAxes)

    nyr = 4
    center = np.zeros(2)
    nx_mhd = mhd_config['nx'][0]
    ny_mhd = mhd_config['ny'][0]
    dx_mhd = mhd_config['dx'][0]
    dy_mhd = mhd_config['dy'][0]
    COLORS = palettable.tableau.Tableau_10.mpl_colors
    lenx = nx_mhd * dx_mhd * 0.5 * lnorm
    for iy in range(nyr):
        center[0] = 0
        center[1] = (ny_mhd * (iy + 0.5) / nyr) * dy_mhd * lnorm
        leny = (ny_mhd/nyr) * dy_mhd * lnorm - 1
        plot_box(center, lenx, leny, axs[2], COLORS[iy])

    rect_cbar = np.copy(rect)
    rect_cbar[0] = rect[0] + 0.02
    rect_cbar[2] = 0.015
    rect_cbar[1] += rect_cbar[3] * 0.25
    rect_cbar[3] *= 0.5
    cbar_ax = fig.add_axes(rect_cbar)
    cbar2 = fig.colorbar(ims[2], cax=cbar_ax, extend='both')
    cbar2.ax.tick_params(labelsize=16)

    # Local spectrum
    with open('config/spectrum_config_10ta.json', 'r') as file_handler:
        config = json.load(file_handler)
    sde_run_config = config[run_name]
    nreduce = sde_run_config["nreduce"]
    nx = mhd_config["nx"][0] // nreduce
    ny = mhd_config["ny"][0] // nreduce

    fname = '../data/' + run_name + '/fp_local_' + tframe_str + '_sum.dat'
    fdata = np.fromfile(fname)
    npp = fdata.shape[0] // (nx * ny)
    fdata = fdata.reshape([ny, nx, npp])
    print("Total number of particles: %d" % np.sum(fdata))
    print("data size: %d %d %d (C-order)" % fdata.shape)
    xs, xe = nx//4, nx*3//4
    fdata_bin = fdata[:, xs:xe, :].reshape(nyr, -1, nx, npp)
    fdata_r = np.sum(np.sum(fdata_bin, axis=1), axis=1)

    p0 = 0.1  # initial particle momentum
    pmom = np.logspace(-2, 1, npp + 1) / p0
    pmom_mid = (pmom[:-1] + pmom[1:]) * 0.5
    dpmom = np.diff(pmom)
    elog = pmom**2
    elog_mid = (elog[:-1] + elog[1:]) * 0.5

    rect[0] = rect0[0] + rect0[2] + hgap
    axs[3] = fig.add_axes(rect)
    axs[3].tick_params(labelsize=16)
    for iy in range(nyr):
        fmom = fdata_r[iy] * pmom_mid / dpmom  # f(p)p^3
        fene = fmom / pmom_mid*2  # f(p)p or f(e)
        fene[fene == 0] = np.nan
        axs[3].loglog(elog_mid, fene, linewidth=2, color=COLORS[iy])
    axs[3].set_xlim(1E-1, 1E4)
    axs[3].set_ylim(1E-2, 5E7)

    axs[1].tick_params(axis='y', labelleft=False)
    axs[0].set_ylabel(r'$y/$Mm', fontsize=20)
    axs[1].set_xlabel(r'$x/$Mm', fontsize=20)
    axs[2].set_xlabel(r'$x/$Mm', fontsize=20)
    axs[2].set_ylabel(r'$y/$Mm', fontsize=20)
    axs[3].set_xlabel(r'$\varepsilon/$keV', fontsize=20)

    tva = mhd_config["dt_out"][0] * plot_config["tframe"] / mhd_config["ly"][0]
    title = r'$t = ' + "{:10.1f}".format(tva) + r'\tau_A$'
    title += ' (frame: %d)' % plot_config["tframe"]
    fig.suptitle(title, fontsize=24)

    fdir = '../img/nrho_absj_spectra/' + run_name + '/'
    mkdir_p(fdir)
    fname = fdir + 'nrho_absj_specta_' + tframe_str + '.jpg'
    fig.savefig(fname)

    if show_plot:
        plt.show()
    else:
        plt.close()


def get_cmd_args():
    """Get command line arguments
    """
    default_mhd_run = 'S1E5_beta01_bg00'
    default_mhd_run_dir = ('/net/scratch3/xiaocanli/mhd/guide_field_scaling/' +
                           default_mhd_run + '/')
    default_sde_run = 'p000_b000_00001_100'
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
    parser.add_argument('--tstart', action="store", default='0', type=int,
                        help='starting time frame')
    parser.add_argument('--tend', action="store", default='10', type=int,
                        help='ending time frame')
    parser.add_argument('--nmin', action="store", default='0.1', type=float,
                        help='minimum density')
    parser.add_argument('--nmax', action="store", default='10.0', type=float,
                        help='maximum density')
    parser.add_argument('--rhessi18', action="store_true", default=False,
                        help='whether to plot diagram for RHESSI18')
    return parser.parse_args()


def analysis_single_frames(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    if args.rhessi18:
        rhessi18(plot_config, mhd_config, show_plot=True)


def process_input(plot_config, mhd_config, args, tframe):
    """process one time frame"""
    plot_config["tframe"] = tframe
    if args.rhessi18:
        rhessi18(plot_config, mhd_config, show_plot=False)


def analysis_multi_frames(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    tframes = range(plot_config["tmin"], plot_config["tmax"] + 1)
    if args.time_loop:
        for tframe in tframes:
            print("Time frame: %d" % tframe)
            plot_config["tframe"] = tframe
            if args.rhessi18:
                rhessi18(plot_config, mhd_config, show_plot=False)
    else:
        ncores = multiprocessing.cpu_count()
        ncores = 8
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
    plot_config["tmax"] = sde_run_config["tmax"]
    plot_config["tmin"] = sde_run_config["tmin"]
    plot_config["eband"] = args.ene_band
    plot_config["run_name"] = run_name
    plot_config["tframe"] = args.tframe
    plot_config["nmin"] = args.nmin
    plot_config["nmax"] = args.nmax
    plot_config["mhd_run"] = args.mhd_run
    plot_config["mhd_run_dir"] = args.mhd_run_dir
    if args.multi_frames:
        analysis_multi_frames(plot_config, mhd_config, args)
    else:
        analysis_single_frames(plot_config, mhd_config, args)


if __name__ == "__main__":
    main()

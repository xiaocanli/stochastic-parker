#!/usr/bin/env python3
"""
Analysis procedures for stochastic integration
"""
from __future__ import print_function

import argparse
import itertools
import json
import multiprocessing

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from matplotlib.colors import LogNorm

from util import load_mhd_config, mkdir_p

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]


def figure_text(eband):
    """get figure text for different energy band

    Args:
        eband: energy band
    """
    if eband == 0:
        eneh = "{%0.1f}" % (2**eband * 0.75)**2
        fig_text = r'$\varepsilon/\varepsilon_0 <' + eneh + '$'
    elif eband == 4:
        enel = "{%0.1f}" % (2**(eband - 2) * 1.5)**2
        fig_text = r'$\varepsilon/\varepsilon_0 >' + enel + '$'
    else:
        enel = "{%0.1f}" % (2**(eband - 1) * 0.75)**2
        eneh = "{%0.1f}" % (2**(eband - 1) * 1.5)**2
        fig_text = r'$' + enel + r'\leq \varepsilon/\varepsilon_0 <' + eneh + '$'

    return fig_text


def spatial_distribution_multi(plot_config, mhd_config, show_plot=True):
    """Plot spatial distributions for multiple energy band

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    fig = plt.figure(figsize=[15, 15])
    rect0 = [0.07, 0.53, 0.42, 0.42]
    rect = np.copy(rect0)
    hgap, vgap = 0.02, 0.02
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
    vmin, vmax = 1, 2E1
    sizes = [mhd_config["xmin"][0], mhd_config["xmax"][0],
             mhd_config["ymin"][0], mhd_config["ymax"][0]]
    for i, j in itertools.product(range(2), range(2)):
        rect[0] = rect0[0] + (hgap + rect0[2]) * i
        rect[1] = rect0[1] - (vgap + rect0[3]) * j
        eband = j * 2 + i + 1
        fband = dists[eband, :]
        fband = np.reshape(fband, (ny, nx))
        print("min, max and mean of the data: %f %f %f %f" %
              (np.min(fband), np.max(fband),
               np.mean(fband), np.std(fband)))
        axs[eband-1] = fig.add_axes(rect)
        ims[eband-1] = axs[eband-1].imshow(fband, cmap=plt.cm.viridis,
                                           aspect='auto', origin='lower',
                                           extent=sizes,
                                           norm=LogNorm(vmin=vmin, vmax=vmax),
                                           interpolation='bicubic')
        axs[eband - 1].tick_params(labelsize=16)
        fig_text = figure_text(eband)
        axs[eband-1].text(0.02, 0.9, fig_text, color='k', fontsize=20,
                          bbox=dict(facecolor='none', alpha=1.0,
                                    edgecolor='none', pad=10.0),
                          horizontalalignment='left', verticalalignment='bottom',
                          transform=axs[eband-1].transAxes)

    axs[2].set_xlabel(r'$x$', fontsize=20)
    axs[3].set_xlabel(r'$x$', fontsize=20)
    axs[0].set_ylabel(r'$y$', fontsize=20)
    axs[2].set_ylabel(r'$y$', fontsize=20)

    axs[0].tick_params(axis='x', labelbottom='off')
    axs[1].tick_params(axis='x', labelbottom='off')
    axs[1].tick_params(axis='y', labelleft='off')
    axs[3].tick_params(axis='y', labelleft='off')

    rect[0] = rect0[0] + rect0[2] * 2 + 0.03
    rect[2] = 0.02
    rect[3] = rect0[3] * 2 + vgap
    cbar_ax = fig.add_axes(rect)
    cbar1 = fig.colorbar(ims[0], cax=cbar_ax)
    cbar1.ax.tick_params(labelsize=16)

    tva = mhd_config["dt_out"][0] * plot_config["tframe"] / mhd_config["ly"][0]
    title = r'$t = ' + "{:10.1f}".format(tva) + r'\tau_A$'
    title += ' (frame: %d)' % plot_config["tframe"]
    fig.suptitle(title, fontsize=24)

    fdir = '../img/nrho/' + run_name + '/'
    mkdir_p(fdir)
    fname = fdir + 'nrho_bands_' + tframe_str + '.jpg'
    fig.savefig(fname, dpi=200)

    if show_plot:
        plt.show()
    else:
        plt.close()


def spatial_distribution_band(plot_config, mhd_config, show_plot=True):
    """Plot particle spatial distribution for one energy band

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    tframe_str = str(plot_config["tframe"]).zfill(4)
    eband = plot_config["eband"]
    run_name = plot_config["run_name"]
    fname = '../data/' + run_name + '/fxy-' + tframe_str + '_sum.dat'
    fdata = np.fromfile(fname)
    nx, = plot_config["nx"]
    ny, = plot_config["ny"]
    print("data size: %d %d" % (nx, ny))
    dists = fdata.reshape([fdata.shape[0] // (nx * ny), -1])
    print("Total number of particles: %d" % np.sum(dists))

    fband = dists[eband, :]
    print("particle number in band %d: %d" % (eband, fband.size))
    fband = np.reshape(fband, (ny, nx))
    print("min, max, mean, and std: %f %f %f %f" %
          (np.min(fband), np.max(fband), np.mean(fband), np.std(fband)))

    fig = plt.figure(figsize=[12, 12])
    rect = [0.07, 0.1, 0.86, 0.86]
    ax = fig.add_axes(rect)
    vmin, vmax = 1, 2E1
    sizes = [mhd_config["xmin"][0], mhd_config["xmax"][0],
             mhd_config["ymin"][0], mhd_config["ymax"][0]]
    img = ax.imshow(fband, cmap=plt.cm.viridis, aspect='auto',
                    origin='lower', extent=sizes,
                    # vmin=vmin, vmax=vmax,
                    norm=LogNorm(vmin=vmin, vmax=vmax),
                    interpolation='bicubic')
    ax.tick_params(labelsize=16)
    fig_text = figure_text(eband)
    ax.text(0.02, 0.9, fig_text, color='k', fontsize=20,
            bbox=dict(facecolor='none', alpha=1.0, edgecolor='none', pad=10.0),
            horizontalalignment='left', verticalalignment='bottom',
            transform=ax.transAxes)

    ax.set_xlabel(r'$x$', fontsize=20)
    ax.set_ylabel(r'$y$', fontsize=20)

    rect[0] += rect[2] + 0.01
    rect[2] = 0.02
    cbar_ax = fig.add_axes(rect)
    cbar1 = fig.colorbar(img, cax=cbar_ax)
    cbar1.ax.tick_params(labelsize=16)

    fdir = '../img/nrho/' + run_name + '/'
    mkdir_p(fdir)
    fname = fdir + 'nrho_' + tframe_str + '_' + str(eband) + '.jpg'
    fig.savefig(fname, dpi=200)

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
    parser.add_argument('--multi_bands', action="store_true", default=False,
                        help='whether to plot multiple energy bands')
    parser.add_argument('--tstart', action="store", default='0', type=int,
                        help='starting time frame')
    parser.add_argument('--tend', action="store", default='10', type=int,
                        help='ending time frame')
    return parser.parse_args()


def analysis_single_frames(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    if args.multi_bands:
        spatial_distribution_multi(plot_config, mhd_config, show_plot=True)
    else:
        spatial_distribution_band(plot_config, mhd_config, show_plot=True)


def process_input(plot_config, mhd_config, args, tframe):
    """process one time frame"""
    plot_config["tframe"] = tframe
    if args.multi_bands:
        spatial_distribution_multi(plot_config, mhd_config, show_plot=False)
    else:
        spatial_distribution_band(plot_config, mhd_config, show_plot=False)


def analysis_multi_frames(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    tframes = range(plot_config["tmin"], plot_config["tmax"] + 1)
    ncores = multiprocessing.cpu_count()
    Parallel(n_jobs=ncores)(delayed(process_input)(plot_config, mhd_config,
                                                   args, tframe)
                            for tframe in tframes)


def main():
    """business logic for when running this module as the primary one!"""
    args = get_cmd_args()
    with open('config/spectrum_config.json', 'r') as file_handler:
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
    if args.multi_frames:
        analysis_multi_frames(plot_config, mhd_config, args)
    else:
        analysis_single_frames(plot_config, mhd_config, args)


if __name__ == "__main__":
    main()

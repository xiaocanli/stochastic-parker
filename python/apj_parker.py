#!/usr/bin/env python3
"""
Analysis procedures for ApJ paper
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

import particle_trajectory as traj
from sde_util import load_mhd_config, mkdir_p
from spatial_distribution import figure_text

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]


def read_az(mhd_run_dir, nx, ny, tframe):
    """Read the out-of-plane current density

    Args:
        mhd_run_dir: MHD run directory
        nx, ny: data dimensions
        tframe: time frame
    """
    az = np.zeros((nx, ny))
    fname = mhd_run_dir + 'data/Ay.gda'
    az = np.memmap(fname, dtype='float32', mode='r',
                   offset=nx*ny*tframe*4, shape=(nx, ny))
    return az


def spatial_distribution_band(plot_config, mhd_config, show_plot=True):
    """Plot particle spatial distribution for one energy band

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    fig = plt.figure(figsize=[12, 8])
    rect = [0.07, 0.12, 0.25, 0.8]
    axs = []
    mpl.rcParams['contour.negative_linestyle'] = 'solid'
    for tframe in plot_config["tframes"]:
        tframe_str = str(tframe).zfill(4)
        eband = plot_config["eband"]
        run_name = plot_config["run_name"]
        fname = '../data/' + run_name + '/fxy-' + tframe_str + '_sum.dat'
        fdata = np.fromfile(fname)
        nx, = plot_config["nx"]
        ny, = plot_config["ny"]
        nx_mhd = mhd_config["nx"][0]
        ny_mhd = mhd_config["ny"][0]
        az = read_az(plot_config["mhd_run_dir"], nx_mhd, ny_mhd, tframe)
        xgrid = np.linspace(mhd_config["xmin"], mhd_config["xmax"],
                            mhd_config["nx"])
        ygrid = np.linspace(mhd_config["ymin"], mhd_config["ymax"],
                            mhd_config["ny"])
        print("data size: %d %d" % (nx, ny))
        dists = fdata.reshape([fdata.shape[0] // (nx * ny), -1])
        print("Total number of particles: %d" % np.sum(dists))

        fband = dists[eband, :]
        print("particle number in band %d: %d" % (eband, fband.size))
        fband = np.reshape(fband, (ny, nx))
        print("min, max, mean, and std: %f %f %f %f" %
              (np.min(fband), np.max(fband), np.mean(fband), np.std(fband)))

        ax = fig.add_axes(rect)
        axs.append(ax)
        vmin, vmax = 0.1, 1E1
        lx = mhd_config["xmax"][0] - mhd_config["xmin"][0]
        sizes = [mhd_config["xmin"][0] + lx*0.5, mhd_config["xmax"][0],
                 mhd_config["ymin"][0], mhd_config["ymax"][0]]
        img = ax.imshow(fband[:, nx//2:], cmap=plt.cm.viridis, aspect='auto',
                        origin='lower', extent=sizes,
                        # vmin=vmin, vmax=vmax,
                        norm=LogNorm(vmin=vmin, vmax=vmax),
                        interpolation='bicubic')
        level_az = np.linspace(np.min(az[nx_mhd//2:, :]),
                               np.max(az[nx_mhd//2:, :]), 15)
        ax.contour(xgrid[nx_mhd//2:], ygrid,
                   az[nx_mhd//2:, :].T, colors='black',
                   linewidths=0.5, levels=level_az)
        ax.tick_params(labelsize=16)
        rect[0] += rect[2] + 0.03
        ax.set_xlabel(r'$x$', fontsize=20)
        tva = mhd_config["dt_out"][0] * tframe / mhd_config["ly"][0]
        title = r'$t = ' + "{:10.1f}".format(tva) + r'\tau_A$'
        plt.title(title, fontsize=24)

    axs[0].set_ylabel(r'$y$', fontsize=20)
    axs[1].tick_params(axis='y', labelleft='off')
    axs[2].tick_params(axis='y', labelleft='off')

    rect[0] -= 0.01
    rect[2] = 0.02
    cbar_ax = fig.add_axes(rect)
    cbar1 = fig.colorbar(img, cax=cbar_ax)
    cbar1.ax.tick_params(labelsize=16)
    fig_text = figure_text(eband)
    title = r"$n($" + fig_text + r"$)$"
    cbar1.ax.set_ylabel(title, fontsize=24)

    fdir = '../img/apj_parker/'
    mkdir_p(fdir)
    tstr = ""
    for tframe in plot_config["tframes"]:
        tstr += str(tframe) + "_"
    fname = fdir + 'nrho_' + tstr + str(eband) + '.jpg'
    fig.savefig(fname, dpi=200)

    if show_plot:
        plt.show()
    else:
        plt.close()


def ptl_traj(mhd_run_info, sde_run, mhd_config, config):
    """Plot particle trajectory

    Args:
        mhd_run_info: information of the MHD run
        sde_run: SDE run name
        mhd_config: MHD configuration
        config: additional configuration
    """
    fpath = "../data/traj_data/"
    opath = "../img/ptl_traj/" + mhd_run_info["run_name"] + '/' + sde_run + '/'
    mkdir_p(opath)
    fname = (fpath + mhd_run_info["run_name"] + '/' + sde_run +
             "/tracked_particle_points_" +
             str(config["mpi_rank"]).zfill(4) + ".dat")
    nptl = traj.get_nptl(fname)
    ptl = traj.read_ptl_traj(fname, nptl, config["ptl_index"])
    ptl = np.array([list(item) for item in ptl])
    ptl = traj.boundary_cross(ptl)
    xptl = ptl[:-1, 0]
    yptl = ptl[:-1, 1]
    pptl = ptl[:-1, 2]
    tptl = ptl[:-1, 4]
    ptl_range = [np.min(xptl), np.max(xptl), np.min(yptl), np.max(yptl)]

    fig = plt.figure(figsize=[20, 10])
    rect = [0.05, 0.12, 0.15, 0.8]
    hgap = 0.06

    # tframe = config["tframes"][0]
    # mhd_time = mhd_config["dt_out"][0] * tframe
    # # find the closed time point
    # tindex, _ = traj.find_nearest(tptl, mhd_time)
    # mhd_box, rho = traj.read_mhd_data(mhd_run_info, tframe)
    # prange, mhd_data = traj.adjust_mhd_data(mhd_box, rho, ptl_range)
    tindices = []

    for (i, tframe) in enumerate(config["tframes"]):
        mhd_time = mhd_config["dt_out"][0] * tframe
        tva = mhd_time / mhd_config["ly"][0] # in Alfven crossing time
        title = r'$t = ' + "{:10.1f}".format(tva) + r'\tau_A$'

        # find the closed time point
        tindex, _ = traj.find_nearest(tptl, mhd_time)
        tindices.append(tindex)
        mhd_box, rho = traj.read_mhd_data(mhd_run_info, tframe)
        prange, mhd_data = traj.adjust_mhd_data(mhd_box, rho, ptl_range)

        ax = fig.add_axes(rect)
        if i > 0:
            ax.tick_params(axis='y', labelleft='off')
        ylabeling = False if i > 0 else True
        img = ax.imshow(mhd_data, extent=prange,
                        vmin=0.5, vmax=3.0,
                        cmap=plt.cm.viridis,
                        aspect='auto', origin='lower')
        traj.plot_particle_2d(('x', xptl), ('y', yptl), ax, tindex, ylabeling)
        ax.set_xlim(prange[:2])
        ax.set_ylim(prange[2:])
        plt.title(title, fontsize=24)
        rect[0] += rect[2] + 0.015

    rect[2] = 0.015
    cbar_ax = fig.add_axes(rect)
    cbar = fig.colorbar(img, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=16)
    cbar.ax.set_xlabel(r'$\rho$', fontsize=24)
    cbar.ax.xaxis.set_label_position('top')

    tindices = np.asarray(tindices)
    rect[0] += rect[2] + 0.07
    rect[2] = 0.17
    ax1 = fig.add_axes(rect)
    traj.plot_particle_2d(('x', xptl), ('p', pptl), ax1, tindices[0])
    for tindex in tindices[1:]:
        ax1.plot(xptl[tindex], pptl[tindex], marker='o',
                 markersize=10, color="red")

    rect[0] += rect[2] + 0.015
    ax2 = fig.add_axes(rect)
    ax2.tick_params(axis='y', labelleft='off')
    tptl /= mhd_config["ly"][0] # in Alfven crossing time
    traj.plot_particle_2d(('t', tptl), ('p', pptl), ax2,
                          tindices[0], ylabeling=False)
    for tindex in tindices[1:]:
        ax2.plot(tptl[tindex], pptl[tindex], marker='o',
                 markersize=10, color="red")


    out_dir = '../img/apj_parker/'
    mkdir_p(out_dir)
    fig_name = (out_dir + "ptl_" + str(config["mpi_rank"]) +
                "_"  + str(config["ptl_index"]) + ".jpg")
    fig.savefig(fig_name, dpi=200)

    plt.show()
    # plt.close()


def get_cmd_args():
    """Get command line arguments
    """
    default_mhd_run = 'S1E5_beta01_bg00'
    default_mhd_run_dir = ('/net/scratch3/xiaocanli/mhd/guide_field_scaling/' +
                           default_mhd_run + '/')
    default_sde_run = 'p133_b000_0008_001_l'
    parser = argparse.ArgumentParser(description='Spatial distribution')
    parser.add_argument('--mhd_run', action="store",
                        default=default_mhd_run, help='MHD run name')
    parser.add_argument('--mhd_run_dir', action="store",
                        default=default_mhd_run_dir, help='MHD run directory')
    parser.add_argument('--mhd_run_type', action="store", default="reconnection",
                        help='MHD run type')
    parser.add_argument('--mhd_code', action="store", default="Athena",
                        help='MHD code')
    parser.add_argument('--config_name', action="store",
                        default="athinput.reconnection",
                        help='MHD configuration filename')
    parser.add_argument('--sde_run', action="store",
                        default=default_sde_run, help='SDE run name')
    parser.add_argument('--tframe', action="store", default='100', type=int,
                        help='Time frame')
    parser.add_argument('--ene_band', action="store", default='2', type=int,
                        help='energy band')
    parser.add_argument('--multi_frames', action="store_true", default=False,
                        help='whether to analyze multiple frames')
    parser.add_argument('--tstart', action="store", default='0', type=int,
                        help='starting time frame')
    parser.add_argument('--tend', action="store", default='10', type=int,
                        help='ending time frame')
    parser.add_argument('--spatial_dist', action="store_true", default=False,
                        help='whether to plot spatial distribution')
    parser.add_argument('--ptl_traj', action="store_true", default=False,
                        help='whether to plot particle trajectories')
    return parser.parse_args()


def analysis_single_frame(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    if args.spatial_dist:
        spatial_distribution_band(plot_config, mhd_config, show_plot=True)
    elif args.ptl_traj:
        mhd_run_info = traj.get_mhd_info(args)
        config = {}
        config["mpi_rank"] = 1
        config["ptl_index"] = 99
        config["tframes"] = [150, 155, 186]
        ptl_traj(mhd_run_info, args.sde_run, mhd_config, config)


def process_input(plot_config, mhd_config, args, tframe):
    """process one time frame"""
    plot_config["tframe"] = tframe
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
    with open('config/spectrum_config_bg.json', 'r') as file_handler:
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
    plot_config["mhd_run_dir"] = args.mhd_run_dir
    tframes = [150, 175, 200]
    plot_config["tframes"] = tframes
    if args.multi_frames:
        analysis_multi_frames(plot_config, mhd_config, args)
    else:
        analysis_single_frame(plot_config, mhd_config, args)


if __name__ == "__main__":
    main()

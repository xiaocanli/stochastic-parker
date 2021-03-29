#!/usr/bin/env python3
"""
Analysis procedures for plotting particle trajectory
"""
from __future__ import print_function

import argparse
import math
import multiprocessing
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed

from sde_util import mkdir_p, load_mhd_config, find_nearest

sys.path.insert(0, '/global/homes/x/xiaocan/Git/mhd_analysis_sli')
import mhd_data

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]


def get_nptl(fname):
    """get the number of tracked particles

    Args:
        fname: the filename of the tracer file
    """
    with open(fname, 'r') as file_h:
        tmp = np.memmap(file_h, dtype='int32', mode='r',
                        offset=0, shape=(1), order='F')
        nptl = tmp[0]
    return nptl


def read_ptl_traj(fname, nptl, ptl_index):
    """Read particle trajectories data

    Args:
        fname: the filename of the tracer file
        nptl: number of particles
        ptl_index: particle index
    """
    with open(fname, 'r') as file_h:
        offset = 4
        nsteps_tracked = np.memmap(file_h, dtype='int32', mode='r',
                                   offset=offset, shape=(nptl), order='F')
        ptl_type = np.dtype([('x', np.float64),
                             ('y', np.float64),
                             ('p', np.float64),
                             ('weight', np.float64),
                             ('t', np.float64),
                             ('dt', np.float64),
                             ('split_time', np.int32),
                             ('count_flag', np.int32),
                             ('tag', np.int32),
                             ('nsteps_tracking', np.int32)])
        ptl_size = ptl_type.itemsize
        offset += nptl * 4 + np.sum(nsteps_tracked[:ptl_index]) * ptl_size
        ptl = np.memmap(file_h, dtype=ptl_type, mode='r', offset=offset,
                        shape=(nsteps_tracked[ptl_index]), order='F')
    return ptl


def boundary_cross(ptl):
    """Shift particle if it crosses boundaries

    Args:
        ptl: particle information
    """
    deltax = np.diff(ptl[:-1, 0])
    deltay = np.diff(ptl[:-1, 1])
    if np.any(np.abs(deltax) > 0.5):
        xc_index = np.array(np.where(np.abs(deltax) > 0.5)) + 1
        sizex = xc_index.shape
        if sizex[1] > 1:
            xc_index = np.squeeze(xc_index)
            for xci in xc_index:
                shift = -2 if deltax[xci-1] > 0 else 2
                ptl[xci:, 0] += shift
        else:
            xci = xc_index[0, 0]
            shift = -2 if deltax[xci-1] > 0 else 2
            ptl[xci:, 0] += shift
    if np.any(np.abs(deltay) > 0.5):
        yc_index = np.array(np.where(np.abs(deltay) > 0.5)) + 1
        sizey = yc_index.shape
        if sizey[1] > 1:
            yc_index = np.squeeze(yc_index)
            for yci in yc_index:
                shift = -2 if deltay[yci-1] > 0 else 2
                ptl[yci:, 1] += shift
        else:
            yci = yc_index[0, 0]
            shift = -2 if deltay[yci-1] > 0 else 2
            ptl[yci:, 1] += shift

    return ptl


def read_mhd_data(mhd_run_info, tframe):
    """Read MHD data

    Args:
        mhd_run_info: information of the MHD run
        tframe: time frame
    """
    xmesh, ymesh, data = mhd_data.read_fields_data(mhd_run_info, tframe)
    ny, nx = xmesh.shape
    mhd_box = [xmesh[0, 0], xmesh[0, -1], ymesh[0, 0], ymesh[-1, 0], nx, ny]
    rho = data[:, :, 0].T
    return (mhd_box, rho)


def adjust_mhd_data(mhd_box, mhd_field, ptl_range):
    """Adjust MHD data based on particle range

    Args:
        mhd_box: MHD box range
        mhd_field: MHD simulation field
        ptl_range: particle range
    """
    lenx = mhd_box[1] - mhd_box[0]
    leny = mhd_box[3] - mhd_box[2]
    dxm = lenx / mhd_box[4]
    dym = leny / mhd_box[5]
    prange = np.copy(ptl_range)
    prange[0] = math.floor(ptl_range[0]) + 0.1
    prange[1] = math.ceil(ptl_range[1]) - 0.1
    ycenter = (prange[2] + prange[3]) * 0.5
    prange[2] = ycenter - leny * 0.5
    prange[3] = ycenter + leny * 0.5
    ixl = int(math.floor((prange[0] - mhd_box[0]) / dxm))
    ixr = int(math.ceil((prange[1] - mhd_box[0]) / dxm))
    iyb = int(math.floor((prange[2] - mhd_box[2]) / dym))
    iyt = int(math.ceil((prange[3] - mhd_box[2]) / dym))
    fdata = np.zeros((iyt-iyb+1, ixr-ixl+1))
    mdata = np.copy(mhd_field)
    mdata = np.roll(mdata, -iyb, axis=0)
    mdata = np.roll(mdata, -ixl, axis=1)
    fdata = mdata[0:iyt-iyb+1, 0:ixr-ixl+1]
    return (prange, fdata)


def plot_particle_2d(xdata, ydata, ax, tindex, ylabeling=True):
    """Plot particle trajectory in a 2D plane

    Args:
        xdata: data along the x-direction (pair)
        ydata: data along the y-direction (pair)
        ax: plot axis handler
        tindex: time index corresponding to current MHD time frame
    """
    ax.plot(xdata[1], ydata[1], color='k', linewidth=2)
    ax.plot(xdata[1][tindex], ydata[1][tindex],
            marker='o', markersize=10, color="red")
    ax.tick_params(labelsize=20)
    ax.set_xlabel(r'$' + xdata[0] + '$', fontsize=24)
    if ylabeling:
        ax.set_ylabel(r'$' + ydata[0] + '$', fontsize=20)


def plot_ptl_traj(mhd_run_info, sde_run, tframe, mhd_config, mpi_size):
    """Plot particle trajectories

    Args:
        mhd_run_info: information of the MHD run
        sde_run: SDE run name
        tframe: time frame
    """
    mhd_box, rho = read_mhd_data(mhd_run_info, tframe)

    mhd_time = mhd_config["dt_out"][0] * tframe
    tva = mhd_time / mhd_config["ly"][0] # in Alfven crossing time
    title = r'$t = ' + "{:10.1f}".format(tva) + r'\tau_A$'
    title += ' (frame: %d)' % tframe

    fpath = "../data/traj_data/"
    opath = "../img/ptl_traj/" + mhd_run_info["run_name"] + '/' + sde_run + '/'
    mkdir_p(opath)
    for mpi_rank in range(mpi_size):
        fname = (fpath + mhd_run_info["run_name"] + '/' + sde_run +
                 "/tracked_particle_points_" + str(mpi_rank).zfill(4) + ".dat")
        nptl = get_nptl(fname)
        odir1 = opath + "mpi_rank_" + str(mpi_rank) + "/"
        mkdir_p(odir1)
        for ptl_index in range(nptl-2, nptl):
            ptl = read_ptl_traj(fname, nptl, ptl_index)
            ptl = np.array([list(item) for item in ptl])
            ptl = boundary_cross(ptl)
            xptl = ptl[:-1, 0]
            yptl = ptl[:-1, 1]
            pptl = ptl[:-1, 2]
            tptl = ptl[:-1, 4]
            ptl_range = [np.min(xptl), np.max(xptl), np.min(yptl), np.max(yptl)]

            # find the closed time point
            tindex, _ = find_nearest(tptl, mhd_time)

            prange, mhd_data = adjust_mhd_data(mhd_box, rho, ptl_range)

            fig = plt.figure(figsize=[20, 10])
            rect = [0.05, 0.12, 0.2, 0.8]
            hgap = 0.06
            ax = fig.add_axes(rect)
            img = ax.imshow(mhd_data, extent=prange,
                            vmin=0.5, vmax=3.0,
                            cmap=plt.cm.viridis,
                            aspect='auto', origin='lower')
            plot_particle_2d(('x', xptl), ('y', yptl), ax, tindex)
            ax.set_xlim(prange[:2])
            ax.set_ylim(prange[2:])

            rect[0] += rect[2] + 0.01
            rect[2] = 0.015
            cbar_ax = fig.add_axes(rect)
            cbar = fig.colorbar(img, cax=cbar_ax)
            cbar.ax.tick_params(labelsize=16)
            cbar.ax.set_xlabel(r'$\rho$', fontsize=24)
            cbar.ax.xaxis.set_label_position('top')

            rect[0] += rect[2] + hgap
            rect[2] = 0.29
            ax1 = fig.add_axes(rect)
            plot_particle_2d(('x', xptl), ('p', pptl), ax1, tindex)

            rect[0] += rect[2] + hgap
            ax2 = fig.add_axes(rect)
            plot_particle_2d(('t', tptl), ('p', pptl), ax2, tindex)

            fig.suptitle(title, fontsize=24)

            out_dir = odir1 + "ptl_" + str(ptl_index) + "/"
            mkdir_p(out_dir)
            fig_name = out_dir + "tframe_" + str(tframe) + ".jpg"
            fig.savefig(fig_name, dpi=200)

            # plt.show()
            plt.close()


def get_cmd_args():
    """Get command line arguments """
    default_run_name = 'S1E5_beta01_bg00'
    default_run_dir = ('/net/scratch3/xiaocanli/mhd/guide_field_scaling/' +
                       default_run_name + '/')
    default_sde_run = 'p000_b000_00001_100'
    parser = argparse.ArgumentParser(description='Plot particle trajectory')
    parser.add_argument('--mhd_run_dir', action="store", default=default_run_dir,
                        help='MHD run directory')
    parser.add_argument('--mhd_run', action="store", default=default_run_name,
                        help='MHD run name')
    parser.add_argument('--mhd_run_type', action="store", default="reconnection",
                        help='MHD run type')
    parser.add_argument('--mhd_code', action="store", default="Athena",
                        help='MHD code')
    parser.add_argument('--config_name', action="store",
                        default="athinput.reconnection",
                        help='MHD configuration filename')
    parser.add_argument('--sde_run', action="store",
                        default=default_sde_run, help='SDE run name')
    parser.add_argument('--multi_frames', action="store_true", default=False,
                        help='whether to analyze multiple frames')
    parser.add_argument('--tstart', action="store", default='0', type=int,
                        help='starting time frame')
    parser.add_argument('--tend', action="store", default='10', type=int,
                        help='ending time frame')
    return parser.parse_args()


def get_mhd_info(args):
    """Get MHD run information
    """
    mhd_run_info = {}
    mhd_run_info["run_name"] = args.mhd_run
    mhd_run_info["run_dir"] = args.mhd_run_dir
    mhd_run_info["run_type"] = args.mhd_run_type
    mhd_run_info["mhd_code"] = args.mhd_code
    mhd_run_info["config_name"] = mhd_run_info["run_dir"] + args.config_name
    return mhd_run_info


def process_input(mhd_run_info, mhd_config, args, tframe):
    """process one time frame"""
    plot_ptl_traj(mhd_run_info, args.sde_run, tframe, mhd_config, 8)


def analysis_multi_frames(mhd_run_info, mhd_config, args):
    """Analysis for multiple time frames
    """
    tframes = range(args.tstart, args.tend + 1)
    ncores = multiprocessing.cpu_count()
    ncores = 4
    Parallel(n_jobs=ncores)(delayed(process_input)(mhd_run_info,
                                                   mhd_config,
                                                   args, tframe)
                            for tframe in tframes)


def main():
    """business logic for when running this module as the primary one!"""
    args = get_cmd_args()
    mhd_run_info = get_mhd_info(args)
    mhd_config = load_mhd_config(args.mhd_run_dir)
    if args.multi_frames:
        analysis_multi_frames(mhd_run_info, mhd_config, args)
    else:
        plot_ptl_traj(mhd_run_info, args.sde_run, 200, mhd_config, 8)


if __name__ == "__main__":
    main()

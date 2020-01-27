#!/usr/bin/env python3
"""
Analysis procedures for analyzing the compression
"""
from __future__ import print_function

import argparse
import json
import math
import multiprocessing
import sys

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from scipy import signal
from scipy.ndimage.filters import gaussian_filter, median_filter
from scipy.optimize import curve_fit

from sde_util import mkdir_p, load_mhd_config, find_nearest, div0

sys.path.insert(0, '/users/xiaocanli/Git/mhd_analysis_sli')
import mhd_data

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

def gauss(x, *p):
    """Gaussian function
    """
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))


def read_mhd_vel(mhd_run_info, tframe):
    """Read MHD velocity data

    Args:
        mhd_run_info: information of the MHD run
        tframe: time frame
    """
    xmesh, ymesh, data = mhd_data.read_fields_data(mhd_run_info, tframe)
    ny, nx = xmesh.shape
    mhd_box = [xmesh[0, 0], xmesh[0, -1], ymesh[0, 0], ymesh[-1, 0], nx, ny]
    vx = data[:, :, 2].T
    vy = data[:, :, 3].T
    vz = data[:, :, 4].T
    return (mhd_box, vx, vy, vz)


def compression_dist(mhd_run_info, tframe):
    """Calculate the distribution of the divergence of velocity field

    Args:
        mhd_run_info: information of the MHD run
        tframe: time frame
    """
    mhd_box, vx, vy, _ = read_mhd_vel(mhd_run_info, tframe)
    dx_mhd = (mhd_box[1] - mhd_box[0]) / mhd_box[4]
    dy_mhd = (mhd_box[3] - mhd_box[2]) / mhd_box[5]
    divv = (np.gradient(vx, axis=1) / dx_mhd +
            np.gradient(vy, axis=0) / dy_mhd)
    nbins = 200
    vbins = np.linspace(-500, 300, nbins + 1)
    fdivv, bins_edges = np.histogram(divv, bins=vbins)
    vbins = (vbins[:-1] + vbins[1:]) * 0.5
    fdist = np.zeros((2, nbins))
    fdist[0, :] = vbins
    fdist[1, :] = fdivv
    # fdir = "../data/divv_dist/" + mhd_run_info["run_name"] + "/"
    # mkdir_p(fdir)
    # fname = fdir + "divv_dist_" + str(tframe) + ".dat"
    # fdist.tofile(fname)
    plt.semilogy(vbins, fdivv)
    plt.show()


def rebin(a, shape):
    """https://stackoverflow.com/a/8090605/2561161
    """
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)


def reduced_compression(mhd_run_info, tframe, nreduce):
    """Compression with reduced size

    Args:
        mhd_run_info: information of the MHD run
        tframe: time frame
        nreduce: reduce factor
    """
    mhd_box, vx, vy, _ = read_mhd_vel(mhd_run_info, tframe)
    dx_mhd = (mhd_box[1] - mhd_box[0]) / mhd_box[4]
    dy_mhd = (mhd_box[3] - mhd_box[2]) / mhd_box[5]
    divv = (np.gradient(vx, axis=1) / dx_mhd +
            np.gradient(vy, axis=0) / dy_mhd)
    ny_mhd, nx_mhd = vx.shape
    reduced_shape = [ny_mhd//nreduce, nx_mhd//nreduce]
    divv = np.asarray(rebin(divv, reduced_shape))
    vx = np.asarray(rebin(vx, reduced_shape))
    vy = np.asarray(rebin(vy, reduced_shape))
    return (vx, vy, divv)


def power_index_dist(plot_config, mhd_config, mhd_run_info, show_plot=True):
    """Plot particle spatial distribution for one energy band

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
        mhd_run_info: information of the MHD run
    """
    tframe = plot_config["tframe"]
    tframe_str = str(tframe).zfill(4)
    run_name = plot_config["run_name"]
    fname = '../data/' + run_name + '/fxy-' + tframe_str + '_sum.dat'
    fdata = np.fromfile(fname)
    nx, = plot_config["nx"]
    ny, = plot_config["ny"]
    print("data size: %d %d" % (nx, ny))
    dists = fdata.reshape([fdata.shape[0] // (nx * ny), -1])
    print("Total number of particles: %d" % np.sum(dists))

    ftot = np.sum(dists, axis=0).reshape(ny, nx)
    sigma = 5
    # ftot = gaussian_filter(ftot, sigma)
    ftot = median_filter(ftot, sigma)
    # ng = 3
    # kernel = np.ones((ng,ng)) / float(ng*ng)
    # ftot = signal.convolve2d(ftot, kernel, mode='same')
    nreduce = plot_config["nreduce"]
    vx, vy, divv = reduced_compression(mhd_run_info, tframe, nreduce)
    dx = mhd_config["lx"] * nreduce / mhd_config["nx"]
    dy = mhd_config["ly"] * nreduce / mhd_config["ny"]
    fadv = (vx * np.gradient(ftot, axis=1) / dx +
            vy * np.gradient(ftot, axis=0) / dy)
    kappa = 0.003
    fdiff = kappa * (np.gradient(np.gradient(ftot, axis=1) / dx, axis=1) / dx +
                     np.gradient(np.gradient(ftot, axis=0) / dy, axis=0) / dy)
    facc = divv * ftot / 3
    alpha = div0(fadv - fdiff, facc)
    print("Mina and max: %f, %f" % (np.min(alpha), np.max(alpha)))
    alpha_selected = alpha[np.abs(alpha) < 30]
    nbins = 100
    abins = np.linspace(-30, 30, nbins + 1)
    sigma = 3
    alpha = gaussian_filter(alpha, sigma)
    # falpha, bins_edges = np.histogram(alpha[128:, 380:388], bins=abins)
    falpha, bins_edges = np.histogram(alpha, bins=abins)
    fdist = np.zeros((2, nbins))
    fdist[0, :] = (abins[1:] + abins[:-1]) * 0.5
    fdist[1, :] = falpha
    # fdir = "../data/pindex_dist/" + mhd_run_info["run_name"] + "/"
    # mkdir_p(fdir)
    # fname = fdir + "pindex_dist_" + str(tframe) + ".dat"
    # fdist.tofile(fname)
    plt.plot(abins[:-1], falpha)
    plt.show()


def power_index_dist2(plot_config, mhd_config, mhd_run_info, show_plot=True):
    """Plot particle spatial distribution for one energy band

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
        mhd_run_info: information of the MHD run
    """
    tframe = plot_config["tframe"]
    tframe_str = str(tframe).zfill(4)
    run_name = plot_config["run_name"]
    fname = '../data/' + run_name + '/fxy-' + tframe_str + '_sum.dat'
    fdata = np.fromfile(fname)
    nx, = plot_config["nx"]
    ny, = plot_config["ny"]
    print("data size: %d %d" % (nx, ny))
    dists = fdata.reshape([fdata.shape[0] // (nx * ny), -1])
    print("Total number of particles: %d" % np.sum(dists))

    ftot = np.sum(dists, axis=0).reshape(ny, nx)
    xs, xe, ys, ye = nx//2, nx, 128, ny
    ftot_sum = np.sum(ftot[ys:ye, xs:xe], axis=0)
    dx, = mhd_config["lx"] / nx
    dy, = mhd_config["ly"] / ny
    xgrid = np.linspace(dx * xs, dx * xe, xe - xs)
    ygrid = np.linspace(dy * ys, dy * ye, ye - ys)
    # p0 is the initial guess for the fitting coefficients
    p0 = [np.max(ftot_sum), (xs + xe) * dx * 0.5, 1.]
    coeff, var_matrix = curve_fit(gauss, xgrid, ftot_sum, p0=p0)
    ftot_fit = gauss(xgrid, *coeff)

    # ftot_sum = median_filter(ftot_sum, 5)
    # ng = 5
    # kernel = np.ones(ng) / float(ng)
    # ftot_sum = signal.convolve(ftot_sum, kernel, 'same')
    # dftot = np.gradient(np.exp(ftot_sum))
    # plt.plot(xgrid, ftot_fit)
    plt.plot(xgrid, ftot_sum)
    plt.show()


def plot_alpha_distribution(mhd_run_info, tstart, tend):
    """Plot the distribution power index alpha

    Args:
        mhd_run_info: information of the MHD run
        log_bin: whether to use logarithm scale
        show_plot: whether to show plot
    """
    fig = plt.figure(figsize=[7, 5])
    rect = [0.13, 0.13, 0.7, 0.8]
    ax = fig.add_axes(rect)
    fdir = "../data/pindex_dist/" + mhd_run_info["run_name"] + "/"
    nframes = tend - tstart + 1
    for tframe in range(tstart, tend + 1):
        fname = fdir + "pindex_dist_" + str(tframe) + ".dat"
        data = np.fromfile(fname)
        dsize, = data.shape
        abins = data[:dsize//2]
        fdist = data[dsize//2:]
        color = plt.cm.jet((tframe - tstart) / float(nframes), 1)
        ax.plot(abins, fdist, color=color, linewidth=2)

    mhd_config = mhd_data.read_mhd_config(mhd_run_info["config_name"],
                                          mhd_run_info["mhd_code"])
    tmin = tstart * mhd_config.dt_out / (mhd_config.xmax - mhd_config.xmin)
    tmax = tend * mhd_config.dt_out / (mhd_config.xmax - mhd_config.xmin)
    rect[0] += rect[2] + 0.01
    rect[2] = 0.04
    cax = fig.add_axes(rect)
    sm = plt.cm.ScalarMappable(cmap=plt.cm.jet,
                               norm=plt.Normalize(vmin=tmin, vmax=tmax))
    # fake up the array of the scalar mappable. Urgh...
    sm._A = []
    cbar = fig.colorbar(sm, cax=cax)
    cbar.set_label(r'$t/\tau_A$', fontsize=20)
    cbar.ax.tick_params(labelsize=16)

    ax.tick_params(labelsize=16)
    ax.set_xlabel(r'$\alpha$', fontsize=20)
    ax.set_ylabel(r'$f(\alpha)$', fontsize=20)
    fpath = '../img/pindex_dist/'
    mkdir_p(fpath)
    plt.savefig(fpath + 'pindex_dist_' + mhd_run_info["run_name"] + '.pdf')
    plt.show()


def get_cmd_args():
    """Get command line arguments """
    default_run_name = 'S1E5_beta01_bg00'
    default_run_dir = ('/net/scratch3/xiaocanli/mhd/guide_field_scaling/' +
                       default_run_name + '/')
    default_sde_run = 'p000_b000_0003_100'
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
    parser.add_argument('--tframe', action="store", default='150', type=int,
                        help= 'Time frame')
    parser.add_argument('--multi_frames', action="store_true", default=False,
                        help='whether to analyze multiple frames')
    parser.add_argument('--tstart', action="store", default='0', type=int,
                        help='starting time frame')
    parser.add_argument('--tend', action="store", default='10', type=int,
                        help='ending time frame')
    parser.add_argument('--calc_alpha', action="store_true", default=False,
                        help='whether to calculate the distribution of alpha')
    parser.add_argument('--plot_alpha', action="store_true", default=False,
                        help='whether to plot the distribution of alpha')
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


def process_input(plot_config, mhd_config, mhd_run_info, args, tframe):
    """process one time frame"""
    plot_config["tframe"] = tframe
    if args.calc_alpha:
        power_index_dist(plot_config, mhd_config, mhd_run_info, show_plot=False)


def analysis_multi_frames(plot_config, mhd_config, mhd_run_info, args):
    """Analysis for multiple time frames
    """
    tframes = range(args.tstart, args.tend + 1)
    ncores = multiprocessing.cpu_count()
    ncores = 8
    Parallel(n_jobs=ncores)(delayed(process_input)(plot_config, mhd_config,
                                                   mhd_run_info, args, tframe)
                            for tframe in tframes)


def main():
    """business logic for when running this module as the primary one!"""
    args = get_cmd_args()
    mhd_run_info = get_mhd_info(args)
    with open('config/spectrum_config_10ta.json', 'r') as file_handler:
        config = json.load(file_handler)
    mhd_config = load_mhd_config(args.mhd_run_dir)
    plot_config = {}
    run_name = args.mhd_run + "/" + args.sde_run
    sde_run_config = config[run_name]
    nreduce = sde_run_config["nreduce"]
    plot_config["nx"] = mhd_config["nx"] // nreduce
    plot_config["ny"] = mhd_config["ny"] // nreduce
    plot_config["nreduce"] = nreduce
    plot_config["tmax"] = sde_run_config["tmax"]
    plot_config["tmin"] = sde_run_config["tmin"]
    plot_config["run_name"] = run_name
    plot_config["tframe"] = args.tframe
    if args.multi_frames:
        analysis_multi_frames(plot_config, mhd_config,
                              mhd_run_info, args)
    else:
        # compression_dist(mhd_run_info, args.tframe)
        if args.calc_alpha:
            # power_index_dist(plot_config, mhd_config,
            #                  mhd_run_info, show_plot=True)
            power_index_dist2(plot_config, mhd_config,
                              mhd_run_info, show_plot=True)
        if args.plot_alpha:
            plot_alpha_distribution(mhd_run_info, args.tstart, args.tend)


if __name__ == "__main__":
    main()

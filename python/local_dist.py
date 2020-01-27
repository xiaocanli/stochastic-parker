#!/usr/bin/env python3
"""
Analysis procedures for local distribution
"""
from __future__ import print_function

import argparse
import itertools
import json
import multiprocessing
import sys

import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from matplotlib.colors import LogNorm

import particle_trajectory as traj
from sde_util import load_mhd_config, mkdir_p

sys.path.insert(0, '/users/xiaocanli/Git/mhd_analysis_sli')
import mhd_data

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

FONT = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 24}


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


def find_nearest(array, value):
    """Find nearest value in an array
    """
    idx = (np.abs(array-value)).argmin()
    return (idx, array[idx])


def reduce_local_dist(plot_config, mhd_config, show_plot=True):
    """Reduce local particle spatial distribution

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    run_name = plot_config["run_name"]
    run_type = plot_config["run_type"]
    tframe = plot_config["tframe"]
    # Normalization parameters depend on runs
    if run_type == "Fan-Early":
        rho0 = 1.2E10  # cm^-3
        L0 = 5.0E9 # in cm
        nrx, nry, nrp = 4, 2, 4
    elif run_type == "Bin-Fan":
        rho0 = 1.0E9  # cm^-3
        L0 = 6.2E9 # in cm
        nrx, nry, nrp = 1, 1, 4
    elif run_type == "Harris_UMN":
        rho0 = 1.0E9  # cm^-3
        L0 = 2.5E9 # in cm
        nrx, nry, nrp = 4, 4, 1

    # Dimensions for local spectra
    nx, = plot_config["nx"]
    ny, = plot_config["ny"]
    nxr = nx // nrx
    nyr = ny // nry
    xstart = nxr//2

    # Read and reduce the local energy spectra
    tframe_str = str(tframe).zfill(4)
    fname = '../data/' + run_name + '/fp_local_' + tframe_str + '_sum.dat'
    fdata = np.fromfile(fname)
    npp = fdata.shape[0] // (nx * ny)
    fdata = fdata.reshape([ny, nx, npp])
    print("Total number of particles: %d" % np.sum(fdata))
    print("data size: %d %d %d (C-order)" % fdata.shape)
    dists = fdata.reshape(ny//nry, nry, nx//nrx, nrx, npp//nrp, nrp)
    dists_r = np.sum(np.sum(np.sum(dists, axis=5), axis=3), axis=1)
    print("reduced data size: %d %d %d (C-order)" % dists_r.shape)

    # Momentum bins (simulation unit)
    pmin = 1E-2
    pmax = 1E1
    p0 = 1E-1
    pmom = np.logspace(math.log10(pmin), math.log10(pmax), npp + 1)
    pmom_mid = (pmom[:-1] + pmom[1:]) * 0.5
    dpmom = np.diff(pmom)
    pmom_r = pmom[::nrp]
    pmom_mid_r = (pmom_r[:-1] + pmom_r[1:]) * 0.5
    dpmom_r = np.diff(pmom_r)
    fmom_r = dists_r * pmom_mid_r / dpmom_r  # f(p)p^3
    fene_r = fmom_r / pmom_mid_r*2  # f(p)p or f(e)

    # Normalized energy bins
    ene0 = 1.0  # 1keV
    ene_shift = 1.0  # parameter to adjust the energy of the injected particles
    elog_r = ene0 * pmom_r**2 / p0**2  # in keV
    elog_r *= ene_shift  # shift the injection energy
    elog_mid_r = (elog_r[:-1] + elog_r[1:]) * 0.5
    delog_r = np.diff(elog_r)

    # Normalize the spectrum
    # We assume the initially inject particles are about 5% of all particles
    nptl_tot_real = L0**2 * rho0
    ene_cut = 10 * ene0 # about 1% of the simulated particles
    cutoff, _ = find_nearest(pmom, math.sqrt(ene_cut) * p0)
    cutoff_r = cutoff // nrp
    nptl_tot = np.sum(fdata)
    nptl_above_cutoff = np.sum(fdata[:, :, cutoff:])
    fnorm = nptl_tot_real * 0.05 / nptl_tot
    fnorm *= nxr * nyr / L0**2
    print("Assuming nonthermal start from %f keV" %
          ((pmom[cutoff]**2 * ene_shift) / p0**2))
    print("Total Number of particles: %d" % nptl_tot)
    print("Number of particles above break: %d" % nptl_above_cutoff)
    print("Non-thermal fraction: %f" % (nptl_above_cutoff/nptl_tot))

    if plot_config["check_dist"]:
        # Plot the spatial distribution of the high-energy electrons
        L0_Mm = L0 / 1E8  # to Mm
        dists_r *= fnorm
        dist_2d = np.sum(dists_r[:, :, cutoff_r+1:], axis=2)
        if run_type == "Fan-Early":
            rect = [0.12, 0.10, 0.7, 0.85]
            fig = plt.figure(figsize=[7, 10])
            ax1 = fig.add_axes(rect)
            extent_box = [L0_Mm*0.5, L0_Mm, 0, L0_Mm]
            vmin, vmax = 0, 1E9
        elif run_type == "Bin-Fan":
            rect = [0.12, 0.12, 0.7, 0.85]
            fig = plt.figure(figsize=[8, 7])
            ax1 = fig.add_axes(rect)
            extent_box = [0, L0_Mm, 0, L0_Mm]
            vmin, vmax = 0, 1E9
        if run_type == "Harris_UMN":
            rect = [0.12, 0.10, 0.68, 0.85]
            fig = plt.figure(figsize=[7, 10])
            ax1 = fig.add_axes(rect)
            extent_box = [L0_Mm*0.5, L0_Mm, 0, L0_Mm]
            vmin, vmax = 0, 2E7
        img = ax1.imshow(dist_2d[:, xstart:], cmap=plt.cm.inferno,
                         vmin=vmin, vmax=vmax,
                         extent=extent_box,
                         aspect='auto', origin='lower',
                         interpolation='bicubic')
        ecut = elog_r[cutoff_r+1]
        ecut_s = "{%0.1f}" % ecut
        label1 = r'$\varepsilon > ' + ecut_s + '$ keV'
        ax1.text(0.05, 0.95, label1, color='white', fontsize=24,
                 horizontalalignment='left', verticalalignment='center',
                 transform=ax1.transAxes)
        ax1.set_xlabel(r'$x$/Mm', fontdict=FONT, fontsize=24)
        ax1.set_ylabel(r'$y$/Mm', fontdict=FONT, fontsize=24)
        ax1.tick_params(labelsize=20)
        rect[0] += rect[2] + 0.02
        rect[2] = 0.04
        cbar_ax = fig.add_axes(rect)
        cbar = fig.colorbar(img, cax=cbar_ax)
        cbar.ax.tick_params(labelsize=16)
        cbar.ax.set_ylabel(r'$\rho$ (cm$^{-3}$)', fontsize=24)
        fdir = '../img/fp_local/' + run_name + '/'
        mkdir_p(fdir)
        fname = fdir + "fp_local_reduced_" + str(tframe) + ".jpg"
        fig.savefig(fname, dpi=200)

    # Save the reduced spectrum
    fdir = '../data/' + run_name + '/reduced/'
    mkdir_p(fdir)
    fname = fdir + "fe_local_reduced_" + str(tframe) + ".dat"
    fene_r *= fnorm  # normalize to real number density
    fene_r *= 0.5 * p0**2 / (ene0 * ene_shift)  # normalize to keV^-1
    fene_r[:, xstart:, :].tofile(fname)

    if plot_config["check_dist"]:
        # Check the saved local spectrum
        fdata = np.fromfile(fname)
        fdata = fdata.reshape([nyr, nxr-xstart, -1])
        fene_r = fdata * delog_r[np.newaxis, np.newaxis, :]
        dist_2d = np.sum(fene_r[:, :, cutoff_r+1:], axis=2)
        rect = [0.12, 0.10, 0.68, 0.85]
        fig = plt.figure(figsize=[7, 10])
        ax1 = fig.add_axes(rect)
        img = ax1.imshow(dist_2d, cmap=plt.cm.inferno,
                         vmin=vmin, vmax=vmax,
                         extent=extent_box,
                         aspect='auto', origin='lower',
                         interpolation='bicubic')

        if show_plot:
            plt.show()
        else:
            plt.close('all')


def reduce_mhd_data(plot_config, mhd_config, mhd_run_info):
    """Reduce MHD data size

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
        mhd_run_info: information of the MHD run
    """
    run_type = plot_config["run_type"]
    if run_type == "Fan-Early":
        rho0 = 1.2E10  # cm^-3
        b0 = 50 # Gauss
        T0 = 6.0E6  # K
        beta0 = 0.1
        nrx, nry = 64, 32
    if run_type == "Harris_UMN":
        rho0 = 1.0E9  # cm^-3
        b0 = 50 # Gauss
        T0 = 1.0E6  # K
        beta0 = 0.1
        nrx, nry = 64, 64
    tframe = plot_config["tframe"]
    xmesh, ymesh, data = mhd_data.read_fields_data(mhd_run_info, tframe)
    ny, nx = xmesh.shape
    mhd_box = [xmesh[0, 0], xmesh[0, -1], ymesh[0, 0], ymesh[-1, 0], nx, ny]
    rho = data[:, :, 0].T
    pre = data[:, :, 1].T
    bx = data[:, :, 5].T
    by = data[:, :, 6].T
    bz = data[:, :, 7].T
    rho = rho.reshape(ny//nry, nry, nx//nrx, nrx)
    pre = pre.reshape(ny//nry, nry, nx//nrx, nrx)
    bx = bx.reshape(ny//nry, nry, nx//nrx, nrx)
    by = by.reshape(ny//nry, nry, nx//nrx, nrx)
    bz = bz.reshape(ny//nry, nry, nx//nrx, nrx)
    rho_r = np.mean(np.mean(rho, axis=3), axis=1)
    pre_r = np.mean(np.mean(pre, axis=3), axis=1)
    bx_r = np.mean(np.mean(bx, axis=3), axis=1)
    by_r = np.mean(np.mean(by, axis=3), axis=1)
    bz_r = np.mean(np.mean(bz, axis=3), axis=1)
    absB_r = np.sqrt(bx_r**2 + by_r**2 + bz_r**2)
    T_r = pre_r / rho_r / (beta0 * 0.5)
    rho_r *= rho0
    bx *= b0
    by *= b0
    bz *= b0
    absB_r *= b0
    T_r *= T0

    nxr = nx // nrx
    nyr = ny // nry

    fdir = plot_config["mhd_run_dir"] + "data_reduced/"
    mkdir_p(fdir)
    fname = fdir + 'rho_' + str(tframe) + '.dat'
    rho_r[:, nxr//2:].tofile(fname)
    fname = fdir + 'T_' + str(tframe) + '.dat'
    T_r[:, nxr//2:].tofile(fname)
    fname = fdir + 'bx_' + str(tframe) + '.dat'
    bx_r[:, nxr//2:].tofile(fname)
    fname = fdir + 'by_' + str(tframe) + '.dat'
    by_r[:, nxr//2:].tofile(fname)
    fname = fdir + 'bz_' + str(tframe) + '.dat'
    bz_r[:, nxr//2:].tofile(fname)
    fname = fdir + 'absB_' + str(tframe) + '.dat'
    absB_r[:, nxr//2:].tofile(fname)


def calc_va(b0, n0, verbose=False):
    """Calculate the Alfven speed in m/s

    Args:
        b0: magnetic field strength in Gauss
        n0: particle number density in cm^-3
    """
    pmass = 1.6726219E-27  # in kilogram
    mu0 = 4 * math.pi * 1E-7
    va = b0 * 1E-4 / math.sqrt(mu0 * n0 * 1E6 * pmass)
    print("The Alfven speed is %f km/s" % (va/1E3))
    return va


def plot_reduced_mhd(plot_config, mhd_config, mhd_run_info, show_plot=True):
    """Plot reduced MHD data
    """
    run_type = plot_config["run_type"]
    if run_type == "Fan-Early":
        rho0 = 1.2E10  # cm^-3
        b0 = 50 # Gauss
        T0 = 6.0E6  # K
        L0 = 5.0E9 # in cm
        nrx, nry = 64, 32
    elif run_type == "Harris_UMN":
        rho0 = 1.0E9  # cm^-3
        b0 = 50 # Gauss
        T0 = 1.0E6  # K
        L0 = 2.5E9 # in cm
        nrx, nry = 64, 64
    va = calc_va(b0, rho0) # in m/s
    time_norm = L0 / (va * 1E2)
    rho_norm = 1.0E9  # cm^-3
    bnorm = 50 # Gauss
    Tnorm = 1.0E6  # K
    nx, = mhd_config["nx"]
    ny, = mhd_config["ny"]
    nxr = nx // nrx
    nyr = ny // nry
    fdir = plot_config["mhd_run_dir"] + "data_reduced/"
    tframe = plot_config["tframe"]
    L0_Mm = L0 / 1E8
    extent_box = [0, L0_Mm*0.5, 0, L0_Mm]

    fig = plt.figure(figsize=[7, 5])
    rect = [0.09, 0.13, 0.28, 0.8]
    hgap, vgap = 0.02, 0.02

    fname = fdir + 'rho_' + str(tframe) + '.dat'
    fdata = np.fromfile(fname, dtype=np.float32)
    fdata = fdata.reshape((nyr, nxr//2))
    ax = fig.add_axes(rect)
    img = ax.imshow(fdata/rho_norm, extent=extent_box, cmap=plt.cm.plasma,
                    aspect='auto', origin='lower',
                    interpolation='bicubic',
                    vmin=1, vmax=10)
    ax.tick_params(labelsize=12)
    ax.set_xlabel(r'$x$/Mm', fontsize=16)
    ax.set_ylabel(r'$y$/Mm', fontsize=16)
    ax.tick_params(bottom=True, top=False, left=True, right=True)
    ax.tick_params(axis='x', which='minor', direction='in', top=True)
    ax.tick_params(axis='x', which='major', direction='in', top=True)
    ax.tick_params(axis='y', which='minor', direction='in')
    ax.tick_params(axis='y', which='major', direction='in')
    rect_cbar = np.copy(rect)
    rect_cbar[0] += 0.02
    rect_cbar[2] = 0.015
    rect_cbar[1] = 0.4
    rect_cbar[3] = 0.3
    cbar_ax = fig.add_axes(rect_cbar)
    cbar = fig.colorbar(img, cax=cbar_ax)
    cbar.set_label(r'$\rho (10^9\text{cm}^{-3})$', color='w', fontsize=12)
    cbar.ax.yaxis.set_tick_params(color='w')
    cbar.outline.set_edgecolor('w')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')
    cbar.ax.tick_params(labelsize=12, color='w')

    rect[0] += rect[2] + hgap
    fname = fdir + 'T_' + str(tframe) + '.dat'
    fdata = np.fromfile(fname, dtype=np.float32)
    fdata = fdata.reshape((nyr, nxr//2))
    ax = fig.add_axes(rect)
    img = ax.imshow(fdata/Tnorm, extent=extent_box, cmap=plt.cm.viridis,
                    aspect='auto', origin='lower',
                    interpolation='bicubic',
                    vmin=0.1, vmax=10)
    ax.tick_params(labelsize=12)
    ax.tick_params(axis='y', labelleft=False)
    ax.set_xlabel(r'$x$/Mm', fontsize=16)
    time = mhd_config["dt_out"][0] * time_norm * tframe
    tname = "Time: %0.1f s" % time
    plt.title(tname, fontsize=16)
    ax.tick_params(bottom=True, top=False, left=True, right=True)
    ax.tick_params(axis='x', which='minor', direction='in', top=True)
    ax.tick_params(axis='x', which='major', direction='in', top=True)
    ax.tick_params(axis='y', which='minor', direction='in')
    ax.tick_params(axis='y', which='major', direction='in')
    rect_cbar = np.copy(rect)
    rect_cbar[0] += 0.02
    rect_cbar[2] = 0.015
    rect_cbar[1] = 0.4
    rect_cbar[3] = 0.3
    cbar_ax = fig.add_axes(rect_cbar)
    cbar = fig.colorbar(img, cax=cbar_ax)
    cbar.set_label(r'$T (10^6\text{K})$', color='w', fontsize=12)
    cbar.ax.yaxis.set_tick_params(color='w')
    cbar.outline.set_edgecolor('w')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')
    cbar.ax.tick_params(labelsize=12, color='w')

    rect[0] += rect[2] + hgap
    fname = fdir + 'absB_' + str(tframe) + '.dat'
    fdata = np.fromfile(fname, dtype=np.float32)
    fdata = fdata.reshape((nyr, nxr//2))
    ax = fig.add_axes(rect)
    img = ax.imshow(fdata/bnorm, extent=extent_box, cmap=plt.cm.plasma,
                    aspect='auto', origin='lower',
                    interpolation='bicubic',
                    vmin=0.0, vmax=2.0)
    ax.tick_params(labelsize=12)
    ax.tick_params(axis='y', labelleft=False)
    ax.set_xlabel(r'$x$/Mm', fontsize=16)
    ax.tick_params(bottom=True, top=False, left=True, right=True)
    ax.tick_params(axis='x', which='minor', direction='in', top=True)
    ax.tick_params(axis='x', which='major', direction='in', top=True)
    ax.tick_params(axis='y', which='minor', direction='in')
    ax.tick_params(axis='y', which='major', direction='in')
    rect_cbar = np.copy(rect)
    rect_cbar[0] += 0.02
    rect_cbar[2] = 0.015
    rect_cbar[1] = 0.4
    rect_cbar[3] = 0.3
    cbar_ax = fig.add_axes(rect_cbar)
    cbar = fig.colorbar(img, cax=cbar_ax)
    cbar.set_label(r'$B (\text{50 Gauss})$', color='w', fontsize=12)
    cbar.ax.yaxis.set_tick_params(color='w')
    cbar.outline.set_edgecolor('w')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')
    cbar.ax.tick_params(labelsize=12, color='w')

    img_dir = '../img/' + mhd_run_info["run_name"] + '/mhd_fields_reduced/'
    mkdir_p(img_dir)
    fname = img_dir + 'mhd_fields_' + str(tframe) + '.jpg'
    fig.savefig(fname, dpi=400)
    if not show_plot:
        plt.close()
    if show_plot:
        plt.show()


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


def get_cmd_args():
    """Get command line arguments
    """
    default_mhd_run = 'S1E5_beta01_bg00'
    default_mhd_run_dir = ('/net/scratch3/xiaocanli/mhd/guide_field_scaling/' +
                           default_mhd_run + '/')
    default_sde_run = 'p000_b000_00047_100_l'
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
    parser.add_argument('--tframe', action="store", default='200', type=int,
                        help='Time frame')
    parser.add_argument('--multi_frames', action="store_true", default=False,
                        help='whether to analyze multiple frames')
    parser.add_argument('--time_loop', action="store_true", default=False,
                        help='whether to loop over time instead of using joblib')
    parser.add_argument('--tstart', action="store", default='0', type=int,
                        help='starting time frame')
    parser.add_argument('--tend', action="store", default='10', type=int,
                        help='ending time frame')
    parser.add_argument('--local_dist', action="store_true", default=False,
                        help='whether to plot spatial distribution')
    parser.add_argument('--check_dist', action="store_true", default=False,
                        help='whether to check local spatial distribution')
    parser.add_argument('--reduce_mhd', action="store_true", default=False,
                        help='whether to reduce MHD data size')
    parser.add_argument('--run_type', action="store", default="Fan-Early",
                        help='What kind of run')
    parser.add_argument('--plot_rmhd', action="store_true", default=False,
                        help='whether to plot reduced MHD data')
    return parser.parse_args()


def analysis_single_frame(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    mhd_run_info = get_mhd_info(args)
    if args.local_dist:
        reduce_local_dist(plot_config, mhd_config, show_plot=True)
    elif args.reduce_mhd:
        reduce_mhd_data(plot_config, mhd_config, mhd_run_info)
    elif args.plot_rmhd:
        plot_reduced_mhd(plot_config, mhd_config, mhd_run_info)


def process_input(plot_config, mhd_config, args, tframe):
    """process one time frame"""
    plot_config["tframe"] = tframe
    mhd_run_info = get_mhd_info(args)
    if args.local_dist:
        reduce_local_dist(plot_config, mhd_config, show_plot=False)
    elif args.reduce_mhd:
        reduce_mhd_data(plot_config, mhd_config, mhd_run_info)
    elif args.plot_rmhd:
        plot_reduced_mhd(plot_config, mhd_config, mhd_run_info, show_plot=False)


def analysis_multi_frames(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    tframes = range(plot_config["tmin"], plot_config["tmax"] + 1)
    mhd_run_info = get_mhd_info(args)
    if args.time_loop:
        if args.local_dist:
            for tframe in tframes:
                print("Time frame: %d" % tframe)
                plot_config["tframe"] = tframe
                reduce_local_dist(plot_config, mhd_config, show_plot=False)
        elif args.plot_rmhd:
            for tframe in tframes:
                print("Time frame: %d" % tframe)
                plot_config["tframe"] = tframe
                plot_reduced_mhd(plot_config, mhd_config,
                                 mhd_run_info, show_plot=False)
    else:
        ncores = multiprocessing.cpu_count()
        ncores = 8
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
    plot_config["run_name"] = run_name
    plot_config["tframe"] = args.tframe
    plot_config["mhd_run_dir"] = args.mhd_run_dir
    plot_config["run_type"] = args.run_type
    plot_config["check_dist"] = args.check_dist
    tframes = [150, 175, 200]
    plot_config["tframes"] = tframes
    if args.multi_frames:
        analysis_multi_frames(plot_config, mhd_config, args)
    else:
        analysis_single_frame(plot_config, mhd_config, args)


if __name__ == "__main__":
    main()

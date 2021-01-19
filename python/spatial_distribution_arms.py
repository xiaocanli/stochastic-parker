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

from evtk.hl import gridToVTK
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from matplotlib.colors import LogNorm, SymLogNorm
from scipy.interpolate import interp2d

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

    # for iband in range(nbands):
    #     print("Band %d" % iband)
    #     fname = '../data/' + run_name + '/dist_' + tframe_str + '_' + str(iband)
    #     gridToVTK(fname, p3d, t3d, r3d, pointData = {"pressure" : fdata[iband]})

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
    nmin, nmax = 10, 2E4

    rr, tt = np.meshgrid(rpos_r, tpos_r)
    xx = rr * np.cos(tt)
    yy = rr * np.sin(tt)
    rpos_new = np.linspace(rpos_r.min(), rpos_r.max(), npr*nreduce)
    tpos_new = np.linspace(tpos_r.min(), tpos_r.max(), ntr*nreduce)
    rr_new, tt_new = np.meshgrid(rpos_new, tpos_new)
    xx_new = rr_new * np.cos(tt_new)
    yy_new = rr_new * np.sin(tt_new)

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
        f = interp2d(rpos_r, tpos_r, fdata_sum, kind='cubic')
        fdata_new = f(rpos_new, tpos_new)
        fdata_new[fdata_new<0] = 1E-5
        im = ax.pcolormesh(xx_new, yy_new, fdata_new,
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
    rpos_new = np.linspace(rpos_r.min(), rpos_r.max(), npr*nreduce)
    ppos_new = np.linspace(ppos_r.min(), ppos_r.max(), ntr*nreduce)
    rr_new, pp_new = np.meshgrid(rpos_new, ppos_new)
    xx_new = rr_new * np.sin(pp_new)
    yy_new = rr_new * np.cos(pp_new)

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
        f = interp2d(rpos_r, ppos_r, fdata_sum, kind='cubic')
        fdata_new = f(rpos_new, ppos_new)
        fdata_new[fdata_new<0] = 1E-5
        im = ax.pcolormesh(xx_new, yy_new, fdata_new,
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


def plot_dist_2d_cuts(plot_config, mhd_config, show_plot=True):
    """Plot 2D cuts of spatial distribution

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
    ng = 2  # number of ghost cells

    mhd_run_dir = plot_config["mhd_run_dir"]
    fname = mhd_run_dir + "bin_data/xpos.dat"
    rpos = np.fromfile(fname)
    fname = mhd_run_dir + "bin_data/ypos.dat"
    tpos = np.fromfile(fname) + math.pi * 0.5
    fname = mhd_run_dir + "bin_data/zpos.dat"
    ppos = np.fromfile(fname)

    nx = mhd_config["nx"][0] + ng*2
    ny = mhd_config["ny"][0] + ng*2
    nz = mhd_config["nz"][0] + ng*2
    lx = rpos.max() - rpos.min()
    ly = tpos.max() - tpos.min()
    lz = ppos.max() - ppos.min()
    dr = lx / nx
    dt = ly / ny
    dp = lz / nz
    fpath = plot_config["mhd_run_dir"] + 'bin_data/'
    fname = fpath + 'mhd_data_' + str(tframe).zfill(4)
    mhd_fields = np.fromfile(fname, dtype=np.float32)
    mhd_fields = mhd_fields.reshape((nz, ny, nx, 8))

    p3d, t3d, r3d = np.meshgrid(ppos, tpos, rpos, indexing='ij')
    divv = (np.gradient(mhd_fields[:, :, :, 0], dr, axis=2) +
            (2.0 * mhd_fields[:, :, :, 0] +
             np.gradient(mhd_fields[:, :, :, 1], dt, axis=1)) / r3d +
            (np.cos(t3d) * mhd_fields[:, :, :, 1] +
             np.gradient(mhd_fields[:, :, :, 2], dp, axis=0)) / (r3d * np.sin(t3d)))

    rpos_r = rpos[ng:-ng:nreduce]
    tpos_r = tpos[ng:-ng:nreduce]
    ppos_r = ppos[ng:-ng:nreduce]

    run_name = plot_config["run_name"]
    tframe_str = str(tframe).zfill(4)
    fname = '../data/' + run_name + '/fxy-' + tframe_str + '_sum.dat'
    fdata = np.fromfile(fname).reshape([-1, npr, ntr, nrr])
    p3d, t3d, r3d = np.meshgrid(ppos_r, tpos_r, rpos_r, indexing='ij')
    fdata /= r3d**2 * np.sin(t3d)
    nbands, _, _, _ = fdata.shape

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
    rpos_new = np.linspace(rpos_r.min(), rpos_r.max(), npr*nreduce)
    tpos_new = np.linspace(tpos_r.min(), tpos_r.max(), ntr*nreduce)
    rr_new, tt_new = np.meshgrid(rpos_new, tpos_new)
    xx_new = rr_new * np.cos(tt_new)
    yy_new = rr_new * np.sin(tt_new)
    rr_o, tt_o = np.meshgrid(rpos, tpos)
    xx_o = rr_o * np.cos(tt_o)
    yy_o = rr_o * np.sin(tt_o)

    fnorms = [1, 4, 16, 64]

    for iphi in range(npr):
        print("phi %d of %d" % (iphi, npr))
        fig = plt.figure(figsize=[20, 10])
        rect0 = [0.05, 0.55, 0.27, 0.40]
        hgap, vgap = 0.03, 0.03
        bands, bande = plot_config["high_bands"]
        for iband, eband in enumerate(range(bands, bande+1)):
            row = iband // 2
            col = iband % 2
            rect = np.copy(rect0)
            rect[0] = rect0[0] + col * (rect0[2] + hgap)
            rect[1] = rect0[1] - row * (rect0[3] + vgap)
            ax = fig.add_axes(rect)
            fdata_cut = fdata[eband, iphi, :, :] * fnorms[iband]
            f = interp2d(rpos_r, tpos_r, fdata_cut, kind='cubic')
            fdata_new = f(rpos_new, tpos_new)
            fdata_new[fdata_new<=0] = 1E-5
            im = ax.pcolormesh(xx_new, yy_new, fdata_new,
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

        rect = np.copy(rect0)
        rect[0] += rect0[2] * 2 + hgap + 0.08
        rect[1] -= (rect[3] - vgap) * 0.5
        ax = fig.add_axes(rect)
        im = ax.pcolormesh(xx_o, yy_o, divv[iphi*nreduce, :, :],
                           cmap=plt.cm.seismic,
                           vmin=-1E5, vmax=1E5)
        ax.set_xlabel(r'$x/R_e$', fontsize=16)
        ax.set_ylabel(r'$y/R_e$', fontsize=16)
        ax.tick_params(labelsize=12)
        fig_text = r'$\nabla\cdot\boldsymbol{V}$'
        ax.text(0.05, 0.85, fig_text, color='k', fontsize=16,
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)
        rect_cbar = np.copy(rect)
        rect_cbar[1] -= 0.1
        rect_cbar[3] = 0.02
        cbar_ax = fig.add_axes(rect_cbar)
        cbar1 = fig.colorbar(im, cax=cbar_ax, extend='both',
                             orientation="horizontal")
        cbar1.ax.tick_params(labelsize=16)

        L0 = 6.9634E8  # m
        v0 = 1E3       # m/s
        t0 = L0 / v0
        tva = mhd_config["dt_out"][0] * t0
        title = 'Frame: %d' % plot_config["tframe"]
        fig.suptitle(title, fontsize=24)

        fdir = '../img/nrho/' + run_name + '/'
        mkdir_p(fdir)
        fname = fdir + 'nrho_bands_phi_' + tframe_str + '_' + str(iphi) + '.jpg'
        fig.savefig(fname)
        # plt.show()
        plt.close()

    rr, pp = np.meshgrid(rpos_r, ppos_r)
    xx = rr * np.sin(pp)
    yy = rr * np.cos(pp)
    rpos_new = np.linspace(rpos_r.min(), rpos_r.max(), npr*nreduce)
    ppos_new = np.linspace(ppos_r.min(), ppos_r.max(), ntr*nreduce)
    rr_new, pp_new = np.meshgrid(rpos_new, ppos_new)
    xx_new = rr_new * np.sin(pp_new)
    yy_new = rr_new * np.cos(pp_new)
    rr_o, pp_o = np.meshgrid(rpos, ppos)
    xx_o = rr_o * np.sin(pp_o)
    yy_o = rr_o * np.cos(pp_o)

    fnorms = [1, 4, 16, 64]

    for itheta in range(ntr):
        print("theta %d of %d" % (itheta, ntr))
        fig = plt.figure(figsize=[20, 10])
        rect0 = [0.05, 0.55, 0.27, 0.40]
        hgap, vgap = 0.03, 0.03
        bands, bande = plot_config["high_bands"]
        for iband, eband in enumerate(range(bands, bande+1)):
            row = iband // 2
            col = iband % 2
            rect = np.copy(rect0)
            rect[0] = rect0[0] + col * (rect0[2] + hgap)
            rect[1] = rect0[1] - row * (rect0[3] + vgap)
            ax = fig.add_axes(rect)
            fdata_cut = fdata[eband, :, itheta, :] * fnorms[iband]
            f = interp2d(rpos_r, ppos_r, fdata_cut, kind='cubic')
            fdata_new = f(rpos_new, ppos_new)
            fdata_new[fdata_new<0] = 1e-5
            im = ax.pcolormesh(xx_new, yy_new, fdata_new,
                               norm=LogNorm(vmin=nmin, vmax=nmax))
            ax.tick_params(labelsize=12)
            if row == 1:
                ax.set_xlabel(r'$x/r_e$', fontsize=16)
            else:
                ax.tick_params(axis='x', labelleft=False)
            if col == 0:
                ax.set_ylabel(r'$z/r_e$', fontsize=16)
            else:
                ax.tick_params(axis='y', labelleft=False)
            norm = fnorms[iband]
            e0 = plot_config["e0"]
            enel = "{%0.0f}" % (ebins_edges[eband] * e0)
            eneh = "{%0.0f}" % (ebins_edges[eband+1] * e0)
            ftext = (r'$' + enel + r'$kev' + r'$' +
                     r'< \varepsilon <' +
                     eneh + r'$kev')
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

        rect = np.copy(rect0)
        rect[0] += rect0[2] * 2 + hgap + 0.08
        rect[1] -= (rect[3] - vgap) * 0.5
        ax = fig.add_axes(rect)
        im = ax.pcolormesh(xx_o, yy_o, divv[:, itheta*nreduce, :],
                           cmap=plt.cm.seismic,
                           vmin=-1E5, vmax=1E5)
        ax.set_xlabel(r'$x/R_e$', fontsize=16)
        ax.set_ylabel(r'$y/R_e$', fontsize=16)
        ax.tick_params(labelsize=12)
        fig_text = r'$\nabla\cdot\boldsymbol{V}$'
        ax.text(0.05, 0.85, fig_text, color='k', fontsize=16,
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)

        l0 = 6.9634e8  # m
        v0 = 1e3       # m/s
        t0 = l0 / v0
        tva = mhd_config["dt_out"][0] * t0
        title = 'frame: %d' % plot_config["tframe"]
        fig.suptitle(title, fontsize=24)

        fdir = '../img/nrho/' + run_name + '/'
        mkdir_p(fdir)
        fname = fdir + 'nrho_bands_theta_' + tframe_str + '_' + str(itheta) + '.jpg'
        fig.savefig(fname)
        plt.close()


def get_cmd_args():
    """get command line arguments
    """
    default_mhd_run = 'arms_test'
    default_mhd_run_dir = ('/net/scratch3/xiaocan/mhd/' + default_mhd_run + '/')
    default_sde_run = 'test_run'
    parser = argparse.ArgumentParser(description='spatial distribution')
    parser.add_argument('--mhd_run', action="store",
                        default=default_mhd_run, help='mhd run name')
    parser.add_argument('--mhd_run_dir', action="store",
                        default=default_mhd_run_dir, help='mhd run directory')
    parser.add_argument('--sde_run', action="store",
                        default=default_sde_run, help='sde run name')
    parser.add_argument('--tframe', action="store", default='100', type=int,
                        help='time frame')
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
                        help='whether to plot 2d spatial distributions')
    parser.add_argument('--dist_2d_cuts', action="store_true", default=False,
                        help='whether to plot 2d cuts of the spatial distributions')
    return parser.parse_args()


def analysis_single_frames(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    if args.dist_2d:
        plot_dist_2d(plot_config, mhd_config, show_plot=True)
    elif args.dist_2d_cuts:
        plot_dist_2d_cuts(plot_config, mhd_config, show_plot=True)


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

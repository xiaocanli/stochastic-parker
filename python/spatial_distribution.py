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

from sde_util import load_mhd_config, mkdir_p

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]


def figure_text(eband, e0=1.0):
    """get figure text for different energy band

    Args:
        eband: energy band
        e0: initial energy
    """
    if eband == 0:
        eneh = "{%0.0f}" % ((2**eband * 0.75)**2 * e0)
        fig_text = r'$\varepsilon <' + eneh + r'$keV'
    elif eband == 5:
        enel = "{%0.0f}" % ((2**(eband - 2) * 1.5)**2 * e0)
        fig_text = r'$\varepsilon >' + enel + r'$keV'
    else:
        enel = "{%0.0f}" % ((2**(eband - 1) * 0.75)**2 * e0)
        eneh = "{%0.0f}" % ((2**(eband - 1) * 1.5)**2 * e0)
        fig_text = (r'$' + enel + r'$keV' + r'$'
                    r'< \varepsilon <' +
                    eneh + r'$keV')

    return fig_text


def spatial_distribution_multi(plot_config, mhd_config, show_plot=True):
    """Plot spatial distributions for multiple energy band

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    mpl.rc('text', usetex=True)
    fig = plt.figure(figsize=[10, 10])
    rect0 = [0.08, 0.52, 0.40, 0.40]
    rect = np.copy(rect0)
    hgap, vgap = 0.04, 0.02
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
    for i, j in itertools.product(range(2), range(2)):
        rect[0] = rect0[0] + (hgap + rect0[2]) * i
        rect[1] = rect0[1] - (vgap + rect0[3]) * j
        eband = j * 2 + i + 1
        fband = dists[eband, :]
        fnorm = 4**(eband - 1)
        fband = np.reshape(fband, (ny, nx)) * fnorm
        fband += vmin * 0.01
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
        fig_text = str(fnorm) + r'$n($' + figure_text(eband) + r'$)$'
        color = 'k' if eband == 0 else 'w'
        axs[eband-1].text(0.02, 0.9, fig_text, color=color, fontsize=20,
                          bbox=dict(facecolor='none', alpha=1.0,
                                    edgecolor='none', pad=10.0),
                          horizontalalignment='left', verticalalignment='bottom',
                          transform=axs[eband-1].transAxes)

    axs[2].set_xlabel(r'$x$', fontsize=20)
    axs[3].set_xlabel(r'$x$', fontsize=20)
    axs[0].set_ylabel(r'$y$', fontsize=20)
    axs[2].set_ylabel(r'$y$', fontsize=20)

    axs[0].tick_params(axis='x', labelbottom=False)
    axs[1].tick_params(axis='x', labelbottom=False)
    axs[1].tick_params(axis='y', labelleft=False)
    axs[3].tick_params(axis='y', labelleft=False)

    rect[0] = rect0[0] + rect0[2] * 2 + hgap + 0.01
    rect[2] = 0.015
    rect[3] = rect0[3] * 2 + vgap
    rect[1] += rect[3] * 0.25
    rect[3] *= 0.5
    cbar_ax = fig.add_axes(rect)
    cbar1 = fig.colorbar(ims[0], cax=cbar_ax, extend='both')
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
    print("particle number in band %d: %d" % (eband, np.sum(fband)))
    fband = np.reshape(fband, (ny, nx))
    print("min, max, mean, and std: %f %f %f %f" %
          (np.min(fband), np.max(fband), np.mean(fband), np.std(fband)))

    fig = plt.figure(figsize=[7, 6.5])
    rect = [0.1, 0.12, 0.78, 0.84]
    ax = fig.add_axes(rect)
    vmin = plot_config['nmin']
    vmax = plot_config['nmax']
    sizes = [mhd_config["xmin"][0], mhd_config["xmax"][0],
             mhd_config["ymin"][0], mhd_config["ymax"][0]]
    img = ax.imshow(fband + 0.1*vmin, cmap=plt.cm.hot, aspect='auto',
                    origin='lower', extent=sizes,
                    # vmin=vmin, vmax=vmax,
                    norm=LogNorm(vmin=vmin, vmax=vmax),
                    interpolation='bicubic')
    ax.tick_params(labelsize=12)
    fig_text = figure_text(eband)
    ax.text(0.02, 0.9, fig_text, color='k', fontsize=20,
            bbox=dict(facecolor='none', alpha=1.0, edgecolor='none', pad=10.0),
            horizontalalignment='left', verticalalignment='bottom',
            transform=ax.transAxes)

    ax.set_xlabel(r'$x$', fontsize=16)
    ax.set_ylabel(r'$y$', fontsize=16)

    rect[0] += rect[2] + 0.01
    rect[2] = 0.02
    rect[1] += rect[3] * 0.25
    rect[3] *= 0.5
    cbar_ax = fig.add_axes(rect)
    cbar1 = fig.colorbar(img, cax=cbar_ax, extend='both')
    cbar1.ax.tick_params(labelsize=16)

    fdir = '../img/nrho/' + run_name + '/'
    mkdir_p(fdir)
    fname = fdir + 'nrho_' + tframe_str + '_' + str(eband) + '.jpg'
    fig.savefig(fname, dpi=200)

    if show_plot:
        plt.show()
    else:
        plt.close()


def spatial_distribution_tri_frames(plot_config, mhd_config, show_plot=True):
    """Plot particle spatial distribution for one energy band for three frames

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    # tframes = [100, 125, 160]
    tframes = [147, 190, 201]
    fig = plt.figure(figsize=[6, 4])
    rect0 = [0.1, 0.12, 0.24, 0.8]
    rect = np.copy(rect0)
    hgap = 0.02
    # L0 = 75.0  # in Mm
    L0 = 200.0  # in Mm
    for iframe, tframe in enumerate(tframes):
        tframe_str = str(tframe).zfill(4)
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
        print("particle number in band %d: %d" % (eband, np.sum(fband)))
        fband = np.reshape(fband, (ny, nx))
        fband = fband[:, nx//4:nx*3//4]
        print("min, max, mean, and std: %f %f %f %f" %
              (np.min(fband), np.max(fband), np.mean(fband), np.std(fband)))

        rect[0] = rect0[0] + (rect[2] + hgap) * iframe
        ax = fig.add_axes(rect)
        vmin = plot_config['nmin']
        vmax = plot_config['nmax']
        sizes = [0.5*mhd_config["xmin"][0], 0.5*mhd_config["xmax"][0],
                 mhd_config["ymin"][0], mhd_config["ymax"][0]]
        sizes = np.asarray(sizes) * L0
        img = ax.imshow(fband + 0.1*vmin, cmap=plt.cm.viridis, aspect='auto',
                        origin='lower', extent=sizes,
                        # vmin=vmin, vmax=vmax,
                        norm=LogNorm(vmin=vmin, vmax=vmax),
                        interpolation='bicubic')
        ax.tick_params(labelsize=10)
        if iframe > 0:
            ax.tick_params(axis='y', labelleft=False)
        else:
            ax.set_ylabel(r'$y$ (Mm)', fontsize=12)
        ax.set_xlabel(r'$x$ (Mm)', fontsize=12)
        mhd_time = mhd_config["dt_out"][0] * tframe
        tva = mhd_time * L0
        title = r'$t = ' + "{:10.1f}".format(tva) + r'\text{ s}$'
        plt.title(title, fontsize=12)

    rect[0] += rect[2] + 0.01
    rect[2] = 0.02
    rect[1] += rect[3] * 0.125
    rect[3] *= 0.75
    cbar_ax = fig.add_axes(rect)
    cbar1 = fig.colorbar(img, cax=cbar_ax, extend='both')
    cbar1.ax.tick_params(labelsize=10)
    fig_text = r'$n($' + figure_text(eband) + r'$)$'
    cbar1.ax.set_ylabel(fig_text, fontsize=12)

    fdir = '../img/nrho/' + run_name + '/'
    mkdir_p(fdir)
    fname = fdir + 'nrho_tri_frames_' + str(eband) + '.pdf'
    fig.savefig(fname, dpi=200)

    if show_plot:
        plt.show()
    else:
        plt.close()


def spatial_distribution_mhd_field(plot_config, mhd_config, show_plot=True):
    """Plot spatial distribution for one energy band and corresponding MHD fields

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    # tframes = [100, 125, 160]
    fig = plt.figure(figsize=[4, 4])
    rect0 = [0.13, 0.12, 0.4, 0.8]
    rect = np.copy(rect0)
    hgap = 0.03
    # L0 = 75.0  # in Mm
    L0 = 200.0  # in Mm

    tframe = plot_config["tframe"]
    tframe_str = str(tframe).zfill(4)
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
    print("particle number in band %d: %d" % (eband, np.sum(fband)))
    fband = np.reshape(fband, (ny, nx))
    fband = fband[:, nx//4:nx*3//4]
    print("min, max, mean, and std: %f %f %f %f" %
          (np.min(fband), np.max(fband), np.mean(fband), np.std(fband)))

    nx = mhd_config["nx"][0] + 4
    ny = mhd_config["ny"][0] + 4
    dx = mhd_config["dx"][0]
    dy = mhd_config["dy"][0]
    fpath = plot_config["mhd_run_dir"] + 'bin_data/'
    fname = fpath + 'mhd_data_' + str(tframe).zfill(4)
    mhd_fields = np.fromfile(fname, dtype=np.float32)
    mhd_fields = mhd_fields.reshape((ny, nx, 8))
    bx = mhd_fields[:, :, 4]
    by = mhd_fields[:, :, 5]
    jz = np.gradient(by, dx, axis=1) - np.gradient(bx, dy, axis=0)
    rect = np.copy(rect0)
    ax = fig.add_axes(rect)
    sizes = [0.5*mhd_config["xmin"][0], 0.5*mhd_config["xmax"][0],
             mhd_config["ymin"][0], mhd_config["ymax"][0]]
    sizes = np.asarray(sizes) * L0
    img = ax.imshow(jz[:, nx//4:nx*3//4], extent=sizes, cmap=plt.cm.seismic,
                    vmin=-500, vmax=500,
                    aspect='auto', origin='lower', interpolation='bicubic')
    ax.tick_params(labelsize=10)
    ax.set_xlabel(r'$x$/Mm', fontsize=12)
    ax.set_ylabel(r'$y$ (Mm)', fontsize=12)

    rect_cbar = np.copy(rect)
    rect_cbar[0] += 0.02
    rect_cbar[2] = 0.015
    rect_cbar[1] = 0.4
    rect_cbar[3] = 0.3
    cbar_ax = fig.add_axes(rect_cbar)
    cbar = fig.colorbar(img, cax=cbar_ax, extend="both")
    cbar.ax.tick_params(labelsize=10, color='k')
    cbar.ax.yaxis.set_tick_params(color='k')
    cbar.outline.set_edgecolor('k')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='k')
    fig_text = r'$j_z$'
    # cbar.ax.set_ylabel(fig_text, fontsize=12, color='w')
    ax.text(0.05, 0.95, fig_text, color='k', fontsize=12,
            horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    rect[0] += rect[2] + hgap
    ax = fig.add_axes(rect)
    vmin = plot_config['nmin']
    vmax = plot_config['nmax']
    sizes = [0.5*mhd_config["xmin"][0], 0.5*mhd_config["xmax"][0],
             mhd_config["ymin"][0], mhd_config["ymax"][0]]
    sizes = np.asarray(sizes) * L0
    img = ax.imshow(fband + 0.1*vmin, cmap=plt.cm.viridis, aspect='auto',
                    origin='lower', extent=sizes,
                    # vmin=vmin, vmax=vmax,
                    norm=LogNorm(vmin=vmin, vmax=vmax),
                    interpolation='bicubic')
    ax.tick_params(labelsize=10)
    ax.set_xlabel(r'$x$ (Mm)', fontsize=12)
    ax.tick_params(axis='y', labelleft=False)
    mhd_time = mhd_config["dt_out"][0] * tframe
    tva = mhd_time * L0
    title = r'$t = ' + "{:10.1f}".format(tva) + r'\text{ s}$'
    plt.suptitle(title, fontsize=12)

    rect_cbar = np.copy(rect)
    rect_cbar[0] += 0.02
    rect_cbar[2] = 0.015
    rect_cbar[1] = 0.4
    rect_cbar[3] = 0.3
    cbar_ax = fig.add_axes(rect_cbar)
    cbar = fig.colorbar(img, cax=cbar_ax, extend="both")
    cbar.ax.tick_params(labelsize=10, color='w')
    cbar.ax.yaxis.set_tick_params(color='w')
    cbar.outline.set_edgecolor('w')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')
    fig_text = r'$n($' + figure_text(eband) + r'$)$'
    # cbar.ax.set_ylabel(fig_text, fontsize=12, color='w')
    ax.text(0.05, 0.95, fig_text, color='white', fontsize=12,
            horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)


    fdir = '../img/nrho/' + run_name + '/'
    mkdir_p(fdir)
    fname = fdir + 'nrho_jz_' + str(eband) + '.pdf'
    # fig.savefig(fname, dpi=200)

    if show_plot:
        plt.show()
    else:
        plt.close()


def spatial_distribution_mhd_field_multi(plot_config, mhd_config, show_plot=True):
    """Plot spatial distribution for multi energy bands and corresponding MHD fields

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    # L0 = 75.0  # in Mm
    L0 = 200.0  # in Mm

    tframe = plot_config["tframe"]

    nx = mhd_config["nx"][0] + 4
    ny = mhd_config["ny"][0] + 4
    dx = mhd_config["dx"][0]
    dy = mhd_config["dy"][0]
    fpath = plot_config["mhd_run_dir"] + 'bin_data/'
    fname = fpath + 'mhd_data_' + str(tframe).zfill(4)
    mhd_fields = np.fromfile(fname, dtype=np.float32)
    mhd_fields = mhd_fields.reshape((ny, nx, 8))
    bx = mhd_fields[:, :, 4]
    by = mhd_fields[:, :, 5]
    jz = np.gradient(by, dx, axis=1) - np.gradient(bx, dy, axis=0)

    fig = plt.figure(figsize=[24, 6])
    rect0 = [0.04, 0.1, 0.145, 0.83]
    rect = np.copy(rect0)
    hgap = 0.015
    ax = fig.add_axes(rect)
    sizes = [mhd_config["xmin"][0], mhd_config["xmax"][0],
             mhd_config["ymin"][0], mhd_config["ymax"][0]]
    sizes = np.asarray(sizes) * L0
    # fdata = [:, nx//4:nx*3//4]
    fdata = jz
    img = ax.imshow(fdata, extent=sizes, cmap=plt.cm.seismic,
                    vmin=-500, vmax=500,
                    aspect='auto', origin='lower', interpolation='bicubic')
    ax.tick_params(labelsize=12)
    ax.set_xlabel(r'$x$/Mm', fontsize=16)
    ax.set_ylabel(r'$y$ (Mm)', fontsize=16)

    rect_cbar = np.copy(rect)
    rect_cbar[0] += 0.01
    rect_cbar[2] = 0.005
    rect_cbar[1] = 0.4
    rect_cbar[3] = 0.3
    cbar_ax = fig.add_axes(rect_cbar)
    cbar = fig.colorbar(img, cax=cbar_ax, extend="both")
    cbar.ax.tick_params(labelsize=10, color='k')
    cbar.ax.yaxis.set_tick_params(color='k')
    cbar.outline.set_edgecolor('k')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='k')
    fig_text = r'$j_z$'
    # cbar.ax.set_ylabel(fig_text, fontsize=12, color='w')
    ax.text(0.05, 0.95, fig_text, color='k', fontsize=16,
            horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    fnorms = [1, 2, 4, 32]

    tframe_str = str(tframe).zfill(4)
    run_name = plot_config["run_name"]
    fname = '../data/' + run_name + '/fxy-' + tframe_str + '_sum.dat'
    fdata = np.fromfile(fname)
    nx, = plot_config["nx"]
    ny, = plot_config["ny"]
    nbands = int(fdata[0])
    pbins_edges = fdata[1:nbands+2]
    pinit = 0.1
    pbins_edges /= pinit
    ebins_edges = pbins_edges**2
    print("data size: %d %d" % (nx, ny))
    ndp, = fdata.shape
    dists = fdata[nbands+2:].reshape([(ndp - nbands - 2) // (nx * ny), -1])
    print("Total number of particles: %d" % np.sum(dists))

    fband_ycuts = []
    bands, bande = plot_config["high_bands"]
    for iband, eband in enumerate(range(bands, bande+1)):
        fband = dists[eband, :]
        print("particle number in band %d: %d" % (eband, np.sum(fband)))
        print("min, max, mean, and std: %f %f %f %f" %
              (np.min(fband), np.max(fband), np.mean(fband), np.std(fband)))
        fband *= fnorms[iband]
        fband = np.reshape(fband, (ny, nx))
        # fband = fband[:, nx//4:nx*3//4]
        _, nxl = fband.shape
        fdata_cut = np.mean(fband[:, nxl//2-3:nxl//2+4], axis=1)
        fband_ycuts.append(fdata_cut)
        rect[0] += rect[2] + hgap
        ax = fig.add_axes(rect)
        vmin = plot_config['nmin']
        vmax = plot_config['nmax']
        img = ax.imshow(fband + 0.1*vmin, cmap=plt.cm.viridis,
                        aspect='auto', origin='lower', extent=sizes,
                        # vmin=vmin, vmax=vmax,
                        norm=LogNorm(vmin=vmin, vmax=vmax),
                        interpolation='bicubic')
        ax.tick_params(labelsize=12)
        ax.set_xlabel(r'$x$ (Mm)', fontsize=16)
        ax.tick_params(axis='y', labelleft=False)

        rect_cbar = np.copy(rect)
        rect_cbar[0] += 0.01
        rect_cbar[2] = 0.005
        rect_cbar[1] = 0.4
        rect_cbar[3] = 0.3
        cbar_ax = fig.add_axes(rect_cbar)
        cbar = fig.colorbar(img, cax=cbar_ax, extend="both")
        cbar.ax.tick_params(labelsize=10, color='w')
        cbar.ax.yaxis.set_tick_params(color='w')
        cbar.outline.set_edgecolor('w')
        plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')
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
        # cbar.ax.set_ylabel(fig_text, fontsize=12, color='w')
        ax.text(0.05, 0.95, fig_text, color='white', fontsize=12,
                horizontalalignment='left', verticalalignment='center',
                transform=ax.transAxes)

    rect[0] += rect[2] + hgap
    ax = fig.add_axes(rect)
    ygrid = np.linspace(sizes[2], sizes[3], ny)
    for iband, fband_y in enumerate(fband_ycuts):
        label = r'$' + str(iband + 1) + '$'
        ax.semilogx(fband_y, ygrid, nonposx="mask", label=label)
    ax.legend(loc=1, prop={'size': 16}, ncol=1,
              shadow=False, fancybox=False, frameon=False)
    ax.grid()
    ax.set_xlim([1E-1, 1E3])
    ax.set_xticks(np.logspace(-1, 3, num=5))
    ax.set_xlabel('Density', fontsize=16)
    ax.tick_params(axis='y', labelleft=False)
    # fname = '../data/' + run_name + '/fp-' + str(tframe).zfill(4) + '_sum.dat'
    # data = np.fromfile(fname, dtype=np.float64)
    # dsz, = data.shape
    # nbins = dsz // 2
    # pinit = 0.1
    # pmom = data[0:nbins] / pinit
    # elog = pmom**2
    # elog *= plot_config["e0"]  # to keV
    # fname = '../data/' + run_name + '/fp_local_' + tframe_str + '_sum.dat'
    # fdata = np.fromfile(fname).reshape([ny, nx, -1])
    # fpy = np.sum(fdata, axis=1)
    # Xn, Yn = np.meshgrid(elog, ygrid)
    # img = ax.pcolormesh(Xn, Yn, fpy, cmap=plt.cm.plasma,
    #                     norm=LogNorm(vmin=vmin, vmax=vmax))
    # ax.set_xlim([1E0, 1E3])
    # ax.set_xscale('log')

    ax.set_ylim(sizes[2:])
    ax.tick_params(labelsize=12)

    mhd_time = mhd_config["dt_out"][0] * tframe
    tva = mhd_time * L0
    title = r'$t = ' + "{:10.1f}".format(tva) + r'\text{ s}$'
    plt.suptitle(title, fontsize=20)


    fdir = '../img/nrho/' + run_name + '/'
    mkdir_p(fdir)
    fname = fdir + 'nrho_jz_' + str(tframe) + '.jpg'
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
    parser.add_argument('--time_loop', action="store_true", default=False,
                        help='whether to use a time loop to analyze multiple frames')
    parser.add_argument('--multi_bands', action="store_true", default=False,
                        help='whether to plot multiple energy bands')
    parser.add_argument('--tri_frames', action="store_true", default=False,
                        help='whether to plot three time frames')
    parser.add_argument('--with_mhd', action="store_true", default=False,
                        help='whether to plot with MHD fields')
    parser.add_argument('--tstart', action="store", default='0', type=int,
                        help='starting time frame')
    parser.add_argument('--tend', action="store", default='10', type=int,
                        help='ending time frame')
    parser.add_argument('--nmin', action="store", default='0.1', type=float,
                        help='minimum density')
    parser.add_argument('--nmax', action="store", default='10.0', type=float,
                        help='maximum density')
    return parser.parse_args()


def analysis_single_frames(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    if args.tri_frames:
        spatial_distribution_tri_frames(plot_config, mhd_config, show_plot=True)
    if args.with_mhd:
        if args.multi_bands:
            spatial_distribution_mhd_field_multi(plot_config, mhd_config, show_plot=True)
        else:
            spatial_distribution_mhd_field(plot_config, mhd_config, show_plot=True)
    else:
        if args.multi_bands:
            spatial_distribution_multi(plot_config, mhd_config, show_plot=True)
        else:
            spatial_distribution_band(plot_config, mhd_config, show_plot=True)


def process_input(plot_config, mhd_config, args, tframe):
    """process one time frame"""
    plot_config["tframe"] = tframe
    if args.with_mhd:
        if args.multi_bands:
            spatial_distribution_mhd_field_multi(plot_config, mhd_config, show_plot=False)
        else:
            spatial_distribution_mhd_field(plot_config, mhd_config, show_plot=False)
    else:
        if args.multi_bands:
            spatial_distribution_multi(plot_config, mhd_config, show_plot=False)
        else:
            spatial_distribution_band(plot_config, mhd_config, show_plot=False)


def analysis_multi_frames(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    tframes = range(plot_config["tmin"], plot_config["tmax"] + 1)
    if args.time_loop:
        for tframe in tframes:
            print("Time frame: %d" % tframe)
            plot_config["tframe"] = tframe
            if args.with_mhd:
                if args.multi_bands:
                    spatial_distribution_mhd_field_multi(plot_config, mhd_config, show_plot=False)
                else:
                    spatial_distribution_mhd_field(plot_config, mhd_config, show_plot=False)
            else:
                if args.multi_bands:
                    spatial_distribution_multi(plot_config, mhd_config, show_plot=False)
                else:
                    spatial_distribution_band(plot_config, mhd_config, show_plot=False)
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
    plot_config["nmin"] = args.nmin
    plot_config["nmax"] = args.nmax
    plot_config["mhd_run_dir"] = args.mhd_run_dir
    if args.multi_frames:
        analysis_multi_frames(plot_config, mhd_config, args)
    else:
        analysis_single_frames(plot_config, mhd_config, args)


if __name__ == "__main__":
    main()

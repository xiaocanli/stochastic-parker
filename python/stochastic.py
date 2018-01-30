"""
Analysis procedures for stochastic integration
"""
import argparse
import collections
import itertools
import math
import multiprocessing
import os.path
import struct
import sys

import json
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from matplotlib import rc
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D

sys.path.insert(0, '/net/scratch3/xiaocanli/mhd/python/')
from mhd_analysis_2d import read_fields_data
from momentum_energy_distributions import particle_momentum_energy_distributions
from shell_functions import *

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

font = {'family' : 'serif',
        #'color'  : 'darkred',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 24,
        }

def sum_particle_distribution(run_name, ct, mpi_size, nbins, p, elog):
    """
    """
    fp = np.zeros(nbins)
    fe = np.zeros(nbins)
    for i in range(0, mpi_size):
        fname = '../data/' + run_name + '/fp-' + str(ct).zfill(4) + '_' + \
                str(i).zfill(4) + '.dat'
        data = np.fromfile(fname, dtype=np.float64)
        fp += data[nbins:] / np.gradient(p)
        fe += data[nbins:] / np.gradient(elog)
    fname = '../data/' + run_name + '/fp-' + str(ct).zfill(4) + '_sum.dat'
    fp.tofile(fname)
    fname = '../data/' + run_name + '/fe-' + str(ct).zfill(4) + '_sum.dat'
    fe.tofile(fname)


def plot_spatial_distributions(ct, nx, ny, run_name):
    """
    """
    fig = plt.figure(figsize=[10, 10])
    xs0, ys0 = 0.1, 0.55
    w1, h1 = 0.4, 0.4
    hgap, vgap = 0.02, 0.02
    fname = '../data/' + run_name + '/fxy-' + str(ct).zfill(4) + '.dat'
    data = np.genfromtxt(fname)

    axs = [None] * 4
    ims = [None] * 4
    for i, j in itertools.product(range(2), range(2)):
        xs = xs0 + (hgap + w1) * i
        ys = ys0 - (vgap + h1) * j
        n = j * 2 + i
        f0 = data[:, n]
        f0 = np.reshape(f0, (nx, ny))
        print("min, max and mean of the data: %f %f %f %f" %
              (np.min(f0), np.max(f0), np.mean(f0), np.std(f0)))
        axs[n] = fig.add_axes([xs, ys, w1, h1])
        vmin, vmax = 0.01, 1.0
        ims[n] = axs[n].imshow(f0, cmap=plt.cm.Blues, aspect='auto',
                origin='lower',
                # vmin = vmin, vmax = vmax,
                norm=LogNorm(vmin=vmin, vmax=vmax),
                interpolation='bicubic')
        axs[n].tick_params(labelsize=16)
        el = "{%0.1f}" % (0.5 + n)**2
        eh = "{%0.1f}" % (1.5 + n)**2
        title = r'$' + el + r'\leq \varepsilon/\varepsilon_0 <' + eh + '$'
        axs[n].text(0.02, 0.9, title, color='w', fontsize=20, 
                bbox=dict(facecolor='none', alpha=1.0, edgecolor='none', pad=10.0),
                horizontalalignment='left', verticalalignment='bottom',
                transform = axs[n].transAxes)

    axs[2].set_xlabel(r'$x$', fontsize=20)
    axs[3].set_xlabel(r'$x$', fontsize=20)
    axs[0].set_ylabel(r'$y$', fontsize=20)
    axs[2].set_ylabel(r'$y$', fontsize=20)

    axs[0].tick_params(axis='x', labelbottom='off')
    axs[1].tick_params(axis='x', labelbottom='off')
    axs[1].tick_params(axis='y', labelleft='off')
    axs[3].tick_params(axis='y', labelleft='off')

    cbar_ax = fig.add_axes([xs+w1+0.01, ys, 0.02, h1*2+vgap])
    cbar1 = fig.colorbar(ims[0], cax=cbar_ax)
    cbar1.ax.tick_params(labelsize=16)

    fdir = '../img/nrho/' + run_name + '/'
    mkdir_p(fdir)
    fname = fdir + 'nrho_' + str(ct) + '.jpg'
    fig.savefig(fname, dpi=200)

    plt.close()
    # plt.show()


def plot_spatial_distributions_high_energy(ct, nx, ny, run_name, ene_bin=0):
    """Plot particle spatial distribution

    Args:
        ct: time frame
        nx, ny: the dimension of the distributions
        run_name: name of the simulation
        ene_bin: energy bin to plot
    """
    fdata = np.zeros(nx * ny * 4)
    fname = '../data/' + run_name + '/fxy-' + str(ct).zfill(4) + '_sum.dat'
    fdata = np.fromfile(fname)
        
    sz, = fdata.shape
    dists = fdata.reshape([5, sz/5])
    print(np.sum(dists))

    f0 = dists[ene_bin, :]
    print(f0.size, nx, ny)
    f0 = np.reshape(f0, (ny, nx))
    print("nx, ny: %d %d" % (nx, ny))
    print("min, max, mean, and std: %f %f %f %f" %
          (np.min(f0), np.max(f0), np.mean(f0), np.std(f0)))

    fig = plt.figure(figsize=[12, 12])
    xs, ys = 0.07, 0.1
    w1, h1 = 0.86, 0.86
    hgap, vgap = 0.02, 0.02
    ax = fig.add_axes([xs, ys, w1, h1])
    vmin, vmax = 1, 2E1
    ylim1 = 0.0
    ylim2 = 2.0
    xmin, xmax = 0, 2
    # ax.plot(np.sum(f0, axis=0), color='k')
    im = ax.imshow(f0, cmap=plt.cm.viridis, aspect='auto',
            origin='lower', extent=[xmin, xmax, ylim1, ylim2],
            # vmin=vmin, vmax=vmax,
            norm=LogNorm(vmin=vmin, vmax=vmax),
            interpolation='bicubic')
    ax.tick_params(labelsize=16)
    if ene_bin == 0:
        eh = "{%0.1f}" % (2**ene_bin * 0.75)**2
        title = r'$\varepsilon/\varepsilon_0 <' + eh + '$'
    elif ene_bin == 4:
        el = "{%0.1f}" % (2**(ene_bin - 2) * 1.5)**2
        title = r'$\varepsilon/\varepsilon_0 >' + el + '$'
    else:
        el = "{%0.1f}" % (2**(ene_bin - 1) * 0.75)**2
        eh = "{%0.1f}" % (2**(ene_bin - 1) * 1.5)**2
        title = r'$' + el + r'\leq \varepsilon/\varepsilon_0 <' + eh + '$'
    ax.text(0.02, 0.9, title, color='k', fontsize=20, 
            bbox=dict(facecolor='none', alpha=1.0, edgecolor='none', pad=10.0),
            horizontalalignment='left', verticalalignment='bottom',
            transform = ax.transAxes)

    ax.set_xlabel(r'$x$', fontsize=20)
    ax.set_ylabel(r'$y$', fontsize=20)

    cbar_ax = fig.add_axes([xs+w1+0.01, ys, 0.02, h1])
    cbar1 = fig.colorbar(im, cax=cbar_ax)
    cbar1.ax.tick_params(labelsize=16)

    fdir = '../img/nrho/' + run_name + '/'
    mkdir_p(fdir)
    fname = fdir + 'nrho_' + str(ct) + '_' + str(ene_bin) + '.jpg'
    fig.savefig(fname, dpi=200)

    # plt.close()
    # plt.show()


def spatial_distributions(ct, nx, ny, run_name):
    """
    """
    fig = plt.figure(figsize=[12, 6])
    xs0, ys0 = 0.07, 0.15
    w1, h1 = 0.4, 0.8
    hgap, vgap = 0.05, 0.02
    fname = '../data/' + run_name + '/fxy-' + str(ct).zfill(4) + '.dat'
    data = np.genfromtxt(fname)

    axs = [None] * 2
    ims = [None] * 2
    ns = [1, 2]
    xmin, xmax = 0, 2
    ymin, ymax = 0, 2
    for i in range(2):
        xs = xs0 + (hgap + w1) * i
        n = i
        f0 = data[:, ns[i]]
        f0 = np.reshape(f0, (nx, ny))
        f0[-1, :] = f0[0, :]
        f0 = np.roll(f0, nx/2, axis=0)
        print("min, max and mean of the data: %f %f %f %f" %
              (np.min(f0), np.max(f0), np.mean(f0), np.std(f0)))
        axs[n] = fig.add_axes([xs, ys0, w1, h1])
        vmin, vmax = 0.01, 1.0
        ims[n] = axs[n].imshow(f0, cmap=plt.cm.Blues, aspect='auto',
                origin='lower', extent=[xmin, xmax, ymin, ymax],
                # vmin = vmin, vmax = vmax,
                norm = LogNorm(vmin=vmin, vmax=vmax),
                interpolation='bicubic')
        axs[n].tick_params(labelsize=16)
        el = "{%0.1f}" % (0.5 + ns[i])**2
        eh = "{%0.1f}" % (1.5 + ns[i])**2
        title = r'$' + el + r'\leq \varepsilon/\varepsilon_0 <' + eh + '$'
        axs[n].text(0.02, 0.9, title, color='k', fontsize=20, 
                bbox=dict(facecolor='none', alpha=1.0, edgecolor='none', pad=10.0),
                horizontalalignment='left', verticalalignment='bottom',
                transform = axs[n].transAxes)
        axs[i].set_xlabel(r'$x$', fontsize=20)

    axs[0].set_ylabel(r'$y$', fontsize=20)
    axs[1].tick_params(axis='y', labelleft='off')

    cbar_ax = fig.add_axes([xs+w1+0.01, ys0, 0.02, h1])
    cbar1 = fig.colorbar(ims[0], cax=cbar_ax)
    cbar1.ax.tick_params(labelsize=16)

    fdir = '../img/nrho/' + run_name + '/'
    mkdir_p(fdir)
    fname = fdir + 'nrho_' + str(ct) + '.jpg'
    fig.savefig(fname, dpi=200)

    # plt.close()
    plt.show()


def plot_ptl_traj(run_dir):
    """
    """
    fname = '../data/tracked_particle_points_0007.dat'
    fh = open(fname, 'r')
    offset = 0
    tmp = np.memmap(fh, dtype='int32', mode='r', offset=offset,
                     shape=(1), order='F')
    nptl = tmp[0]
    offset = 4
    nsteps_tracked = np.memmap(fh, dtype='int32', mode='r', offset=offset,
                               shape=(nptl), order='F')
    ptl_type = np.dtype([('x', np.float64), ('y', np.float64),
                         ('p', np.float64), ('weight', np.float64),
                         ('t', np.float64), ('dt', np.float64),
                         ('split_time', np.int32), ('count_flag', np.int32),
                         ('tag', np.int32), ('nsteps_tracking', np.int32)])
    ptl_size = ptl_type.itemsize
    iptl = 99
    offset += nptl * 4 + np.sum(nsteps_tracked[:iptl]) * ptl_size
    ptl = np.memmap(fh, dtype=ptl_type, mode='r', offset=offset,
                    shape=(nsteps_tracked[iptl]), order='F')
    fh.close()

    ptl = np.array([list(item) for item in ptl])
    dx = np.diff(ptl[:-1, 0])
    dy = np.diff(ptl[:-1, 1])
    if np.any(np.abs(dx) > 0.5):
        xc_index = np.array(np.where(np.abs(dx) > 0.5)) + 1
        sz = xc_index.shape
        if sz[1] > 1:
            xc_index = np.squeeze(xc_index)
            for i in xc_index:
                shift = -2 if dx[i-1] > 0 else 2
                ptl[i:, 0] += shift
        else:
            i = xc_index[0, 0]
            shift = -2 if dx[i-1] > 0 else 2
            ptl[i:, 0] += shift
    if np.any(np.abs(dy) > 0.5):
        yc_index = np.array(np.where(np.abs(dy) > 0.5)) + 1
        sz = yc_index.shape
        if sz[1] > 1:
            yc_index = np.squeeze(yc_index)
            for i in yc_index:
                shift = -2 if dy[i-1] > 0 else 2
                ptl[i:, 1] += shift
        else:
            i = yc_index[0, 0]
            shift = -2 if dy[i-1] > 0 else 2
            ptl[i:, 1] += shift
    dp = np.gradient(ptl[:, 2])
    x = ptl[:-1, 0]
    y = ptl[:-1, 1]
    p = ptl[:-1, 2]

    xx, yy, data = read_fields_data(run_dir, 200, 'reconnection')
    xmin = np.min(xx)
    ymin = np.min(yy)
    xmax = np.max(xx)
    ymax = np.max(yy)

    fig = plt.figure(figsize=[20, 10])
    xs, ys = 0.1, 0.15
    w1, h1 = 0.4, 0.8
    gap = 0.05
    ax = fig.add_axes([xs, ys, w1, h1])
    rho = data[:, :, 0].T
    lx = xmax - xmin
    ly = ymax - ymin
    for i in range(3):
        xmin1 = xmin + (i - 1) * lx
        xmax1 = xmin1 + lx
        for j in range(3):
            ymin1 = ymin + (j - 1) * ly
            ymax1 = ymin1 + ly
            im = ax.imshow(rho, extent=[xmin1, xmax1, ymin1, ymax1],
                           aspect='auto', origin='lower')
    # plt.scatter(ptl[:-1:10, 0], ptl[:-1:10, 1], c=dp[:-1:10])
    ax.plot(x, y, color='k')
    xc = (np.min(x) + np.max(x)) * 0.5
    yc = (np.min(y) + np.max(y)) * 0.5
    ax.set_xlim([xc - 1.0, xc + 1.0])
    ax.set_ylim([yc - 1.0, yc + 1.0])

    xs += w1 + gap
    ax1 = fig.add_axes([xs, ys, w1, h1])
    ax1.plot(x, p, color='k', linewidth=2)

    plt.show()


def load_mhd_config(run_dir):
    """Load MHD simulation configuration
    """
    fname = run_dir + 'bin_data/mhd_config.dat'
    mtype = np.dtype([('dx', float), ('dy', float), ('dz', float),
                      ('xmin', float), ('ymin', float), ('zmin', float),
                      ('xmax', float), ('ymax', float), ('zmax', float),
                      ('lx', float), ('ly', float), ('lz', float),
                      ('dt_out', float),
                      ('nx', np.int32), ('ny', np.int32), ('nz', np.int32),
                      ('nxs', np.int32), ('nys', np.int32), ('nzs', np.int32),
                      ('topox', np.int32), ('topoy', np.int32), ('topoz', np.int32),
                      ('nvar', np.int32),
                      ('bcx', np.int32), ('bcy', np.int32), ('bcz', np.int32)])
    mhd_config = np.fromfile(fname, dtype=mtype)
    return mhd_config


def get_cmd_args():
    """Get command line arguments
    """
    default_run_name = 'shock_interaction/p000_b000_001'
    default_run_dir = '/net/scratch3/xiaocanli/mhd/shock/' 
    parser = argparse.ArgumentParser(description='Stochastic integration of Parker equation')
    parser.add_argument('--run_dir', action="store", default=default_run_dir,
                        help='run directory')
    parser.add_argument('--run_name', action="store", default=default_run_name,
                        help='run name')
    parser.add_argument('--multi_frames', action="store_true", default=False,
                        help='whether analyzing multiple frames')
    parser.add_argument('--ene_bin', action="store", default='1', type=int,
                        help='Energy bin')
    parser.add_argument('--edist', action="store_true", default=False,
                        help='whether to plot energy distribution')
    parser.add_argument('--sdist', action="store_true", default=False,
                        help='whether to plot spatial distribution')
    parser.add_argument('--multi', action="store_true", default=False,
                        help='whether to plot multiple frames')
    parser.add_argument('--power_test', action="store_true", default=False,
                        help='whether this is a test for power-law fitting')
    return parser.parse_args()


if __name__ == "__main__":
    args = get_cmd_args()
    run_name = args.run_name
    run_dir = args.run_dir
    with open('config/spectrum_config.json', 'r') as f:
        config = json.load(f)
    plot_config = config[run_name]
    tmax = plot_config["tmax"]
    tmin = plot_config["tmin"]
    mhd_config = load_mhd_config(run_dir)
    nreduce = plot_config["nreduce"]
    nx = mhd_config["nx"] / nreduce
    ny = mhd_config["ny"] / nreduce
    ebin = args.ene_bin
    cts = range(tmin, tmax)
    def processInput(job_id):
        print(job_id)
        ct = job_id
        if args.sdist:
            plot_spatial_distributions_high_energy(ct, nx, ny, run_name, ene_bin=ebin)
            plt.close()
    ncores = multiprocessing.cpu_count()
    if args.multi:
        Parallel(n_jobs=ncores)(delayed(processInput)(ct) for ct in cts)
    # spatial_distributions(ct, nx, ny, run_name)
    if args.edist:
        particle_momentum_energy_distributions(plot_config, args.power_test)
    if args.sdist and (not args.multi):
        plot_spatial_distributions_high_energy(tmax, nx, ny, run_name, ene_bin=ebin)
        plt.show()
    # plot_ptl_traj(run_dir)

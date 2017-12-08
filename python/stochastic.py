"""
Analysis procedures for stochastic integration
"""
import collections
import itertools
import math
import multiprocessing
import os.path
import struct
import sys

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
        # print i
        fname = '../data/' + run_name + '/fp-' + str(ct).zfill(4) + '_' + \
                str(i).zfill(4) + '.dat'
        data = np.fromfile(fname, dtype=np.float64)
        fp += data[nbins:] / np.gradient(p)
        fe += data[nbins:] / np.gradient(elog)
    fname = '../data/' + run_name + '/fp-' + str(ct).zfill(4) + '_sum.dat'
    fp.tofile(fname)
    fname = '../data/' + run_name + '/fe-' + str(ct).zfill(4) + '_sum.dat'
    fe.tofile(fname)


def plot_particle_distributions(ct, ax1, ax2, run_name, show_plot=True,
        is_power=True, **kwargs):
    """
    """
    fname = '../data/' + run_name + '/fp-' + str(ct).zfill(4) + '_0000.dat'
    # fname = '../data/fp-' + str(ct).zfill(4) + '_0000.dat'
    data = np.fromfile(fname, dtype=np.float64)
    sz, = data.shape
    nbins = sz / 2
    print nbins
    p0 = 0.1
    p = data[0:nbins] / p0
    elog = p**2
    # fp = data[nbins:] / np.gradient(p)
    # fe = data[nbins:] / np.gradient(elog)
    # for i in range(1, 128):
    #     print i
    #     fname = '../data/' + run_name + '/fp-' + str(ct).zfill(4) + '_' + \
    #             str(i).zfill(4) + '.dat'
    #     # fname = '../data/fp-' + str(ct).zfill(4) + '_' + \
    #     #         str(i).zfill(4) + '.dat'
    #     data = np.fromfile(fname, dtype=np.float64)
    #     fp += data[nbins:] / np.gradient(p)
    #     fe += data[nbins:] / np.gradient(elog)
    fname = '../data/' + run_name + '/fp-' + str(ct).zfill(4) + '_sum.dat'
    fp = np.fromfile(fname, dtype=np.float64)
    fname = '../data/' + run_name + '/fe-' + str(ct).zfill(4) + '_sum.dat'
    fe = np.fromfile(fname, dtype=np.float64)

    if 'color' in kwargs:
        color = kwargs['color']
    else:
        color = 'b'
    ax1.loglog(p, fp, linewidth=2, color=color)
    ax1.set_xlim([1E-1, 1E2])
    ax1.set_ylim([1E-3, 1E8])
    # ax1.set_xlim([3E-1, 2E1])
    # ax1.set_ylim([1E-1, 1E8])
    # ax1.set_xlim([3E-1, 1E1])
    # ax1.set_ylim([1E-1, 1E8])

    if is_power:
        pindex = -3.9
        sindex, eindex = 45, 127
        pratio = 2E7
        # pindex = -5.4
        # sindex, eindex = 50, 75
        # pratio = 2E7
        fpower = p[sindex:eindex]**pindex * pratio
        ax1.loglog(p[sindex:eindex], fpower, linewidth=2, color='k')
        power_index = "{%0.1f}" % pindex
        tname = r'$\sim p^{' + power_index + '}$'
        ax1.text(0.72, 0.7, tname, color='black', fontsize=24,
                horizontalalignment='left', verticalalignment='center',
                transform = ax1.transAxes)

    ax1.set_xlabel(r'$p/p_0$', fontdict=font, fontsize=24)
    ax1.set_ylabel(r'$f(p)$', fontdict=font, fontsize=24)
    ax1.tick_params(labelsize=20)
    
    fe /= np.sqrt(elog)
    ax2.loglog(elog, fe, linewidth=2, color=color)
    # ax2.set_xlim([1E-1, 1E3])
    # ax2.set_ylim([1E-5, 1E8])
    # ax2.set_xlim([1E-1, 2E2])
    # ax2.set_ylim([1E-1, 1E8])
    # ax2.set_xlim([1E-1, 1E3])
    # ax2.set_ylim([1E-2, 1E6])
    # ax2.set_xlim([1E-1, 2E2])
    # ax2.set_ylim([1E-2, 2E6])
    ax2.set_xlim([1E-1, 2E1])
    ax2.set_ylim([1E-2, 5E6])

    if is_power:
        # pindex = -2.95
        # pratio = 1E7
        # pindex = -4.0
        # pratio = 1E7
        pindex = -5.8
        pratio = 5E6
        power_index = "{%0.1f}" % pindex
        # pindex = -1.5
        # pratio = 1E4
        # power_index = "{%0.1f}" % pindex
        fpower = elog[sindex:eindex]**pindex * pratio
        ax2.loglog(elog[sindex:eindex], fpower, linewidth=2, color='k')
        tname = r'$\sim \varepsilon^{' + power_index + '}$'
        ax2.text(0.7, 0.7, tname, color='black', fontsize=24,
                horizontalalignment='left', verticalalignment='center',
                transform = ax2.transAxes)
        # print 'Fraction of power-law particles:', \
        #         np.sum(data[sindex:eindex, 1]) / np.sum(data[:, 1])
        # print 'Extend in energy of the power-law part:', \
        #         elog[eindex] / elog[sindex]


    ax2.set_xlabel(r'$\varepsilon/\varepsilon_0$', fontdict=font, fontsize=24)
    ax2.set_ylabel(r'$f(\varepsilon)$', fontdict=font, fontsize=24)
    ax2.tick_params(labelsize=20)
    
    if show_plot:
        plt.show()
    else:
        pass


def particle_distributions(ct, run_name, is_multi=False):
    """
    """
    xs, ys = 0.15, 0.15
    w1, h1 = 0.8, 0.8
    fig1 = plt.figure(figsize=[7, 5])
    ax1 = fig1.add_axes([xs, ys, w1, h1])
    fig2 = plt.figure(figsize=[7, 5])
    ax2 = fig2.add_axes([xs, ys, w1, h1])

    fdir = '../img/spectrum/' + run_name + '/'
    mkdir_p(fdir)

    ntp = 201
    if is_multi:
        for ct in range(0, ntp):
            color = plt.cm.jet(ct/float(ntp), 1)
            kwargs = {"color": color}
            is_power = False if ct < ntp - 1 else True
            is_power = False
            plot_particle_distributions(ct, ax1, ax2, run_name, show_plot=False,
                                        is_power=is_power, **kwargs)
        fp_name = 'fp_time_1.eps'
        fe_name = 'fe_time_1.eps'
    else:
        plot_particle_distributions(ct, ax1, ax2, run_name)
        fp_name = 'fp_' + str(ct) + '.eps'
        fe_name = 'fe_' + str(ct) + '.eps'

    fig1.savefig(fdir + fp_name)
    fig2.savefig(fdir + fe_name)

    plt.show()


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
        print "min, max and mean of the data:", np.min(f0), np.max(f0), \
                np.mean(f0), np.std(f0)
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


def plot_spatial_distributions_high_energy(ct, nx, ny, run_name):
    """
    """
    fig = plt.figure(figsize=[12, 3])
    xs, ys = 0.07, 0.19
    w1, h1 = 0.8, 0.76
    hgap, vgap = 0.02, 0.02
    ax = fig.add_axes([xs, ys, w1, h1])
    # fname = '../data/' + run_name + '/fxy-' + str(ct).zfill(4) + '.dat'
    # data = np.genfromtxt(fname)
    fdata = np.zeros(nx * ny * 4)
    mpi_size = 16
    for mpi_rank in range(mpi_size):
        fname = '../data/' + run_name + '/fxy-' + str(ct).zfill(4) + '_0001.dat'
        data = np.fromfile(fname)
        fdata += data
        
    sz, = fdata.shape
    dists = fdata.reshape([4, sz/4])

    n = 0
    f0 = dists[n, :]
    f0 = np.reshape(f0, (nx, ny))
    print "min, max and mean of the data:", np.min(f0), np.max(f0), \
            np.mean(f0), np.std(f0)
    ax = fig.add_axes([xs, ys, w1, h1])
    vmin, vmax = 1, 1E3
    x1, x2 = 0, nx/2
    ylim1 = 0.0
    ylim2 = 1.0
    xmin, xmax = 0, 2
    f0 = f0.T
    im = ax.imshow(f0[x1:x2], cmap=plt.cm.jet, aspect='auto',
            origin='lower', extent=[xmin, xmax, ylim1, ylim2],
            # vmin=vmin, vmax=vmax,
            norm=LogNorm(vmin=vmin, vmax=vmax),
            interpolation='bicubic')
    ax.tick_params(labelsize=16)
    el = "{%0.1f}" % (0.5 + n)**2
    eh = "{%0.1f}" % (1.5 + n)**2
    title = r'$' + el + r'\leq \varepsilon/\varepsilon_0 <' + eh + '$'
    ax.text(0.02, 0.8, title, color='k', fontsize=20, 
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
    fname = fdir + 'nrho_' + str(ct) + '.jpg'
    fig.savefig(fname, dpi=200)

    # plt.close()
    plt.show()


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
        print "min, max and mean of the data:", np.min(f0), np.max(f0), \
                np.mean(f0), np.std(f0)
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
    fname = '../data/tracked_particle_points_0000.dat'
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


if __name__ == "__main__":
    cmdargs = sys.argv
    if (len(cmdargs) > 1):
        ct = int(cmdargs[1])
    else:
        ct = 10
    # run_name = 'S16E5_beta001_bg0/p133_b000'
    # run_name = 'S1E4_beta001_bg05/p133_b100'
    # run_name = 'S9E4_beta001_bg00/p133_b000'
    # run_name = 'athena_S1E5_beta01/p000_b000'
    # run_name = 'athena_S1E5_beta01/p133_b000'
    run_name = 'athena_S1E5_beta01/p133_b100'
    # run_name = 'athena_S1E5_beta01_bg10/p000_b000'
    # run_name = 'athena_S1E5_beta01_bg10/p133_b000'
    # run_name = 'athena_S1E5_beta01_bg10/p133_b100'
    # run_name = ''
    # run_path = '/net/scratch3/xiaocanli/mhd/athena4.2/bin/S1E5_vl_hlld_beta01/' 
    run_path = '/net/scratch3/xiaocanli/mhd/athena4.2/bin/S1E5_vl_hlld_beta01_bg02/' 
    # nx = ny = 1536 / 4
    # nx = ny = 2048 / 4
    # cts = range(201)
    # mpi_size = 128
    # fname = '../data/' + run_name + '/fp-' + str(ct).zfill(4) + '_0000.dat'
    # data = np.fromfile(fname, dtype=np.float64)
    # sz, = data.shape
    # nbins = sz / 2
    # print nbins
    # p0 = 0.1
    # p = data[0:nbins] / p0
    # elog = p**2
    def processInput(job_id):
        print job_id
        ct = job_id
        # plot_spatial_distributions(job_id, nx, ny, run_name)
        # plot_spatial_distributions_high_energy(job_id, nx, ny, run_name)
        sum_particle_distribution(run_name, ct, mpi_size, nbins, p, elog)
    ncores = multiprocessing.cpu_count()
    # Parallel(n_jobs=ncores)(delayed(processInput)(ct) for ct in cts)
    # spatial_distributions(ct, nx, ny, run_name)
    # particle_distributions(ct, run_name, is_multi=True)
    # plot_spatial_distributions_high_energy(ct, nx, ny, run_name)
    plot_ptl_traj(run_path)

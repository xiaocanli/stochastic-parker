"""
Analysis procedures for stochastic integration
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LogNorm
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import math
import os.path
import struct
import collections
import sys
from shell_functions import *
import itertools
import multiprocessing
from joblib import Parallel, delayed

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

font = {'family' : 'serif',
        #'color'  : 'darkred',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 24,
        }


def plot_particle_distributions(ct, ax1, ax2, show_plot=True, is_power=True,
        **kwargs):
    """
    """
    fname = '../data/fp-' + str(ct).zfill(4) + '.dat'
    data = np.genfromtxt(fname)
    p = data[:, 0]
    elog = p**2
    fp = data[:, 1] / np.gradient(p)
    fe = data[:, 1] / np.gradient(elog)

    if 'color' in kwargs:
        color = kwargs['color']
    else:
        color = 'b'
    ax1.loglog(p, fp, linewidth=2, color=color)
    ax1.set_xlim([3E-1, 1E1])
    # ax1.set_ylim([1E1, 1E7])

    if is_power:
        pindex = -4.0
        sindex, eindex = 70, 95
        pratio = 1E6
        fpower = elog[sindex:eindex]**pindex * pratio
        ax1.loglog(elog[sindex:eindex], fpower, linewidth=2, color='k')
        power_index = "{%0.1f}" % pindex
        tname = r'$\sim p^{' + power_index + '}$'
        ax1.text(0.6, 0.7, tname, color='black', fontsize=24,
                horizontalalignment='left', verticalalignment='center',
                transform = ax1.transAxes)

    ax1.set_xlabel(r'$p/p_0$', fontdict=font, fontsize=24)
    ax1.set_ylabel(r'$f(p)$', fontdict=font, fontsize=24)
    ax1.tick_params(labelsize=20)
    
    ax2.loglog(elog, fe, linewidth=2, color=color)
    ax2.set_xlim([1E-1, 1E2])
    # ax2.set_ylim([1E1, 1E7])

    if is_power:
        pindex = -2.5
        sindex, eindex = 75, 127
        pratio = 5E5
        fpower = elog[sindex:eindex]**pindex * pratio
        ax2.loglog(elog[sindex:eindex], fpower, linewidth=2, color='k')
        power_index = "{%0.1f}" % pindex
        tname = r'$\sim \varepsilon^{' + power_index + '}$'
        ax2.text(0.6, 0.7, tname, color='black', fontsize=24,
                horizontalalignment='left', verticalalignment='center',
                transform = ax2.transAxes)
        print 'Fraction of power-law particles:', \
                np.sum(data[sindex:eindex, 1]) / np.sum(data[:, 1])
        print 'Extend in energy of the power-law part:', \
                elog[eindex] / elog[sindex]


    ax2.set_xlabel(r'$\varepsilon/\varepsilon_0$', fontdict=font, fontsize=24)
    ax2.set_ylabel(r'$f(\varepsilon)$', fontdict=font, fontsize=24)
    ax2.tick_params(labelsize=20)
    
    if show_plot:
        plt.show()
    else:
        pass


def particle_distributions(ct, is_multi=False):
    """
    """
    xs, ys = 0.15, 0.15
    w1, h1 = 0.8, 0.8
    fig1 = plt.figure(figsize=[7, 5])
    ax1 = fig1.add_axes([xs, ys, w1, h1])
    fig2 = plt.figure(figsize=[7, 5])
    ax2 = fig2.add_axes([xs, ys, w1, h1])

    fdir = '../img/spectrum/'
    mkdir_p(fdir)

    ntp = 109
    if is_multi:
        for ct in range(ntp):
            color = plt.cm.jet(ct/float(ntp), 1)
            kwargs = {"color": color}
            is_power = False if ct < ntp - 1 else True
            plot_particle_distributions(ct, ax1, ax2, show_plot=False,
                                        is_power=is_power, **kwargs)
        fp_name = 'fp_time.eps'
        fe_name = 'fe_time.eps'
    else:
        plot_particle_distributions(ct, ax1, ax2)
        fp_name = 'fp_' + str(ct) + '.eps'
        fe_name = 'fe_' + str(ct) + '.eps'

    fig1.savefig(fdir + fp_name)
    fig2.savefig(fdir + fe_name)

    plt.show()


def plot_spatial_distributions(ct, nx, ny):
    """
    """
    fig = plt.figure(figsize=[10, 10])
    xs0, ys0 = 0.1, 0.55
    w1, h1 = 0.4, 0.4
    hgap, vgap = 0.02, 0.02
    fname = '../data/fxy-' + str(ct).zfill(4) + '.dat'
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
        ims[n] = axs[n].imshow(f0, cmap=plt.cm.jet, aspect='auto',
                origin='lower',
                vmin = vmin, vmax = vmax,
                # norm=LogNorm(vmin=vmin, vmax=vmax),
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

    fdir = '../img/nrho/'
    mkdir_p(fdir)
    fname = fdir + 'nrho_' + str(ct) + '.jpg'
    fig.savefig(fname, dpi=300)

    plt.close()
    # plt.show()


if __name__ == "__main__":
    cmdargs = sys.argv
    if (len(cmdargs) > 1):
        ct = int(cmdargs[1])
    else:
        ct = 10
    # plot_momentrum_distribution(ct)
    # plot_spatial_distributions(ct, 1024, 1024)
    # particle_distributions(ct)
    particle_distributions(ct, is_multi=True)
    cts = range(108)
    def processInput(job_id):
        print job_id
        plot_spatial_distributions(job_id, 1024, 1024)
    ncores = multiprocessing.cpu_count()
    # Parallel(n_jobs=ncores)(delayed(processInput)(ct) for ct in cts)

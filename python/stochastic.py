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

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

font = {'family' : 'serif',
        #'color'  : 'darkred',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 24,
        }


def plot_momentrum_distribution(ct):
    """
    """
    fig = plt.figure(figsize=[7, 5])
    xs, ys = 0.15, 0.15
    w1, h1 = 0.8, 0.8
    ax = fig.add_axes([xs, ys, w1, h1])
    data = np.genfromtxt('../data/fp-0000.dat')
    p = data[:, 0]
    f = data[:, 1] / np.gradient(p)
    f /= p
    elog = p**2
    ax.loglog(elog, f, linewidth=2, color='r')
    fname = '../data/fp-' + str(ct).zfill(4) + '.dat'
    data = np.genfromtxt(fname)
    f = data[:, 1] / np.gradient(p)
    f /= p
    ax.loglog(elog, f, linewidth=2, color='b')
    ax.set_xlim([1E-1, 1E2])
    # ax.set_ylim([1E1, 1E7])

    pindex = -4.0
    sindex, eindex = 200, 800
    pratio = 1E4
    fpower = elog[sindex:eindex]**pindex * pratio
    ax.loglog(elog[sindex:eindex], fpower, linewidth=2, color='k')

    ax.set_xlabel(r'$E/E_0$', fontdict=font, fontsize=24)
    ax.set_ylabel(r'$f(E)$', fontdict=font, fontsize=24)
    ax.tick_params(labelsize=20)
    
    mkdir_p('../img/')
    fig.savefig('../img/fe.eps')

    plt.show()


def plot_spatial_distributions(ct, nx, ny):
    """
    """
    fig = plt.figure(figsize=[7, 5])
    xs, ys = 0.15, 0.15
    w1, h1 = 0.7, 0.8
    ax = fig.add_axes([xs, ys, w1, h1])
    fname = '../data/fxy-' + str(ct).zfill(4) + '.dat'
    data = np.genfromtxt(fname)
    f0 = data[:, 0]
    f0 = np.reshape(f0, (nx, ny))
    im1 = ax.imshow(f0, cmap=plt.cm.jet, aspect='auto',
            origin='lower', interpolation='bicubic')
    ax.tick_params(labelsize=16)
    ax.set_xlabel(r'$x$', fontsize=20)
    ax.set_ylabel(r'$y$', fontsize=20)
    cbar_ax = fig.add_axes([xs+w1+0.01, ys, 0.02, h1])
    cbar1 = fig.colorbar(im1, cax=cbar_ax)
    cbar1.ax.tick_params(labelsize=16)

    plt.show()


if __name__ == "__main__":
    cmdargs = sys.argv
    ct = int(cmdargs[1])
    # plot_momentrum_distribution(ct)
    plot_spatial_distributions(ct, 1024, 1024)

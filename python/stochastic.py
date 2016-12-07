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


def plot_momentrum_distribution():
    """
    """
    fig = plt.figure(figsize=[7, 5])
    xs, ys = 0.15, 0.15
    w1, h1 = 0.8, 0.8
    ax = fig.add_axes([xs, ys, w1, h1])
    data = np.genfromtxt('../data/fp-0000.dat')
    p = data[:, 0]
    f = data[:, 1] / np.gradient(p)
    ax.loglog(p, f, linewidth=2, color='r')
    data = np.genfromtxt('../data/fp-0200.dat')
    f = data[:, 1] / np.gradient(p)
    ax.loglog(p, f, linewidth=2, color='b')
    ax.set_xlim([1E-2, 1E0])
    ax.set_ylim([1E2, 1E7])

    ax.set_xlabel(r'$p$', fontdict=font, fontsize=24)
    ax.set_ylabel(r'$f(p)$', fontdict=font, fontsize=24)
    ax.tick_params(labelsize=20)
    
    mkdir_p('../img/')
    fig.savefig('../img/fp.eps')

    plt.show()


if __name__ == "__main__":
    plot_momentrum_distribution()

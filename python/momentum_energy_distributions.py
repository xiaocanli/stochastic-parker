"""
Functions for plotting particle momentum and energy distributions
"""
import argparse
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
from scipy.optimize import curve_fit

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


def func_power(x, a, b):
    """Function for fitting with power-law expression.
    """
    return b * np.power(x, a)


def plot_momentum_distribution(p, fp, ax, plot_config, **kwargs):
    """Plot particle momentum distribution

    Args:
        p: momentum bins
        fp: momentum distribution
        ax: plotting axis
        plot_config: plotting configuration
    """
    if 'color' in kwargs:
        color = kwargs['color']
    else:
        color = 'b'
    ax.loglog(p, fp, linewidth=2, color=color)
    ax.set_xlim(plot_config["xlim_p"])
    ax.set_ylim(plot_config["ylim_p"])
    npow = len(plot_config["power_low"])
    power_test = kwargs["power_test"]
    xlim_log = np.log10(plot_config["xlim_p"])
    ylim_log = np.log10(plot_config["ylim_p"])
    rangex = xlim_log[1] - xlim_log[0]
    rangey = ylim_log[1] - ylim_log[0]
    nbin, = fp.shape

    if kwargs["plot_power"]:
        for i in range(npow):
            sindex = plot_config["power_low"][i]
            if power_test:
                eindex = nbin - 1
            else:
                eindex = plot_config["power_high"][i]
            prange = p[sindex:eindex]
            popt, pcov = curve_fit(func_power, prange, fp[sindex:eindex])
            pindex = popt[0]
            if power_test:
                pconst = popt[1]
            else:
                pconst = popt[1] * 3
            fpower = func_power(prange, pindex, pconst)
            ax.loglog(prange, fpower, linewidth=2, color='k')
            ax.plot([p[sindex], p[sindex]], ax.get_ylim(),
                    color='k', linestyle='--')
            power_index = "{%0.2f}" % pindex
            tname = r'$\sim p^{' + power_index + '}$'
            midx = math.log10(prange[(eindex - sindex)//3])
            midy = math.log10(fpower[(eindex - sindex)//3]) + 0.5  # shift by one order
            xtext = (midx - xlim_log[0]) / rangex
            ytext = (midy - ylim_log[0]) / rangey
            ax.text(xtext, ytext, tname, color='black', fontsize=24,
                    horizontalalignment='left', verticalalignment='center',
                    transform = ax.transAxes)

    ax.set_xlabel(r'$p/p_0$', fontdict=font, fontsize=24)
    ax.set_ylabel(r'$f(p)p^3$', fontdict=font, fontsize=24)
    ax.tick_params(labelsize=20)


def plot_energy_distribution(elog, fe, ax, plot_config, **kwargs):
    """Plot particle energy distribution

    Args:
        elog: energy bins
        fe: energy distribution
        ax: plotting axis
        plot_config: plotting configuration
    """
    if 'color' in kwargs:
        color = kwargs['color']
    else:
        color = 'b'
    ax.loglog(elog, fe, linewidth=2, color=color)
    ax.set_xlim(plot_config["xlim_e"])
    ax.set_ylim(plot_config["ylim_e"])
    npow = len(plot_config["power_low"])
    power_test = kwargs["power_test"]
    xlim_log = np.log10(plot_config["xlim_e"])
    ylim_log = np.log10(plot_config["ylim_e"])
    rangex = xlim_log[1] - xlim_log[0]
    rangey = ylim_log[1] - ylim_log[0]
    nbin, = fe.shape

    if kwargs["plot_power"]:
        for i in range(npow):
            sindex = plot_config["power_low"][i]
            if power_test:
                eindex = nbin - 1
            else:
                eindex = plot_config["power_high"][i]
            erange = elog[sindex:eindex]
            popt, pcov = curve_fit(func_power, erange, fe[sindex:eindex])
            pindex = popt[0]
            if power_test:
                pconst = popt[1]
            else:
                pconst = popt[1] * 5
            fpower = func_power(erange, pindex, pconst)
            ax.loglog(erange, fpower, linewidth=2, color='k')
            ax.plot([elog[sindex], elog[sindex]], ax.get_ylim(),
                    color='k', linestyle='--')
            power_index = "{%0.2f}" % pindex
            tname = r'$\sim \varepsilon^{' + power_index + '}$'
            midx = math.log10(erange[(eindex - sindex)//3])
            midy = math.log10(fpower[(eindex - sindex)//3]) + 0.5  # shift by one order
            xtext = (midx - xlim_log[0]) / rangex
            ytext = (midy - ylim_log[0]) / rangey
            ax.text(xtext, ytext, tname, color='black', fontsize=24,
                    horizontalalignment='left', verticalalignment='center',
                    transform = ax.transAxes)

    ax.set_xlabel(r'$\varepsilon/\varepsilon_0$', fontdict=font, fontsize=24)
    ax.set_ylabel(r'$f(p)p$', fontdict=font, fontsize=24)
    ax.tick_params(labelsize=20)

    if kwargs["show_plot"]:
        plt.show()
    else:
        pass


def plot_particle_distributions(ct, ax1, ax2, plot_config, **kwargs):
    """Plot particle energy and momentum distributions for one time frame

    Args:
        ct: time frame
        ax1: figure axis for momentum distribution
        ax2: figure axis for energy distribution
        plot_config: plotting configuration
    """
    run_name = plot_config["run_name"]
    fname = '../data/' + run_name + '/fp-' + str(ct).zfill(4) + '_sum.dat'
    data = np.fromfile(fname, dtype=np.float64)
    sz, = data.shape
    nbins = sz // 2
    p0 = 0.1
    p = data[0:nbins] / p0
    elog = p**2
    fp = data[nbins:] * p / np.gradient(p)
    fe = fp / p**2

    plot_momentum_distribution(p, fp, ax1, plot_config, **kwargs)
    plot_energy_distribution(elog, fe, ax2, plot_config, **kwargs)


def particle_momentum_energy_distributions(plot_config, power_test=True):
    """Plot momentum and energy distributions

    Args:
        plot_config: plotting configuration
        power_test: test for power-law fitting
    """
    xs, ys = 0.15, 0.15
    w1, h1 = 0.8, 0.8
    fig1 = plt.figure(figsize=[7, 5])
    ax1 = fig1.add_axes([xs, ys, w1, h1])
    fig2 = plt.figure(figsize=[7, 5])
    ax2 = fig2.add_axes([xs, ys, w1, h1])

    fdir = '../img/spectrum/' + plot_config["run_name"] + '/'
    mkdir_p(fdir)

    tmax = plot_config["tmax"]
    tmin = plot_config["tmin"]
    ntp = tmax - tmin
    kwargs = {"color": 'r', "show_plot": False, "plot_power": False,
              "power_test": power_test}
    for ct in range(tmin, tmax):
        color = plt.cm.jet((ct-tmin)/float(ntp), 1)
        kwargs["color"] = color
        kwargs["plot_power"] = False if ct < tmax - 1 else True
        plot_particle_distributions(ct, ax1, ax2, plot_config, **kwargs)
    fp_name = 'fp_time_1.eps'
    fe_name = 'fe_time_1.eps'

    fig1.savefig(fdir + fp_name)
    fig2.savefig(fdir + fe_name)

    plt.show()


if __name__ == "__main__":
    pass

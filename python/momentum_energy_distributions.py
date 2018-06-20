#!/usr/bin/env python3
"""
Functions for analyzing particle momentum and energy distributions
"""
from __future__ import print_function

import argparse
import json
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import palettable
from scipy.optimize import curve_fit

from sde_util import mkdir_p

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]
COLORS = palettable.colorbrewer.qualitative.Set1_9.mpl_colors

FONT = {'family' : 'serif',
        'color'  : 'black',
        'weight' : 'normal',
        'size'   : 24}


def func_power(xvar, pindex, const):
    """Function for fitting with power-law expression.
    """
    return const * np.power(xvar, pindex)


def plot_momentum_distribution(pmom, fmom, ax, plot_config, **kwargs):
    """Plot particle momentum distribution

    Args:
        pmom: momentum bins
        fmom momentum distribution
        ax: plotting axis
        plot_config: plotting configuration
    """
    if 'color' in kwargs:
        color = kwargs['color']
    else:
        color = 'b'
    ax.loglog(pmom, fmom, linewidth=2, color=color)
    ax.set_xlim(plot_config["xlim_p"])
    ax.set_ylim(plot_config["ylim_p"])
    npow = len(plot_config["power_low"])
    power_test = kwargs["power_test"]
    xlim_log = np.log10(plot_config["xlim_p"])
    ylim_log = np.log10(plot_config["ylim_p"])
    rangex = xlim_log[1] - xlim_log[0]
    rangey = ylim_log[1] - ylim_log[0]
    nbin, = fmom.shape

    if kwargs["plot_power"]:
        for i in range(npow):
            sindex = plot_config["power_low"][i]
            if power_test:
                eindex = nbin - 1
            else:
                eindex = plot_config["power_high"][i]
            prange = pmom[sindex:eindex]
            popt, _ = curve_fit(func_power, prange, fmom[sindex:eindex])
            pindex = popt[0]
            if power_test:
                pconst = popt[1]
            else:
                pconst = popt[1] * 3
            fpower = func_power(prange, pindex, pconst)
            ax.loglog(prange, fpower, linewidth=2, color='k')
            ax.plot([pmom[sindex], pmom[sindex]], ax.get_ylim(),
                    color='k', linestyle='--')
            power_index = "{%0.2f}" % pindex
            tname = r'$\sim p^{' + power_index + '}$'
            midx = math.log10(prange[(eindex - sindex)//3])
            midy = math.log10(fpower[(eindex - sindex)//3]) + 0.5
            xtext = (midx - xlim_log[0]) / rangex
            ytext = (midy - ylim_log[0]) / rangey
            ax.text(xtext, ytext, tname, color='black', fontsize=24,
                    horizontalalignment='left', verticalalignment='center',
                    transform=ax.transAxes)

    ax.set_xlabel(r'$p/p_0$', fontdict=FONT, fontsize=24)
    ax.set_ylabel(r'$f(p)p^3$', fontdict=FONT, fontsize=24)
    ax.tick_params(labelsize=20)


def plot_energy_distribution(elog, fene, ax, plot_config, **kwargs):
    """Plot particle energy distribution

    Args:
        elog: energy bins
        fene: energy distribution
        ax: plotting axis
        plot_config: plotting configuration
    """
    if 'color' in kwargs:
        color = kwargs['color']
    else:
        color = 'b'
    ax.loglog(elog, fene, linewidth=2, color=color)
    ax.set_xlim(plot_config["xlim_e"])
    ax.set_ylim(plot_config["ylim_e"])
    npow = len(plot_config["power_low"])
    power_test = kwargs["power_test"]
    xlim_log = np.log10(plot_config["xlim_e"])
    ylim_log = np.log10(plot_config["ylim_e"])
    rangex = xlim_log[1] - xlim_log[0]
    rangey = ylim_log[1] - ylim_log[0]
    nbin, = fene.shape

    if kwargs["plot_power"]:
        for i in range(npow):
            sindex = plot_config["power_low"][i]
            if power_test:
                eindex = nbin - 1
            else:
                eindex = plot_config["power_high"][i]
            erange = elog[sindex:eindex]
            popt, _ = curve_fit(func_power, erange, fene[sindex:eindex])
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
            midy = math.log10(fpower[(eindex - sindex)//3]) + 0.5
            xtext = (midx - xlim_log[0]) / rangex
            ytext = (midy - ylim_log[0]) / rangey
            ax.text(xtext, ytext, tname, color='black', fontsize=24,
                    horizontalalignment='left', verticalalignment='center',
                    transform=ax.transAxes)

    ax.set_xlabel(r'$\varepsilon/\varepsilon_0$', fontdict=FONT, fontsize=24)
    ax.set_ylabel(r'$f(p)p$', fontdict=FONT, fontsize=24)
    ax.tick_params(labelsize=20)

    if kwargs["show_plot"]:
        plt.show()
    else:
        pass


def plot_particle_distributions(tframe, ax1, ax2, plot_config, **kwargs):
    """Plot particle energy and momentum distributions for one time frame

    Args:
        tframe time frame
        ax1: figure axis for momentum distribution
        ax2: figure axis for energy distribution
        plot_config: plotting configuration
    """
    run_name = plot_config["run_name"]
    fname = '../data/' + run_name + '/fp-' + str(tframe).zfill(4) + '_sum.dat'
    data = np.fromfile(fname, dtype=np.float64)
    dsz, = data.shape
    nbins = dsz // 2
    pinit = 0.1
    pmom = data[0:nbins] / pinit
    elog = pmom**2
    fmom = data[nbins:] * pmom / np.gradient(pmom)
    fene = fmom / pmom**2

    plot_momentum_distribution(pmom, fmom, ax1, plot_config, **kwargs)
    plot_energy_distribution(elog, fene, ax2, plot_config, **kwargs)


def momentum_energy_distributions(plot_config, power_test=True):
    """Plot momentum and energy distributions

    Args:
        plot_config: plotting configuration
        power_test: test for power-law fitting
    """
    rect = [0.15, 0.15, 0.8, 0.8]
    fig1 = plt.figure(figsize=[7, 5])
    ax1 = fig1.add_axes(rect)
    fig2 = plt.figure(figsize=[7, 5])
    ax2 = fig2.add_axes(rect)

    fdir = '../img/spectrum/' + plot_config["run_name"] + '/'
    mkdir_p(fdir)

    tmax = plot_config["tmax"]
    tmin = plot_config["tmin"]
    kwargs = {"color": 'r', "show_plot": False,
              "plot_power": False, "power_test": power_test}
    for tframe in range(tmin, tmax):
        color = plt.cm.jet((tframe-tmin)/float(tmax-tmin), 1)
        kwargs["color"] = color
        kwargs["plot_power"] = False if tframe < tmax - 1 else True
        plot_particle_distributions(tframe, ax1, ax2, plot_config, **kwargs)
    fp_name = 'fp_time_1.eps'
    fe_name = 'fe_time_1.eps'

    fig1.savefig(fdir + fp_name)
    fig2.savefig(fdir + fe_name)

    plt.show()


def energy_distributions_bg(config, sde_run):
    """Plot energy distributions for runs with different guide field

    Args:
        config: plotting configurations
        sde_run: stochastic integration with specific diffusion coefficient
    """
    rect = [0.15, 0.15, 0.8, 0.8]
    fig1 = plt.figure(figsize=[7, 5])
    ax = fig1.add_axes(rect)

    fdir = '../img/spectrum/spectra_bg/'
    mkdir_p(fdir)

    mhd_runs = ["S1E5_beta01_bg00",
                "S1E5_beta01_bg02",
                "S1E5_beta01_bg05",
                "S1E5_beta01_bg10"]

    run_labels = [r"$B_g=0.0$", r"$B_g=0.2$", r"$B_g=0.5$", r"$B_g=1.0$"]
    line_styles = ['--', '-.', ':']
    ax.set_prop_cycle('color', COLORS)

    plots = []
    power_plots = []
    for index, mhd_run in enumerate(mhd_runs):
        run_name = mhd_run + "/" + sde_run
        plot_config = config[run_name]
        tframe = str(plot_config["tmax"]).zfill(4)
        fname = '../data/' + run_name + '/fp-' + tframe + '_sum.dat'
        data = np.fromfile(fname, dtype=np.float64)
        dsz, = data.shape
        nbins = dsz // 2
        pinit = 0.1
        pmom = data[0:nbins] / pinit
        elog = pmom**2
        fmom = data[nbins:] * pmom / np.gradient(pmom)
        fene = fmom / pmom**2
        plot1, = ax.loglog(elog, fene, linewidth=2, label=run_labels[index])
        plots.append(plot1)
        npow = len(plot_config["power_low"])
        if index != 1:
            for i in range(npow):
                sindex = plot_config["power_low"][i]
                eindex = plot_config["power_high"][i]
                erange = elog[sindex:eindex]
                popt, _ = curve_fit(func_power, erange, fene[sindex:eindex])
                pindex = popt[0]
                pconst = popt[1] * 3
                fpower = func_power(erange, pindex, pconst)
                power_index = "{%0.2f}" % pindex
                tname = r'$\sim\varepsilon^{' + power_index + '}$'
                plot2, = ax.loglog(erange, fpower, linewidth=2,
                                   color=COLORS[index],
                                   linestyle=line_styles[i],
                                   label=tname)
                power_plots.append(plot2)

    leg1 = ax.legend(handles=plots, loc=1, prop={'size': 20}, ncol=1,
                     shadow=False, fancybox=False, frameon=False)
    ax.add_artist(leg1)
    ax.legend(handles=power_plots, loc=3, prop={'size': 20}, ncol=1,
              shadow=False, fancybox=False, frameon=False)
    ax.set_xlim([1E-1, 2E2])
    ax.set_ylim([1E-2, 1E7])
    ax.set_yticks(np.logspace(-1, 7, num=5))
    ax.set_xlabel(r'$\varepsilon/\varepsilon_0$', fontdict=FONT, fontsize=24)
    ax.set_ylabel(r'$f(\varepsilon)$', fontdict=FONT, fontsize=24)
    ax.tick_params(labelsize=20)
    fe_name = "fe_" + sde_run + ".eps"
    fig1.savefig(fdir + fe_name)
    plt.show()


def get_cmd_args():
    """Get command line arguments
    """
    default_mhd_run = 'S1E5_beta01_bg00'
    default_sde_run = 'p000_b000_001_100'
    parser = argparse.ArgumentParser(description='Momentum and energy distributions')
    parser.add_argument('--mhd_run', action="store",
                        default=default_mhd_run, help='MHD run name')
    parser.add_argument('--sde_run', action="store",
                        default=default_sde_run, help='SDE run name')
    parser.add_argument('--power_test', action="store_true", default=False,
                        help='whether this is a test for power-law fitting')
    parser.add_argument('--multi_bg', action="store_true", default=False,
                        help='whether plotting for MHD runs with different guide field')
    return parser.parse_args()


def main():
    """business logic for when running this module as the primary one!"""
    args = get_cmd_args()
    mhd_run = args.mhd_run
    sde_run = args.sde_run
    run_name = mhd_run + "/" + sde_run
    power_test = args.power_test
    with open('config/spectrum_config_10ta.json', 'r') as file_handler:
        config = json.load(file_handler)
    plot_config = config[run_name]
    if args.multi_bg:
        energy_distributions_bg(config, sde_run)
    else:
        momentum_energy_distributions(plot_config, power_test)


if __name__ == "__main__":
    main()

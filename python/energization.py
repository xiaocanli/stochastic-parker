#!/usr/bin/env python3
"""
Functions for analyzing particle energization
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

from sde_util import load_mhd_config, mkdir_p

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


def div0(a, b):
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0]

    From: http://stackoverflow.com/a/35696047/2561161

    """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        c[~np.isfinite(c)] = 0  # -inf inf NaN
    return c


def find_nearest(array,value):
    """Find in the array the element closest to the given value
    """
    idx = (np.abs(array-value)).argmin()
    return (idx, array[idx])


def energization(run_name, tframe):
    """Plot energization terms for a single run
    """
    nframes = tframe
    fname = '../data/' + run_name + '/fdpdt-' + str(tframe).zfill(4) + '_sum.dat'
    data = np.fromfile(fname, dtype=np.float64)
    dsz, = data.shape
    nbins = dsz // 3
    pinit = 0.1
    pmom = data[0:nbins] / pinit
    elog = pmom**2
    fnptl = data[nbins:nbins*2]
    fdpdt = data[nbins*2:nbins*3]

    fname = '../data/' + run_name + '/fp-' + str(tframe).zfill(4) + '_sum.dat'
    data = np.fromfile(fname, dtype=np.float64)
    ntot = data[nbins:]

    for tframe in range(40, nframes):
        fname = '../data/' + run_name + '/fdpdt-' + str(tframe).zfill(4) + '_sum.dat'
        data = np.fromfile(fname, dtype=np.float64)
        fnptl += data[nbins:nbins*2]
        fdpdt += data[nbins*2:nbins*3]

    fdpdt = div0(fdpdt, fnptl)  # divided by particle number

    rect = [0.15, 0.15, 0.8, 0.8]
    fig1 = plt.figure(figsize=[7, 5])
    ax1 = fig1.add_axes(rect)

    unit_index, _ = find_nearest(pmom, 1.5)
    neg_index = np.argmax(fdpdt[unit_index:] <= 0)
    neg_index += unit_index
    print("Starting and ending p range: %f, %f" %
          (pmom[unit_index], pmom[neg_index]))
    popt = np.polyfit(np.log10(pmom[unit_index:neg_index]),
                      np.log10(fdpdt[unit_index:neg_index]), 1)
    # popt[0] = 1.17
    # popt[1] *= 1.05
    pindex = popt[0]
    print("Index and norm: %f %f" % (popt[0], 10**popt[1]))
    fpower = 10**popt[1] * pmom**popt[0]

    power_index = "{%0.2f}" % pindex
    tname = r'$\sim p^{' + power_index + '}$'

    ax1.loglog(pmom, fdpdt, linewidth=2, color=COLORS[0])
    ax1.loglog(pmom, fpower, linewidth=2, color='k', label=tname)
    # ax1.loglog(pmom, ntot, linewidth=2, color=COLORS[0])
    ax1.set_xlim([1.0, 100])
    # ax1.set_ylim([-5, 5])
    ax1.tick_params(labelsize=16)
    ax1.set_xlabel(r'$p/p_0$', fontdict=FONT, fontsize=20)
    ax1.set_ylabel(r'$dp/dt$', fontdict=FONT, fontsize=20)
    leg1 = ax1.legend(loc=4, prop={'size': 20}, ncol=1,
                      shadow=False, fancybox=False, frameon=False)

    fdir = '../img/dpdt/' + run_name + '/'
    mkdir_p(fdir)

    plt.show()


def energization_multi(run_config):
    """Plot energization terms for a multiple runs

    Args:
        run_config: configurations for these runs, including
            tmin, tmax, mhd_run, sde_runs, run_labels
    """
    nrun = len(run_config["sde_runs"])
    rect = [0.15, 0.15, 0.8, 0.8]
    fig1 = plt.figure(figsize=[7, 5])
    ax1 = fig1.add_axes(rect)
    COLORS = palettable.tableau.Tableau_10.mpl_colors
    ax1.set_prop_cycle('color', COLORS)

    p1s = []
    p2s = []
    power_norms = []

    for norm, sde_run, run_label in zip(run_config["norms"],
                                        run_config["sde_runs"],
                                        run_config["run_labels"]):
        run_name = run_config["mhd_run"] + "/" + sde_run
        fname = ('../data/' + run_name + '/fdpdt-' +
                 str(run_config["tmax"]).zfill(4) + '_sum.dat')
        data = np.fromfile(fname, dtype=np.float64)
        dsz, = data.shape
        nbins = dsz // 3
        pinit = 0.1
        pmom = data[0:nbins] / pinit
        elog = pmom**2
        fnptl = data[nbins:nbins*2]
        fdpdt = data[nbins*2:nbins*3]

        for tframe in range(run_config["tmin"], run_config["tmax"]):
            fname = ('../data/' + run_name + '/fdpdt-' +
                     str(tframe).zfill(4) + '_sum.dat')
            data = np.fromfile(fname, dtype=np.float64)
            fnptl += data[nbins:nbins*2]
            fdpdt += data[nbins*2:nbins*3]

        fdpdt = div0(fdpdt, fnptl)  # divided by particle number

        if sde_run != "p000_b000_001_100":
            unit_index, _ = find_nearest(pmom, 2.5)
        else:
            unit_index, _ = find_nearest(pmom, 1.05)
        neg_index = np.argmax(fdpdt[unit_index:] <= 0)
        neg_index += unit_index
        print("Starting and ending p range: %f, %f" %
              (pmom[unit_index], pmom[neg_index]))
        popt = np.polyfit(np.log10(pmom[unit_index:neg_index]),
                          np.log10(fdpdt[unit_index:neg_index]), 1)
        pindex = popt[0]
        power_norms.append(10**popt[1]*norm)
        print("Index and norm: %f %f" % (popt[0], 10**popt[1]*norm))
        fpower = 10**popt[1] * pmom**popt[0]
        power_index = "{%0.2f}" % pindex
        tname = r'$\sim p^{' + power_index + '}$'

        p1, = ax1.loglog(pmom, fdpdt*norm, linewidth=2, label=run_label)
        p2, = ax1.loglog(pmom[:neg_index], fpower[:neg_index]*norm*2,
                         linewidth=2, color=p1.get_color(), linestyle='--',
                         label=tname)
        p1s.append(p1)
        p2s.append(p2)

    ax1.set_xlim([1.0, 100])
    ax1.set_ylim([1E-4, 1E3])
    ax1.tick_params(labelsize=16)
    ax1.set_xlabel(r'$p/p_0$', fontdict=FONT, fontsize=20)
    ax1.set_ylabel(r'$\left<dp/dt\right>$', fontdict=FONT, fontsize=20)
    leg1 = ax1.legend(handles=p1s, loc=4, prop={'size': 20}, ncol=2,
                      shadow=False, fancybox=False, frameon=False)
    ax1.add_artist(leg1)
    leg2 = ax1.legend(handles=p2s, loc=2, prop={'size': 20}, ncol=2,
                      shadow=False, fancybox=False, frameon=False)

    fdir = '../img/dpdt/'
    mkdir_p(fdir)

    fname = fdir + 'fdpdt_' + run_config["mhd_run"] + '.pdf'
    fig1.savefig(fname)

    rect = [0.15, 0.15, 0.8, 0.8]
    fig1 = plt.figure(figsize=[7, 5])
    ax1 = fig1.add_axes(rect)
    ax1.set_prop_cycle('color', COLORS)
    power_norms = np.asarray(power_norms)
    norms = np.asarray(run_config["norms"])
    ax1.scatter(norms, power_norms, c=COLORS[:nrun], s=100)
    # ax1.semilogx(norms[::-1], power_norms[::-1], linewidth=2,
    #              color='k')
    norms_new = np.logspace(-1, 2)
    pindex = 1.36
    fpower = power_norms[-1] * norms_new**pindex
    power_index = "{%0.2f}" % pindex
    tname = r'$\sim v_A^{' + power_index + '}$'
    ax1.semilogx(norms_new, fpower, linewidth=2, color='k', label=tname)
    ax1.tick_params(labelsize=16)
    ax1.set_xlim([0.5, 50])
    ax1.set_ylim([-0.1, 0.7])
    ax1.set_xlabel(r'$v_A/v_{A0}$', fontdict=FONT, fontsize=20)
    ax1.set_ylabel('Acceleration rate normalization', fontdict=FONT, fontsize=20)
    leg1 = ax1.legend(loc=2, prop={'size': 20}, ncol=1,
                      shadow=False, fancybox=False, frameon=False)
    fname = fdir + 'power_norm_' + run_config["mhd_run"] + '.pdf'
    fig1.savefig(fname)

    plt.show()


def energization_time_multi(run_config):
    """Plot energization terms for a multiple runs

    Args:
        run_config: configurations for these runs, including
            tmin, tmax, mhd_run, sde_runs, run_labels
    """
    mhd_config = load_mhd_config(run_config["mhd_run_dir"])
    COLORS = palettable.tableau.Tableau_10.mpl_colors
    tinterval = 10
    tstart = 100
    tmin = math.floor(run_config["tmin"] / tinterval) * tinterval
    tmax = math.ceil((run_config["tmax"] + 1) / tinterval) * tinterval
    # tframes = np.arange(tmin + tinterval//2, tmax + tinterval//2, tinterval)
    tframes = np.arange(tmin, tmax, tinterval)
    nframes, = tframes.shape
    pindices = np.zeros(nframes)
    y0s = np.zeros(nframes)
    start_frame = tstart // tinterval - tmin // tinterval
    tva = mhd_config["dt_out"][0] * tframes / mhd_config["ly"][0]

    rect = [0.15, 0.15, 0.8, 0.8]
    fig1 = plt.figure(figsize=[7, 5])
    ax1 = fig1.add_axes(rect)
    ax1.set_prop_cycle('color', COLORS)
    fig2 = plt.figure(figsize=[7, 5])
    ax2 = fig2.add_axes(rect)
    ax2.set_prop_cycle('color', COLORS)

    p1s = []
    p2s = []

    for norm, sde_run, run_label in zip(run_config["norms"],
                                        run_config["sde_runs"],
                                        run_config["run_labels"]):
        run_name = run_config["mhd_run"] + "/" + sde_run
        fname = ('../data/' + run_name + '/fdpdt-' +
                 str(run_config["tmax"]).zfill(4) + '_sum.dat')
        data = np.fromfile(fname, dtype=np.float64)
        dsz, = data.shape
        nbins = dsz // 3
        pinit = 0.1
        pmom = data[0:nbins] / pinit
        elog = pmom**2
        fnptl_time = np.zeros((nframes, nbins))
        fdpdt_time = np.zeros((nframes, nbins))

        for tframe in range(tstart, run_config["tmax"]+1):
            fname = ('../data/' + run_name + '/fdpdt-' +
                     str(tframe).zfill(4) + '_sum.dat')
            data = np.fromfile(fname, dtype=np.float64)
            iframe = tframe // tinterval - tmin // tinterval
            fnptl_time[iframe, :] += data[nbins:nbins*2]
            fdpdt_time[iframe, :] += data[nbins*2:nbins*3]

        for iframe in range(start_frame, nframes):
            fnptl = fnptl_time[iframe, :]
            fdpdt = fdpdt_time[iframe, :]
            fdpdt = div0(fdpdt, fnptl)  # divided by particle number

            if sde_run != "p000_b000_001_100":
                unit_index, _ = find_nearest(pmom, 2.5)
            else:
                unit_index, _ = find_nearest(pmom, 1.5)
            neg_index = np.argmax(fdpdt[unit_index:] <= 0)
            neg_index += unit_index
            print("Starting and ending p range: %f, %f" %
                  (pmom[unit_index], pmom[neg_index]))
            popt = np.polyfit(np.log10(pmom[unit_index:neg_index]),
                              np.log10(fdpdt[unit_index:neg_index]), 1)
            pindex = popt[0]
            pindices[iframe] = pindex
            y0s[iframe] = 10**popt[1]*norm
            print("Index and norm: %f %f" % (popt[0], 10**popt[1]*norm))

        p1, = ax1.plot(tva, pindices, linewidth=2, label=run_label)
        p2, = ax2.semilogy(tva, y0s, linewidth=2, label=run_label)
        p1s.append(p1)
        p2s.append(p2)

    ax1.plot([tva[start_frame], tva[-1]], [1.1, 1.1], color='k',
             linestyle='--', linewidth=0.5)
    ax1.set_xlim([tva[start_frame], tva[-1]])
    ax1.set_ylim([0, 2])
    ax1.tick_params(labelsize=16)
    ax1.set_xlabel(r'$t/\tau_A$', fontdict=FONT, fontsize=20)
    ax1.set_ylabel('Acceleration rate index', fontdict=FONT, fontsize=20)

    ax2.set_xlim([tva[start_frame], tva[-1]])
    ax2.set_ylim([1E-3, 1E1])
    ax2.tick_params(labelsize=16)
    ax2.set_xlabel(r'$t/\tau_A$', fontdict=FONT, fontsize=20)
    ax2.set_ylabel('Acceleration rate normalization', fontdict=FONT, fontsize=20)
    leg1 = ax1.legend(handles=p1s, loc=4, prop={'size': 20}, ncol=2,
                      shadow=False, fancybox=False, frameon=False)
    leg2 = ax2.legend(handles=p2s, loc=1, prop={'size': 20}, ncol=2,
                      shadow=False, fancybox=False, frameon=False)

    fdir = '../img/dpdt/'
    mkdir_p(fdir)

    fname = fdir + 'pindex_time_' + run_config["mhd_run"] + '.pdf'
    fig1.savefig(fname)

    fname = fdir + 'pnorm_time_' + run_config["mhd_run"] + '.pdf'
    fig2.savefig(fname)

    plt.show()


def get_cmd_args():
    """Get command line arguments
    """
    default_mhd_run = 'S1E5_beta01_bg00'
    default_mhd_run_dir = ('/net/scratch3/xiaocanli/mhd/guide_field_scaling/' +
                           default_mhd_run + '/')
    default_sde_run = 'p000_b000_001_100'
    parser = argparse.ArgumentParser(description='Momentum and energy distributions')
    parser.add_argument('--mhd_run', action="store",
                        default=default_mhd_run, help='MHD run name')
    parser.add_argument('--mhd_run_dir', action="store",
                        default=default_mhd_run_dir, help='MHD run directory')
    parser.add_argument('--sde_run', action="store",
                        default=default_sde_run, help='SDE run name')
    parser.add_argument('--tframe', action="store", default='100', type=int,
                        help='current time frame')
    parser.add_argument('--multi_runs', action="store_true", default=False,
                        help='whether doing analysis for multiple runs')
    return parser.parse_args()


def main():
    """business logic for when running this module as the primary one!"""
    args = get_cmd_args()
    mhd_run = args.mhd_run
    sde_run = args.sde_run
    run_name = mhd_run + "/" + sde_run
    with open('config/spectrum_config_10ta.json', 'r') as file_handler:
        config = json.load(file_handler)
    if args.multi_runs:
        sde_runs = ["p000_b000_00003_100",
                    "p000_b000_0001_100",
                    "p000_b000_0003_100",
                    "p000_b000_001_100"]
        norms = [100.0/3, 10.0, 10.0/3, 1]
        run_labels = [r"$\kappa=0.0003\kappa_0$",
                      r"$\kappa=0.001\kappa_0$",
                      r"$\kappa=0.003\kappa_0$",
                      r"$\kappa=0.01\kappa_0$"]
        run_config = {"tmin": 40,
                      "tmax": args.tframe,
                      "mhd_run": mhd_run,
                      "mhd_run_dir": args.mhd_run_dir,
                      "sde_runs": sde_runs,
                      "norms": norms,
                      "run_labels": run_labels}
        # energization_multi(run_config)
        energization_time_multi(run_config)
    else:
        energization(run_name, args.tframe)


if __name__ == "__main__":
    main()

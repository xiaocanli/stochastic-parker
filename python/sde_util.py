"""
Utility functions
"""
import errno
import math
import os

import matplotlib as mpl
import matplotlib.pylab as plt
import numpy as np
from matplotlib import rc
from scipy.constants import physical_constants

mpl.rcParams['contour.negative_linestyle'] = 'solid'
rc('font', **{'family': 'serif', 'serif': ['Palatino']})
rc('text', usetex=True)
mpl.rcParams["text.latex.preamble"] = \
        (r"\usepackage{amsmath, bm}" +
         r"\DeclareMathAlphabet{\mathsfit}{\encodingdefault}{\sfdefault}{m}{sl}" +
         r"\SetMathAlphabet{\mathsfit}{bold}{\encodingdefault}{\sfdefault}{bx}{sl}" +
         r"\newcommand{\tensorsym}[1]{\bm{\mathsfit{#1}}}")

def mkdir_p(path):
    """Create directory recursively
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def load_mhd_config(run_dir):
    """Load MHD simulation configuration

    The configuration might be different from that for the original MHD
    simulation after reorganizing the MHD data.

    Arguments:
        run_dir (string): MHD run directory

    Returns:
        mhd_config (numpy ndarray): MHD configuration
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


def get_mhd_info(mhd_run_dir, config_name):
    """Get MHD information from the configuration file

    Arguments:
        mhd_run_dir (string): MHD run directory
        config_name (string): MHD simulation configuration file name
    """
    with open(mhd_run_dir + "/" + config_name) as f:
        contents = f.readlines()
    f.close()
    mhd_info = {}
    for line in contents:
        if "<" in line and ">" in line and "<" == line[0]:
            block_name = line[1:line.find(">")]
            mhd_info[block_name] = {}
        else:
            if line[0] != "#" and "=" in line:
                line_splits = line.split("=")
                tail = line_splits[1].split("\n")
                data = tail[0].split("#")
                ltmp = line_splits[0].strip()
                try:
                    mhd_info[block_name][ltmp] = float(data[0])
                except ValueError:
                    mhd_info[block_name][ltmp] = data[0].strip()
    return mhd_info


def find_nearest(array, value):
    """Find nearest value in an array
    """
    idx = (np.abs(array-value)).argmin()
    return (idx, array[idx])


def div0(a, b):
    """ ignore / 0, div0( [-1, 0, 1], 0 ) -> [0, 0, 0]

    From: http://stackoverflow.com/a/35696047/2561161

    """
    with np.errstate(divide='ignore', invalid='ignore'):
        c = np.true_divide(a, b)
        c[~np.isfinite(c)] = 0  # -inf inf NaN
    return c


def calc_va_cgs(ni, b0):
    """Calculate the Alfven speed

    Args:
        ni: density (cm^-3)
        b0: magnetic field (Gauss)
    """
    mi = 1.6726219E-24  # proton mass in gram
    va = b0 / math.sqrt(4 * math.pi * ni * mi)  # CGS in cm/s
    return va


def read_momentum_dist(run_name, tframe=0):
    """read momentum distribution

    Arguments:
        run_name (string): SDE run name
        tframe (integer, optional): time frame

    Returns:
        pbins (1D numpy array): momentum bins
        pbins_edge (1D numpy array): momentum bins edges
        fp (1D numpy array): momentum distribution
    """
    fname = '../data/' + run_name + '/fp-' + str(tframe).zfill(4) + '_sum.dat'
    data = np.fromfile(fname, dtype=np.float64)
    dsz, = data.shape
    nbins = dsz // 2
    pbins_edge = data[:nbins]
    pbins = 0.5 * (pbins_edge[1:] + pbins_edge[:-1])
    fp = data[nbins:-1] / np.diff(pbins_edge)
    fp /= pbins**2
    return (pbins, pbins_edge, fp)


def get_energy_spectrum(sde_run_config, tframe=0):
    """get energy spectrum from momentum distribution

    Arguments:
        sde_run_config (dict): SDE run configuration
        tframe (integer, optional): time frame

    Returns:
        ebins_kev (1D numpy array): energy bins
        fe (1D numpy array): energy distribution
    """
    run_name = sde_run_config["run_name"]
    pbins, pbins_edge, fp = read_momentum_dist(run_name, tframe)
    dnptl_bins = fp * pbins**2
    dpbins = np.diff(pbins_edge)  # the bins sizes
    dnptl_bins *= dpbins  # number of particles in each bin
    ebins_kev, debins_kev = get_ebins_kev(sde_run_config)
    fe = dnptl_bins / debins_kev
    return (ebins_kev, fe)


def get_ebins_kev(sde_run_config):
    """get energy bins in kev

    Arguments:
        sde_run_config (dict): SDE run configuration
    """
    if "species" in sde_run_config:
        species = sde_run_config["species"]
    else:
        species = "electron"
    mtext = species + " mass energy equivalent in MeV"
    rest_ene_kev = physical_constants[mtext][0] * 1E3
    e0_kev = sde_run_config["e0"]
    gamma = e0_kev / rest_ene_kev + 1.0
    p0_mc = math.sqrt(gamma**2 - 1)  # p0 in terms of mc
    tmin = sde_run_config["tmin"]
    run_name = sde_run_config["run_name"]
    pbins, pbins_edge, _ = read_momentum_dist(run_name, tmin)
    pinit = 0.1  # commonly used values
    pbins_mc = pbins * p0_mc / pinit
    pbins_edge_mc = pbins_edge * p0_mc / pinit
    gamma_bins = np.sqrt(pbins_mc**2 + 1)
    gamma_bins_edge = np.sqrt(pbins_edge_mc**2 + 1)
    ebins_kev = (gamma_bins - 1) * rest_ene_kev
    ebins_edge_kev = (gamma_bins_edge - 1) * rest_ene_kev
    debins_kev = np.diff(ebins_edge_kev)

    return (ebins_kev, debins_kev)


def plot_energy_spectrum(elog, fene, ax, plot_config, **kwargs):
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
    fene[fene == 0.0] = np.nan
    ax.loglog(elog, fene, linewidth=1, color=color)
    ax.set_xlim(plot_config["xlim_e"])
    ax.set_ylim(plot_config["ylim_e"])
    npow = len(plot_config["power_low"])
    xlim_log = np.log10(plot_config["xlim_e"])
    ylim_log = np.log10(plot_config["ylim_e"])
    rangex = xlim_log[1] - xlim_log[0]
    rangey = ylim_log[1] - ylim_log[0]
    nbin, = fene.shape

    if kwargs["plot_power"]:
        for i in range(npow):
            sindex = plot_config["power_low"][i]
            eindex = plot_config["power_high"][i]
            erange = elog[sindex:eindex]
            popt, _ = curve_fit(func_power, erange, fene[sindex:eindex])
            pindex = popt[0]
            pconst = popt[1] * 2
            fpower = func_power(erange, pindex, pconst)
            eindex1, _ = find_nearest(elog, 150.0)
            erange1 = elog[sindex:eindex1 + 1]
            fpower1 = func_power(erange1, -3.6, pconst)
            fpower1 *= fpower[0] / fpower1[0]
            erange2 = elog[eindex1:]
            fpower2 = func_power(erange2, -5.8, pconst)
            fpower2 *= fpower1[-1] / fpower2[0]
            ax.loglog(erange, fpower, linewidth=2, color='k')
            ax.loglog(erange1, fpower1, linewidth=1, color='k', linestyle="-")
            ax.loglog(erange2, fpower2, linewidth=1, color='k', linestyle="-")
            ax.plot([elog[sindex], elog[sindex]],
                    ax.get_ylim(),
                    color='k',
                    linestyle='--')
            power_index = "{%0.2f}" % pindex
            tname = r'$\sim \varepsilon^{' + power_index + '}$'
            midx = math.log10(erange[(eindex - sindex) // 3])
            midy = math.log10(fpower[(eindex - sindex) // 3]) + 0.5
            xtext = (midx - xlim_log[0]) / rangex
            ytext = (midy - ylim_log[0]) / rangey
            ax.text(xtext,
                    ytext,
                    tname,
                    color='black',
                    fontsize=20,
                    horizontalalignment='left',
                    verticalalignment='center',
                    transform=ax.transAxes)


if __name__ == "__main__":
    pass

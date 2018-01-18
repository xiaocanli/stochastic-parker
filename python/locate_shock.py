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
from struct import calcsize, unpack

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from matplotlib import rc
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.mplot3d import Axes3D

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


def fread(f, fmt):
    u = unpack(fmt, f.read(calcsize(fmt)))
    if (len(fmt)==1):
        u = u[0]
    return u 


def read_header(fh):
    """Read file header
    """
    header = {}
    header['nx4']  = fread(fh, 4*"i")
    header['time'] = fread(fh, 1*"f")
    header['bbox'] = fread(fh, 4*"f")
    return header


def read_fields_data_sli(run_dir, tframe):
    """Reader fields data at time frame tframe

    The MHD simulation is done by Shengtai Li's MHD code
    """
    print("Time frame = %d" % tframe)
    ext = "%04d" % tframe
    filename = "%s%s" % (run_dir + "bin_out", ext)
    f = open(filename)
    header = read_header(f)
    nvar, nx, ny, nprocs = header['nx4']
    bbox = header['bbox']

    lx = bbox[1] - bbox[0]
    ly = bbox[3] - bbox[2]
    dx = lx / nx
    dy = ly / ny
    x = bbox[0] + (np.arange(0, nx) + 0.5) * dx
    y = bbox[2] + (np.arange(0, ny) + 0.5) * dy
    xx, yy = np.meshgrid(x, y)
    data = np.zeros((nx, ny, nvar), dtype=np.float32, order='F')

    offset = 36  # Header size (5 float and 4 integer)
    n4 = np.zeros(4, dtype='int32')
    print("Reading data ...")
    for i in range(0, nprocs):
        n4 = np.memmap(filename, dtype='int32', mode='r',
                       offset=offset, shape=(4))
        ix  = n4[0]   # starting x-pos
        iy  = n4[1]   # starting y-pos
        nx1 = n4[2]   # number cell in x
        ny1 = n4[3]   # number cell in y
        offset += 16
        data[ix:ix+nx1, iy:iy+ny1, :] = np.memmap(filename, dtype='float32',
                mode='r', offset=offset, shape=(nx1, ny1, nvar), order='F')
        offset += nx1 * ny1 * nvar * 4

    print "  done."
    f.close()
    return (xx, yy, data)


def locate_shock_position(run_dir, tframe):
    """
    Locate shock position in a 2D simulation
    """
    # data is saved as rho, p, vx, vy, bx, by, ...
    xx, yy, data = read_fields_data_sli(run_dir, tframe)
    x = xx[0, :]
    y = yy[:, 0]
    xmin, xmax = x[0], x[-1]
    ymin, ymax = y[0], y[-1]
    rho = data[:, :, 0].T
    p = data[:, :, 1].T
    vx = data[:, :, 2].T
    vy = data[:, :, 3].T
    bx = data[:, :, 4].T
    by = data[:, :, 5].T
    fdata = vx
    divy_fdata, divx_fdata = np.gradient(fdata)
    ixs_div_min = np.argmax(np.abs(divx_fdata), axis=1)

    fig = plt.figure(figsize=[12, 10])
    xs, ys = 0.08, 0.08
    w1, h1 = 0.8, 0.86
    gap = 0.02
    ax = fig.add_axes([xs, ys, w1, h1])
    im = ax.imshow(divx_fdata, extent=[xmin, xmax, ymin, ymax],
                   aspect='auto', origin='lower')
    ax.plot(x[ixs_div_min], y, color='w')
    xs += w1 + gap
    cbar_ax = fig.add_axes([xs, ys, 0.03, h1])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=16)
    ax.tick_params(labelsize=16)
    ax.set_xlabel(r'$x$', fontsize=20)
    ax.set_ylabel(r'$y$', fontsize=20)

    fdir = '../img/shock_loc/'
    mkdir_p(fdir)
    fname = fdir + 'shock_loc_' + str(tframe) + '.jpg'
    fig.savefig(fname, dpi=200)

    plt.close()
    # plt.show()


if __name__ == "__main__":
    run_dir = '/net/scratch3/xiaocanli/mhd/shock/bin_data/'
    tframe = 40
    # locate_shock_position(run_dir, tframe)
    def processInput(job_id):
        print(job_id)
        tframe = job_id
        locate_shock_position(run_dir, tframe)
    ncores = multiprocessing.cpu_count()
    ncores = 10
    cts = range(127)
    Parallel(n_jobs=ncores)(delayed(processInput)(ct) for ct in cts)

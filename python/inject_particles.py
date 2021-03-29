#!/usr/bin/env python3
"""
Analysis procedures to test how to inject particles
"""
from __future__ import print_function

import argparse
import itertools
import json
import math
import multiprocessing

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from matplotlib.colors import LogNorm
from scipy import signal
from scipy.interpolate import interp1d

from sde_util import load_mhd_config, mkdir_p

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}"


def inject_at_jz(plot_config, mhd_config, show_plot=True):
    """inject particles where current density is large

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    L0 = 200.0  # in Mm

    tframe = plot_config["tframe"]
    nx = mhd_config["nx"][0] + 4
    ny = mhd_config["ny"][0] + 4
    dx = mhd_config["dx"][0]
    dy = mhd_config["dy"][0]
    xmin, = mhd_config["xmin"]
    xmax, = mhd_config["xmax"]
    ymin, = mhd_config["ymin"]
    ymax, = mhd_config["ymax"]
    lx_mhd = xmax - xmin
    ly_mhd = ymax - ymin
    fpath = plot_config["mhd_run_dir"] + 'bin_data/'
    fname = fpath + 'mhd_data_' + str(tframe).zfill(4)
    mhd_fields = np.fromfile(fname, dtype=np.float32)
    mhd_fields = mhd_fields.reshape((ny, nx, 8))
    bx = mhd_fields[:, :, 4]
    by = mhd_fields[:, :, 5]
    bz = mhd_fields[:, :, 6]
    jx = np.gradient(bz, dy, axis=0)
    jy = -np.gradient(bz, dx, axis=1)
    jz = np.gradient(by, dx, axis=1) - np.gradient(bx, dy, axis=0)
    absj = np.sqrt(jx**2 + jy**2 + jz**2)

    fig = plt.figure(figsize=[4, 8])
    rect0 = [0.15, 0.1, 0.8, 0.85]
    rect = np.copy(rect0)
    hgap = 0.03
    rect = np.copy(rect0)
    ax = fig.add_axes(rect)
    sizes = [mhd_config["xmin"][0], mhd_config["xmax"][0],
             mhd_config["ymin"][0], mhd_config["ymax"][0]]
    sizes = np.asarray(sizes) * L0
    xs = 0
    fdata = absj[:, xs:(nx-xs)]
    xgrid = np.linspace(sizes[0], sizes[1], nx-xs*2)
    ygrid = np.linspace(sizes[2], sizes[3], ny)
    nx_grid, = xgrid.shape
    ny_grid, = ygrid.shape
    img = ax.imshow(fdata, extent=sizes, cmap=plt.cm.viridis,
                    # norm=LogNorm(vmin=1E2, vmax=1E3),
                    vmin=4E2, vmax=1E3,
                    aspect='auto', origin='lower', interpolation='none')
    # absj_maxy = np.argmax(fdata[:, nx_grid//2])
    # fdata_maxx = np.argmax(fdata, axis=1)
    # cond = (nx_grid-fdata_maxx-fdata_maxx) < 32
    # condy = np.logical_and(ygrid > 10, ygrid < 50)
    # cond = np.logical_and(cond, condy)
    # ax.plot(xgrid[fdata_maxx][cond], ygrid[cond], color='k')
    # ax.plot(xgrid[nx_grid-fdata_maxx-1][cond], ygrid[cond], color='k')
    ax.tick_params(labelsize=12)
    ax.set_xlabel(r'$x$/Mm', fontsize=16)
    ax.set_ylabel(r'$y$ (Mm)', fontsize=16)

    rect_cbar = np.copy(rect)
    rect_cbar[0] += 0.02
    rect_cbar[2] = 0.03
    rect_cbar[1] = 0.4
    rect_cbar[3] = 0.3
    cbar_ax = fig.add_axes(rect_cbar)
    cbar = fig.colorbar(img, cax=cbar_ax, extend="both")
    cbar.ax.tick_params(labelsize=10, color='w')
    cbar.ax.yaxis.set_tick_params(color='w')
    cbar.outline.set_edgecolor('w')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')
    fig_text = r'$|\boldsymbol{j}|$'
    ax.text(0.05, 0.95, fig_text, color='w', fontsize=16,
            horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    fdir = '../data/' + plot_config["mhd_run"] + '/exhaust_boundary/'
    fname = fdir + 'xz_right_' + str(tframe) + '.dat'
    xlist_right, ylist_right = np.fromfile(fname).reshape([2, -1])
    fname = fdir + 'xz_left_' + str(tframe) + '.dat'
    xlist_left, ylist_left = np.fromfile(fname).reshape([2, -1])
    if not np.all(xlist_right > xlist_left):
        iy_max = np.argmax(ylist_right)
        xlist_right = np.copy(xlist_right[:iy_max+1])
        ylist_right = np.copy(ylist_right[:iy_max+1])
        xlist_left = lx_mhd - xlist_right
        ylist_left = np.copy(ylist_right)
    ax.plot(xlist_left*L0, ylist_left*L0, color='w',
            linewidth=1, linestyle='--')
    ax.plot(xlist_right*L0, ylist_right*L0, color='w',
            linewidth=1, linestyle='--')
    ax.set_xlim([sizes[0], sizes[1]])
    ax.set_ylim([sizes[2], sizes[3]])

    # nbins = 100
    # jbins = np.logspace(-1, 3, nbins+1)
    # jbins_mid = 0.5 * (jbins[1:] + jbins[:-1])
    # jdist, _ = np.histogram(absj, bins=jbins)
    # plt.loglog(jbins_mid, jdist)

    fdir = '../img/inject_jz/' + plot_config["mhd_run"] + '/'
    mkdir_p(fdir)
    fname = fdir + 'jz_' + str(tframe) + '.jpg'
    fig.savefig(fname, dpi=200)

    if show_plot:
        plt.show()
    else:
        plt.close()


def middle_step_rk4(x, y, nx, ny, dx, dy, bx, by):
    """Middle step of Runge-Kutta method to trace the magnetic field line

    Args:
        x, y: the coordinates of current point
        nx, ny: the dimensions of the data
        bx, by: the magnetic field arrays
    """
    x0 = (x + 1.5*dx) / dx
    y0 = (y + 1.5*dy) / dy
    ix1 = int(math.floor(x0))
    iy1 = int(math.floor(y0))
    offsetx = x0 - ix1
    offsety = y0 - iy1
    ix2 = ix1 + 1
    iy2 = iy1 + 1
    v1 = (1.0 - offsetx) * (1.0 - offsety)
    v2 = offsetx * (1.0 - offsety)
    v3 = offsetx * offsety
    v4 = (1.0 - offsetx) * offsety
    bx0 = (bx[iy1, ix1] * v1 + bx[iy1, ix2] * v2 +
           bx[iy2, ix2] * v3 + bx[iy2, ix1] * v4)
    by0 = (by[iy1, ix1] * v1 + by[iy1, ix2] * v2 +
           by[iy2, ix2] * v3 + by[iy2, ix1] * v4)
    absB = math.sqrt(bx0**2 + by0**2)
    deltax1 = bx0 / absB
    deltay1 = by0 / absB
    return (deltax1, deltay1)


def trace_field_line(bx, by, x0, y0, lx, ly):
    """Tracer magnetic field line

    Args:
        bx, by: magnetic field
        x0, y0: starting point of field line tracing
        lx, ly: length along x- and y- direction
    """
    ny, nx = bx.shape
    dx = lx / nx
    dy = ly / ny
    lx_mhd = (nx - 4) * dx
    ly_mhd = (ny - 4) * dy
    xmin, xmax = 0, lx_mhd
    ymin, ymax = 0, ly_mhd

    nstep = 0
    deltas = math.sqrt(dx**2 + dy**2) * 0.1
    hds = deltas * 0.5
    total_lengh = 0
    xs = x0
    ys = y0
    xlist = [xs]
    ylist = [ys]
    x, y = xs, ys
    xcond = x >= xmin and x <= xmax
    ycond = y >= ymin and y <= ymax
    while xcond and ycond and total_lengh < 5 * ly:
        deltax1, deltay1 = middle_step_rk4(x, y, nx, ny, dx, dy, bx, by)
        x1 = x + deltax1 * hds
        y1 = y + deltay1 * hds
        deltax2, deltay2 = middle_step_rk4(x1, y1, nx, ny, dx, dy, bx, by)
        x2 = x + deltax2 * hds
        y2 = y + deltay2 * hds
        deltax3, deltay3 = middle_step_rk4(x2, y2, nx, ny, dx, dy, bx, by)
        x3 = x + deltax3 * deltas
        y3 = y + deltay3 * deltas
        deltax4, deltay4 = middle_step_rk4(x3, y3, nx, ny, dx, dy, bx, by)
        x += deltas/6 * (deltax1 + 2*deltax2 + 2*deltax3 + deltax4)
        y += deltas/6 * (deltay1 + 2*deltay2 + 2*deltay3 + deltay4)
        total_lengh += deltas
        xlist.append(x)
        ylist.append(y)
        nstep += 1
        length = math.sqrt((x-xs)**2 + (y-ys)**2)
        if length < dx and nstep > 20:
            break
        xcond = x >= xmin and x <= xmax
        ycond = y >= ymin and y <= ymax

    _, _ = middle_step_rk4(xs, ys, nx, ny, dx, dy, bx, by)
    xlist = np.asarray(xlist)
    ylist = np.asarray(ylist)
    return (xlist, ylist)


def get_exhaust_boundary(plot_config, mhd_config, show_plot=True):
    """Get the exhaust boundary for the reconnection layer

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    tframe = plot_config["tframe"]
    nx = mhd_config["nx"][0] + 4
    ny = mhd_config["ny"][0] + 4
    dx = mhd_config["dx"][0]
    dy = mhd_config["dy"][0]
    fpath = plot_config["mhd_run_dir"] + 'bin_data/'
    fname = fpath + 'mhd_data_' + str(tframe).zfill(4)
    mhd_fields = np.fromfile(fname, dtype=np.float32)
    mhd_fields = mhd_fields.reshape((ny, nx, 8))
    bx = mhd_fields[:, :, 4]
    by = mhd_fields[:, :, 5]
    bz = mhd_fields[:, :, 6]

    xmin, = mhd_config["xmin"]
    xmax, = mhd_config["xmax"]
    ymin, = mhd_config["ymin"]
    ymax, = mhd_config["ymax"]
    nx_mhd, = mhd_config["nx"]
    ny_mhd, = mhd_config["ny"]
    lx_mhd = xmax - xmin
    ly_mhd = ymax - ymin
    dx = lx_mhd / nx_mhd
    dy = ly_mhd / nx_mhd
    xmin_grid = xmin - 2 * dx
    xmax_grid = xmax + 2 * dx
    ymin_grid = ymin - 2 * dy
    ymax_grid = ymax + 2 * dy
    lx_grid = xmax_grid - xmin_grid
    ly_grid = ymax_grid - ymin_grid
    x1, x2 = 0.5 * (xmin + xmax), xmax
    while (x2 - x1) > 0.1 * dx:
        xmid = (x1 + x2) * 0.5
        # print("Starting x-position: %f" % xmid)
        xlist, ylist = trace_field_line(bx, by, xmid, 0, lx_grid, ly_grid)
        if np.any(xlist < 0.5*lx_mhd) and ylist.max() < 0.25*ly_mhd:
            x1 = xmid
        else:
            x2 = xmid

    xlist, ylist = trace_field_line(bx, by, x2, 0, lx_grid, ly_grid)

    fdir = '../data/' + plot_config["mhd_run"] + '/exhaust_boundary/'
    mkdir_p(fdir)

    xz = np.asarray([xlist, ylist])
    fname = fdir + 'xz_right_' + str(tframe) + '.dat'
    xz.tofile(fname)
    xz = np.asarray([lx_mhd - xlist, ylist])
    fname = fdir + 'xz_left_' + str(tframe) + '.dat'
    xz.tofile(fname)

    fig = plt.figure(figsize=[4, 8])
    rect0 = [0.15, 0.1, 0.8, 0.85]
    rect = np.copy(rect0)
    hgap = 0.03
    rect = np.copy(rect0)
    ax = fig.add_axes(rect)
    sizes = [xmin_grid, xmax_grid, ymin_grid, ymax_grid]
    L0 = 200  # in Mm
    sizes = np.asarray(sizes) * L0
    jx = np.gradient(bz, dy, axis=0)
    jy = -np.gradient(bz, dx, axis=1)
    jz = np.gradient(by, dx, axis=1) - np.gradient(bx, dy, axis=0)
    absj = np.sqrt(jx**2 + jy**2 + jz**2)
    img = ax.imshow(absj, extent=sizes, cmap=plt.cm.viridis,
                    norm=LogNorm(vmin=1E1, vmax=1E3),
                    # vmin=0, vmax=1E3,
                    aspect='auto', origin='lower', interpolation='none')
    xlist *= L0
    ylist *= L0
    p1, = ax.plot(xlist, ylist, color='w', linewidth=1)
    ax.plot(lx_mhd*L0 - xlist, ylist, color=p1.get_color(), linewidth=1)
    ax.tick_params(labelsize=12)
    ax.set_xlabel(r'$x$/Mm', fontsize=16)
    ax.set_ylabel(r'$y$ (Mm)', fontsize=16)

    if show_plot:
        plt.show()
    else:
        plt.close()


def calc_deltab(plot_config, mhd_config):
    """Calculate the magnetic fluctuation

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    L0 = 200.0  # in Mm

    nghost = 2
    tframe = plot_config["tframe"]
    nx = mhd_config["nx"][0] + 4
    ny = mhd_config["ny"][0] + 4
    xmin, = mhd_config["xmin"]
    xmax, = mhd_config["xmax"]
    ymin, = mhd_config["ymin"]
    ymax, = mhd_config["ymax"]
    nx_mhd, = mhd_config["nx"]
    ny_mhd, = mhd_config["ny"]
    lx_mhd = xmax - xmin
    ly_mhd = ymax - ymin
    dx = lx_mhd / nx_mhd
    dy = ly_mhd / ny_mhd
    xmin_grid = xmin - nghost * dx
    xmax_grid = xmax + nghost * dx
    ymin_grid = ymin - nghost * dy
    ymax_grid = ymax + nghost * dy
    lx_grid = xmax_grid - xmin_grid
    ly_grid = ymax_grid - ymin_grid
    xmhd = np.linspace(0.5*dx, lx_mhd-0.5*dx, nx_mhd)
    ymhd = np.linspace(0.5*dy, ly_mhd-0.5*dy, ny_mhd)
    xgrid = np.linspace(-(nghost-0.5)*dx, lx_mhd+(nghost-0.5)*dx, nx)
    ygrid = np.linspace(-(nghost-0.5)*dy, ly_mhd+(nghost-0.5)*dy, ny)

    fpath = plot_config["mhd_run_dir"] + 'bin_data/'
    fname = fpath + 'mhd_data_' + str(tframe).zfill(4)
    mhd_fields = np.fromfile(fname, dtype=np.float32)
    mhd_fields = mhd_fields.reshape((ny, nx, 8))
    bx = mhd_fields[:, :, 4]
    by = mhd_fields[:, :, 5]
    bz = mhd_fields[:, :, 6]
    jx = np.gradient(bz, dy, axis=0)
    jy = -np.gradient(bz, dx, axis=1)
    jz = np.gradient(by, dx, axis=1) - np.gradient(bx, dy, axis=0)
    absj = np.sqrt(jx**2 + jy**2 + jz**2)

    fdir = '../data/' + plot_config["mhd_run"] + '/exhaust_boundary/'
    fname = fdir + 'xz_right_' + str(tframe) + '.dat'
    xlist_right, ylist_right = np.fromfile(fname).reshape([2, -1])
    fname = fdir + 'xz_left_' + str(tframe) + '.dat'
    xlist_left, ylist_left = np.fromfile(fname).reshape([2, -1])
    if not np.all(xlist_right > xlist_left):
        iy_max = np.argmax(ylist_right)
        xlist_right = np.copy(xlist_right[:iy_max+1])
        ylist_right = np.copy(ylist_right[:iy_max+1])
        xlist_left = lx_mhd - xlist_right
        ylist_left = np.copy(ylist_right)

    ix_top = len(xlist_right) - np.argmax(xlist_right[::-1] - xlist_left[::-1] > 0.1)
    ix = np.argmin(xlist_right[:ix_top] - xlist_left[:ix_top])
    cond1 = (xlist_right[:ix] - xlist_left[:ix]) < 0.01
    ips = np.argmin(cond1)
    cond2 = (xlist_right[ix:] - xlist_left[ix:]) < 0.04
    ipe = np.argmax(cond2) + ix
    cond = np.concatenate((cond1, cond2), axis=0)
    cond = np.logical_and(cond, ylist_right < ylist_left.max()*0.8)
    ys = ylist_right[cond][0]
    ye = ylist_right[cond][-1]
    fleft = interp1d(ylist_left, xlist_left)
    fright = interp1d(ylist_right, xlist_right)
    xlist_left_mhd = np.zeros(ny_mhd)
    xlist_right_mhd = np.zeros(ny_mhd)
    sep = np.zeros(ny_mhd)
    if ymhd[-1] > ylist_left[-1]:
        iy_max = np.argmax(ymhd > ylist_left.max())
        xlist_left_mhd[:iy_max] = fleft(ymhd[:iy_max])
        xlist_right_mhd[:iy_max] = fright(ymhd[:iy_max])
    else:
        xlist_left_mhd = fleft(ymhd)
        xlist_right_mhd = fright(ymhd)
        iy_max = ny_mhd - 1

    ix_left_mhd = np.floor((xlist_left_mhd - xmhd[0]) / dx + nghost).astype(int)
    ix_right_mhd = np.ceil((xlist_right_mhd - xmhd[0]) / dx + nghost).astype(int)
    # ix_right_mhd = nx - ix_left_mhd
    dist = xlist_right - xlist_left
    cond = np.logical_and(ylist_right > ys, ylist_right < ye)
    y_xpoint = ylist_right[np.argmin(dist[cond]) + np.argmax(ylist_right > ys)]
    ntranx = 20  # Number of transition cells at left and right
    ntrany_bot = 150  # Number of transition cells at top and bottom
    ntrany_top = 0.5 * (ye - y_xpoint) // dy // 2
    if ntrany_top < ntrany_bot:
        ntrany_top = ntrany_bot
    db_y = np.zeros(ny)
    db_y = 0.5 * (np.tanh((ygrid - ys)/(ntrany_bot*dy)) -
                  np.tanh((ygrid - (ye + y_xpoint) * 0.5)/(ntrany_top*dy)))
    db_y -= db_y.min()
    deltab = np.zeros([ny, nx])
    # # plt.plot(ix_left_mhd)
    # plt.plot(ix_left_mhd + ix_right_mhd)
    # plt.show()
    for iy in range(nghost, iy_max+nghost):
        ix1 = ix_left_mhd[iy-nghost]
        ix2 = ix_right_mhd[iy-nghost]
        deltab[iy, :] = 0.5 * db_y[iy] * (np.tanh((xgrid-xgrid[ix1-ntranx])/(ntranx*dx)) -
                                          np.tanh((xgrid-xgrid[ix2+ntranx])/(ntranx*dx)))
        deltab[iy, :] += 1E-2
    for iy in range(nghost):
        deltab[iy, :] = deltab[nghost]
    for iy in range(iy_max+nghost, ny):
        deltab[iy, :] = 1E-2
    # deltab[deltab < 1E-2] = 1E-2
    deltab = deltab.astype(np.float32)

    fname = fpath + 'deltab_' + str(tframe).zfill(4)
    deltab.tofile(fname)


def calc_correlation_length(plot_config, mhd_config):
    """Calculate the turbulence correlation length

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    L0 = 200.0  # in Mm
    lc_min = 2000  # in km

    nghost = 2
    tframe = plot_config["tframe"]
    nx = mhd_config["nx"][0] + nghost*2
    ny = mhd_config["ny"][0] + nghost*2
    xmin, = mhd_config["xmin"]
    xmax, = mhd_config["xmax"]
    ymin, = mhd_config["ymin"]
    ymax, = mhd_config["ymax"]
    nx_mhd, = mhd_config["nx"]
    ny_mhd, = mhd_config["ny"]
    lx_mhd = xmax - xmin
    ly_mhd = ymax - ymin
    dx = lx_mhd / nx_mhd
    dy = ly_mhd / ny_mhd
    xmin_grid = xmin - nghost * dx
    xmax_grid = xmax + nghost * dx
    ymin_grid = ymin - nghost * dy
    ymax_grid = ymax + nghost * dy
    lx_grid = xmax_grid - xmin_grid
    ly_grid = ymax_grid - ymin_grid
    xmhd = np.linspace(0.5*dx, lx_mhd-0.5*dx, nx_mhd)
    ymhd = np.linspace(0.5*dy, ly_mhd-0.5*dy, ny_mhd)
    xgrid = np.linspace(-(nghost-0.5)*dx, lx_mhd+(nghost-0.5)*dx, nx)
    ygrid = np.linspace(-(nghost-0.5)*dy, ly_mhd+(nghost-0.5)*dy, ny)

    fpath = plot_config["mhd_run_dir"] + 'bin_data/'
    fname = fpath + 'mhd_data_' + str(tframe).zfill(4)
    mhd_fields = np.fromfile(fname, dtype=np.float32)
    mhd_fields = mhd_fields.reshape((ny, nx, 8))
    bx = mhd_fields[:, :, 4]
    by = mhd_fields[:, :, 5]
    bz = mhd_fields[:, :, 6]
    absb = np.sqrt(bx**2 + by**2 + bz**2)

    # Right and left exhaust boundaries
    fdir = '../data/' + plot_config["mhd_run"] + '/exhaust_boundary/'
    fname = fdir + 'xz_right_' + str(tframe) + '.dat'
    xlist_right, ylist_right = np.fromfile(fname).reshape([2, -1])
    fname = fdir + 'xz_left_' + str(tframe) + '.dat'
    xlist_left, ylist_left = np.fromfile(fname).reshape([2, -1])
    if not np.all(xlist_right > xlist_left):
        iy_max = np.argmax(ylist_right)
        xlist_right = np.copy(xlist_right[:iy_max+1])
        ylist_right = np.copy(ylist_right[:iy_max+1])
        xlist_left = lx_mhd - xlist_right
        ylist_left = np.copy(ylist_right)

    # Find one point close to the top of the large flux rope
    ix_top = len(xlist_right) - np.argmax((xlist_right[::-1] - xlist_left[::-1]) > 0.1)
    ix = np.argmin(xlist_right[:ix_top] - xlist_left[:ix_top]) # X-point
    cond1 = (xlist_right[:ix] - xlist_left[:ix]) < 0.01
    ips = np.argmin(cond1)  # one point below the X-point and close to looptop
    cond2 = (xlist_right[ix:] - xlist_left[ix:]) < 0.04
    ipe = np.argmax(cond2) + ix # one point where the top exhaust opens up
    cond = np.concatenate((cond1, cond2), axis=0)
    cond = np.logical_and(cond, ylist_right < ylist_right.max()*0.8)

    # Bottom and top of reconnection exhaust
    ys = ylist_right[cond][0]   # Close to looptop
    ye = ylist_right[cond][-1]  # Below the large flux rope

    # Find the top of the fluxrope. If the top is out of domain, then use ny_mhd
    # And get the coordinates for the exhaust and flux rope boundaries.
    fleft = interp1d(ylist_left, xlist_left)
    fright = interp1d(ylist_right, xlist_right)
    xlist_left_mhd = np.zeros(ny_mhd)
    xlist_right_mhd = np.zeros(ny_mhd)
    if ymhd[-1] > ylist_left[-1]:
        iy_top = np.argmax(ymhd > ylist_left.max()) # Top of the large flux rope
        xlist_left_mhd[:iy_top] = fleft(ymhd[:iy_top])
        xlist_right_mhd[:iy_top] = fright(ymhd[:iy_top])
        # The exhaust boundaries intersect with left and right domain boundaries
        if (xlist_right[-1] - xlist_left[-1]) > lx_mhd:
            xlist_left_mhd[iy_top:] = xmin
            xlist_right_mhd[iy_top:] = xmax
            iy_max = ny_mhd - 1
        else:
            iy_max = iy_top
    else:
        xlist_left_mhd = fleft(ymhd)
        xlist_right_mhd = fright(ymhd)
        iy_max = ny_mhd - 1

    # x-indices of the exhaust and flux rope boundaries
    ix_left_mhd = np.floor((xlist_left_mhd - xmhd[0]) / dx + nghost).astype(int)
    ix_right_mhd = np.ceil((xlist_right_mhd - xmhd[0]) / dx + nghost).astype(int)
    # ix_right_mhd = nx - ix_left_mhd

    # separation between the left and right boundaries of the exhaust and flux rope
    iy_looptop = int(ys / dy)
    sep = np.zeros(ny_mhd)
    sep[iy_looptop:] = (xlist_right_mhd[iy_looptop:] -
                        xlist_left_mhd[iy_looptop:]) * L0 * 1E3  # in km
    sep[:iy_looptop] = (xlist_right_mhd[iy_looptop] -
                        xlist_left_mhd[iy_looptop]) * L0 * 1E3

    lc_norm = plot_config["lc_norm"]  # in km
    lc_min = plot_config["lc_min"]  # in km
    lc = np.zeros([ny, nx])
    absb_sqrt = np.sqrt(absb)
    ntranx = 20  # Number of transition cells at left and right
    for iy in range(nghost, iy_max+nghost):
        ix1 = ix_left_mhd[iy-nghost]
        ix2 = ix_right_mhd[iy-nghost]
        lc_b = lc_norm / absb[iy, ix1-4*ntranx]  # background
        if lc_b > sep[iy]:
            lc[iy, :] = (np.tanh((xgrid-xgrid[ix1-ntranx])/(ntranx*dx)) -
                         np.tanh((xgrid-xgrid[ix2+ntranx])/(ntranx*dx)))
            lc[iy, :] = lc[iy, :].max() - lc[iy, :]
            lc[iy, :] *= (lc_b - sep[iy]) / lc[iy, :].max()
            lc[iy, :] += sep[iy]
        else:
            lc[iy, :] = (np.tanh((xgrid-xgrid[ix1-ntranx])/(ntranx*dx)) -
                         np.tanh((xgrid-xgrid[ix2+ntranx])/(ntranx*dx)))
            lc[iy, :] *= (sep[iy] - lc_b) / lc[iy, :].max()
            lc[iy, :] += lc_b
        lc[iy, :ix1-4*ntranx] = lc_norm / absb[iy, :ix1-4*ntranx]
        lc[iy, ix2+4*ntranx:] = lc_norm / absb[iy, ix2+4*ntranx:]
    lc[lc < lc_min] = lc_min  # Set the floor value
    # iy_lf = int(ye / dy)  # behind the large flux rope
    # lc[iy_lf:, :] = lc_norm / absb[iy_lf:, :]
    for iy in range(nghost):
        lc[iy, :] = lc[nghost]
    for iy in range(iy_max+nghost, ny):
        if iy_max < ny_mhd - 1:
            lc[iy, :] = lc_norm / absb[iy, :]
        else:
            lc[iy, :] = lc[iy_max+nghost-1]
    # lc[lc<lc_min] = lc_min

    lc = lc.astype(np.float32)
    fname = fpath + 'lc_' + str(tframe).zfill(4)
    lc.tofile(fname)


def plot_deltab(plot_config, mhd_config, show_plot=True):
    """Plot the magnetic fluctuation

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    L0 = 200.0  # in Mm

    tframe = plot_config["tframe"]
    nx = mhd_config["nx"][0] + 4
    ny = mhd_config["ny"][0] + 4
    dx = mhd_config["dx"][0]
    dy = mhd_config["dy"][0]
    xmin, = mhd_config["xmin"]
    xmax, = mhd_config["xmax"]
    ymin, = mhd_config["ymin"]
    ymax, = mhd_config["ymax"]
    nx_mhd, = mhd_config["nx"]
    ny_mhd, = mhd_config["ny"]
    lx_mhd = xmax - xmin
    ly_mhd = ymax - ymin
    dx = lx_mhd / nx_mhd
    dy = ly_mhd / nx_mhd
    xmin_grid = xmin - 2 * dx
    xmax_grid = xmax + 2 * dx
    ymin_grid = ymin - 2 * dy
    ymax_grid = ymax + 2 * dy
    lx_grid = xmax_grid - xmin_grid
    ly_grid = ymax_grid - ymin_grid
    xmhd = np.linspace(0.5*dx, lx_mhd-0.5*dx, nx_mhd)
    ymhd = np.linspace(0.5*dy, ly_mhd-0.5*dy, ny_mhd)
    xgrid = np.linspace(-1.5*dx, lx_mhd+1.5*dx, nx)
    ygrid = np.linspace(-1.5*dy, ly_mhd+1.5*dy, ny)

    fdir = '../data/' + plot_config["mhd_run"] + '/exhaust_boundary/'
    fname = fdir + 'xz_right_' + str(tframe) + '.dat'
    xlist_right, ylist_right = np.fromfile(fname).reshape([2, -1])
    fname = fdir + 'xz_left_' + str(tframe) + '.dat'
    xlist_left, ylist_left = np.fromfile(fname).reshape([2, -1])
    if not np.all(xlist_right > xlist_left):
        iy_max = np.argmax(ylist_right)
        xlist_right = np.copy(xlist_right[:iy_max+1])
        ylist_right = np.copy(ylist_right[:iy_max+1])
        xlist_left = lx_mhd - xlist_right
        ylist_left = np.copy(ylist_right)

    fpath = plot_config["mhd_run_dir"] + 'bin_data/'
    fname = fpath + 'deltab_' + str(tframe).zfill(4)
    deltab = np.fromfile(fname, dtype=np.float32).reshape([ny, nx])

    fig = plt.figure(figsize=[5.5, 8])
    rect0 = [0.15, 0.1, 0.8, 0.85]
    rect = np.copy(rect0)
    hgap = 0.03
    rect = np.copy(rect0)
    ax = fig.add_axes(rect)
    sizes = [xmin_grid, xmax_grid, ymin_grid, ymax_grid]
    sizes = np.asarray(sizes) * L0

    ax.plot(xlist_left*L0, ylist_left*L0, color='w',
            linewidth=1, linestyle='--')
    ax.plot(xlist_right*L0, ylist_right*L0, color='w',
            linewidth=1, linestyle='--')
    img = ax.imshow(deltab, extent=sizes, cmap=plt.cm.plasma,
                    norm=LogNorm(vmin=1E-3, vmax=3E-2),
                    aspect='auto', origin='lower', interpolation='bicubic')

    rect_cbar = np.copy(rect)
    rect_cbar[0] += 0.02
    rect_cbar[2] = 0.02
    rect_cbar[1] = 0.4
    rect_cbar[3] = 0.3
    cbar_ax = fig.add_axes(rect_cbar)
    cbar = fig.colorbar(img, cax=cbar_ax, extend="both")
    cbar.ax.tick_params(labelsize=10, color='w')
    cbar.ax.yaxis.set_tick_params(color='w')
    cbar.outline.set_edgecolor('w')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')
    fig_text = r'$\delta B^2/B_0^2$'
    cbar_ax.set_ylabel(fig_text, fontsize=16, color='w')
    # ax.text(0.05, 0.95, fig_text, color='white', fontsize=12,
    #         horizontalalignment='left', verticalalignment='center',
    #         transform=ax.transAxes)

    mhd_time = mhd_config["dt_out"][0] * tframe
    tva = mhd_time * L0
    title = r'$t = ' + "{:10.1f}".format(tva) + r'\text{ s}$'
    plt.suptitle(title, fontsize=20)

    ax.set_xlim(sizes[:2])
    ax.set_ylim(sizes[2:])
    ax.tick_params(labelsize=12)
    ax.set_xlabel(r'$x$/Mm', fontsize=16)
    ax.set_ylabel(r'$y$ (Mm)', fontsize=16)

    fdir = '../img/deltab/' + plot_config["mhd_run"] + '/'
    mkdir_p(fdir)
    fname = fdir + 'deltab_' + str(tframe) + '.jpg'
    fig.savefig(fname, dpi=200)

    if show_plot:
        plt.show()
    else:
        plt.close()

    # plt.plot(deltab[:, nx_mhd//2])
    # plt.show()


def plot_correlation_length(plot_config, mhd_config, show_plot=True):
    """Plot the turbulence correlation length

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    L0 = 200.0  # in Mm

    tframe = plot_config["tframe"]
    nx = mhd_config["nx"][0] + 4
    ny = mhd_config["ny"][0] + 4
    dx = mhd_config["dx"][0]
    dy = mhd_config["dy"][0]
    xmin, = mhd_config["xmin"]
    xmax, = mhd_config["xmax"]
    ymin, = mhd_config["ymin"]
    ymax, = mhd_config["ymax"]
    nx_mhd, = mhd_config["nx"]
    ny_mhd, = mhd_config["ny"]
    lx_mhd = xmax - xmin
    ly_mhd = ymax - ymin
    dx = lx_mhd / nx_mhd
    dy = ly_mhd / nx_mhd
    xmin_grid = xmin - 2 * dx
    xmax_grid = xmax + 2 * dx
    ymin_grid = ymin - 2 * dy
    ymax_grid = ymax + 2 * dy
    lx_grid = xmax_grid - xmin_grid
    ly_grid = ymax_grid - ymin_grid
    xmhd = np.linspace(0.5*dx, lx_mhd-0.5*dx, nx_mhd)
    ymhd = np.linspace(0.5*dy, ly_mhd-0.5*dy, ny_mhd)
    xgrid = np.linspace(-1.5*dx, lx_mhd+1.5*dx, nx)
    ygrid = np.linspace(-1.5*dy, ly_mhd+1.5*dy, ny)

    fdir = '../data/' + plot_config["mhd_run"] + '/exhaust_boundary/'
    fname = fdir + 'xz_right_' + str(tframe) + '.dat'
    xlist_right, ylist_right = np.fromfile(fname).reshape([2, -1])
    fname = fdir + 'xz_left_' + str(tframe) + '.dat'
    xlist_left, ylist_left = np.fromfile(fname).reshape([2, -1])
    if not np.all(xlist_right > xlist_left):
        iy_max = np.argmax(ylist_right)
        xlist_right = np.copy(xlist_right[:iy_max+1])
        ylist_right = np.copy(ylist_right[:iy_max+1])
        xlist_left = lx_mhd - xlist_right
        ylist_left = np.copy(ylist_right)

    fpath = plot_config["mhd_run_dir"] + 'bin_data/'
    fname = fpath + 'lc_' + str(tframe).zfill(4)
    lc = np.fromfile(fname, dtype=np.float32).reshape([ny, nx])

    fig = plt.figure(figsize=[5.5, 8])
    rect0 = [0.15, 0.1, 0.8, 0.85]
    rect = np.copy(rect0)
    hgap = 0.03
    rect = np.copy(rect0)
    ax = fig.add_axes(rect)
    sizes = [xmin_grid, xmax_grid, ymin_grid, ymax_grid]
    sizes = np.asarray(sizes) * L0

    ax.plot(xlist_left*L0, ylist_left*L0, color='w',
            linewidth=1, linestyle='-')
    ax.plot(xlist_right*L0, ylist_right*L0, color='w',
            linewidth=1, linestyle='-')
    img = ax.imshow(lc, extent=sizes, cmap=plt.cm.viridis,
                    norm=LogNorm(vmin=1E2, vmax=5E4),
                    aspect='auto', origin='lower', interpolation='bicubic')

    rect_cbar = np.copy(rect)
    rect_cbar[0] += 0.02
    rect_cbar[2] = 0.02
    rect_cbar[1] = 0.4
    rect_cbar[3] = 0.3
    cbar_ax = fig.add_axes(rect_cbar)
    cbar = fig.colorbar(img, cax=cbar_ax, extend="both")
    cbar.ax.tick_params(labelsize=10, color='w')
    cbar.ax.yaxis.set_tick_params(color='w')
    cbar.outline.set_edgecolor('w')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')
    fig_text = r'$L_c/\text{km}$'
    cbar_ax.set_ylabel(fig_text, fontsize=16, color='w')
    # ax.text(0.05, 0.95, fig_text, color='white', fontsize=12,
    #         horizontalalignment='left', verticalalignment='center',
    #         transform=ax.transAxes)

    mhd_time = mhd_config["dt_out"][0] * tframe
    tva = mhd_time * L0
    title = r'$t = ' + "{:10.1f}".format(tva) + r'\text{ s}$'
    plt.suptitle(title, fontsize=20)

    ax.set_xlim(sizes[:2])
    ax.set_ylim(sizes[2:])
    ax.tick_params(labelsize=12)
    ax.set_xlabel(r'$x$/Mm', fontsize=16)
    ax.set_ylabel(r'$y$ (Mm)', fontsize=16)

    fdir = '../img/lc/' + plot_config["mhd_run"] + '/'
    mkdir_p(fdir)
    fname = fdir + 'lc_' + str(tframe) + '.jpg'
    fig.savefig(fname, dpi=200)

    if show_plot:
        plt.show()
    else:
        plt.close()


def calc_shear_params(plot_config, mhd_config, show_plot=True):
    """calculate flow shear parameters

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    L0 = 200.0  # in Mm

    tframe = plot_config["tframe"]
    nx = mhd_config["nx"][0] + 4
    ny = mhd_config["ny"][0] + 4
    dx = mhd_config["dx"][0]
    dy = mhd_config["dy"][0]
    fpath = plot_config["mhd_run_dir"] + 'bin_data/'
    fname = fpath + 'mhd_data_' + str(tframe).zfill(4)
    mhd_fields = np.fromfile(fname, dtype=np.float32)
    mhd_fields = mhd_fields.reshape((ny, nx, 8))
    vx = mhd_fields[:, :, 0]
    vy = mhd_fields[:, :, 1]
    dvx_dx = np.gradient(vx, axis=1)
    dvx_dy = np.gradient(vx, axis=0)
    dvy_dx = np.gradient(vy, axis=1)
    dvy_dy = np.gradient(vy, axis=0)
    divv = dvx_dx + dvy_dy
    gshear = ((2*(dvx_dy+dvy_dx)**2 + 4*(dvx_dx**2+dvy_dy**2)) / 30 -
              (2*divv**2 / 45))
    print("Min and Max of shear parameter: %f, %f" % (gshear.min(), gshear.max()))
    print("Min and Max of compression: %f, %f" % (divv.min(), divv.max()))

    fig = plt.figure(figsize=[8, 8])
    rect0 = [0.1, 0.1, 0.4, 0.85]
    rect = np.copy(rect0)
    hgap = 0.03
    rect = np.copy(rect0)
    ax = fig.add_axes(rect)
    sizes = [mhd_config["xmin"][0], 0.25*mhd_config["xmax"][0],
             mhd_config["ymin"][0], mhd_config["ymax"][0]]
    sizes = np.asarray(sizes) * L0
    xs = 3*nx//8
    fdata = gshear[:, xs:(nx-xs)]
    xgrid = np.linspace(sizes[0], sizes[1], nx-xs*2)
    ygrid = np.linspace(sizes[2], sizes[3], ny)
    nx_grid, = xgrid.shape
    ny_grid, = ygrid.shape
    img = ax.imshow(fdata, extent=sizes, cmap=plt.cm.viridis,
                    # norm=LogNorm(vmin=1E2, vmax=1E3),
                    vmin=0, vmax=1E-3,
                    aspect='auto', origin='lower', interpolation='none')
    ax.tick_params(labelsize=12)
    ax.set_xlabel(r'$x$/Mm', fontsize=16)
    ax.set_ylabel(r'$y$ (Mm)', fontsize=16)

    rect_cbar = np.copy(rect)
    rect_cbar[0] += 0.02
    rect_cbar[2] = 0.015
    rect_cbar[1] = 0.4
    rect_cbar[3] = 0.3
    cbar_ax = fig.add_axes(rect_cbar)
    cbar = fig.colorbar(img, cax=cbar_ax, extend="both")
    cbar.ax.tick_params(labelsize=10, color='w')
    cbar.ax.yaxis.set_tick_params(color='w')
    cbar.outline.set_edgecolor('w')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')
    fig_text = r'$\Gamma$'
    ax.text(0.05, 0.95, fig_text, color='w', fontsize=16,
            horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    rect[0] += rect[2] + hgap
    ax = fig.add_axes(rect)
    fdata = divv[:, xs:(nx-xs)]
    img = ax.imshow(fdata, extent=sizes, cmap=plt.cm.seismic,
                    # norm=LogNorm(vmin=1E2, vmax=1E3),
                    vmin=-1E-2, vmax=1E-2,
                    aspect='auto', origin='lower', interpolation='none')
    ax.tick_params(labelsize=12)
    ax.set_xlabel(r'$x$/Mm', fontsize=16)
    ax.tick_params(axis='y', labelleft=False)

    rect_cbar = np.copy(rect)
    rect_cbar[0] += 0.02
    rect_cbar[2] = 0.015
    rect_cbar[1] = 0.4
    rect_cbar[3] = 0.3
    cbar_ax = fig.add_axes(rect_cbar)
    cbar = fig.colorbar(img, cax=cbar_ax, extend="both")
    cbar.ax.tick_params(labelsize=10, color='k')
    cbar.ax.yaxis.set_tick_params(color='k')
    cbar.outline.set_edgecolor('k')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='k')
    fig_text = r'$\Delta\cdot\boldsymbol{V}$'
    ax.text(0.05, 0.95, fig_text, color='k', fontsize=16,
            horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    fdir = '../img/shear_params/' + plot_config["mhd_run"] + '/'
    mkdir_p(fdir)
    # fname = fdir + 'shear_params_' + str(tframe) + '.jpg'
    # fig.savefig(fname, dpi=200)

    if show_plot:
        plt.show()
    else:
        plt.close()


def get_cmd_args():
    """Get command line arguments
    """
    default_mhd_run = 'run_fr++_20200310_noamr'
    default_mhd_run_dir = ('/net/scratch3/xiaocan/mhd/fluxropes/sep_10th_2017/' +
                           default_mhd_run + '/')
    parser = argparse.ArgumentParser(description='Spatial distribution')
    parser.add_argument('--mhd_run', action="store",
                        default=default_mhd_run, help='MHD run name')
    parser.add_argument('--mhd_run_dir', action="store",
                        default=default_mhd_run_dir, help='MHD run directory')
    parser.add_argument('--tframe', action="store", default='100', type=int,
                        help='Time frame')
    parser.add_argument('--multi_frames', action="store_true", default=False,
                        help='whether to analyze multiple frames')
    parser.add_argument('--time_loop', action="store_true", default=False,
                        help='whether to use a time loop to analyze multiple frames')
    parser.add_argument('--tstart', action="store", default='0', type=int,
                        help='starting time frame')
    parser.add_argument('--tend', action="store", default='10', type=int,
                        help='ending time frame')
    parser.add_argument('--inject_at_jz', action="store_true", default=False,
                        help='Inject particles where jz is large')
    parser.add_argument('--calc_deltab', action="store_true", default=False,
                        help='Calculate the magnetic fluctuation')
    parser.add_argument('--calc_lc', action="store_true", default=False,
                        help='Calculate the correlation length')
    parser.add_argument('--plot_deltab', action="store_true", default=False,
                        help='Plot the magnetic fluctuation')
    parser.add_argument('--plot_lc', action="store_true", default=False,
                        help='Plot the correlation length')
    parser.add_argument('--exhaust_boundary', action="store_true", default=False,
                        help='whether to get the boundary of reconnection exhaust')
    parser.add_argument('--shear_params', action="store_true", default=False,
                        help='whether to shear parameters')
    parser.add_argument('--lc_norm', action="store", default='5000.0', type=float,
                        help='Normalization for the turbulence correlation length (km)')
    parser.add_argument('--lc_min', action="store", default='500.0', type=float,
                        help='Minimum of the turbulence correlation length (km)')
    return parser.parse_args()


def analysis_single_frames(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    if args.inject_at_jz:
        inject_at_jz(plot_config, mhd_config, show_plot=True)
    elif args.calc_deltab:
        calc_deltab(plot_config, mhd_config)
    elif args.calc_lc:
        calc_correlation_length(plot_config, mhd_config)
    elif args.plot_deltab:
        plot_deltab(plot_config, mhd_config, show_plot=True)
    elif args.plot_lc:
        plot_correlation_length(plot_config, mhd_config, show_plot=True)
    elif args.exhaust_boundary:
        get_exhaust_boundary(plot_config, mhd_config, show_plot=True)
    if args.shear_params:
        calc_shear_params(plot_config, mhd_config, show_plot=True)


def process_input(plot_config, mhd_config, args, tframe):
    plot_config["tframe"] = tframe
    print("Time frame: %d" % tframe)
    if args.exhaust_boundary:
        get_exhaust_boundary(plot_config, mhd_config, show_plot=False)
    elif args.calc_deltab:
        calc_deltab(plot_config, mhd_config)
    elif args.calc_lc:
        calc_correlation_length(plot_config, mhd_config)


def analysis_multi_frames(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    tframes = range(plot_config["tstart"], plot_config["tend"] + 1)
    if args.time_loop:
        for tframe in tframes:
            print("Time frame: %d" % tframe)
            plot_config["tframe"] = tframe
            if args.inject_at_jz:
                inject_at_jz(plot_config, mhd_config, show_plot=False)
            elif args.plot_deltab:
                plot_deltab(plot_config, mhd_config, show_plot=False)
            elif args.calc_lc:
                calc_correlation_length(plot_config, mhd_config)
            elif args.plot_lc:
                plot_correlation_length(plot_config, mhd_config, show_plot=False)
    else:
        ncores = multiprocessing.cpu_count()
        # ncores = 18
        Parallel(n_jobs=ncores)(delayed(process_input)(plot_config, mhd_config,
                                                       args, tframe)
                                for tframe in tframes)


def main():
    """business logic for when running this module as the primary one!"""
    args = get_cmd_args()
    mhd_config = load_mhd_config(args.mhd_run_dir)
    plot_config = {}
    plot_config["tframe"] = args.tframe
    plot_config["mhd_run"] = args.mhd_run
    plot_config["mhd_run_dir"] = args.mhd_run_dir
    plot_config["tstart"] = args.tstart
    plot_config["tend"] = args.tend
    plot_config["lc_norm"] = args.lc_norm
    plot_config["lc_min"] = args.lc_min
    if args.multi_frames:
        analysis_multi_frames(plot_config, mhd_config, args)
    else:
        analysis_single_frames(plot_config, mhd_config, args)


if __name__ == "__main__":
    main()

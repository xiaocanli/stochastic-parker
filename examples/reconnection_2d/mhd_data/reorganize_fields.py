"""
Reorganize MHD fields for Parker's transport equation.
Basically, we add two ghost cells on each side of the domain.
"""
from __future__ import print_function

import argparse
import math
import multiprocessing

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from joblib import Parallel, delayed

import mhd_data
from util import mkdir_p

plt.style.use("seaborn-deep")
mpl.rc('text', usetex=True)
mpl.rcParams["text.latex.preamble"] = \
        (r"\usepackage{amsmath, bm}" +
         r"\DeclareMathAlphabet{\mathsfit}{\encodingdefault}{\sfdefault}{m}{sl}" +
         r"\SetMathAlphabet{\mathsfit}{bold}{\encodingdefault}{\sfdefault}{bx}{sl}" +
         r"\newcommand{\tensorsym}[1]{\bm{\mathsfit{#1}}}")

def save_mhd_fields_with_ghost(mhd_run_info, data_range, tframe, gfilter=False):
    """Save MHD fields with ghost cells

    The 8 fields include vx, vy, vz, v, bx, by, bz, b.
    rho and pressure are saved as separated files.
    """
    if mhd_run_info["mhd_code"] == "Athena++":
        xmesh, _, ftmp = mhd_data.read_fields_data(mhd_run_info, tframe)
        ny, nx = ftmp["rho"][0].shape
        athenapp_fields = ['rho', 'press', 'vel1', 'vel2', 'vel3', 'Bcc1', 'Bcc2', 'Bcc3']
        shape = [nx, ny, len(athenapp_fields)]
        fdata = np.zeros(shape)
        for ivar, var in enumerate(athenapp_fields):
            fdata[:, :, ivar] = ftmp[var][0].T
    else:
        xmesh, _, fdata = mhd_data.read_fields_data(mhd_run_info, tframe)
        ny, nx = xmesh.shape
    ixs = math.floor(data_range[0] * nx)
    ixe = math.ceil(data_range[1] * nx)
    iys = math.floor(data_range[2] * ny)
    iye = math.ceil(data_range[3] * ny)
    nxs = ixe - ixs
    nys = iye - iys
    nxsh = nxs
    if mhd_run_info["xmirror"]:
        nxs *= 2
    nysh = nys
    if mhd_run_info["ymirror"]:
        nys *= 2
    mhd_fields = np.zeros((nys+4, nxs+4, 8), dtype=np.float32)
    rho = np.zeros((nys+4, nxs+4), dtype=np.float32)
    pre = np.zeros((nys+4, nxs+4), dtype=np.float32)
    if mhd_run_info["with_z_component"]:
        # fdata: rho, p, vx, vy, vz, bx, by, bz
        absb = np.sqrt(np.sum(fdata[:, :, 5:8]**2, axis=2))
        if mhd_run_info["mhd_code"] == "Athena++" and mhd_run_info["output_type"] == "Fluxrope.out2":
            mhd_fields[2:nysh+2, nxsh+2:nxs+2, 0] = fdata[ixs:ixe, iys:iye, 2].T
            mhd_fields[2:nysh+2, nxsh+2:nxs+2, 1] = fdata[ixs:ixe, iys:iye, 3].T
            mhd_fields[2:nysh+2, nxsh+2:nxs+2, 2] = fdata[ixs:ixe, iys:iye, 4].T
            mhd_fields[2:nysh+2, nxsh+2:nxs+2, 3] = fdata[ixs:ixe, iys:iye, 0].T
            mhd_fields[2:nysh+2, nxsh+2:nxs+2, 4] = fdata[ixs:ixe, iys:iye, 5].T
            mhd_fields[2:nysh+2, nxsh+2:nxs+2, 5] = fdata[ixs:ixe, iys:iye, 6].T
            mhd_fields[2:nysh+2, nxsh+2:nxs+2, 6] = fdata[ixs:ixe, iys:iye, 7].T
            mhd_fields[2:nysh+2, nxsh+2:nxs+2, 7] = absb[ixs:ixe, iys:iye].T
            rho[2:nysh+2, nxsh+2:nxs+2] = fdata[ixs:ixe, iys:iye, 0].T
            pre[2:nysh+2, nxsh+2:nxs+2] = fdata[ixs:ixe, iys:iye, 1].T
        else:
            mhd_fields[2:nysh+2, 2:nxsh+2, 0] = fdata[ixs:ixe, iys:iye, 2].T
            mhd_fields[2:nysh+2, 2:nxsh+2, 1] = fdata[ixs:ixe, iys:iye, 3].T
            mhd_fields[2:nysh+2, 2:nxsh+2, 2] = fdata[ixs:ixe, iys:iye, 4].T
            mhd_fields[2:nysh+2, 2:nxsh+2, 3] = fdata[ixs:ixe, iys:iye, 0].T
            mhd_fields[2:nysh+2, 2:nxsh+2, 4] = fdata[ixs:ixe, iys:iye, 5].T
            mhd_fields[2:nysh+2, 2:nxsh+2, 5] = fdata[ixs:ixe, iys:iye, 6].T
            mhd_fields[2:nysh+2, 2:nxsh+2, 6] = fdata[ixs:ixe, iys:iye, 7].T
            mhd_fields[2:nysh+2, 2:nxsh+2, 7] = absb[ixs:ixe, iys:iye].T
            rho[2:nys+2, 2:nxsh+2] = fdata[ixs:ixe, iys:iye, 0].T
            pre[2:nys+2, 2:nxsh+2] = fdata[ixs:ixe, iys:iye, 1].T
    else:
        # fdata: rho, p, vx, vy, bx, by (and maybe Az)
        absb = np.sqrt(np.sum(fdata[:, :, 4:6]**2, axis=2))
        if mhd_run_info["mhd_code"] == "Athena++" and mhd_run_info["output_type"] == "Fluxrope.out2":
            mhd_fields[2:nysh+2, nxsh+2:nxs+2, 0] = fdata[ixs:ixe, iys:iye, 2].T
            mhd_fields[2:nysh+2, nxsh+2:nxs+2, 1] = fdata[ixs:ixe, iys:iye, 3].T
            mhd_fields[2:nysh+2, nxsh+2:nxs+2, 3] = fdata[ixs:ixe, iys:iye, 0].T
            mhd_fields[2:nysh+2, nxsh+2:nxs+2, 4] = fdata[ixs:ixe, iys:iye, 4].T
            mhd_fields[2:nysh+2, nxsh+2:nxs+2, 5] = fdata[ixs:ixe, iys:iye, 5].T
            mhd_fields[2:nysh+2, nxsh+2:nxs+2, 7] = absb[ixs:ixe, iys:iye].T
            rho[2:nysh+2, nxsh+2:nxs+2] = fdata[ixs:ixe, iys:iye, 0].T
            pre[2:nysh+2, nxsh+2:nxs+2] = fdata[ixs:ixe, iys:iye, 1].T
        else:
            mhd_fields[2:nysh+2, 2:nxsh+2, 0] = fdata[ixs:ixe, iys:iye, 2].T
            mhd_fields[2:nysh+2, 2:nxsh+2, 1] = fdata[ixs:ixe, iys:iye, 3].T
            mhd_fields[2:nysh+2, 2:nxsh+2, 3] = fdata[ixs:ixe, iys:iye, 0].T
            mhd_fields[2:nysh+2, 2:nxsh+2, 4] = fdata[ixs:ixe, iys:iye, 4].T
            mhd_fields[2:nysh+2, 2:nxsh+2, 5] = fdata[ixs:ixe, iys:iye, 5].T
            mhd_fields[2:nysh+2, 2:nxsh+2, 7] = absb[ixs:ixe, iys:iye].T
            rho[2:nysh+2, 2:nxsh+2] = fdata[ixs:ixe, iys:iye, 0].T
            pre[2:nysh+2, 2:nxsh+2] = fdata[ixs:ixe, iys:iye, 1].T

    if mhd_run_info["xmirror"]:
        if mhd_run_info["mhd_code"] == "Athena++" and mhd_run_info["output_type"] == "Fluxrope.out2":
            mhd_fields[:, 2:nxsh+2, 0] = -mhd_fields[:, nxs+1:nxsh+1:-1, 0]
            mhd_fields[:, 2:nxsh+2, 1:4] = mhd_fields[:, nxs+1:nxsh+1:-1, 1:4]
            mhd_fields[:, 2:nxsh+2, 4] = mhd_fields[:, nxs+1:nxsh+1:-1, 4]
            mhd_fields[:, 2:nxsh+2, 5] = -mhd_fields[:, nxs+1:nxsh+1:-1, 5]
            mhd_fields[:, 2:nxsh+2, 6:] = mhd_fields[:, nxs+1:nxsh+1:-1, 6:]
            rho[:, 2:nxsh+2] = rho[:, nxs+1:nxsh+1:-1]
            pre[:, 2:nxsh+2] = pre[:, nxs+1:nxsh+1:-1]
        else:
            mhd_fields[:, nxsh+2:nxs+2, 0] = -mhd_fields[:, nxsh+1:1:-1, 0]
            mhd_fields[:, nxsh+2:nxs+2, 1:4] = mhd_fields[:, nxsh+1:1:-1, 1:4]
            mhd_fields[:, nxsh+2:nxs+2, 4] = mhd_fields[:, nxsh+1:1:-1, 4]
            mhd_fields[:, nxsh+2:nxs+2, 5] = -mhd_fields[:, nxsh+1:1:-1, 5]
            mhd_fields[:, nxsh+2:nxs+2, 6:] = mhd_fields[:, nxsh+1:1:-1, 6:]
            rho[:, nxsh+2:nxs+2] = rho[:, nxsh+1:1:-1]
            pre[:, nxsh+2:nxs+2] = pre[:, nxsh+1:1:-1]
    if mhd_run_info["ymirror"]:
        mhd_fields[nysh+2:nys+2, :, 0] = mhd_fields[nysh+1:1:-1, :, 0]
        mhd_fields[nysh+2:nys+2, :, 1] = -mhd_fields[nysh+1:1:-1, :, 1]
        mhd_fields[nysh+2:nys+2, :, 2:4] = mhd_fields[nysh+1:1:-1, :, 2:4]
        mhd_fields[nysh+2:nys+2, :, 4] = -mhd_fields[nysh+1:1:-1, :, 4]
        mhd_fields[nysh+2:nys+2, :, 5] = mhd_fields[nysh+1:1:-1, :, 5]
        mhd_fields[nysh+2:nys+2, :, 6:] = mhd_fields[nysh+1:1:-1, :, 6:]
        rho[nysh+2:nys+2, :] = rho[nysh+1:1:-1, :]
        pre[nysh+2:nys+2, :] = pre[nysh+1:1:-1, :]

    if mhd_run_info["boundary"] == 0:
        mhd_fields[0:2, :, :] = mhd_fields[nys-1:nys+1, :, :]
        mhd_fields[:, 0:2, :] = mhd_fields[:, nxs-1:nxs+1, :]
        mhd_fields[nys+2:, :, :] = mhd_fields[3:5, :, :]
        mhd_fields[:, nxs+2:, :] = mhd_fields[:, 3:5, :]
        mhd_fields[0:2, 0:2, :] = mhd_fields[nys-1:nys+1, nxs-1:nxs+1, :]
        mhd_fields[nys+2:, nxs+2:, :] = mhd_fields[3:5, 3:5, :]
        mhd_fields[0:2, nxs+2:, :] = mhd_fields[nys-1:nys+1, 3:5, :]
        mhd_fields[nys+2:, 0:2, :] = mhd_fields[3:5, nxs-1:nxs+1, :]
        rho[0:2, :] = rho[nys-1:nys+1, :]
        rho[:, 0:2] = rho[:, nxs-1:nxs+1]
        rho[nys+2:, :] = rho[3:5, :]
        rho[:, nxs+2:] = rho[:, 3:5]
        rho[0:2, 0:2] = rho[nys-1:nys+1, nxs-1:nxs+1]
        rho[nys+2:, nxs+2:] = rho[3:5, 3:5]
        rho[0:2, nxs+2:] = rho[nys-1:nys+1, 3:5]
        rho[nys+2:, 0:2] = rho[3:5, nxs-1:nxs+1]
        pre[0:2, :] = pre[nys-1:nys+1, :]
        pre[:, 0:2] = pre[:, nxs-1:nxs+1]
        pre[nys+2:, :] = pre[3:5, :]
        pre[:, nxs+2:] = pre[:, 3:5]
        pre[0:2, 0:2] = pre[nys-1:nys+1, nxs-1:nxs+1]
        pre[nys+2:, nxs+2:] = pre[3:5, 3:5]
        pre[0:2, nxs+2:] = pre[nys-1:nys+1, 3:5]
        pre[nys+2:, 0:2] = pre[3:5, nxs-1:nxs+1]
    elif mhd_run_info["boundary"] == 1:
        mhd_fields[0:2, :, :] = mhd_fields[3:1:-1, :, :]
        mhd_fields[:, 0:2, :] = mhd_fields[:, 3:1:-1, :]
        mhd_fields[nys+2:, :, :] = mhd_fields[nys+1:nys-1:-1, :, :]
        mhd_fields[:, nxs+2:, :] = mhd_fields[:, nxs+1:nxs-1:-1, :]
        mhd_fields[0:2, 0:2, :] = mhd_fields[3:1:-1, 3:1:-1, :]
        mhd_fields[nys+2:, nxs+2:, :] = mhd_fields[nys+1:nys-1:-1, nxs+1:nxs-1:-1, :]
        mhd_fields[0:2, nxs+2:, :] = mhd_fields[3:1:-1, nxs+1:nxs-1:-1, :]
        mhd_fields[nys+2:, 0:2, :] = mhd_fields[nys+1:nys-1:-1, 3:1:-1, :]
        rho[0:2, :] = rho[3:1:-1, :]
        rho[:, 0:2] = rho[:, 3:1:-1]
        rho[nys+2:, :] = rho[nys+1:nys-1:-1, :]
        rho[:, nxs+2:] = rho[:, nxs+1:nxs-1:-1]
        rho[0:2, 0:2] = rho[3:1:-1, 3:1:-1]
        rho[nys+2:, nxs+2:] = rho[nys+1:nys-1:-1, nxs+1:nxs-1:-1]
        rho[0:2, nxs+2:] = rho[3:1:-1, nxs+1:nxs-1:-1]
        rho[nys+2:, 0:2] = rho[nys+1:nys-1:-1, 3:1:-1]
        pre[0:2, :] = pre[3:1:-1, :]
        pre[:, 0:2] = pre[:, 3:1:-1]
        pre[nys+2:, :] = pre[nys+1:nys-1:-1, :]
        pre[:, nxs+2:] = pre[:, nxs+1:nxs-1:-1]
        pre[0:2, 0:2] = pre[3:1:-1, 3:1:-1]
        pre[nys+2:, nxs+2:] = pre[nys+1:nys-1:-1, nxs+1:nxs-1:-1]
        pre[0:2, nxs+2:] = pre[3:1:-1, nxs+1:nxs-1:-1]
        pre[nys+2:, 0:2] = pre[nys+1:nys-1:-1, 3:1:-1]

    fpath = mhd_run_info["run_dir"] + 'bin_data/'
    mkdir_p(fpath)
    fname = fpath + 'mhd_data_' + str(tframe).zfill(4)
    print("Field shape: ", mhd_fields.shape)
    print("Whether it is column-major: ", np.isfortran(mhd_fields))
    if gfilter:
        for i in range(8):
            mhd_fields[:, :, i] = gaussian_filter(mhd_fields[:, :, i], sigma=3)
    mhd_fields.tofile(fname)

    fname = fpath + 'rho_' + str(tframe).zfill(4)
    rho.tofile(fname)
    fname = fpath + 'pre_' + str(tframe).zfill(4)
    pre.tofile(fname)

    if tframe == 0:
        save_mhd_config(mhd_run_info, nxs, nys)


def save_mhd_config(mhd_run_info, nxs, nys):
    """Save MHD configuration
    """
    mhd_config = mhd_data.read_mhd_config(mhd_run_info["config_name"],
                                          mhd_run_info["mhd_code"])
    double_data = np.zeros(13)
    int_data = np.zeros(14, dtype=np.int32)
    nx_mhd = nxs
    ny_mhd = nys
    nz_mhd = mhd_config.nz
    lx_mhd = (mhd_config.xmax - mhd_config.xmin) * nxs / mhd_config.nx
    ly_mhd = (mhd_config.ymax - mhd_config.ymin) * nys / mhd_config.ny
    lz_mhd = mhd_config.zmax - mhd_config.zmin
    mpi_sizex = 16
    mpi_sizey = 16
    mpi_sizez = 1
    if nx_mhd > 0:
        double_data[0] = lx_mhd / nx_mhd
    else:
        double_data[0] = 0.0
    if ny_mhd > 0:
        double_data[1] = ly_mhd / ny_mhd
    else:
        double_data[1] = 0.0
    if nz_mhd > 0:
        double_data[2] = lz_mhd / nz_mhd
    else:
        double_data[2] = 0.0
    double_data[3] = 0.0
    double_data[4] = 0.0
    double_data[5] = 0.0
    double_data[6] = lx_mhd
    double_data[7] = ly_mhd
    double_data[8] = lz_mhd
    double_data[9] = lx_mhd
    double_data[10] = ly_mhd
    double_data[11] = lz_mhd
    double_data[12] = mhd_config.dt_out
    int_data[0] = nx_mhd
    int_data[1] = ny_mhd
    int_data[2] = nz_mhd
    int_data[3] = nx_mhd // mpi_sizex
    int_data[4] = ny_mhd // mpi_sizey
    int_data[5] = nz_mhd // mpi_sizez
    int_data[6] = mpi_sizex
    int_data[7] = mpi_sizey
    int_data[8] = mpi_sizez
    int_data[9] = 9
    int_data[10] = 0 # Periodic boundary condition as default

    fpath = mhd_run_info["run_dir"] + 'bin_data/'
    mkdir_p(fpath)
    fname = fpath + 'mhd_config.dat'
    double_data.tofile(fname)
    with open(fname, 'a') as file_handler:
        int_data.tofile(file_handler)


def test_reorganized_mhd_data(mhd_run_info, mhd_config, data_range, tframe, show_plot=True):
    """Test reorganized MHD fields
    """
    fpath = mhd_run_info["run_dir"] + 'bin_data/'
    fname = fpath + 'mhd_data_' + str(tframe).zfill(4)
    ixs = math.floor(data_range[0] * mhd_config.nx)
    ixe = math.ceil(data_range[1] * mhd_config.nx)
    iys = math.floor(data_range[2] * mhd_config.ny)
    iye = math.ceil(data_range[3] * mhd_config.ny)
    nxs = ixe - ixs
    nys = iye - iys
    nxsh = nxs
    nysh = nys
    if mhd_run_info["xmirror"]:
        nxs *= 2
    if mhd_run_info["ymirror"]:
        nys *= 2
    nx = nxs + 4
    ny = nys + 4
    mhd_fields = np.fromfile(fname, dtype=np.float32)
    mhd_fields = mhd_fields.reshape((ny, nx, 8))
    # plt.plot(mhd_fields[0, :, 5])
    dx = (mhd_config.xmax - mhd_config.xmin) / mhd_config.nx
    dy = (mhd_config.ymax - mhd_config.ymin) / mhd_config.ny
    vx = mhd_fields[:, :, 0]
    vy = mhd_fields[:, :, 1]
    rho = mhd_fields[:, :, 3]
    bx = mhd_fields[:, :, 4]
    by = mhd_fields[:, :, 5]
    bz = mhd_fields[:, :, 6]
    absB = np.sqrt(bx**2 + by**2 + bz**2)
    print("Min and Max of magnetic field: %f, %f" % (absB.min(), absB.max()))
    va = absB / np.sqrt(rho)
    divv = np.gradient(vx, axis=1) / dx + np.gradient(vy, axis=0) / dy
    divv_max = 10
    jz = np.gradient(by, dx, axis=1) - np.gradient(bx, dy, axis=0)
    fig = plt.figure(figsize=[10, 8])
    rect = [0.10, 0.12, 0.75, 0.8]
    ax = fig.add_axes(rect)
    sizes = [mhd_config.xmin, mhd_config.xmax,
             mhd_config.ymin, mhd_config.ymax]
    # img = plt.imshow(divv, origin='lower', cmap=plt.cm.seismic,
    #                  vmin=-divv_max, vmax=divv_max)
    # img = plt.imshow(vx, origin='lower', cmap=plt.cm.seismic,
    #                  vmin=-0.1, vmax=0.1)
    # img = plt.imshow(by, origin='lower', cmap=plt.cm.seismic)
    # img = plt.imshow(absB, origin='lower', cmap=plt.cm.viridis)
    img = ax.imshow(jz, extent=sizes, cmap=plt.cm.seismic,
                    vmin=-500, vmax=500,
                    aspect='auto', origin='lower', interpolation='none')
    # img = ax.imshow(rho, extent=sizes, cmap=plt.cm.seismic,
    #                 vmin=0, vmax=4,
    #                 aspect='auto', origin='lower', interpolation='none')
    # img = ax.imshow(va, extent=sizes, cmap=plt.cm.seismic,
    #                 vmin=0, vmax=2,
    #                 aspect='auto', origin='lower', interpolation='none')
    rect[0] += rect[2] + 0.02
    rect[2] = 0.03
    cbar_ax = fig.add_axes(rect)
    cbar = fig.colorbar(img, cax=cbar_ax)
    ax.tick_params(labelsize=16)
    cbar.ax.tick_params(labelsize=16)
    ax.set_xlabel(r'$x$', fontsize=20)
    ax.set_ylabel(r'$y$', fontsize=20)
    cbar.ax.set_ylabel(r'$j_z$', fontsize=20)
    tva = mhd_config.dt_out * tframe / (mhd_config.ymax - mhd_config.ymin)
    title = r'$t = ' + "{:10.2f}".format(tva) + r'\tau_A$'
    title += ' (frame: %d)' % tframe
    fig.suptitle(title, fontsize=24)
    fpath = '../img/jz/' + mhd_run_info["run_name"] + '/'
    mkdir_p(fpath)
    plt.savefig(fpath + 'jz_' + str(tframe) + '.jpg', dpi=200)
    if show_plot:
        plt.show()
    else:
        plt.close()


def mhd_fields_tri(mhd_run_info, mhd_config, data_range, tframe, show_plot=True):
    """Plot three frames of MHD fields
    """
    fpath = mhd_run_info["run_dir"] + 'bin_data/'
    ixs = math.floor(data_range[0] * mhd_config.nx)
    ixe = math.ceil(data_range[1] * mhd_config.nx)
    iys = math.floor(data_range[2] * mhd_config.ny)
    iye = math.ceil(data_range[3] * mhd_config.ny)
    L0 = 200  # in Mm
    sizes = [-mhd_config.xmax*data_range[1], mhd_config.xmax*data_range[1],
             mhd_config.ymin*data_range[2], mhd_config.ymax*data_range[3]]
    dx = (mhd_config.xmax - mhd_config.xmin) / mhd_config.nx
    dy = (mhd_config.ymax - mhd_config.ymin) / mhd_config.ny
    sizes = np.asarray(sizes) * L0
    nxs = ixe - ixs
    nys = iye - iys
    nxsh = nxs
    nysh = nys
    if mhd_run_info["xmirror"]:
        nxs *= 2
    if mhd_run_info["ymirror"]:
        nys *= 2
    nx = nxs + 4
    ny = nys + 4
    fig = plt.figure(figsize=[16, 8])
    rect = [0.05, 0.12, 0.2, 0.8]
    tframes = [120, 146, 189, 200]
    for iframe, tframe in enumerate(tframes):
        ax = fig.add_axes(rect)
        fname = fpath + 'mhd_data_' + str(tframe).zfill(4)
        mhd_fields = np.fromfile(fname, dtype=np.float32)
        mhd_fields = mhd_fields.reshape((ny, nx, 8))
        bx = mhd_fields[:, :, 4]
        by = mhd_fields[:, :, 5]
        jz = np.gradient(by, dx, axis=1) - np.gradient(bx, dy, axis=0)
        img = ax.imshow(jz, extent=sizes, cmap=plt.cm.seismic,
                        vmin=-500, vmax=500,
                        aspect='auto', origin='lower', interpolation='none')
        ax.tick_params(labelsize=16)
        ax.set_xlabel(r'$x$/Mm', fontsize=20)
        if iframe == 0:
            ax.set_ylabel(r'$y$/Mm', fontsize=20)
        else:
            ax.tick_params(axis='y', labelleft=False)
        tva = mhd_config.dt_out * tframe * L0
        title = r'$t = ' + "{:10.1f}".format(tva) + r's$'
        ax.set_title(title, fontsize=24)
        rect[0] += rect[2] + 0.02
    rect[0] -= rect[2] + 0.02
    rect[0] += rect[2] + 0.02
    rect[2] = 0.01
    cbar_ax = fig.add_axes(rect)
    cbar = fig.colorbar(img, cax=cbar_ax)
    cbar.ax.tick_params(labelsize=16)
    cbar.ax.set_ylabel(r'$j_z$', fontsize=20)
    fpath = '../img/jz/' + mhd_run_info["run_name"] + '/'
    mkdir_p(fpath)
    plt.savefig(fpath + 'jz_tri.jpg', dpi=200)
    if show_plot:
        plt.show()
    else:
        plt.close()


def plot_vy_temp(mhd_run_info, mhd_config, data_range, tframe, show_plot=True):
    """Plot one frame of Vy and temperature
    """
    fpath = mhd_run_info["run_dir"] + 'bin_data/'
    ixs = math.floor(data_range[0] * mhd_config.nx)
    ixe = math.ceil(data_range[1] * mhd_config.nx)
    iys = math.floor(data_range[2] * mhd_config.ny)
    iye = math.ceil(data_range[3] * mhd_config.ny)
    sizes = [-mhd_config.xmax*data_range[1], mhd_config.xmax*data_range[1],
             mhd_config.ymin*data_range[2], mhd_config.ymax*data_range[3]]
    dx = (mhd_config.xmax - mhd_config.xmin) / mhd_config.nx
    dy = (mhd_config.ymax - mhd_config.ymin) / mhd_config.ny
    nxs = ixe - ixs
    nys = iye - iys
    nxsh = nxs
    nysh = nys
    if mhd_run_info["xmirror"]:
        nxs *= 2
    if mhd_run_info["ymirror"]:
        nys *= 2
    nx = nxs + 4
    ny = nys + 4
    fig = plt.figure(figsize=[4, 4])
    rect0 = [0.13, 0.12, 0.4, 0.85]
    hgap = 0.03
    rect = np.copy(rect0)
    ax = fig.add_axes(rect)
    fname = fpath + 'mhd_data_' + str(tframe).zfill(4)
    mhd_fields = np.fromfile(fname, dtype=np.float32)
    mhd_fields = mhd_fields.reshape((ny, nx, 8))
    vy = mhd_fields[:, :, 1]
    img = ax.imshow(vy, extent=sizes, cmap=plt.cm.seismic,
                    vmin=-1, vmax=1,
                    aspect='auto', origin='lower', interpolation='bicubic')
    ax.tick_params(bottom=True, top=False, left=True, right=False)
    ax.tick_params(axis='x', which='minor', direction='in')
    ax.tick_params(axis='x', which='major', direction='in')
    ax.tick_params(axis='y', which='minor', direction='in')
    ax.tick_params(axis='y', which='major', direction='in')
    ax.set_xlim([-0.3, 0.3])
    ax.tick_params(labelsize=10)
    ax.set_xlabel(r'$x$', fontsize=12)
    ax.set_ylabel(r'$y$', fontsize=12)
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
    fig_text = r'$V_y/V_A$'
    # cbar.ax.set_ylabel(fig_text, fontsize=12, color='w')
    ax.text(0.05, 0.95, fig_text, color='k', fontsize=12,
            horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    fname = fpath + 'pre_' + str(tframe).zfill(4)
    pre = np.fromfile(fname, dtype=np.float32).reshape([ny, nx])
    fname = fpath + 'rho_' + str(tframe).zfill(4)
    rho = np.fromfile(fname, dtype=np.float32).reshape([ny, nx])
    temp0 = 0.05
    temp = pre / rho / temp0
    rect = np.copy(rect0)
    rect[0] += rect[2] + hgap
    ax = fig.add_axes(rect)
    img = ax.imshow(temp, extent=sizes, cmap=plt.cm.plasma,
                    vmin=0, vmax=10,
                    aspect='auto', origin='lower', interpolation='bicubic')
    ax.tick_params(bottom=True, top=False, left=True, right=False)
    ax.tick_params(axis='x', which='minor', direction='in')
    ax.tick_params(axis='x', which='major', direction='in')
    ax.tick_params(axis='y', which='minor', direction='in')
    ax.tick_params(axis='y', which='major', direction='in')
    ax.set_xlim([-0.3, 0.3])
    ax.tick_params(labelsize=10)
    ax.set_xlabel(r'$x$', fontsize=12)
    ax.tick_params(axis='y', labelleft=False)
    rect_cbar = np.copy(rect)
    rect_cbar[0] += 0.02
    rect_cbar[2] = 0.015
    rect_cbar[1] = 0.4
    rect_cbar[3] = 0.3
    cbar_ax = fig.add_axes(rect_cbar)
    cbar = fig.colorbar(img, cax=cbar_ax, extend="max")
    cbar.ax.tick_params(labelsize=10, color='w')
    cbar.ax.yaxis.set_tick_params(color='w')
    cbar.outline.set_edgecolor('w')
    plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='w')
    fig_text = r'$T/T_0$'
    # cbar.ax.set_ylabel(fig_text, fontsize=12, color='w')
    ax.text(0.05, 0.95, fig_text, color='w', fontsize=12,
            horizontalalignment='left', verticalalignment='center',
            transform=ax.transAxes)

    fpath = '../img/vy_temp/' + mhd_run_info["run_name"] + '/'
    mkdir_p(fpath)
    fname = fpath + 'vy_temp_' + str(tframe) + '.pdf'
    plt.savefig(fname, dpi=200)
    if show_plot:
        plt.show()
    else:
        plt.close()


def get_cmd_args():
    """Get command line arguments """
    default_run_name = 'S1E5_beta01_bg00'
    default_run_dir = ('/net/scratch3/xiaocanli/mhd/guide_field_scaling/' +
                       default_run_name + '/')
    parser = argparse.ArgumentParser(description='Reorganize MHD fields')
    parser.add_argument('--run_dir', action="store", default=default_run_dir,
                        help='MHD run directory')
    parser.add_argument('--run_name', action="store", default=default_run_name,
                        help='MHD run name')
    parser.add_argument('--output_type', action="store", default="reconnection",
                        help='MHD output type')
    parser.add_argument('--mhd_code', action="store", default="Athena",
                        help='MHD code')
    parser.add_argument('--with_z_component', action="store_true", default=False,
                        help='whether vz and bz are in original MHD fields')
    parser.add_argument('--config_name', action="store",
                        default="athinput.reconnection",
                        help='MHD configuration filename')
    parser.add_argument('--boundary', action="store", default='0', type=int,
                        help='Boundary condition: 0 for periodic, 1 for open')
    parser.add_argument('--tframe', action="store", default='30', type=int,
                        help='Time frame')
    parser.add_argument('--tstart', action="store", default='0', type=int,
                        help='Starting time frame')
    parser.add_argument('--tend', action="store", default='10', type=int,
                        help='Ending time frame')
    parser.add_argument('--multi_frames', action="store_true", default=False,
                        help='whether to analyze multiple frames')
    parser.add_argument('--check_data', action="store_true", default=False,
                        help='whether to check reorganized data')
    parser.add_argument('--rx_min', action="store", default='0.0', type=float,
                        help='Minimum of x range ([0, 1))')
    parser.add_argument('--rx_max', action="store", default='1.0', type=float,
                        help='Maximum of x range ((0, 1])')
    parser.add_argument('--ry_min', action="store", default='0.0', type=float,
                        help='Minimum of y range ([0, 1))')
    parser.add_argument('--ry_max', action="store", default='1.0', type=float,
                        help='Maximum of y range ((0, 1])')
    parser.add_argument('--xmirror', action="store_true", default=False,
                        help='whether only simulating half domain along x')
    parser.add_argument('--ymirror', action="store_true", default=False,
                        help='whether only simulating half domain along y')
    parser.add_argument('--organize_mhd_field', action="store_true", default=False,
                        help='whether to organize MHD fields')
    parser.add_argument('--save_mhd_config', action="store_true", default=False,
                        help='whether to save original MHD configuration')
    parser.add_argument('--mhd_fields_tri', action="store_true", default=False,
                        help='whether to plot 3 frames of MHD fields')
    parser.add_argument('--vy_temp', action="store_true", default=False,
                        help='whether to plot Vy and temperature')
    parser.add_argument('--gaussian_filter', action="store_true", default=False,
                        help='whether to apply Gaussian filter to the data')
    return parser.parse_args()


def process_input(mhd_run_info, mhd_config, data_range, tframe, args):
    """process one time frame"""
    if args.organize_mhd_field:
        save_mhd_fields_with_ghost(mhd_run_info, data_range, tframe,
                                   args.gaussian_filter)
    elif args.check_data:
        test_reorganized_mhd_data(mhd_run_info, mhd_config, data_range,
                                  tframe, show_plot=False)


def main():
    """business logic for when running this module as the primary one!"""
    args = get_cmd_args()
    mhd_run_info = {}
    mhd_run_info["run_name"] = args.run_name
    mhd_run_info["run_dir"] = args.run_dir
    mhd_run_info["output_type"] = args.output_type
    mhd_run_info["mhd_code"] = args.mhd_code
    mhd_run_info["config_name"] = mhd_run_info["run_dir"] + args.config_name
    mhd_run_info["boundary"] = args.boundary
    mhd_run_info["xmirror"] = args.xmirror
    mhd_run_info["ymirror"] = args.ymirror
    mhd_run_info["with_z_component"] = args.with_z_component
    mhd_config = mhd_data.read_mhd_config(mhd_run_info["config_name"],
                                          mhd_run_info["mhd_code"])
    if args.save_mhd_config:
        mhd_data.save_mhd_config(mhd_config, mhd_run_info["run_dir"])
    data_range = [args.rx_min, args.rx_max, args.ry_min, args.ry_max]
    if not args.multi_frames:
        if args.organize_mhd_field:
            save_mhd_fields_with_ghost(mhd_run_info, data_range, args.tframe)
        elif args.check_data:
            test_reorganized_mhd_data(mhd_run_info, mhd_config,
                                      data_range, args.tframe)
        elif args.mhd_fields_tri:
            mhd_fields_tri(mhd_run_info, mhd_config, data_range,
                           args.tframe, show_plot=True)
        elif args.vy_temp:
            plot_vy_temp(mhd_run_info, mhd_config, data_range,
                         args.tframe, show_plot=True)
    else:
        tframes = range(args.tstart, args.tend + 1)
        ncores = multiprocessing.cpu_count()
        ncores = 8
        Parallel(n_jobs=ncores)(delayed(process_input)(mhd_run_info,
                                                       mhd_config,
                                                       data_range,
                                                       tframe, args)
                                for tframe in tframes)


if __name__ == "__main__":
    main()

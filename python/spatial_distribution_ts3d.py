#!/usr/bin/env python3
"""
Spatial distribution of energetic particles for 3D TS simulations
"""
from __future__ import print_function

import argparse
import itertools
import json
import math
import multiprocessing

from evtk.hl import gridToVTK
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed
from matplotlib.colors import LogNorm, SymLogNorm
from scipy.interpolate import interp2d

from sde_util import load_mhd_config, mkdir_p

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]

def convert_dist_vtk(plot_config, mhd_config, show_plot=True):
    """convert distributions to VTK format

    Args:
        plot_config: plot configuration in dictionary
        mhd_config: MHD simulation configuration
    """
    tframe = plot_config["tframe"]
    nreduce = plot_config["nreduce"]
    nx_mhd, = mhd_config["nx"]  # r
    ny_mhd, = mhd_config["ny"]  # theta
    nz_mhd, = mhd_config["nz"]  # phi
    nxr = nx_mhd // nreduce
    nyr = ny_mhd // nreduce
    nzr = nz_mhd // nreduce

    mhd_run_dir = plot_config["mhd_run_dir"]
    fname = mhd_run_dir + "bin_data/xpos.dat"
    xpos = np.fromfile(fname)
    fname = mhd_run_dir + "bin_data/ypos.dat"
    ypos = np.fromfile(fname)
    fname = mhd_run_dir + "bin_data/zpos.dat"
    zpos = np.fromfile(fname)

    xpos_r = xpos[::nreduce]
    ypos_r = ypos[::nreduce]
    zpos_r = zpos[::nreduce]

    run_name = plot_config["run_name"]
    tframe_str = str(tframe).zfill(4)
    fname = '../data/' + run_name + '/fxy-' + tframe_str + '_sum.dat'
    fdata = np.fromfile(fname)
    nbands = int(fdata[0])
    parray_bands = fdata[1:nbands+2]
    fdata = fdata[nbands+2:].reshape([-1, nzr, nyr, nxr])
    z3d, y3d, x3d = np.meshgrid(zpos_r, ypos_r, xpos_r, indexing='ij')
    nbands, _, _, _ = fdata.shape

    for iband in range(nbands):
        print("Band %d" % iband)
        fname = '../data/' + run_name + '/dist_' + tframe_str + '_' + str(iband)
        gridToVTK(fname, z3d, y3d, x3d, pointData = {"density" : fdata[iband]})


def get_cmd_args():
    """get command line arguments
    """
    default_mhd_run = 'ts_3d_test'
    default_mhd_run_dir = ('/net/scratch3/xiaocan/mhd/' + default_mhd_run + '/')
    default_sde_run = 'p000_b000_0004_001'
    parser = argparse.ArgumentParser(description='spatial distribution')
    parser.add_argument('--mhd_run', action="store",
                        default=default_mhd_run, help='mhd run name')
    parser.add_argument('--mhd_run_dir', action="store",
                        default=default_mhd_run_dir, help='mhd run directory')
    parser.add_argument('--sde_run', action="store",
                        default=default_sde_run, help='sde run name')
    parser.add_argument('--tframe', action="store", default='100', type=int,
                        help='time frame')
    parser.add_argument('--ene_band', action="store", default='2', type=int,
                        help='energy band')
    parser.add_argument('--multi_frames', action="store_true", default=False,
                        help='whether to analyze multiple frames')
    parser.add_argument('--time_loop', action="store_true", default=False,
                        help='whether to use a time loop to analyze multiple frames')
    parser.add_argument('--multi_bands', action="store_true", default=False,
                        help='whether to plot multiple energy bands')
    parser.add_argument('--tstart', action="store", default='0', type=int,
                        help='starting time frame')
    parser.add_argument('--tend', action="store", default='10', type=int,
                        help='ending time frame')
    parser.add_argument('--dist_vtk', action="store_true", default=False,
                        help='whether to convert distribution data to VTK')
    return parser.parse_args()


def analysis_single_frames(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    if args.dist_vtk:
        convert_dist_vtk(plot_config, mhd_config, show_plot=True)


def process_input(plot_config, mhd_config, args, tframe):
    """process one time frame"""
    plot_config["tframe"] = tframe
    if args.dist_vtk:
        convert_dist_vtk(plot_config, mhd_config, show_plot=False)


def analysis_multi_frames(plot_config, mhd_config, args):
    """Analysis for multiple time frames
    """
    tframes = range(plot_config["tmin"], plot_config["tmax"] + 1)
    if args.time_loop:
        for tframe in tframes:
            print("Time frame: %d" % tframe)
            plot_config["tframe"] = tframe
            if args.dist_vtk:
                convert_dist_vtk(plot_config, mhd_config, show_plot=False)
    else:
        ncores = multiprocessing.cpu_count()
        ncores = 16
        Parallel(n_jobs=ncores)(delayed(process_input)(plot_config, mhd_config,
                                                       args, tframe)
                                for tframe in tframes)


def main():
    """business logic for when running this module as the primary one!"""
    args = get_cmd_args()
    with open('config/spectrum_config_10ta.json', 'r') as file_handler:
        config = json.load(file_handler)
    mhd_config = load_mhd_config(args.mhd_run_dir)
    plot_config = {}
    run_name = args.mhd_run + "/" + args.sde_run
    sde_run_config = config[run_name]
    nreduce = sde_run_config["nreduce"]
    plot_config["nx"] = mhd_config["nx"] // nreduce
    plot_config["ny"] = mhd_config["ny"] // nreduce
    plot_config["nz"] = mhd_config["nz"] // nreduce
    plot_config["nreduce"] = nreduce
    plot_config["tmax"] = sde_run_config["tmax"]
    plot_config["tmin"] = sde_run_config["tmin"]
    if "e0" in sde_run_config:
        plot_config["e0"] = sde_run_config["e0"]
    else:
        plot_config["e0"] = 1.0
    if "high_bands" in sde_run_config:
        plot_config["high_bands"] = sde_run_config["high_bands"]
    else:
        plot_config["high_bands"] = [2, 5]
    plot_config["eband"] = args.ene_band
    plot_config["run_name"] = run_name
    plot_config["tframe"] = args.tframe
    plot_config["mhd_run_dir"] = args.mhd_run_dir
    if args.multi_frames:
        analysis_multi_frames(plot_config, mhd_config, args)
    else:
        analysis_single_frames(plot_config, mhd_config, args)


if __name__ == "__main__":
    main()

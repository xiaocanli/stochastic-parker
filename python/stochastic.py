#!/usr/bin/env python3
"""
Analysis procedures for stochastic integration
"""
from __future__ import print_function

import argparse
import json
import multiprocessing

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from joblib import Parallel, delayed

from momentum_energy_distributions import momentum_energy_distributions
from sde_util import mkdir_p

mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)
mpl.rcParams['text.latex.preamble'] = [r"\usepackage{amsmath}"]


def get_cmd_args():
    """Get command line arguments
    """
    default_run_name = 'shock_interaction/p000_b000_001'
    default_run_dir = '/net/scratch3/xiaocanli/mhd/shock/'
    parser = argparse.ArgumentParser(description='Stochastic integration of Parker equation')
    parser.add_argument('--run_dir', action="store", default=default_run_dir,
                        help='run directory')
    parser.add_argument('--run_name', action="store", default=default_run_name,
                        help='run name')
    parser.add_argument('--multi_frames', action="store_true", default=False,
                        help='whether analyzing multiple frames')
    parser.add_argument('--ene_bin', action="store", default='1', type=int,
                        help='Energy bin')
    parser.add_argument('--edist', action="store_true", default=False,
                        help='whether to plot energy distribution')
    parser.add_argument('--sdist', action="store_true", default=False,
                        help='whether to plot spatial distribution')
    parser.add_argument('--multi', action="store_true", default=False,
                        help='whether to plot multiple frames')
    parser.add_argument('--power_test', action="store_true", default=False,
                        help='whether this is a test for power-law fitting')
    return parser.parse_args()


def main():
    """business logic for when running this module as the primary one!"""
    args = get_cmd_args()
    run_name = args.run_name
    run_dir = args.run_dir
    # def processInput(job_id):
    #     print(job_id)
    #     ct = job_id
    #     if args.sdist:
    #         plot_spatial_distributions_high_energy(ct, nx, ny, run_name, ene_bin=ebin)
    #         plt.close()
    # ncores = multiprocessing.cpu_count()
    # if args.multi:
    #     Parallel(n_jobs=ncores)(delayed(processInput)(ct) for ct in cts)
    # # spatial_distributions(ct, nx, ny, run_name)
    # if args.edist:
    #     momentum_energy_distributions(plot_config, args.power_test)
    # if args.sdist and (not args.multi):
    #     plot_spatial_distributions_high_energy(tmax, nx, ny, run_name, ene_bin=ebin)
    #     plt.show()


if __name__ == "__main__":
    main()

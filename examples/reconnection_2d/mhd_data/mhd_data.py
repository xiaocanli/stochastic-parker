#!/bin/env python3
"""
Module for MHD data and configuration
"""
from __future__ import print_function

import numpy as np

import athena_data
import athenapp_data
from util import mkdir_p

def read_fields_data(mhd_run_info, tframe):
    """Read MHD data
    """
    run_dir = mhd_run_info["run_dir"]
    mhd_code = mhd_run_info["mhd_code"]
    if mhd_code == "Athena":
        output_type = mhd_run_info["output_type"]
        xmesh, ymesh, fdata = athena_data.read_fields_data_athena(run_dir,
                                                                  tframe,
                                                                  output_type)
        return (xmesh, ymesh, fdata)
    elif mhd_code == "Athena++":
        output_type = mhd_run_info["output_type"]
        xmesh, ymesh, fdata = athenapp_data.read_fields_data_athenapp(run_dir,
                                                                      tframe,
                                                                      output_type)
        return (xmesh, ymesh, fdata)
    else:
        print("MHD code is not supported for now")


def read_mhd_config(config_filename, mhd_code):
    """Read MHD simulation configuration
    """
    if mhd_code == "Athena":
        return athena_data.read_mhd_config_athena(config_filename)
    elif mhd_code == "Athena++":
        return athenapp_data.read_mhd_config_athena(config_filename)
    else:
        print("MHD code is not supported for now")


def save_mhd_config(mhd_config, run_directory):
    """Save MHD configuration

    Note that the configuration structure is different from mhd_configuration
    above. It is used in Fortran code.
    """
    double_data = np.zeros(13)
    int_data = np.zeros(14, dtype=np.int32)
    lx_mhd = mhd_config.xmax - mhd_config.xmin
    ly_mhd = mhd_config.ymax - mhd_config.ymin
    lz_mhd = mhd_config.zmax - mhd_config.zmin
    nx_mhd = mhd_config.nx
    ny_mhd = mhd_config.ny
    nz_mhd = mhd_config.nz
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
    double_data[3] = mhd_config.xmin
    double_data[4] = mhd_config.ymin
    double_data[5] = mhd_config.zmin
    double_data[6] = mhd_config.xmax
    double_data[7] = mhd_config.ymax
    double_data[8] = mhd_config.zmax
    double_data[9] = lx_mhd
    double_data[10] = ly_mhd
    double_data[11] = lz_mhd
    double_data[12] = mhd_config.dt_out
    int_data[0] = mhd_config.nx
    int_data[1] = mhd_config.ny
    int_data[2] = mhd_config.nz
    int_data[3] = mhd_config.nx // mhd_config.mpi_sizex
    int_data[4] = mhd_config.ny // mhd_config.mpi_sizey
    int_data[5] = mhd_config.nz // mhd_config.mpi_sizez
    int_data[6] = mhd_config.mpi_sizex
    int_data[7] = mhd_config.mpi_sizey
    int_data[8] = mhd_config.mpi_sizez
    int_data[9] = 9
    int_data[10] = 0 # Periodic boundary condition as default

    fpath = run_directory + 'bin_data/'
    mkdir_p(fpath)
    fname = fpath + 'mhd_config.dat'
    double_data.tofile(fname)
    with open(fname, 'a') as file_handler:
        int_data.tofile(file_handler)


def main():
    """business logic for when running this module as the primary one!"""
    pass


if __name__ == "__main__":
    main()

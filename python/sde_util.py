"""
Utility functions
"""
import errno
import os

import numpy as np


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


if __name__ == "__main__":
    pass

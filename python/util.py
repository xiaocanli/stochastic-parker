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


if __name__ == "__main__":
    pass

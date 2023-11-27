Installation
=============

This package uses a `Multiple stream Mersenne Twister
PRNG <http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_en.html>`__
to generated the required random numbers. This package uses
`FLAP <https://github.com/szaghi/FLAP>`__ to deal with comand line
arguments. The FLAP pacakge is recommended to install using the
`FoBiS.py <https://github.com/szaghi/FoBiS>`__, which is a building system for Fortran projects. Currently, the code can only be run on CPUs.

Requirments
-----------

-  Programs that are commonly installed on HPC clusters
   
   -  **Git**: for cloning the code from GitHub
   -  **CMake**: 3.9.0 or higher
   -  **GCC** or **Intel** Fortran compilers.
   -  **MPI**: OpenMPI or MPICH. It should work with the others.
   -  **HDF5**: parallel version.
   .. tip::
      On Perlmutter@NERSC, you can load the modules by ``module load cpu cmake cray-hdf5-parallel``

-  `FoBiS.py <https://github.com/szaghi/FoBiS>`__: a building system for Fortran projects

   -  The installation wiki:
      https://github.com/szaghi/FoBiS/wiki/Install.
   -  It is recommended to use PyPI to install it:
      ``pip install FoBiS.py --user``
   -  After installing, make sure that ``FoBiS.py`` is in your ``$PATH`` by including
      
      .. code:: sh

         export PATH=$HOME/.local/bin/:$PATH

      in your ``.bashrc``.

-  `Multiple stream Mersenne Twister
   PRNG <http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_en.html>`__: for generating pseudorandom numbers parallelly.

      This program module splits single long period random number series from Mersenne Twister (MT) into multiple (almost) independent streams.

   -  Download:
      http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_f90.tar.gz
   -  Load the same packages as listed above.

      .. code-block:: bash

         tar zxvf mt_stream_f90.tar.gz
         cd mt_stream_f90-1.11
      
      Depending on the system type, you might need to modify ``FC`` in ``Makefile`` by hand. The default is `ifort`, which is for the Intel compilers. For `gfortran`, please use the script ``config/gfortran_mt_stream.sh`` to make the relevant changes in ``Makefile``. First, copy the script to ``mt_stream_f90-1.11``. Then,

      .. code:: sh

         ./gfortran_mt_stream.sh

      .. note::
         If you are using newer versions of `gfortran`, there will be compiling errors. Therefore, we need to update the source code slightly using ``config/fix_mt_stream.sh``. First, copy the script to ``mt_stream_f90-1.11``. Then,

         .. code:: sh

            ./fix_mt_stream.sh

   -  You can then compile the code by
  
      .. code:: sh

         make

   - After compiling the code, you need to set the environment variable ``MT_STREAM`` to the installation directory. For example,

      .. code:: sh

         export MT_STREAM=$HOME/local/mt_stream_f90-1.11

     This environment variable will be used to compile the GPAT code.

     .. tip::
        You can put this into your ``.bashrc``.

Download
--------

.. code:: sh

   git clone https://github.com/xiaocanli/stochastic-parker 

Install
-------

In the directory ``stochastic-parker``,

.. code:: sh

   mkdir build
   cd build
   cmake ..
   make
   make install

To turn on OpenMP parallelization, please use ``cmake -DUSE_OPENMP="On" ..``. To turn on ``AVX512`` for the KNL nodes, please also include ``-DUSE_AVX512="On"``.

.. tip::
 On Perlmutter@NERSC, you can load the modules by ``module load cpu cmake cray-hdf5-parallel python``

.. note::
   Please use the same version of ``python`` as when installing ``FoBiS.py``. Otherwise, ``FLAP`` will not be compiled correctly.

.. note::
   The code is not carefully optimized for the KNL nodes. Please use it in caution.

You can create a soft link of the executable ``stochastic-mhd.exec`` in the scratch filesystem for running the code.
Install Guide
=============

This package uses a `Multiple stream Mersenne Twister
PRNG <http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_en.html>`__
to generated the required random numbers. This package uses
`FLAP <https://github.com/szaghi/FLAP>`__ to deal with comand line
arguments. The FLAP pacakge is recommended to install using the
`FoBiS <https://github.com/szaghi/FoBiS>`__, which is a building system
for Fortran projects.

Requirments
-----------

-  Programs that are commonly installed on HPC clusters
   
   -  **Git**: for cloning the code from GitHub
   -  **CMake**: 3.9.0 or higher
   -  **GCC** or **Intel** Fortran compilers.
   -  **MPI**: OpenMPI or MPICH. It should work with the others.
   -  **HDF5**: parallel version.
-  `FoBiS.py <https://github.com/szaghi/FoBiS>`__: a building system for Fortran projects

   -  The installation wiki:
      https://github.com/szaghi/FoBiS/wiki/Install.
   -  It is recommended to use PyPI to install it:
      ``pip install FoBiS.py --user``
   -  After installing, make sure that ``FoBiS.py`` is in your ``$PATH``

-  `Multiple stream Mersenne Twister
   PRNG <http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_en.html>`__: for generating pseudorandom numbers parallelly.

      This program module splits single long period random number series from Mersenne Twister (MT) into multiple (almost) independent streams.

   -  Download:
      http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_f90.tar.gz
   -  Load the same packages as listed above.

      .. code:: sh

         tar zxvf mt_stream_f90.tar.gz
         cd mt_stream_f90-1.11
         make

   -  Its ``Makefile`` uses Intel compiler as default. To use GCC
      compilers, you have to modify the ``Makefile``.
   -  After installation, set the environment variable ``MT_STREAM`` to
      the installation directory. For example,

      .. code:: sh

         export MT_STREAM=$HOME/local/mt_stream_f90-1.11

      .. note::

       On Cori@NERSC, please set ``MT_STREAM_KNL`` or ``MT_STREAM_HSW``, depending on whether KNL nodes (``USE_AVX512`` is used during compiling ``stochastic-parker``) or Haswell nodes are used.

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

.. note::
   The code is not carefully optimized for the KNL nodes. Please use it in caution.

You can link the executable ``stochastic-mhd.exec`` to the scratch filesystem for running the code.
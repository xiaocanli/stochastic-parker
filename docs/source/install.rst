Install Guide
=============

Overview
--------

This package uses a `Multiple stream Mersenne Twister
PRNG <http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_en.html>`__
to generated the required random numbers. This package uses
`FLAP <https://github.com/szaghi/FLAP>`__ to deal with comand line
arguments. The FLAP pacakge is recommended to install using the
`FoBiS <https://github.com/szaghi/FoBiS>`__, which is a building system
for Fortran projects.

Requirments
-----------

-  **Git**. The newer the better.
-  **CMake** 3.9.0 or higher
-  **GCC** or **Intel** compilers. Not sure if it works for different
   versions.
-  **MPI**. OpenMPI or MPICH. Not fully tested.
-  `FoBiS <https://github.com/szaghi/FoBiS>`__

   -  The installation wiki:
      https://github.com/szaghi/FoBiS/wiki/Install.
   -  It is recommended to use PyPI to install it:
      ``pip install FoBiS.py --user``

-  `Multiple stream Mersenne Twister
   PRNG <http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_en.html>`__

   -  Download:
      http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_f90.tar.gz
   -  Load the same packages as listed above.

   .. code:: sh

      $ tar zxvf mt_stream_f90.tar.gz
      $ cd mt_stream_f90-1.11
      $ make

   -  Its ``Makefile`` uses Intel compiler as default. To use GCC
      compilers, you have to modify the ``Makefile``.
   -  After installation, set the environment variable ``MT_STREAM`` to
      the installation directory. For example,

   .. code:: sh

      export MT_STREAM=$HOME/local/mt_stream_f90-1.11

Download
--------

.. code:: sh

   $ git clone https://github.com/xiaocanli/stochastic-parker 

Install
-------

-  In the directory ``stochastic-parker``,

.. code:: sh

   $ mkdir build
   $ cd build
   $ cmake ..
   $ make
   $ make install

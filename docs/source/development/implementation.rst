Implementation Details
======================

This section provides an overview of the code implementation details. The code solves particle transport equations using stochastic differential equations (SDEs). To solve the SDEs, the code uses pseudo particles to sample the particle distributions. The code is designed to be modular and scalable. It is parallelized using MPI and OpenMP and can handle large-scale 3D simulations. The code is organized into several modules, each of which is responsible for a specific aspect of the simulation. See the next section for more details. In general, the code needs to read MHD data, interpolate it at particle positions, evolve the particles, and produce diagnostics.

- Initialization Phase:

  - Set up the simulation parameters (``simulation_setup_module``, ``mhd_config``, and ``mpi_module``).
  - Read the MHD data from files (``mhd_data_parallel``)
  - Initialize the particle data (``particle_module``).
  - Initialize diagnostics (``diagnostics``).

- Main Simulation Loop:

  - Inject pseudo particles (``particle_module``).
  - Interpolate MHD data to particle positions (``mhd_data_parallel``).
  - Evolve particles using SDEs and MHD fields (``particle_module``, ``random_number_generator``).

- Diagnostics and Data Handling:
  
  - Produce diagnostics of the simulation (``diagnostics``).
  - Handle data read/write operations (``hdf5_io``, ``mpi_io``).

Module Structure
----------------

The code is organized into several modules to improve modularity and maintainability. The modules are saved in ``src/modules`` directory. These modules are used by ``src/stochastic_mhd.f90`` to perform the simulation. The main modules include:

.. note:: 
    This section needs to be updated to include more details about each of these modules.

- ``particle_module``: This module contains particle data and the methods to inject, remove, and push particles. All the core algorithms of SDE are implemented in this module.
- ``mhd_data_parallel``: This module contains the data structures and methods to initialize and read MHD simulation data. It also contains the methods to calculate the gradients of these fields and to interpolate the MHD data to the particle positions. Additionally, this module also includes the methods for spatially dependent turbulence correlation length and turbulence amplitude.
- ``simulation_setup_module``: This module contains the methods to set up the simulation parameters, such as the simulation MPI topology, the MHD field configuration, particle boundary conditions, and neighbors of local MPI domains.
- ``hdf5_io``: This module contains the methods to create, open, close, read and write data in HDF5 format. It can deal with data up to 5D and with/without parallel I/O.
- ``mhd_config``: This module contains the methods to read the MHD configuration file and to set up the MHD simulation parameters, including the grid and time stamps if needed.
- ``constants``: This module contains constants used in the simulation.
- ``diagnostics``: This module contains the methods to diagnose the simulations, such as the particle distribution functions, particle data, escaped particle data, and methods for quick checks as the simulation proceeds.
- ``mpi_io``: This module contains the methods to read and write data in parallel using MPI-IO.
- ``mpi_module``: This module contains the methods to set up the MPI topology and to communicate between different MPI domains.
- ``random_number_generator``: This module contains the methods to generate random numbers using the Mersenne Twister algorithm.
- ``read_confg``: This module contains the methods to read the simulation configuration file.
- ``acc_region_surface``: This module contains the methods for the surface to separate the acceleration region and the non-acceleration region. It is useful when we want to isolate certain particle acceleration regions, for example, near the shock front.
- ``mhd_data_sli``: This module contains the methods to read the MHD simulation data from Shengtai Li's MHD code. NOTE: This module is not used in the current version of the code and will be removed in the future.
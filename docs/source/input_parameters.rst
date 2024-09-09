Input Parameters
================

The input parameters are mostly in script ``diffusion*.sh`` (e.g., ``diffusion_reconnection.sh`` for the 2D reconnection example). A few others are in ``conf*.dat`` (e.g., ``conf.dat`` for the 2D reconnection example). The input parameters are explained below.

diffusion.sh
------------

- ``MACHINE``: either ``Frontera`` or ``Perlmutter``. ``Perlmutter`` is the default.
- ``quota_hour``: the maximum wall time in hours. The code dump restart files if the simulation is not finished when the wall time is reached (30 mins before that in the code). The restart files are saved in the ``restart`` sub-directory under the diagnostic directory. These are data files ``particles_*.h5`` for particle data, ``mtss.mpi_rank.thread_id`` for the state of the random number generators, ``particle_module_stare_*.h5`` for the state of the particle module, and ``latest_restart`` for the last time frame before duming restart files.
- ``restart_flag``: whether the simulation is a restart from previous simulations. ``.false.`` is the default, which means it is not a restart. ``.true.`` means it is a restart. When it is a restart, the code will read the particle data and ``random_number_generater`` from the previous simulation and continue the simulation from the last time step.
- ``conf``: the name of the configuration file. ``conf.dat`` is the default.
- ``mhd_config_filename``: the name of the MHD configuration file. ``mhd_config.dat`` is the default.
- ``focused_transport``: whether to solve the focused transport equation. ``.false.`` is the default, which means solving the Parker transport equation.
- ``mpi_size``: the total number of MPI processes.
- ``ntasks_per_node``: the number of MPI processes per node.
- ``size_mpi_sub``: the size of an MPI sub-communicator. ``1`` is the default, which means there is only one MPI communicator (MPI_COMM_WORLD). If not ``1``, the total number of MPI processes must be divisible by ``size_mpi_sub``. The MPI communications are done within each sub-communicator. In this way, we don't need to divide the domain into too many sub-domains, which is not efficient for MPI communications. To get the global properties (e.g., global particle spectrum), we will use the cross-communicators.
- ``nptl``: the total number of particles for one MPI process.
- ``nptl_max``: the maximum number of particles for one MPI process. If the number of particles exceeds ``nptl_max``, we will stop injecting more particles or splitting particles. Typically, we can choose up to 10 million.
- ``ts``: the start time step.
- ``te``: the end time step.
- ``tmax_mhd``: the maximum time frame for the MHD fields. Typically, ``tmax_mhd`` is larger than ``te``, so the transport modeling will have enough MHD fields to use. When ``tmax_mhd`` is smaller than ``te``, there are two cases: 1) when ``single_time_frame=1`` (discussed below), the code will keep running until ``te``; 2) when ``single_time_frame=0``, the code will stop reading MHD fields after ``tmax_mhd`` and use the last MHD fields for the rest of the simulation.
- ``single_time_frame``: whether to use only one-time frame. ``0`` is the default, which means using multiple time frames (from ``ts`` to ``te``). ``1`` means using only one time frame (``ts``). Then, the simulation will run ``te-ts+1`` time steps using the same MHD fields.
- ``time_interp``: whether to interpolate in-between time frames. ``1`` is the default, which means interpolation. ``0`` means no interpolation.
- ``dist_flag``: the initial particle distribution. ``0`` is the default, which means Maxwellian. ``1`` means delta function. ``2`` means power-law.
- ``power_index``: the power-law index for the initial momentum distribution if ``dist_flag=2``.
- ``split_flag``: whether to split particles when they reach certain momentum. ``0`` is the default, which means no splitting. ``1`` means splitting.
- ``split_ratio``: the momentum increase ratio for particle splitting. Assuming the particle's initial momentum is ``p0``, the particle will be splitted into two particles once the particle momentum reaches ``p0*split_ratio``. The rul of thumb is to set ``split_ratio`` to be around :math:`2^{1/\gamma}`, where :math:`\gamma` is the power-law index for the particle momentum distribution.
- ``pmin_split``: the minimum momentum (in ``p0``) to start splitting particles. Instead of starting from ``p0``, we can start from a higher momentum to avoid splitting too many particles at low momentum. It is particularly important when injecting Maxwellian particles, where the number of particles is large at low momentum.
- ``local_dist``: whether to diagnose local particle distribution.
- ``dump_escaped_dist``: whether to dump the escaped particle distributions, either ``.true`` or ``.false``. If ``.false.``, the particles will be removed once they escape from the boundaries. Otherwise, the escaped particles will be accumulated in one MHD data interval. The distributions of the escaped particles are stored in files ``escaped_dist_*.h5``. They have the same grid sizes as the local particle distributions but with one dimension less. For example, if we have a 3D simulation, the escaped particle distributions will be 2D for the low-x, high-x, low-y, high-y, low-z, and high-z boundaries.

.. warning::
    The escaped particle distributions are only dumped when ``local_dist=.true.``. When ``local_dist=.false.``, the escaped particle distributions will not be dumped event if ``dump_escaped_dist=.true.``.

- ``dump_escaped``: whether to dump the escaped particles. ``dump_escaped_dist`` has to be ``.true.`` for this to work. If ``dump_escaped_dist=.false.``, the raw data of the escaped particles will be dumped. Otherwise, only the distributions of the escaped particles will be dumped. Note that the raw data of the escaped particles can be large. It is not recommended to dump the raw data of the escaped particles unless you want to analyze the raw data in detail.
- ``track_particle_flag``: whether to track particles. ``.false.`` is the default, which means not tracking particles. ``.true.`` means tracking particles. When tracking particles, we will select ``nptl_selected`` particles with the highest energy, rerun the simulation, and output the trajectories of these particles.
- ``particle_tags_file``: the file name for the particle tags. It is used to select the particles to track. The file ``tags_selected_*.h5`` has one dataset ``tags``. The dataset is a 1D array of integers. These tags can be obtained from analyzing the particle data. For example, we can select the particles with the highest energy or particles in certain regions and output their tags to the file ``tags_selected_*.h5``. The file is located in the same directory as the diagnostic directory.
- ``nsteps_interval``: the steps interval to track particles. If it is ``1``, we will track particles at every time step. If it is ``10``, we will track particles at every 10 steps. Particles can be evolves for millions of time steps. It is not recommended to track particles at every time step unless you want to analyze the particle trajectories in great detail.
- ``inject_new_ptl``: whether to inject new particles at every step. The values are ``.true.`` or ``.false.``. ``.true.`` means injecting new particles at every MHD step. ``.false.`` means only injecting particles at the beginning.
- ``inject_same_nptl``: whether to inject the same number of particles every step. ``.true.`` is the default, which means injecting the same number of particles. ``.false.`` means injecting a different number of particles every step. Then, the number of injected particles depends on how particles are injected (see below). ``inject_same_nptl`` only works when ``inject_new_ptl=.true.``.
- ``tmax_to_inject``: the maximum time frame to inject particles. It is used to stop injecting particles after some time. It is useful, for example, when we want to run simulations for a long time but only inject particles at the beginning. Note that ``tmax_to_inject`` should be larger than ``ts``. Note that it only works when ``inject_new_ptl=.true.``.
- ``inject_part_box``: whether to inject particles in part of the simulation box. ``.false.`` is the default, which means injecting particles in the whole box. ``.true.`` means injecting particles in part of the box.
- ``ptl_xmin``: the minimum x coordinate for particle injection. It is the absolute coordinate.
- ``ptl_xmax``: the maximum x coordinate for particle injection.
- ``ptl_ymin``: the minimum y coordinate for particle injection.
- ``ptl_ymax``: the maximum y coordinate for particle injection.
- ``ptl_zmin``: the minimum z coordinate for particle injection.
- ``ptl_zmax``: the maximum z coordinate for particle injection.
- ``inject_large_jz``: whether to inject particles where jz is large. ``.false.`` is the default, which means not injecting particles where jz is large. ``.true.`` means injecting particles where jz is large.
- ``jz_min``: the minimum jz for injection when ``inject_large_jz=.true.``. It is to make sure that the particles are only injected in regions with large jz. It is useful for reconnection simulations, where we might want to inject particles in the reconnection region. Depending on the simulation, we might need to adjust this value because each simulation has a different jz range or different normalizations. Note that ``inject_large_jz=.true.``, ``inject_large_absj=.true.``, ``inject_large_db2=.true.``,  ``inject_large_divv=.true.``, and ``inject_large_rho=.true.`` are exclusive. Only one of them can be ``.true.``.
- ``ncells_large_jz_norm``: the normalization for the number of cells with large jz. It is used to determine the number of particles to inject. The number of particles to inject is ``nptl*ncells(jz>jz_min)/ncells_large_jz_norm``, where ``ncells(jz>jz_min)`` is the number of cells with jz larger than ``jz_min``.
- ``inject_large_absj``: whether to inject particles where :math:`|j|` is large. ``.false.`` is the default, which means not injecting particles where :math:`|j|` is large. ``.true.`` means injecting particles where :math:`|j|` is large.
- ``absj_min``: the minimum :math:`|j|` for injection when ``inject_large_absj=.true.``. It is to make sure that the particles are only injected in regions with large :math:`|j|`.
- ``ncells_large_absj_norm``: the normalization for the number of cells with large :math:`|j|`. It is used to determine the number of particles to inject. The number of particles to inject is ``nptl*ncells(|j|>absj_min)/ncells_large_absj_norm``, where ``ncells(|j|>absj_min)`` is the number of cells with :math:`|j|` larger than absj_min.
- ``inject_large_db2``: whether to inject particles where turbulence amplitude (db2) is large. ``.false.`` is the default, which means not injecting particles where db2 is large. ``.true.`` means injecting particles where db2 is large. When ``inject_large_db2=.true.``, the code needs to read an additional file ``deltab*.dat`` to get the turbulence amplitude.
- ``db2_min``: the minimum db2 for injection when ``inject_large_db2=.true.``. It is to make sure that the particles are only injected in regions with large db2. It is useful for reconnection simulations, where we might want to inject particles in the reconnection region, where reconnection-driven turbulence can be intense.
- ``ncells_large_db2_norm``: the normalization for the number of cells with large db2. It is used to determine the number of particles to inject. The number of particles to inject is ``nptl*ncells(db2>db2_min)/ncells_large_db2_norm``, where ``ncells(db2>db2_min)`` is the number of cells with db2 larger than ``db2_min``.

.. note:: 
    The functionality of ``inject_large_db2`` is not fully tested. It is not recommended to use it for now. Additionally, we need to understand spatially dependent turbulence amplitude better.

- ``inject_large_divv``: whether to inject particles where flow compression divv is negatively large. ``.false.`` is the default, which means not injecting particles where divv is negatively large. ``.true.`` means injecting particles where divv is negatively large.
- ``divv_min``: the minimum divv for injection when ``inject_large_divv=.true.``. This is to make sure that the particles are only injected in regions with large volumes. It is useful for reconnection or shock simulations, where we might want to inject particles near regions with strong flow compression.
- ``ncells_large_divv_norm``: the normalization for the number of cells with large ``divv``. It is used to determine the number of particles to inject. The number of particles to inject is ``nptl*ncells(|divv|>|divv_min|)/ncells_large_divv_norm``, where ``ncells(|div|>|divv_min|)`` is the number of cells with ``|divv|`` larger than ``|divv_min|``.
- ``inject_large_rho``: whether to inject particles where density is large. ``.false.`` is the default, which means not injecting particles where density is large. ``.true.`` means injecting particles where density is large.
- ``rho_min``: the minimum density for injection when ``inject_large_rho=.true.``. It is to make sure that the particles are only injected in regions with large density.
- ``ncells_large_rho_norm``: the normalization for the number of cells with large density. It is used to determine the number of particles to inject. The number of particles to inject is ``nptl*ncells(rho>rho_min)/ncells_large_rho_norm``, where ``ncells(rho>rho_min)`` is the number of cells with density larger than ``rho_min``.
- ``dpp_wave``: whether to include momentum diffusion due to wave scattering. ``0`` is the default, which means not including momentum diffusion due to wave scattering. ``1`` means including momentum diffusion due to wave scattering.
- ``dpp_shear``: whether to include momentum diffusion due to flow shear. ``0`` is the default, which means not including momentum diffusion due to flow shear. ``1`` means including momentum diffusion due to flow shear.
- ``weak_scattering``: whether particle scattering is in the weak-scattering regime. ``1`` is the default, which means the weak-scattering regime. ``0`` means the strong-scattering regime.
- ``deltab_flag``: whether to have spatially dependent turbulence amplitude. ``0`` is the default, which means that the turbulence amplitude is spatially uniform. ``1`` means having spatially dependent turbulence amplitude. When ``deltab_flag=1``, the code needs to read an additional file ``deltab*.dat`` to get the turbulence amplitude.
- ``correlation_flag``: whether to have spatially dependent turbulence correlation length. ``0`` is the default, which means that the turbulence correlation length is spatially uniform. ``1`` means having spatially dependent turbulence correlation length. When ``correlation_flag=1``, the code needs to read an additional file ``lc*.dat`` to get the turbulence correlation length.

.. note:: 
    The functionalities of ``deltab_flag`` and ``corrlation_flag`` are not fully tested. It is not recommended to use them for now. Additionally, we need to understand spatially dependent turbulence amplitude and correlation length better.

- ``ndim_field``: the dimension of the field. The values can be ``1``, ``2``, or ``3``. ``1`` means 1D simulation, but it is not fully tested.
- ``drift_param1``: the parameter 1 for particle drift. It is used to determine the drift velocity. See the Theory section for details.
- ``drift_param2``: the parameter 2 for particle drift. It is used to determine the drift velocity. See the Theory section for details.
- ``charge``: the charge of the particle in unit charge. ``-1`` is the default, which means electron.
- ``spherical_coord``: whether the grid is spherical. ``0`` is the default, which means the grid is Cartesian. ``1`` means the grid is spherical.
- ``uniform_grid``: whether the grid is uniform. ``1`` is the default, which means the grid is uniform. ``0`` means the grid is non-uniform. Then, we need to the coordinates of the grid points in files ``xpos.dat``, ``ypos.dat``, and ``zpos.dat``, located in the same directory as the MHD configuration file.
- ``check_drift_2d``: whether to check particle drift in 2D simulations. ``0`` is the default, which means not checking particle drift in 2D simulations. ``1`` means checking particle drift in 2D simulations. It is useful for 2D simulations, where we can check how much particles drift along the out-of-plane direction.
- ``particle_data_dump``: whether to dump particle data. ``0`` is the default, which means not dumping particle data. ``1`` means dumping particle data. When dumping particle data, the code will output the particle data at every output time step. The particle data are stored in files ``particles_*.h5``. Since particle data can be large, it is not recommended to dump particle data unless you want to analyze the particle data.
- ``include_3rd_dim``: whether to include transport along the 3rd-dim in 2D simulations. ``0`` is the default, which means not including transport along the 3rd-dim in 2D simulations. ``1`` means including transport along the 3rd-dim in 2D simulations. It is useful for 2D simulations, where we can check how much particles transport along the out-of-plane direction.
- ``acc_by_surface``: whether the acceleration region is separated by a surface. ``0`` is the default, which means the acceleration region is not separated by a surface. ``1`` means the acceleration region is separated by a surface. Then, we need to specify the surface file name and the normal direction of the surface to get the 2D surfaces that separate the acceleration region. The surface file name is specified by ``surface_filename1`` and ``surface_filename2``. The normal direction of the surface is specified by ``surface_norm1`` and ``surface_norm2``. The surface normal direction can be ``+x``, ``-x``, ``+y``, ``-y``, ``+z``, or ``-z``. It is useful when we want to selectively turn on/off particle acceleration in certain regions. For example, we can turn on particle acceleration only in the reconnection region or termination shock region in flare simulations. To get the 2D surfaces separating different acceleration regions, we need to write our own scripts to look into the MHD simulation data.

.. note:: 
    The functionality of ``acc_by_surface`` is not fully tested. It is not recommended to use it for now. Example scripts to get the 2D surfaces separating different acceleration regions will be provided later.

- ``surface2_existed``: whether the second surface exists. ``0`` is the default, which means the second surface does not exist. Then, we only have two regions. ``1`` means the second surface exists. Then, we will have acceleration regions separated by these two surfaces.
- ``varying_dt_mhd``: whether the time interval for MHD fields is varying. ``0`` is the default, which means the time interval for MHD fields is uniform. ``1`` means the time interval for MHD fields is varying. It is useful when the MHD simulation has varying time intervals. For example, the MHD simulation might have a large time interval at the beginning and a small time interval later. Then, we can use ``varying_dt_mhd=1`` to use the varying time interval for MHD fields. When ``varying_dt_mhd=1``, the code needs to read an additional file ``time_stamps.dat`` to get the time stamps for each MHD frame. The file is located in the same directory as the MHD data files.

Then, the script will modify the configuration file ``conf.dat``. The parameters in ``conf.dat`` are explained below. Additionally, a few other parameters are modified in this script for more flexibility.

- ``tau0_scattering``: the scattering time for initial particles. It is only used for momentum diffusion due to wave scattering. It is not used for Parker transport. The parameters are calculated based on the initial particle momentum and turbulence properties in ``sde.py``.
- ``duu0``: the normalization for pitch-angle diffusion coefficient. It is only used in the focused transport equation. The parameters are calculated based on the initial particle momentum and turbulence properties in ``sde.py``.
- ``particle_v0``: the particle speed/velocity normalization. It is only used in the focused transport equation. The parameters are calculated based on the initial particle momentum and turbulence properties in ``sde.py``.
- ``dir_mhd_data``: the directory for MHD simulation data.
- ``diagnostic_directory``: the directory for diagnostics data.

conf.dat
--------

- ``b0``: initial magnetic field strength (deprecated).
- ``p0``: initial particle momentum. Its value is arbitrary. 0.1 is typically used so that the particle momentum is not too small or too large. Note that ``p0`` corresponds to particles with the input diffusion coefficients.
- ``pmin``: the minimum particle momentum. It is used when injecting particles and when calculating the global particle spectrum. It is typically set to ``1E-2``.
- ``pmax``: the maximum particle momentum. It is used when injecting particles and when calculating the global particle spectrum. It is typically set to ``1E1``.
- ``momentum_dependency``: whether the diffusion coefficients depend on particle momentum. ``1`` is the default, which means the diffusion coefficients depend on particle momentum. ``0`` means the diffusion coefficients do not depend on particle momentum.
- ``pindex``: the power-law index for the momentum dependency of the diffusion coefficients. It is only used when ``momentum_dependency=1``. It is typically set to ``3-5/3=4/3=1.3333333``, where ``5/3`` is the turbulence spectral slope for the Kolmogorov spectrum. It can be modified in ``difffusion.sh`` when using different turbulence models.
- ``mag_dependency``: whether the diffusion coefficients depend on magnetic field strength. ``1`` is the default, which means the diffusion coefficients depend on magnetic field strength. ``0`` means the diffusion coefficients do not depend on magnetic field strength.
- ``kpara0``: the normalization for the parallel diffusion coefficient. It is calculated based on the initial particle momentum, magnetic field, and turbulence properties in ``sde.py``.
- ``kret``: the ratio of the perpendicular diffusion coefficient to the parallel diffusion coefficient. It is typically set to less than ``0.1``.
- ``dt_min``: the minimum time step allowed to avoid infinite time step.
- ``dt_min_rel``: the minimum relative time step w.r.t. one field time interval. ``dt_min`` is set to ``dt_min_rel`` times the time interval for MHD fields if the latter is larger than ``dt_min``.

.. note::
    The time step is adaptive. It is calculated based on the particle momentum, magnetic field, pitch angle, and diffusion coefficients. The rule of thumb for ``dt_min_rel`` is ``1E-6``. For MHD simulations with lower resolution, it can be up to ``1E-4``. For MHD simulations with very high resolutions, a large ``dt_min_rel`` might lead to wrong results in these high-resolution simulations, while a small ``dt_min_rel`` might lead to a long simulation time. We suggest doing a convergence test to get the optimal value. 

- ``dt_max_rel``: the maximum relative time step w.r.t. one field time interval to avoid a time step too large, which could cause the particles to jump over multiple grid cells.
- ``npp_global``: the number of momentum bins for the global particle spectrum.
- ``nmu_global``: the number of pitch-angle bins for global particle distributions.
- ``dump_interval1``: the interval to dump local particle distributions. It is only used when ``local_dist=1`` in ``diffusion.sh``.
- ``pmin1``: the minimum particle momentum for local particle distributions.
- ``pmax1``: the maximum particle momentum for local particle distributions.
- ``npbins1``: the number of momentum bins for local particle distributions.
- ``nmu1``: the number of pitch-angle bins for local particle distributions.
- ``rx1``: reduced factor along the x direction for local particle distributions. For every ``rx1`` grid cell along the x direction, we will have one bin for local particle distributions.
- ``ry1``: reduced factor along the y direction for local particle distributions. For every ``ry1`` grid cell along the y direction, we will have one bin for local particle distributions.
- ``rz1``: reduced factor along the z-direction for local particle distributions. For every ``rz1`` grid cell along the z direction, we will have one bin for local particle distributions.

.. note::
    The other three local distributions are similar. We can adjust the number of bins and reduce factors to get different distributions. For example, we can get a distribution with higher resolution in the momentum space and lower resolution in the pitch-angle space by increasing ``npbins`` and decreasing ``nmu``. Or we can get distributions with higher momentum resolution but coarse spatial resolution by increasing ``rx``, ``ry``, and ``rz``.
    
.. note::
    We only dump local distributions every few MHD output intervals. When ``dump_interval`` is larger than the number of MHD outputs, it will not dump the distribution. In this way, we don't have to dump all four kinds of local distributions.

- ``acc_region_flag``: whether to turn on particle acceleration in certain regions. ``0`` is the default, which means turning on particle acceleration in the entire region. ``1`` means turning on particle acceleration in certain regions. When ``acc_region_flag=1``, we need to specify the acceleration region. The acceleration region is specified by ``acc_xmin``, ``acc_xmax``, ``acc_ymin``, ``acc_ymax``, ``acc_zmin``, and ``acc_zmax``. These are the relative values from 0 to 1. The acceleration region is a box with the minimum coordinate (``acc_xmin``, ``acc_ymin``, ``acc_zmin``) and the maximum coordinate (``acc_xmax``, ``acc_ymax``, ``acc_zmax``). It is useful when we want to selectively turn on/off particle acceleration in certain regions. For example, we can turn on particle acceleration only in the reconnection region or termination shock region in flare simulations. If we set ``acc_xmax`` or ``acc_ymax`` or ``acc_zmax`` to negative values, the acceleration in the entire simulation domain will be turned off.
- ``pbcx``: the boundary condition for particles along the x direction. ``0`` is the default, which means periodic boundary condition. ``1`` means open boundary condition.
- ``pbcy``: the boundary condition for particles along the y direction. ``0`` is the default, which means periodic boundary condition. ``1`` means open boundary condition.
- ``pbcz``: the boundary condition for particles along the z direction. ``0`` is the default, which means periodic boundary condition. ``1`` means open boundary condition.

.. note:: 
    Additional boundary conditions should be included in the future, such as reflecting boundary condition.

- ``mpi_sizex``: the number of MPI processes along the x direction. It is default to ``1`` when ``size_mpi_sub=1``. Otherwise, ``mpi_sizex*mpi_sizey*mpi_sizez`` should be equal to ``size_mpi_sub``.
- ``mpi_sizey``: the number of MPI processes along the y direction. It is default to ``1`` when ``size_mpi_sub=1``. Otherwise, ``mpi_sizex*mpi_sizey*mpi_sizez`` should be equal to ``size_mpi_sub``.
- ``mpi_sizez``: the number of MPI processes along the z direction. It is default to ``1`` when ``size_mpi_sub=1``. Otherwise, ``mpi_sizex*mpi_sizey*mpi_sizez`` should be equal to ``size_mpi_sub``.

.. note:: 
    When ``size_mpi_sub>1`` in ``diffusion.sh``. ``mpi_sizex*mpi_sizey*mpi_sizez`` should be equal to ``size_mpi_sub``. Otherwise, the code will stop.

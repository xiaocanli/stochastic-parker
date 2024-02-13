Particle Module
===============

This section provides an overview of the particle module, including the particle data structure, the methods for injecting, removing, and pushing particles.

Particle Data Structure
-----------------------

.. code-block:: fortran

    type particle_type
        integer(i1)  :: split_times     !< Particle splitting times
        integer(i1)  :: count_flag      !< Only count particle when it is 1
        integer(i4)  :: origin          !< The origin MPI rank of the particle
        integer(i4)  :: nsteps_tracked  !< # of tracked steps in MHD data interval
        integer(i4)  :: nsteps_pushed   !< # of steps have been pushed
        integer(i4)  :: tag_injected    !< Particle tag 1 when injected
        integer(i4)  :: tag_splitted    !< Particle tag 2 after splitting
        real(dp) :: x, y, z, p          !< Position and momentum
        real(dp) :: v, mu               !< Velocity and cosine of pitch-angle
        real(dp) :: weight, t, dt       !< Particle weight, time and time step
    end type particle_type

- ``x``, ``y``, ``z``: position of the particle in the simulation domain. The positions are in the global coordinate system. For spherical coordinates, ``x`` is the radial distance :math:`r`, ``y`` is the polar angle :math:`\theta`, and ``z`` is the azimuthal angle :math:`\phi`.
- ``p``: momentum magnitude of the particle.
- ``v``: velocity magnitude of the particle. We typically assume non-relativistic particles, so that :math:`v \sim p`.
- ``mu``: cosine of the pitch angle of the particle.
- ``weight``: weight of the particle. The weight is initialized to be 1.0. When a particle reaches certain momentum, it will be splitted into two identical particles with half of the weight. Due to stochastic nature of the process, these particles will evolve differently after the splitting.
- ``t``: time of the particle. The time is initialized to be the time when particle is injected.
- ``dt``: time step of the particle. The time step is adaptive and is updated at each time step.
- ``split_times``: number of times the particle has been splitted.
- ``count_flag``: flag to count the particle. If the flag is 1, the particle will be counted when accumulating particle distributions. When the flag is negative (-1 to -6), it means that the particle has escaped from one of the six boundaries of the simulation domain (for 2D, -1 to -4 for the four boundaries). When the flag is 0, the particle has negative momentum and will be removed from the simulation.
- ``origin``: the origin MPI rank of the particle. Combined with particle ``tag_injected`` and ``tag_splitted``, it can be used to identify the particle for particle tracking.
- ``nsteps_tracked``: number of steps the particle has been tracked in the MHD data interval. It will be reset to 0 at the start of the MHD data interval.
- ``nsteps_pushed``: number of steps the particle has been pushed since they are injected. It is useful when tracking particles. It will accummulate to ``nsteps_interval`` and then reset to 0 once the tracked particle is recorded.
- ``tag_injected``: tag of the particle when it is injected. The tag is unique on each MPI rank.
- ``tag_splitted``: tag of the particle after it is splitted. The tag is unique on each MPI rank.

.. note:: 

    The tag of a particle only unique on each MPI rank. When combined with the origin MPI rank, the particle can be uniquely identified in the simulation. Since ``tag_splitted`` is a integer4, the maximum number of particles that one particle can be splitted into is 2^31. If the number of particles exceeds this limit, it means that we probably need to use a larger ``split_ratio``.

Particle Tracking
-----------------

As described above, we can track particles using the combination of ``tag_injected``, ``tag_splitted``, and ``origin``. The particle tracking is useful for diagnosing how the particles are accelerated and transported. To get the particle tracking information, we need to runt the simulation twice. The first run is a regular run. We need to run the simulation to a certain time frame and dump the particles (also included in the restart files). Then, we need to use a post-processing tool (e.g., a Python script) to select the particles we want to track and save the particle information to a file. The second run is a tracking run. We need to set ``track_particle_flag = .true.`` in the input file and set the ``particle_tags_file`` to the file that contains the particle information. The tracking run is similar to the regular run except that we need to load the ``particle_tags_file`` and use the information to identify the particles we want to track.

- When a particle is injected, the particle is assigned with an ``origin`` MPI rank and a ``tag_injected``, which starts from 0 and increases by 1 for each particle injected on each MPI rank. When the particle is splitted, the ``tag_splitted`` is assigned to the particle. ``tag_splitted`` starts from 1. When one particle is splitted into two particles, each particle has a half of the original weight. One of the two particles has the same ``tag_splitted`` as the original particle, and the other has the ``tag_splitted=tag_splitted_original+2**(split_times-1)``. The ``split_times`` is the number of times the particle has been splitted. The ``tag_splitted`` is used to identify the particle after it is splitted. The combination of ``tag_injected``, ``tag_splitted``, and ``origin`` can be used to uniquely identify the particle in the simulation.
- The first time we run the simulation, we need to set ``track_particle_flag = .false.`` in the input file. The particle information will be saved in the restart files. We can use the restart files to get the particle information. We can use a Python script to select the particles we want to track. For example, we can select the partilces with the highest energies form the particle file using the following script.

.. code-block:: python
    
    tframe = 200
    run_name = sde_run_config["run_name"]
    fname = "../data/" + run_name + "/restart/particles_" + str(tframe).zfill(4) + ".h5"
    with h5py.File(fname, 'r') as fh:
        p = fh["p"][:]
        tag_injected = fh["tag_injected"][:]
        tag_splitted = fh["tag_splitted"][:]
        origin = fh["origin"][:]
        split_times = fh["split_times"][:]
    sort_index = np.argsort(p)
    psort = p[sort_index]
    tag_injected_sort = tag_injected[sort_index]
    tag_splitted_sort = tag_splitted[sort_index]
    origin_sort = origin[sort_index]
    split_times_sort = split_times[sort_index]
    nptl_selected = 1000
    # select the highest-energy particles
    origin_s = origin_sort[-nptl_selected:]
    tag_injected_s = tag_injected_sort[-nptl_selected:]
    tag_splitted_s = tag_splitted_sort[-nptl_selected:]
    split_times_s = split_times_sort[-nptl_selected:]
    # trace back to the injected particle
    nsplit_max = np.max(split_times_s)
    tags_all = np.zeros([nptl_selected, nsplit_max+2], dtype=np.int32)
    tags_all[:, 0] = origin_s  # the first col is the origin
    tags_all[:, 1] = tag_injected_s  # the second col is tag_injected
    for iptl in range(nptl_selected):
        nsplit = split_times_s[iptl]
        if nsplit > 0:
            tags_all[iptl, nsplit+1] = tag_splitted_s[iptl]
            for isplit in range(nsplit, 1, -1):
                if tags_all[iptl, isplit+1] > 2**(isplit-1):
                    tags_all[iptl, isplit] = tags_all[iptl, isplit+1] - 2**(isplit-1)
                else:
                    tags_all[iptl, isplit] = tags_all[iptl, isplit+1]
    # sort by the origin, tag_injected, and then all tag_splitted
    ind = np.lexsort(tags_all[:, ::-1].T)
    tags_all_sorted = tags_all[ind]
    # save the selected particles to a file
    odir = "../data/" + run_name + "/particle_tracking/"
    mkdir_p(odir)
    fname = odir + "tags_selected_01.h5"
    with h5py.File(fname, 'w') as fh:
        fh.create_dataset("tags", tags_all_sorted.shape, data=tags_all_sorted)

- The second time we run the simulation, we need to set ``track_particle_flag = .true.`` in the input file. We need to set the ``particle_tags_file`` to the file that contains the particle information. The tracking run is similar to the regular run except that we need to load the ``particle_tags_file`` and use the information to identify the particles we want to track.
- The simulation will output tracked particle information every MHD intervel. The files are saved in ``particle_tracking`` under the diagnostic directory. We can get all the particle information using a script like the following.

.. code:: python

    fdir = "../data/" + run_name + "/particle_tracking/"
    ts, te = 51, 200
    nframes = te - ts + 1
    fname = fdir + "particles_tracked_" + str(ts).zfill(4) + ".h5"
    with h5py.File(fname, "r") as fh:
        x = fh["x"][:, :]
    nptl_tracking, nsteps_tracked = x.shape
    t = np.zeros([nptl_tracking, nsteps_tracked*nframes])
    x = np.zeros([nptl_tracking, nsteps_tracked*nframes])
    y = np.zeros([nptl_tracking, nsteps_tracked*nframes])
    p = np.zeros([nptl_tracking, nsteps_tracked*nframes])
    tag_splitted = np.zeros([nptl_tracking, nsteps_tracked*nframes])
    for tframe in range(ts, te+1):
        fname = fdir + "particles_tracked_" + str(tframe).zfill(4) + ".h5"
        it1 = (tframe - ts) * nsteps_tracked
        it2 = it1 + nsteps_tracked
        with h5py.File(fname, "r") as fh:
            t[:, it1:it2] = fh["t"][:, :]
            x[:, it1:it2] = fh["x"][:, :]
            y[:, it1:it2] = fh["y"][:, :]
            p[:, it1:it2] = fh["p"][:, :]
            tag_splitted[:, it1:it2] = fh["tag_splitted"][:, :]

To plot one particle trajectory, we can use, for example,

.. code:: python

    iptl = 100
    ind = t[iptl, :] > 0
    plt.plot(x[iptl, ind], y[iptl, ind])
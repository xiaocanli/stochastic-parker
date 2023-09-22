!*******************************************************************************
!< Module of particle data and methods to inject, remove and push particles
!*******************************************************************************
module particle_module
    use constants, only: fp, dp
    use simulation_setup_module, only: ndim_field
    use mhd_config_module, only: uniform_grid_flag, spherical_coord_flag
    use mhd_data_parallel, only: nfields, ngrads
    use mpi_module
    use omp_lib
    use hdf5
    implicit none
    private
    save
    public init_particles, free_particles, inject_particles_spatial_uniform, &
        read_particle_params, particle_mover, split_particle, &
        set_particle_datatype_mpi, free_particle_datatype_mpi, &
        select_particles_tracking, init_particle_tracking, &
        free_particle_tracking, init_tracked_particle_points, &
        free_tracked_particle_points, negative_particle_tags, &
        save_tracked_particle_points, record_tracked_particle_init, &
        inject_particles_at_shock, inject_particles_at_large_jz, &
        inject_particles_at_large_db2, inject_particles_at_large_divv, &
        set_dpp_params, set_duu_params, set_flags_params, set_drift_parameters, &
        set_flag_check_drift_2d, get_interp_paramters, &
        get_interp_paramters_spherical

    public particle_type, ptls, escaped_ptls, &
        nptl_current, nptl_escaped, nptl_escaped_max, nptl_max, &
        spherical_coord_flag, leak, leak_negp, nptl_split, &
        pmin, pmax, COUNT_FLAG_INBOX, COUNT_FLAG_OTHERS

    type particle_type
        real(dp) :: x, y, z, p      !< Position and momentum
        real(dp) :: v, mu           !< Velocity and cosine of pitch-angle
        real(dp) :: weight, t, dt   !< Particle weight, time and time step
        integer  :: split_times     !< Particle splitting times
        integer  :: count_flag      !< Only count particle when it is 1
        integer  :: tag             !< Particle tag
        integer  :: nsteps_tracking !< Total particle tracking steps
    end type particle_type

    real(dp) :: pmin  !< Minimum particle momentum
    real(dp) :: pmax  !< Maximum particle momentum

    integer, parameter :: COUNT_FLAG_INBOX  = 1  !< For in-box particles
    integer, parameter :: COUNT_FLAG_ESCAPE = -1 !< For escaped particles
    integer, parameter :: COUNT_FLAG_OTHERS = 0  !< For other particles

    integer :: particle_datatype_mpi
    type(particle_type), allocatable, dimension(:) :: ptls
    type(particle_type), allocatable, dimension(:, :) :: senders
    type(particle_type), allocatable, dimension(:, :) :: recvers
    integer, allocatable, dimension(:) :: nsenders, nrecvers
    !dir$ attributes align:128 :: ptls

    integer :: nptl_current         !< Number of particles currently in the box
    integer :: nptl_old             !< Number of particles without receivers
    integer :: nptl_max             !< Maximum number of particles allowed
    integer :: nptl_split           !< Number of particles from splitting
    integer :: nptl_inject          !< Number of injected particles
    integer :: tag_max              !< Maximum particle tag
    real(dp) :: leak                !< Leaking particles from boundary considering weight
    real(dp) :: leak_negp           !< Leaking particles with negative momentum

    real(dp) :: kpara0              !< kpara for particles with momentum p0
    real(dp) :: kret                !< The ratio of kpara to kperp
    integer :: momentum_dependency  !< kappa dependency on particle momentum
    integer :: mag_dependency       !< kappa dependency on magnetic field
    integer :: acc_region_flag      !< flag for whether to turn on acceleration in certain region
    real(dp) :: pindex              !< power index for the momentum dependency
    real(dp) :: p0    !< the standard deviation of the Gaussian distribution of momentum
    real(dp) :: b0    !< Initial magnetic field strength
    type kappa_type
        real(dp) :: knorm0                    !< normalization for spatial diffusion coefficient
        real(dp) :: kpara, kperp              !< Parallel and perpendicular kappa
        real(dp) :: skpara, skperp            !< Square root of 2*kpara and 2*kperp
        real(dp) :: skpara_perp               !< Square root of 2*(kpara - kperp)
        real(dp) :: kxx, kyy, kzz             !< Diagonal components of kappa tensor
        real(dp) :: kxy, kxz, kyz             !< Other components of kappa tensor
        real(dp) :: dkxx_dx, dkyy_dy, dkzz_dz !< Gradients of diagonal components
        real(dp) :: dkxy_dx, dkxy_dy          !< Gradients of other components
        real(dp) :: dkxz_dx, dkxz_dz
        real(dp) :: dkyz_dy, dkyz_dz
    end type kappa_type !< 21 variables

    real(dp), dimension(6) :: acc_region  !< from 0 to 1 (xmin, xmax, ymin, ymax, zmin, zmax)

    real(dp) :: dt_min      !< Minimum time step
    real(dp) :: dt_max      !< Maximum time step
    real(dp) :: dt_min_rel  !< Minimum time step w.r.t. one field time interval
    real(dp) :: dt_max_rel  !< Maximum time step w.r.t. one field time interval

    !< Momentum diffusion
    logical :: dpp_wave_flag, dpp_shear_flag
    !< Whether particle scattering is weak (tau*Omega >> 1)
    logical :: weak_scattering
    !< The scattering time for initial particles (kpara = v^2\tau/3)
    real(dp) :: tau0

    !< Pitch-angle diffusion
    real(dp) :: duu0

    !< Other flags and parameters
    logical :: deltab_flag, correlation_flag

    !< Particle drift
    real(dp) :: drift1 ! ev_ABL_0/p_0c
    real(dp) :: drift2 ! mev_ABL_0/p_0^2
    integer :: pcharge ! Particle charge in the unit of of e
    logical :: check_drift_2d ! Whether to check drift in 2D simulations

    !< Escaped particles
    type(particle_type), allocatable, dimension(:) :: escaped_ptls
    integer :: nptl_escaped, nptl_escaped_max
    !dir$ attributes align:128 :: escaped_ptls

    !< Whether to include transport along the 3rd dimension in 2D simulations
    logical :: include_3rd_dim_in2d_flag

    !< Whether the acceleration region is separated by a time-varying surface
    logical :: acc_by_surface_flag

    !< Particle tracking
    integer, allocatable, dimension(:) :: nsteps_tracked_ptls, noffsets_tracked_ptls
    integer, allocatable, dimension(:) :: tags_selected_ptls
    type(particle_type), allocatable, dimension(:) :: ptl_traj_points
    integer :: nsteps_tracked_tot   !< Total tracking steps for all tracked particles

    contains

    !---------------------------------------------------------------------------
    !< Initialize particle data
    !< Args:
    !<  nptl_max_allowed: the maximum number of particles allowed
    !---------------------------------------------------------------------------
    subroutine init_particles(nptl_max_allowed)
        implicit none
        integer, intent(in) :: nptl_max_allowed
        nptl_max = nptl_max_allowed
        allocate(ptls(nptl_max))
        ptls%x = 0.0
        ptls%y = 0.0
        ptls%z = 0.0
        ptls%p = 0.0
        ptls%v = 0.0
        ptls%mu = 0.0
        ptls%weight = 0.0
        ptls%t = 0.0
        ptls%dt = 0.0
        ptls%split_times = 0
        ptls%count_flag = COUNT_FLAG_OTHERS
        ptls%tag = 0
        ptls%nsteps_tracking = 0
        nptl_current = 0     ! No particle initially
        nptl_split = 0
        tag_max = 0

        !< Particles crossing domain boundaries
        allocate(senders(nptl_max / 10, ndim_field*2))
        allocate(recvers(nptl_max / 10, ndim_field*2))
        allocate(nsenders(ndim_field*2))
        allocate(nrecvers(ndim_field*2))
        senders%x = 0.0
        senders%y = 0.0
        senders%z = 0.0
        senders%p = 0.0
        senders%v = 0.0
        senders%mu = 0.0
        senders%weight = 0.0
        senders%t = 0.0
        senders%dt = 0.0
        senders%split_times = 0
        senders%count_flag = COUNT_FLAG_OTHERS
        senders%tag = 0
        senders%nsteps_tracking = 0

        recvers%x = 0.0
        recvers%y = 0.0
        recvers%z = 0.0
        recvers%p = 0.0
        recvers%v = 0.0
        recvers%mu = 0.0
        recvers%weight = 0.0
        recvers%t = 0.0
        recvers%dt = 0.0
        recvers%split_times = 0
        recvers%count_flag = COUNT_FLAG_OTHERS
        recvers%tag = 0
        recvers%nsteps_tracking = 0

        nsenders = 0
        nrecvers = 0

        !< Leaked particles
        leak = 0.0_dp
        leak_negp = 0.0_dp
    end subroutine init_particles

    !---------------------------------------------------------------------------
    !< Free particle data
    !---------------------------------------------------------------------------
    subroutine free_particles
        implicit none
        deallocate(ptls, senders, recvers)
        deallocate(nsenders, nrecvers)
    end subroutine free_particles

    !---------------------------------------------------------------------------
    !< Set parameters for momentum diffusion
    !< Args:
    !<  dpp_wave_int(integer): flag for momentum diffusion due to wave scattering
    !<  dpp_shear_int(integer): flag for momentum diffusion due to flow shear
    !<  weak_scattering_flag(integer): flag for weak-scattering regime
    !<  tau0_scattering: the scattering time for initial particles
    !---------------------------------------------------------------------------
    subroutine set_dpp_params(dpp_wave_int, dpp_shear_int, &
            weak_scattering_flag, tau0_scattering)
        implicit none
        integer, intent(in) :: dpp_wave_int, dpp_shear_int, weak_scattering_flag
        real(dp), intent(in) :: tau0_scattering
        dpp_wave_flag = .false.
        dpp_shear_flag = .false.
        weak_scattering = .false.
        if (dpp_wave_int == 1) dpp_wave_flag = .true.
        if (dpp_shear_int == 1) dpp_shear_flag = .true.
        if (weak_scattering_flag == 1) weak_scattering = .true.
        tau0 = tau0_scattering
    end subroutine set_dpp_params

    !---------------------------------------------------------------------------
    !< Set parameters for pitch-angle diffusion
    !< Args:
    !<  duu_init: duu for particles with p0
    !---------------------------------------------------------------------------
    subroutine set_duu_params(duu_init)
        implicit none
        real(dp), intent(in) :: duu_init
        duu0 = duu_init
    end subroutine set_duu_params

    !---------------------------------------------------------------------------
    !< Set other flags and parameters
    !< Args:
    !<  deltab_flag_int(integer): flag for magnetic fluctuation
    !<  correlation_flag_int(integer): flag for turbulence correlation length
    !<  include_3rd_dim(integer): flag for whether to include 3rd-dim transport
    !<  acc_by_surface(integer): flag for whether to separate the acceleration
    !<      region by a time-varying surface
    !---------------------------------------------------------------------------
    subroutine set_flags_params(deltab_flag_int, correlation_flag_int, &
            include_3rd_dim, acc_by_surface)
        implicit none
        integer, intent(in) :: deltab_flag_int, correlation_flag_int, &
            include_3rd_dim, acc_by_surface
        deltab_flag = .false.
        correlation_flag = .false.
        include_3rd_dim_in2d_flag = .false.
        acc_by_surface_flag = .false.
        if (deltab_flag_int == 1) deltab_flag = .true.
        if (correlation_flag_int == 1) correlation_flag = .true.
        if (include_3rd_dim == 1) include_3rd_dim_in2d_flag = .true.
        if (acc_by_surface == 1) acc_by_surface_flag = .true.
    end subroutine set_flags_params

    !---------------------------------------------------------------------------
    !< Set drift parameters for particle drift in 3D
    !< Args:
    !<  drift_param1: ev_ABL_0/pc
    !<  drift_param2: mev_ABL_0/p^2
    !<  charge: particle charge in the unit of e
    !---------------------------------------------------------------------------
    subroutine set_drift_parameters(drift_param1, drift_param2, charge)
        implicit none
        real(dp), intent(in) :: drift_param1, drift_param2
        integer, intent(in) :: charge
        drift1 = drift_param1
        drift2 = drift_param2
        pcharge = charge
    end subroutine set_drift_parameters

    !---------------------------------------------------------------------------
    !< Set the flag for checking particle drift in 2D simulations
    !< Args:
    !<  check_drift_2d_flag(integer): flag for checking drift in 2D
    !---------------------------------------------------------------------------
    subroutine set_flag_check_drift_2d(check_drift_2d_flag)
        implicit none
        integer, intent(in) :: check_drift_2d_flag
        check_drift_2d = .false.
        if (check_drift_2d_flag == 1) check_drift_2d = .true.
    end subroutine set_flag_check_drift_2d

    !---------------------------------------------------------------------------
    !< Set MPI datatype for particle type
    !---------------------------------------------------------------------------
    subroutine set_particle_datatype_mpi
        implicit none
        integer :: oldtypes(0:1), blockcounts(0:1)
        integer :: offsets(0:1), extent
        ! Setup description of the 9 MPI_DOUBLE fields.
        offsets(0) = 0
        oldtypes(0) = MPI_DOUBLE_PRECISION
        blockcounts(0) = 9
        ! Setup description of the 4 MPI_INTEGER fields.
        call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent, ierr)
        offsets(1) = blockcounts(0) * extent
        oldtypes(1) = MPI_INTEGER
        blockcounts(1) = 4
        ! Define structured type and commit it.
        call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, &
            particle_datatype_mpi, ierr)
        call MPI_TYPE_COMMIT(particle_datatype_mpi, ierr)
    end subroutine set_particle_datatype_mpi

    !---------------------------------------------------------------------------
    !< Free MPI datatype for particle type
    !---------------------------------------------------------------------------
    subroutine free_particle_datatype_mpi
        implicit none
        call MPI_TYPE_FREE(particle_datatype_mpi, ierr)
    end subroutine free_particle_datatype_mpi

    !---------------------------------------------------------------------------
    !< Inject one particle
    !< Args:
    !<  xpos, ypos, zpos: particle position
    !<  nptl_current: current number of particles
    !<  dist_flag: 0 for Maxwellian. 1 for delta function. 2 for power-law
    !<  particle_v0: particle velocity at p0 in the normalized velocity
    !<  mu: the cosine of the particle pitch-angle
    !<  ct_mhd: MHD simulation time frame, starting from 1
    !<  dt: the time interval
    !<  power_index: power-law index if dist_flag==2
    !---------------------------------------------------------------------------
    subroutine inject_one_particle(xpos, ypos, zpos, nptl_current, dist_flag, &
            particle_v0, mu, ct_mhd, dt, power_index)
        use mhd_config_module, only: tstamps_mhd
        use random_number_generator, only: unif_01, two_normals
        implicit none
        real(dp), intent(in) :: xpos, ypos, zpos, dt, power_index
        real(dp), intent(in) :: particle_v0, mu
        integer, intent(in) :: nptl_current, dist_flag, ct_mhd
        real(dp) :: r01, norm, fxp, ptmp, ftest
        ptls(nptl_current)%x = xpos
        ptls(nptl_current)%y = ypos
        ptls(nptl_current)%z = zpos
        if (dist_flag == 0) then
            ftest = 1.0
            fxp = 0.5
            do while (ftest > fxp)
                ptmp = (unif_01(0) * (pmax - pmin) + pmin) / p0  ! Need to normalized
                fxp = ptmp**2 * exp(-0.5*ptmp**2)
                ftest = unif_01(0) * 0.75  ! The maximum value is about 0.736
            enddo
            ptls(nptl_current)%p = ptmp * p0
        else if (dist_flag == 1) then
            ptls(nptl_current)%p = p0
        else if (dist_flag == 2) then
            r01 = unif_01(0)
            if (int(power_index) == 1) then
                ptls(nptl_current)%p = (pmax / p0)**r01 * p0
            else
                norm = pmax**(-power_index + 1) - p0**(-power_index + 1)
                ptls(nptl_current)%p = &
                    (r01 * norm + p0**(-power_index + 1))**(1.0 / (-power_index + 1))
            endif
        endif
        ptls(nptl_current)%v = particle_v0 * ptls(nptl_current)%p / p0
        ptls(nptl_current)%mu = mu
        ptls(nptl_current)%weight = 1.0d0
        ptls(nptl_current)%t = tstamps_mhd(ct_mhd)
        ptls(nptl_current)%dt = dt
        ptls(nptl_current)%split_times = 0
        ptls(nptl_current)%count_flag = COUNT_FLAG_INBOX
        tag_max = tag_max + 1
        ptls(nptl_current)%tag = tag_max
        ptls(nptl_current)%nsteps_tracking = 0
    end subroutine inject_one_particle

    !---------------------------------------------------------------------------
    !< Inject particles which are spatially uniform
    !< Args:
    !<  nptl: number of particles to be injected
    !<  dt: the time interval
    !<  dist_flag: 0 for Maxwellian. 1 for delta function. 2 for power-law
    !<  particle_v0: particle velocity at p0 in the normalized velocity
    !<  ct_mhd: MHD simulation time frame, starting from 1
    !<  part_box: box to inject particles
    !<  power_index: power-law index if dist_flag==2
    !---------------------------------------------------------------------------
    subroutine inject_particles_spatial_uniform(nptl, dt, dist_flag, &
            particle_v0, ct_mhd, part_box, power_index)
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: get_ncells_large_jz
        use random_number_generator, only: unif_01
        implicit none
        integer, intent(in) :: nptl, dist_flag, ct_mhd
        real(dp), intent(in) :: dt, power_index, particle_v0
        real(dp), intent(in), dimension(6) :: part_box
        integer :: i, imod2
        real(dp) :: xmin, ymin, zmin, xmax, ymax, zmax
        real(dp) :: xmin_box, ymin_box, zmin_box
        real(dp) :: xmax_box, ymax_box, zmax_box
        real(dp) :: jz_min, xtmp, ytmp, ztmp, mu_tmp
        integer :: ncells_large_jz, ncells_large_jz_g
        logical :: inbox

        xmin_box = part_box(1)
        ymin_box = part_box(2)
        zmin_box = part_box(3)
        xmax_box = part_box(4)
        ymax_box = part_box(5)
        zmax_box = part_box(6)
        xmin = fconfig%xmin
        ymin = fconfig%ymin
        zmin = fconfig%zmin
        xmax = fconfig%xmax
        ymax = fconfig%ymax
        zmax = fconfig%zmax

        jz_min = -1.0 ! set negative so all cells with counted if the part_box
        ncells_large_jz = get_ncells_large_jz(jz_min, spherical_coord_flag, part_box)
        call MPI_ALLREDUCE(ncells_large_jz, ncells_large_jz_g, 1, &
            MPI_INTEGER, MPI_SUM, mpi_sub_comm, ierr)
        call MPI_BARRIER(mpi_sub_comm, ierr)
        nptl_inject = int(nptl * mpi_sub_size * dble(ncells_large_jz) / dble(ncells_large_jz_g))

        do i = 1, nptl_inject
            nptl_current = nptl_current + 1
            if (nptl_current > nptl_max) nptl_current = nptl_max
            inbox = .false.
            do while (.not. inbox)
                xtmp = unif_01(0) * (xmax - xmin) + xmin
                ytmp = unif_01(0) * (ymax - ymin) + ymin
                ztmp = unif_01(0) * (zmax - zmin) + zmin
                inbox = xtmp > xmin_box .and. xtmp < xmax_box .and. &
                        ytmp > ymin_box .and. ytmp < ymax_box .and. &
                        ztmp > zmin_box .and. ztmp < zmax_box
            enddo
            mu_tmp = 2.0d0 * unif_01(0) - 1.0d0
            call inject_one_particle(xtmp, ytmp, ztmp, nptl_current, &
                dist_flag, particle_v0, mu_tmp, ct_mhd, dt, power_index)
        enddo

        if (mpi_rank == master) then
            write(*, "(A)") "Finished injecting particles"
        endif
    end subroutine inject_particles_spatial_uniform

    !---------------------------------------------------------------------------
    !< Inject particles at shock
    !< Args:
    !<  nptl: number of particles to be injected
    !<  dt: the time interval
    !<  dist_flag: 0 for Maxwellian. 1 for delta function. 2 for power-law
    !<  particle_v0: particle velocity at p0 in the normalized velocity
    !<  ct_mhd: MHD simulation time frame, starting from 1
    !<  power_index: power-law index if dist_flag==2
    !---------------------------------------------------------------------------
    subroutine inject_particles_at_shock(nptl, dt, dist_flag, particle_v0, &
            ct_mhd, power_index)
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config, tstamps_mhd
        use mhd_data_parallel, only: interp_shock_location
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: nptl, dist_flag, ct_mhd
        real(dp), intent(in) :: dt, power_index, particle_v0
        integer :: i, iy, iz
        real(dp) :: xmin, ymin, xmax, ymax, zmin, zmax
        real(dp) :: ry, rz, dpy, dpz, shock_xpos
        real(dp) :: r01, norm, fxp, ptmp, ftest
        integer, dimension(2) :: pos
        real(dp), dimension(4) :: weights

        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax
        zmin = fconfig%zmin
        zmax = fconfig%zmax

        nptl_inject = nptl
        do i = 1, nptl
            nptl_current = nptl_current + 1
            if (nptl_current > nptl_max) nptl_current = nptl_max
            ptls(nptl_current)%y = unif_01(0)*(ymax-ymin) + ymin
            dpy = ptls(nptl_current)%y / mhd_config%dy
            iy = floor(dpy)
            ptls(nptl_current)%z = unif_01(0)*(zmax-zmin) + zmin
            dpz = ptls(nptl_current)%z / mhd_config%dz
            iz = floor(dpz)
            !< We assume that particles are inject at the earlier shock location
            pos = (/ iy, iz /)
            ry = dpy - iy
            rz = dpy - iy
            weights(1) = (1 - ry) * (1 - rz)
            weights(2) = ry * (1 - rz)
            weights(3) = (1 - ry) * rz
            weights(4) = ry * rz
            shock_xpos = interp_shock_location(pos, weights, 0.0d0) + 2  ! Two ghost cells
            ptls(nptl_current)%x = shock_xpos * (xmax - xmin) / fconfig%nxg
            if (dist_flag == 0) then
                ftest = 1.0
                fxp = 0.5
                do while (ftest > fxp)
                    ptmp = (unif_01(0) * (pmax - pmin) + pmin) / p0  ! Need to normalized
                    fxp = ptmp**2 * exp(-0.5*ptmp**2)
                    ftest = unif_01(0) * 0.75  ! The maximum value is about 0.736
                enddo
                ptls(nptl_current)%p = ptmp * p0
            else if (dist_flag == 1) then
                ptls(nptl_current)%p = p0
            else if (dist_flag == 2) then
                r01 = unif_01(0)
                if (int(power_index) == 1) then
                    ptls(nptl_current)%p = (pmax / p0)**r01 * p0
                else
                    norm = pmax**(-power_index + 1) - p0**(-power_index + 1)
                    ptls(nptl_current)%p = &
                        (r01 * norm + p0**(-power_index + 1))**(1.0 / (-power_index + 1))
                endif
            endif
            ptls(nptl_current)%v = particle_v0 * ptls(nptl_current)%p / p0
            ptls(nptl_current)%mu = 2.0d0 * unif_01(0) - 1.0d0
            ptls(nptl_current)%weight = 1.0d0
            ptls(nptl_current)%t = tstamps_mhd(ct_mhd)
            ptls(nptl_current)%dt = dt
            ptls(nptl_current)%split_times = 0
            ptls(nptl_current)%count_flag = COUNT_FLAG_INBOX
            tag_max = tag_max + 1
            ptls(nptl_current)%tag = tag_max
            ptls(nptl_current)%nsteps_tracking = 0
        enddo
        if (mpi_rank == master) then
            write(*, "(A)") "Finished injecting particles at the shock"
        endif
    end subroutine inject_particles_at_shock

    !---------------------------------------------------------------------------
    !< Gets the corner indices and weights for linear interpolation
    !< Args:
    !<  px, py, pz: corner position in real values
    !<  pos: indices of the cells where the particle is
    !<  weights: for linear interpolation
    !---------------------------------------------------------------------------
    subroutine get_interp_paramters(px, py, pz, pos, weights)
        implicit none
        real(dp), intent(in) :: px, py, pz
        real(dp), dimension(8), intent(out) :: weights
        integer, dimension(3), intent(out) :: pos
        real(dp) :: rx, ry, rz, rx1, ry1, rz1

        if (ndim_field == 1) then
            pos = (/ floor(px)+1,  1, 1 /)
            ry = 0.0
            rz = 0.0
        else if (ndim_field == 2) then
            pos = (/ floor(px)+1,  floor(py)+1, 1 /)
            ry = py - pos(2) + 1
            rz = 0.0
        else
            pos = (/ floor(px)+1,  floor(py)+1, floor(pz)+1 /)
            ry = py - pos(2) + 1
            rz = pz - pos(3) + 1
        endif
        rx = px - pos(1) + 1
        rx1 = 1.0 - rx
        ry1 = 1.0 - ry
        rz1 = 1.0 - rz
        weights(1) = rx1 * ry1 * rz1
        weights(2) = rx * ry1 * rz1
        weights(3) = rx1 * ry * rz1
        weights(4) = rx * ry * rz1
        weights(5) = rx1 * ry1 * rz
        weights(6) = rx * ry1 * rz
        weights(7) = rx1 * ry * rz
        weights(8) = rx * ry * rz
    end subroutine get_interp_paramters

    !---------------------------------------------------------------------------
    !< Binary search. We are searching for the cell not the exact grid point.
    !< Modified from https://rosettacode.org/wiki/Binary_search#Fortran
    !---------------------------------------------------------------------------
    recursive function binarySearch_R (a, value) result (bsresult)
        real(dp), intent(in) :: a(:), value
        integer :: bsresult, mid

        mid = size(a)/2 + 1
        if (size(a) == 0) then
            bsresult = 0
        else if (a(mid) > value) then
            bsresult = binarySearch_R(a(:mid-1), value)
        else if (a(mid) < value) then
            bsresult = binarySearch_R(a(mid+1:), value)
            bsresult = mid + bsresult
        else
            bsresult = mid
        end if
    end function binarySearch_R

    !---------------------------------------------------------------------------
    !< Gets the corner indices and weights for linear interpolation in spherical
    !< coordinates
    !< Args:
    !<  x, y, z: corner position in real values
    !<  pos: indices of the cells where the particle is
    !<  weights: for interpolation
    !---------------------------------------------------------------------------
    subroutine get_interp_paramters_spherical(x, y, z, pos, weights)
        use constants, only: pi
        use mhd_data_parallel, only: xpos_local, ypos_local, zpos_local
        implicit none
        real(dp), intent(in) :: x, y, z  ! r, theta, phi
        real(dp), dimension(8), intent(out) :: weights
        integer, dimension(3), intent(out) :: pos
        real(dp) :: rx, ry, rz, rx1, ry1, rz1, dv
        real(dp) :: ctheta1, ctheta2, ctheta
        real(dp) :: x3_1, x3_2, x3, one_third

        ! The lbound of the position array is -1, so -2
        pos(1) = binarySearch_R(xpos_local, x) - 2
        if (ndim_field > 1) then
            pos(2) = binarySearch_R(ypos_local, y) - 2
            if (ndim_field == 3) then
                pos(3) = binarySearch_R(zpos_local, z) - 2
            else
                pos(3) = 1
            endif
        else
            pos(2) = 1
        endif
        ! We assume here theta is from 0 to pi
        ctheta1 = cos(ypos_local(pos(2)))
        ctheta2 = cos(ypos_local(pos(2)+1))
        ctheta = cos(y)
        one_third = 1.0_dp / 3.0_dp
        x3 = x**3 * one_third
        x3_1 = xpos_local(pos(1))**3 * one_third
        x3_2 = xpos_local(pos(1)+1)**3 * one_third

        rx = x3 - x3_1
        rx1 = x3_2 - x3
        if (ndim_field > 1) then
            ry = ctheta1 - ctheta
            ry1 = ctheta - ctheta2
            if (ndim_field == 3) then
                rz = z - zpos_local(pos(3))
                rz1 = zpos_local(pos(3)+1) - z
            else
                rz = 0.0
                rz1 = 1.0
            endif
        else
            ry = 0.0
            ry1 = 1.0
        endif

        ! Volume of the cell
        dv = (x3_2 - x3_1) * (ry + ry1) * (rz + rz1)

        weights(1) = rx1 * ry1 * rz1
        weights(2) = rx * ry1 * rz1
        weights(3) = rx1 * ry * rz1
        weights(4) = rx * ry * rz1
        weights(5) = rx1 * ry1 * rz
        weights(6) = rx * ry1 * rz
        weights(7) = rx1 * ry * rz
        weights(8) = rx * ry * rz

        weights = weights / dv
    end subroutine get_interp_paramters_spherical

    !---------------------------------------------------------------------------
    !< Inject particles where current density jz is large
    !< Note that this might not work if each MPI rank only handle part of the
    !< MHD simulation, because jz is always close to 0 in some part of the MHD
    !< simulations. Therefore, be smart on how to participate the simulation.
    !< Args:
    !<  nptl: number of particles to be injected
    !<  dt: the time interval
    !<  dist_flag: 0 for Maxwellian. 1 for delta function. 2 for power-law
    !<  particle_v0: particle velocity at p0 in the normalized velocity
    !<  ct_mhd: MHD simulation time frame, starting from 1
    !<  inject_same_nptl: whether to inject the same of particles every step
    !<  jz_min: the minimum jz
    !<  ncells_large_jz_norm ! Normalization for the number of cells with large jz
    !<  part_box: box to inject particles
    !<  power_index: power-law index if dist_flag==2
    !---------------------------------------------------------------------------
    subroutine inject_particles_at_large_jz(nptl, dt, dist_flag, particle_v0, &
            ct_mhd, inject_same_nptl, jz_min, ncells_large_jz_norm, part_box, &
            power_index)
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: interp_fields
        use mhd_data_parallel, only: get_ncells_large_jz
        use random_number_generator, only: unif_01
        implicit none
        integer, intent(in) :: nptl, dist_flag, ct_mhd
        integer, intent(in) :: inject_same_nptl, ncells_large_jz_norm
        real(dp), intent(in) :: dt, jz_min, power_index, particle_v0
        real(dp), intent(in), dimension(6) :: part_box
        real(dp) :: xmin, ymin, zmin, xmax, ymax, zmax
        real(dp) :: xmin_box, ymin_box, zmin_box
        real(dp) :: xmax_box, ymax_box, zmax_box
        real(dp) :: xtmp, ytmp, ztmp, px, py, pz
        real(dp) :: dxm, dym, dzm, dby_dx, dbx_dy
        real(dp) :: rt, jz, by, mu_tmp
        real(dp), dimension(2) :: rands
        real(dp) :: r01, norm
        real(dp), dimension(nfields+ngrads) :: fields !< Fields at particle position
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights
        integer :: i, imod2
        integer :: ncells_large_jz, ncells_large_jz_g
        !dir$ attributes align:256 :: fields

        xmin_box = part_box(1)
        ymin_box = part_box(2)
        zmin_box = part_box(3)
        xmax_box = part_box(4)
        ymax_box = part_box(5)
        zmax_box = part_box(6)
        xmin = fconfig%xmin
        ymin = fconfig%ymin
        zmin = fconfig%zmin
        xmax = fconfig%xmax
        ymax = fconfig%ymax
        zmax = fconfig%zmax
        dxm = mhd_config%dx
        dym = mhd_config%dy
        dzm = mhd_config%dz

        xtmp = xmin_box
        ytmp = ymin_box
        ztmp = zmin_box
        px = 0.0_dp
        py = 0.0_dp
        pz = 0.0_dp

        !< When each MPI only reads part of the MHD data, the number of cells
        !< with larger jz in the local domain will be different from mpi_rank
        !< to mpi_rank. That's why we need to redistribute the number of particles
        !< to inject.
        ncells_large_jz = get_ncells_large_jz(jz_min, spherical_coord_flag, part_box)
        if (inject_same_nptl == 1) then
            call MPI_ALLREDUCE(ncells_large_jz, ncells_large_jz_g, 1, &
                MPI_INTEGER, MPI_SUM, mpi_sub_comm, ierr)
            nptl_inject = int(nptl * mpi_sub_size * &
                (dble(ncells_large_jz) / dble(ncells_large_jz_g)))
        else
            nptl_inject = int(nptl * mpi_sub_size * &
                (dble(ncells_large_jz) / dble(ncells_large_jz_norm)))
        endif

        do i = 1, nptl_inject
            nptl_current = nptl_current + 1
            if (nptl_current > nptl_max) nptl_current = nptl_max
            jz = -2.0_dp
            do while (jz < jz_min)
                xtmp = unif_01(0) * (xmax - xmin) + xmin
                ytmp = unif_01(0) * (ymax - ymin) + ymin
                ztmp = unif_01(0) * (zmax - zmin) + zmin
                if (xtmp >= xmin_box .and. xtmp <= xmax_box .and. &
                    ytmp >= ymin_box .and. ytmp <= ymax_box .and. &
                    ztmp >= zmin_box .and. ztmp <= zmax_box) then
                    px = (xtmp - xmin) / dxm
                    py = (ytmp - ymin) / dym
                    pz = (ztmp - zmin) / dzm
                    rt = 0.0_dp
                    if (spherical_coord_flag) then
                        call get_interp_paramters_spherical(xtmp, ytmp, ztmp, pos, weights)
                    else
                        call get_interp_paramters(px, py, pz, pos, weights)
                    endif
                    call interp_fields(pos, weights, rt, fields)
                    dbx_dy = fields(nfields+14)
                    dby_dx = fields(nfields+16)
                    if (spherical_coord_flag) then
                        by = fields(6)
                        jz = abs(dby_dx + (by - dbx_dy) / xtmp)
                    else
                        jz = abs(dby_dx - dbx_dy)
                    endif
                else ! not in part_box
                    jz = -3.0_dp
                endif
            enddo
            mu_tmp = 2.0d0 * unif_01(0) - 1.0d0
            call inject_one_particle(xtmp, ytmp, ztmp, nptl_current, &
                dist_flag, particle_v0, mu_tmp, ct_mhd, dt, power_index)
        enddo
        if (mpi_rank == master) then
            write(*, "(A)") "Finished injecting particles where jz is large"
        endif
    end subroutine inject_particles_at_large_jz

    !---------------------------------------------------------------------------
    !< Inject particles where the turbulence variance is large
    !< Note that this might not work if each MPI rank only handle part of the
    !< MHD simulation, because jz is always close to 0 in some part of the MHD
    !< simulations. Therefore, be smart on how to participate the simulation.
    !< Args:
    !<  nptl: number of particles to be injected
    !<  dt: the time interval
    !<  dist_flag: 0 for Maxwellian. 1 for delta function. 2 for power-law
    !<  particle_v0: particle velocity at p0 in the normalized velocity
    !<  ct_mhd: MHD simulation time frame, starting from 1
    !<  inject_same_nptl: whether to inject the same of particles every step
    !<  db2_min: the minimum db2
    !<  ncells_large_db2_norm ! Normalization for the number of cells with large db2
    !<  part_box: box to inject particles
    !<  power_index: power-law index if dist_flag==2
    !---------------------------------------------------------------------------
    subroutine inject_particles_at_large_db2(nptl, dt, dist_flag, particle_v0, &
            ct_mhd, inject_same_nptl, db2_min, ncells_large_db2_norm, part_box, &
            power_index)
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: interp_magnetic_fluctuation
        use mhd_data_parallel, only: get_ncells_large_db2
        use random_number_generator, only: unif_01
        implicit none
        integer, intent(in) :: nptl, dist_flag, ct_mhd
        integer, intent(in) :: inject_same_nptl, ncells_large_db2_norm
        real(dp), intent(in) :: dt, db2_min, power_index, particle_v0
        real(dp), intent(in), dimension(6) :: part_box
        real(dp) :: xmin, ymin, zmin, xmax, ymax, zmax
        real(dp) :: xmin_box, ymin_box, zmin_box
        real(dp) :: xmax_box, ymax_box, zmax_box
        real(dp) :: xtmp, ytmp, ztmp, px, py, pz
        real(dp) :: dxm, dym, dzm
        real(dp) :: rt, db2, by, mu_tmp
        real(dp), dimension(2) :: rands
        real(dp) :: r01, norm
        real(dp), dimension(4) :: db2_array
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights
        integer :: i, imod2
        integer :: ncells_large_db2, ncells_large_db2_g
        !dir$ attributes align:32 :: db2_array

        xmin_box = part_box(1)
        ymin_box = part_box(2)
        zmin_box = part_box(3)
        xmax_box = part_box(4)
        ymax_box = part_box(5)
        zmax_box = part_box(6)
        xmin = fconfig%xmin
        ymin = fconfig%ymin
        zmin = fconfig%zmin
        xmax = fconfig%xmax
        ymax = fconfig%ymax
        zmax = fconfig%zmax
        dxm = mhd_config%dx
        dym = mhd_config%dy
        dzm = mhd_config%dz

        xtmp = xmin_box
        ytmp = ymin_box
        ztmp = zmin_box
        px = 0.0_dp
        py = 0.0_dp
        pz = 0.0_dp

        !< When each MPI only reads part of the MHD data, the number of cells
        !< with larger db2 in the local domain will be different from mpi_rank
        !< to mpi_rank. That's why we need to redistribute the number of particles
        !< to inject.
        ncells_large_db2 = get_ncells_large_db2(db2_min, spherical_coord_flag, part_box)
        if (inject_same_nptl == 1) then
            call MPI_ALLREDUCE(ncells_large_db2, ncells_large_db2_g, 1, &
                MPI_INTEGER, MPI_SUM, mpi_sub_comm, ierr)
            nptl_inject = int(nptl * mpi_sub_size * &
                (dble(ncells_large_db2) / dble(ncells_large_db2_g)))
        else
            nptl_inject = int(nptl * mpi_sub_size * &
                (dble(ncells_large_db2) / dble(ncells_large_db2_norm)))
        endif

        do i = 1, nptl_inject
            nptl_current = nptl_current + 1
            if (nptl_current > nptl_max) nptl_current = nptl_max
            db2 = -2.0_dp
            do while (db2 < db2_min)
                xtmp = unif_01(0) * (xmax - xmin) + xmin
                ytmp = unif_01(0) * (ymax - ymin) + ymin
                ztmp = unif_01(0) * (zmax - zmin) + zmin
                if (xtmp >= xmin_box .and. xtmp <= xmax_box .and. &
                    ytmp >= ymin_box .and. ytmp <= ymax_box .and. &
                    ztmp >= zmin_box .and. ztmp <= zmax_box) then
                    px = (xtmp - xmin) / dxm
                    py = (ytmp - ymin) / dym
                    pz = (ztmp - zmin) / dzm
                    rt = 0.0_dp
                    if (spherical_coord_flag) then
                        call get_interp_paramters_spherical(xtmp, ytmp, ztmp, pos, weights)
                    else
                        call get_interp_paramters(px, py, pz, pos, weights)
                    endif
                    call interp_magnetic_fluctuation(pos, weights, rt, db2_array)
                    db2 = db2_array(1)
                else ! not in part_box
                    db2 = -3.0_dp
                endif
            enddo
            mu_tmp = 2.0d0 * unif_01(0) - 1.0d0
            call inject_one_particle(xtmp, ytmp, ztmp, nptl_current, &
                dist_flag, particle_v0, mu_tmp, ct_mhd, dt, power_index)
        enddo
        if (mpi_rank == master) then
            write(*, "(A)") "Finished injecting particles where db2 is large"
        endif
    end subroutine inject_particles_at_large_db2

    !---------------------------------------------------------------------------
    !< Inject particles where the turbulence variance is large
    !< Note that this might not work if each MPI rank only handle part of the
    !< MHD simulation, because jz is always close to 0 in some part of the MHD
    !< simulations. Therefore, be smart on how to participate the simulation.
    !< Args:
    !<  nptl: number of particles to be injected
    !<  dt: the time interval
    !<  dist_flag: 0 for Maxwellian. 1 for delta function. 2 for power-law
    !<  particle_v0: particle velocity at p0 in the normalized velocity
    !<  ct_mhd: MHD simulation time frame, starting from 1
    !<  inject_same_nptl: whether to inject the same of particles every step
    !<  divv_min: the minimum divv
    !<  ncells_large_divv_norm ! Normalization for the number of cells with large divv
    !<  part_box: box to inject particles
    !<  power_index: power-law index if dist_flag==2
    !---------------------------------------------------------------------------
    subroutine inject_particles_at_large_divv(nptl, dt, dist_flag, particle_v0, &
            ct_mhd, inject_same_nptl, divv_min, ncells_large_divv_norm, part_box, &
            power_index)
        use constants, only: pi
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: interp_fields
        use mhd_data_parallel, only: get_ncells_large_divv
        use random_number_generator, only: unif_01
        implicit none
        integer, intent(in) :: nptl, dist_flag, ct_mhd
        integer, intent(in) :: inject_same_nptl, ncells_large_divv_norm
        real(dp), intent(in) :: dt, divv_min, power_index, particle_v0
        real(dp), intent(in), dimension(6) :: part_box
        real(dp) :: xmin, ymin, zmin, xmax, ymax, zmax
        real(dp) :: xmin_box, ymin_box, zmin_box
        real(dp) :: xmax_box, ymax_box, zmax_box
        real(dp) :: xtmp, ytmp, ztmp, px, py, pz
        real(dp) :: dxm, dym, dzm
        real(dp) :: rt, divv, by, mu_tmp
        real(dp), dimension(2) :: rands
        real(dp) :: r01, norm
        real(dp), dimension(nfields+ngrads) :: fields !< Fields at particle position
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights
        real(dp) :: ctheta, istheta
        integer :: i, imod2
        integer :: ncells_large_divv, ncells_large_divv_g
        !dir$ attributes align:256 :: fields

        xmin_box = part_box(1)
        ymin_box = part_box(2)
        zmin_box = part_box(3)
        xmax_box = part_box(4)
        ymax_box = part_box(5)
        zmax_box = part_box(6)
        xmin = fconfig%xmin
        ymin = fconfig%ymin
        zmin = fconfig%zmin
        xmax = fconfig%xmax
        ymax = fconfig%ymax
        zmax = fconfig%zmax
        dxm = mhd_config%dx
        dym = mhd_config%dy
        dzm = mhd_config%dz

        xtmp = xmin_box
        ytmp = ymin_box
        ztmp = zmin_box
        px = 0.0_dp
        py = 0.0_dp
        pz = 0.0_dp

        !< When each MPI only reads part of the MHD data, the number of cells
        !< with larger divv in the local domain will be different from mpi_rank
        !< to mpi_rank. That's why we need to redistribute the number of particles
        !< to inject.
        ncells_large_divv = get_ncells_large_divv(divv_min, spherical_coord_flag, part_box)
        if (inject_same_nptl == 1) then
            call MPI_ALLREDUCE(ncells_large_divv, ncells_large_divv_g, 1, &
                MPI_INTEGER, MPI_SUM, mpi_sub_comm, ierr)
            nptl_inject = int(nptl * mpi_sub_size * &
                (dble(ncells_large_divv) / dble(ncells_large_divv_g)))
        else
            nptl_inject = int(nptl * mpi_sub_size * &
                (dble(ncells_large_divv) / dble(ncells_large_divv_norm)))
        endif

        do i = 1, nptl_inject
            nptl_current = nptl_current + 1
            if (nptl_current > nptl_max) nptl_current = nptl_max
            divv = 2.0_dp
            do while (-divv < divv_min)
                xtmp = unif_01(0) * (xmax - xmin) + xmin
                ytmp = unif_01(0) * (ymax - ymin) + ymin
                ztmp = unif_01(0) * (zmax - zmin) + zmin
                if (xtmp >= xmin_box .and. xtmp <= xmax_box .and. &
                    ytmp >= ymin_box .and. ytmp <= ymax_box .and. &
                    ztmp >= zmin_box .and. ztmp <= zmax_box) then
                    px = (xtmp - xmin) / dxm
                    py = (ytmp - ymin) / dym
                    pz = (ztmp - zmin) / dzm
                    rt = 0.0_dp
                    if (spherical_coord_flag) then
                        call get_interp_paramters_spherical(xtmp, ytmp, ztmp, pos, weights)
                    else
                        call get_interp_paramters(px, py, pz, pos, weights)
                    endif
                    call interp_fields(pos, weights, rt, fields)
                    if (spherical_coord_flag) then
                        divv = fields(nfields+1) + 2.0 * fields(1) / xtmp
                        if (ndim_field > 1) then
                            ctheta = cos(ytmp)
                            istheta = 1.0 / sin(ytmp)
                            divv = divv + (fields(nfields+5) + &
                                           fields(2)*ctheta*istheta) / xtmp
                            if (ndim_field > 2) then
                                divv = divv + fields(3)*istheta / xtmp
                            endif
                        endif
                    else
                        divv = fields(nfields+1)
                        if (ndim_field > 1) then
                            divv = divv + fields(nfields+5)
                            if (ndim_field > 2) then
                                divv = divv + fields(nfields+9)
                            endif
                        endif
                    endif
                else ! not in part_box
                    divv = 3.0_dp
                endif
            enddo
            mu_tmp = 2.0d0 * unif_01(0) - 1.0d0
            call inject_one_particle(xtmp, ytmp, ztmp, nptl_current, &
                dist_flag, particle_v0, mu_tmp, ct_mhd, dt, power_index)
        enddo
        if (mpi_rank == master) then
            write(*, "(A)") "Finished injecting particles where divv is negatively large"
        endif
    end subroutine inject_particles_at_large_divv

    !---------------------------------------------------------------------------
    !< Particle mover in one cycle
    !< Args:
    !<  t0: the starting time
    !<  dtf: the time interval of the MHD fields
    !<  track_particle_flag: whether to track particles, 0 for no, 1 for yes
    !<  nptl_selected: number of selected particles
    !<  nsteps_interval: save particle points every nsteps_interval
    !<  num_fine_steps: number of fine time steps
    !<  focused_transport: whether to the Focused Transport equation
    !---------------------------------------------------------------------------
    subroutine particle_mover_one_cycle(t0, dtf, track_particle_flag, &
            nptl_selected, nsteps_interval, num_fine_steps, focused_transport)
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: interp_fields, interp_magnetic_fluctuation, &
            interp_correlation_length
        use acc_region_surface, only: interp_acc_surface
        implicit none
        real(dp), intent(in) :: t0, dtf
        integer, intent(in) :: track_particle_flag, nptl_selected
        integer, intent(in) :: nsteps_interval, num_fine_steps
        logical, intent(in) :: focused_transport
        real(dp) :: dxm, dym, dzm, xmin, xmax, ymin, ymax, zmin, zmax
        real(dp) :: xmin1, xmax1, ymin1, ymax1, zmin1, zmax1
        real(dp) :: deltax, deltay, deltaz, deltap, deltav, deltamu
        real(dp) :: dt_target, dt_fine
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights
        real(dp) :: px, py, pz, rt
        integer :: i, tracking_step, offset, step, thread_id
        type(particle_type) :: ptl
        type(kappa_type) :: kappa
        real(dp) :: surface_height1, surface_height2
        real(dp), dimension(nfields+ngrads) :: fields
        real(dp), dimension(4) :: db2, lc
        !dir$ attributes align:256 :: fields
        !dir$ attributes align:32 :: db2
        !dir$ attributes align:32 :: lc

        dt_fine = dtf / num_fine_steps
        dxm = mhd_config%dx
        dym = mhd_config%dy
        dzm = mhd_config%dz
        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax
        zmin = fconfig%zmin
        zmax = fconfig%zmax
        xmin1 = xmin - dxm * 0.5
        xmax1 = xmax + dxm * 0.5
        ymin1 = ymin - dym * 0.5
        ymax1 = ymax + dym * 0.5
        zmin1 = zmin - dzm * 0.5
        zmax1 = zmax + dzm * 0.5

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ptl, kappa, &
        !$OMP& fields, db2, lc, surface_height1, surface_height2, &
        !$OMP& deltax, deltay, deltaz, deltap, deltav, deltamu, &
        !$OMP& dt_target, pos, weights, px, py, pz, rt, &
        !$OMP& tracking_step, offset, step, thread_id)
        thread_id = 0
#if (defined USE_OPENMP)
        thread_id = OMP_GET_THREAD_NUM()
#endif
        !$OMP DO
        do i = nptl_old + 1, nptl_current
            ptl = ptls(i)
            deltax = 0.0
            deltay = 0.0
            deltaz = 0.0
            deltap = 0.0
            deltav = 0.0
            deltamu = 0.0

            ! The targeted time depends on the time stamp of the particle.
            step = ceiling((ptl%t - t0) / dt_fine)
            if (step <= 0) then
                dt_target = dt_fine
            else
                dt_target = step * dt_fine
            endif
            if (dt_target > dtf) then
                dt_target = dtf
            endif

            ! Safe check
            if (ptl%p < 0.0 .and. ptl%count_flag == COUNT_FLAG_INBOX) then
                ptl%count_flag = COUNT_FLAG_OTHERS
                !$OMP ATOMIC UPDATE
                leak_negp = leak_negp + ptl%weight
            else
                call particle_boundary_condition(ptl, xmin1, xmax1, &
                    ymin1, ymax1, zmin1, zmax1)
            endif
            if (ptl%count_flag /= COUNT_FLAG_INBOX) then
                ptls(i) = ptl
                cycle
            endif

            ! Loop over fine time steps between each MHD time interval
            ! The loop might be terminated early when dt_target is reached.
            do while (dt_target < (dtf + dt_fine*0.1)) ! 0.1 is for safe comparison
                if (ptl%count_flag /= COUNT_FLAG_INBOX) then
                    exit
                endif
                do while ((ptl%t - t0) < dt_target .and. ptl%count_flag == COUNT_FLAG_INBOX)
                    ! Check if particles are leaked or out of the local domain
                    if (ptl%p < 0.0) then
                        ptl%count_flag = COUNT_FLAG_OTHERS
                        !$OMP ATOMIC UPDATE
                        leak_negp = leak_negp + ptl%weight
                    else
                        call particle_boundary_condition(ptl, xmin1, xmax1, &
                            ymin1, ymax1, zmin1, zmax1)
                    endif
                    if (ptl%count_flag /= COUNT_FLAG_INBOX) then
                        exit
                    endif

                    ! Field interpolation parameters
                    px = (ptl%x-xmin) / dxm
                    py = (ptl%y-ymin) / dym
                    pz = (ptl%z-zmin) / dzm
                    rt = (ptl%t - t0) / dtf
                    if (spherical_coord_flag) then
                        call get_interp_paramters_spherical(ptl%x, ptl%y, ptl%z, pos, weights)
                    else
                        call get_interp_paramters(px, py, pz, pos, weights)
                    endif
                    call interp_fields(pos, weights, rt, fields)
                    if (deltab_flag) then
                        call interp_magnetic_fluctuation(pos, weights, rt, db2)
                    endif
                    if (correlation_flag) then
                        call interp_correlation_length(pos, weights, rt, lc)
                    endif
                    call calc_spatial_diffusion_coefficients(ptl, focused_transport, &
                        fields, db2, lc, kappa)
                    call set_time_step(t0, dt_target, fields, kappa, ptl)
                    if (focused_transport) then
                        if (ndim_field == 1) then
                            call push_particle_1d_ft(thread_id, rt, ptl, fields, db2, &
                                lc, kappa, deltax, deltap, deltav, deltamu)
                        else if (ndim_field == 2) then
                            if (include_3rd_dim_in2d_flag) then
                                call push_particle_2d_include_3rd_ft(thread_id, rt, ptl, &
                                    fields, db2, lc, kappa, deltax, deltay, deltaz, &
                                    deltap, deltav, deltamu)
                            else
                                call push_particle_2d_ft(thread_id, rt, ptl, fields, db2, &
                                    lc, kappa, deltax, deltay, deltap, deltav, deltamu)
                            endif
                        else
                            if (acc_by_surface_flag) then
                                call interp_acc_surface(pos, weights, rt, &
                                    surface_height1, surface_height2)
                            endif
                            call push_particle_3d_ft(thread_id, rt, surface_height1, &
                                surface_height2, ptl, fields, db2, lc, kappa, &
                                deltax, deltay, deltaz, deltap, deltav, deltamu)
                        endif
                    else
                        if (ndim_field == 1) then
                            call push_particle_1d(thread_id, rt, ptl, fields, db2, &
                                lc, kappa, deltax, deltap)
                        else if (ndim_field == 2) then
                            if (include_3rd_dim_in2d_flag) then
                                call push_particle_2d_include_3rd(thread_id, rt, ptl, &
                                    fields, db2, lc, kappa, deltax, deltay, deltaz, deltap)
                            else
                                call push_particle_2d(thread_id, rt, ptl, fields, db2, &
                                    lc, kappa, deltax, deltay, deltap)
                            endif
                        else
                            if (acc_by_surface_flag) then
                                call interp_acc_surface(pos, weights, rt, &
                                    surface_height1, surface_height2)
                            endif
                            call push_particle_3d(thread_id, rt, surface_height1, &
                                surface_height2, ptl, fields, db2, lc, kappa, &
                                deltax, deltay, deltaz, deltap)
                        endif
                    endif

                    ! Number of particle tracking steps
                    ptl%nsteps_tracking = ptl%nsteps_tracking + 1

                    ! Track particles
                    if (track_particle_flag == 1) then
                        if (ptl%tag < 0 .and. &
                            mod(ptl%nsteps_tracking, nsteps_interval) == 0) then
                            tracking_step = ptl%nsteps_tracking / nsteps_interval
                            offset = noffsets_tracked_ptls(-ptl%tag)
                            ptl_traj_points(offset + tracking_step + 1) = ptl
                        endif
                    endif
                enddo ! while loop inside a fine time step

                ! Make sure ptl%t reach the targeted time exactly
                if ((ptl%t - t0) > dt_target .and. ptl%count_flag == COUNT_FLAG_INBOX) then
                    ptl%x = ptl%x - deltax
                    ptl%y = ptl%y - deltay
                    ptl%z = ptl%z - deltaz
                    ptl%p = ptl%p - deltap
                    ptl%v = ptl%v - deltav
                    ptl%mu = ptl%mu - deltamu
                    ptl%t = ptl%t - ptl%dt
                    ptl%dt = t0 + dt_target - ptl%t
                    if (ptl%dt > 0) then
                        ptl%nsteps_tracking = ptl%nsteps_tracking - 1

                        px = (ptl%x-xmin) / dxm
                        py = (ptl%y-ymin) / dym
                        pz = (ptl%z-zmin) / dzm
                        rt = (ptl%t - t0) / dtf
                        if (spherical_coord_flag) then
                            call get_interp_paramters_spherical(&
                                ptl%x, ptl%y, ptl%z, pos, weights)
                        else
                            call get_interp_paramters(px, py, pz, pos, weights)
                        endif
                        call interp_fields(pos, weights, rt, fields)
                        if (deltab_flag) then
                            call interp_magnetic_fluctuation(pos, weights, rt, db2)
                        endif
                        if (correlation_flag) then
                            call interp_correlation_length(pos, weights, rt, lc)
                        endif
                        call calc_spatial_diffusion_coefficients(ptl, focused_transport, &
                            fields, db2, lc, kappa)
                        if (focused_transport) then
                            if (ndim_field == 1) then
                                call push_particle_1d_ft(thread_id, rt, ptl, fields, db2, &
                                    lc, kappa, deltax, deltap, deltav, deltamu)
                            else if (ndim_field == 2) then
                                if (include_3rd_dim_in2d_flag) then
                                    call push_particle_2d_include_3rd_ft(thread_id, rt, ptl, &
                                        fields, db2, lc, kappa, deltax, deltay, deltaz, deltap, &
                                        deltav, deltamu)
                                else
                                    call push_particle_2d_ft(thread_id, rt, ptl, fields, db2, &
                                        lc, kappa, deltax, deltay, deltap, deltav, deltamu)
                                endif
                            else
                                if (acc_by_surface_flag) then
                                    call interp_acc_surface(pos, weights, rt, surface_height1, surface_height2)
                                endif
                                call push_particle_3d_ft(thread_id, rt, surface_height1, surface_height2, &
                                    ptl, fields, db2, lc, kappa, deltax, deltay, deltaz, deltap, &
                                    deltav, deltamu)
                            endif
                        else
                            if (ndim_field == 1) then
                                call push_particle_1d(thread_id, rt, ptl, fields, db2, &
                                    lc, kappa, deltax, deltap)
                            else if (ndim_field == 2) then
                                if (include_3rd_dim_in2d_flag) then
                                    call push_particle_2d_include_3rd(thread_id, rt, ptl, &
                                        fields, db2, lc, kappa, deltax, deltay, deltaz, deltap)
                                else
                                    call push_particle_2d(thread_id, rt, ptl, fields, db2, &
                                        lc, kappa, deltax, deltay, deltap)
                                endif
                            else
                                if (acc_by_surface_flag) then
                                    call interp_acc_surface(pos, weights, rt, surface_height1, surface_height2)
                                endif
                                call push_particle_3d(thread_id, rt, surface_height1, surface_height2, &
                                    ptl, fields, db2, lc, kappa, deltax, deltay, deltaz, deltap)
                            endif
                        endif
                        ptl%nsteps_tracking = ptl%nsteps_tracking + 1

                        ! Track particles
                        if (track_particle_flag == 1) then
                            if (ptl%tag < 0 .and. &
                                mod(ptl%nsteps_tracking, nsteps_interval) == 0) then
                                tracking_step = ptl%nsteps_tracking / nsteps_interval
                                offset = noffsets_tracked_ptls(-ptl%tag)
                                ptl_traj_points(offset + tracking_step + 1) = ptl
                            endif
                        endif
                    endif
                    if (ptl%p < 0.0) then
                        ptl%count_flag = COUNT_FLAG_OTHERS
                        !$OMP ATOMIC UPDATE
                        leak_negp = leak_negp + ptl%weight
                    else
                        call particle_boundary_condition(ptl, xmin1, xmax1, &
                            ymin1, ymax1, zmin1, zmax1)
                    endif
                endif

                dt_target = dt_target + dt_fine
            enddo ! Fine time step loop
            ptls(i) = ptl
        enddo ! Loop over particles
        !$OMP END PARALLEL
    end subroutine particle_mover_one_cycle

    !---------------------------------------------------------------------------
    !< Move particles using the MHD simulation data as background fields
    !< Args:
    !<  focused_transport: whether to the Focused Transport equation
    !<  track_particle_flag: whether to track particles, 0 for no, 1 for yes
    !<  nptl_selected: number of selected particles
    !<  nsteps_interval: save particle points every nsteps_interval
    !<  mhd_tframe: MHD time frame, starting from 1
    !<  num_fine_steps: number of fine time steps
    !<  dump_escaped: whether to dump escaped particles
    !---------------------------------------------------------------------------
    subroutine particle_mover(focused_transport, track_particle_flag, &
            nptl_selected, nsteps_interval, mhd_tframe, num_fine_steps, &
            dump_escaped)
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config, tstamps_mhd
        implicit none
        logical, intent(in) :: focused_transport
        integer, intent(in) :: track_particle_flag, nptl_selected, nsteps_interval
        integer, intent(in) :: mhd_tframe, num_fine_steps
        logical, intent(in) :: dump_escaped
        integer :: i, local_flag, global_flag, ncycle
        logical :: all_particles_in_box
        real(dp) :: t0, dtf, xmin, xmax, ymin, ymax, zmin, zmax
        type(particle_type) :: ptl
        all_particles_in_box = .false.
        nptl_old = 0

        t0 = tstamps_mhd(mhd_tframe)
        dtf = tstamps_mhd(mhd_tframe+1) - t0
        call set_dt_min_max(dtf)

        ncycle = 0
        local_flag = 0
        global_flag = 0

        do while (.not. all_particles_in_box)
            ncycle = ncycle + 1
            nsenders = 0
            nrecvers = 0
            if (nptl_old < nptl_current) then
                call particle_mover_one_cycle(t0, dtf, track_particle_flag, &
                    nptl_selected, nsteps_interval, num_fine_steps, &
                    focused_transport)
            endif
            call remove_particles(dump_escaped)
            call send_recv_particles
            call add_neighbor_particles  ! Also update nptl_old, nptl_current
            if (sum(nrecvers) > 0) then
                local_flag = sum(nrecvers)
            else
                local_flag = 0
            endif
            global_flag = 0
            call MPI_ALLREDUCE(local_flag, global_flag, 1, MPI_INTEGER, &
                MPI_SUM, mpi_sub_comm, ierr)
            if (global_flag > 0) then
                all_particles_in_box = .false.
            else
                all_particles_in_box = .true.
            endif
        enddo
        if (mpi_rank == master) then
            write(*, "(A, I0)") "Number of cycles: ", ncycle
        endif

        !< Send particles in ghost cells to neighbors.
        !< We do not need to push particles after they are sent to neighbors.
        !< To do this, we can reduce number of particles crossing local domain
        !< domain boundaries in the next step.
        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax
        zmin = fconfig%zmin
        zmax = fconfig%zmax
        nsenders = 0
        nrecvers = 0
        do i = 1, nptl_current
            ptl = ptls(i)
            if (ptl%p < 0.0) then
                ptl%count_flag = COUNT_FLAG_OTHERS
                leak_negp = leak_negp + ptl%weight
            else
                call particle_boundary_condition(ptl, xmin, xmax, ymin, ymax, zmin, zmax)
            endif
            ptls(i) = ptl
        enddo
        call remove_particles(dump_escaped)
        call send_recv_particles
        call add_neighbor_particles
    end subroutine particle_mover

    !---------------------------------------------------------------------------
    !< Particle boundary conditions
    !< Args:
    !<  plt: particle structure
    !<  xmin, xmax: min and max along the x-direction
    !<  ymin, ymax: min and max along the y-direction
    !<  zmin, zmax: min and max along the z-direction
    !---------------------------------------------------------------------------
    subroutine particle_boundary_condition(ptl, xmin, xmax, ymin, ymax, zmin, zmax)
        use simulation_setup_module, only: neighbors, mpi_ix, mpi_iy, mpi_iz, &
            mpi_sizex, mpi_sizey, mpi_sizez
        use mhd_config_module, only: mhd_config
        implicit none
        real(dp), intent(in) :: xmin, xmax, ymin, ymax, zmin, zmax
        type(particle_type), intent(inout) :: ptl

        if (ptl%x < xmin .and. ptl%count_flag == COUNT_FLAG_INBOX) then
            if (neighbors(1) < 0) then
                !$OMP ATOMIC UPDATE
                leak = leak + ptl%weight
                ptl%count_flag = COUNT_FLAG_ESCAPE
            else if (neighbors(1) == mpi_sub_rank) then
                ptl%x = ptl%x - xmin + xmax
            else if (neighbors(1) == mpi_sub_rank + mpi_sizex - 1) then
                ptl%x = ptl%x - mhd_config%xmin + mhd_config%xmax
                !$OMP CRITICAL
                nsenders(1) = nsenders(1) + 1
                senders(nsenders(1), 1) = ptl
                !$OMP END CRITICAL
                ptl%count_flag = COUNT_FLAG_OTHERS
            else
                !$OMP CRITICAL
                nsenders(1) = nsenders(1) + 1
                senders(nsenders(1), 1) = ptl
                !$OMP END CRITICAL
                ptl%count_flag = COUNT_FLAG_OTHERS
            endif
        else if (ptl%x > xmax .and. ptl%count_flag == COUNT_FLAG_INBOX) then
            if (neighbors(2) < 0) then
                !$OMP ATOMIC UPDATE
                leak = leak + ptl%weight
                ptl%count_flag = COUNT_FLAG_ESCAPE
            else if (neighbors(2) == mpi_sub_rank) then
                ptl%x = ptl%x - xmax + xmin
            else if (neighbors(2) == mpi_sub_rank - mpi_sizex + 1) then !< simulation boundary
                ptl%x = ptl%x - mhd_config%xmax + mhd_config%xmin
                !$OMP CRITICAL
                nsenders(2) = nsenders(2) + 1
                senders(nsenders(2), 2) = ptl
                !$OMP END CRITICAL
                ptl%count_flag = COUNT_FLAG_OTHERS
            else
                !$OMP CRITICAL
                nsenders(2) = nsenders(2) + 1
                senders(nsenders(2), 2) = ptl
                !$OMP END CRITICAL
                ptl%count_flag = COUNT_FLAG_OTHERS
            endif
        endif

        if (ndim_field > 1) then
            !< We need to make sure the count_flag is not set to 0.
            !< Otherwise, we met send the particles to two different neighbors.
            if (ptl%y < ymin .and. ptl%count_flag == COUNT_FLAG_INBOX) then
                if (neighbors(3) < 0) then
                    !$OMP ATOMIC UPDATE
                    leak = leak + ptl%weight
                    ptl%count_flag = COUNT_FLAG_ESCAPE
                else if (neighbors(3) == mpi_sub_rank) then
                    ptl%y = ptl%y - ymin + ymax
                else if (neighbors(3) == mpi_sub_rank + (mpi_sizey - 1) * mpi_sizex) then
                    ptl%y = ptl%y - mhd_config%ymin + mhd_config%ymax
                    !$OMP CRITICAL
                    nsenders(3) = nsenders(3) + 1
                    senders(nsenders(3), 3) = ptl
                    !$OMP END CRITICAL
                    ptl%count_flag = COUNT_FLAG_OTHERS
                else
                    !$OMP CRITICAL
                    nsenders(3) = nsenders(3) + 1
                    senders(nsenders(3), 3) = ptl
                    !$OMP END CRITICAL
                    ptl%count_flag = COUNT_FLAG_OTHERS
                endif
            else if (ptl%y > ymax .and. ptl%count_flag == COUNT_FLAG_INBOX) then
                if (neighbors(4) < 0) then
                    !$OMP ATOMIC UPDATE
                    leak = leak + ptl%weight
                    ptl%count_flag = COUNT_FLAG_ESCAPE
                else if (neighbors(4) == mpi_sub_rank) then
                    ptl%y = ptl%y - ymax + ymin
                else if (neighbors(4) == mpi_sub_rank - (mpi_sizey - 1) * mpi_sizex) then
                    ptl%y = ptl%y - mhd_config%ymax + mhd_config%ymin
                    !$OMP CRITICAL
                    nsenders(4) = nsenders(4) + 1
                    senders(nsenders(4), 4) = ptl
                    !$OMP END CRITICAL
                    ptl%count_flag = COUNT_FLAG_OTHERS
                else
                    !$OMP CRITICAL
                    nsenders(4) = nsenders(4) + 1
                    senders(nsenders(4), 4) = ptl
                    !$OMP END CRITICAL
                    ptl%count_flag = COUNT_FLAG_OTHERS
                endif
            endif
        endif

        if (ndim_field == 3 .or. (ndim_field == 2 .and. include_3rd_dim_in2d_flag)) then
            if (ptl%z < zmin .and. ptl%count_flag == COUNT_FLAG_INBOX) then
                if (neighbors(5) < 0) then
                    !$OMP ATOMIC UPDATE
                    leak = leak + ptl%weight
                    ptl%count_flag = COUNT_FLAG_ESCAPE
                else if (neighbors(5) == mpi_sub_rank) then
                    ptl%z = ptl%z - zmin + zmax
                else if (neighbors(5) == mpi_sub_rank + (mpi_sizez - 1) * mpi_sizey * mpi_sizex) then
                    ptl%z = ptl%z - mhd_config%zmin + mhd_config%zmax
                    !$OMP CRITICAL
                    nsenders(5) = nsenders(5) + 1
                    senders(nsenders(5), 5) = ptl
                    !$OMP END CRITICAL
                    ptl%count_flag = COUNT_FLAG_OTHERS
                else
                    !$OMP CRITICAL
                    nsenders(5) = nsenders(5) + 1
                    senders(nsenders(5), 5) = ptl
                    !$OMP END CRITICAL
                    ptl%count_flag = COUNT_FLAG_OTHERS
                endif
            else if (ptl%z > zmax .and. ptl%count_flag == COUNT_FLAG_INBOX) then
                if (neighbors(6) < 0) then
                    !$OMP ATOMIC UPDATE
                    leak = leak + ptl%weight
                    ptl%count_flag = COUNT_FLAG_ESCAPE
                else if (neighbors(6) == mpi_sub_rank) then
                    ptl%z = ptl%z - zmax + zmin
                else if (neighbors(6) == mpi_sub_rank - (mpi_sizez - 1) * mpi_sizey * mpi_sizex) then
                    ptl%z = ptl%z - mhd_config%zmax + mhd_config%zmin
                    !$OMP CRITICAL
                    nsenders(6) = nsenders(6) + 1
                    senders(nsenders(6), 6) = ptl
                    !$OMP END CRITICAL
                    ptl%count_flag = COUNT_FLAG_OTHERS
                else
                    !$OMP CRITICAL
                    nsenders(6) = nsenders(6) + 1
                    senders(nsenders(6), 6) = ptl
                    !$OMP END CRITICAL
                    ptl%count_flag = COUNT_FLAG_OTHERS
                endif
            endif
        endif
    end subroutine particle_boundary_condition

    !---------------------------------------------------------------------------
    !< Send and receiver particles from one neighbor
    !< This assumes MPI size along this direction is 1 or an even number >= 2.
    !< Args:
    !<  send_id, recv_id: sender and receiver ID
    !<  mpi_direc: MPI rank along one direction
    !---------------------------------------------------------------------------
    subroutine send_recv_particle_one_neighbor(send_id, recv_id, mpi_direc)
        use simulation_setup_module, only: neighbors
        implicit none
        integer, intent(in) :: send_id, recv_id, mpi_direc
        integer :: nsend, nrecv
        nrecv = 0
        nsend = nsenders(send_id)
        if (neighbors(send_id) /= mpi_sub_rank .and. neighbors(send_id) >= 0) then
            call MPI_SEND(nsend, 1, MPI_INTEGER, neighbors(send_id), &
                mpi_sub_rank, mpi_sub_comm, ierr)
        endif
        if (neighbors(recv_id) /= mpi_sub_rank .and. neighbors(recv_id) >= 0) then
            call MPI_RECV(nrecv, 1, MPI_INTEGER, &
                neighbors(recv_id), neighbors(recv_id), mpi_sub_comm, status, ierr)
        endif
        nrecvers(recv_id) = nrecv
        if (mod(mpi_direc, 2) == 0) then
            if (nsend > 0) then
                call MPI_SEND(senders(1:nsend, send_id), nsend, &
                    particle_datatype_mpi, neighbors(send_id), mpi_sub_rank, &
                    mpi_sub_comm, ierr)
            endif
            if (nrecv > 0) then
                call MPI_RECV(recvers(1:nrecv, recv_id), nrecv, &
                    particle_datatype_mpi, neighbors(recv_id), &
                    neighbors(recv_id), mpi_sub_comm, status, ierr)
            endif
        else
            if (nrecv > 0) then
                call MPI_RECV(recvers(1:nrecv, recv_id), nrecv, &
                    particle_datatype_mpi, neighbors(recv_id), &
                    neighbors(recv_id), mpi_sub_comm, status, ierr)
            endif
            if (nsend > 0) then
                call MPI_SEND(senders(1:nsend, send_id), nsend, &
                    particle_datatype_mpi, neighbors(send_id), mpi_sub_rank, &
                    mpi_sub_comm, ierr)
            endif
        endif
    end subroutine send_recv_particle_one_neighbor

    !---------------------------------------------------------------------------
    !< Send particles to neighbors
    !---------------------------------------------------------------------------
    subroutine send_recv_particles
        use simulation_setup_module, only: mpi_ix, mpi_iy, mpi_iz
        implicit none
        call send_recv_particle_one_neighbor(1, 2, mpi_ix) !< Right -> Left
        call send_recv_particle_one_neighbor(2, 1, mpi_ix) !< Left -> Right
        if (ndim_field > 1) then
            call send_recv_particle_one_neighbor(3, 4, mpi_iy) !< Top -> Bottom
            call send_recv_particle_one_neighbor(4, 3, mpi_iy) !< Bottom -> Top
            if (ndim_field == 3) then
                call send_recv_particle_one_neighbor(5, 6, mpi_iz) !< Front -> Back
                call send_recv_particle_one_neighbor(6, 5, mpi_iz) !< Back -> Front
            endif
        endif
    end subroutine send_recv_particles

    !---------------------------------------------------------------------------
    !< Calculate the spatial diffusion coefficients
    !< Args:
    !<  ptl: particle structure
    !<  focused_transport: whether to the focused transport equation
    !<  fields: fields and their gradients at particle position
    !<  db2: turbulence variance and its gradients at particle position
    !<  lc: turbulence correlation length and its gradients at particle position
    !<  kappa: kappa and related variables (For focused transport, kappa is
    !<         actually the perpendicular diffusion coefficients
    !---------------------------------------------------------------------------
    subroutine calc_spatial_diffusion_coefficients(ptl, focused_transport, &
            fields, db2, lc, kappa)
        implicit none
        type(particle_type), intent(in) :: ptl
        logical, intent(in) :: focused_transport
        real(dp), dimension(*), intent(in) :: fields
        real(dp), dimension(*), intent(in) :: db2, lc
        type(kappa_type), intent(out) :: kappa
        real(dp) :: knorm
        real(dp) :: bx, by, bz, b, ib1, ib2, ib3, ib4
        real(dp) :: dbx_dx, dby_dx, dbz_dx
        real(dp) :: dbx_dy, dby_dy, dbz_dy
        real(dp) :: dbx_dz, dby_dz, dbz_dz
        real(dp) :: db_dx, db_dy, db_dz
        real(dp) :: dkdx, dkdy, dkdz
        real(dp) :: kpp

        bx = fields(5)
        by = fields(6)
        bz = fields(7)
        b = dsqrt(bx**2 + by**2 + bz**2)
        if (b == 0) then
            ib1 = 1.0
        else
            ib1 = 1.0_dp / b
        endif
        ib2 = ib1 * ib1
        ib3 = ib1 * ib2
        ib4 = ib2 * ib2

        kappa%knorm0 = 1.0_dp
        if (mag_dependency == 1) then
            kappa%knorm0 = kappa%knorm0 * ib1**(1./3.)
        endif

        ! Magnetic fluctuation dB^2/B^2
        if (deltab_flag) then
            kappa%knorm0 = kappa%knorm0 / db2(1)
        endif

        ! Turbulence correlation length
        ! Make sure that lc is non-zero in the data file!!!
        if (correlation_flag) then
            kappa%knorm0 = kappa%knorm0 * lc(1)**(2./3.)
        endif

        ! Momentum-dependent kappa
        if (momentum_dependency == 1) then
            knorm = kappa%knorm0 * (ptl%p / p0)**pindex
        else
            knorm = kappa%knorm0
        endif

        kappa%kpara = kpara0 * knorm
        kappa%kperp = kappa%kpara * kret

        kappa%skpara = dsqrt(2.0 * kappa%kpara)
        kappa%skperp = dsqrt(2.0 * kappa%kperp)
        kappa%skpara_perp = dsqrt(2.0 * (kappa%kpara - kappa%kperp))

        if (ndim_field == 1) then
            dkdx = 0.0_dp
            if (mag_dependency == 1) then
                dkdx = -db_dx * ib1 / 3
            endif
            if (deltab_flag) then
                dkdx = dkdx - db2(2) / db2(1)
            endif
            if (correlation_flag) then
                dkdx = dkdx + 2 * lc(2) / (3 * lc(1))
            endif
            if (focused_transport) then
                kappa%dkxx_dx = kappa%kperp * dkdx
                if (spherical_coord_flag) then
                    kappa%kxx = kappa%kperp
                endif
            else
                kappa%dkxx_dx = kappa%kpara * dkdx
                if (spherical_coord_flag) then
                    kappa%kxx = kappa%kpara
                endif
            endif
        else if (ndim_field == 2) then
            if (include_3rd_dim_in2d_flag) then
                dbx_dx = fields(nfields+13)
                dbx_dy = fields(nfields+14)
                dbx_dz = 0.0_dp
                dby_dx = fields(nfields+16)
                dby_dy = fields(nfields+17)
                dby_dz = 0.0_dp
                dbz_dx = fields(nfields+19)
                dbz_dy = fields(nfields+20)
                dbz_dz = 0.0_dp
                db_dx = fields(nfields+22)
                db_dy = fields(nfields+23)
                db_dz = 0.0_dp
                dkdx = 0.0_dp
                dkdy = 0.0_dp
                dkdz = 0.0_dp
                if (mag_dependency == 1) then
                    dkdx = -db_dx * ib1 / 3
                    dkdy = -db_dy * ib1 / 3
                endif
                if (deltab_flag) then
                    dkdx = dkdx - db2(2) / db2(1)
                    dkdy = dkdy - db2(3) / db2(1)
                endif
                if (correlation_flag) then
                    dkdx = dkdx + 2 * lc(2) / (3 * lc(1))
                    dkdy = dkdy + 2 * lc(3) / (3 * lc(1))
                endif
                if (focused_transport) then
                    kpp = kappa%kpara - kappa%kperp
                else
                    kpp = -kappa%kperp
                endif
                kappa%dkxx_dx = kappa%kperp*dkdx + kpp*dkdx*bx**2*ib2 + &
                    2.0*kpp*bx*(dbx_dx*b-bx*db_dx)*ib3
                kappa%dkyy_dy = kappa%kperp*dkdy + kpp*dkdy*by**2*ib2 + &
                    2.0*kpp*by*(dby_dy*b-by*db_dy)*ib3
                kappa%dkzz_dz = kappa%kperp*dkdz + kpp*dkdz*bz**2*ib2 + &
                    2.0*kpp*bz*(dbz_dz*b-bz*db_dz)*ib3
                kappa%dkxy_dx = kpp*dkdx*bx*by*ib2 + kpp * &
                    ((dbx_dx*by+bx*dby_dx)*ib2 - 2.0*bx*by*db_dx*ib3)
                kappa%dkxy_dy = kpp*dkdy*bx*by*ib2 + kpp * &
                    ((dbx_dy*by+bx*dby_dy)*ib2 - 2.0*bx*by*db_dy*ib3)
                kappa%dkxz_dx = kpp*dkdx*bx*bz*ib2 + kpp * &
                    ((dbx_dx*bz+bx*dbz_dx)*ib2 - 2.0*bx*bz*db_dx*ib3)
                kappa%dkxz_dz = kpp*dkdz*bx*bz*ib2 + kpp * &
                    ((dbx_dz*bz+bx*dbz_dz)*ib2 - 2.0*bx*bz*db_dz*ib3)
                kappa%dkyz_dy = kpp*dkdy*by*bz*ib2 + kpp * &
                    ((dby_dy*bz+by*dbz_dy)*ib2 - 2.0*by*bz*db_dy*ib3)
                kappa%dkyz_dz = kpp*dkdz*by*bz*ib2 + kpp * &
                    ((dby_dz*bz+by*dbz_dz)*ib2 - 2.0*by*bz*db_dz*ib3)
                kappa%kxx = kappa%kperp + kpp * bx * bx * ib2
                kappa%kyy = kappa%kperp + kpp * by * by * ib2
                kappa%kzz = kappa%kperp + kpp * bz * bz * ib2
                kappa%kxy = kpp * bx * by * ib2
                kappa%kxz = kpp * bx * bz * ib2
                kappa%kyz = kpp * by * bz * ib2
            else
                dbx_dx = fields(nfields+13)
                dbx_dy = fields(nfields+14)
                dby_dx = fields(nfields+16)
                dby_dy = fields(nfields+17)
                db_dx = fields(nfields+22)
                db_dy = fields(nfields+23)
                dkdx = 0.0_dp
                dkdy = 0.0_dp
                if (mag_dependency == 1) then
                    dkdx = -db_dx * ib1 / 3
                    dkdy = -db_dy * ib1 / 3
                endif
                if (deltab_flag) then
                    dkdx = dkdx - db2(2) / db2(1)
                    dkdy = dkdy - db2(3) / db2(1)
                endif
                if (correlation_flag) then
                    dkdx = dkdx + 2 * lc(2) / (3 * lc(1))
                    dkdy = dkdy + 2 * lc(3) / (3 * lc(1))
                endif
                if (focused_transport) then
                    kpp = kappa%kpara - kappa%kperp
                else
                    kpp = -kappa%kperp
                endif
                kappa%dkxx_dx = kappa%kperp*dkdx + kpp*dkdx*bx**2*ib2 + &
                    2.0*kpp*bx*(dbx_dx*b-bx*db_dx)*ib3
                kappa%dkyy_dy = kappa%kperp*dkdy + kpp*dkdy*by**2*ib2 + &
                    2.0*kpp*by*(dby_dy*b-by*db_dy)*ib3
                kappa%dkxy_dx = kpp*dkdx*bx*by*ib2 + kpp * &
                    ((dbx_dx*by+bx*dby_dx)*ib2 - 2.0*bx*by*db_dx*ib3)
                kappa%dkxy_dy = kpp*dkdy*bx*by*ib2 + kpp * &
                    ((dbx_dy*by+bx*dby_dy)*ib2 - 2.0*bx*by*db_dy*ib3)
                kappa%kxx = kappa%kperp + kpp * bx * bx * ib2
                kappa%kyy = kappa%kperp + kpp * by * by * ib2
                kappa%kxy = kpp * bx * by * ib2
            endif
        else
            dbx_dx = fields(nfields+13)
            dbx_dy = fields(nfields+14)
            dbx_dz = fields(nfields+15)
            dby_dx = fields(nfields+16)
            dby_dy = fields(nfields+17)
            dby_dz = fields(nfields+18)
            dbz_dx = fields(nfields+19)
            dbz_dy = fields(nfields+20)
            dbz_dz = fields(nfields+21)
            db_dx = fields(nfields+22)
            db_dy = fields(nfields+23)
            db_dz = fields(nfields+24)
            dkdx = 0.0_dp
            dkdy = 0.0_dp
            dkdz = 0.0_dp
            if (mag_dependency == 1) then
                dkdx = -db_dx * ib1 / 3
                dkdy = -db_dy * ib1 / 3
                dkdz = -db_dz * ib1 / 3
            endif
            if (deltab_flag) then
                dkdx = dkdx - db2(2) / db2(1)
                dkdy = dkdy - db2(3) / db2(1)
                dkdz = dkdz - db2(4) / db2(1)
            endif
            if (correlation_flag) then
                dkdx = dkdx + 2 * lc(2) / (3 * lc(1))
                dkdy = dkdy + 2 * lc(3) / (3 * lc(1))
                dkdz = dkdz + 2 * lc(4) / (3 * lc(1))
            endif
            if (focused_transport) then
                kpp = kappa%kpara - kappa%kperp
            else
                kpp = -kappa%kperp
            endif
            kappa%dkxx_dx = kappa%kperp*dkdx + kpp*dkdx*bx**2*ib2 + &
                2.0*kpp*bx*(dbx_dx*b-bx*db_dx)*ib3
            kappa%dkyy_dy = kappa%kperp*dkdy + kpp*dkdy*by**2*ib2 + &
                2.0*kpp*by*(dby_dy*b-by*db_dy)*ib3
            kappa%dkzz_dz = kappa%kperp*dkdz + kpp*dkdz*bz**2*ib2 + &
                2.0*kpp*bz*(dbz_dz*b-bz*db_dz)*ib3
            kappa%dkxy_dx = kpp*dkdx*bx*by*ib2 + kpp * &
                ((dbx_dx*by+bx*dby_dx)*ib2 - 2.0*bx*by*db_dx*ib3)
            kappa%dkxy_dy = kpp*dkdy*bx*by*ib2 + kpp * &
                ((dbx_dy*by+bx*dby_dy)*ib2 - 2.0*bx*by*db_dy*ib3)
            kappa%dkxz_dx = kpp*dkdx*bx*bz*ib2 + kpp * &
                ((dbx_dx*bz+bx*dbz_dx)*ib2 - 2.0*bx*bz*db_dx*ib3)
            kappa%dkxz_dz = kpp*dkdz*bx*bz*ib2 + kpp * &
                ((dbx_dz*bz+bx*dbz_dz)*ib2 - 2.0*bx*bz*db_dz*ib3)
            kappa%dkyz_dy = kpp*dkdy*by*bz*ib2 + kpp * &
                ((dby_dy*bz+by*dbz_dy)*ib2 - 2.0*by*bz*db_dy*ib3)
            kappa%dkyz_dz = kpp*dkdz*by*bz*ib2 + kpp * &
                ((dby_dz*bz+by*dbz_dz)*ib2 - 2.0*by*bz*db_dz*ib3)
            kappa%kxx = kappa%kperp + kpp * bx * bx * ib2
            kappa%kyy = kappa%kperp + kpp * by * by * ib2
            kappa%kzz = kappa%kperp + kpp * bz * bz * ib2
            kappa%kxy = kpp * bx * by * ib2
            kappa%kxz = kpp * bx * bz * ib2
            kappa%kyz = kpp * by * bz * ib2
        endif
    end subroutine calc_spatial_diffusion_coefficients

    !---------------------------------------------------------------------------
    !< Read particle parameters including the diffusion coefficients
    !< Args:
    !<  conf_file: configuration file name
    !---------------------------------------------------------------------------
    subroutine read_particle_params(conf_file)
        use read_config, only: get_variable
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: dt_mhd_max
        implicit none
        character(*), intent(in) :: conf_file
        logical :: condx, condy, condz
        real(fp) :: temp
        integer :: fh

        if (mpi_rank == master) then
            fh = 10
            open(unit=fh, file='config/'//trim(conf_file), status='old')
            b0 = get_variable(fh, 'b0', '=')
            p0 = get_variable(fh, 'p0', '=')
            pmin = get_variable(fh, 'pmin', '=')
            pmax = get_variable(fh, 'pmax', '=')
            temp = get_variable(fh, 'momentum_dependency', '=')
            momentum_dependency = int(temp)
            if (momentum_dependency == 1) then
                pindex = get_variable(fh, 'pindex', '=')
            endif
            temp = get_variable(fh, 'mag_dependency', '=')
            mag_dependency = int(temp)
            kpara0 = get_variable(fh, 'kpara0', '=')
            kret = get_variable(fh, 'kret', '=')
            dt_min_rel = get_variable(fh, 'dt_min_rel', '=')
            dt_max_rel = get_variable(fh, 'dt_max_rel', '=')

            acc_region_flag = 0  !< In default, particles can be accelerated everywhere
            temp = get_variable(fh, 'acc_region_flag', '=')
            acc_region_flag = int(temp)
            acc_region(1) = get_variable(fh, 'acc_xmin', '=')
            acc_region(2) = get_variable(fh, 'acc_xmax', '=')
            acc_region(3) = get_variable(fh, 'acc_ymin', '=')
            acc_region(4) = get_variable(fh, 'acc_ymax', '=')
            acc_region(5) = get_variable(fh, 'acc_zmin', '=')
            acc_region(6) = get_variable(fh, 'acc_zmax', '=')
            close(fh)
            !< echo the information
            print *, "---------------------------------------------------"
            write(*, "(A)") " Particle parameters including diffusion coefficients"
            write(*, "(A, E13.6E2)") " The standard deviation of momentum distribution is ", p0
            write(*, "(A, E13.6E2, E13.6E2)") &
                " Minimum and maximum particle momentum", pmin, pmax
            write(*, "(A, E13.6E2)") " Initial magnetic field strength", b0
            write(*, "(A,I0)") " kappa dependency on p: ", momentum_dependency
            if (momentum_dependency == 1) then
                write(*, "(A,E13.6E2)") " kappa ~ p^", pindex
            else
                write(*, "(A)") " kappa doesn't depend on particle momentum"
            endif
            if (mag_dependency == 1) then
                write(*, "(A)") " kappa ~ 1/B"
            else
                write(*, "(A)") " kappa doesn't depend on B"
            endif
            write(*, "(A,E13.6E2)") " kpara0 = ", kpara0
            write(*, "(A,E13.6E2)") " kperp / kpara = ", kret
            write(*, "(A,E13.6E2)") " Minimum relative time step = ", dt_min_rel
            write(*, "(A,E13.6E2)") " Maximum relative time step = ", dt_max_rel
            if (acc_region_flag == 1) then
                write(*, "(A)") " Particle can only be accelerated in part of the box"
                write(*, "(A,F13.6)") " acc_xmin = ", acc_region(1)
                write(*, "(A,F13.6)") " acc_xmax = ", acc_region(2)
                write(*, "(A,F13.6)") " acc_ymin = ", acc_region(3)
                write(*, "(A,F13.6)") " acc_ymax = ", acc_region(4)
                write(*, "(A,F13.6)") " acc_zmin = ", acc_region(5)
                write(*, "(A,F13.6)") " acc_zmax = ", acc_region(6)
            endif
            print *, "---------------------------------------------------"
        endif
        call MPI_BCAST(momentum_dependency, 1, MPI_INTEGER, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_BCAST(mag_dependency, 1, MPI_INTEGER, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_BCAST(pindex, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(p0, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(pmin, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(pmax, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(b0, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(kpara0, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(kret, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(dt_min_rel, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(dt_max_rel, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(acc_region_flag, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)

        ! These might change if MHD time interval is varying
        call set_dt_min_max(dt_mhd_max)

        if (acc_region_flag == 1) then
            call MPI_BCAST(acc_region, 6, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        else
            acc_region(1) = 0.0
            acc_region(2) = 1.0
            acc_region(3) = 0.0
            acc_region(4) = 1.0
            acc_region(5) = 0.0
            acc_region(6) = 1.0
        endif
    end subroutine read_particle_params

    !---------------------------------------------------------------------------
    !< Determine the time step. X, Y, Z are for r, theta, phi in
    !< spherical coordinates
    !< Args;
    !<  t0: the initial time for current particle
    !<  dtf: the time interval between fine diagnostics
    !<  fields: fields and their gradients at particle position
    !<  kappa: kappa and related variables
    !<  ptl: particle structure
    !---------------------------------------------------------------------------
    subroutine set_time_step(t0, dtf, fields, kappa, ptl)
        use constants, only: pi
        use mhd_config_module, only: mhd_config
        implicit none
        real(dp), intent(in) :: t0, dtf
        real(dp), dimension(*), intent(in) :: fields
        type(kappa_type), intent(in) :: kappa
        type(particle_type), intent(inout) :: ptl
        real(dp) :: tmp30, tmp40, bx, by, bz, b, ib, ib2, ib3
        real(dp) :: vx, vy, vz, dxm, dym, dzm, dt1, dt2, dt3
        real(dp) :: vdx, vdy, vdz, vdp
        real(dp) :: dbx_dy, dbx_dz, dby_dx, dby_dz, dbz_dx, dbz_dy
        real(dp) :: db_dx, db_dy, db_dz, bxn, byn, bzn, bxyn, ibxyn
        real(dp) :: a1, b1, c1, Qpp, Qpm, Qmp, Qmm, ctheta, istheta, ir, ir2
        real(dp) :: qtmp1, qtmp2, atmp
        real(dp) :: gbr, gbt, gbp, p11, p12, p13, p22, p23, p33

        if (ndim_field == 1) then
            vx = fields(1)
            dxm = mhd_config%dx
            tmp30 = kappa%skpara
            if (spherical_coord_flag) then
                tmp40 = abs(vx + kappa%dkxx_dx + 2.0*kappa%kxx/ptl%x)
            else
                tmp40 = abs(vx + kappa%dkxx_dx)
            endif
            if (tmp40 .ne. 0.0d0) then
                if (tmp30 > 0) then
                    ptl%dt = min(dxm/(10.0*tmp40), (tmp30/tmp40)**2) * 0.5d0
                else
                    ptl%dt = dxm/(10.0*tmp40)
                endif
            else
                ptl%dt = dt_min
            endif
        else if (ndim_field == 2) then
            vx = fields(1)
            vy = fields(2)
            bx = fields(5)
            by = fields(6)
            bz = fields(7)
            b = dsqrt(bx**2 + by**2 + bz**2)
            if (b == 0.0d0) then
                ib = 0.0_dp
            else
                ib = 1.0_dp / b
            endif
            bxn = bx * ib
            byn = by * ib
            dxm = mhd_config%dx
            dym = mhd_config%dy

            if (spherical_coord_flag) then
                ctheta = cos(ptl%y)
                istheta = 1.0 / sin(ptl%y)
                ir = 1.0 / ptl%x
                ir2 = 1.0 / ptl%x**2
                a1 = kappa%kxx
                b1 = kappa%kxy * ir
                c1 = kappa%kyy * ir2
                atmp = dsqrt((a1-c1)**2 + 4*b1**2)
                Qpp = atmp + (a1 + c1)
                Qmp = atmp - (a1 + c1)
                Qpm = atmp + (a1 - c1)
                Qmm = atmp - (a1 - c1)
                qtmp1 = dsqrt(-Qmp/(Qmm**2+4*b1**2))
                qtmp2 = dsqrt(Qpp/(Qpm**2+4*b1**2))
            endif

            ! X-direction
            if (spherical_coord_flag) then
                tmp30 = abs(-Qmm * qtmp1) + abs(Qpm * qtmp2)
                tmp40 = abs(vx + kappa%dkxx_dx + &
                    (2.0*kappa%kxx + kappa%dkxy_dy + kappa%kxy*ctheta*istheta)*ir)
            else
                tmp30 = kappa%skperp + kappa%skpara_perp * abs(bxn)
                tmp40 = abs(vx + kappa%dkxx_dx + kappa%dkxy_dy)
            endif
            if (tmp40 .ne. 0.0d0) then
                if (tmp30 > 0) then
                    ! dt1 = min(dxm/(10.0*tmp40), (tmp30/tmp40)**2) * 0.5d0
                    ! dt1 = min((tmp30/tmp40)**2, dxm**2/tmp30**2) * 0.3
                    dt1 = min(2 * kappa%kxx/tmp40**2, 0.5 * dxm**2/kappa%kxx) * 0.3
                    ! dt1 = min(dxm**2/kappa%kxx, dxm/tmp40) * 0.3
                else
                    dt1 = dxm/(10.0*tmp40)
                endif
            else
                dt1 = dt_min
            endif

            ! Y-direction
            if (spherical_coord_flag) then
                tmp30 = abs(2*b1*qtmp1) + abs(2*b1*qtmp2)
                tmp40 = abs((vy+kappa%dkxy_dx)*ir + &
                    (kappa%kxy + kappa%dkyy_dy + kappa%kyy*ctheta*istheta)*ir2)
            else
                tmp30 = kappa%skperp + kappa%skpara_perp * abs(byn)
                tmp40 = abs(vy + kappa%dkxy_dx + kappa%dkyy_dy)
            endif
            if (tmp40 .ne. 0.0d0) then
                if (tmp30 > 0) then
                    ! dt2 = min(dym/(10.0*tmp40), (tmp30/tmp40)**2) * 0.5d0
                    ! dt2 = min((tmp30/tmp40)**2, dym**2/tmp30**2) * 0.3
                    dt2 = min(2 * kappa%kyy/tmp40**2, 0.5 * dym**2/kappa%kyy) * 0.3
                    ! dt2 = min(dym**2/kappa%kyy, dym/tmp40) * 0.3
                else
                    dt2 = dym/(10.0*tmp40)
                endif
            else
                dt2 = dt_min
            endif
            ptl%dt = min(dt1, dt2)
        else
            vx = fields(1)
            vy = fields(2)
            vz = fields(3)
            bx = fields(5)
            by = fields(6)
            bz = fields(7)
            b = dsqrt(bx**2 + by**2 + bz**2)
            if (b == 0.0d0) then
                ib = 0.0_dp
            else
                ib = 1.0_dp / b
            endif
            bxn = bx * ib
            byn = by * ib
            bzn = bz * ib
            bxyn = dsqrt(bxn**2 + byn**2)
            if (bxyn == 0.0d0) then
                ibxyn = 0.0d0
            else
                ibxyn = 1.0_dp / bxyn
            endif
            dxm = mhd_config%dx
            dym = mhd_config%dy
            dzm = mhd_config%dz

            if (spherical_coord_flag) then
                ctheta = cos(ptl%y)
                istheta = 1.0 / sin(ptl%y)
                ir = 1.0 / ptl%x
                ir2 = ir**2
                p11 = dsqrt(2.0*(kappa%kxx*kappa%kyz**2 + &
                                 kappa%kyy*kappa%kxz**2 + &
                                 kappa%kzz*kappa%kxy**2 - &
                                 2.0*kappa%kxy*kappa%kxz*kappa%kyz - &
                                 kappa%kxx*kappa%kyy*kappa%kzz) / &
                                 (kappa%kyz**2 - kappa%kyy*kappa%kzz))
                p12 = (kappa%kxz*kappa%kyz - kappa%kxy*kappa%kzz) * &
                    dsqrt(2.0*(kappa%kyy - (kappa%kyz**2/kappa%kzz))) / &
                    (kappa%kyz**2 - kappa%kyy*kappa%kzz)
                p13 = dsqrt(2.0 / kappa%kzz) * kappa%kxz
                p22 = dsqrt(2.0 * (kappa%kyy - kappa%kyz**2/kappa%kzz)) * ir
                p23 = dsqrt(2.0 / kappa%kzz) * kappa%kyz * ir
                p33 = dsqrt(2.0 * kappa%kzz) * istheta * ir
            endif

            ! Drift velocity
            dbx_dy = fields(nfields+14)
            dbx_dz = fields(nfields+15)
            dby_dx = fields(nfields+16)
            dby_dz = fields(nfields+18)
            dbz_dx = fields(nfields+19)
            dbz_dy = fields(nfields+20)
            db_dx = fields(nfields+22)
            db_dy = fields(nfields+23)
            db_dz = fields(nfields+24)
            ib2 = ib * ib
            ib3 = ib * ib2
            vdp = 1.0 / (3 * pcharge) / dsqrt((drift1*p0/ptl%p)**2 + (drift2*p0**2/ptl%p**2)**2)
            if (spherical_coord_flag) then
                gbr = db_dx
                gbt = db_dy * ir
                gbp = db_dz * ir * istheta
                vdx = vdp * ((dbz_dy + bz*ctheta*istheta - dby_dz*istheta)*ir*ib2 - &
                             (gbt*bz - gbp*by)*2.0*ib3)
                vdy = vdp * ((dbx_dz*istheta*ir - dbz_dx - bz*ir)*ib2 - &
                             (gbp*bx - gbr*bz)*2.0*ib3)
                vdz = vdp * ((dby_dx + by*ir - dbx_dy*ir)*ib2 - &
                             (gbr*by - gbt*bx)*2.0*ib3)
            else
                vdx = vdp * ((dbz_dy-dby_dz)*ib2 - 2*(bz*db_dy-by*db_dz)*ib3)
                vdy = vdp * ((dbx_dz-dbz_dx)*ib2 - 2*(bx*db_dz-bz*db_dx)*ib3)
                vdz = vdp * ((dby_dx-dbx_dy)*ib2 - 2*(by*db_dx-bx*db_dy)*ib3)
            endif

            ! X-direction
            if (spherical_coord_flag) then
                tmp30 = abs(p11) + abs(p12) + abs(p13)
                tmp40 = abs(vx + vdx + kappa%dkxx_dx + &
                    (2.0*kappa%kxx + kappa%dkxy_dy + &
                     kappa%kxy*ctheta*istheta + &
                     kappa%dkxz_dz*istheta)*ir)
            else
                tmp30 = abs(bxn)*kappa%skpara + &
                    abs(bxn*bzn)*kappa%skperp*ibxyn + abs(byn)*kappa%skperp*ibxyn
                tmp40 = abs(vx + vdx + kappa%dkxx_dx + kappa%dkxy_dy + kappa%dkxz_dz)
            endif
            if (tmp40 .ne. 0.0d0) then
                if (tmp30 > 0) then
                    dt1 = min(dxm/(10.0*tmp40), (tmp30/tmp40)**2) * 0.5d0
                else
                    dt1 = dxm/(10.0*tmp40)
                endif
            else
                dt1 = dt_min
            endif

            ! Y-direction
            if (spherical_coord_flag) then
                tmp30 = abs(p22) + abs(p33)
                tmp40 = abs((vy + vdy + kappa%dkxy_dx)*ir + &
                    (kappa%kxy + kappa%dkyy_dy + &
                     kappa%kyy*ctheta*istheta + kappa%dkyz_dz*istheta)*ir2)
            else
                tmp30 = abs(byn)*kappa%skpara + &
                    abs(byn*bzn)*kappa%skperp*ibxyn + abs(bxn)*kappa%skperp*ibxyn
                tmp40 = abs(vy + vdy + kappa%dkxy_dx + kappa%dkyy_dy + kappa%dkyz_dz)
            endif
            if (tmp40 .ne. 0.0d0) then
                if (tmp30 > 0) then
                    dt2 = min(dym/(10.0*tmp40), (tmp30/tmp40)**2) * 0.5d0
                else
                    dt2 = dym/(10.0*tmp40)
                endif
            else
                dt2 = dt_min
            endif

            ! Z-direction
            if (spherical_coord_flag) then
                tmp30 = abs(p33)
                tmp40 = abs((vz + vdz + kappa%dkxz_dx)*istheta*ir + &
                    (kappa%kxz + kappa%dkyz_dy + kappa%dkzz_dz*istheta)*istheta*ir2)
            else
                tmp30 = abs(bzn)*kappa%skpara + bxyn*kappa%skperp
                tmp40 = abs(vz + vdz + kappa%dkxz_dx + kappa%dkyz_dy + kappa%dkzz_dz)
            endif
            if (tmp40 .ne. 0.0d0) then
                if (tmp30 > 0) then
                    dt3 = min(dym/(10.0*tmp40), (tmp30/tmp40)**2) * 0.5d0
                else
                    dt3 = dym/(10.0*tmp40)
                endif
            else
                dt3 = dt_min
            endif

            ptl%dt = min(dt1, dt2, dt3)
        endif

        !!< Make sure the time step is not too large. Adding dt_min to make
        !!< sure to exit the while where this routine is called
        !if ((ptl%t + ptl%dt - t0) > dtf) then
        !    ! ptl%dt = t0 + dtf - ptl%t + dt_min
        !    ptl%dt = t0 + dtf - ptl%t
        !endif

        !< Make sure the time step is not too small
        if (ptl%dt .lt. dt_min) then
            ptl%dt = dt_min
        endif

        !< Make sure the time step is not too large
        if (ptl%dt .gt. dt_max) then
            ptl%dt = dt_max
        endif
    end subroutine set_time_step

    !---------------------------------------------------------------------------
    !< Check whether particle in the in the acceleration region (rectangular box)
    !< Args:
    !<  ptl: one particle
    !---------------------------------------------------------------------------
    function particle_in_acceleration_region(ptl) result (ptl_in_region)
        use mhd_config_module, only: mhd_config
        implicit none
        type(particle_type), intent(in) :: ptl
        logical :: ptl_in_region, inx, iny, inz
        real(dp) :: xnorm, ynorm, znorm  !< 0 - 1
        ptl_in_region = .true.
        inx = .true.
        iny = .true.
        inz = .true.
        xnorm = (ptl%x - mhd_config%xmin) / mhd_config%lx
        inx = (xnorm >= acc_region(1)) .and. (xnorm <= acc_region(2))
        if (ndim_field > 1) then
            ynorm = (ptl%y - mhd_config%ymin) / mhd_config%ly
            iny = (ynorm >= acc_region(3)) .and. (ynorm <= acc_region(4))
        endif
        if (ndim_field == 3) then
            znorm = (ptl%z - mhd_config%zmin) / mhd_config%lz
            inz = (znorm >= acc_region(5)) .and. (znorm <= acc_region(6))
        endif
        ptl_in_region = inx .and. iny .and. inz
    end function particle_in_acceleration_region

    !---------------------------------------------------------------------------
    !< Push particle for a single step for a 1D simulation
    !< Args:
    !<  thread_id: thread ID staring from 0
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  ptl: particle structure
    !<  fields: fields and their gradients at particle position
    !<  db2: turbulence variance and its gradients at particle position
    !<  lc: turbulence correlation length and its gradients at particle position
    !<  kappa: kappa and related variables
    !<  deltax, deltap: the change of x and p in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_1d(thread_id, rt, ptl, fields, db2, &
            lc, kappa, deltax, deltap)
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: interp_fields, &
            interp_magnetic_fluctuation, interp_correlation_length
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        real(dp), intent(in) :: rt
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        real(dp), dimension(*), intent(inout) :: db2, lc
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltap
        real(dp) :: xtmp
        real(dp) :: sdt, dvx_dx, divv, gshear
        real(dp) :: b, ib, bx, by, bz, vx, px, rt1
        real(dp) :: xmin, xmax, dxm
        reaL(dp) :: xmin1, xmax1, dxmh
        real(dp) :: skpara, skpara1
        real(dp) :: ran1, sqrt3
        real(dp) :: rho, va ! Plasma density and Alfven speed
        real(dp) :: rands(2)
        real(dp) :: sigmaxx, sigmayy, sigmazz ! shear tensor
        real(dp) :: bbsigma ! b_ib_jsigma_ij
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights

        vx = fields(1)
        b = fields(8)
        dvx_dx = fields(nfields+1)
        xmin = fconfig%xmin
        xmax = fconfig%xmax
        dxm = mhd_config%dx
        dxmh = 0.5 * dxm

        !< The field data has two ghost cells, so the particles can cross the
        !< boundary without causing segment fault errors
        xmin1 = xmin - dxmh
        xmax1 = xmax + dxmh

        deltax = 0.0
        deltap = 0.0

        sdt = dsqrt(ptl%dt)
        if (spherical_coord_flag) then
            deltax = (vx + kappa%dkxx_dx + 2.0*kappa%kxx/ptl%x) * ptl%dt
        else
            deltax = (vx + kappa%dkxx_dx)*ptl%dt
        endif
        xtmp = ptl%x + deltax + kappa%skpara*sdt
        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3

        !< We originally tried to decrease the time step when xtmp or ytmp are out-of-bound,
        !< but decreasing the time step does not necessarily make the moving distance smaller.
        !< Therefore, we switch between first-order and second-order method.
        if (xtmp < xmin1 .or. xtmp > xmax1) then
            !< First-order method
            deltax = deltax + ran1*kappa%skpara*sdt
        else
            !< Second-order method. It requires xtmp and ytmp are in the local domain.
            px = (xtmp - xmin) / dxm
            if (spherical_coord_flag) then
                call get_interp_paramters_spherical(ptl%x, 0.0_dp, 0.0_dp, pos, weights)
            else
                call get_interp_paramters(px, 0.0_dp, 0.0_dp, pos, weights)
            endif
            call interp_fields(pos, weights, rt, fields)
            if (deltab_flag) then
                call interp_magnetic_fluctuation(pos, weights, rt, db2)
            endif
            if (correlation_flag) then
                call interp_correlation_length(pos, weights, rt, lc)
            endif

            skpara = kappa%skpara

            call calc_spatial_diffusion_coefficients(ptl, .false., fields, &
                db2, lc, kappa)

            skpara1 = kappa%skpara ! Diffusion coefficient at predicted position

            deltax = deltax + ran1*skpara*sdt + &
                (skpara1-skpara)*(ran1*ran1-1.0)*sdt/2.0
        endif

        ptl%x = ptl%x + deltax
        ptl%t = ptl%t + ptl%dt

        ! Momentum
        if (spherical_coord_flag) then
            divv = dvx_dx + 2.0*vx/ptl%x
        else
            divv = dvx_dx
        endif
        deltap = -ptl%p * divv * ptl%dt / 3.0d0
        ! Momentum diffusion due to wave scattering
        if (dpp_wave_flag) then
            rho = fields(4)
            va = b / dsqrt(rho)
            if (momentum_dependency == 1) then
                deltap = deltap + (8*ptl%p / (27*kappa%kpara)) * va**2 * ptl%dt
            else
                deltap = deltap + (4*ptl%p / (9*kappa%kpara)) * va**2 * ptl%dt
            endif
            ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
            deltap = deltap + ran1 * va * ptl%p * dsqrt(2/(9*kappa%kpara)) * sdt
        endif

        ! Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            sigmaxx = dvx_dx - divv / 3
            sigmayy = -divv / 3
            sigmazz = -divv / 3
            if (weak_scattering) then
                bx = fields(5)
                by = fields(6)
                bz = fields(7)
                b = dsqrt(bx**2 + by**2 + bz**2)
                if (b == 0) then
                    ib = 0.0
                else
                    ib = 1.0 / b
                endif
                bbsigma = sigmaxx * bx**2 + sigmayy * by**2 + sigmazz * bz**2
                bbsigma = bbsigma * ib * ib
                gshear = bbsigma**2 / 5
            else
                gshear = 2 * (sigmaxx**2 + sigmayy**2 + sigmazz**2) / 15
            endif
            if (gshear > 0) then
                deltap = deltap + (2 + pindex) * gshear * tau0 * kappa%knorm0 * &
                    ptl%p**(pindex-1) * p0**(2.0-pindex) * ptl%dt
                if (.not. dpp_wave_flag) then
                    ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
                endif
                deltap = deltap + dsqrt(2 * gshear * tau0 * kappa%knorm0 * &
                    ptl%p**pindex * p0**(2.0-pindex)) * ran1 * sdt
            endif
        endif

        if (acc_region_flag == 1) then
            if (particle_in_acceleration_region(ptl)) then
                ptl%p = ptl%p + deltap
            else
                deltap = 0.0
            endif
        else
            ptl%p = ptl%p + deltap
        endif
        if (ptl%p < 0.5 * p0) then
            ptl%p = ptl%p - deltap
            deltap = 0.5 * p0 - ptl%p
            ptl%p = 0.5 * p0
        endif
    end subroutine push_particle_1d

    !---------------------------------------------------------------------------
    !< Push particle for a single step for a 1D simulation for focused transport
    !< Args:
    !<  thread_id: thread ID staring from 0
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  ptl: particle structure
    !<  fields: fields and their gradients at particle position
    !<  db2: turbulence variance and its gradients at particle position
    !<  lc: turbulence correlation length and its gradients at particle position
    !<  kappa: kappa and related variables
    !<  deltax: the change of x in this step
    !<  deltap: the change of p in this step
    !<  deltav: the change of v in this step
    !<  deltamu: the change of mu in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_1d_ft(thread_id, rt, ptl, fields, db2, &
            lc, kappa, deltax, deltap, deltav, deltamu)
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: interp_fields, &
            interp_magnetic_fluctuation, interp_correlation_length
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        real(dp), intent(in) :: rt
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        real(dp), dimension(*), intent(inout) :: db2, lc
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltap, deltav, deltamu
        real(dp) :: xtmp
        real(dp) :: sdt, dvx_dx, dvy_dx, dvz_dx
        real(dp) :: divv, bb_gradv, bv_gradv, gshear
        real(dp) :: b, ib, bx, by, bz, rt1
        real(dp) :: vx, vy, vz, cot_theta, ir
        real(dp) :: mu2, acc_rate, duu, dmu_dt, duu_du, dtmp
        real(dp) :: dp_dpp, div_bnorm, db_dx, h0, turb_gamma
        real(dp) :: xmin, xmax, dxm
        reaL(dp) :: xmin1, xmax1, dxmh
        real(dp) :: skpara, skpara1
        real(dp) :: ran1, sqrt3
        real(dp) :: rho, va ! Plasma density and Alfven speed
        real(dp) :: rands(2)
        real(dp) :: sigmaxx, sigmayy, sigmazz ! shear tensor
        real(dp) :: bbsigma ! b_ib_jsigma_ij
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights

        vx = fields(1)
        vy = fields(2)
        vz = fields(3)
        bx = fields(5)
        by = fields(6)
        bz = fields(7)
        b = dsqrt(bx**2 + by**2 + bz**2)
        dvx_dx = fields(nfields+1)
        dvy_dx = fields(nfields+4)
        dvz_dx = fields(nfields+7)
        xmin = fconfig%xmin
        xmax = fconfig%xmax
        dxm = mhd_config%dx
        dxmh = 0.5 * dxm

        !< The field data has two ghost cells, so the particles can cross the
        !< boundary without causing segment fault errors
        xmin1 = xmin - dxmh
        xmax1 = xmax + dxmh

        deltax = 0.0
        deltap = 0.0

        sdt = dsqrt(ptl%dt)
        if (spherical_coord_flag) then
            deltax = (vx + ptl%v * ptl%mu + kappa%dkxx_dx + 2.0*kappa%kxx/ptl%x) * ptl%dt
        else
            deltax = (vx + ptl%v * ptl%mu + kappa%dkxx_dx)*ptl%dt
        endif
        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltax = deltax + ran1*kappa%skperp*sdt

        ptl%x = ptl%x + deltax
        ptl%t = ptl%t + ptl%dt

        ! Momentum and velocity
        if (spherical_coord_flag) then
            ir = 1.0d0 / ptl%x
            cot_theta = 1.0d0 / tan(ptl%y)
            divv = dvx_dx + 2.0 * vx * ir
            bb_gradv = (bx*bx*dvx_dx - bx*(by*vy + bz*vz) * ir + &
                        by*bx*dvy_dx + by*(by*vx - cot_theta*bz*vz) * ir + &
                        bz*bx*dvz_dx + bz**2 * (vx + cot_theta*vy) * ir) / (b*b)
            bv_gradv = (bx*vx*dvx_dx - bx*(vy*vy + vz*vz) * ir + &
                        by*vx*dvy_dx + by*(vy*vx - cot_theta*vz*vz) * ir + &
                        bz*vx*dvz_dx + bz*vz * (vx + cot_theta*vy) * ir) / (b*b)
        else
            divv = dvx_dx
            bb_gradv = (bx*bx*dvx_dx + by*bx*dvy_dx + bz*bx*dvz_dx) / (b*b)
            bv_gradv = (bx*vx*dvx_dx + by*vx*dvy_dx + bz*vx*dvz_dx) / b
        endif
        mu2 = ptl%mu * ptl%mu
        acc_rate = -(0.5 *(1 - mu2) * divv + &
                     0.5 *(3 * mu2 - 1) * bb_gradv + &
                     ptl%mu * bv_gradv / ptl%v)
        deltap = ptl%p * acc_rate * ptl%dt
        deltav = ptl%v * acc_rate * ptl%dt

        ! Momentum diffusion due to wave scattering
        if (dpp_wave_flag) then
            rho = fields(4)
            va = b / dsqrt(rho)
            dp_dpp = 0.0d0
            if (momentum_dependency == 1) then
                dp_dpp = (8*ptl%p / (27*kappa%kpara)) * va**2 * ptl%dt
            else
                dp_dpp = (4*ptl%p / (9*kappa%kpara)) * va**2 * ptl%dt
            endif
            ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
            dp_dpp = dp_dpp + ran1 * va * ptl%p * dsqrt(2/(9*kappa%kpara)) * sdt
            deltap = deltap + dp_dpp
            deltav = deltav + ptl%v * dp_dpp / ptl%p
        endif

        ! Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            sigmaxx = dvx_dx - divv / 3
            sigmayy = -divv / 3
            sigmazz = -divv / 3
            if (weak_scattering) then
                if (b == 0) then
                    ib = 0.0
                else
                    ib = 1.0 / b
                endif
                bbsigma = sigmaxx * bx**2 + sigmayy * by**2 + sigmazz * bz**2
                bbsigma = bbsigma * ib * ib
                gshear = bbsigma**2 / 5
            else
                gshear = 2 * (sigmaxx**2 + sigmayy**2 + sigmazz**2) / 15
            endif
            dp_dpp = 0.0d0
            if (gshear > 0) then
                dp_dpp = (2 + pindex) * gshear * tau0 * kappa%knorm0 * &
                    ptl%p**(pindex-1) * p0**(2.0-pindex) * ptl%dt
                if (.not. dpp_wave_flag) then
                    ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
                endif
                dp_dpp = dsqrt(2 * gshear * tau0 * kappa%knorm0 * &
                    ptl%p**pindex * p0**(2.0-pindex)) * ran1 * sdt
            endif
            deltap = deltap + dp_dpp
            deltav = deltav + ptl%v * dp_dpp / ptl%p
        endif

        ! Pitch-angle evolution
        db_dx = fields(nfields+22)
        div_bnorm = -bx * db_dx / (b*b)
        dmu_dt = ptl%v * div_bnorm + ptl%mu * divv - &
            3 * ptl%mu * bb_gradv - 2 * bv_gradv / ptl%v
        dmu_dt = dmu_dt * (1-mu2) * 0.5
        h0 = 0.2  ! Enabling pitch-angle scattering around mu=0
        turb_gamma = 3.0 - pindex  ! Turbulence spectral index
        dtmp = abs(ptl%mu)**(turb_gamma-1) + h0
        duu = duu0 * (1-mu2) * dtmp
        if (ptl%mu .gt. 0.0d0) then
            duu_du = duu0 * (-2*ptl%mu * dtmp + &
                             (1-mu2) * abs(ptl%mu)**(turb_gamma-2))
        else if (ptl%mu .lt. 0.0d0) then
            duu_du = duu0 * (-2*ptl%mu * dtmp - &
                             (1-mu2) * abs(ptl%mu)**(turb_gamma-2))
        else
            duu_du = 0.0d0
        endif

        if (momentum_dependency == 1) then
            dtmp = (ptl%p / p0)**(turb_gamma-1)
            duu = duu * dtmp
            duu_du = duu_du * dtmp
        endif
        deltamu = (dmu_dt + duu_du) * ptl%dt
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltamu = deltamu + ran1*dsqrt(2.0*duu)*sdt
        ptl%mu = ptl%mu + deltamu

        if (acc_region_flag == 1) then
            if (particle_in_acceleration_region(ptl)) then
                ptl%p = ptl%p + deltap
                ptl%v = ptl%v + deltav
            else
                deltap = 0.0
                deltav = 0.0
            endif
        else
            ptl%p = ptl%p + deltap
            ptl%v = ptl%v + deltav
        endif

        if (ptl%p < 0.5 * p0) then
            ptl%v = ptl%v - deltav
            deltav = ptl%v * 0.5 * p0 / ptl%p - ptl%v
            ptl%v = ptl%v + deltav
            ptl%p = ptl%p - deltap
            deltap = 0.5 * p0 - ptl%p
            ptl%p = 0.5 * p0
        endif
    end subroutine push_particle_1d_ft

    !---------------------------------------------------------------------------
    !< Push particle for a single step for a 2D simulation
    !< Args:
    !<  thread_id: thread ID staring from 0
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  ptl: particle structure
    !<  fields: fields and their gradients at particle position
    !<  db2: turbulence variance and its gradients at particle position
    !<  lc: turbulence correlation length and its gradients at particle position
    !<  kappa: kappa and related variables
    !<  deltax, deltay, deltap: the change of x, y and p in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_2d(thread_id, rt, ptl, fields, db2, &
            lc, kappa, deltax, deltay, deltap)
        use constants, only: pi
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: interp_fields, &
            interp_magnetic_fluctuation, interp_correlation_length
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        real(dp), intent(in) :: rt
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        real(dp), dimension(*), intent(inout) :: db2, lc
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltay, deltap
        real(dp) :: sdt, dvx_dx, dvx_dy, dvy_dx, dvy_dy, divv, gshear
        real(dp) :: bx, by, bz, b, vx, vy, px, py, pz, rt1, ib
        real(dp) :: bx1, by1, btot1, ib1
        real(dp) :: dbx_dy, dby_dx, db_dx, db_dy, dbz_dx, dbz_dy
        real(dp) :: vdx, vdy, vdz, vdp, ib2, ib3, gbr, gbt
        real(dp) :: xmin, ymin, xmax, ymax, dxm, dym
        reaL(dp) :: xmin1, ymin1, xmax1, ymax1, dxmh, dymh
        real(dp) :: skperp, skpara_perp, skperp1, skpara_perp1
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rho, va ! Plasma density and Alfven speed
        real(dp) :: vflow ! Alfven speed + background flow
        real(dp) :: rands(2)
        real(dp) :: a1, b1, c1, Qpp, Qpm, Qmp, Qmm
        real(dp) :: ctheta, istheta, cot_theta, ir, ir2
        real(dp) :: qtmp1, qtmp2, atmp
        real(dp) :: a1_1, b1_1, c1_1, Qpp_1, Qpm_1, Qmp_1, Qmm_1
        real(dp) :: qtmp1_1, qtmp2_1
        real(dp) :: deltaz
        real(dp) :: sigmaxx, sigmayy, sigmazz, sigmaxy ! shear tensor
        real(dp) :: bbsigma ! b_ib_jsigma_ij
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights

        vx = fields(1)
        vy = fields(2)
        bx = fields(5)
        by = fields(6)
        bz = fields(7)
        b = dsqrt(bx**2 + by**2 + bz**2)
        dvx_dx = fields(nfields+1)
        dvy_dy = fields(nfields+5)
        xmin = fconfig%xmin
        ymin = fconfig%ymin
        xmax = fconfig%xmax
        ymax = fconfig%ymax
        dxm = mhd_config%dx
        dym = mhd_config%dy
        dxmh = 0.5 * dxm
        dymh = 0.5 * dym

        !< The field data has two ghost cells, so the particles can cross the
        !< boundary without causing segment fault errors
        xmin1 = xmin - dxmh
        ymin1 = ymin - dymh
        xmax1 = xmax + dxmh
        ymax1 = ymax + dymh

        deltax = 0.0
        deltay = 0.0
        deltap = 0.0

        sdt = dsqrt(ptl%dt)
        if (b == 0) then
            ib = 0.0
        else
            ib = 1.0 / b
        endif

        if (spherical_coord_flag) then
            ctheta = cos(ptl%y)
            istheta = 1.0 / sin(ptl%y)
            cot_theta = ctheta * istheta
            ir = 1.0 / ptl%x
            ir2 = 1.0 / ptl%x**2
            a1 = kappa%kxx
            b1 = kappa%kxy * ir
            c1 = kappa%kyy * ir2
            atmp = dsqrt((a1-c1)**2 + 4*b1**2)
            Qpp = atmp + (a1 + c1)
            Qmp = atmp - (a1 + c1)
            Qpm = atmp + (a1 - c1)
            Qmm = atmp - (a1 - c1)
            qtmp1 = dsqrt(-Qmp/(Qmm**2+4*b1**2))
            qtmp2 = dsqrt(Qpp/(Qpm**2+4*b1**2))
        endif

        !< When Bz is non-zero, there are in-plane drifts
        dbz_dx = fields(nfields+19)
        dbz_dy = fields(nfields+20)
        db_dx = fields(nfields+22)
        db_dy = fields(nfields+23)
        ib2 = ib * ib
        ib3 = ib * ib2
        vdp = 1.0 / (3 * pcharge) / dsqrt((drift1*p0/ptl%p)**2 + (drift2*p0**2/ptl%p**2)**2)
        if (spherical_coord_flag) then
            gbr = db_dx
            gbt = db_dy * ir
            vdx = vdp * ((dbz_dy + cot_theta * bz) * ir * ib2 - &
                         gbt * bz * 2.0 * ib3)
            vdy = vdp * (-(dbz_dx + bz * ir) * ib2 + gbr * bz * 2.0 * ib3)
            deltax = vdx * ptl%dt
            deltay = vdy * ir * ptl%dt
        else
            vdx = vdp * (dbz_dy * ib2 - 2.0 * bz * db_dy * ib3)
            vdy = vdp * (-dbz_dx * ib2 + 2.0 * bz * db_dx * ib3)
            deltax = vdx * ptl%dt
            deltay = vdy * ptl%dt
        endif

        !< Check particle drift along the out-of-plane direction
        if (check_drift_2d) then
            dbx_dy = fields(nfields+14)
            dby_dx = fields(nfields+16)
            if (spherical_coord_flag) then
                vdz = vdp * ((dby_dx + by*ir - dbx_dy*ir)*ib2 - &
                             (gbr*by - gbt*bx*ir)*2.0*ib3)
                deltaz = vdz * istheta * ir * ptl%dt
            else
                vdz = vdp * ((dby_dx-dbx_dy)*ib2 - 2*(by*db_dx-bx*db_dy)*ib3)
                deltaz = vdz * ptl%dt
            endif
        else
            deltaz = 0.0
        endif

        if (spherical_coord_flag) then
            deltax = deltax + (vx + kappa%dkxx_dx + &
                (2.0*kappa%kxx + kappa%dkxy_dy + kappa%kxy*cot_theta)*ir)*ptl%dt
            deltay = deltay + ((vy + kappa%dkxy_dx)*ir + &
                (kappa%kxy + kappa%dkyy_dy + kappa%kyy*cot_theta)*ir2)*ptl%dt
        else
            deltax = deltax + (vx + kappa%dkxx_dx + kappa%dkxy_dy)*ptl%dt
            deltay = deltay + (vy + kappa%dkxy_dx + kappa%dkyy_dy)*ptl%dt
        endif
        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran2 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran3 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3

        ! rands = two_normals(thread_id)
        ! ran1 = rands(1)
        ! ran2 = rands(2)
        ! rands = two_normals(thread_id)
        ! ran3 = rands(1)

        !< First-order method
        if (spherical_coord_flag) then
            deltax = deltax + (-Qmm*qtmp1*ran1 + Qpm*qtmp2*ran2)*sdt
            deltay = deltay + (2*b1*qtmp1*ran1 + 2*b1*qtmp2*ran2)*sdt
        else
            deltax = deltax + ran1*kappa%skperp*sdt + ran3*kappa%skpara_perp*sdt*bx*ib
            deltay = deltay + ran2*kappa%skperp*sdt + ran3*kappa%skpara_perp*sdt*by*ib
        endif

        ptl%x = ptl%x + deltax
        ptl%y = ptl%y + deltay
        ptl%z = ptl%z + deltaz
        ptl%t = ptl%t + ptl%dt

        ! Momentum
        if (spherical_coord_flag) then
            divv = dvx_dx + (2.0*vx + dvy_dy + vy*ctheta*istheta)*ir
        else
            divv = dvx_dx + dvy_dy
        endif
        deltap = -ptl%p * divv * ptl%dt / 3.0d0
        ! Momentum diffusion due to wave scattering
        if (dpp_wave_flag) then
            rho = fields(4)
            va = b / dsqrt(rho)
            vx = vx + sign(va * bx * ib, vx)
            vy = vy + sign(va * by * ib, vy)
            vflow = dsqrt(vx**2 + vy**2)
            if (momentum_dependency == 1) then
                deltap = deltap + (8*ptl%p / (27*kappa%kpara)) * va**2 * ptl%dt
            else
                deltap = deltap + (4*ptl%p / (9*kappa%kpara)) * va**2 * ptl%dt
            endif
            ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
            ! rands = two_normals(thread_id)
            ! ran1 = rands(1)
            deltap = deltap + ran1 * va * ptl%p * dsqrt(2/(9*kappa%kpara)) * sdt
        endif

        ! Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            dvx_dy = fields(nfields+2)
            dvy_dx = fields(nfields+4)
            sigmaxx = dvx_dx - divv / 3
            sigmayy = dvy_dy - divv / 3
            sigmazz = -divv / 3
            sigmaxy = (dvx_dy + dvy_dx) / 2
            if (weak_scattering) then
                bbsigma = sigmaxx * bx**2 + sigmayy * by**2 + &
                    sigmazz * bz**2 + 2.0 * sigmaxy * bx * by
                bbsigma = bbsigma * ib * ib
                gshear = bbsigma**2 / 5
            else
                gshear = 2 * (sigmaxx**2 + sigmayy**2 + sigmazz**2 + 2 * sigmaxy**2) / 15
            endif
            if (gshear > 0) then
                deltap = deltap + (2 + pindex) * gshear * tau0 * kappa%knorm0 * &
                    ptl%p**(pindex-1) * p0**(2.0-pindex) * ptl%dt
                if (.not. dpp_wave_flag) then
                    ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
                endif
                deltap = deltap + dsqrt(2 * gshear * tau0 * kappa%knorm0 * &
                    ptl%p**pindex * p0**(2.0-pindex)) * ran1 * sdt
            endif
        endif

        if (acc_region_flag == 1) then
            if (particle_in_acceleration_region(ptl)) then
                ptl%p = ptl%p + deltap
            else
                deltap = 0.0
            endif
        else
            ptl%p = ptl%p + deltap
        endif
        if (ptl%p < 0.5 * p0) then
            ptl%p = ptl%p - deltap
            deltap = 0.5 * p0 - ptl%p
            ptl%p = 0.5 * p0
        endif
    end subroutine push_particle_2d

    !---------------------------------------------------------------------------
    !< Push particle for a single step for a 2D simulation for focused transport
    !< Args:
    !<  thread_id: thread ID staring from 0
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  ptl: particle structure
    !<  fields: fields and their gradients at particle position
    !<  db2: turbulence variance and its gradients at particle position
    !<  lc: turbulence correlation length and its gradients at particle position
    !<  kappa: kappa and related variables
    !<  deltax: the change of x in this step
    !<  deltay: the change of y in this step
    !<  deltap: the change of p in this step
    !<  deltav: the change of v in this step
    !<  deltamu: the change of mu in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_2d_ft(thread_id, rt, ptl, fields, db2, &
            lc, kappa, deltax, deltay, deltap, deltav, deltamu)
        use constants, only: pi
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: interp_fields, &
            interp_magnetic_fluctuation, interp_correlation_length
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        real(dp), intent(in) :: rt
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        real(dp), dimension(*), intent(inout) :: db2, lc
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltay, deltap
        real(dp), intent(out) :: deltav, deltamu
        real(dp) :: xtmp, ytmp
        real(dp) :: sdt, divv, bb_gradv, bv_gradv, gshear
        real(dp) :: dvx_dx, dvx_dy, dvy_dx, dvy_dy, dvz_dx, dvz_dy
        real(dp) :: bx, by, bz, b, ib, rt1, bxn, byn, bzn, ibxyn
        real(dp) :: vx, vy, vz, vbx, vby
        real(dp) :: bx1, by1, btot1, ib1
        real(dp) :: dbx_dx, dbx_dy, dby_dx, dby_dy, dbz_dx, dbz_dy, db_dx, db_dy
        real(dp) :: vdx, vdy, vdz, vdp, ib2, ib3, gbr, gbt
        real(dp) :: xmin, ymin, xmax, ymax, dxm, dym
        reaL(dp) :: xmin1, ymin1, xmax1, ymax1, dxmh, dymh
        real(dp) :: skperp, skpara_perp, skperp1, skpara_perp1
        real(dp) :: mu2, muf1, muf2, kx, ky, kz, bdot_curvb, acc_rate
        real(dp) :: dp_dpp, div_bnorm, h0, turb_gamma
        real(dp) :: duu, dmu_dt, duu_du, dtmp
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rho, va ! Plasma density and Alfven speed
        real(dp) :: vflow ! Alfven speed + background flow
        real(dp) :: rands(2)
        real(dp) :: a1, b1, c1, Qpp, Qpm, Qmp, Qmm
        real(dp) :: ctheta, istheta, cot_theta, ir, ir2
        real(dp) :: qtmp1, qtmp2, atmp
        real(dp) :: a1_1, b1_1, c1_1, Qpp_1, Qpm_1, Qmp_1, Qmm_1
        real(dp) :: qtmp1_1, qtmp2_1
        real(dp) :: deltaz
        real(dp) :: sigmaxx, sigmayy, sigmazz, sigmaxy ! shear tensor
        real(dp) :: bbsigma ! b_ib_jsigma_ij
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights

        vx = fields(1)
        vy = fields(2)
        vz = fields(3)
        bx = fields(5)
        by = fields(6)
        bz = fields(7)
        b = dsqrt(bx**2 + by**2 + bz**2)
        dvx_dx = fields(nfields+1)
        dvx_dy = fields(nfields+2)
        dvy_dx = fields(nfields+4)
        dvy_dy = fields(nfields+5)
        dvz_dx = fields(nfields+7)
        dvz_dy = fields(nfields+8)
        xmin = fconfig%xmin
        ymin = fconfig%ymin
        xmax = fconfig%xmax
        ymax = fconfig%ymax
        dxm = mhd_config%dx
        dym = mhd_config%dy
        dxmh = 0.5 * dxm
        dymh = 0.5 * dym

        !< The field data has two ghost cells, so the particles can cross the
        !< boundary without causing segment fault errors
        xmin1 = xmin - dxmh
        ymin1 = ymin - dymh
        xmax1 = xmax + dxmh
        ymax1 = ymax + dymh

        deltax = 0.0
        deltay = 0.0
        deltap = 0.0
        deltav = 0.0
        deltamu = 0.0

        sdt = dsqrt(ptl%dt)
        if (b == 0) then
            ib = 0.0
        else
            ib = 1.0 / b
        endif

        if (spherical_coord_flag) then
            ctheta = cos(ptl%y)
            istheta = 1.0 / sin(ptl%y)
            cot_theta = ctheta * istheta
            ir = 1.0 / ptl%x
            ir2 = 1.0 / ptl%x**2
            a1 = kappa%kxx
            b1 = kappa%kxy * ir
            c1 = kappa%kyy * ir2
            atmp = dsqrt((a1-c1)**2 + 4*b1**2)
            Qpp = atmp + (a1 + c1)
            Qmp = atmp - (a1 + c1)
            Qpm = atmp + (a1 - c1)
            Qmm = atmp - (a1 - c1)
            qtmp1 = dsqrt(-Qmp/(Qmm**2+4*b1**2))
            qtmp2 = dsqrt(Qpp/(Qpm**2+4*b1**2))
        endif

        !< When Bz is non-zero, there are in-plane drifts
        dbx_dx = fields(nfields+13)
        dbx_dy = fields(nfields+14)
        dby_dx = fields(nfields+16)
        dby_dy = fields(nfields+17)
        dbz_dx = fields(nfields+19)
        dbz_dy = fields(nfields+20)
        db_dx = fields(nfields+22)
        db_dy = fields(nfields+23)
        ib2 = ib * ib
        ib3 = ib * ib2
        vdp = 1.0 / pcharge / dsqrt((drift1*p0/ptl%p)**2 + (drift2*p0**2/ptl%p**2)**2)
        mu2 = ptl%mu**2
        muf1 = 0.5 * (1.0 - mu2)
        muf2 = 0.5 * (3.0*mu2 - 1.0)
        if (spherical_coord_flag) then
            gbr = db_dx
            gbt = db_dy * ir
            kx = bx * dbx_dx + by * dbx_dy * ir - (by**2 + bz**2) * ir
            ky = bx * dby_dx + by * dby_dy * ir + (bx*by - cot_theta*bz**2) * ir
            kz = bx * dbz_dx + by * dbz_dy * ir + (bx*bz + cot_theta*by*bz) * ir
            bdot_curvb = bx * (dbz_dy + cot_theta * bz) * ir + &
                         by * (-dbz_dx - bz * ir) + &
                         bz * (dby_dx + by * ir - dbx_dy * ir)
            vdx = vdp * (muf1 * (-bz * gbt) * ib2 + &
                         mu2 * (by * kz - bz * ky) * ib3 + &
                         muf1 * bx * bdot_curvb * ib3)
            vdy = vdp * (muf1 * (bz * gbr) * ib2 + &
                         mu2 * (bz * kx - bx * kz) * ib3 + &
                         muf1 * by * bdot_curvb * ib3)
            deltax = vdx * ptl%dt
            deltay = vdy * ir * ptl%dt
        else
            kx = bx * dbx_dx + by * dbx_dy
            ky = bx * dby_dx + by * dby_dy
            kz = bx * dbz_dx + by * dbz_dy
            bdot_curvb = bx * dbz_dy - by * dbz_dx + bz * (dby_dx - dbx_dy)
            vdx = vdp * (muf1 * (-bz * db_dy) * ib2 + &
                         mu2 * (by * kz - bz * ky) * ib3 + &
                         muf1 * bx * bdot_curvb * ib3)
            vdy = vdp * (muf1 * ( bz * db_dx) * ib2 + &
                         mu2 * (bz * kx - bx * kz) * ib3 + &
                         muf1 * by * bdot_curvb * ib3)
            deltax = vdx * ptl%dt
            deltay = vdy * ptl%dt
        endif

        !< Check particle drift along the out-of-plane direction
        if (check_drift_2d) then
            if (spherical_coord_flag) then
                vdz = vdp * (muf1 * (bx * gbt - by * gbr) * ib2 + &
                             mu2 * (bx * ky - by * kx) * ib3 + &
                             muf1 * bz * bdot_curvb * ib3)
                deltaz = vdz * istheta * ir * ptl%dt
            else
                vdz = vdp * (muf1 * (bx * db_dy - by * db_dx) * ib2 + &
                             mu2 * (bx * ky - by * kx) * ib3 + &
                             muf1 * bz * bdot_curvb * ib3)
                deltaz = vdz * ptl%dt
            endif
        else
            deltaz = 0.0
        endif

        vbx = ptl%v * ptl%mu * ib
        vby = vbx * by ! particle velocity along the magnetic field
        vbx = vbx * bx
        if (spherical_coord_flag) then
            deltax = deltax + (vx + vbx + kappa%dkxx_dx + &
                (2.0*kappa%kxx + kappa%dkxy_dy + kappa%kxy*cot_theta)*ir)*ptl%dt
            deltay = deltay + ((vy + vby + kappa%dkxy_dx)*ir + &
                (kappa%kxy + kappa%dkyy_dy + kappa%kyy*cot_theta)*ir2)*ptl%dt
        else
            deltax = deltax + (vx + vbx + kappa%dkxx_dx + kappa%dkxy_dy)*ptl%dt
            deltay = deltay + (vy + vby + kappa%dkxy_dx + kappa%dkyy_dy)*ptl%dt
        endif
        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran2 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3

        ! rands = two_normals(thread_id)
        ! ran1 = rands(1)
        ! ran2 = rands(2)
        ! rands = two_normals(thread_id)
        ! ran3 = rands(1)

        !< First-order method
        if (spherical_coord_flag) then
            deltax = deltax + (-Qmm*qtmp1*ran1 + Qpm*qtmp2*ran2)*sdt
            deltay = deltay + (2*b1*qtmp1*ran1 + 2*b1*qtmp2*ran2)*sdt
        else
            bxn = bx * ib
            byn = by * ib
            bzn = bz * ib
            ibxyn = 1.0d0 / dsqrt(bxn**2 + byn**2)
            deltax = deltax + kappa%skperp * ibxyn * sdt * (-bxn*bzn*ran1 - by*ran2)
            deltay = deltay + kappa%skperp * ibxyn * sdt * (-byn*bzn*ran1 + bx*ran2)
        endif

        ptl%x = ptl%x + deltax
        ptl%y = ptl%y + deltay
        ptl%z = ptl%z + deltaz
        ptl%t = ptl%t + ptl%dt

        ! Momentum and velocity
        if (spherical_coord_flag) then
            divv = dvx_dx + (2.0*vx + dvy_dy + vy*cot_theta)*ir
            bb_gradv = (bx*bx*dvx_dx + bx*by*dvx_dy*ir - bx*(by*vy + bz*vz)*ir + &
                        by*bx*dvy_dx + by*by*dvy_dy*ir + by*(by*vx - cot_theta*bz*vz)*ir + &
                        bz*bx*dvz_dx + bz*by*dvz_dy*ir + bz**2*(vx + cot_theta*vy)*ir) * ib2
            bv_gradv = (bx*vx*dvx_dx + bx*vy*dvx_dy*ir - bx*(vy*vy + vz*vz)*ir + &
                        by*vx*dvy_dx + by*vy*dvy_dy*ir + by*(vy*vx - cot_theta*vz*vz)*ir + &
                        bz*vx*dvz_dx + bz*vy*dvz_dy*ir + bz*vz*(vx + cot_theta*vy)*ir) * ib
        else
            divv = dvx_dx + dvy_dy
            bb_gradv = (bx * (bx * dvx_dx + by * dvx_dy) + &
                        by * (bx * dvy_dx + by * dvy_dy) + &
                        bz * (bx * dvz_dx + by * dvz_dy)) * ib2
            bv_gradv = (bx * (vx * dvx_dx + vy * dvx_dy) + &
                        by * (vx * dvy_dx + vy * dvy_dy) + &
                        bz * (vx * dvz_dx + vy * dvz_dy)) * ib
        endif
        acc_rate = -(muf1 * divv + muf2 * bb_gradv + ptl%mu * bv_gradv / ptl%v)
        deltap = ptl%p * acc_rate * ptl%dt
        deltav = ptl%v * acc_rate * ptl%dt

        ! Momentum diffusion due to wave scattering
        if (dpp_wave_flag) then
            rho = fields(4)
            va = b / dsqrt(rho)
            vx = vx + sign(va * bx * ib, vx)
            vy = vy + sign(va * by * ib, vy)
            vflow = dsqrt(vx**2 + vy**2)
            dp_dpp = 0.0d0
            if (momentum_dependency == 1) then
                dp_dpp = (8*ptl%p / (27*kappa%kpara)) * va**2 * ptl%dt
            else
                dp_dpp = (4*ptl%p / (9*kappa%kpara)) * va**2 * ptl%dt
            endif
            ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
            ! rands = two_normals(thread_id)
            ! ran1 = rands(1)
            dp_dpp = dp_dpp + ran1 * va * ptl%p * dsqrt(2/(9*kappa%kpara)) * sdt
            deltap = deltap + dp_dpp
            deltav = deltav + ptl%v * dp_dpp / ptl%p
        endif

        ! Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            dvx_dy = fields(nfields+2)
            dvy_dx = fields(nfields+4)
            sigmaxx = dvx_dx - divv / 3
            sigmayy = dvy_dy - divv / 3
            sigmazz = -divv / 3
            sigmaxy = (dvx_dy + dvy_dx) / 2
            if (weak_scattering) then
                bbsigma = sigmaxx * bx**2 + sigmayy * by**2 + &
                    sigmazz * bz**2 + 2.0 * sigmaxy * bx * by
                bbsigma = bbsigma * ib * ib
                gshear = bbsigma**2 / 5
            else
                gshear = 2 * (sigmaxx**2 + sigmayy**2 + sigmazz**2 + 2 * sigmaxy**2) / 15
            endif
            dp_dpp = 0.0d0
            if (gshear > 0) then
                dp_dpp = (2 + pindex) * gshear * tau0 * kappa%knorm0 * &
                    ptl%p**(pindex-1) * p0**(2.0-pindex) * ptl%dt
                if (.not. dpp_wave_flag) then
                    ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
                endif
                dp_dpp = dsqrt(2 * gshear * tau0 * kappa%knorm0 * &
                    ptl%p**pindex * p0**(2.0-pindex)) * ran1 * sdt
            endif
            deltap = deltap + dp_dpp
            deltav = deltav + ptl%v * dp_dpp / ptl%p
        endif

        ! Pitch-angle evolution
        div_bnorm = -(bx * db_dx + by * db_dy * ir) * ib2
        dmu_dt = ptl%v * div_bnorm + ptl%mu * divv - &
            3 * ptl%mu * bb_gradv - 2 * bv_gradv / ptl%v
        dmu_dt = dmu_dt * muf1
        h0 = 0.2  ! Enabling pitch-angle scattering around mu=0
        turb_gamma = 3.0 - pindex  ! Turbulence spectral index
        dtmp = abs(ptl%mu)**(turb_gamma-1) + h0
        duu = duu0 * (1-mu2) * dtmp
        if (ptl%mu .gt. 0.0d0) then
            duu_du = duu0 * (-2*ptl%mu * dtmp + &
                             (1-mu2) * abs(ptl%mu)**(turb_gamma-2))
        else if (ptl%mu .lt. 0.0d0) then
            duu_du = duu0 * (-2*ptl%mu * dtmp - &
                             (1-mu2) * abs(ptl%mu)**(turb_gamma-2))
        else
            duu_du = 0.0d0
        endif

        if (momentum_dependency == 1) then
            dtmp = (ptl%p / p0)**(turb_gamma-1)
            duu = duu * dtmp
            duu_du = duu_du * dtmp
        endif
        deltamu = (dmu_dt + duu_du) * ptl%dt
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltamu = deltamu + ran1*dsqrt(2.0*duu)*sdt
        ptl%mu = ptl%mu + deltamu

        if (acc_region_flag == 1) then
            if (particle_in_acceleration_region(ptl)) then
                ptl%p = ptl%p + deltap
                ptl%v = ptl%v + deltav
            else
                deltap = 0.0
                deltav = 0.0
            endif
        else
            ptl%p = ptl%p + deltap
            ptl%v = ptl%v + deltav
        endif

        if (ptl%p < 0.5 * p0) then
            ptl%v = ptl%v - deltav
            deltav = ptl%v * 0.5 * p0 / ptl%p - ptl%v
            ptl%v = ptl%v + deltav
            ptl%p = ptl%p - deltap
            deltap = 0.5 * p0 - ptl%p
            ptl%p = 0.5 * p0
        endif
    end subroutine push_particle_2d_ft

    !---------------------------------------------------------------------------
    !< Push particle for a single step for a 2D simulation but include transport
    !< along the 3rd dimension
    !< Args:
    !<  thread_id: thread ID staring from 0
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  ptl: particle structure
    !<  fields: fields and their gradients at particle position
    !<  db2: turbulence variance and its gradients at particle position
    !<  lc: turbulence correlation length and its gradients at particle position
    !<  kappa: kappa and related variables
    !<  deltax, deltay, deltaz, deltap: the change of x, y, z and p in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_2d_include_3rd(thread_id, rt, ptl, fields, db2, &
            lc, kappa, deltax, deltay, deltaz, deltap)
        use constants, only: pi
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: interp_fields, &
            interp_magnetic_fluctuation, interp_correlation_length
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        real(dp), intent(in) :: rt
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        real(dp), dimension(*), intent(inout) :: db2, lc
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltay, deltaz, deltap
        real(dp) :: rt1, sdt
        real(dp) :: dvx_dx, dvx_dy, dvx_dz
        real(dp) :: dvy_dx, dvy_dy, dvy_dz
        real(dp) :: dvz_dx, dvz_dy, dvz_dz
        real(dp) :: divv, gshear
        real(dp) :: vx, vy, vz
        real(dp) :: bx, by, bz, b, ib, ib2, ib3
        real(dp) :: bxn, byn, bzn, bxyn, ibxyn
        real(dp) :: px, py, p
        real(dp) :: dbx_dy, dbx_dz, dby_dx, dby_dz, dbz_dx, dbz_dy
        real(dp) :: db_dx, db_dy, db_dz
        real(dp) :: vdx, vdy, vdz, vdp
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rho, va ! Plasma density and Alfven speed
        real(dp) :: rands(2)
        real(dp) :: ctheta, istheta, ir, ir2
        real(dp) :: gbr, gbt, gbp, p11, p12, p13, p22, p23, p33
        real(dp) :: sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, sigmayz ! shear tensor
        real(dp) :: bbsigma ! b_ib_jsigma_ij
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights

        vx = fields(1)
        vy = fields(2)
        vz = fields(3)
        bx = fields(5)
        by = fields(6)
        bz = fields(7)
        b = dsqrt(bx**2 + by**2 + bz**2)
        if (b == 0.0d0) then
            ib = 0.0_dp
        else
            ib = 1.0_dp / b
        endif
        bxn = bx * ib
        byn = by * ib
        bzn = bz * ib
        bxyn = dsqrt(bxn**2 + byn**2)
        if (bxyn == 0.0d0) then
            ibxyn = 0.0d0
        else
            ibxyn = 1.0_dp / bxyn
        endif
        dvx_dx = fields(nfields+1)
        dvy_dy = fields(nfields+5)
        dvz_dz = 0.0_dp

        if (spherical_coord_flag) then
            ctheta = cos(ptl%y)
            istheta = 1.0 / sin(ptl%y)
            ir = 1.0 / ptl%x
            ir2 = 1.0 / ptl%x**2
            p11 = dsqrt(2.0*(kappa%kxx*kappa%kyz**2 + &
                             kappa%kyy*kappa%kxz**2 + &
                             kappa%kzz*kappa%kxy**2 - &
                             2.0*kappa%kxy*kappa%kxz*kappa%kyz - &
                             kappa%kxx*kappa%kyy*kappa%kzz) / &
                             (kappa%kyz**2 - kappa%kyy*kappa%kzz))
            p12 = (kappa%kxz*kappa%kyz - kappa%kxy*kappa%kzz) * &
                dsqrt(2.0*(kappa%kyy - (kappa%kyz**2/kappa%kzz))) / &
                (kappa%kyz**2 - kappa%kyy*kappa%kzz)
            p13 = dsqrt(2.0 / kappa%kzz) * kappa%kxz
            p22 = dsqrt(2.0 * (kappa%kyy - kappa%kyz**2/kappa%kzz)) * ir
            p23 = dsqrt(2.0 / kappa%kzz) * kappa%kyz * ir
            p33 = dsqrt(2.0 * kappa%kzz) * istheta * ir
        endif

        ! Drift velocity
        dbx_dy = fields(nfields+14)
        dbx_dz = 0.0_dp
        dby_dx = fields(nfields+16)
        dby_dz = 0.0_dp
        dbz_dx = fields(nfields+19)
        dbz_dy = fields(nfields+20)
        db_dx = fields(nfields+22)
        db_dy = fields(nfields+23)
        db_dz = 0.0_dp
        ib2 = ib * ib
        ib3 = ib * ib2
        vdp = 1.0 / (3 * pcharge) / dsqrt((drift1*p0/ptl%p)**2 + (drift2*p0**2/ptl%p**2)**2)
        if (spherical_coord_flag) then
            gbr = db_dx
            gbt = db_dy * ir
            gbp = db_dz * ir * istheta
            vdx = vdp * ((dbz_dy + bz*ctheta*istheta - dby_dz*istheta)*ir*ib2 - &
                         (gbt*bz - gbp*by)*2.0*ib3)
            vdy = vdp * ((dbx_dz*istheta*ir - dbz_dx - bz*ir)*ib2 - &
                         (gbp*bx - gbr*bz)*2.0*ib3)
            vdz = vdp * ((dby_dx + by*ir - dbx_dy*ir)*ib2 - &
                         (gbr*by - gbt*bx)*2.0*ib3)
        else
            vdx = vdp * ((dbz_dy-dby_dz)*ib2 - 2*(bz*db_dy-by*db_dz)*ib3)
            vdy = vdp * ((dbx_dz-dbz_dx)*ib2 - 2*(bx*db_dz-bz*db_dx)*ib3)
            vdz = vdp * ((dby_dx-dbx_dy)*ib2 - 2*(by*db_dx-bx*db_dy)*ib3)
        endif

        kappa%dkxz_dz = 0.0_dp
        kappa%dkyz_dz = 0.0_dp
        kappa%dkzz_dz = 0.0_dp
        if (spherical_coord_flag) then
            deltax = (vx + vdx + kappa%dkxx_dx + &
                (2.0*kappa%kxx + kappa%dkxy_dy + &
                kappa%kxy*ctheta*istheta + kappa%dkxz_dz*istheta)*ir) * ptl%dt
            deltay = ((vy + vdy + kappa%dkxy_dx)*ir + &
                (kappa%kxy + kappa%dkyy_dy + kappa%kyy*ctheta*istheta + &
                kappa%dkyz_dz*istheta)*ir2) * ptl%dt
            deltaz = ((vz + vdz + kappa%dkxz_dx)*istheta*ir + &
                (kappa%kxz + kappa%dkyz_dy + kappa%dkzz_dz*istheta)*istheta*ir2) * ptl%dt
        else
            deltax = (vx + vdx + kappa%dkxx_dx + kappa%dkxy_dy + kappa%dkxz_dz) * ptl%dt
            deltay = (vy + vdy + kappa%dkxy_dx + kappa%dkyy_dy + kappa%dkyz_dz) * ptl%dt
            deltaz = (vz + vdz + kappa%dkxz_dx + kappa%dkyz_dy + kappa%dkzz_dz) * ptl%dt
        endif

        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran2 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran3 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3

        sdt = dsqrt(ptl%dt)

        if (spherical_coord_flag) then
            deltax = deltax + (p11*ran1 + p12*ran2 + p13*ran3)*sdt
            deltay = deltay + (p22*ran2 + p23*ran3)*sdt
            deltaz = deltaz + p33*ran3*sdt
        else
            deltax = deltax + &
                (bxn*kappa%skpara*ran1 - bxn*bzn*kappa%skperp*ibxyn*ran2 - &
                 byn*kappa%skperp*ibxyn*ran3)*sdt
            deltay = deltay + &
                (byn*kappa%skpara*ran1 - byn*bzn*kappa%skperp*ibxyn*ran2 + &
                 bxn*kappa%skperp*ibxyn*ran3)*sdt
            deltaz = deltaz + (bzn*kappa%skpara*ran1 + bxyn*kappa%skperp*ran2)*sdt
        endif

        ptl%x = ptl%x + deltax
        ptl%y = ptl%y + deltay
        ptl%z = ptl%z + deltaz
        ptl%t = ptl%t + ptl%dt

        ! Momentum
        if (spherical_coord_flag) then
            divv = dvx_dx + &
                (2.0*vx + dvy_dy + vy*ctheta*istheta + dvz_dz*istheta)*ir
        else
            divv = dvx_dx + dvy_dy + dvz_dz
        endif
        deltap = -ptl%p * divv * ptl%dt / 3.0d0
        ! Momentum diffusion due to wave scattering
        if (dpp_wave_flag) then
            rho = fields(4)
            va = b / dsqrt(rho)
            if (momentum_dependency == 1) then
                deltap = deltap + (8*ptl%p / (27*kappa%kpara)) * va**2 * ptl%dt
            else
                deltap = deltap + (4*ptl%p / (9*kappa%kpara)) * va**2 * ptl%dt
            endif
            ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
            deltap = deltap + ran1 * va * ptl%p * dsqrt(2/(9*kappa%kpara)) * sdt
        endif

        ! Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            dvx_dy = fields(nfields+2)
            dvx_dz = 0.0_dp
            dvy_dx = fields(nfields+4)
            dvy_dz = 0.0_dp
            dvz_dx = fields(nfields+7)
            dvz_dy = fields(nfields+8)
            sigmaxx = dvx_dx - divv / 3
            sigmayy = dvy_dy - divv / 3
            sigmazz = dvz_dz - divv / 3
            sigmaxy = (dvx_dy + dvy_dx) / 2
            sigmaxz = (dvx_dz + dvz_dx) / 2
            sigmayz = (dvy_dz + dvz_dy) / 2
            if (weak_scattering) then
                bbsigma = sigmaxx * bx**2 + sigmayy * by**2 + sigmazz * bz**2 + &
                    2.0 * (sigmaxy * bx * by + sigmaxz * bx * bz + sigmayz * by * bz)
                bbsigma = bbsigma * ib * ib
                gshear = bbsigma**2 / 5
            else
                gshear = 2 * (sigmaxx**2 + sigmayy**2 + sigmazz**2 + &
                    2 * (sigmaxy**2 + sigmaxz**2 + sigmayz**2)) / 15
            endif
            deltap = deltap + (2 + pindex) * gshear * tau0 * kappa%knorm0 * &
                ptl%p**(pindex-1) * p0**(2.0-pindex) * ptl%dt
            if (.not. dpp_wave_flag) then
                ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
            endif
            deltap = deltap + dsqrt(2 * gshear * tau0 * kappa%knorm0 * &
                ptl%p**pindex * p0**(2.0-pindex)) * ran1 * sdt
        endif

        if (acc_region_flag == 1) then
            if (particle_in_acceleration_region(ptl)) then
                ptl%p = ptl%p + deltap
            else
                deltap = 0.0
            endif
        else
            ptl%p = ptl%p + deltap
        endif
        if (ptl%p < 0.5 * p0) then
            ptl%p = ptl%p - deltap
            deltap = 0.5 * p0 - ptl%p
            ptl%p = 0.5 * p0
        endif
    end subroutine push_particle_2d_include_3rd

    !---------------------------------------------------------------------------
    !< Push particle for a single step for a 2D simulation but include transport
    !< along the 3rd dimension when solving focused transport equation
    !< Args:
    !<  thread_id: thread ID staring from 0
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  ptl: particle structure
    !<  fields: fields and their gradients at particle position
    !<  db2: turbulence variance and its gradients at particle position
    !<  lc: turbulence correlation length and its gradients at particle position
    !<  kappa: kappa and related variables
    !<  deltax: the change of x in this step
    !<  deltay: the change of y in this step
    !<  deltaz: the change of z in this step
    !<  deltap: the change of p in this step
    !<  deltav: the change of v in this step
    !<  deltamu: the change of mu in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_2d_include_3rd_ft(thread_id, rt, ptl, fields, db2, &
            lc, kappa, deltax, deltay, deltaz, deltap, deltav, deltamu)
        use constants, only: pi
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: interp_fields, &
            interp_magnetic_fluctuation, interp_correlation_length
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        real(dp), intent(in) :: rt
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        real(dp), dimension(*), intent(inout) :: db2, lc
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltay, deltaz, deltap, deltav, deltamu
        real(dp) :: rt1, sdt
        real(dp) :: dvx_dx, dvx_dy
        real(dp) :: dvy_dx, dvy_dy
        real(dp) :: dvz_dx, dvz_dy
        real(dp) :: divv, bb_gradv, bv_gradv, gshear
        real(dp) :: vx, vy, vz, vbx, vby, vbz
        real(dp) :: bx, by, bz, b, ib, ib2, ib3
        real(dp) :: bxn, byn, bzn, bxyn, ibxyn
        real(dp) :: px, py, p
        real(dp) :: dbx_dx, dbx_dy
        real(dp) :: dby_dx, dby_dy
        real(dp) :: dbz_dx, dbz_dy
        real(dp) :: db_dx, db_dy
        real(dp) :: vdx, vdy, vdz, vdp
        real(dp) :: mu2, muf1, muf2, kx, ky, kz, bdot_curvb, acc_rate
        real(dp) :: dp_dpp, div_bnorm, h0, turb_gamma
        real(dp) :: duu, dmu_dt, duu_du, dtmp
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rho, va ! Plasma density and Alfven speed
        real(dp) :: rands(2)
        real(dp) :: ctheta, istheta, cot_theta, ir, ir2
        real(dp) :: gbr, gbt, p11, p12, p13, p22, p23, p33
        real(dp) :: sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, sigmayz ! shear tensor
        real(dp) :: bbsigma ! b_ib_jsigma_ij
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights

        vx = fields(1)
        vy = fields(2)
        vz = fields(3)
        bx = fields(5)
        by = fields(6)
        bz = fields(7)
        b = dsqrt(bx**2 + by**2 + bz**2)
        if (b == 0.0d0) then
            ib = 0.0_dp
        else
            ib = 1.0_dp / b
        endif
        bxn = bx * ib
        byn = by * ib
        bzn = bz * ib
        bxyn = dsqrt(bxn**2 + byn**2)
        if (bxyn == 0.0d0) then
            ibxyn = 0.0d0
        else
            ibxyn = 1.0_dp / bxyn
        endif

        if (spherical_coord_flag) then
            ctheta = cos(ptl%y)
            istheta = 1.0 / sin(ptl%y)
            ir = 1.0 / ptl%x
            ir2 = 1.0 / ptl%x**2
            p11 = dsqrt(2.0*(kappa%kxx*kappa%kyz**2 + &
                             kappa%kyy*kappa%kxz**2 + &
                             kappa%kzz*kappa%kxy**2 - &
                             2.0*kappa%kxy*kappa%kxz*kappa%kyz - &
                             kappa%kxx*kappa%kyy*kappa%kzz) / &
                             (kappa%kyz**2 - kappa%kyy*kappa%kzz))
            p12 = (kappa%kxz*kappa%kyz - kappa%kxy*kappa%kzz) * &
                dsqrt(2.0*(kappa%kyy - (kappa%kyz**2/kappa%kzz))) / &
                (kappa%kyz**2 - kappa%kyy*kappa%kzz)
            p13 = dsqrt(2.0 / kappa%kzz) * kappa%kxz
            p22 = dsqrt(2.0 * (kappa%kyy - kappa%kyz**2/kappa%kzz)) * ir
            p23 = dsqrt(2.0 / kappa%kzz) * kappa%kyz * ir
            p33 = dsqrt(2.0 * kappa%kzz) * istheta * ir
        endif

        ! Drift velocity
        dbx_dx = fields(nfields+13)
        dbx_dy = fields(nfields+14)
        dby_dx = fields(nfields+16)
        dby_dy = fields(nfields+17)
        dbz_dx = fields(nfields+19)
        dbz_dy = fields(nfields+20)
        db_dx = fields(nfields+22)
        db_dy = fields(nfields+23)
        ib2 = ib * ib
        ib3 = ib * ib2
        vdp = 1.0 / pcharge / dsqrt((drift1*p0/ptl%p)**2 + (drift2*p0**2/ptl%p**2)**2)
        mu2 = ptl%mu**2
        muf1 = 0.5 * (1.0 - mu2)
        muf2 = 0.5 * (3.0*mu2 - 1.0)
        if (spherical_coord_flag) then
            gbr = db_dx
            gbt = db_dy * ir
            kx = bx * dbx_dx + ir * (by * dbx_dy  - (by**2 + bz**2))
            ky = bx * dby_dx + ir * (by * dby_dy + (bx*by - cot_theta*bz**2))
            kz = bx * dbz_dx + ir * (by * dbz_dy + (bx*bz + cot_theta*by*bz))
            bdot_curvb = bx * (dbz_dy + cot_theta * bz) * ir + &
                         by * (-dbz_dx - bz * ir) + &
                         bz * (dby_dx + by * ir - dbx_dy * ir)
            vdx = vdp * (muf1 * (-bz * gbt) * ib2 + &
                         mu2 * (by * kz - bz * ky) * ib3 + &
                         muf1 * bx * bdot_curvb * ib3)
            vdy = vdp * (muf1 * (bz * gbr) * ib2 + &
                         mu2 * (bz * kx - bx * kz) * ib3 + &
                         muf1 * by * bdot_curvb * ib3)
            vdz = vdp * (muf1 * (bx * gbt - by * gbr) * ib2 + &
                         mu2 * (bx * ky - by * kx) * ib3 + &
                         muf1 * bz * bdot_curvb * ib3)
        else
            kx = bx * dbx_dx + by * dbx_dy
            ky = bx * dby_dx + by * dby_dy
            kz = bx * dbz_dx + by * dbz_dy
            bdot_curvb = bx * (dbz_dy) + &
                         by * (-dbz_dx) + &
                         bz * (dby_dx - dbx_dy)
            vdx = vdp * (muf1 * (-bz * db_dy) * ib2 + &
                         mu2 * (by * kz - bz * ky) * ib3 + &
                         muf1 * bx * bdot_curvb * ib3)
            vdy = vdp * (muf1 * (bz * db_dx) * ib2 + &
                         mu2 * (bz * kx - bx * kz) * ib3 + &
                         muf1 * by * bdot_curvb * ib3)
            vdz = vdp * (muf1 * (bx * db_dy - by * db_dx) * ib2 + &
                         mu2 * (bx * ky - by * kx) * ib3 + &
                         muf1 * bz * bdot_curvb * ib3)
        endif

        vbx = ptl%v * ptl%mu * ib
        vby = vbx * by ! particle velocity along the magnetic field
        vbz = vbx * bz
        vbx = vbx * bx
        if (spherical_coord_flag) then
            deltax = (vx + vbx + vdx + kappa%dkxx_dx + &
                (2.0*kappa%kxx + kappa%dkxy_dy + kappa%kxy*cot_theta)*ir) * ptl%dt
            deltay = ((vy + vby + vdy + kappa%dkxy_dx)*ir + &
                (kappa%kxy + kappa%dkyy_dy + kappa%kyy*cot_theta)*ir2) * ptl%dt
            deltaz = ((vz + vbz + vdz + kappa%dkxz_dx)*istheta*ir + &
                (kappa%kxz + kappa%dkyz_dy)*istheta*ir2) * ptl%dt
        else
            deltax = (vx + vbx + vdx + kappa%dkxx_dx + kappa%dkxy_dy) * ptl%dt
            deltay = (vy + vby + vdy + kappa%dkxy_dx + kappa%dkyy_dy) * ptl%dt
            deltaz = (vz + vbz + vdz + kappa%dkxz_dx + kappa%dkyz_dy) * ptl%dt
        endif

        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran2 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran3 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3

        sdt = dsqrt(ptl%dt)

        if (spherical_coord_flag) then
            deltax = deltax + (p11*ran1 + p12*ran2 + p13*ran3)*sdt
            deltay = deltay + (p22*ran2 + p23*ran3)*sdt
            deltaz = deltaz + p33*ran3*sdt
        else
            deltax = deltax + (-bxn*bzn*kappa%skperp*ibxyn*ran1 - &
                               byn*kappa%skperp*ibxyn*ran2)*sdt
            deltay = deltay + (-byn*bzn*kappa%skperp*ibxyn*ran1 + &
                               bxn*kappa%skperp*ibxyn*ran2)*sdt
            deltaz = deltaz + bxyn*kappa%skperp*ran1*sdt
        endif

        ptl%x = ptl%x + deltax
        ptl%y = ptl%y + deltay
        ptl%z = ptl%z + deltaz
        ptl%t = ptl%t + ptl%dt

        ! Momentum and velocity
        dvx_dx = fields(nfields+1)
        dvx_dy = fields(nfields+2)
        dvy_dx = fields(nfields+4)
        dvy_dy = fields(nfields+5)
        dvz_dx = fields(nfields+7)
        dvz_dy = fields(nfields+8)
        if (spherical_coord_flag) then
            divv = dvx_dx + (2.0*vx + dvy_dy + vy*cot_theta)*ir
            bb_gradv = (bx*bx*dvx_dx + bx*by*dvx_dy*ir - bx*(by*vy + bz*vz)*ir + &
                        by*bx*dvy_dx + by*by*dvy_dy*ir + by*(by*vx - cot_theta*bz*vz)*ir + &
                        bz*bx*dvz_dx + bz*by*dvz_dy*ir + bz**2*(vx + cot_theta*vy)*ir) * ib2
            bv_gradv = (bx*vx*dvx_dx + bx*vy*dvx_dy*ir - bx*(vy*vy + vz*vz)*ir + &
                        by*vx*dvy_dx + by*vy*dvy_dy*ir + by*(vy*vx - cot_theta*vz*vz)*ir + &
                        bz*vx*dvz_dx + bz*vy*dvz_dy*ir + bz*vz*(vx + cot_theta*vy)*ir) * ib
        else
            divv = dvx_dx + dvy_dy
            bb_gradv = (bx * (bx * dvx_dx + by * dvx_dy) + &
                        by * (bx * dvy_dx + by * dvy_dy) + &
                        bz * (bx * dvz_dx + by * dvz_dy)) * ib2
            bv_gradv = (bx * (vx * dvx_dx + vy * dvx_dy) + &
                        by * (vx * dvy_dx + vy * dvy_dy) + &
                        bz * (vx * dvz_dx + vy * dvz_dy)) * ib
        endif
        acc_rate = -(muf1 * divv + muf2 * bb_gradv + ptl%mu * bv_gradv / ptl%v)
        deltap = ptl%p * acc_rate * ptl%dt
        deltav = ptl%v * acc_rate * ptl%dt

        ! Momentum diffusion due to wave scattering
        if (dpp_wave_flag) then
            rho = fields(4)
            va = b / dsqrt(rho)
            if (momentum_dependency == 1) then
                dp_dpp = (8*ptl%p / (27*kappa%kpara)) * va**2 * ptl%dt
            else
                dp_dpp = (4*ptl%p / (9*kappa%kpara)) * va**2 * ptl%dt
            endif
            ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
            dp_dpp = dp_dpp + ran1 * va * ptl%p * dsqrt(2/(9*kappa%kpara)) * sdt
            deltap = deltap + dp_dpp
            deltav = deltav + ptl%v * dp_dpp / ptl%p
        endif

        ! Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            dvx_dy = fields(nfields+2)
            dvy_dx = fields(nfields+4)
            dvz_dx = fields(nfields+7)
            dvz_dy = fields(nfields+8)
            sigmaxx = dvx_dx - divv / 3
            sigmayy = dvy_dy - divv / 3
            sigmazz = -divv / 3
            sigmaxy = (dvx_dy + dvy_dx) / 2
            sigmaxz = dvz_dx / 2
            sigmayz = dvz_dy / 2
            if (weak_scattering) then
                bbsigma = sigmaxx * bx**2 + sigmayy * by**2 + sigmazz * bz**2 + &
                    2.0 * (sigmaxy * bx * by + sigmaxz * bx * bz + sigmayz * by * bz)
                bbsigma = bbsigma * ib * ib
                gshear = bbsigma**2 / 5
            else
                gshear = 2 * (sigmaxx**2 + sigmayy**2 + sigmazz**2 + &
                    2 * (sigmaxy**2 + sigmaxz**2 + sigmayz**2)) / 15
            endif
            dp_dpp = 0.0d0
            if (gshear > 0) then
                dp_dpp = (2 + pindex) * gshear * tau0 * kappa%knorm0 * &
                    ptl%p**(pindex-1) * p0**(2.0-pindex) * ptl%dt
                if (.not. dpp_wave_flag) then
                    ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
                endif
                dp_dpp = dsqrt(2 * gshear * tau0 * kappa%knorm0 * &
                    ptl%p**pindex * p0**(2.0-pindex)) * ran1 * sdt
            endif
            deltap = deltap + dp_dpp
            deltav = deltav + ptl%v * dp_dpp / ptl%p
        endif

        ! Pitch-angle evolution
        div_bnorm = -(bx * db_dx + by * db_dy * ir) * ib2
        dmu_dt = ptl%v * div_bnorm + ptl%mu * divv - &
            3 * ptl%mu * bb_gradv - 2 * bv_gradv / ptl%v
        dmu_dt = dmu_dt * muf1
        h0 = 0.2  ! Enabling pitch-angle scattering around mu=0
        turb_gamma = 3.0 - pindex  ! Turbulence spectral index
        dtmp = abs(ptl%mu)**(turb_gamma-1) + h0
        duu = duu0 * (1-mu2) * dtmp
        if (ptl%mu .gt. 0.0d0) then
            duu_du = duu0 * (-2*ptl%mu * dtmp + &
                             (1-mu2) * abs(ptl%mu)**(turb_gamma-2))
        else if (ptl%mu .lt. 0.0d0) then
            duu_du = duu0 * (-2*ptl%mu * dtmp - &
                             (1-mu2) * abs(ptl%mu)**(turb_gamma-2))
        else
            duu_du = 0.0d0
        endif

        if (momentum_dependency == 1) then
            dtmp = (ptl%p / p0)**(turb_gamma-1)
            duu = duu * dtmp
            duu_du = duu_du * dtmp
        endif
        deltamu = (dmu_dt + duu_du) * ptl%dt
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltamu = deltamu + ran1*dsqrt(2.0*duu)*sdt
        ptl%mu = ptl%mu + deltamu

        if (acc_region_flag == 1) then
            if (particle_in_acceleration_region(ptl)) then
                ptl%p = ptl%p + deltap
                ptl%v = ptl%v + deltav
            else
                deltap = 0.0
                deltav = 0.0
            endif
        else
            ptl%p = ptl%p + deltap
            ptl%v = ptl%v + deltav
        endif

        if (ptl%p < 0.5 * p0) then
            ptl%v = ptl%v - deltav
            deltav = ptl%v * 0.5 * p0 / ptl%p - ptl%v
            ptl%v = ptl%v + deltav
            ptl%p = ptl%p - deltap
            deltap = 0.5 * p0 - ptl%p
            ptl%p = 0.5 * p0
        endif
    end subroutine push_particle_2d_include_3rd_ft

    !---------------------------------------------------------------------------
    !< Push particle for a single step for a 3D simulation
    !< Args:
    !<  thread_id: thread ID staring from 0
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  surface_height1: the height of surface 1 separating the acceleration region
    !<  surface_height2: the height of surface 2 separating the acceleration region
    !<  ptl: particle structure
    !<  fields: fields and their gradients at particle position
    !<  db2: turbulence variance and its gradients at particle position
    !<  lc: turbulence correlation length and its gradients at particle position
    !<  kappa: kappa and related variables
    !<  deltax, deltay, deltaz, deltap: the change of x, y, z and p in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_3d(thread_id, rt, surface_height1, surface_height2, &
            ptl, fields, db2, lc, kappa, deltax, deltay, deltaz, deltap)
        use constants, only: pi
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: interp_fields, &
            interp_magnetic_fluctuation, interp_correlation_length
        use acc_region_surface, only: check_above_acc_surface
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        real(dp), intent(in) :: rt, surface_height1, surface_height2
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        real(dp), dimension(*), intent(inout) :: db2, lc
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltay, deltaz, deltap
        real(dp) :: rt1, sdt
        real(dp) :: dvx_dx, dvx_dy, dvx_dz
        real(dp) :: dvy_dx, dvy_dy, dvy_dz
        real(dp) :: dvz_dx, dvz_dy, dvz_dz
        real(dp) :: divv, gshear
        real(dp) :: vx, vy, vz
        real(dp) :: bx, by, bz, b, ib, ib2, ib3
        real(dp) :: bxn, byn, bzn, bxyn, ibxyn
        real(dp) :: px, py, pz
        real(dp) :: dbx_dy, dbx_dz, dby_dx, dby_dz, dbz_dx, dbz_dy
        real(dp) :: db_dx, db_dy, db_dz
        real(dp) :: vdx, vdy, vdz, vdp
        real(dp) :: xmin, ymin, zmin, xmax, ymax, zmax, dxm, dym, dzm
        reaL(dp) :: xmin1, ymin1, zmin1, xmax1, ymax1, zmax1, dxmh, dymh, dzmh
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rho, va ! Plasma density and Alfven speed
        real(dp) :: rands(2)
        real(dp) :: ctheta, istheta, ir, ir2
        real(dp) :: gbr, gbt, gbp, p11, p12, p13, p22, p23, p33
        real(dp) :: sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, sigmayz ! shear tensor
        real(dp) :: bbsigma ! b_ib_jsigma_ij
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights
        logical :: in_acc_region

        vx = fields(1)
        vy = fields(2)
        vz = fields(3)
        bx = fields(5)
        by = fields(6)
        bz = fields(7)
        b = dsqrt(bx**2 + by**2 + bz**2)
        if (b == 0.0d0) then
            ib = 0.0_dp
        else
            ib = 1.0_dp / b
        endif
        bxn = bx * ib
        byn = by * ib
        bzn = bz * ib
        bxyn = dsqrt(bxn**2 + byn**2)
        if (bxyn == 0.0d0) then
            ibxyn = 0.0d0
        else
            ibxyn = 1.0_dp / bxyn
        endif
        dvx_dx = fields(nfields+1)
        dvy_dy = fields(nfields+5)
        dvz_dz = fields(nfields+9)
        xmin = fconfig%xmin
        ymin = fconfig%ymin
        zmin = fconfig%zmin
        xmax = fconfig%xmax
        ymax = fconfig%ymax
        zmax = fconfig%zmax
        dxm = mhd_config%dx
        dym = mhd_config%dy
        dzm = mhd_config%dz
        dxmh = 0.5 * dxm
        dymh = 0.5 * dym
        dzmh = 0.5 * dzm

        if (spherical_coord_flag) then
            ctheta = cos(ptl%y)
            istheta = 1.0 / sin(ptl%y)
            ir = 1.0 / ptl%x
            ir2 = 1.0 / ptl%x**2
            p11 = dsqrt(2.0*(kappa%kxx*kappa%kyz**2 + &
                             kappa%kyy*kappa%kxz**2 + &
                             kappa%kzz*kappa%kxy**2 - &
                             2.0*kappa%kxy*kappa%kxz*kappa%kyz - &
                             kappa%kxx*kappa%kyy*kappa%kzz) / &
                             (kappa%kyz**2 - kappa%kyy*kappa%kzz))
            p12 = (kappa%kxz*kappa%kyz - kappa%kxy*kappa%kzz) * &
                dsqrt(2.0*(kappa%kyy - (kappa%kyz**2/kappa%kzz))) / &
                (kappa%kyz**2 - kappa%kyy*kappa%kzz)
            p13 = dsqrt(2.0 / kappa%kzz) * kappa%kxz
            p22 = dsqrt(2.0 * (kappa%kyy - kappa%kyz**2/kappa%kzz)) * ir
            p23 = dsqrt(2.0 / kappa%kzz) * kappa%kyz * ir
            p33 = dsqrt(2.0 * kappa%kzz) * istheta * ir
        endif

        ! Drift velocity
        dbx_dy = fields(nfields+14)
        dbx_dz = fields(nfields+15)
        dby_dx = fields(nfields+16)
        dby_dz = fields(nfields+18)
        dbz_dx = fields(nfields+19)
        dbz_dy = fields(nfields+20)
        db_dx = fields(nfields+22)
        db_dy = fields(nfields+23)
        db_dz = fields(nfields+24)
        ib2 = ib * ib
        ib3 = ib * ib2
        vdp = 1.0 / (3 * pcharge) / dsqrt((drift1*p0/ptl%p)**2 + (drift2*p0**2/ptl%p**2)**2)
        if (spherical_coord_flag) then
            gbr = db_dx
            gbt = db_dy * ir
            gbp = db_dz * ir * istheta
            vdx = vdp * ((dbz_dy + bz*ctheta*istheta - dby_dz*istheta)*ir*ib2 - &
                         (gbt*bz - gbp*by)*2.0*ib3)
            vdy = vdp * ((dbx_dz*istheta*ir - dbz_dx - bz*ir)*ib2 - &
                         (gbp*bx - gbr*bz)*2.0*ib3)
            vdz = vdp * ((dby_dx + by*ir - dbx_dy*ir)*ib2 - &
                         (gbr*by - gbt*bx)*2.0*ib3)
        else
            vdx = vdp * ((dbz_dy-dby_dz)*ib2 - 2*(bz*db_dy-by*db_dz)*ib3)
            vdy = vdp * ((dbx_dz-dbz_dx)*ib2 - 2*(bx*db_dz-bz*db_dx)*ib3)
            vdz = vdp * ((dby_dx-dbx_dy)*ib2 - 2*(by*db_dx-bx*db_dy)*ib3)
        endif

        !< The field data has two ghost cells, so the particles can cross the
        !< boundary without causing segment fault errors
        xmin1 = xmin - dxmh
        ymin1 = ymin - dymh
        xmax1 = xmax + dxmh
        ymax1 = ymax + dymh

        if (spherical_coord_flag) then
            deltax = (vx + vdx + kappa%dkxx_dx + &
                (2.0*kappa%kxx + kappa%dkxy_dy + &
                kappa%kxy*ctheta*istheta + kappa%dkxz_dz*istheta)*ir) * ptl%dt
            deltay = ((vy + vdy + kappa%dkxy_dx)*ir + &
                (kappa%kxy + kappa%dkyy_dy + kappa%kyy*ctheta*istheta + &
                kappa%dkyz_dz*istheta)*ir2) * ptl%dt
            deltaz = ((vz + vdz + kappa%dkxz_dx)*istheta*ir + &
                (kappa%kxz + kappa%dkyz_dy + kappa%dkzz_dz*istheta)*istheta*ir2) * ptl%dt
        else
            deltax = (vx + vdx + kappa%dkxx_dx + kappa%dkxy_dy + kappa%dkxz_dz) * ptl%dt
            deltay = (vy + vdy + kappa%dkxy_dx + kappa%dkyy_dy + kappa%dkyz_dz) * ptl%dt
            deltaz = (vz + vdz + kappa%dkxz_dx + kappa%dkyz_dy + kappa%dkzz_dz) * ptl%dt
        endif

        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran2 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran3 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3

        sdt = dsqrt(ptl%dt)

        if (spherical_coord_flag) then
            deltax = deltax + (p11*ran1 + p12*ran2 + p13*ran3)*sdt
            deltay = deltay + (p22*ran2 + p23*ran3)*sdt
            deltaz = deltaz + p33*ran3*sdt
        else
            deltax = deltax + &
                (bxn*kappa%skpara*ran1 - bxn*bzn*kappa%skperp*ibxyn*ran2 - &
                 byn*kappa%skperp*ibxyn*ran3)*sdt
            deltay = deltay + &
                (byn*kappa%skpara*ran1 - byn*bzn*kappa%skperp*ibxyn*ran2 + &
                 bxn*kappa%skperp*ibxyn*ran3)*sdt
            deltaz = deltaz + (bzn*kappa%skpara*ran1 + bxyn*kappa%skperp*ran2)*sdt
        endif

        ! Momentum
        if (spherical_coord_flag) then
            divv = dvx_dx + &
                (2.0*vx + dvy_dy + vy*ctheta*istheta + dvz_dz*istheta)*ir
        else
            divv = dvx_dx + dvy_dy + dvz_dz
        endif
        deltap = -ptl%p * divv * ptl%dt / 3.0d0
        ! Momentum diffusion due to wave scattering
        if (dpp_wave_flag) then
            rho = fields(4)
            va = b / dsqrt(rho)
            if (momentum_dependency == 1) then
                deltap = deltap + (8*ptl%p / (27*kappa%kpara)) * va**2 * ptl%dt
            else
                deltap = deltap + (4*ptl%p / (9*kappa%kpara)) * va**2 * ptl%dt
            endif
            ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
            deltap = deltap + ran1 * va * ptl%p * dsqrt(2/(9*kappa%kpara)) * sdt
        endif

        ! Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            dvx_dy = fields(nfields+2)
            dvx_dz = fields(nfields+3)
            dvy_dx = fields(nfields+4)
            dvy_dz = fields(nfields+6)
            dvz_dx = fields(nfields+7)
            dvz_dy = fields(nfields+8)
            sigmaxx = dvx_dx - divv / 3
            sigmayy = dvy_dy - divv / 3
            sigmazz = dvz_dz - divv / 3
            sigmaxy = (dvx_dy + dvy_dx) / 2
            sigmaxz = (dvx_dz + dvz_dx) / 2
            sigmayz = (dvy_dz + dvz_dy) / 2
            if (weak_scattering) then
                bbsigma = sigmaxx * bx**2 + sigmayy * by**2 + sigmazz * bz**2 + &
                    2.0 * (sigmaxy * bx * by + sigmaxz * bx * bz + sigmayz * by * bz)
                bbsigma = bbsigma * ib * ib
                gshear = bbsigma**2 / 5
            else
                gshear = 2 * (sigmaxx**2 + sigmayy**2 + sigmazz**2 + &
                    2 * (sigmaxy**2 + sigmaxz**2 + sigmayz**2)) / 15
            endif
            deltap = deltap + (2 + pindex) * gshear * tau0 * kappa%knorm0 * &
                ptl%p**(pindex-1) * p0**(2.0-pindex) * ptl%dt
            if (.not. dpp_wave_flag) then
                ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
            endif
            deltap = deltap + dsqrt(2 * gshear * tau0 * kappa%knorm0 * &
                ptl%p**pindex * p0**(2.0-pindex)) * ran1 * sdt
        endif

        if (acc_region_flag == 1) then
            in_acc_region = particle_in_acceleration_region(ptl)
            if (acc_by_surface_flag) then
                in_acc_region = in_acc_region .and. &
                    check_above_acc_surface(ptl%x, ptl%y, ptl%z, surface_height1, surface_height2)
            endif
            if (in_acc_region) then
                ptl%p = ptl%p + deltap
            else
                deltap = 0.0
            endif
        else
            ptl%p = ptl%p + deltap
        endif
        if (ptl%p < 0.5 * p0) then
            ptl%p = ptl%p - deltap
            deltap = 0.5 * p0 - ptl%p
            ptl%p = 0.5 * p0
        endif

        ptl%x = ptl%x + deltax
        ptl%y = ptl%y + deltay
        ptl%z = ptl%z + deltaz
        ptl%t = ptl%t + ptl%dt
    end subroutine push_particle_3d

    !---------------------------------------------------------------------------
    !< Push particle for a single step for a 3D simulation for focused transport
    !< Args:
    !<  thread_id: thread ID staring from 0
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  surface_height1: the height of surface 1 separating the acceleration region
    !<  surface_height2: the height of surface 2 separating the acceleration region
    !<  ptl: particle structure
    !<  fields: fields and their gradients at particle position
    !<  db2: turbulence variance and its gradients at particle position
    !<  lc: turbulence correlation length and its gradients at particle position
    !<  kappa: kappa and related variables
    !<  deltax: the change of x in this step
    !<  deltay: the change of y in this step
    !<  deltaz: the change of z in this step
    !<  deltap: the change of p in this step
    !<  deltav: the change of v in this step
    !<  deltamu: the change of mu in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_3d_ft(thread_id, rt, surface_height1, surface_height2, &
            ptl, fields, db2, lc, kappa, deltax, deltay, deltaz, deltap, &
            deltav, deltamu)
        use constants, only: pi
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: interp_fields, &
            interp_magnetic_fluctuation, interp_correlation_length
        use acc_region_surface, only: check_above_acc_surface
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        real(dp), intent(in) :: rt, surface_height1, surface_height2
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        real(dp), dimension(*), intent(inout) :: db2, lc
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltay, deltaz, deltap, deltav, deltamu
        real(dp) :: rt1, sdt
        real(dp) :: dvx_dx, dvx_dy, dvx_dz
        real(dp) :: dvy_dx, dvy_dy, dvy_dz
        real(dp) :: dvz_dx, dvz_dy, dvz_dz
        real(dp) :: divv, bb_gradv, bv_gradv, gshear
        real(dp) :: vx, vy, vz, vbx, vby, vbz
        real(dp) :: bx, by, bz, b, ib, ib2, ib3
        real(dp) :: bxn, byn, bzn, bxyn, ibxyn
        real(dp) :: px, py, pz
        real(dp) :: dbx_dx, dbx_dy, dbx_dz
        real(dp) :: dby_dx, dby_dy, dby_dz
        real(dp) :: dbz_dx, dbz_dy, dbz_dz
        real(dp) :: db_dx, db_dy, db_dz
        real(dp) :: vdx, vdy, vdz, vdp
        real(dp) :: xmin, ymin, zmin, xmax, ymax, zmax, dxm, dym, dzm
        reaL(dp) :: xmin1, ymin1, zmin1, xmax1, ymax1, zmax1, dxmh, dymh, dzmh
        real(dp) :: mu2, muf1, muf2, kx, ky, kz, bdot_curvb, acc_rate
        real(dp) :: dp_dpp, div_bnorm, h0, turb_gamma
        real(dp) :: duu, dmu_dt, duu_du, dtmp
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rho, va ! Plasma density and Alfven speed
        real(dp) :: rands(2)
        real(dp) :: ctheta, istheta, cot_theta, ir, ir2
        real(dp) :: gbr, gbt, gbp, p11, p12, p13, p22, p23, p33
        real(dp) :: sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, sigmayz ! shear tensor
        real(dp) :: bbsigma ! b_ib_jsigma_ij
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights
        logical :: in_acc_region

        vx = fields(1)
        vy = fields(2)
        vz = fields(3)
        bx = fields(5)
        by = fields(6)
        bz = fields(7)
        b = dsqrt(bx**2 + by**2 + bz**2)
        if (b == 0.0d0) then
            ib = 0.0_dp
        else
            ib = 1.0_dp / b
        endif
        bxn = bx * ib
        byn = by * ib
        bzn = bz * ib
        bxyn = dsqrt(bxn**2 + byn**2)
        if (bxyn == 0.0d0) then
            ibxyn = 0.0d0
        else
            ibxyn = 1.0_dp / bxyn
        endif
        xmin = fconfig%xmin
        ymin = fconfig%ymin
        zmin = fconfig%zmin
        xmax = fconfig%xmax
        ymax = fconfig%ymax
        zmax = fconfig%zmax
        dxm = mhd_config%dx
        dym = mhd_config%dy
        dzm = mhd_config%dz
        dxmh = 0.5 * dxm
        dymh = 0.5 * dym
        dzmh = 0.5 * dzm

        if (spherical_coord_flag) then
            ctheta = cos(ptl%y)
            istheta = 1.0 / sin(ptl%y)
            cot_theta = ctheta * istheta
            ir = 1.0 / ptl%x
            ir2 = 1.0 / ptl%x**2
            p11 = dsqrt(2.0*(kappa%kxx*kappa%kyz**2 + &
                             kappa%kyy*kappa%kxz**2 + &
                             kappa%kzz*kappa%kxy**2 - &
                             2.0*kappa%kxy*kappa%kxz*kappa%kyz - &
                             kappa%kxx*kappa%kyy*kappa%kzz) / &
                             (kappa%kyz**2 - kappa%kyy*kappa%kzz))
            p12 = (kappa%kxz*kappa%kyz - kappa%kxy*kappa%kzz) * &
                dsqrt(2.0*(kappa%kyy - (kappa%kyz**2/kappa%kzz))) / &
                (kappa%kyz**2 - kappa%kyy*kappa%kzz)
            p13 = dsqrt(2.0 / kappa%kzz) * kappa%kxz
            p22 = dsqrt(2.0 * (kappa%kyy - kappa%kyz**2/kappa%kzz)) * ir
            p23 = dsqrt(2.0 / kappa%kzz) * kappa%kyz * ir
            p33 = dsqrt(2.0 * kappa%kzz) * istheta * ir
        endif

        ! Drift velocity
        dbx_dx = fields(nfields+13)
        dbx_dy = fields(nfields+14)
        dbx_dz = fields(nfields+15)
        dby_dx = fields(nfields+16)
        dby_dy = fields(nfields+17)
        dby_dz = fields(nfields+18)
        dbz_dx = fields(nfields+19)
        dbz_dy = fields(nfields+20)
        dbz_dz = fields(nfields+21)
        db_dx = fields(nfields+22)
        db_dy = fields(nfields+23)
        db_dz = fields(nfields+24)
        ib2 = ib * ib
        ib3 = ib * ib2
        vdp = 1.0 / pcharge / dsqrt((drift1*p0/ptl%p)**2 + (drift2*p0**2/ptl%p**2)**2)
        mu2 = ptl%mu**2
        muf1 = 0.5 * (1.0 - mu2)
        muf2 = 0.5 * (3.0*mu2 - 1.0)
        if (spherical_coord_flag) then
            gbr = db_dx
            gbt = db_dy * ir
            gbp = db_dz * ir * istheta
            kx = bx * dbx_dx + ir * (by * dbx_dy + &
                                     bz * dbx_dz * istheta - &
                                    (by**2 + bz**2))
            ky = bx * dby_dx + ir * (by * dby_dy + &
                                     bz * dby_dz * istheta + &
                                     (bx*by - cot_theta*bz**2))
            kz = bx * dbz_dx + ir * (by * dbz_dy + &
                                     bz * dbz_dz * istheta + &
                                     (bx*bz + cot_theta*by*bz))
            bdot_curvb = bx * (dbz_dy + cot_theta * bz - dby_dz * istheta) * ir + &
                         by * (dbx_dz * istheta * ir - dbz_dx - bz * ir) + &
                         bz * (dby_dx + by * ir - dbx_dy * ir)
            vdx = vdp * (muf1 * (by * gbp - bz * gbt) * ib2 + &
                         mu2 * (by * kz - bz * ky) * ib3 + &
                         muf1 * bx * bdot_curvb * ib3)
            vdy = vdp * (muf1 * (bz * gbr - bx * gbp) * ib2 + &
                         mu2 * (bz * kx - bx * kz) * ib3 + &
                         muf1 * by * bdot_curvb * ib3)
            vdz = vdp * (muf1 * (bx * gbt - by * gbr) * ib2 + &
                         mu2 * (bx * ky - by * kx) * ib3 + &
                         muf1 * bz * bdot_curvb * ib3)
        else
            kx = bx * dbx_dx + by * dbx_dy + bz * dbx_dz
            ky = bx * dby_dx + by * dby_dy + bz * dby_dz
            kz = bx * dbz_dx + by * dbz_dy + bz * dbz_dz
            bdot_curvb = bx * (dbz_dy - dby_dz) + &
                         by * (dbx_dz - dbz_dx) + &
                         bz * (dby_dx - dbx_dy)
            vdx = vdp * (muf1 * (by * db_dz - bz * db_dy) * ib2 + &
                         mu2 * (by * kz - bz * ky) * ib3 + &
                         muf1 * bx * bdot_curvb * ib3)
            vdy = vdp * (muf1 * (bz * db_dx - bx * db_dz) * ib2 + &
                         mu2 * (bz * kx - bx * kz) * ib3 + &
                         muf1 * by * bdot_curvb * ib3)
            vdz = vdp * (muf1 * (bx * db_dy - by * db_dx) * ib2 + &
                         mu2 * (bx * ky - by * kx) * ib3 + &
                         muf1 * bz * bdot_curvb * ib3)
        endif

        !< The field data has two ghost cells, so the particles can cross the
        !< boundary without causing segment fault errors
        xmin1 = xmin - dxmh
        ymin1 = ymin - dymh
        xmax1 = xmax + dxmh
        ymax1 = ymax + dymh

        vbx = ptl%v * ptl%mu * ib
        vby = vbx * by ! particle velocity along the magnetic field
        vbz = vbx * bz
        vbx = vbx * bx
        if (spherical_coord_flag) then
            deltax = (vx + vbx + vdx + kappa%dkxx_dx + &
                (2.0*kappa%kxx + kappa%dkxy_dy + &
                kappa%kxy*ctheta*istheta + kappa%dkxz_dz*istheta)*ir) * ptl%dt
            deltay = ((vy + vby + vdy + kappa%dkxy_dx)*ir + &
                (kappa%kxy + kappa%dkyy_dy + kappa%kyy*ctheta*istheta + &
                kappa%dkyz_dz*istheta)*ir2) * ptl%dt
            deltaz = ((vz + vbz + vdz + kappa%dkxz_dx)*istheta*ir + &
                (kappa%kxz + kappa%dkyz_dy + kappa%dkzz_dz*istheta)*istheta*ir2) * ptl%dt
        else
            deltax = (vx + vbx + vdx + kappa%dkxx_dx + kappa%dkxy_dy + kappa%dkxz_dz) * ptl%dt
            deltay = (vy + vby + vdy + kappa%dkxy_dx + kappa%dkyy_dy + kappa%dkyz_dz) * ptl%dt
            deltaz = (vz + vbz + vdz + kappa%dkxz_dx + kappa%dkyz_dy + kappa%dkzz_dz) * ptl%dt
        endif

        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran2 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran3 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3

        sdt = dsqrt(ptl%dt)

        if (spherical_coord_flag) then
            deltax = deltax + (p11*ran1 + p12*ran2 + p13*ran3)*sdt
            deltay = deltay + (p22*ran2 + p23*ran3)*sdt
            deltaz = deltaz + p33*ran3*sdt
        else
            deltax = deltax + (-bxn*bzn*kappa%skperp*ibxyn*ran1 - &
                               byn*kappa%skperp*ibxyn*ran2)*sdt
            deltay = deltay + (-byn*bzn*kappa%skperp*ibxyn*ran1 + &
                               bxn*kappa%skperp*ibxyn*ran2)*sdt
            deltaz = deltaz + bxyn*kappa%skperp*ran1*sdt
        endif

        ! Momentum and velocity
        dvx_dx = fields(nfields+1)
        dvx_dy = fields(nfields+2)
        dvx_dz = fields(nfields+3)
        dvy_dx = fields(nfields+4)
        dvy_dy = fields(nfields+5)
        dvy_dz = fields(nfields+6)
        dvz_dx = fields(nfields+7)
        dvz_dy = fields(nfields+8)
        dvz_dz = fields(nfields+9)
        if (spherical_coord_flag) then
            divv = dvx_dx + (2.0*vx + dvy_dy + vy*cot_theta + dvz_dz*istheta)*ir
            dtmp = istheta * ir
            bb_gradv = (bx*bx*dvx_dx + bx*by*dvx_dy*ir + bx*bz*dvx_dz*dtmp - bx*(by*vy + bz*vz)*ir + &
                        by*bx*dvy_dx + by*by*dvy_dy*ir + by*bz*dvy_dz*dtmp + by*(by*vx - cot_theta*bz*vz)*ir + &
                        bz*bx*dvz_dx + bz*by*dvz_dy*ir + bz*bz*dvz_dz*dtmp + bz**2*(vx + cot_theta*vy)*ir) * ib2
            bv_gradv = (bx*vx*dvx_dx + bx*vy*dvx_dy*ir + bx*vz*dvx_dz*dtmp - bx*(vy*vy + vz*vz)*ir + &
                        by*vx*dvy_dx + by*vy*dvy_dy*ir + by*vz*dvy_dz*dtmp + by*(vy*vx - cot_theta*vz*vz)*ir + &
                        bz*vx*dvz_dx + bz*vy*dvz_dy*ir + bz*vz*dvz_dz*dtmp + bz*vz*(vx + cot_theta*vy)*ir) * ib
        else
            divv = dvx_dx + dvy_dy + dvz_dz
            bb_gradv = (bx * (bx * dvx_dx + by * dvx_dy + bz * dvx_dz) + &
                        by * (bx * dvy_dx + by * dvy_dy + bz * dvy_dz) + &
                        bz * (bx * dvz_dx + by * dvz_dy + bz * dvz_dz)) * ib2
            bv_gradv = (bx * (vx * dvx_dx + vy * dvx_dy + vz * dvx_dz) + &
                        by * (vx * dvy_dx + vy * dvy_dy + vz * dvy_dz) + &
                        bz * (vx * dvz_dx + vy * dvz_dy + vz * dvz_dz)) * ib
        endif
        acc_rate = -(muf1 * divv + muf2 * bb_gradv + ptl%mu * bv_gradv / ptl%v)
        deltap = ptl%p * acc_rate * ptl%dt
        deltav = ptl%v * acc_rate * ptl%dt

        ! Momentum diffusion due to wave scattering
        if (dpp_wave_flag) then
            rho = fields(4)
            va = b / dsqrt(rho)
            if (momentum_dependency == 1) then
                dp_dpp = (8*ptl%p / (27*kappa%kpara)) * va**2 * ptl%dt
            else
                dp_dpp = (4*ptl%p / (9*kappa%kpara)) * va**2 * ptl%dt
            endif
            ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
            dp_dpp = dp_dpp + ran1 * va * ptl%p * dsqrt(2/(9*kappa%kpara)) * sdt
            deltap = deltap + dp_dpp
            deltav = deltav + ptl%v * dp_dpp / ptl%p
        endif

        ! Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            sigmaxx = dvx_dx - divv / 3
            sigmayy = dvy_dy - divv / 3
            sigmazz = dvz_dz - divv / 3
            sigmaxy = (dvx_dy + dvy_dx) / 2
            sigmaxz = (dvx_dz + dvz_dx) / 2
            sigmayz = (dvy_dz + dvz_dy) / 2
            if (weak_scattering) then
                bbsigma = sigmaxx * bx**2 + sigmayy * by**2 + sigmazz * bz**2 + &
                    2.0 * (sigmaxy * bx * by + sigmaxz * bx * bz + sigmayz * by * bz)
                bbsigma = bbsigma * ib * ib
                gshear = bbsigma**2 / 5
            else
                gshear = 2 * (sigmaxx**2 + sigmayy**2 + sigmazz**2 + &
                    2 * (sigmaxy**2 + sigmaxz**2 + sigmayz**2)) / 15
            endif
            dp_dpp = 0.0d0
            if (gshear > 0) then
                dp_dpp = (2 + pindex) * gshear * tau0 * kappa%knorm0 * &
                    ptl%p**(pindex-1) * p0**(2.0-pindex) * ptl%dt
                if (.not. dpp_wave_flag) then
                    ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
                endif
                dp_dpp = dsqrt(2 * gshear * tau0 * kappa%knorm0 * &
                    ptl%p**pindex * p0**(2.0-pindex)) * ran1 * sdt
            endif
            deltap = deltap + dp_dpp
            deltav = deltav + ptl%v * dp_dpp / ptl%p
        endif

        ! Pitch-angle evolution
        div_bnorm = -(bx * db_dx + by * db_dy * ir + bz * db_dz * ir * istheta) * ib2
        dmu_dt = ptl%v * div_bnorm + ptl%mu * divv - &
            3 * ptl%mu * bb_gradv - 2 * bv_gradv / ptl%v
        dmu_dt = dmu_dt * muf1
        h0 = 0.2  ! Enabling pitch-angle scattering around mu=0
        turb_gamma = 3.0 - pindex  ! Turbulence spectral index
        dtmp = abs(ptl%mu)**(turb_gamma-1) + h0
        duu = duu0 * (1-mu2) * dtmp
        if (ptl%mu .gt. 0.0d0) then
            duu_du = duu0 * (-2*ptl%mu * dtmp + &
                             (1-mu2) * abs(ptl%mu)**(turb_gamma-2))
        else if (ptl%mu .lt. 0.0d0) then
            duu_du = duu0 * (-2*ptl%mu * dtmp - &
                             (1-mu2) * abs(ptl%mu)**(turb_gamma-2))
        else
            duu_du = 0.0d0
        endif

        if (momentum_dependency == 1) then
            dtmp = (ptl%p / p0)**(turb_gamma-1)
            duu = duu * dtmp
            duu_du = duu_du * dtmp
        endif
        deltamu = (dmu_dt + duu_du) * ptl%dt
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltamu = deltamu + ran1*dsqrt(2.0*duu)*sdt
        ptl%mu = ptl%mu + deltamu

        if (acc_region_flag == 1) then
            in_acc_region = particle_in_acceleration_region(ptl)
            if (acc_by_surface_flag) then
                in_acc_region = in_acc_region .and. &
                    check_above_acc_surface(ptl%x, ptl%y, ptl%z, surface_height1, surface_height2)
            endif
            if (in_acc_region) then
                ptl%p = ptl%p + deltap
                ptl%v = ptl%v + deltav
            else
                deltap = 0.0
                deltav = 0.0
            endif
        else
            ptl%p = ptl%p + deltap
            ptl%v = ptl%v + deltav
        endif
        if (ptl%p < 0.5 * p0) then
            ptl%v = ptl%v - deltav
            deltav = ptl%v * 0.5 * p0 / ptl%p - ptl%v
            ptl%v = ptl%v + deltav
            ptl%p = ptl%p - deltap
            deltap = 0.5 * p0 - ptl%p
            ptl%p = 0.5 * p0
        endif

        ptl%x = ptl%x + deltax
        ptl%y = ptl%y + deltay
        ptl%z = ptl%z + deltaz
        ptl%t = ptl%t + ptl%dt
    end subroutine push_particle_3d_ft

    !---------------------------------------------------------------------------
    !< Resize escaped particle array
    !---------------------------------------------------------------------------
    subroutine resize_escaped_particles
        implicit none
        type(particle_type), allocatable, dimension(:) :: escaped_ptls_tmp
        integer :: i, nmax
        allocate(escaped_ptls_tmp(nptl_escaped_max))
        escaped_ptls_tmp = escaped_ptls
        ! Reallocate the escaped particle data array
        nmax = max(int(1.25 * nptl_escaped_max), nptl_escaped_max + nptl_current)
        nptl_escaped_max = nmax
        deallocate(escaped_ptls)
        allocate(escaped_ptls(nptl_escaped_max))
        escaped_ptls%x = 0.0
        escaped_ptls%y = 0.0
        escaped_ptls%z = 0.0
        escaped_ptls%p = 0.0
        escaped_ptls%v = 0.0
        escaped_ptls%mu = 0.0
        escaped_ptls%weight = 0.0
        escaped_ptls%t = 0.0
        escaped_ptls%dt = 0.0
        escaped_ptls%split_times = 0
        escaped_ptls%count_flag = COUNT_FLAG_OTHERS
        escaped_ptls%tag = 0
        escaped_ptls%nsteps_tracking = 0
        escaped_ptls(:nptl_escaped) = escaped_ptls_tmp(:nptl_escaped)
        deallocate(escaped_ptls_tmp)
    end subroutine resize_escaped_particles

    !---------------------------------------------------------------------------
    !< Remove particles from simulation if their count_flags are 0.
    !< Args:
    !<  dump_escaped: whether to dump escaped particles
    !---------------------------------------------------------------------------
    subroutine remove_particles(dump_escaped)
        implicit none
        logical, intent(in) :: dump_escaped
        integer :: i, nremoved
        type(particle_type) :: ptl1

        ! Resize the escaped particle data array if necessary
        if (dump_escaped) then
            if ((nptl_escaped + nptl_current) > nptl_escaped_max) then
                call resize_escaped_particles
            endif
        endif

        if (nptl_current > 0) then
            nremoved = 0
            i = 1
            do while (i <= nptl_current)
                if ((nptl_current - i) == (nremoved - 1)) then
                    exit
                endif
                if (ptls(i)%count_flag == COUNT_FLAG_INBOX) then
                    i = i + 1
                else
                    if (ptls(i)%count_flag == COUNT_FLAG_ESCAPE) then
                        ! Copy escaped particles
                        nptl_escaped = nptl_escaped + 1
                        if (dump_escaped) then
                            escaped_ptls(nptl_escaped) = ptls(i)
                        endif
                    endif
                    ptl1 = ptls(nptl_current-nremoved)
                    ptls(nptl_current-nremoved) = ptls(i)
                    ptls(i) = ptl1
                    nremoved = nremoved + 1
                endif
            enddo
            nptl_current = nptl_current - nremoved
        endif
    end subroutine remove_particles

    !---------------------------------------------------------------------------
    !< Add particles from neighbors
    !---------------------------------------------------------------------------
    subroutine add_neighbor_particles
        implicit none
        integer :: i, j, nrecv
        nptl_old = nptl_current
        do i = 1, ndim_field*2
            nrecv = nrecvers(i)
            if (nrecv > 0) then
                ptls(nptl_current+1:nptl_current+nrecv) = recvers(1:nrecv, i)
                nptl_current = nptl_current + nrecv
            endif
        enddo
    end subroutine add_neighbor_particles

    !---------------------------------------------------------------------------
    !< When a particle get to a certain energy, split it to two particles
    !< to increase statistics at high energy
    !< Args:
    !<  split_ratio (> 1.0): momentum increase ratio for particle splitting
    !<  pmin_split (> 1.0): the minimum momentum (in terms in p0) to start
    !<                      splitting particles
    !---------------------------------------------------------------------------
    subroutine split_particle(split_ratio, pmin_split)
        implicit none
        real(dp), intent(in) :: split_ratio, pmin_split
        integer :: i, nptl
        real(dp) :: p_threshold
        type(particle_type) :: ptl, ptl_new
        nptl = nptl_current
        do i = 1, nptl
            ptl = ptls(i)
            p_threshold = pmin_split * p0 * split_ratio**ptl%split_times
            if (ptl%p > p_threshold .and. ptl%p <= pmax) then
                nptl_current = nptl_current + 1
                if (nptl_current > nptl_max) then
                    nptl_current = nptl_max
                    return
                endif
                nptl_split = nptl_split + 1
                ptl%weight = 0.5**(1.0 + ptl%split_times)
                ptl%split_times = ptl%split_times + 1
                ptls(nptl_current) = ptl
                tag_max = tag_max + 1
                ptls(nptl_current)%tag = tag_max
                ptls(i) = ptl
            endif
        enddo
    end subroutine split_particle

    !---------------------------------------------------------------------------
    !< quicksort particles w.r.t their energy/momentum
    !---------------------------------------------------------------------------
    recursive subroutine quicksort_particle(ptls, first, last)
        implicit none
        type(particle_type), dimension(:) :: ptls
        type(particle_type) :: ptl_tmp
        real(dp) :: p_pivot
        integer :: first, last
        integer :: i, j

        p_pivot = ptls((first+last) / 2)%p
        i = first
        j = last
        do
            do while (ptls(i)%p < p_pivot)
                i = i + 1
            end do
            do while (p_pivot < ptls(j)%p)
                j = j - 1
            end do
            if (i >= j) exit
            ptl_tmp = ptls(i)
            ptls(i) = ptls(j)
            ptls(j) = ptl_tmp
            i = i + 1
            j = j - 1
        end do
        if (first < i-1) call quicksort_particle(ptls, first, i-1)
        if (j+1 < last)  call quicksort_particle(ptls, j+1, last)
    end subroutine quicksort_particle

    !---------------------------------------------------------------------------
    !< Initialize particle tracking
    !< Args:
    !<  nptl_selected: number of selected particles
    !---------------------------------------------------------------------------
    subroutine init_particle_tracking(nptl_selected)
        implicit none
        integer, intent(in) :: nptl_selected
        allocate(nsteps_tracked_ptls(nptl_selected))
        allocate(noffsets_tracked_ptls(nptl_selected))
        allocate(tags_selected_ptls(nptl_selected))
        nsteps_tracked_ptls = 0
        noffsets_tracked_ptls = 0
        tags_selected_ptls = 0
    end subroutine init_particle_tracking

    !---------------------------------------------------------------------------
    !< Free particle tracking
    !---------------------------------------------------------------------------
    subroutine free_particle_tracking
        implicit none
        deallocate(nsteps_tracked_ptls, noffsets_tracked_ptls)
        deallocate(tags_selected_ptls)
    end subroutine free_particle_tracking

    !---------------------------------------------------------------------------
    !< Select particles for particle tracking
    !< Args:
    !<  nptl: number of initial particles
    !<  nptl_selected: number of selected particles
    !<  nsteps_interval: save particle points every nsteps_interval
    !---------------------------------------------------------------------------
    subroutine select_particles_tracking(nptl, nptl_selected, nsteps_interval)
        implicit none
        integer, intent(in) :: nptl, nptl_selected, nsteps_interval
        integer :: iptl, i
        call quicksort_particle(ptls, 1, nptl_current)

        !< Select high-energy particles
        iptl = nptl_current-nptl_selected+1
        tags_selected_ptls = ptls(iptl:nptl_current)%tag
        nsteps_tracked_ptls = ptls(iptl:nptl_current)%nsteps_tracking

        !< Only part of the trajectory points are saved
        nsteps_tracked_ptls = (nsteps_tracked_ptls + nsteps_interval - 1) / nsteps_interval
        nsteps_tracked_ptls = nsteps_tracked_ptls + 1  ! Include the initial point

        do i = 2, nptl_selected
            noffsets_tracked_ptls(i) = noffsets_tracked_ptls(i-1) + &
                nsteps_tracked_ptls(i-1)
        enddo
        nsteps_tracked_tot = noffsets_tracked_ptls(nptl_selected) + &
                             nsteps_tracked_ptls(nptl_selected)

        !< Set particles to their initial values and set nptl_current to 0
        ptls%x = 0.0
        ptls%y = 0.0
        ptls%z = 0.0
        ptls%p = 0.0
        ptls%v = 0.0
        ptls%mu = 0.0
        ptls%weight = 0.0
        ptls%t = 0.0
        ptls%dt = 0.0
        ptls%split_times = 0
        ptls%count_flag = 0
        ptls%tag = 0
        ptls%nsteps_tracking = 0
        nptl_current = 0
        tag_max = 0
    end subroutine select_particles_tracking

    !---------------------------------------------------------------------------
    !< Initialize tracked particle points
    !< Args:
    !<  nptl_selected: number of selected particles
    !---------------------------------------------------------------------------
    subroutine init_tracked_particle_points(nptl_selected)
        implicit none
        integer, intent(in) :: nptl_selected
        allocate(ptl_traj_points(nsteps_tracked_tot))
        ptl_traj_points%x = 0.0
        ptl_traj_points%y = 0.0
        ptl_traj_points%z = 0.0
        ptl_traj_points%p = 0.0
        ptl_traj_points%v = 0.0
        ptl_traj_points%mu = 0.0
        ptl_traj_points%weight = 0.0
        ptl_traj_points%t = 0.0
        ptl_traj_points%dt = 0.0
        ptl_traj_points%split_times = 0
        ptl_traj_points%count_flag = 0
        ptl_traj_points%tag = 0
        ptl_traj_points%nsteps_tracking = 0
    end subroutine init_tracked_particle_points

    !---------------------------------------------------------------------------
    !< Make the flags of tracked particles negative, so they can be easily
    !< identified.
    !< Args:
    !<  nptl_selected: number of selected particles
    !---------------------------------------------------------------------------
    subroutine negative_particle_tags(nptl_selected)
        implicit none
        integer, intent(in) :: nptl_selected
        integer :: iptl, itag
        do iptl = nptl_current - nptl_inject + 1, nptl_current
            do itag = 1, nptl_selected
                if (ptls(iptl)%tag == tags_selected_ptls(itag)) then
                    ptls(iptl)%tag = -itag
                endif
            enddo
        enddo
    end subroutine negative_particle_tags

    !---------------------------------------------------------------------------
    !< Save the initial information of tracked particles
    !< Args:
    !<  nptl_selected: number of selected particles
    !---------------------------------------------------------------------------
    subroutine record_tracked_particle_init(nptl_selected)
        implicit none
        integer, intent(in) :: nptl_selected
        integer :: iptl, offset
        do iptl = nptl_current - nptl_inject + 1, nptl_current
            if (ptls(iptl)%tag < 0) then
                offset = noffsets_tracked_ptls(-ptls(iptl)%tag)
                ptl_traj_points(offset + 1) = ptls(iptl)
            endif
        enddo
    end subroutine record_tracked_particle_init

    !---------------------------------------------------------------------------
    !< Free tracked particle points
    !---------------------------------------------------------------------------
    subroutine free_tracked_particle_points
        implicit none
        deallocate(ptl_traj_points)
    end subroutine free_tracked_particle_points

    !---------------------------------------------------------------------------
    !< Save tracked particle points
    !< Args:
    !<  nptl_selected: number of selected particles
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine save_tracked_particle_points(nptl_selected, file_path)
        implicit none
        integer, intent(in) :: nptl_selected
        character(*), intent(in) :: file_path
        character(len=4) :: mrank
        character(len=128) :: fname
        integer :: fh, pos1
        write (mrank,'(i4.4)') mpi_rank
        fh = 41
        fname = trim(file_path)//'tracked_particle_points_'//mrank//'.dat'
        open(unit=fh, file=fname, access='stream', status='unknown', &
            form='unformatted', action='write')
        pos1 = 1
        write(fh, pos=pos1) nptl_selected
        pos1 = pos1 + sizeof(fp)
        write(fh, pos=pos1) nsteps_tracked_ptls
        pos1 = pos1 + nptl_selected * sizeof(fp)
        write(fh, pos=pos1) ptl_traj_points
        close(fh)
    end subroutine save_tracked_particle_points

    !---------------------------------------------------------------------------
    !< set minimum and maximum particle step
    !< Args:
    !<  dtf: time interval of the MHD fields
    !---------------------------------------------------------------------------
    subroutine set_dt_min_max(dtf)
        implicit none
        real(dp), intent(in) :: dtf
        dt_min = dt_min_rel * dtf
        dt_max = dt_max_rel * dtf
    end subroutine set_dt_min_max
end module particle_module

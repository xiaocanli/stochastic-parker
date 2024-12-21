!*******************************************************************************
!< Module of particle data and methods to inject, remove and push particles
!*******************************************************************************
module particle_module
    use constants, only: i1, i4, i8, sp, dp
    use simulation_setup_module, only: ndim_field
    use mhd_config_module, only: uniform_grid_flag, spherical_coord_flag
    use mhd_data_parallel, only: nfields, ngrads
    use mhd_data_parallel, only: xpos_local, ypos_local, zpos_local
    use mpi_module
    use omp_lib
    use hdf5
    implicit none
    private
    save
    public init_particles, free_particles, inject_particles_spatial_uniform, &
        read_particle_params, particle_mover, split_particle, &
        set_particle_datatype_mpi, free_particle_datatype_mpi, &
        init_particle_tracking, free_particle_tracking, dump_tracked_particles, &
        inject_particles_at_shock, inject_particles_at_large_jz, &
        inject_particles_at_large_absj, inject_particles_at_large_db2, &
        inject_particles_at_large_divv, inject_particles_at_large_rho, &
        set_dpp_params, set_duu_params, set_flags_params, set_drift_parameters, &
        set_flag_check_drift_2d, get_interp_paramters, &
        get_interp_paramters_spherical, read_particles, &
        save_particle_module_state, read_particle_module_state, &
        binarySearch_R

    public particle_type, ptls, escaped_ptls, &
        nptl_current, nptl_escaped, nptl_escaped_max, nptl_max, &
        spherical_coord_flag, leak, leak_negp, nptl_split, &
        pmin, pmax
    public COUNT_FLAG_INBOX, COUNT_FLAG_OTHERS, &
        COUNT_FLAG_ESCAPE_LX, COUNT_FLAG_ESCAPE_HX, &
        COUNT_FLAG_ESCAPE_LY, COUNT_FLAG_ESCAPE_HY, &
        COUNT_FLAG_ESCAPE_LZ, COUNT_FLAG_ESCAPE_HZ

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
        real(dp) :: padding             !< padding for safety
    end type particle_type

    real(dp) :: pmin  !< Minimum particle momentum
    real(dp) :: pmax  !< Maximum particle momentum

    integer, parameter :: COUNT_FLAG_INBOX  = 1  !< For in-box particles
    integer, parameter :: COUNT_FLAG_ESCAPE_LX = -1 !< Escaped from low-x boundary
    integer, parameter :: COUNT_FLAG_ESCAPE_HX = -2 !< Escaped from high-x boundary
    integer, parameter :: COUNT_FLAG_ESCAPE_LY = -3 !< Escaped from low-y boundary
    integer, parameter :: COUNT_FLAG_ESCAPE_HY = -4 !< Escaped from high-y boundary
    integer, parameter :: COUNT_FLAG_ESCAPE_LZ = -5 !< Escaped from low-z boundary
    integer, parameter :: COUNT_FLAG_ESCAPE_HZ = -6 !< Escaped from high-z boundary
    integer, parameter :: COUNT_FLAG_OTHERS = 0  !< For other particles

    integer :: particle_datatype_mpi
    type(particle_type), allocatable, dimension(:) :: ptls
    type(particle_type), allocatable, dimension(:, :) :: senders
    type(particle_type), allocatable, dimension(:, :) :: recvers
    integer, allocatable, dimension(:) :: nsenders, nrecvers
    !dir$ attributes align:128 :: ptls
    !dir$ attributes align:128 :: senders
    !dir$ attributes align:128 :: recvers

    integer :: nptl_current     !< Number of particles currently in the box
    integer :: nptl_old         !< Number of particles without receivers
    integer :: nptl_max         !< Maximum number of particles allowed
    integer :: nptl_split       !< Number of particles from splitting
    integer :: nptl_inject      !< Number of injected particles
    real(dp) :: leak            !< Leaking particles from boundary considering weight
    real(dp) :: leak_negp       !< Leaking particles with negative momentum

    ! tagging
    integer(i4) :: tag_max          !< Maximum particle tag

    real(dp) :: kpara0              !< kpara for particles with momentum p0
    real(dp) :: kret                !< The ratio of kpara to kperp
    integer :: momentum_dependency  !< kappa dependency on particle momentum
    integer :: mag_dependency       !< kappa dependency on magnetic field
    integer :: acc_region_flag      !< flag for whether to turn on acceleration in certain region
    real(dp) :: gamma_turb          !< turbulence spectral index
    real(dp) :: pindex              !< power index for the momentum dependency
    real(dp) :: p0    !< the standard deviation of the Gaussian distribution of momentum
    real(dp) :: b0    !< Initial magnetic field strength
    type kappa_type
        real(dp) :: knorm_para, knorm_perp    !< normalizations for spatial diffusion coefficient
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

    !< Pitch-angle
    real(dp), parameter :: mu_max=0.99  !< Avoiding parallel streaming
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
    logical :: track_particle_flag = .false. ! if init_particle_tracking is called, it will be .true.
    integer :: nptl_tracking   ! Number of particles to track
    integer :: split_times_max ! Maximum number of split times among these particles
    integer :: nsteps_tracking_max  ! Maximum number tracking steps backed dt_min_rel
    integer(i4), allocatable, dimension(:, :) :: tags_tracking
    type(particle_type), allocatable, dimension(:, :) :: particles_tracked

    interface read_ptl_element
        module procedure &
            read_integer1_element, read_integer4_element, &
            read_integer8_element, read_double_element
    end interface read_ptl_element

    interface write_tracked_ptl_element
        module procedure &
            write_integer1_tracked_element, write_integer4_tracked_element, &
            write_integer8_tracked_element, write_double_tracked_element
    end interface write_tracked_ptl_element

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
        ptls%origin = 0
        ptls%nsteps_tracked = 0
        ptls%nsteps_pushed = 0
        ptls%tag_injected = 0
        ptls%tag_splitted = 0
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
        senders%origin = 0
        senders%nsteps_tracked = 0
        senders%nsteps_pushed = 0
        senders%tag_injected = 0
        senders%tag_splitted = 0

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
        recvers%origin = 0
        recvers%nsteps_tracked = 0
        recvers%nsteps_pushed = 0
        recvers%tag_injected = 0
        recvers%tag_splitted = 0

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
        integer :: oldtypes(0:2), blockcounts(0:2)
        integer(KIND=MPI_ADDRESS_KIND) :: offsets(0:2)
        integer(i8) :: lb, extent
        ! Setup description of the 2 MPI_INTEGER1 fields.
        offsets(0) = 0
        oldtypes(0) = MPI_INTEGER1
        blockcounts(0) = 2
        ! Setup description of the 5 MPI_INTEGER4 fields.
        call MPI_TYPE_GET_EXTENT(MPI_INTEGER1, lb, extent, ierr)
        offsets(1) = blockcounts(0) * extent + offsets(0)
        oldtypes(1) = MPI_INTEGER4
        blockcounts(1) = 5
        ! Setup description of the 9 MPI_DOUBLE_PRECISION fields + 1 padding
        call MPI_TYPE_GET_EXTENT(MPI_INTEGER4, lb, extent, ierr)
        offsets(2) = blockcounts(1) * extent + offsets(1)
        oldtypes(2) = MPI_DOUBLE_PRECISION
        blockcounts(2) = 10
        ! Define structured type and commit it.
        call MPI_TYPE_CREATE_STRUCT(3, blockcounts, offsets, oldtypes, &
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
        integer, intent(in) :: nptl_current
        integer, intent(in) :: dist_flag, ct_mhd
        real(dp) :: r01, norm, fxp, ptmp, ftest, dt_mhd
        integer :: iptl_lo, iptl_hi
        ptls(nptl_current)%x = xpos
        ptls(nptl_current)%y = ypos
        ptls(nptl_current)%z = zpos
        if (dist_flag == 0) then
            ftest = 1.0
            fxp = 0.5
            do while (ftest > fxp)
                ptmp = (unif_01(0) * (pmax - pmin) + pmin) / p0  ! Need to normalized
                fxp = ptmp**2 * exp(-ptmp**2)
                ftest = unif_01(0) * 0.37  ! The maximum value is about 0.3679
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
        dt_mhd = tstamps_mhd(ct_mhd+1) - tstamps_mhd(ct_mhd)
        ptls(nptl_current)%t = tstamps_mhd(ct_mhd) + unif_01(0) * dt_mhd
        ptls(nptl_current)%dt = dt
        ptls(nptl_current)%split_times = 0
        ptls(nptl_current)%count_flag = COUNT_FLAG_INBOX
        ptls(nptl_current)%origin = mpi_rank
        ptls(nptl_current)%nsteps_tracked = 0
        ptls(nptl_current)%nsteps_pushed = 0
        ptls(nptl_current)%tag_injected = tag_max
        tag_max = tag_max + 1
        ptls(nptl_current)%tag_splitted = 1
        if (track_particle_flag .and. &
            is_particle_selected(ptls(nptl_current), iptl_lo, iptl_hi)) then
            ptls(nptl_current)%nsteps_tracked = 1
            ptls(nptl_current)%tag_injected = -ptls(nptl_current)%tag_injected
            ptls(nptl_current)%tag_splitted = -1
            particles_tracked(1, iptl_lo:iptl_hi) = ptls(nptl_current)
        endif
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
        integer, intent(in) :: nptl
        integer, intent(in) :: dist_flag, ct_mhd
        real(dp), intent(in) :: dt, power_index, particle_v0
        real(dp), intent(in), dimension(6) :: part_box
        integer :: i
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

        if (fconfig%nxg == fconfig%nxf .and. &
            fconfig%nyg == fconfig%nyf .and. &
            fconfig%nzg == fconfig%nzf) then
            ! each MPI rank reads all the data
            nptl_inject = nptl
            do i = 1, nptl_inject
                nptl_current = nptl_current + 1
                if (nptl_current > nptl_max) nptl_current = nptl_max
                inbox = .false.
                xtmp = unif_01(0) * (xmax_box - xmin_box) + xmin_box
                ytmp = unif_01(0) * (ymax_box - ymin_box) + ymin_box
                ztmp = unif_01(0) * (zmax_box - zmin_box) + zmin_box
                mu_tmp = mu_max * (2.0d0 * unif_01(0) - 1.0d0)
                call inject_one_particle(xtmp, ytmp, ztmp, nptl_current, &
                    dist_flag, particle_v0, mu_tmp, ct_mhd, dt, power_index)
            enddo
        else
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
                mu_tmp = mu_max * (2.0d0 * unif_01(0) - 1.0d0)
                call inject_one_particle(xtmp, ytmp, ztmp, nptl_current, &
                    dist_flag, particle_v0, mu_tmp, ct_mhd, dt, power_index)
            enddo
        endif


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
        integer, intent(in) :: nptl
        integer, intent(in) :: dist_flag, ct_mhd
        real(dp), intent(in) :: dt, power_index, particle_v0
        integer :: i
        integer :: iy, iz
        real(dp) :: xmin, ymin, xmax, ymax, zmin, zmax
        real(dp) :: ry, rz, dpy, dpz, shock_xpos
        real(dp) :: r01, norm, fxp, ptmp, ftest
        integer :: iptl_lo, iptl_hi
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
            ptls(nptl_current)%mu = mu_max * (2.0d0 * unif_01(0) - 1.0d0)
            ptls(nptl_current)%weight = 1.0d0
            ptls(nptl_current)%t = tstamps_mhd(ct_mhd)
            ptls(nptl_current)%dt = dt
            ptls(nptl_current)%split_times = 0
            ptls(nptl_current)%count_flag = COUNT_FLAG_INBOX
            ptls(nptl_current)%origin = mpi_rank
            ptls(nptl_current)%nsteps_tracked = 0
            ptls(nptl_current)%nsteps_pushed = 0
            ptls(nptl_current)%tag_injected = tag_max
            tag_max = tag_max + 1
            ptls(nptl_current)%tag_splitted = 1
            if (track_particle_flag .and. &
                is_particle_selected(ptls(nptl_current), iptl_lo, iptl_hi)) then
                ptls(nptl_current)%nsteps_tracked = 1
                ptls(nptl_current)%tag_injected = -ptls(nptl_current)%tag_injected
                ptls(nptl_current)%tag_splitted = -1
                particles_tracked(1, iptl_lo:iptl_hi) = ptls(nptl_current)
            endif
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
        integer, intent(in) :: nptl
        integer, intent(in) :: dist_flag, ct_mhd
        integer, intent(in) :: ncells_large_jz_norm
        logical, intent(in) :: inject_same_nptl
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
        integer :: i, ix, iy, iz
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
        if (inject_same_nptl) then
            call MPI_ALLREDUCE(ncells_large_jz, ncells_large_jz_g, 1, &
                MPI_INTEGER, MPI_SUM, mpi_sub_comm, ierr)
            ! The global sum may be a little different from nptl
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
                    rt = 0.0_dp
                    if (spherical_coord_flag) then
                        call get_interp_paramters_spherical(xtmp, ytmp, ztmp, pos, weights)
                    else
                        if (uniform_grid_flag) then
                            px = (xtmp - xmin) / dxm
                            py = (ytmp - ymin) / dym
                            pz = (ztmp - zmin) / dzm
                        else
                            ix = binarySearch_R(xpos_local, xtmp) - 2
                            iy = binarySearch_R(ypos_local, ytmp) - 2
                            iz = binarySearch_R(zpos_local, ztmp) - 2
                            px = (xtmp - xpos_local(ix)) / (xpos_local(ix+1) - xpos_local(ix))
                            py = (ytmp - ypos_local(iy)) / (ypos_local(iy+1) - ypos_local(iy))
                            pz = (ztmp - zpos_local(iz)) / (zpos_local(iz+1) - zpos_local(iz))
                        endif
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
            mu_tmp = mu_max * (2.0d0 * unif_01(0) - 1.0d0)
            call inject_one_particle(xtmp, ytmp, ztmp, nptl_current, &
                dist_flag, particle_v0, mu_tmp, ct_mhd, dt, power_index)
        enddo
        if (mpi_rank == master) then
            write(*, "(A)") "Finished injecting particles where jz is large"
        endif
    end subroutine inject_particles_at_large_jz

    !---------------------------------------------------------------------------
    !< Inject particles where current density absJ is large
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
    !<  absj_min: the minimum absj
    !<  ncells_large_absj_norm ! Normalization for the number of cells with large absj
    !<  part_box: box to inject particles
    !<  power_index: power-law index if dist_flag==2
    !---------------------------------------------------------------------------
    subroutine inject_particles_at_large_absj(nptl, dt, dist_flag, particle_v0, &
            ct_mhd, inject_same_nptl, absj_min, ncells_large_absj_norm, part_box, &
            power_index)
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: interp_fields
        use mhd_data_parallel, only: get_ncells_large_absj
        use random_number_generator, only: unif_01
        implicit none
        integer, intent(in) :: nptl
        integer, intent(in) :: dist_flag, ct_mhd
        integer, intent(in) :: ncells_large_absj_norm
        logical, intent(in) :: inject_same_nptl
        real(dp), intent(in) :: dt, absj_min, power_index, particle_v0
        real(dp), intent(in), dimension(6) :: part_box
        real(dp) :: xmin, ymin, zmin, xmax, ymax, zmax
        real(dp) :: xmin_box, ymin_box, zmin_box
        real(dp) :: xmax_box, ymax_box, zmax_box
        real(dp) :: xtmp, ytmp, ztmp, px, py, pz
        real(dp) :: dxm, dym, dzm
        real(dp) :: dbx_dy, dbx_dz, dby_dx, dby_dz, dbz_dx, dbz_dy
        real(dp) :: rt, absj, mu_tmp
        real(dp) :: bx, by, bz
        real(dp) :: ctheta, stheta, r
        real(dp), dimension(2) :: rands
        real(dp) :: r01, norm
        real(dp), dimension(nfields+ngrads) :: fields !< Fields at particle position
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights
        integer :: i, ix, iy, iz
        integer :: ncells_large_absj, ncells_large_absj_g
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
        !< with larger absj in the local domain will be different from mpi_rank
        !< to mpi_rank. That's why we need to redistribute the number of particles
        !< to inject.
        ncells_large_absj = get_ncells_large_absj(absj_min, spherical_coord_flag, part_box)
        if (inject_same_nptl) then
            call MPI_ALLREDUCE(ncells_large_absj, ncells_large_absj_g, 1, &
                MPI_INTEGER, MPI_SUM, mpi_sub_comm, ierr)
            ! The global sum may be a little different from nptl
            nptl_inject = int(nptl * mpi_sub_size * &
                (dble(ncells_large_absj) / dble(ncells_large_absj_g)))
        else
            nptl_inject = int(nptl * mpi_sub_size * &
                (dble(ncells_large_absj) / dble(ncells_large_absj_norm)))
        endif

        do i = 1, nptl_inject
            nptl_current = nptl_current + 1
            if (nptl_current > nptl_max) nptl_current = nptl_max
            absj = -2.0_dp
            do while (absj < absj_min)
                xtmp = unif_01(0) * (xmax - xmin) + xmin
                ytmp = unif_01(0) * (ymax - ymin) + ymin
                ztmp = unif_01(0) * (zmax - zmin) + zmin
                if (xtmp >= xmin_box .and. xtmp <= xmax_box .and. &
                    ytmp >= ymin_box .and. ytmp <= ymax_box .and. &
                    ztmp >= zmin_box .and. ztmp <= zmax_box) then
                    rt = 0.0_dp
                    if (spherical_coord_flag) then
                        call get_interp_paramters_spherical(xtmp, ytmp, ztmp, pos, weights)
                    else
                        if (uniform_grid_flag) then
                            px = (xtmp - xmin) / dxm
                            py = (ytmp - ymin) / dym
                            pz = (ztmp - zmin) / dzm
                        else
                            ix = binarySearch_R(xpos_local, xtmp) - 2
                            iy = binarySearch_R(ypos_local, ytmp) - 2
                            iz = binarySearch_R(zpos_local, ztmp) - 2
                            px = (xtmp - xpos_local(ix)) / (xpos_local(ix+1) - xpos_local(ix))
                            py = (ytmp - ypos_local(iy)) / (ypos_local(iy+1) - ypos_local(iy))
                            pz = (ztmp - zpos_local(iz)) / (zpos_local(iz+1) - zpos_local(iz))
                        endif
                        call get_interp_paramters(px, py, pz, pos, weights)
                    endif
                    call interp_fields(pos, weights, rt, fields)
                    dbx_dy = fields(nfields+14)
                    dbx_dz = fields(nfields+15)
                    dby_dx = fields(nfields+16)
                    dby_dz = fields(nfields+18)
                    dbz_dx = fields(nfields+19)
                    dbz_dy = fields(nfields+20)
                    if (spherical_coord_flag) then
                        bx = fields(5) ! r
                        by = fields(6) ! theta
                        bz = fields(7) ! phi
                        ctheta = cos(ytmp)
                        stheta = sin(ytmp)
                        r = xtmp
                        absj = sqrt(((ctheta*bz + stheta*dbz_dy - dby_dz) / (r*stheta))**2 + &
                                    ((dbx_dz/stheta - bz - r*dbz_dx) / r)**2 + &
                                    ((by + r*dby_dx - dbx_dy) / r)**2)
                    else
                        absj = sqrt((dby_dz - dbz_dy)**2 + &
                                    (dbz_dx - dbx_dz)**2 + &
                                    (dbx_dy - dby_dx)**2)
                    endif
                else ! not in part_box
                    absj = -3.0_dp
                endif
            enddo
            mu_tmp = mu_max * (2.0d0 * unif_01(0) - 1.0d0)
            call inject_one_particle(xtmp, ytmp, ztmp, nptl_current, &
                dist_flag, particle_v0, mu_tmp, ct_mhd, dt, power_index)
        enddo
        if (mpi_rank == master) then
            write(*, "(A)") "Finished injecting particles where absj is large"
        endif
    end subroutine inject_particles_at_large_absj

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
        integer, intent(in) :: nptl
        integer, intent(in) :: dist_flag, ct_mhd
        integer, intent(in) :: ncells_large_db2_norm
        logical, intent(in) :: inject_same_nptl
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
        real(dp), dimension(4) :: db2_slab, db2_2d
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights
        integer :: i, ix, iy, iz
        integer :: ncells_large_db2, ncells_large_db2_g
        !dir$ attributes align:32 :: db2_slab, db2_2d

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
        if (inject_same_nptl) then
            call MPI_ALLREDUCE(ncells_large_db2, ncells_large_db2_g, 1, &
                MPI_INTEGER, MPI_SUM, mpi_sub_comm, ierr)
            ! The global sum may be a little different from nptl
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
                    rt = 0.0_dp
                    if (spherical_coord_flag) then
                        call get_interp_paramters_spherical(xtmp, ytmp, ztmp, pos, weights)
                    else
                        if (uniform_grid_flag) then
                            px = (xtmp - xmin) / dxm
                            py = (ytmp - ymin) / dym
                            pz = (ztmp - zmin) / dzm
                        else
                            ix = binarySearch_R(xpos_local, xtmp) - 2
                            iy = binarySearch_R(ypos_local, ytmp) - 2
                            iz = binarySearch_R(zpos_local, ztmp) - 2
                            px = (xtmp - xpos_local(ix)) / (xpos_local(ix+1) - xpos_local(ix))
                            py = (ytmp - ypos_local(iy)) / (ypos_local(iy+1) - ypos_local(iy))
                            pz = (ztmp - zpos_local(iz)) / (zpos_local(iz+1) - zpos_local(iz))
                        endif
                        call get_interp_paramters(px, py, pz, pos, weights)
                    endif
                    call interp_magnetic_fluctuation(pos, weights, rt, db2_slab, db2_2d)
                    db2 = db2_slab(1)
                else ! not in part_box
                    db2 = -3.0_dp
                endif
            enddo
            mu_tmp = mu_max * (2.0d0 * unif_01(0) - 1.0d0)
            call inject_one_particle(xtmp, ytmp, ztmp, nptl_current, &
                dist_flag, particle_v0, mu_tmp, ct_mhd, dt, power_index)
        enddo
        if (mpi_rank == master) then
            write(*, "(A)") "Finished injecting particles where db2 is large"
        endif
    end subroutine inject_particles_at_large_db2

    !---------------------------------------------------------------------------
    !< Inject particles where the compression is negatively large
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
        integer, intent(in) :: nptl
        integer, intent(in) :: dist_flag, ct_mhd
        integer, intent(in) :: ncells_large_divv_norm
        logical, intent(in) :: inject_same_nptl
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
        integer :: i, ix, iy, iz
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
        if (inject_same_nptl) then
            call MPI_ALLREDUCE(ncells_large_divv, ncells_large_divv_g, 1, &
                MPI_INTEGER, MPI_SUM, mpi_sub_comm, ierr)
            ! The global sum may be a little different from nptl
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
                    rt = 0.0_dp
                    if (spherical_coord_flag) then
                        call get_interp_paramters_spherical(xtmp, ytmp, ztmp, pos, weights)
                    else
                        if (uniform_grid_flag) then
                            px = (xtmp - xmin) / dxm
                            py = (ytmp - ymin) / dym
                            pz = (ztmp - zmin) / dzm
                        else
                            ix = binarySearch_R(xpos_local, xtmp) - 2
                            iy = binarySearch_R(ypos_local, ytmp) - 2
                            iz = binarySearch_R(zpos_local, ztmp) - 2
                            px = (xtmp - xpos_local(ix)) / (xpos_local(ix+1) - xpos_local(ix))
                            py = (ytmp - ypos_local(iy)) / (ypos_local(iy+1) - ypos_local(iy))
                            pz = (ztmp - zpos_local(iz)) / (zpos_local(iz+1) - zpos_local(iz))
                        endif
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
            mu_tmp = mu_max * (2.0d0 * unif_01(0) - 1.0d0)
            call inject_one_particle(xtmp, ytmp, ztmp, nptl_current, &
                dist_flag, particle_v0, mu_tmp, ct_mhd, dt, power_index)
        enddo
        if (mpi_rank == master) then
            write(*, "(A)") "Finished injecting particles where divv is negatively large"
        endif
    end subroutine inject_particles_at_large_divv

    !---------------------------------------------------------------------------
    !< Inject particles where plasma density rho is large
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
    !<  rho_min: the minimum rho
    !<  ncells_large_rho_norm ! Normalization for the number of cells with large rho
    !<  part_box: box to inject particles
    !<  power_index: power-law index if dist_flag==2
    !---------------------------------------------------------------------------
    subroutine inject_particles_at_large_rho(nptl, dt, dist_flag, particle_v0, &
            ct_mhd, inject_same_nptl, rho_min, ncells_large_rho_norm, part_box, &
            power_index)
        use constants, only: pi
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: interp_fields
        use mhd_data_parallel, only: get_ncells_large_rho
        use random_number_generator, only: unif_01
        implicit none
        integer, intent(in) :: nptl
        integer, intent(in) :: dist_flag, ct_mhd
        integer, intent(in) :: ncells_large_rho_norm
        logical, intent(in) :: inject_same_nptl
        real(dp), intent(in) :: dt, rho_min, power_index, particle_v0
        real(dp), intent(in), dimension(6) :: part_box
        real(dp) :: xmin, ymin, zmin, xmax, ymax, zmax
        real(dp) :: xmin_box, ymin_box, zmin_box
        real(dp) :: xmax_box, ymax_box, zmax_box
        real(dp) :: xtmp, ytmp, ztmp, px, py, pz
        real(dp) :: dxm, dym, dzm
        real(dp) :: rt, rho, mu_tmp
        real(dp), dimension(2) :: rands
        real(dp) :: r01, norm
        real(dp), dimension(nfields+ngrads) :: fields !< Fields at particle position
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights
        integer :: i, ix, iy, iz
        integer :: ncells_large_rho, ncells_large_rho_g
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
        !< with larger rho in the local domain will be different from mpi_rank
        !< to mpi_rank. That's why we need to redistribute the number of particles
        !< to inject.
        ncells_large_rho = get_ncells_large_rho(rho_min, part_box)
        if (inject_same_nptl) then
            call MPI_ALLREDUCE(ncells_large_rho, ncells_large_rho_g, 1, &
                MPI_INTEGER, MPI_SUM, mpi_sub_comm, ierr)
            ! The global sum may be a little different from nptl
            nptl_inject = int(nptl * mpi_sub_size * &
                (dble(ncells_large_rho) / dble(ncells_large_rho_g)))
        else
            nptl_inject = int(nptl * mpi_sub_size * &
                (dble(ncells_large_rho) / dble(ncells_large_rho_norm)))
        endif

        do i = 1, nptl_inject
            nptl_current = nptl_current + 1
            if (nptl_current > nptl_max) nptl_current = nptl_max
            rho = 0.0_dp
            do while (rho < rho_min)
                xtmp = unif_01(0) * (xmax - xmin) + xmin
                ytmp = unif_01(0) * (ymax - ymin) + ymin
                ztmp = unif_01(0) * (zmax - zmin) + zmin
                if (xtmp >= xmin_box .and. xtmp <= xmax_box .and. &
                    ytmp >= ymin_box .and. ytmp <= ymax_box .and. &
                    ztmp >= zmin_box .and. ztmp <= zmax_box) then
                    rt = 0.0_dp
                    if (spherical_coord_flag) then
                        call get_interp_paramters_spherical(xtmp, ytmp, ztmp, pos, weights)
                    else
                        if (uniform_grid_flag) then
                            px = (xtmp - xmin) / dxm
                            py = (ytmp - ymin) / dym
                            pz = (ztmp - zmin) / dzm
                        else
                            ix = binarySearch_R(xpos_local, xtmp) - 2
                            iy = binarySearch_R(ypos_local, ytmp) - 2
                            iz = binarySearch_R(zpos_local, ztmp) - 2
                            px = (xtmp - xpos_local(ix)) / (xpos_local(ix+1) - xpos_local(ix))
                            py = (ytmp - ypos_local(iy)) / (ypos_local(iy+1) - ypos_local(iy))
                            pz = (ztmp - zpos_local(iz)) / (zpos_local(iz+1) - zpos_local(iz))
                        endif
                        call get_interp_paramters(px, py, pz, pos, weights)
                    endif
                    call interp_fields(pos, weights, rt, fields)
                    rho = fields(4)
                else ! not in part_box
                    rho = 0.0_dp
                endif
            enddo
            mu_tmp = mu_max * (2.0d0 * unif_01(0) - 1.0d0)
            call inject_one_particle(xtmp, ytmp, ztmp, nptl_current, &
                dist_flag, particle_v0, mu_tmp, ct_mhd, dt, power_index)
        enddo
        if (mpi_rank == master) then
            write(*, "(A)") "Finished injecting particles where rho is large"
        endif
    end subroutine inject_particles_at_large_rho

    !---------------------------------------------------------------------------
    !< Particle mover in one cycle
    !< Args:
    !<  t0: the starting time
    !<  dtf: the time interval of the MHD fields
    !<  nsteps_interval: save particle points every nsteps_interval
    !<  num_fine_steps: number of fine time steps
    !<  focused_transport: whether to the Focused Transport equation
    !<  nlgc: whether to use NLGC to evaluate kperp
    !<  kperp_kpara: ratio between kperp and kpara for initial particles and normalizations
    !---------------------------------------------------------------------------
    subroutine particle_mover_one_cycle(t0, dtf, nsteps_interval, &
            num_fine_steps, focused_transport, nlgc, kperp_kpara)
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: interp_fields, interp_magnetic_fluctuation, &
            interp_correlation_length
        use acc_region_surface, only: interp_acc_surface
        implicit none
        real(dp), intent(in) :: t0, dtf, kperp_kpara
        integer, intent(in) :: nsteps_interval, num_fine_steps
        logical, intent(in) :: focused_transport, nlgc
        real(dp) :: dxm, dym, dzm, xmin, xmax, ymin, ymax, zmin, zmax
        real(dp) :: dxl, dxu, dyl, dyu, dzl, dzu ! boundary cell sizes if non-uniform grid
        real(dp) :: xmin1, xmax1, ymin1, ymax1, zmin1, zmax1
        real(dp) :: deltax, deltay, deltaz, deltap, deltav, deltamu
        real(dp) :: dt_target, dt_fine
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights
        real(dp) :: px, py, pz, rt, dt_old
        integer :: i, nx, ny, nz, ix, iy, iz
        integer :: step, thread_id
        integer :: iptl_lo, iptl_hi
        type(particle_type) :: ptl
        type(kappa_type) :: kappa
        real(dp) :: surface_height1, surface_height2
        real(dp), dimension(nfields+ngrads) :: fields
        real(dp), dimension(4) :: db2_slab, db2_2d, lc_slab, lc_2d
        !dir$ attributes align:256 :: fields
        !dir$ attributes align:32 :: db2_slab
        !dir$ attributes align:32 :: db2_2d
        !dir$ attributes align:32 :: lc_slab
        !dir$ attributes align:32 :: lc_2d

        dt_fine = dtf / num_fine_steps
        dxm = mhd_config%dx
        dym = mhd_config%dy
        dzm = mhd_config%dz
        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz
        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax
        zmin = fconfig%zmin
        zmax = fconfig%zmax
        if (uniform_grid_flag) then
            xmin1 = xmin - dxm * 0.5
            xmax1 = xmax + dxm * 0.5
            ymin1 = ymin - dym * 0.5
            ymax1 = ymax + dym * 0.5
            zmin1 = zmin - dzm * 0.5
            zmax1 = zmax + dzm * 0.5
        else
            dxl = xpos_local(2) - xpos_local(1)
            dyl = ypos_local(2) - ypos_local(1)
            dzl = zpos_local(2) - zpos_local(1)
            dxu = xpos_local(nx+1) - xpos_local(nx)
            dyu = ypos_local(ny+1) - ypos_local(ny)
            dzu = zpos_local(nz+1) - zpos_local(nz)
            xmin1 = xmin - dxl * 0.5
            xmax1 = xmax + dxu * 0.5
            ymin1 = ymin - dyl * 0.5
            ymax1 = ymax + dyu * 0.5
            zmin1 = zmin - dzl * 0.5
            zmax1 = zmax + dzu * 0.5
        endif

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(ptl, kappa, &
        !$OMP& fields, db2_slab, db2_2d, lc_slab, lc_2d, &
        !$OMP& surface_height1, surface_height2, &
        !$OMP& deltax, deltay, deltaz, deltap, deltav, deltamu, &
        !$OMP& dt_target, pos, weights, ix, iy, iz, px, py, pz, rt, &
        !$OMP& step, thread_id, iptl_lo, iptl_hi)
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
                    if (uniform_grid_flag) then
                        px = (ptl%x - xmin) / dxm
                        py = (ptl%y - ymin) / dym
                        pz = (ptl%z - zmin) / dzm
                    else
                        ix = binarySearch_R(xpos_local, ptl%x) - 2
                        iy = binarySearch_R(ypos_local, ptl%y) - 2
                        iz = binarySearch_R(zpos_local, ptl%z) - 2
                        px = (ptl%x - xpos_local(ix)) / (xpos_local(ix+1) - xpos_local(ix))
                        py = (ptl%y - ypos_local(iy)) / (ypos_local(iy+1) - ypos_local(iy))
                        pz = (ptl%z - zpos_local(iz)) / (zpos_local(iz+1) - zpos_local(iz))
                    endif
                    rt = (ptl%t - t0) / dtf
                    if (spherical_coord_flag) then
                        call get_interp_paramters_spherical(ptl%x, ptl%y, ptl%z, pos, weights)
                    else
                        call get_interp_paramters(px, py, pz, pos, weights)
                    endif
                    call interp_fields(pos, weights, rt, fields)
                    if (deltab_flag) then
                        call interp_magnetic_fluctuation(pos, weights, rt, db2_slab, db2_2d)
                    endif
                    if (correlation_flag) then
                        call interp_correlation_length(pos, weights, rt, lc_slab, lc_2d)
                    endif
                    if (nlgc) then
                        call calc_spatial_diffusion_coefficients_nlgc(ptl, focused_transport, &
                            kperp_kpara, fields, db2_slab, db2_2d, lc_slab, lc_2d, kappa)
                    else
                        call calc_spatial_diffusion_coefficients(ptl, focused_transport, &
                            fields, db2_slab, lc_slab, kappa)
                    endif
                    if (focused_transport) then
                        if (ndim_field == 1) then
                            call push_particle_1d_ft(thread_id, rt, ptl, fields, &
                                db2_slab, lc_slab, kappa, .false., deltax, deltap, deltav, deltamu)
                        else if (ndim_field == 2) then
                            if (include_3rd_dim_in2d_flag) then
                                call push_particle_2d_include_3rd_ft(thread_id, rt, ptl, &
                                    fields, db2_slab, lc_slab, kappa, .false., &
                                    deltax, deltay, deltaz, deltap, deltav, deltamu)
                            else
                                call push_particle_2d_ft(thread_id, rt, ptl, fields, &
                                    db2_slab, lc_slab, kappa, .false., &
                                    deltax, deltay, deltap, deltav, deltamu)
                            endif
                        else
                            if (acc_by_surface_flag) then
                                call interp_acc_surface(pos, weights, rt, &
                                    surface_height1, surface_height2)
                            endif
                            call push_particle_3d_ft(thread_id, rt, surface_height1, &
                                surface_height2, ptl, fields, db2_slab, lc_slab, kappa, &
                                .false., deltax, deltay, deltaz, deltap, deltav, deltamu)
                        endif
                    else
                        if (ndim_field == 1) then
                            call push_particle_1d(thread_id, rt, ptl, fields, &
                                kappa, .false., deltax, deltap)
                        else if (ndim_field == 2) then
                            if (include_3rd_dim_in2d_flag) then
                                call push_particle_2d_include_3rd(thread_id, rt, ptl, &
                                    fields, kappa, .false., deltax, deltay, deltaz, deltap)
                            else
                                call push_particle_2d(thread_id, rt, ptl, fields, &
                                    kappa, .false., deltax, deltay, deltap)
                            endif
                        else
                            if (acc_by_surface_flag) then
                                call interp_acc_surface(pos, weights, rt, &
                                    surface_height1, surface_height2)
                            endif
                            call push_particle_3d(thread_id, rt, surface_height1, &
                                surface_height2, ptl, fields, kappa, .false., &
                                deltax, deltay, deltaz, deltap)
                        endif
                    endif

                    ! Number of steps the particle has been pushed
                    ptl%nsteps_pushed = mod(ptl%nsteps_pushed + 1, nsteps_interval)

                    ! Track particles
                    if (track_particle_flag) then
                        if (ptl%tag_splitted < 0 .and. ptl%nsteps_pushed == 0) then
                            call locate_particle(ptl, iptl_lo, iptl_hi)
                            ptl%nsteps_tracked = ptl%nsteps_tracked + 1
                            particles_tracked(ptl%nsteps_tracked, iptl_lo:iptl_hi) = ptl
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
                    dt_old = ptl%dt  ! save it for quick_check
                    ptl%dt = t0 + dt_target - ptl%t
                    if (ptl%dt > 0) then
                        if (ptl%tag_splitted < 0 .and. ptl%nsteps_pushed == 0) then
                            ! back one step if necessary
                            ptl%nsteps_tracked = ptl%nsteps_tracked - 1
                            ptl%nsteps_pushed = nsteps_interval - 2
                        else
                            ptl%nsteps_pushed = ptl%nsteps_pushed - 1
                        endif

                        if (uniform_grid_flag) then
                            px = (ptl%x - xmin) / dxm
                            py = (ptl%y - ymin) / dym
                            pz = (ptl%z - zmin) / dzm
                        else
                            ix = binarySearch_R(xpos_local, ptl%x) - 2
                            iy = binarySearch_R(ypos_local, ptl%y) - 2
                            iz = binarySearch_R(zpos_local, ptl%z) - 2
                            px = (ptl%x - xpos_local(ix)) / (xpos_local(ix+1) - xpos_local(ix))
                            py = (ptl%y - ypos_local(iy)) / (ypos_local(iy+1) - ypos_local(iy))
                            pz = (ptl%z - zpos_local(iz)) / (zpos_local(iz+1) - zpos_local(iz))
                        endif
                        rt = (ptl%t - t0) / dtf
                        if (spherical_coord_flag) then
                            call get_interp_paramters_spherical(&
                                ptl%x, ptl%y, ptl%z, pos, weights)
                        else
                            call get_interp_paramters(px, py, pz, pos, weights)
                        endif
                        call interp_fields(pos, weights, rt, fields)
                        if (deltab_flag) then
                            call interp_magnetic_fluctuation(pos, weights, rt, db2_slab, db2_2d)
                        endif
                        if (correlation_flag) then
                            call interp_correlation_length(pos, weights, rt, lc_slab, lc_2d)
                        endif
                        if (nlgc) then
                            call calc_spatial_diffusion_coefficients_nlgc(ptl, focused_transport, &
                                kperp_kpara, fields, db2_slab, db2_2d, lc_slab, lc_2d, kappa)
                        else
                            call calc_spatial_diffusion_coefficients(ptl, focused_transport, &
                                fields, db2_slab, lc_slab, kappa)
                        endif
                        if (focused_transport) then
                            if (ndim_field == 1) then
                                call push_particle_1d_ft(thread_id, rt, ptl, fields, &
                                    db2_slab, lc_slab, kappa, .true., deltax, deltap, deltav, deltamu)
                            else if (ndim_field == 2) then
                                if (include_3rd_dim_in2d_flag) then
                                    call push_particle_2d_include_3rd_ft(thread_id, rt, ptl, &
                                        fields, db2_slab, lc_slab, kappa, .true., &
                                        deltax, deltay, deltaz, deltap, deltav, deltamu)
                                else
                                    call push_particle_2d_ft(thread_id, rt, ptl, fields, &
                                        db2_slab, lc_slab, kappa, .true., &
                                        deltax, deltay, deltap, deltav, deltamu)
                                endif
                            else
                                if (acc_by_surface_flag) then
                                    call interp_acc_surface(pos, weights, rt, &
                                        surface_height1, surface_height2)
                                endif
                                call push_particle_3d_ft(thread_id, rt, surface_height1, surface_height2, &
                                    ptl, fields, db2_slab, lc_slab, kappa, .true., &
                                    deltax, deltay, deltaz, deltap, deltav, deltamu)
                            endif
                        else
                            if (ndim_field == 1) then
                                call push_particle_1d(thread_id, rt, ptl, fields, &
                                    kappa, .true., deltax, deltap)
                            else if (ndim_field == 2) then
                                if (include_3rd_dim_in2d_flag) then
                                    call push_particle_2d_include_3rd(thread_id, rt, ptl, &
                                        fields, kappa, .true., deltax, deltay, deltaz, deltap)
                                else
                                    call push_particle_2d(thread_id, rt, ptl, fields, &
                                        kappa, .true., deltax, deltay, deltap)
                                endif
                            else
                                if (acc_by_surface_flag) then
                                    call interp_acc_surface(pos, weights, rt, &
                                        surface_height1, surface_height2)
                                endif
                                call push_particle_3d(thread_id, rt, surface_height1, surface_height2, &
                                    ptl, fields, kappa, .true., deltax, deltay, deltaz, deltap)
                            endif
                        endif
                        ptl%nsteps_pushed = mod(ptl%nsteps_pushed + 1, nsteps_interval)

                        ! Track particles
                        if (track_particle_flag) then
                            if (ptl%tag_splitted < 0 .and. ptl%nsteps_pushed == 0) then
                                call locate_particle(ptl, iptl_lo, iptl_hi)
                                ptl%nsteps_tracked = ptl%nsteps_tracked + 1
                                particles_tracked(ptl%nsteps_tracked, iptl_lo:iptl_hi) = ptl
                            endif
                        endif

                    endif  ! if (ptl%dt > 0)

                    ptl%dt = dt_old  ! change it back for quick_check

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
    !<  nlgc: whether to use NLGC to evaluate kperp
    !<  kperp_kpara: ratio between kperp and kpara for initial particles and normalizations
    !<  nsteps_interval: save particle points every nsteps_interval
    !<  mhd_tframe: MHD time frame, starting from 1
    !<  num_fine_steps: number of fine time steps
    !<  dump_escaped_dist: whether to dump distributions of the escaped particles
    !---------------------------------------------------------------------------
    subroutine particle_mover(focused_transport, nlgc, kperp_kpara, &
            nsteps_interval, mhd_tframe, num_fine_steps, dump_escaped_dist)
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config, tstamps_mhd
        implicit none
        logical, intent(in) :: focused_transport, nlgc
        integer, intent(in) :: nsteps_interval
        integer, intent(in) :: mhd_tframe, num_fine_steps
        logical, intent(in) :: dump_escaped_dist
        real(dp), intent(in) :: kperp_kpara
        real(dp) :: dxm, dym, dzm, xmin, xmax, ymin, ymax, zmin, zmax
        real(dp) :: dxl, dxu, dyl, dyu, dzl, dzu ! boundary cell sizes if non-uniform grid
        real(dp) :: xmin1, xmax1, ymin1, ymax1, zmin1, zmax1
        integer :: local_flag, global_flag, ncycle, iptl
        integer :: nx, ny, nz
        logical :: all_particles_in_box
        real(dp) :: t0, dtf
        type(particle_type) :: ptl

        all_particles_in_box = .false.
        nptl_old = 0

        t0 = tstamps_mhd(mhd_tframe)
        dtf = tstamps_mhd(mhd_tframe+1) - t0
        call set_dt_min_max(dtf)

        dxm = mhd_config%dx
        dym = mhd_config%dy
        dzm = mhd_config%dz
        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz
        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax
        zmin = fconfig%zmin
        zmax = fconfig%zmax
        if (uniform_grid_flag) then
            xmin1 = xmin - dxm * 0.5
            xmax1 = xmax + dxm * 0.5
            ymin1 = ymin - dym * 0.5
            ymax1 = ymax + dym * 0.5
            zmin1 = zmin - dzm * 0.5
            zmax1 = zmax + dzm * 0.5
        else
            dxl = xpos_local(2) - xpos_local(1)
            dyl = ypos_local(2) - ypos_local(1)
            dzl = zpos_local(2) - zpos_local(1)
            dxu = xpos_local(nx+1) - xpos_local(nx)
            dyu = ypos_local(ny+1) - ypos_local(ny)
            dzu = zpos_local(nz+1) - ypos_local(nz)
            xmin1 = xmin - dxl * 0.5
            xmax1 = xmax + dxu * 0.5
            ymin1 = ymin - dyl * 0.5
            ymax1 = ymax + dyu * 0.5
            zmin1 = zmin - dzl * 0.5
            zmax1 = zmax + dzu * 0.5
        endif

        ncycle = 0
        local_flag = 0
        global_flag = 0

        ! Reset nsteps_tracked to 1 at the start of the MHD interval
        ! It is set to 1 instead of 0 because nsteps_tracked could be set to 1
        ! when injecting particles at the beginning of the MHD interval
        ptls(:nptl_current)%nsteps_tracked = 1

        do while (.not. all_particles_in_box)
            ncycle = ncycle + 1
            nsenders = 0
            nrecvers = 0
            if (nptl_old < nptl_current) then
                call particle_mover_one_cycle(t0, dtf, nsteps_interval, &
                    num_fine_steps, focused_transport, nlgc, kperp_kpara)
            else
                do iptl = 1, nptl_current
                    ptl = ptls(iptl)
                    if (ptl%p < 0.0 .and. ptl%count_flag /= COUNT_FLAG_INBOX) then
                        ptl%count_flag = COUNT_FLAG_OTHERS
                        leak_negp = leak_negp + ptl%weight
                    else
                        call particle_boundary_condition(ptl, xmin1, xmax1, &
                            ymin1, ymax1, zmin1, zmax1)
                    endif
                enddo
            endif
            call remove_particles(dump_escaped_dist)
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
        nsenders = 0
        nrecvers = 0
        do iptl = 1, nptl_current
            ptl = ptls(iptl)
            if (ptl%p < 0.0 .and. ptl%count_flag /= COUNT_FLAG_INBOX) then
                ptl%count_flag = COUNT_FLAG_OTHERS
                leak_negp = leak_negp + ptl%weight
            else
                call particle_boundary_condition(ptl, xmin, xmax, ymin, ymax, zmin, zmax)
            endif
            ptls(iptl) = ptl
        enddo
        call remove_particles(dump_escaped_dist)
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
                ptl%count_flag = COUNT_FLAG_ESCAPE_LX
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
                ptl%count_flag = COUNT_FLAG_ESCAPE_HX
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
                    ptl%count_flag = COUNT_FLAG_ESCAPE_LY
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
                    ptl%count_flag = COUNT_FLAG_ESCAPE_HY
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
                    ptl%count_flag = COUNT_FLAG_ESCAPE_LZ
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
                    ptl%count_flag = COUNT_FLAG_ESCAPE_HZ
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
        integer :: nsend, nrecv, iptl
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
    !<  db2_slab: turbulence variance and its gradients at particle position
    !<  lc_slab: turbulence correlation length and its gradients at particle position
    !<  kappa: kappa and related variables (For focused transport, kappa is
    !<         actually the perpendicular diffusion coefficients
    !---------------------------------------------------------------------------
    subroutine calc_spatial_diffusion_coefficients(ptl, focused_transport, &
            fields, db2_slab, lc_slab, kappa)
        implicit none
        type(particle_type), intent(in) :: ptl
        logical, intent(in) :: focused_transport
        real(dp), dimension(*), intent(in) :: fields
        real(dp), dimension(*), intent(in) :: db2_slab
        real(dp), dimension(*), intent(in) :: lc_slab
        type(kappa_type), intent(out) :: kappa
        real(dp) :: knorm, knorm_para, knorm_perp
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
        if (b < EPSILON(b)) then
            ib1 = 1.0
        else
            ib1 = 1.0_dp / b
        endif
        ib2 = ib1 * ib1
        ib3 = ib1 * ib2
        ib4 = ib2 * ib2

        kappa%knorm_para = 1.0_dp
        kappa%knorm_perp = 1.0_dp
        if (mag_dependency == 1) then
            kappa%knorm_para = kappa%knorm_para * b**(gamma_turb - 2.0_dp)
        endif

        ! Magnetic fluctuation dB^2/B^2
        if (deltab_flag) then
            kappa%knorm_para = kappa%knorm_para / db2_slab(1)
        endif

        ! Turbulence correlation length
        ! Make sure that lc is non-zero in the data file!!!
        if (correlation_flag) then
            kappa%knorm_para = kappa%knorm_para * lc_slab(1)**(gamma_turb - 1.0_dp)
        endif

        ! Momentum-dependent kappa
        if (momentum_dependency == 1) then
            knorm = kappa%knorm_para * (ptl%p / p0)**pindex
        else
            knorm = kappa%knorm_para
        endif

        kappa%knorm_perp = kappa%knorm_para
        kappa%kpara = kpara0 * knorm
        kappa%kperp = kappa%kpara * kret

        kappa%skpara = dsqrt(2.0 * kappa%kpara)
        kappa%skperp = dsqrt(2.0 * kappa%kperp)
        kappa%skpara_perp = dsqrt(2.0 * (kappa%kpara - kappa%kperp))

        if (ndim_field == 1) then
            dkdx = 0.0_dp
            if (mag_dependency == 1) then
                dkdx = db_dx * ib1 * (gamma_turb - 2.0_dp)
            endif
            if (deltab_flag) then
                dkdx = dkdx - db2_slab(2) / db2_slab(1)
            endif
            if (correlation_flag) then
                dkdx = dkdx + (gamma_turb - 1.0_dp) * lc_slab(2) / lc_slab(1)
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
                    dkdx = db_dx * ib1 * (gamma_turb - 2.0_dp)
                    dkdy = db_dy * ib1 * (gamma_turb - 2.0_dp)
                endif
                if (deltab_flag) then
                    dkdx = dkdx - db2_slab(2) / db2_slab(1)
                    dkdy = dkdy - db2_slab(3) / db2_slab(1)
                endif
                if (correlation_flag) then
                    dkdx = dkdx + (gamma_turb - 1.0_dp) * lc_slab(2) / lc_slab(1)
                    dkdy = dkdy + (gamma_turb - 1.0_dp) * lc_slab(3) / lc_slab(1)
                endif
                if (focused_transport) then
                    kpp = -kappa%kperp
                else
                    kpp = kappa%kpara - kappa%kperp
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
                    dkdx = db_dx * ib1 * (gamma_turb - 2.0_dp)
                    dkdy = db_dy * ib1 * (gamma_turb - 2.0_dp)
                endif
                if (deltab_flag) then
                    dkdx = dkdx - db2_slab(2) / db2_slab(1)
                    dkdy = dkdy - db2_slab(3) / db2_slab(1)
                endif
                if (correlation_flag) then
                    dkdx = dkdx + (gamma_turb - 1.0_dp) * lc_slab(2) / lc_slab(1)
                    dkdy = dkdy + (gamma_turb - 1.0_dp) * lc_slab(3) / lc_slab(1)
                endif
                if (focused_transport) then
                    kpp = -kappa%kperp
                else
                    kpp = kappa%kpara - kappa%kperp
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
                dkdx = db_dx * (gamma_turb - 2.0_dp)
                dkdy = db_dy * (gamma_turb - 2.0_dp)
                dkdz = db_dz * (gamma_turb - 2.0_dp)
            endif
            if (deltab_flag) then
                dkdx = dkdx - db2_slab(2) / db2_slab(1)
                dkdy = dkdy - db2_slab(3) / db2_slab(1)
                dkdz = dkdz - db2_slab(4) / db2_slab(1)
            endif
            if (correlation_flag) then
                dkdx = dkdx + (gamma_turb - 1.0_dp) * lc_slab(2) / lc_slab(1)
                dkdy = dkdy + (gamma_turb - 1.0_dp) * lc_slab(3) / lc_slab(1)
                dkdz = dkdz + (gamma_turb - 1.0_dp) * lc_slab(4) / lc_slab(1)
            endif
            if (focused_transport) then
                kpp = -kappa%kperp
            else
                kpp = kappa%kpara - kappa%kperp
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
   !< Calculate the spatial diffusion coefficients
   !< Args:
   !<  ptl: particle structure
   !<  focused_transport: whether to the focused transport equation
   !<  kperp_kpara: ratio between kperp and kpara for initial particles and normalizations
   !<  fields: fields and their gradients at particle position
   !<  db2_slab, db2_2d: turbulence variance and its gradients at particle position
   !<  lc_slab, lc_2d: turbulence correlation length and its gradients at particle position
   !<  kappa: kappa and related variables (For focused transport, kappa is
   !<         actually the perpendicular diffusion coefficients
   !---------------------------------------------------------------------------
   subroutine calc_spatial_diffusion_coefficients_nlgc(ptl, focused_transport, &
           kperp_kpara, fields, db2_slab, db2_2d, lc_slab, lc_2d, kappa)
       implicit none
       type(particle_type), intent(in) :: ptl
       logical, intent(in) :: focused_transport
       real(dp), dimension(*), intent(in) :: fields
       real(dp), dimension(*), intent(in) :: db2_slab, db2_2d
       real(dp), dimension(*), intent(in) :: lc_slab, lc_2d
       real(dp), intent(in) :: kperp_kpara
       type(kappa_type), intent(out) :: kappa
       real(dp) :: knorm_para, knorm_perp
       real(dp) :: bx, by, bz, b, ib1, ib2, ib3, ib4
       real(dp) :: dbx_dx, dby_dx, dbz_dx
       real(dp) :: dbx_dy, dby_dy, dbz_dy
       real(dp) :: dbx_dz, dby_dz, dbz_dz
       real(dp) :: db_dx, db_dy, db_dz
       real(dp) :: dkpara_dx, dkpara_dy, dkpara_dz
       real(dp) :: dkperp_dx, dkperp_dy, dkperp_dz
       real(dp) :: kpp

       bx = fields(5)
       by = fields(6)
       bz = fields(7)
       b = dsqrt(bx**2 + by**2 + bz**2)
       if (b < EPSILON(b)) then
           ib1 = 1.0
       else
           ib1 = 1.0_dp / b
       endif
       ib2 = ib1 * ib1
       ib3 = ib1 * ib2
       ib4 = ib2 * ib2

       kappa%knorm_para = 1.0_dp
       kappa%knorm_perp = 1.0_dp
       if (mag_dependency == 1) then
           kappa%knorm_para = kappa%knorm_para * b**(gamma_turb - 2.0_dp)
           kappa%knorm_perp = kappa%knorm_perp * b**((gamma_turb - 2.0_dp) / 3.0_dp)
       endif

       ! Magnetic fluctuation dB^2/B^2
       if (deltab_flag) then
           kappa%knorm_para = kappa%knorm_para / db2_slab(1)
           kappa%knorm_perp = kappa%knorm_perp * db2_slab(1)**(-1.0_dp/3.0_dp) * &
               db2_2d(1)**(2.0_dp/3.0_dp)
       endif

       ! Turbulence correlation length
       ! Make sure that lc is non-zero in the data file!!!
       if (correlation_flag) then
           kappa%knorm_para = kappa%knorm_para * lc_slab(1)**(gamma_turb - 1.0_dp)
           kappa%knorm_perp = kappa%knorm_perp * lc_slab(1)**((gamma_turb - 1.0_dp) / 3.0_dp) * &
               lc_2d(1)**(2.0_dp/3.0_dp)
       endif

       ! Momentum-dependent kappa
       if (momentum_dependency == 1) then
           knorm_para = kappa%knorm_para * (ptl%p / p0)**pindex
           knorm_perp = kappa%knorm_perp * (ptl%p / p0)**((5.0_dp-gamma_turb)/3.0_dp)
       else
           knorm_para = kappa%knorm_para
           knorm_perp = kappa%knorm_perp
       endif

       kappa%kpara = kpara0 * knorm_para
       kappa%kperp = kpara0 * kperp_kpara * knorm_perp * ptl%mu**2

       kappa%skpara = dsqrt(2.0 * kappa%kpara)
       kappa%skperp = dsqrt(2.0 * kappa%kperp)
       kappa%skpara_perp = dsqrt(2.0 * (kappa%kpara - kappa%kperp))

       if (ndim_field == 1) then
           dkpara_dx = 0.0_dp
           dkperp_dx = 0.0_dp
           if (mag_dependency == 1) then
               dkpara_dx = db_dx * ib1 * (gamma_turb - 2.0_dp)
               dkperp_dx = db_dx * ib1 * (gamma_turb - 2.0_dp) / 3.0_dp
           endif
           if (deltab_flag) then
               dkpara_dx = dkpara_dx - db2_slab(2) / db2_slab(1)
               dkperp_dx = dkperp_dx - db2_slab(2) / db2_slab(1) / 3.0_dp + &
                   2.0_dp * db2_2d(2) / db2_2d(1) / 3.0_dp
           endif
           if (correlation_flag) then
               dkpara_dx = dkpara_dx + (gamma_turb - 1.0_dp) * lc_slab(2) / lc_slab(1)
               dkperp_dx = dkperp_dx + (gamma_turb - 1.0_dp) * lc_slab(2) / lc_slab(1) / 3.0_dp + &
                   2.0_dp * lc_2d(2) / lc_2d(1) / 3.0_dp
           endif
           if (focused_transport) then
               kappa%dkxx_dx = kappa%kperp * dkperp_dx
               if (spherical_coord_flag) then
                   kappa%kxx = kappa%kperp
               endif
           else
               kappa%dkxx_dx = kappa%kpara * dkpara_dx
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
               dkpara_dx = 0.0_dp
               dkpara_dy = 0.0_dp
               dkpara_dz = 0.0_dp
               dkperp_dx = 0.0_dp
               dkperp_dy = 0.0_dp
               dkperp_dz = 0.0_dp
               if (mag_dependency == 1) then
                   dkpara_dx = db_dx * ib1 * (gamma_turb - 2.0_dp)
                   dkpara_dy = db_dy * ib1 * (gamma_turb - 2.0_dp)
                   dkperp_dx = db_dx * ib1 * (gamma_turb - 2.0_dp) / 3.0_dp
                   dkperp_dy = db_dy * ib1 * (gamma_turb - 2.0_dp) / 3.0_dp
               endif
               if (deltab_flag) then
                   dkpara_dx = dkpara_dx - db2_slab(2) / db2_slab(1)
                   dkpara_dy = dkpara_dy - db2_slab(3) / db2_slab(1)
                   dkperp_dx = dkperp_dx - db2_slab(2) / db2_slab(1) / 3.0_dp + &
                       2.0 * db2_2d(2) / db2_2d(1) / 3.0_dp
                   dkperp_dy = dkperp_dy - db2_slab(3) / db2_slab(1) / 3.0_dp + &
                       2.0 * db2_2d(3) / db2_2d(1) / 3.0_dp
               endif
               if (correlation_flag) then
                   dkpara_dx = dkpara_dx + (gamma_turb - 1.0_dp) * lc_slab(2) / lc_slab(1)
                   dkpara_dy = dkpara_dy + (gamma_turb - 1.0_dp) * lc_slab(3) / lc_slab(1)
                   dkperp_dx = dkperp_dx + (gamma_turb - 1.0_dp) * lc_slab(2) / lc_slab(1) / 3.0_dp + &
                       2.0_dp * lc_2d(2) / lc_2d(1) / 3.0_dp
                   dkperp_dy = dkperp_dy + (gamma_turb - 1.0_dp) * lc_slab(3) / lc_slab(1) / 3.0_dp + &
                       2.0_dp * lc_2d(3) / lc_2d(1) / 3.0_dp
               endif
               if (focused_transport) then
                   kpp = -kappa%kperp
               else
                   kpp = kappa%kpara - kappa%kperp
               endif
               kappa%dkxx_dx = kappa%kperp*dkperp_dx + &
                   (kappa%kpara*dkpara_dx - kappa%kperp*dkperp_dx)*bx**2*ib2 + &
                   2.0*kpp*bx*(dbx_dx*b-bx*db_dx)*ib3
               kappa%dkyy_dy = kappa%kperp*dkperp_dy + &
                   (kappa%kpara*dkpara_dy - kappa%kperp*dkperp_dy)*by**2*ib2 + &
                   2.0*kpp*by*(dby_dy*b-by*db_dy)*ib3
               kappa%dkzz_dz = kappa%kperp*dkperp_dz + &
                   (kappa%kpara*dkpara_dz - kappa%kperp*dkperp_dz)*bz**2*ib2 + &
                   2.0*kpp*bz*(dbz_dz*b-bz*db_dz)*ib3
               kappa%dkxy_dx = (kappa%kpara*dkpara_dx - kappa%kperp*dkperp_dx)*bx*by*ib2 + kpp * &
                   ((dbx_dx*by+bx*dby_dx)*ib2 - 2.0*bx*by*db_dx*ib3)
               kappa%dkxy_dy = (kappa%kpara*dkpara_dy - kappa%kperp*dkperp_dy)*bx*by*ib2 + kpp * &
                   ((dbx_dy*by+bx*dby_dy)*ib2 - 2.0*bx*by*db_dy*ib3)
               kappa%dkxz_dx = (kappa%kpara*dkpara_dx - kappa%kperp*dkperp_dx)*bx*bz*ib2 + kpp * &
                   ((dbx_dx*bz+bx*dbz_dx)*ib2 - 2.0*bx*bz*db_dx*ib3)
               kappa%dkxz_dz = (kappa%kpara*dkpara_dz - kappa%kperp*dkperp_dz)*bx*bz*ib2 + kpp * &
                   ((dbx_dz*bz+bx*dbz_dz)*ib2 - 2.0*bx*bz*db_dz*ib3)
               kappa%dkyz_dy = (kappa%kpara*dkpara_dy - kappa%kperp*dkperp_dy)*by*bz*ib2 + kpp * &
                   ((dby_dy*bz+by*dbz_dy)*ib2 - 2.0*by*bz*db_dy*ib3)
               kappa%dkyz_dz = (kappa%kpara*dkpara_dz - kappa%kperp*dkperp_dz)*by*bz*ib2 + kpp * &
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
               dkpara_dx = 0.0_dp
               dkpara_dy = 0.0_dp
               dkperp_dx = 0.0_dp
               dkperp_dy = 0.0_dp
               if (mag_dependency == 1) then
                   dkpara_dx = db_dx * ib1 * (gamma_turb - 2.0_dp)
                   dkpara_dy = db_dy * ib1 * (gamma_turb - 2.0_dp)
                   dkperp_dx = db_dx * ib1 * (gamma_turb - 2.0_dp) / 3.0_dp
                   dkperp_dy = db_dy * ib1 * (gamma_turb - 2.0_dp) / 3.0_dp
               endif
               if (deltab_flag) then
                   dkpara_dx = dkpara_dx - db2_slab(2) / db2_slab(1)
                   dkpara_dy = dkpara_dy - db2_slab(3) / db2_slab(1)
                   dkperp_dx = dkperp_dx - db2_slab(2) / db2_slab(1) / 3.0_dp + &
                       2.0 * db2_2d(2) / db2_2d(1) / 3.0_dp
                   dkperp_dy = dkperp_dy - db2_slab(3) / db2_slab(1) / 3.0_dp + &
                       2.0 * db2_2d(3) / db2_2d(1) / 3.0_dp
               endif
               if (correlation_flag) then
                   dkpara_dx = dkpara_dx + (gamma_turb - 1.0_dp) * lc_slab(2) / lc_slab(1)
                   dkpara_dy = dkpara_dy + (gamma_turb - 1.0_dp) * lc_slab(3) / lc_slab(1)
                   dkperp_dx = dkperp_dx + (gamma_turb - 1.0_dp) * lc_slab(2) / lc_slab(1) / 3.0_dp + &
                       2.0_dp * lc_2d(2) / lc_2d(1) / 3.0_dp
                   dkperp_dy = dkperp_dy + (gamma_turb - 1.0_dp) * lc_slab(3) / lc_slab(1) / 3.0_dp + &
                       2.0_dp * lc_2d(3) / lc_2d(1) / 3.0_dp
               endif
               if (focused_transport) then
                   kpp = -kappa%kperp
               else
                   kpp = kappa%kpara - kappa%kperp
               endif
               kappa%dkxx_dx = kappa%kperp*dkperp_dx + &
                   (kappa%kpara*dkpara_dx - kappa%kperp*dkperp_dx)*bx**2*ib2 + &
                   2.0*kpp*bx*(dbx_dx*b-bx*db_dx)*ib3
               kappa%dkyy_dy = kappa%kperp*dkperp_dy + &
                   (kappa%kpara*dkpara_dy - kappa%kperp*dkperp_dy)*by**2*ib2 + &
                   2.0*kpp*by*(dby_dy*b-by*db_dy)*ib3
               kappa%dkxy_dx = (kappa%kpara*dkpara_dx - kappa%kperp*dkperp_dx)*bx*by*ib2 + kpp * &
                   ((dbx_dx*by+bx*dby_dx)*ib2 - 2.0*bx*by*db_dx*ib3)
               kappa%dkxy_dy = (kappa%kpara*dkpara_dy - kappa%kperp*dkperp_dy)*bx*by*ib2 + kpp * &
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
           dkpara_dx = 0.0_dp
           dkpara_dy = 0.0_dp
           dkpara_dz = 0.0_dp
           dkperp_dx = 0.0_dp
           dkperp_dy = 0.0_dp
           dkperp_dz = 0.0_dp
           if (mag_dependency == 1) then
               dkpara_dx = db_dx * ib1 * (gamma_turb - 2.0_dp)
               dkpara_dy = db_dy * ib1 * (gamma_turb - 2.0_dp)
               dkpara_dz = db_dz * ib1 * (gamma_turb - 2.0_dp)
               dkperp_dx = db_dx * ib1 * (gamma_turb - 2.0_dp) / 3.0_dp
               dkperp_dy = db_dy * ib1 * (gamma_turb - 2.0_dp) / 3.0_dp
               dkperp_dz = db_dz * ib1 * (gamma_turb - 2.0_dp) / 3.0_dp
           endif
           if (deltab_flag) then
               dkpara_dx = dkpara_dx - db2_slab(2) / db2_slab(1)
               dkpara_dy = dkpara_dy - db2_slab(3) / db2_slab(1)
               dkpara_dz = dkpara_dz - db2_slab(4) / db2_slab(1)
               dkperp_dx = dkperp_dx - db2_slab(2) / db2_slab(1) / 3.0_dp + &
                   2.0 * db2_2d(2) / db2_2d(1) / 3.0_dp
               dkperp_dy = dkperp_dy - db2_slab(3) / db2_slab(1) / 3.0_dp + &
                   2.0 * db2_2d(3) / db2_2d(1) / 3.0_dp
               dkperp_dz = dkperp_dz - db2_slab(4) / db2_slab(1) / 3.0_dp + &
                   2.0 * db2_2d(4) / db2_2d(1) / 3.0_dp
           endif
           if (correlation_flag) then
               dkpara_dx = dkpara_dx + (gamma_turb - 1.0_dp) * lc_slab(2) / lc_slab(1)
               dkpara_dy = dkpara_dy + (gamma_turb - 1.0_dp) * lc_slab(3) / lc_slab(1)
               dkpara_dz = dkpara_dz + (gamma_turb - 1.0_dp) * lc_slab(4) / lc_slab(1)
               dkperp_dx = dkperp_dx + (gamma_turb - 1.0_dp) * lc_slab(2) / lc_slab(1) / 3.0_dp + &
                   2.0_dp * lc_2d(2) / lc_2d(1) / 3.0_dp
               dkperp_dy = dkperp_dy + (gamma_turb - 1.0_dp) * lc_slab(3) / lc_slab(1) / 3.0_dp + &
                   2.0_dp * lc_2d(3) / lc_2d(1) / 3.0_dp
               dkperp_dz = dkperp_dz + (gamma_turb - 1.0_dp) * lc_slab(4) / lc_slab(1) / 3.0_dp + &
                   2.0_dp * lc_2d(4) / lc_2d(1) / 3.0_dp
           endif
           if (focused_transport) then
               kpp = -kappa%kperp
           else
               kpp = kappa%kpara - kappa%kperp
           endif
           kappa%dkxx_dx = kappa%kperp*dkperp_dx + &
               (kappa%kpara*dkpara_dx - kappa%kperp*dkperp_dx)*bx**2*ib2 + &
               2.0*kpp*bx*(dbx_dx*b-bx*db_dx)*ib3
           kappa%dkyy_dy = kappa%kperp*dkperp_dy + &
               (kappa%kpara*dkpara_dy - kappa%kperp*dkperp_dy)*by**2*ib2 + &
               2.0*kpp*by*(dby_dy*b-by*db_dy)*ib3
           kappa%dkzz_dz = kappa%kperp*dkperp_dz + &
               (kappa%kpara*dkpara_dz - kappa%kperp*dkperp_dz)*bz**2*ib2 + &
               2.0*kpp*bz*(dbz_dz*b-bz*db_dz)*ib3
           kappa%dkxy_dx = (kappa%kpara*dkpara_dx - kappa%kperp*dkperp_dx)*bx*by*ib2 + kpp * &
               ((dbx_dx*by+bx*dby_dx)*ib2 - 2.0*bx*by*db_dx*ib3)
           kappa%dkxy_dy = (kappa%kpara*dkpara_dy - kappa%kperp*dkperp_dy)*bx*by*ib2 + kpp * &
               ((dbx_dy*by+bx*dby_dy)*ib2 - 2.0*bx*by*db_dy*ib3)
           kappa%dkxz_dx = (kappa%kpara*dkpara_dx - kappa%kperp*dkperp_dx)*bx*bz*ib2 + kpp * &
               ((dbx_dx*bz+bx*dbz_dx)*ib2 - 2.0*bx*bz*db_dx*ib3)
           kappa%dkxz_dz = (kappa%kpara*dkpara_dz - kappa%kperp*dkperp_dz)*bx*bz*ib2 + kpp * &
               ((dbx_dz*bz+bx*dbz_dz)*ib2 - 2.0*bx*bz*db_dz*ib3)
           kappa%dkyz_dy = (kappa%kpara*dkpara_dy - kappa%kperp*dkperp_dy)*by*bz*ib2 + kpp * &
               ((dby_dy*bz+by*dbz_dy)*ib2 - 2.0*by*bz*db_dy*ib3)
           kappa%dkyz_dz = (kappa%kpara*dkpara_dz - kappa%kperp*dkperp_dz)*by*bz*ib2 + kpp * &
               ((dby_dz*bz+by*dbz_dz)*ib2 - 2.0*by*bz*db_dz*ib3)
           kappa%kxx = kappa%kperp + kpp * bx * bx * ib2
           kappa%kyy = kappa%kperp + kpp * by * by * ib2
           kappa%kzz = kappa%kperp + kpp * bz * bz * ib2
           kappa%kxy = kpp * bx * by * ib2
           kappa%kxz = kpp * bx * bz * ib2
           kappa%kyz = kpp * by * bz * ib2
       endif
   end subroutine calc_spatial_diffusion_coefficients_nlgc

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
        real(sp) :: temp
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
            gamma_turb = get_variable(fh, 'gamma_turb', '=')
            pindex = 3.0_dp - gamma_turb
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
        call MPI_BCAST(gamma_turb, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
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
    !< Calculate momentum diffusion due to wave scattering
    !< Args:
    !<  rho: plasma density
    !<  b: magnetic field strength
    !<  kpara: parallel diffusion coefficient
    !<  ptl: one particle
    !<  dp_dt: 1st-order momentum change rate
    !<  dpp: momentum diffusion
    !---------------------------------------------------------------------------
    subroutine calc_dpp_wave_scattering(rho, b, kpara, ptl, dp_dt, dpp)
        use mhd_config_module, only: mhd_config
        implicit none
        real(dp), intent(in) :: rho, b, kpara
        type(particle_type), intent(in) :: ptl
        real(dp), intent(inout) :: dp_dt, dpp
        real(dp) :: va ! Alfven speed

        va = b / dsqrt(rho)
        if (momentum_dependency == 1) then
            dp_dt = dp_dt + (8*ptl%p / (27*kpara)) * va**2
        else
            dp_dt = dp_dt + (4*ptl%p / (9*kpara)) * va**2
        endif
        dpp = dpp + (ptl%p * va)**2 / (9*kpara)
    end subroutine calc_dpp_wave_scattering

    !---------------------------------------------------------------------------
    !< Calculate momentum diffusion due to flow shear
    !< Args:
    !<  b: magnetic field strength
    !<  bx, by, bz: magnetic field components
    !<  knorm_para: normalization for the parallel diffusion coefficients
    !<  sigmaxx, sigmayy, sigmazz: shear tensor diagonal components
    !<  sigmaxy, sigmaxz, sigmayz: shear tensor other components
    !<  ptl: one particle
    !<  dp_dt: 1st-order momentum change rate
    !<  dpp: momentum diffusion
    !---------------------------------------------------------------------------
    subroutine calc_dpp_flow_shear(b, bx, by, bz, knorm_para, &
            sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, sigmayz, ptl, dp_dt, dpp)
        use mhd_config_module, only: mhd_config
        implicit none
        real(dp), intent(in) :: b, bx, by, bz
        real(dp), intent(in) :: knorm_para
        real(dp), intent(in) :: sigmaxx, sigmayy, sigmazz
        real(dp), intent(in) :: sigmaxy, sigmaxz, sigmayz
        type(particle_type), intent(in) :: ptl
        real(dp), intent(inout) :: dp_dt, dpp
        real(dp) :: ib, bbsigma, gshear

        if (weak_scattering) then
            if (b < EPSILON(b)) then
                ib = 0.0
            else
                ib = 1.0 / b
            endif
            bbsigma = sigmaxx * bx**2 + sigmayy * by**2 + sigmazz * bz**2 + &
                2.0 * (sigmaxy * bx * by + sigmaxz * bx * bz + sigmayz * by * bz)
            bbsigma = bbsigma * ib * ib
            gshear = bbsigma**2 / 5
        else
            gshear = 2 * (sigmaxx**2 + sigmayy**2 + sigmazz**2 + &
                2 * (sigmaxy**2 + sigmaxz**2 + sigmayz**2)) / 15
        endif
        if (gshear > 0.0d0) then
            dp_dt = dp_dt + (2 + pindex) * gshear * tau0 * knorm_para * &
                ptl%p**(pindex-1) * p0**(2.0-pindex)
            dpp = dpp + gshear * tau0 * knorm_para * &
                ptl%p**pindex * p0**(2.0-pindex)
        endif
    end subroutine calc_dpp_flow_shear

    !---------------------------------------------------------------------------
    !< Push particle for a single step for a 1D simulation
    !< Args:
    !<  thread_id: thread ID staring from 0
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  ptl: particle structure
    !<  fields: fields and their gradients at particle position
    !<  kappa: kappa and related variables
    !<  fixed_dt: whether to fix the time step
    !<  deltax, deltap: the change of x and p in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_1d(thread_id, rt, ptl, fields, kappa, fixed_dt, deltax, deltap)
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        logical, intent(in) :: fixed_dt
        real(dp), intent(in) :: rt
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltap
        real(dp) :: sdt, dvx_dx, divv
        real(dp) :: b, ib, bx, by, bz, vx, px, rt1
        real(dp) :: dxm, rho, ran1, sqrt3
        real(dp) :: rands(2)
        real(dp) :: sigmaxx, sigmayy, sigmazz ! shear tensor
        real(dp) :: dx_dt, dp_dt, dpp
        integer :: ix
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights

        vx = fields(1)
        bx = fields(5)
        by = fields(6)
        bz = fields(7)
        b = dsqrt(bx**2 + by**2 + bz**2)
        dvx_dx = fields(nfields+1)
        if (uniform_grid_flag) then
            dxm = mhd_config%dx
        else
            ix = binarySearch_R(xpos_local, ptl%x) - 2
            dxm = xpos_local(ix+1) - xpos_local(ix)
        endif

        deltax = 0.0d0
        deltap = 0.0d0

        !< Calculate the 1st-order and 2nd-order rates
        if (spherical_coord_flag) then
            dx_dt = vx + kappa%dkxx_dx + 2.0*kappa%kxx/ptl%x
            divv = dvx_dx + 2.0*vx/ptl%x
        else
            dx_dt = vx + kappa%dkxx_dx
            divv = dvx_dx
        endif

        dp_dt = -ptl%p * divv / 3.0d0

        !< Momentum diffusion due to wave scattering
        dpp = 0.0d0
        if (dpp_wave_flag) then
            rho = fields(4)
            call calc_dpp_wave_scattering(rho, b, kappa%kpara, ptl, dp_dt, dpp)
        endif

        !< Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            sigmaxx = dvx_dx - divv / 3
            sigmayy = -divv / 3
            sigmazz = -divv / 3
            call calc_dpp_flow_shear(b, bx, by, bz, kappa%knorm_para, &
                sigmaxx, sigmayy, sigmazz, 0.0d0, 0.0d0, 0.0d0, ptl, dp_dt, dpp)
        endif

        !< Set the time step
        if (.not. fixed_dt) then
            if ((dx_dt .ne. 0.0d0) .and. (dp_dt .ne. 0.0d0)) then
                if (kappa%skperp > 0.0_dp) then
                    ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                 (kappa%skperp/dx_dt)**2, &
                                 0.1*ptl%p/abs(dp_dt))
                else
                    ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                 (kappa%skpara/dx_dt)**2, &
                                 0.1*ptl%p/abs(dp_dt))
                endif
            else
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too small
            if (ptl%dt .lt. dt_min) then
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too large
            if (ptl%dt .gt. dt_max) then
                ptl%dt = dt_max
            endif
        endif

        !< Update the particle
        sdt = dsqrt(ptl%dt)
        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltax = dx_dt * ptl%dt + ran1*kappa%skpara*sdt
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltap = dp_dt * ptl%dt + ran1*dsqrt(2*dpp)*sdt

        ptl%x = ptl%x + deltax
        ptl%t = ptl%t + ptl%dt

        if (acc_region_flag == 1) then
            if (particle_in_acceleration_region(ptl)) then
                ptl%p = ptl%p + deltap
            else
                deltap = 0.0
            endif
        else
            ptl%p = ptl%p + deltap
        endif

        !< Make sure the momentum is not too small,
        !< so it might not deviate from the diffusion limit.
        if (ptl%p < 0.25 * p0) then
            ptl%p = ptl%p - deltap
            deltap = 0.25 * p0 - ptl%p
            ptl%p = 0.25 * p0
        endif
    end subroutine push_particle_1d

    !---------------------------------------------------------------------------
    !< Calculate pitch-angle diffusion coefficient
    !---------------------------------------------------------------------------
    subroutine calc_duu(ptl, b, db2_slab, lc_slab, div_bnorm, divv, &
            bb_gradv, bv_gradv, mu2, dmu_dt, duu, duu_du)
        implicit none
        type(particle_type), intent(in) :: ptl
        real(dp), intent(in) :: div_bnorm, divv, bb_gradv, bv_gradv, mu2
        real(dp), intent(in) :: b, db2_slab, lc_slab
        real(dp), intent(out) :: dmu_dt, duu, duu_du
        real(dp) :: dtmp, h0, duu_norm
        dmu_dt = ptl%v * div_bnorm + ptl%mu * divv - &
            3 * ptl%mu * bb_gradv - 2 * bv_gradv / ptl%v
        dmu_dt = dmu_dt * (1-mu2) * 0.5
        h0 = 0.2  ! Enabling pitch-angle scattering around mu=0
        dtmp = abs(ptl%mu)**(gamma_turb-1) + h0
        duu = duu0 * (1-mu2) * dtmp
        if (ptl%mu .gt. 0.0d0) then
            duu_du = duu0 * (-2*ptl%mu * dtmp + &
                             (1-mu2) * abs(ptl%mu)**(gamma_turb-2))
        else if (ptl%mu .lt. 0.0d0) then
            duu_du = duu0 * (-2*ptl%mu * dtmp - &
                             (1-mu2) * abs(ptl%mu)**(gamma_turb-2))
        else
            duu_du = 0.0d0
        endif
        duu_norm = 1.0_dp
        if (mag_dependency == 1) then
            duu_norm = duu_norm * b**(2.0_dp - gamma_turb)
        endif
        if (deltab_flag) then
            duu_norm = duu_norm * db2_slab
        endif
        if (correlation_flag) then
            duu_norm = duu_norm * lc_slab**(1.0 - gamma_turb)
        endif
        if (momentum_dependency == 1) then
            duu_norm = duu_norm * (ptl%p / p0)**(gamma_turb-1)
        endif
        duu_du = duu_du * duu_norm
        duu = duu * duu_norm
        dmu_dt = dmu_dt + duu_du
    end subroutine calc_duu

    !---------------------------------------------------------------------------
    !< Push particle for a single step for a 1D simulation for focused transport
    !< Args:
    !<  thread_id: thread ID staring from 0
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  ptl: particle structure
    !<  fields: fields and their gradients at particle position
    !<  db2_slab: turbulence variance for the slab component
    !<  lc_slab: turbulence correlation length for the slab component
    !<  kappa: kappa and related variables
    !<  fixed_dt: whether to fix the time step
    !<  deltax: the change of x in this step
    !<  deltap: the change of p in this step
    !<  deltav: the change of v in this step
    !<  deltamu: the change of mu in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_1d_ft(thread_id, rt, ptl, fields, db2_slab, lc_slab, &
            kappa, fixed_dt, deltax, deltap, deltav, deltamu)
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        logical, intent(in) :: fixed_dt
        real(dp), intent(in) :: rt
        real(dp), dimension(*), intent(in) :: db2_slab, lc_slab
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltap, deltav, deltamu
        real(dp) :: xtmp
        real(dp) :: sdt, dvx_dx, dvy_dx, dvz_dx
        real(dp) :: divv, bb_gradv, bv_gradv
        real(dp) :: b, ib, bx, by, bz, rt1
        real(dp) :: vx, vy, vz, cot_theta, ir
        real(dp) :: mu2, acc_rate, duu, dmu_dt, duu_du, dtmp
        real(dp) :: dp_dpp, div_bnorm, db_dx, h0
        real(dp) :: dxm, rho, ran1, sqrt3
        real(dp) :: rands(2)
        real(dp) :: sigmaxx, sigmayy, sigmazz ! shear tensor
        real(dp) :: dx_dt, dp_dt, dv_dt, dpp
        integer :: ix
        real(dp), dimension(8) :: weights
        integer, dimension(3) :: pos

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
        if (uniform_grid_flag) then
            dxm = mhd_config%dx
        else
            ix = binarySearch_R(xpos_local, ptl%x) - 2
            dxm = xpos_local(ix+1) - xpos_local(ix)
        endif

        deltax = 0.0d0
        deltap = 0.0d0
        deltav = 0.0d0
        deltamu = 0.0d0

        !< Calculate the 1st-order and 2nd-order rates
        if (spherical_coord_flag) then
            ir = 1.0d0 / ptl%x
            dv_dt = vx + ptl%v * ptl%mu + kappa%dkxx_dx + 2.0 * kappa%kxx * ir
            cot_theta = 1.0d0 / tan(ptl%y)
            divv = dvx_dx + 2.0 * vx * ir
            bb_gradv = (bx*bx*dvx_dx - bx*(by*vy + bz*vz) * ir + &
                        by*bx*dvy_dx + by*(by*vx - cot_theta*bz*vz) * ir + &
                        bz*bx*dvz_dx + bz**2 * (vx + cot_theta*vy) * ir) / (b*b)
            bv_gradv = (bx*vx*dvx_dx - bx*(vy*vy + vz*vz) * ir + &
                        by*vx*dvy_dx + by*(vy*vx - cot_theta*vz*vz) * ir + &
                        bz*vx*dvz_dx + bz*vz * (vx + cot_theta*vy) * ir) / (b*b)
        else
            dv_dt = vx + ptl%v * ptl%mu + kappa%dkxx_dx
            divv = dvx_dx
            bb_gradv = (bx*bx*dvx_dx + by*bx*dvy_dx + bz*bx*dvz_dx) / (b*b)
            bv_gradv = (bx*vx*dvx_dx + by*vx*dvy_dx + bz*vx*dvz_dx) / b
        endif
        mu2 = ptl%mu * ptl%mu
        acc_rate = -(0.5 *(1 - mu2) * divv + &
                     0.5 *(3 * mu2 - 1) * bb_gradv + &
                     ptl%mu * bv_gradv / ptl%v)
        dp_dt = ptl%p * acc_rate

        ! Momentum diffusion due to wave scattering
        dpp = 0.0d0
        if (dpp_wave_flag) then
            rho = fields(4)
            call calc_dpp_wave_scattering(rho, b, kappa%kpara, ptl, dp_dt, dpp)
        endif

        ! Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            sigmaxx = dvx_dx - divv / 3
            sigmayy = -divv / 3
            sigmazz = -divv / 3
            call calc_dpp_flow_shear(b, bx, by, bz, kappa%knorm_para, &
                sigmaxx, sigmayy, sigmazz, 0.0d0, 0.0d0, 0.0d0, ptl, dp_dt, dpp)
        endif

        ! Pitch-angle evolution
        db_dx = fields(nfields+22)
        div_bnorm = -bx * db_dx / (b*b)
        call calc_duu(ptl, b, db2_slab(1), lc_slab(1), div_bnorm, &
            divv, bb_gradv, bv_gradv, mu2, dmu_dt, duu, duu_du)

        !< Set the time step
        if (.not. fixed_dt) then
            if ((dx_dt .ne. 0.0d0) .and. &
                (dp_dt .ne. 0.0d0) .and. &
                (dmu_dt .ne. 0.0d0)) then
                if (kappa%skperp > 0.0_dp) then
                    ptl%dt = min((0.5*dxm/kappa%skperp)**2, &
                                 (kappa%skperp/dx_dt)**2, &
                                 0.1*ptl%p/abs(dp_dt), &
                                 0.1/abs(dmu_dt), 2.0*duu/dmu_dt**2)
                else
                    ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                 (kappa%skpara/dx_dt)**2, &
                                 0.1*ptl%p/abs(dp_dt), &
                                 0.1/abs(dmu_dt), 2.0*duu/dmu_dt**2)
                endif
            else
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too small
            if (ptl%dt .lt. dt_min) then
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too large
            if (ptl%dt .gt. dt_max) then
                ptl%dt = dt_max
            endif
        endif

        !< Update the particle
        sdt = dsqrt(ptl%dt)
        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltax = dx_dt * ptl%dt + ran1*kappa%skpara*sdt
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltap = dp_dt * ptl%dt + ran1*dsqrt(2*dpp)*sdt
        deltav = ptl%v * deltap / ptl%p
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltamu = dmu_dt * ptl%dt + ran1*dsqrt(2*duu)*sdt

        ptl%x = ptl%x + deltax
        ptl%mu = ptl%mu + deltamu
        ptl%t = ptl%t + ptl%dt

        if (ptl%mu > mu_max) then
            deltamu = mu_max - (ptl%mu - deltamu)
            ptl%mu = mu_max
        else if (ptl%mu < -mu_max) then
            deltamu = -mu_max - (ptl%mu - deltamu)
            ptl%mu = -mu_max
        endif

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

        if (ptl%p < 0.25 * p0) then
            ptl%v = ptl%v - deltav
            deltav = ptl%v * 0.25 * p0 / ptl%p - ptl%v
            ptl%v = ptl%v + deltav
            ptl%p = ptl%p - deltap
            deltap = 0.25 * p0 - ptl%p
            ptl%p = 0.25 * p0
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
    !<  kappa: kappa and related variables
    !<  fixed_dt: whether to fix the time step
    !<  deltax, deltay, deltap: the change of x, y and p in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_2d(thread_id, rt, ptl, fields, &
            kappa, fixed_dt, deltax, deltay, deltap)
        use constants, only: pi
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        logical, intent(in) :: fixed_dt
        real(dp), intent(in) :: rt
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltay, deltap
        real(dp) :: sdt, dvx_dx, dvx_dy, dvy_dx, dvy_dy, divv
        real(dp) :: bx, by, bz, b, vx, vy, px, py, pz, rt1, ib
        real(dp) :: bx1, by1, btot1, ib1
        real(dp) :: dbx_dy, dby_dx, db_dx, db_dy, dbz_dx, dbz_dy
        real(dp) :: vdx, vdy, vdz, vdp, ib2, ib3, gbr, gbt
        real(dp) :: dxm, dym
        real(dp) :: skperp, skpara_perp, skperp1, skpara_perp1
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rho ! Plasma density
        real(dp) :: rands(2)
        real(dp) :: a1, b1, c1, Qpp, Qpm, Qmp, Qmm
        real(dp) :: ctheta, istheta, cot_theta, ir, ir2
        real(dp) :: qtmp1, qtmp2, atmp
        real(dp) :: deltaz
        real(dp) :: sigmaxx, sigmayy, sigmazz, sigmaxy ! shear tensor
        real(dp) :: dx_dt, dy_dt, dz_dt, dp_dt, dpp
        integer :: ix, iy
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
        if (uniform_grid_flag) then
            dxm = mhd_config%dx
            dym = mhd_config%dy
        else
            ix = binarySearch_R(xpos_local, ptl%x) - 2
            iy = binarySearch_R(ypos_local, ptl%y) - 2
            dxm = xpos_local(ix+1) - xpos_local(ix)
            dym = ypos_local(iy+1) - ypos_local(iy)
        endif

        deltax = 0.0d0
        deltay = 0.0d0
        deltap = 0.0d0

        if (b < EPSILON(b)) then
            ib = 0.0
        else
            ib = 1.0 / b
        endif

        !< Parameters used below
        if (spherical_coord_flag) then
            ctheta = cos(ptl%y)
            istheta = 1.0 / sin(ptl%y)
            cot_theta = ctheta * istheta
            ir = 1.0 / ptl%x
            ir2 = 1.0 / ptl%x**2
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
        else
            vdx = vdp * (dbz_dy * ib2 - 2.0 * bz * db_dy * ib3)
            vdy = vdp * (-dbz_dx * ib2 + 2.0 * bz * db_dx * ib3)
        endif

        !< Check particle drift along the out-of-plane direction
        if (check_drift_2d) then
            dbx_dy = fields(nfields+14)
            dby_dx = fields(nfields+16)
            if (spherical_coord_flag) then
                vdz = vdp * ((dby_dx + by*ir - dbx_dy*ir)*ib2 - &
                             (gbr*by - gbt*bx*ir)*2.0*ib3)
            else
                vdz = vdp * ((dby_dx-dbx_dy)*ib2 - 2*(by*db_dx-bx*db_dy)*ib3)
            endif
        else
            vdz = 0.0d0
        endif

        !< Calculate the 1st-order and 2nd-order rates
        if (spherical_coord_flag) then
            dx_dt = vx + vdx + kappa%dkxx_dx + &
                (2.0*kappa%kxx + kappa%dkxy_dy + kappa%kxy*cot_theta)*ir
            dy_dt = (vy + vdy + kappa%dkxy_dx)*ir + &
                (kappa%kxy + kappa%dkyy_dy + kappa%kyy*cot_theta)*ir2
            dz_dt = vdz * istheta * ir
            divv = dvx_dx + (2.0*vx + dvy_dy + vy*ctheta*istheta)*ir
        else
            dx_dt = vx + vdx + kappa%dkxx_dx + kappa%dkxy_dy
            dy_dt = vy + vdy + kappa%dkxy_dx + kappa%dkyy_dy
            dz_dt = vdz
            divv = dvx_dx + dvy_dy
        endif

        dp_dt = -ptl%p * divv / 3.0d0

        !< Momentum diffusion due to wave scattering
        dpp = 0.0d0
        if (dpp_wave_flag) then
            rho = fields(4)
            call calc_dpp_wave_scattering(rho, b, kappa%kpara, ptl, dp_dt, dpp)
        endif

        !< Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            dvx_dy = fields(nfields+2)
            dvy_dx = fields(nfields+4)
            sigmaxx = dvx_dx - divv / 3
            sigmayy = dvy_dy - divv / 3
            sigmazz = -divv / 3
            sigmaxy = (dvx_dy + dvy_dx) / 2
            call calc_dpp_flow_shear(b, bx, by, bz, kappa%knorm_para, &
                sigmaxx, sigmayy, sigmazz, sigmaxy, 0.0d0, 0.0d0, ptl, dp_dt, dpp)
        endif

        !< Set the time step
        if (.not. fixed_dt) then
            if ((dx_dt .ne. 0.0d0) .and. &
                (dy_dt .ne. 0.0d0) .and. &
                (dp_dt .ne. 0.0d0)) then
                if (spherical_coord_flag) then
                    if (kappa%skperp > 0.0_dp) then
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym*ptl%x/kappa%skpara)**2, &
                                     (kappa%skperp/dx_dt)**2, &
                                     (kappa%skperp*ir/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt))
                    else
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym*ptl%x/kappa%skpara)**2, &
                                     (kappa%skpara/dx_dt)**2, &
                                     (kappa%skpara*ir/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt))
                    endif
                else
                    if (kappa%skperp > 0.0_dp) then
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym/kappa%skpara)**2, &
                                     (kappa%skperp/dx_dt)**2, &
                                     (kappa%skperp/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt))
                    else
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym/kappa%skpara)**2, &
                                     (kappa%skpara/dx_dt)**2, &
                                     (kappa%skpara/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt))
                    endif
                endif
            else
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too small
            if (ptl%dt .lt. dt_min) then
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too large
            if (ptl%dt .gt. dt_max) then
                ptl%dt = dt_max
            endif
        endif

        !< Update the particle
        sdt = dsqrt(ptl%dt)
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
            a1 = kappa%kxx
            b1 = kappa%kxy * ir
            c1 = kappa%kyy * ir2
            atmp = dsqrt((a1-c1)**2 + 4*b1**2)
            Qpp = atmp + (a1 + c1)
            Qmp = atmp - (a1 + c1)
            Qpm = atmp + (a1 - c1)
            Qmm = atmp - (a1 - c1)
            if (Qmp > 0.0_dp) then
                ! atmp and (a1 + c1) can be very close
                Qmp = 0.0_dp
            endif
            qtmp1 = dsqrt(-Qmp/(Qmm**2+4*b1**2))
            qtmp2 = dsqrt(Qpp/(Qpm**2+4*b1**2))
            deltax = dx_dt * ptl%dt + (-Qmm*qtmp1*ran1 + Qpm*qtmp2*ran2)*sdt
            deltay = dy_dt * ptl%dt + (2*b1*qtmp1*ran1 + 2*b1*qtmp2*ran2)*sdt
        else
            deltax = dx_dt * ptl%dt + ran1*kappa%skperp*sdt + &
                ran3*kappa%skpara_perp*sdt*bx*ib
            deltay = dy_dt * ptl%dt + ran2*kappa%skperp*sdt + &
                ran3*kappa%skpara_perp*sdt*by*ib
        endif
        deltaz = dz_dt * ptl%dt

        ptl%x = ptl%x + deltax
        ptl%y = ptl%y + deltay
        ptl%z = ptl%z + deltaz
        ptl%t = ptl%t + ptl%dt

        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltap = dp_dt * ptl%dt + ran1*dsqrt(2*dpp)*sdt

        if (acc_region_flag == 1) then
            if (particle_in_acceleration_region(ptl)) then
                ptl%p = ptl%p + deltap
            else
                deltap = 0.0
            endif
        else
            ptl%p = ptl%p + deltap
        endif
        if (ptl%p < 0.25 * p0) then
            ptl%p = ptl%p - deltap
            deltap = 0.25 * p0 - ptl%p
            ptl%p = 0.25 * p0
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
    !<  db2_slab: turbulence variance for the slab component
    !<  lc_slab: turbulence correlation length for the slab component
    !<  kappa: kappa and related variables
    !<  fixed_dt: whether to fix the time step
    !<  deltax: the change of x in this step
    !<  deltay: the change of y in this step
    !<  deltap: the change of p in this step
    !<  deltav: the change of v in this step
    !<  deltamu: the change of mu in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_2d_ft(thread_id, rt, ptl, fields, db2_slab, lc_slab, &
            kappa, fixed_dt, deltax, deltay, deltap, deltav, deltamu)
        use constants, only: pi
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        logical, intent(in) :: fixed_dt
        real(dp), intent(in) :: rt
        real(dp), dimension(*), intent(in) :: db2_slab, lc_slab
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltay, deltap
        real(dp), intent(out) :: deltav, deltamu
        real(dp) :: xtmp, ytmp
        real(dp) :: sdt, divv, bb_gradv, bv_gradv
        real(dp) :: dvx_dx, dvx_dy, dvy_dx, dvy_dy, dvz_dx, dvz_dy
        real(dp) :: bx, by, bz, b, ib, rt1, bxn, byn, bzn, ibxyn
        real(dp) :: vx, vy, vz, vbx, vby
        real(dp) :: dbx_dx, dbx_dy, dby_dx, dby_dy, dbz_dx, dbz_dy, db_dx, db_dy
        real(dp) :: vdx, vdy, vdz, vdp, ib2, ib3, gbr, gbt
        real(dp) :: dxm, dym
        real(dp) :: mu2, muf1, muf2, kx, ky, kz, bdot_curvb, acc_rate
        real(dp) :: dp_dpp, div_bnorm, h0
        real(dp) :: duu, dmu_dt, duu_du, dtmp
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rho
        real(dp) :: rands(2)
        real(dp) :: a1, b1, c1, Qpp, Qpm, Qmp, Qmm
        real(dp) :: ctheta, istheta, cot_theta, ir, ir2
        real(dp) :: qtmp1, qtmp2, atmp
        real(dp) :: deltaz
        real(dp) :: sigmaxx, sigmayy, sigmazz, sigmaxy ! shear tensor
        real(dp) :: dx_dt, dy_dt, dz_dt, dp_dt, dv_dt, dpp
        integer :: ix, iy
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights

        vx = fields(1)
        vy = fields(2)
        vz = fields(3)
        bx = fields(5)
        by = fields(6)
        bz = fields(7)
        b = dsqrt(bx**2 + by**2 + bz**2)
        if (uniform_grid_flag) then
            dxm = mhd_config%dx
            dym = mhd_config%dy
        else
            ix = binarySearch_R(xpos_local, ptl%x) - 2
            iy = binarySearch_R(ypos_local, ptl%y) - 2
            dxm = xpos_local(ix+1) - xpos_local(ix)
            dym = ypos_local(iy+1) - ypos_local(iy)
        endif

        deltax = 0.0d0
        deltay = 0.0d0
        deltap = 0.0d0
        deltav = 0.0d0
        deltamu = 0.0d0

        if (b < EPSILON(b)) then
            ib = 0.0
        else
            ib = 1.0 / b
        endif

        !< Parameters used below
        if (spherical_coord_flag) then
            ctheta = cos(ptl%y)
            istheta = 1.0 / sin(ptl%y)
            cot_theta = ctheta * istheta
            ir = 1.0 / ptl%x
            ir2 = 1.0 / ptl%x**2
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
        endif

        !< Check particle drift along the out-of-plane direction
        if (check_drift_2d) then
            if (spherical_coord_flag) then
                vdz = vdp * (muf1 * (bx * gbt - by * gbr) * ib2 + &
                             mu2 * (bx * ky - by * kx) * ib3 + &
                             muf1 * bz * bdot_curvb * ib3)
            else
                vdz = vdp * (muf1 * (bx * db_dy - by * db_dx) * ib2 + &
                             mu2 * (bx * ky - by * kx) * ib3 + &
                             muf1 * bz * bdot_curvb * ib3)
            endif
        else
            vdz = 0.0d0
        endif

        !< Calculate the 1st-order and 2nd-order rates
        vbx = ptl%v * ptl%mu * ib
        vby = vbx * by ! particle velocity along the magnetic field
        vbx = vbx * bx
        dvx_dx = fields(nfields+1)
        dvx_dy = fields(nfields+2)
        dvy_dx = fields(nfields+4)
        dvy_dy = fields(nfields+5)
        dvz_dx = fields(nfields+7)
        dvz_dy = fields(nfields+8)
        if (spherical_coord_flag) then
            dx_dt = vx + vdx + vbx + kappa%dkxx_dx + &
                (2.0*kappa%kxx + kappa%dkxy_dy + kappa%kxy*cot_theta)*ir
            dy_dt = (vy + vdy + vby + kappa%dkxy_dx)*ir + &
                (kappa%kxy + kappa%dkyy_dy + kappa%kyy*cot_theta)*ir2
            dz_dt = vdz * istheta * ir
            divv = dvx_dx + (2.0*vx + dvy_dy + vy*cot_theta)*ir
            bb_gradv = (bx*bx*dvx_dx + bx*by*dvx_dy*ir - bx*(by*vy + bz*vz)*ir + &
                        by*bx*dvy_dx + by*by*dvy_dy*ir + by*(by*vx - cot_theta*bz*vz)*ir + &
                        bz*bx*dvz_dx + bz*by*dvz_dy*ir + bz**2*(vx + cot_theta*vy)*ir) * ib2
            bv_gradv = (bx*vx*dvx_dx + bx*vy*dvx_dy*ir - bx*(vy*vy + vz*vz)*ir + &
                        by*vx*dvy_dx + by*vy*dvy_dy*ir + by*(vy*vx - cot_theta*vz*vz)*ir + &
                        bz*vx*dvz_dx + bz*vy*dvz_dy*ir + bz*vz*(vx + cot_theta*vy)*ir) * ib
        else
            dx_dt = vx + vdx + vbx + kappa%dkxx_dx + kappa%dkxy_dy
            dy_dt = vy + vdy + vby + kappa%dkxy_dx + kappa%dkyy_dy
            dz_dt = vdz
            divv = dvx_dx + dvy_dy
            bb_gradv = (bx * (bx * dvx_dx + by * dvx_dy) + &
                        by * (bx * dvy_dx + by * dvy_dy) + &
                        bz * (bx * dvz_dx + by * dvz_dy)) * ib2
            bv_gradv = (bx * (vx * dvx_dx + vy * dvx_dy) + &
                        by * (vx * dvy_dx + vy * dvy_dy) + &
                        bz * (vx * dvz_dx + vy * dvz_dy)) * ib
        endif
        acc_rate = -(muf1 * divv + muf2 * bb_gradv + ptl%mu * bv_gradv / ptl%v)
        dp_dt = ptl%p * acc_rate

        !< Momentum diffusion due to wave scattering
        dpp = 0.0d0
        if (dpp_wave_flag) then
            rho = fields(4)
            call calc_dpp_wave_scattering(rho, b, kappa%kpara, ptl, dp_dt, dpp)
        endif

        !< Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            dvx_dy = fields(nfields+2)
            dvy_dx = fields(nfields+4)
            sigmaxx = dvx_dx - divv / 3
            sigmayy = dvy_dy - divv / 3
            sigmazz = -divv / 3
            sigmaxy = (dvx_dy + dvy_dx) / 2
            call calc_dpp_flow_shear(b, bx, by, bz, kappa%knorm_para, &
                sigmaxx, sigmayy, sigmazz, sigmaxy, 0.0d0, 0.0d0, ptl, dp_dt, dpp)
        endif

        !< Pitch-angle evolution
        if (spherical_coord_flag) then
            div_bnorm = -(bx * db_dx + by * db_dy * ir) * ib2
        else
            div_bnorm = -(bx * db_dx + by * db_dy) * ib2
        endif
        call calc_duu(ptl, b, db2_slab(1), lc_slab(1), div_bnorm, &
            divv, bb_gradv, bv_gradv, mu2, dmu_dt, duu, duu_du)

        !< Set the time step
        if (.not. fixed_dt) then
            if ((dx_dt .ne. 0.0d0) .and. &
                (dy_dt .ne. 0.0d0) .and. &
                (dp_dt .ne. 0.0d0) .and. &
                (dmu_dt .ne. 0.0d0)) then
                if (spherical_coord_flag) then
                    if (kappa%skperp > 0.0_dp) then
                        ptl%dt = min((0.5*dxm/kappa%skperp)**2, &
                                     (0.5*dym*ptl%x/kappa%skperp)**2, &
                                     (kappa%skperp/dx_dt)**2, &
                                     (kappa%skperp*ir/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt), &
                                     0.1/abs(dmu_dt), &
                                     2.0*duu/dmu_dt**2)
                    else
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym*ptl%x/kappa%skpara)**2, &
                                     (kappa%skpara/dx_dt)**2, &
                                     (kappa%skpara*ir/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt), &
                                     0.1/abs(dmu_dt), &
                                     2.0*duu/dmu_dt**2)
                    endif
                else
                    if (kappa%skperp > 0.0_dp) then
                        ptl%dt = min((0.5*dxm/kappa%skperp)**2, &
                                     (0.5*dym/kappa%skperp)**2, &
                                     (kappa%skperp/dx_dt)**2, &
                                     (kappa%skperp/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt), &
                                     0.1/abs(dmu_dt), &
                                     2.0*duu/dmu_dt**2)
                    else
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym/kappa%skpara)**2, &
                                     (kappa%skpara/dx_dt)**2, &
                                     (kappa%skpara/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt), &
                                     0.1/abs(dmu_dt), &
                                     2.0*duu/dmu_dt**2)
                    endif
                endif
            else
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too small
            if (ptl%dt .lt. dt_min) then
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too large
            if (ptl%dt .gt. dt_max) then
                ptl%dt = dt_max
            endif
        endif

        !< Update the particle
        sdt = dsqrt(ptl%dt)
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
            a1 = kappa%kxx
            b1 = kappa%kxy * ir
            c1 = kappa%kyy * ir2
            atmp = dsqrt((a1-c1)**2 + 4*b1**2)
            Qpp = atmp + (a1 + c1)
            Qmp = atmp - (a1 + c1)
            Qpm = atmp + (a1 - c1)
            Qmm = atmp - (a1 - c1)
            if (Qmp > 0.0_dp) then
                ! atmp and (a1 + c1) can be very close
                Qmp = 0.0_dp
            endif
            qtmp1 = dsqrt(-Qmp/(Qmm**2+4*b1**2))
            qtmp2 = dsqrt(Qpp/(Qpm**2+4*b1**2))
            deltax = dx_dt * ptl%dt + (-Qmm*qtmp1*ran1 + Qpm*qtmp2*ran2)*sdt
            deltay = dy_dt * ptl%dt + (2*b1*qtmp1*ran1 + 2*b1*qtmp2*ran2)*sdt
        else
            bxn = bx * ib
            byn = by * ib
            bzn = bz * ib
            ibxyn = 1.0d0 / dsqrt(bxn**2 + byn**2)
            deltax = dx_dt * ptl%dt + kappa%skperp * ibxyn * sdt * (-bxn*bzn*ran1 - by*ran2)
            deltay = dy_dt * ptl%dt + kappa%skperp * ibxyn * sdt * (-byn*bzn*ran1 + bx*ran2)
        endif
        deltaz = dz_dt * ptl%dt
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltap = dp_dt * ptl%dt + ran1*dsqrt(2*dpp)*sdt
        deltav = ptl%v * deltap / ptl%p
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltamu = dmu_dt * ptl%dt + ran1*dsqrt(2*duu)*sdt

        ptl%x = ptl%x + deltax
        ptl%y = ptl%y + deltay
        ptl%z = ptl%z + deltaz
        ptl%mu = ptl%mu + deltamu
        ptl%t = ptl%t + ptl%dt

        if (ptl%mu > mu_max) then
            deltamu = mu_max - (ptl%mu - deltamu)
            ptl%mu = mu_max
        else if (ptl%mu < -mu_max) then
            deltamu = -mu_max - (ptl%mu - deltamu)
            ptl%mu = -mu_max
        endif

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

        if (ptl%p < 0.25 * p0) then
            ptl%v = ptl%v - deltav
            deltav = ptl%v * 0.25 * p0 / ptl%p - ptl%v
            ptl%v = ptl%v + deltav
            ptl%p = ptl%p - deltap
            deltap = 0.25 * p0 - ptl%p
            ptl%p = 0.25 * p0
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
    !<  kappa: kappa and related variables
    !<  fixed_dt: whether to fix the time step
    !<  deltax, deltay, deltaz, deltap: the change of x, y, z and p in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_2d_include_3rd(thread_id, rt, ptl, fields, &
            kappa, fixed_dt, deltax, deltay, deltaz, deltap)
        use constants, only: pi
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        logical, intent(in) :: fixed_dt
        real(dp), intent(in) :: rt
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltay, deltaz, deltap
        real(dp) :: rt1, sdt
        real(dp) :: dvx_dx, dvx_dy, dvx_dz
        real(dp) :: dvy_dx, dvy_dy, dvy_dz
        real(dp) :: dvz_dx, dvz_dy, dvz_dz
        real(dp) :: divv
        real(dp) :: vx, vy, vz
        real(dp) :: bx, by, bz, b, ib, ib2, ib3
        real(dp) :: bxn, byn, bzn, bxyn, ibxyn
        real(dp) :: px, py, p
        real(dp) :: dbx_dy, dbx_dz, dby_dx, dby_dz, dbz_dx, dbz_dy
        real(dp) :: db_dx, db_dy, db_dz
        real(dp) :: vdx, vdy, vdz, vdp
        real(dp) :: dxm, dym
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rho
        real(dp) :: rands(2)
        real(dp) :: ctheta, istheta, ir, ir2
        real(dp) :: gbr, gbt, gbp, p11, p12, p13, p22, p23, p33
        real(dp) :: sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, sigmayz ! shear tensor
        real(dp) :: dx_dt, dy_dt, dz_dt, dp_dt, dv_dt, dpp
        integer :: ix, iy
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights

        vx = fields(1)
        vy = fields(2)
        vz = fields(3)
        bx = fields(5)
        by = fields(6)
        bz = fields(7)
        b = dsqrt(bx**2 + by**2 + bz**2)
        if (b < EPSILON(b)) then
            ib = 0.0_dp
        else
            ib = 1.0_dp / b
        endif
        bxn = bx * ib
        byn = by * ib
        bzn = bz * ib
        bxyn = dsqrt(bxn**2 + byn**2)
        if (bxyn < EPSILON(bxyn)) then
            ibxyn = 0.0d0
        else
            ibxyn = 1.0_dp / bxyn
        endif
        dvx_dx = fields(nfields+1)
        dvy_dy = fields(nfields+5)
        dvz_dz = 0.0_dp
        if (uniform_grid_flag) then
            dxm = mhd_config%dx
            dym = mhd_config%dy
        else
            ix = binarySearch_R(xpos_local, ptl%x) - 2
            iy = binarySearch_R(ypos_local, ptl%y) - 2
            dxm = xpos_local(ix+1) - xpos_local(ix)
            dym = ypos_local(iy+1) - ypos_local(iy)
        endif

        deltax = 0.0d0
        deltay = 0.0d0
        deltaz = 0.0d0
        deltap = 0.0d0

        !< Parameters used below
        if (spherical_coord_flag) then
            ctheta = cos(ptl%y)
            istheta = 1.0 / sin(ptl%y)
            ir = 1.0 / ptl%x
            ir2 = 1.0 / ptl%x**2
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

        !< Calculate the 1st-order and 2nd-order rates
        kappa%dkxz_dz = 0.0_dp
        kappa%dkyz_dz = 0.0_dp
        kappa%dkzz_dz = 0.0_dp
        if (spherical_coord_flag) then
            dx_dt = vx + vdx + kappa%dkxx_dx + &
                (2.0*kappa%kxx + kappa%dkxy_dy + &
                kappa%kxy*ctheta*istheta + kappa%dkxz_dz*istheta)*ir
            dy_dt = (vy + vdy + kappa%dkxy_dx)*ir + &
                (kappa%kxy + kappa%dkyy_dy + kappa%kyy*ctheta*istheta + &
                kappa%dkyz_dz*istheta)*ir2
            dz_dt = (vz + vdz + kappa%dkxz_dx)*istheta*ir + &
                (kappa%kxz + kappa%dkyz_dy + kappa%dkzz_dz*istheta)*istheta*ir2
            divv = dvx_dx + (2.0*vx + dvy_dy + vy*ctheta*istheta + dvz_dz*istheta)*ir
        else
            dx_dt = vx + vdx + kappa%dkxx_dx + kappa%dkxy_dy + kappa%dkxz_dz
            dy_dt = vy + vdy + kappa%dkxy_dx + kappa%dkyy_dy + kappa%dkyz_dz
            dz_dt = vz + vdz + kappa%dkxz_dx + kappa%dkyz_dy + kappa%dkzz_dz
            divv = dvx_dx + dvy_dy + dvz_dz
        endif
        dp_dt = -ptl%p * divv / 3.0d0

        !< Momentum diffusion due to wave scattering
        dpp = 0.0d0
        if (dpp_wave_flag) then
            rho = fields(4)
            call calc_dpp_wave_scattering(rho, b, kappa%kpara, ptl, dp_dt, dpp)
        endif

        !< Momentum diffusion due to flow shear
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
            call calc_dpp_flow_shear(b, bx, by, bz, kappa%knorm_para, &
                sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, sigmayz, &
                ptl, dp_dt, dpp)
        endif

        !< Set the time step
        if (.not. fixed_dt) then
            if ((dx_dt .ne. 0.0d0) .and. &
                (dy_dt .ne. 0.0d0) .and. &
                (dp_dt .ne. 0.0d0)) then
                if (spherical_coord_flag) then
                    if (kappa%skperp > 0.0_dp) then
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym*ptl%x/kappa%skpara)**2, &
                                     (kappa%skperp/dx_dt)**2, &
                                     (kappa%skperp*ir/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt))
                    else
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym*ptl%x/kappa%skpara)**2, &
                                     (kappa%skpara/dx_dt)**2, &
                                     (kappa%skpara*ir/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt))
                    endif
                else
                    if (kappa%skperp > 0.0_dp) then
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym/kappa%skpara)**2, &
                                     (kappa%skperp/dx_dt)**2, &
                                     (kappa%skperp/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt))
                    else
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym/kappa%skpara)**2, &
                                     (kappa%skpara/dx_dt)**2, &
                                     (kappa%skpara/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt))
                    endif
                endif
            else
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too small
            if (ptl%dt .lt. dt_min) then
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too large
            if (ptl%dt .gt. dt_max) then
                ptl%dt = dt_max
            endif
        endif

        !< Update the particle
        sdt = dsqrt(ptl%dt)
        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran2 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran3 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3

        if (spherical_coord_flag) then
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
            deltax = dx_dt * ptl%dt + (p11*ran1 + p12*ran2 + p13*ran3)*sdt
            deltay = dy_dt * ptl%dt + (p22*ran2 + p23*ran3)*sdt
            deltaz = dz_dt * ptl%dt + p33*ran3*sdt
        else
            deltax = dx_dt * ptl%dt + &
                (bxn*kappa%skpara*ran1 - bxn*bzn*kappa%skperp*ibxyn*ran2 - &
                 byn*kappa%skperp*ibxyn*ran3)*sdt
            deltay = dy_dt * ptl%dt + &
                (byn*kappa%skpara*ran1 - byn*bzn*kappa%skperp*ibxyn*ran2 + &
                 bxn*kappa%skperp*ibxyn*ran3)*sdt
            deltaz = dz_dt * ptl%dt + &
                (bzn*kappa%skpara*ran1 + bxyn*kappa%skperp*ran2)*sdt
        endif

        ptl%x = ptl%x + deltax
        ptl%y = ptl%y + deltay
        ptl%z = ptl%z + deltaz
        ptl%t = ptl%t + ptl%dt

        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltap = dp_dt * ptl%dt + ran1*dsqrt(2*dpp)*sdt

        if (acc_region_flag == 1) then
            if (particle_in_acceleration_region(ptl)) then
                ptl%p = ptl%p + deltap
            else
                deltap = 0.0
            endif
        else
            ptl%p = ptl%p + deltap
        endif
        if (ptl%p < 0.25 * p0) then
            ptl%p = ptl%p - deltap
            deltap = 0.25 * p0 - ptl%p
            ptl%p = 0.25 * p0
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
    !<  db2_slab: turbulence variance for the slab component
    !<  lc_slab: turbulence correlation length for the slab component
    !<  kappa: kappa and related variables
    !<  fixed_dt: whether to fix the time step
    !<  deltax: the change of x in this step
    !<  deltay: the change of y in this step
    !<  deltaz: the change of z in this step
    !<  deltap: the change of p in this step
    !<  deltav: the change of v in this step
    !<  deltamu: the change of mu in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_2d_include_3rd_ft(thread_id, rt, ptl, fields, &
            db2_slab, lc_slab, kappa, fixed_dt, deltax, deltay, deltaz, &
            deltap, deltav, deltamu)
        use constants, only: pi
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        logical, intent(in) :: fixed_dt
        real(dp), intent(in) :: rt
        real(dp), dimension(*), intent(in) :: db2_slab, lc_slab
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltay, deltaz, deltap, deltav, deltamu
        real(dp) :: rt1, sdt
        real(dp) :: dvx_dx, dvx_dy
        real(dp) :: dvy_dx, dvy_dy
        real(dp) :: dvz_dx, dvz_dy
        real(dp) :: divv, bb_gradv, bv_gradv
        real(dp) :: vx, vy, vz, vbx, vby, vbz
        real(dp) :: bx, by, bz, b, ib, ib2, ib3
        real(dp) :: bxn, byn, bzn, bxyn, ibxyn
        real(dp) :: dbx_dx, dbx_dy
        real(dp) :: dby_dx, dby_dy
        real(dp) :: dbz_dx, dbz_dy
        real(dp) :: db_dx, db_dy
        real(dp) :: vdx, vdy, vdz, vdp
        real(dp) :: dxm, dym
        real(dp) :: mu2, muf1, muf2, kx, ky, kz, bdot_curvb, acc_rate
        real(dp) :: dp_dpp, div_bnorm, h0
        real(dp) :: duu, dmu_dt, duu_du, dtmp
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rho
        real(dp) :: rands(2)
        real(dp) :: ctheta, istheta, cot_theta, ir, ir2
        real(dp) :: gbr, gbt, p11, p12, p13, p22, p23, p33
        real(dp) :: sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, sigmayz ! shear tensor
        real(dp) :: dx_dt, dy_dt, dz_dt, dp_dt, dv_dt, dpp
        integer :: ix, iy
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights

        vx = fields(1)
        vy = fields(2)
        vz = fields(3)
        bx = fields(5)
        by = fields(6)
        bz = fields(7)
        b = dsqrt(bx**2 + by**2 + bz**2)
        if (b < EPSILON(b)) then
            ib = 0.0_dp
        else
            ib = 1.0_dp / b
        endif
        bxn = bx * ib
        byn = by * ib
        bzn = bz * ib
        bxyn = dsqrt(bxn**2 + byn**2)
        if (bxyn < EPSILON(bxyn)) then
            ibxyn = 0.0d0
        else
            ibxyn = 1.0_dp / bxyn
        endif
        if (uniform_grid_flag) then
            dxm = mhd_config%dx
            dym = mhd_config%dy
        else
            ix = binarySearch_R(xpos_local, ptl%x) - 2
            iy = binarySearch_R(ypos_local, ptl%y) - 2
            dxm = xpos_local(ix+1) - xpos_local(ix)
            dym = ypos_local(iy+1) - ypos_local(iy)
        endif

        deltax = 0.0d0
        deltay = 0.0d0
        deltaz = 0.0d0
        deltap = 0.0d0
        deltav = 0.0d0
        deltamu = 0.0d0

        !< Parameters used below
        if (spherical_coord_flag) then
            ctheta = cos(ptl%y)
            istheta = 1.0 / sin(ptl%y)
            ir = 1.0 / ptl%x
            ir2 = 1.0 / ptl%x**2
        endif

        !< Drift velocity
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

        !< Calculate the 1st-order and 2nd-order rates
        vbx = ptl%v * ptl%mu * ib
        vby = vbx * by ! particle velocity along the magnetic field
        vbz = vbx * bz
        vbx = vbx * bx
        dvx_dx = fields(nfields+1)
        dvx_dy = fields(nfields+2)
        dvy_dx = fields(nfields+4)
        dvy_dy = fields(nfields+5)
        dvz_dx = fields(nfields+7)
        dvz_dy = fields(nfields+8)
        if (spherical_coord_flag) then
            dx_dt = vx + vbx + vdx + kappa%dkxx_dx + &
                (2.0*kappa%kxx + kappa%dkxy_dy + kappa%kxy*cot_theta)*ir
            dy_dt = (vy + vby + vdy + kappa%dkxy_dx)*ir + &
                (kappa%kxy + kappa%dkyy_dy + kappa%kyy*cot_theta)*ir2
            dz_dt = (vz + vbz + vdz + kappa%dkxz_dx)*istheta*ir + &
                (kappa%kxz + kappa%dkyz_dy)*istheta*ir2
            divv = dvx_dx + (2.0*vx + dvy_dy + vy*cot_theta)*ir
            bb_gradv = (bx*bx*dvx_dx + bx*by*dvx_dy*ir - bx*(by*vy + bz*vz)*ir + &
                        by*bx*dvy_dx + by*by*dvy_dy*ir + by*(by*vx - cot_theta*bz*vz)*ir + &
                        bz*bx*dvz_dx + bz*by*dvz_dy*ir + bz**2*(vx + cot_theta*vy)*ir) * ib2
            bv_gradv = (bx*vx*dvx_dx + bx*vy*dvx_dy*ir - bx*(vy*vy + vz*vz)*ir + &
                        by*vx*dvy_dx + by*vy*dvy_dy*ir + by*(vy*vx - cot_theta*vz*vz)*ir + &
                        bz*vx*dvz_dx + bz*vy*dvz_dy*ir + bz*vz*(vx + cot_theta*vy)*ir) * ib
        else
            dx_dt = vx + vbx + vdx + kappa%dkxx_dx + kappa%dkxy_dy
            dy_dt = vy + vby + vdy + kappa%dkxy_dx + kappa%dkyy_dy
            dz_dt = vz + vbz + vdz + kappa%dkxz_dx + kappa%dkyz_dy
            divv = dvx_dx + dvy_dy
            bb_gradv = (bx * (bx * dvx_dx + by * dvx_dy) + &
                        by * (bx * dvy_dx + by * dvy_dy) + &
                        bz * (bx * dvz_dx + by * dvz_dy)) * ib2
            bv_gradv = (bx * (vx * dvx_dx + vy * dvx_dy) + &
                        by * (vx * dvy_dx + vy * dvy_dy) + &
                        bz * (vx * dvz_dx + vy * dvz_dy)) * ib
        endif
        acc_rate = -(muf1 * divv + muf2 * bb_gradv + ptl%mu * bv_gradv / ptl%v)
        dp_dt = ptl%p * acc_rate

        !< Momentum diffusion due to wave scattering
        dpp = 0.0d0
        if (dpp_wave_flag) then
            rho = fields(4)
            call calc_dpp_wave_scattering(rho, b, kappa%kpara, ptl, dp_dt, dpp)
        endif

        !< Momentum diffusion due to flow shear
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
            call calc_dpp_flow_shear(b, bx, by, bz, kappa%knorm_para, &
                sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, sigmayz, &
                ptl, dp_dt, dpp)
        endif

        !< Pitch-angle evolution
        if (spherical_coord_flag) then
            div_bnorm = -(bx * db_dx + by * db_dy * ir) * ib2
        else
            div_bnorm = -(bx * db_dx + by * db_dy) * ib2
        endif
        call calc_duu(ptl, b, db2_slab(1), lc_slab(1), div_bnorm, &
            divv, bb_gradv, bv_gradv, mu2, dmu_dt, duu, duu_du)

        !< Set the time step
        if (.not. fixed_dt) then
            if ((dx_dt .ne. 0.0d0) .and. &
                (dy_dt .ne. 0.0d0) .and. &
                (dp_dt .ne. 0.0d0) .and. &
                (dmu_dt .ne. 0.0d0)) then
                if (spherical_coord_flag) then
                    if (kappa%skperp > 0.0_dp) then
                        ptl%dt = min((0.5*dxm/kappa%skperp)**2, &
                                     (0.5*dym*ptl%x/kappa%skperp)**2, &
                                     (kappa%skperp/dx_dt)**2, &
                                     (kappa%skperp*ir/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt), &
                                     0.1/abs(dmu_dt), &
                                     2.0*duu/dmu_dt**2)
                    else
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym*ptl%x/kappa%skpara)**2, &
                                     (kappa%skpara/dx_dt)**2, &
                                     (kappa%skpara*ir/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt), &
                                     0.1/abs(dmu_dt), &
                                     2.0*duu/dmu_dt**2)
                    endif
                else
                    if (kappa%skperp > 0.0_dp) then
                        ptl%dt = min((0.5*dxm/kappa%skperp)**2, &
                                     (0.5*dym/kappa%skperp)**2, &
                                     (kappa%skperp/dx_dt)**2, &
                                     (kappa%skperp/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt), &
                                     0.1/abs(dmu_dt), &
                                     2.0*duu/dmu_dt**2)
                    else
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym/kappa%skpara)**2, &
                                     (kappa%skpara/dx_dt)**2, &
                                     (kappa%skpara/dy_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt), &
                                     0.1/abs(dmu_dt), &
                                     2.0*duu/dmu_dt**2)
                    endif
                endif
            else
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too small
            if (ptl%dt .lt. dt_min) then
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too large
            if (ptl%dt .gt. dt_max) then
                ptl%dt = dt_max
            endif
        endif

        !< Update the particle
        sdt = dsqrt(ptl%dt)
        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran2 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran3 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3

        if (spherical_coord_flag) then
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
            deltax = dx_dt * ptl%dt + (p11*ran1 + p12*ran2 + p13*ran3)*sdt
            deltay = dy_dt * ptl%dt + (p22*ran2 + p23*ran3)*sdt
            deltaz = dz_dt * ptl%dt + p33*ran3*sdt
        else
            deltax = dx_dt * ptl%dt + (-bxn*bzn*kappa%skperp*ibxyn*ran1 - &
                byn*kappa%skperp*ibxyn*ran2)*sdt
            deltay = dy_dt * ptl%dt + (-byn*bzn*kappa%skperp*ibxyn*ran1 + &
                bxn*kappa%skperp*ibxyn*ran2)*sdt
            deltaz = dz_dt * ptl%dt + bxyn*kappa%skperp*ran1*sdt
        endif
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltap = dp_dt * ptl%dt + ran1*dsqrt(2*dpp)*sdt
        deltav = ptl%v * deltap / ptl%p
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltamu = dmu_dt * ptl%dt + ran1*dsqrt(2*duu)*sdt

        ptl%x = ptl%x + deltax
        ptl%y = ptl%y + deltay
        ptl%z = ptl%z + deltaz
        ptl%mu = ptl%mu + deltamu
        ptl%t = ptl%t + ptl%dt

        if (ptl%mu > mu_max) then
            deltamu = mu_max - (ptl%mu - deltamu)
            ptl%mu = mu_max
        else if (ptl%mu < -mu_max) then
            deltamu = -mu_max - (ptl%mu - deltamu)
            ptl%mu = -mu_max
        endif

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

        if (ptl%p < 0.25 * p0) then
            ptl%v = ptl%v - deltav
            deltav = ptl%v * 0.25 * p0 / ptl%p - ptl%v
            ptl%v = ptl%v + deltav
            ptl%p = ptl%p - deltap
            deltap = 0.25 * p0 - ptl%p
            ptl%p = 0.25 * p0
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
    !<  kappa: kappa and related variables
    !<  fixed_dt: whether to fix the time step
    !<  deltax, deltay, deltaz, deltap: the change of x, y, z and p in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_3d(thread_id, rt, surface_height1, surface_height2, &
            ptl, fields, kappa, fixed_dt, deltax, deltay, deltaz, deltap)
        use constants, only: pi
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use acc_region_surface, only: check_above_acc_surface
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        logical, intent(in) :: fixed_dt
        real(dp), intent(in) :: rt, surface_height1, surface_height2
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltay, deltaz, deltap
        real(dp) :: rt1, sdt
        real(dp) :: dvx_dx, dvx_dy, dvx_dz
        real(dp) :: dvy_dx, dvy_dy, dvy_dz
        real(dp) :: dvz_dx, dvz_dy, dvz_dz
        real(dp) :: divv
        real(dp) :: vx, vy, vz
        real(dp) :: bx, by, bz, b, ib, ib2, ib3
        real(dp) :: bxn, byn, bzn, bxyn, ibxyn
        real(dp) :: px, py, pz
        real(dp) :: dbx_dy, dbx_dz, dby_dx, dby_dz, dbz_dx, dbz_dy
        real(dp) :: db_dx, db_dy, db_dz
        real(dp) :: vdx, vdy, vdz, vdp
        real(dp) :: dxm, dym, dzm
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rho
        real(dp) :: rands(2)
        real(dp) :: ctheta, istheta, ir, ir2
        real(dp) :: gbr, gbt, gbp, p11, p12, p13, p22, p23, p33
        real(dp) :: sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, sigmayz ! shear tensor
        real(dp) :: dx_dt, dy_dt, dz_dt, dp_dt, dpp
        integer :: ix, iy, iz
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
        if (b < EPSILON(b)) then
            ib = 0.0_dp
        else
            ib = 1.0_dp / b
        endif
        bxn = bx * ib
        byn = by * ib
        bzn = bz * ib
        bxyn = dsqrt(bxn**2 + byn**2)
        if (bxyn < EPSILON(bxyn)) then
            ibxyn = 0.0d0
        else
            ibxyn = 1.0_dp / bxyn
        endif
        dvx_dx = fields(nfields+1)
        dvy_dy = fields(nfields+5)
        dvz_dz = fields(nfields+9)
        if (uniform_grid_flag) then
            dxm = mhd_config%dx
            dym = mhd_config%dy
            dzm = mhd_config%dz
        else
            ix = binarySearch_R(xpos_local, ptl%x) - 2
            iy = binarySearch_R(ypos_local, ptl%y) - 2
            iz = binarySearch_R(zpos_local, ptl%z) - 2
            dxm = xpos_local(ix+1) - xpos_local(ix)
            dym = ypos_local(iy+1) - ypos_local(iy)
            dzm = zpos_local(iz+1) - zpos_local(iz)
        endif

        deltax = 0.0d0
        deltay = 0.0d0
        deltaz = 0.0d0
        deltap = 0.0d0

        !< Parameters used below
        if (spherical_coord_flag) then
            ctheta = cos(ptl%y)
            istheta = 1.0 / sin(ptl%y)
            ir = 1.0 / ptl%x
            ir2 = 1.0 / ptl%x**2
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

        !< Calculate the 1st-order and 2nd-order rates
        if (spherical_coord_flag) then
            dx_dt = vx + vdx + kappa%dkxx_dx + &
                (2.0*kappa%kxx + kappa%dkxy_dy + &
                kappa%kxy*ctheta*istheta + kappa%dkxz_dz*istheta)*ir
            dy_dt = (vy + vdy + kappa%dkxy_dx)*ir + &
                (kappa%kxy + kappa%dkyy_dy + kappa%kyy*ctheta*istheta + &
                kappa%dkyz_dz*istheta)*ir2
            dz_dt = (vz + vdz + kappa%dkxz_dx)*istheta*ir + &
                (kappa%kxz + kappa%dkyz_dy + kappa%dkzz_dz*istheta)*istheta*ir2
            divv = dvx_dx + &
                (2.0*vx + dvy_dy + vy*ctheta*istheta + dvz_dz*istheta)*ir
        else
            dx_dt = vx + vdx + kappa%dkxx_dx + kappa%dkxy_dy + kappa%dkxz_dz
            dy_dt = vy + vdy + kappa%dkxy_dx + kappa%dkyy_dy + kappa%dkyz_dz
            dz_dt = vz + vdz + kappa%dkxz_dx + kappa%dkyz_dy + kappa%dkzz_dz
            divv = dvx_dx + dvy_dy + dvz_dz
        endif
        dp_dt = -ptl%p * divv / 3.0d0

        !< Momentum diffusion due to wave scattering
        dpp = 0.0d0
        if (dpp_wave_flag) then
            rho = fields(4)
            call calc_dpp_wave_scattering(rho, b, kappa%kpara, ptl, dp_dt, dpp)
        endif

        !< Momentum diffusion due to flow shear
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
            call calc_dpp_flow_shear(b, bx, by, bz, kappa%knorm_para, &
                sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, sigmayz, &
                ptl, dp_dt, dpp)
        endif

        !< Set the time step
        if (.not. fixed_dt) then
            if ((dx_dt .ne. 0.0d0) .and. &
                (dy_dt .ne. 0.0d0) .and. &
                (dp_dt .ne. 0.0d0)) then
                if (spherical_coord_flag) then
                    if (kappa%skperp > 0.0_dp) then
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym*ptl%x/kappa%skpara)**2, &
                                     (0.5*dzm*ptl%x*sin(ptl%y)/kappa%skpara)**2, &
                                     (kappa%skperp/dx_dt)**2, &
                                     (kappa%skperp*ir/dy_dt)**2, &
                                     (kappa%skperp*ir*istheta/dz_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt))
                    else
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym*ptl%x/kappa%skpara)**2, &
                                     (0.5*dzm*ptl%x*sin(ptl%y)/kappa%skpara)**2, &
                                     (kappa%skpara/dx_dt)**2, &
                                     (kappa%skpara*ir/dy_dt)**2, &
                                     (kappa%skpara*ir*istheta/dz_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt))
                    endif
                else
                    if (kappa%skperp > 0.0_dp) then
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym/kappa%skpara)**2, &
                                     (0.5*dzm/kappa%skpara)**2, &
                                     (kappa%skperp/dx_dt)**2, &
                                     (kappa%skperp/dy_dt)**2, &
                                     (kappa%skperp/dz_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt))
                    else
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym/kappa%skpara)**2, &
                                     (0.5*dzm/kappa%skpara)**2, &
                                     (kappa%skpara/dx_dt)**2, &
                                     (kappa%skpara/dy_dt)**2, &
                                     (kappa%skpara/dz_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt))
                    endif
                endif
            else
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too small
            if (ptl%dt .lt. dt_min) then
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too large
            if (ptl%dt .gt. dt_max) then
                ptl%dt = dt_max
            endif
        endif

        !< Update the particle
        sdt = dsqrt(ptl%dt)
        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran2 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran3 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3

        if (spherical_coord_flag) then
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
            deltax = dx_dt * ptl%dt + (p11*ran1 + p12*ran2 + p13*ran3)*sdt
            deltay = dy_dt * ptl%dt + (p22*ran2 + p23*ran3)*sdt
            deltaz = dz_dt * ptl%dt + p33*ran3*sdt
        else
            deltax = dx_dt * ptl%dt + &
                (bxn*kappa%skpara*ran1 - bxn*bzn*kappa%skperp*ibxyn*ran2 - &
                 byn*kappa%skperp*ibxyn*ran3)*sdt
            deltay = dy_dt * ptl%dt + &
                (byn*kappa%skpara*ran1 - byn*bzn*kappa%skperp*ibxyn*ran2 + &
                 bxn*kappa%skperp*ibxyn*ran3)*sdt
            deltaz = dz_dt * ptl%dt + &
                (bzn*kappa%skpara*ran1 + bxyn*kappa%skperp*ran2)*sdt
        endif

        ptl%x = ptl%x + deltax
        ptl%y = ptl%y + deltay
        ptl%z = ptl%z + deltaz
        ptl%t = ptl%t + ptl%dt

        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltap = dp_dt * ptl%dt + ran1*dsqrt(2*dpp)*sdt

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
        if (ptl%p < 0.25 * p0) then
            ptl%p = ptl%p - deltap
            deltap = 0.25 * p0 - ptl%p
            ptl%p = 0.25 * p0
        endif
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
    !<  fields: fields and their gradients at particle position
    !<  db2_slab: turbulence variance for the slab component
    !<  kappa: kappa and related variables
    !<  fixed_dt: whether to fix the time step
    !<  deltax: the change of x in this step
    !<  deltay: the change of y in this step
    !<  deltaz: the change of z in this step
    !<  deltap: the change of p in this step
    !<  deltav: the change of v in this step
    !<  deltamu: the change of mu in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_3d_ft(thread_id, rt, surface_height1, surface_height2, &
            ptl, fields, db2_slab, lc_slab, kappa, fixed_dt, deltax, deltay, deltaz, &
            deltap, deltav, deltamu)
        use constants, only: pi
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use acc_region_surface, only: check_above_acc_surface
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: thread_id
        logical, intent(in) :: fixed_dt
        real(dp), intent(in) :: rt, surface_height1, surface_height2
        real(dp), dimension(*), intent(in) :: db2_slab, lc_slab
        type(particle_type), intent(inout) :: ptl
        real(dp), dimension(*), intent(inout) :: fields
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltay, deltaz, deltap, deltav, deltamu
        real(dp) :: rt1, sdt
        real(dp) :: dvx_dx, dvx_dy, dvx_dz
        real(dp) :: dvy_dx, dvy_dy, dvy_dz
        real(dp) :: dvz_dx, dvz_dy, dvz_dz
        real(dp) :: divv, bb_gradv, bv_gradv
        real(dp) :: vx, vy, vz, vbx, vby, vbz
        real(dp) :: bx, by, bz, b, ib, ib2, ib3
        real(dp) :: bxn, byn, bzn, bxyn, ibxyn
        real(dp) :: px, py, pz
        real(dp) :: dbx_dx, dbx_dy, dbx_dz
        real(dp) :: dby_dx, dby_dy, dby_dz
        real(dp) :: dbz_dx, dbz_dy, dbz_dz
        real(dp) :: db_dx, db_dy, db_dz
        real(dp) :: vdx, vdy, vdz, vdp
        real(dp) :: dxm, dym, dzm
        real(dp) :: mu2, muf1, muf2, kx, ky, kz, bdot_curvb, acc_rate
        real(dp) :: dp_dpp, div_bnorm, h0
        real(dp) :: duu, dmu_dt, duu_du, dtmp
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rho
        real(dp) :: rands(2)
        real(dp) :: ctheta, istheta, cot_theta, ir, ir2
        real(dp) :: gbr, gbt, gbp, p11, p12, p13, p22, p23, p33
        real(dp) :: sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, sigmayz ! shear tensor
        real(dp) :: dx_dt, dy_dt, dz_dt, dp_dt, dv_dt, dpp
        integer :: ix, iy, iz
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
        if (b < EPSILON(b)) then
            ib = 0.0_dp
        else
            ib = 1.0_dp / b
        endif
        bxn = bx * ib
        byn = by * ib
        bzn = bz * ib
        bxyn = dsqrt(bxn**2 + byn**2)
        if (bxyn < EPSILON(bxyn)) then
            ibxyn = 0.0d0
        else
            ibxyn = 1.0_dp / bxyn
        endif
        if (uniform_grid_flag) then
            dxm = mhd_config%dx
            dym = mhd_config%dy
            dzm = mhd_config%dz
        else
            ix = binarySearch_R(xpos_local, ptl%x) - 2
            iy = binarySearch_R(ypos_local, ptl%y) - 2
            iz = binarySearch_R(zpos_local, ptl%z) - 2
            dxm = xpos_local(ix+1) - xpos_local(ix)
            dym = ypos_local(iy+1) - ypos_local(iy)
            dzm = zpos_local(iz+1) - zpos_local(iz)
        endif

        deltax = 0.0d0
        deltay = 0.0d0
        deltaz = 0.0d0
        deltap = 0.0d0
        deltav = 0.0d0
        deltamu = 0.0d0

        !< Parameters used below
        if (spherical_coord_flag) then
            ctheta = cos(ptl%y)
            istheta = 1.0 / sin(ptl%y)
            cot_theta = ctheta * istheta
            ir = 1.0 / ptl%x
            ir2 = 1.0 / ptl%x**2
        endif

        !< Drift velocity
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

        !< Calculate the 1st-order and 2nd-order rates
        vbx = ptl%v * ptl%mu * ib
        vby = vbx * by ! particle velocity along the magnetic field
        vbz = vbx * bz
        vbx = vbx * bx
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
            dx_dt = vx + vbx + vdx + kappa%dkxx_dx + &
                (2.0*kappa%kxx + kappa%dkxy_dy + &
                kappa%kxy*ctheta*istheta + kappa%dkxz_dz*istheta)*ir
            dy_dt = (vy + vby + vdy + kappa%dkxy_dx)*ir + &
                (kappa%kxy + kappa%dkyy_dy + kappa%kyy*ctheta*istheta + &
                kappa%dkyz_dz*istheta)*ir2
            dz_dt = (vz + vbz + vdz + kappa%dkxz_dx)*istheta*ir + &
                (kappa%kxz + kappa%dkyz_dy + kappa%dkzz_dz*istheta)*istheta*ir2
            divv = dvx_dx + (2.0*vx + dvy_dy + vy*cot_theta + dvz_dz*istheta)*ir
            dtmp = istheta * ir
            bb_gradv = (bx*bx*dvx_dx + bx*by*dvx_dy*ir + bx*bz*dvx_dz*dtmp - bx*(by*vy + bz*vz)*ir + &
                        by*bx*dvy_dx + by*by*dvy_dy*ir + by*bz*dvy_dz*dtmp + by*(by*vx - cot_theta*bz*vz)*ir + &
                        bz*bx*dvz_dx + bz*by*dvz_dy*ir + bz*bz*dvz_dz*dtmp + bz**2*(vx + cot_theta*vy)*ir) * ib2
            bv_gradv = (bx*vx*dvx_dx + bx*vy*dvx_dy*ir + bx*vz*dvx_dz*dtmp - bx*(vy*vy + vz*vz)*ir + &
                        by*vx*dvy_dx + by*vy*dvy_dy*ir + by*vz*dvy_dz*dtmp + by*(vy*vx - cot_theta*vz*vz)*ir + &
                        bz*vx*dvz_dx + bz*vy*dvz_dy*ir + bz*vz*dvz_dz*dtmp + bz*vz*(vx + cot_theta*vy)*ir) * ib
        else
            dx_dt = vx + vbx + vdx + kappa%dkxx_dx + kappa%dkxy_dy + kappa%dkxz_dz
            dy_dt = vy + vby + vdy + kappa%dkxy_dx + kappa%dkyy_dy + kappa%dkyz_dz
            dz_dt = vz + vbz + vdz + kappa%dkxz_dx + kappa%dkyz_dy + kappa%dkzz_dz
            divv = dvx_dx + dvy_dy + dvz_dz
            bb_gradv = (bx * (bx * dvx_dx + by * dvx_dy + bz * dvx_dz) + &
                        by * (bx * dvy_dx + by * dvy_dy + bz * dvy_dz) + &
                        bz * (bx * dvz_dx + by * dvz_dy + bz * dvz_dz)) * ib2
            bv_gradv = (bx * (vx * dvx_dx + vy * dvx_dy + vz * dvx_dz) + &
                        by * (vx * dvy_dx + vy * dvy_dy + vz * dvy_dz) + &
                        bz * (vx * dvz_dx + vy * dvz_dy + vz * dvz_dz)) * ib
        endif
        acc_rate = -(muf1 * divv + muf2 * bb_gradv + ptl%mu * bv_gradv / ptl%v)
        dp_dt = ptl%p * acc_rate

        !< Momentum diffusion due to wave scattering
        dpp = 0.0d0
        if (dpp_wave_flag) then
            rho = fields(4)
            call calc_dpp_wave_scattering(rho, b, kappa%kpara, ptl, dp_dt, dpp)
        endif

        !< Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            sigmaxx = dvx_dx - divv / 3
            sigmayy = dvy_dy - divv / 3
            sigmazz = dvz_dz - divv / 3
            sigmaxy = (dvx_dy + dvy_dx) / 2
            sigmaxz = (dvx_dz + dvz_dx) / 2
            sigmayz = (dvy_dz + dvz_dy) / 2
            call calc_dpp_flow_shear(b, bx, by, bz, kappa%knorm_para, &
                sigmaxx, sigmayy, sigmazz, sigmaxy, sigmaxz, sigmayz, &
                ptl, dp_dt, dpp)
        endif

        !< Pitch-angle evolution
        if (spherical_coord_flag) then
            div_bnorm = -(bx * db_dx + by * db_dy * ir + bz * db_dz * ir * istheta) * ib2
        else
            div_bnorm = -(bx * db_dx + by * db_dy + bz * db_dz) * ib2
        endif
        call calc_duu(ptl, b, db2_slab(1), lc_slab(1), div_bnorm, &
            divv, bb_gradv, bv_gradv, mu2, dmu_dt, duu, duu_du)

        !< Set the time step
        if (.not. fixed_dt) then
            if ((dx_dt .ne. 0.0d0) .and. &
                (dy_dt .ne. 0.0d0) .and. &
                (dz_dt .ne. 0.0d0) .and. &
                (dp_dt .ne. 0.0d0) .and. &
                (dmu_dt .ne. 0.0d0)) then
                if (spherical_coord_flag) then
                    if (kappa%skperp > 0.0_dp) then
                        ptl%dt = min((0.5*dxm/kappa%skperp)**2, &
                                     (0.5*dym*ptl%x/kappa%skperp)**2, &
                                     (0.5*dzm*ptl%x*sin(ptl%y)/kappa%skperp)**2, &
                                     (kappa%skperp/dx_dt)**2, &
                                     (kappa%skperp*ir/dy_dt)**2, &
                                     (kappa%skperp*ir*istheta/dz_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt), &
                                     0.1/abs(dmu_dt), &
                                     2.0*duu/dmu_dt**2)
                    else
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym*ptl%x/kappa%skpara)**2, &
                                     (0.5*dzm*ptl%x*sin(ptl%y)/kappa%skpara)**2, &
                                     (kappa%skpara/dx_dt)**2, &
                                     (kappa%skpara*ir/dy_dt)**2, &
                                     (kappa%skpara*ir*istheta/dz_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt), &
                                     0.1/abs(dmu_dt), &
                                     2.0*duu/dmu_dt**2)
                    endif
                else
                    if (kappa%skperp > 0.0_dp) then
                        ptl%dt = min((0.5*dxm/kappa%skperp)**2, &
                                     (0.5*dym/kappa%skperp)**2, &
                                     (0.5*dzm/kappa%skperp)**2, &
                                     (kappa%skperp/dx_dt)**2, &
                                     (kappa%skperp/dy_dt)**2, &
                                     (kappa%skperp/dz_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt), &
                                     0.1/abs(dmu_dt), &
                                     2.0*duu/dmu_dt**2)
                    else
                        ptl%dt = min((0.5*dxm/kappa%skpara)**2, &
                                     (0.5*dym/kappa%skpara)**2, &
                                     (0.5*dzm/kappa%skpara)**2, &
                                     (kappa%skpara/dx_dt)**2, &
                                     (kappa%skpara/dy_dt)**2, &
                                     (kappa%skpara/dz_dt)**2, &
                                     0.1*ptl%p/abs(dp_dt), &
                                     0.1/abs(dmu_dt), &
                                     2.0*duu/dmu_dt**2)
                    endif
                endif
            else
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too small
            if (ptl%dt .lt. dt_min) then
                ptl%dt = dt_min
            endif
            !< Make sure the time step is not too large
            if (ptl%dt .gt. dt_max) then
                ptl%dt = dt_max
            endif
        endif

        !< Update the particle
        sdt = dsqrt(ptl%dt)
        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran2 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        ran3 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3

        if (spherical_coord_flag) then
            dtmp = 2.0*(kappa%kxx*kappa%kyz**2 + &
                        kappa%kyy*kappa%kxz**2 + &
                        kappa%kzz*kappa%kxy**2 - &
                        2.0*kappa%kxy*kappa%kxz*kappa%kyz - &
                        kappa%kxx*kappa%kyy*kappa%kzz) / &
                        (kappa%kyz**2 - kappa%kyy*kappa%kzz)
            if (dtmp < 0.0_dp) then ! It can be negatively very small
                dtmp = 0.0_dp
            endif
            p11 = dsqrt(dtmp)
            dtmp = 2.0 * (kappa%kyy - (kappa%kyz**2/kappa%kzz))
            if (dtmp < 0.0_dp) then ! It can be negatively very small
                dtmp = 0.0_dp
            endif
            p12 = (kappa%kxz*kappa%kyz - kappa%kxy*kappa%kzz) * &
                dsqrt(dtmp) / (kappa%kyz**2 - kappa%kyy*kappa%kzz)
            dtmp = 2.0 * (kappa%kyy - kappa%kyz**2/kappa%kzz)
            if (dtmp < 0.0_dp) then ! It can be negatively very small
                dtmp = 0.0_dp
            endif
            p22 = dsqrt(dtmp) * ir
            dtmp = kappa%kzz
            if (kappa%kzz <= 0.0_dp) then ! It can be negatively very small
                p13 = 0.0_dp
                p23 = 0.0_dp
                p33 = 0.0_dp
            else
                p13 = dsqrt(2.0 / kappa%kzz) * kappa%kxz
                p23 = dsqrt(2.0 / kappa%kzz) * kappa%kyz * ir
                p33 = dsqrt(2.0 * kappa%kzz) * istheta * ir
            endif
            deltax = dx_dt * ptl%dt + (p11*ran1 + p12*ran2 + p13*ran3)*sdt
            deltay = dy_dt * ptl%dt + (p22*ran2 + p23*ran3)*sdt
            deltaz = dz_dt * ptl%dt + p33*ran3*sdt
        else
            deltax = dx_dt * ptl%dt + (-bxn*bzn*kappa%skperp*ibxyn*ran1 - &
                byn*kappa%skperp*ibxyn*ran2)*sdt
            deltay = dy_dt * ptl%dt + (-byn*bzn*kappa%skperp*ibxyn*ran1 + &
                bxn*kappa%skperp*ibxyn*ran2)*sdt
            deltaz = dz_dt * ptl%dt + bxyn*kappa%skperp*ran1*sdt
        endif
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltap = dp_dt * ptl%dt + ran1*dsqrt(2*dpp)*sdt
        deltav = ptl%v * deltap / ptl%p
        ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
        deltamu = dmu_dt * ptl%dt + ran1*dsqrt(2*duu)*sdt

        ptl%x = ptl%x + deltax
        ptl%y = ptl%y + deltay
        ptl%z = ptl%z + deltaz
        ptl%mu = ptl%mu + deltamu
        ptl%t = ptl%t + ptl%dt

        if (ptl%mu > mu_max) then
            deltamu = mu_max - (ptl%mu - deltamu)
            ptl%mu = mu_max
        else if (ptl%mu < -mu_max) then
            deltamu = -mu_max - (ptl%mu - deltamu)
            ptl%mu = -mu_max
        endif

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

        if (ptl%p < 0.25 * p0) then
            ptl%v = ptl%v - deltav
            deltav = ptl%v * 0.25 * p0 / ptl%p - ptl%v
            ptl%v = ptl%v + deltav
            ptl%p = ptl%p - deltap
            deltap = 0.25 * p0 - ptl%p
            ptl%p = 0.25 * p0
        endif
    end subroutine push_particle_3d_ft

    !---------------------------------------------------------------------------
    !< Resize escaped particle array
    !---------------------------------------------------------------------------
    subroutine resize_escaped_particles
        implicit none
        type(particle_type), allocatable, dimension(:) :: escaped_ptls_tmp
        integer :: nmax
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
        escaped_ptls%origin = 0
        escaped_ptls%nsteps_tracked = 0
        escaped_ptls%nsteps_pushed = 0
        escaped_ptls%tag_injected = 0
        escaped_ptls%tag_splitted = 0
        escaped_ptls(:nptl_escaped) = escaped_ptls_tmp(:nptl_escaped)
        deallocate(escaped_ptls_tmp)
    end subroutine resize_escaped_particles

    !---------------------------------------------------------------------------
    !< Remove particles from simulation if their count_flags are 0.
    !< Args:
    !<  dump_escaped_dist: whether to dump the distributions of the escaped particles
    !---------------------------------------------------------------------------
    subroutine remove_particles(dump_escaped_dist)
        implicit none
        logical, intent(in) :: dump_escaped_dist
        integer :: i, nremoved
        type(particle_type) :: ptl1

        ! Resize the escaped particle data array if necessary
        if (dump_escaped_dist) then
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
                    if (ptls(i)%count_flag < 0) then
                        ! Copy escaped particles
                        nptl_escaped = nptl_escaped + 1
                        if (dump_escaped_dist) then
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
        integer :: i, j, nrecv, iptl
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
    !<  nsteps_interval: save tracked particle points every nsteps_interval
    !---------------------------------------------------------------------------
    subroutine split_particle(split_ratio, pmin_split, nsteps_interval)
        implicit none
        real(dp), intent(in) :: split_ratio, pmin_split
        integer, intent(in) :: nsteps_interval
        integer :: i, nptl
        integer :: iptl_lo, iptl_hi
        real(dp) :: p_threshold
        type(particle_type) :: ptl
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
                if (ptl%tag_splitted < 0) then ! for tracking particles
                    ptls(nptl_current)%tag_splitted = ptl%tag_splitted - 2**(ptl%split_times-1)
                    ! check the splitted particle
                    if (is_particle_selected(ptls(nptl_current), iptl_lo, iptl_hi)) then
                        if (ptl%nsteps_pushed == 0) then
                            ptls(nptl_current)%nsteps_tracked = ptls(nptl_current)%nsteps_tracked + 1
                            particles_tracked(ptls(nptl_current)%nsteps_tracked, iptl_lo:iptl_hi) = &
                                ptls(nptl_current)
                        endif
                    else
                        ptls(nptl_current)%tag_splitted = -ptls(nptl_current)%tag_splitted
                    endif

                    ! check the original particle
                    if (is_particle_selected(ptl, iptl_lo, iptl_hi)) then
                        if (ptl%nsteps_pushed == 0) then
                            ptl%nsteps_tracked = ptl%nsteps_tracked + 1
                            particles_tracked(ptl%nsteps_tracked, iptl_lo:iptl_hi) = ptl
                        endif
                    else
                        ptl%tag_splitted = -ptl%tag_splitted  ! stop tracking
                    endif
                else
                    ptls(nptl_current)%tag_splitted = ptl%tag_splitted + 2**(ptl%split_times-1)
                endif
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

    !---------------------------------------------------------------------------
    !< Save particle module state for restart
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine save_particle_module_state(iframe, file_path)
        use hdf5_io, only: create_file_h5, write_data_h5, close_file_h5
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        integer, allocatable, dimension(:) :: nptls_current, nptls_split, tags_max
        real(dp), allocatable, dimension(:) :: nptl_leak, nptl_leak_negp
        integer(hsize_t), dimension(1) :: dcount, doffset, dset_dims
        integer(hid_t) :: file_id, dset_id, filespace
        character(len=4) :: ctime
        character(len=256) :: fname
        integer :: error

        if (mpi_rank == master) then
            allocate(nptls_current(mpi_size))
            allocate(nptls_split(mpi_size))
            allocate(tags_max(mpi_size))
            allocate(nptl_leak(mpi_size))
            allocate(nptl_leak_negp(mpi_size))
        endif

        write (ctime,'(i4.4)') iframe
        call MPI_GATHER(nptl_current, 1, MPI_INTEGER, nptls_current, &
            1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_GATHER(nptl_split, 1, MPI_INTEGER, nptls_split, &
            1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_GATHER(tag_max, 1, MPI_INTEGER, tags_max, &
            1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_GATHER(leak, 1, MPI_DOUBLE_PRECISION, nptl_leak, &
            1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_GATHER(leak_negp, 1, MPI_DOUBLE_PRECISION, nptl_leak_negp, &
            1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        if (mpi_rank == master) then
            call h5open_f(error)
            fname = trim(file_path)//'particle_module_state_'//ctime//'.h5'
            call create_file_h5(fname, H5F_ACC_TRUNC_F, file_id, .false., MPI_COMM_WORLD)
            dcount(1) = mpi_size
            doffset(1) = 0
            dset_dims(1) = mpi_size
            call h5screate_simple_f(1, dset_dims, filespace, error)
            call h5dcreate_f(file_id, "nptls_current", H5T_NATIVE_INTEGER, filespace, dset_id, error)
            call write_data_h5(dset_id, dcount, doffset, dset_dims, nptls_current, .false., .false.)
            call h5dclose_f(dset_id, error)
            call h5dcreate_f(file_id, "nptls_split", H5T_NATIVE_INTEGER, filespace, dset_id, error)
            call write_data_h5(dset_id, dcount, doffset, dset_dims, nptls_split, .false., .false.)
            call h5dclose_f(dset_id, error)
            call h5dcreate_f(file_id, "tags_max", H5T_NATIVE_INTEGER, filespace, dset_id, error)
            call write_data_h5(dset_id, dcount, doffset, dset_dims, tags_max, .false., .false.)
            call h5dclose_f(dset_id, error)
            call h5dcreate_f(file_id, "nptl_leak", H5T_NATIVE_DOUBLE, filespace, dset_id, error)
            call write_data_h5(dset_id, dcount, doffset, dset_dims, nptl_leak, .false., .false.)
            call h5dclose_f(dset_id, error)
            call h5dcreate_f(file_id, "nptl_leak_negp", H5T_NATIVE_DOUBLE, filespace, dset_id, error)
            call write_data_h5(dset_id, dcount, doffset, dset_dims, nptl_leak_negp, .false., .false.)
            call h5dclose_f(dset_id, error)
            call h5sclose_f(filespace, error)
            call close_file_h5(file_id)
            call h5close_f(error)

            deallocate(nptls_current, nptls_split, tags_max)
            deallocate(nptl_leak, nptl_leak_negp)
        endif
    end subroutine save_particle_module_state

    !---------------------------------------------------------------------------
    !< Read particle module state for restart
    !< Args:
    !<  iframe: time frame index
    !<  file_path: read data files to this path
    !---------------------------------------------------------------------------
    subroutine read_particle_module_state(iframe, file_path)
        use hdf5_io, only: open_file_h5, read_data_h5, close_file_h5
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        integer, allocatable, dimension(:) :: nptls_current, nptls_split, tags_max
        real(dp), allocatable, dimension(:) :: nptl_leak, nptl_leak_negp
        integer(hsize_t), dimension(1) :: dcount, doffset, dset_dims
        integer(hid_t) :: file_id, dset_id, filespace
        character(len=4) :: ctime
        character(len=256) :: fname
        integer :: error

        if (mpi_rank == master) then
            allocate(nptls_current(mpi_size))
            allocate(nptls_split(mpi_size))
            allocate(tags_max(mpi_size))
            allocate(nptl_leak(mpi_size))
            allocate(nptl_leak_negp(mpi_size))
        endif

        write (ctime,'(i4.4)') iframe
        if (mpi_rank == master) then
            call h5open_f(error)
            fname = trim(file_path)//'particle_module_state_'//ctime//'.h5'
            call open_file_h5(fname, H5F_ACC_RDONLY_F, file_id, .false., MPI_COMM_WORLD)
            dcount(1) = mpi_size
            doffset(1) = 0
            dset_dims(1) = mpi_size
            call h5dopen_f(file_id, "nptls_current", dset_id, error)
            call read_data_h5(dset_id, dcount, doffset, dset_dims, nptls_current, .false., .false.)
            call h5dclose_f(dset_id, error)
            call h5dopen_f(file_id, "nptls_split", dset_id, error)
            call read_data_h5(dset_id, dcount, doffset, dset_dims, nptls_split, .false., .false.)
            call h5dclose_f(dset_id, error)
            call h5dopen_f(file_id, "tags_max", dset_id, error)
            call read_data_h5(dset_id, dcount, doffset, dset_dims, tags_max, .false., .false.)
            call h5dclose_f(dset_id, error)
            call h5dopen_f(file_id, "nptl_leak", dset_id, error)
            call read_data_h5(dset_id, dcount, doffset, dset_dims, nptl_leak, .false., .false.)
            call h5dclose_f(dset_id, error)
            call h5dopen_f(file_id, "nptl_leak_negp", dset_id, error)
            call read_data_h5(dset_id, dcount, doffset, dset_dims, nptl_leak_negp, .false., .false.)
            call h5dclose_f(dset_id, error)
            call close_file_h5(file_id)
            call h5close_f(error)
        endif
        call MPI_SCATTER(nptls_current, 1, MPI_INTEGER, nptl_current, &
            1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_SCATTER(nptls_split, 1, MPI_INTEGER, nptl_split, &
            1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_SCATTER(tags_max, 1, MPI_INTEGER, tag_max, &
            1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_SCATTER(nptl_leak, 1, MPI_DOUBLE_PRECISION, leak, &
            1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_SCATTER(nptl_leak_negp, 1, MPI_DOUBLE_PRECISION, leak_negp, &
            1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)

        if (mpi_rank == master) then
            deallocate(nptls_current, nptls_split, tags_max)
            deallocate(nptl_leak, nptl_leak_negp)
        endif
    end subroutine read_particle_module_state

    !---------------------------------------------------------------------------
    !< Read one double element of the particle data
    !---------------------------------------------------------------------------
    subroutine read_double_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: read_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        real(dp), dimension(:), intent(out) :: fdata
        integer(hid_t) :: dset_id
        integer :: error
        call h5dopen_f(file_id, trim(dset_name), dset_id, error)
        call read_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .true., .true.)
        call h5dclose_f(dset_id, error)
    end subroutine read_double_element

    !---------------------------------------------------------------------------
    !< Read one char-length integer element of the particle data
    !---------------------------------------------------------------------------
    subroutine read_integer1_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: read_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        integer(i1), dimension(:), intent(out) :: fdata
        integer(hid_t) :: dset_id
        integer :: error
        call h5dopen_f(file_id, trim(dset_name), dset_id, error)
        call read_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .true., .true.)
        call h5dclose_f(dset_id, error)
    end subroutine read_integer1_element

    !---------------------------------------------------------------------------
    !< Read one default-length integer element of the particle data
    !---------------------------------------------------------------------------
    subroutine read_integer4_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: read_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        integer(i4), dimension(:), intent(out) :: fdata
        integer(hid_t) :: dset_id
        integer :: error
        call h5dopen_f(file_id, trim(dset_name), dset_id, error)
        call read_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .true., .true.)
        call h5dclose_f(dset_id, error)
    end subroutine read_integer4_element

    !---------------------------------------------------------------------------
    !< Read one long-length integer element of the particle data
    !---------------------------------------------------------------------------
    subroutine read_integer8_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: read_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        integer(i8), dimension(:), intent(out) :: fdata
        integer(hid_t) :: dset_id
        integer :: error
        call h5dopen_f(file_id, trim(dset_name), dset_id, error)
        call read_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .true., .true.)
        call h5dclose_f(dset_id, error)
    end subroutine read_integer8_element

    !---------------------------------------------------------------------------
    !< Read all the particles from file
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine read_particles(iframe, file_path)
        use hdf5_io, only: open_file_h5, close_file_h5, read_data_h5
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        integer(hsize_t), dimension(1) :: dcount, doffset, dset_dims
        integer :: nptl_local
        integer(i8) :: nptl_global, nptl_offset
        integer, allocatable, dimension(:) :: nptls_local
        character(len=128) :: fname
        character(len=4) :: ctime
        integer(hid_t) :: file_id, dset_id, filespace
        integer :: error, i
        logical :: dir_e

        call h5open_f(error)

        ! Read the local number of particles
        write (ctime,'(i4.4)') iframe
        fname = trim(file_path)//'particles_'//ctime//'.h5'
        allocate(nptls_local(mpi_size))
        if (mpi_rank == master) then
            call open_file_h5(fname, H5F_ACC_RDONLY_F, file_id, .false., MPI_COMM_WORLD)
            call h5dopen_f(file_id, "nptl_local", dset_id, error)
            dcount(1) = mpi_size
            doffset(1) = 0
            dset_dims(1) = mpi_size
            call read_data_h5(dset_id, dcount, doffset, dset_dims, &
                nptls_local, .false., .false.)
            call h5dclose_f(dset_id, error)
            call close_file_h5(file_id)
        endif
        call MPI_BCAST(nptls_local, mpi_size, MPI_INTEGER, master, &
            MPI_COMM_WORLD, ierr)
        nptl_current = nptls_local(mpi_rank+1)

        nptl_global = sum(nptls_local + 0_i8)
        nptl_offset = sum(nptls_local(:mpi_rank+1) + 0_i8) - nptl_current
        deallocate(nptls_local)

        ! Read the particle data
        dcount(1) = nptl_current
        doffset(1) = nptl_offset
        dset_dims(1) = nptl_global

        call open_file_h5(fname, H5F_ACC_RDONLY_F, file_id, .true., MPI_COMM_WORLD)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, "x", ptls%x)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, "y", ptls%y)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, "z", ptls%z)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, "p", ptls%p)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, "v", ptls%v)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, "mu", ptls%mu)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, "weight", ptls%weight)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, "t", ptls%t)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, "dt", ptls%dt)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, &
            "split_times", ptls%split_times)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, &
            "count_flag", ptls%count_flag)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, &
            "origin", ptls%origin)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, &
            "nsteps_tracked", ptls%nsteps_tracked)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, &
            "nsteps_pushed", ptls%nsteps_pushed)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, &
            "tag_injected", ptls%tag_injected)
        call read_ptl_element(file_id, dcount, doffset, dset_dims, &
            "tag_splitted", ptls%tag_splitted)

        call close_file_h5(file_id)
        call h5close_f(error)
    end subroutine read_particles

    !---------------------------------------------------------------------------
    !< Initialize particle tracking. This subroutine will read the information of
    !< the selected particles from a file.
    !< Args:
    !<  particle_tags_file: the file containing the particle tags
    !<  nsteps_interval: save particle points every nsteps_interval
    !---------------------------------------------------------------------------
    subroutine init_particle_tracking(particle_tags_file, nsteps_interval)
        use hdf5_io, only: open_file_h5, close_file_h5, read_data_h5
        implicit none
        character(*), intent(in) :: particle_tags_file
        integer, intent(in) :: nsteps_interval
        integer(hid_t) :: file_id, dset_id, filespace
        integer(hsize_t), dimension(2) :: dset_dims, dset_dims_max
        integer(hsize_t), dimension(2) :: dcount, doffset
        integer :: error

        track_particle_flag = .true.

        call h5open_f(error)

        if (mpi_rank == master) then
            call open_file_h5(particle_tags_file, H5F_ACC_RDONLY_F, &
                file_id, .false., MPI_COMM_WORLD)
            call h5dopen_f(file_id, "tags", dset_id, error)
            call h5dget_space_f(dset_id, filespace, error)
            call h5Sget_simple_extent_dims_f(filespace, dset_dims, dset_dims_max, error)
            call h5sclose_f(filespace, error)
            split_times_max = dset_dims(1) - 2 ! the first col is the origin info
            nptl_tracking = dset_dims(2)
        endif
        call MPI_BCAST(nptl_tracking, 1, MPI_INTEGER, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_BCAST(split_times_max, 1, MPI_INTEGER, master, &
            MPI_COMM_WORLD, ierr)

        ! The first col is the origin info. The second col is tag_injected.
        ! The rest are tag_splitted for every split.
        allocate(tags_tracking(split_times_max+2, nptl_tracking))

        if (mpi_rank == master) then
            dcount(1) = split_times_max + 2
            dcount(2) = nptl_tracking
            doffset(1) = 0
            doffset(2) = 0
            call read_data_h5(dset_id, dcount, doffset, dset_dims, &
                tags_tracking, .false., .false.)
            call h5dclose_f(dset_id, error)
            call close_file_h5(file_id)
        endif

        call h5close_f(error)

        call MPI_BCAST(tags_tracking, nptl_tracking*(split_times_max+2), &
            MPI_INTEGER, master, MPI_COMM_WORLD, ierr)

        ! initialize the tracked particles
        ! The worst case is that all the particles are on the same MPI rank
        nsteps_tracking_max = ceiling((1.0 / dt_min_rel) / nsteps_interval) + 1
        allocate(particles_tracked(nsteps_tracking_max, nptl_tracking))
        call reset_tracked_particles
    end subroutine init_particle_tracking

    !---------------------------------------------------------------------------
    !< Reset the tracked particles
    !---------------------------------------------------------------------------
    subroutine reset_tracked_particles
        implicit none
        particles_tracked%x = 0.0
        particles_tracked%y = 0.0
        particles_tracked%z = 0.0
        particles_tracked%p = 0.0
        particles_tracked%v = 0.0
        particles_tracked%mu = 0.0
        particles_tracked%weight = 0.0
        particles_tracked%t = 0.0
        particles_tracked%dt = 0.0
        particles_tracked%split_times = 0
        particles_tracked%count_flag = 0
        particles_tracked%origin = 0
        particles_tracked%nsteps_tracked = 0
        particles_tracked%nsteps_pushed = 0
        particles_tracked%tag_injected = 0
        particles_tracked%tag_splitted = 0
    end subroutine reset_tracked_particles

    !---------------------------------------------------------------------------
    !< Free particle tracking.
    !---------------------------------------------------------------------------
    subroutine free_particle_tracking
        use hdf5_io, only: open_file_h5, close_file_h5, read_data_h5
        implicit none
        deallocate(tags_tracking)
        deallocate(particles_tracked)
    end subroutine free_particle_tracking

    !---------------------------------------------------------------------------
    !< Check if the particle is selected for tracking
    !< Args:
    !<  ptl: particle information
    !<  iptl_lo: the lower index of the particle if tracked
    !<  iptl_hi: the higher index of the particle if tracked
    !---------------------------------------------------------------------------
    function is_particle_selected(ptl, iptl_lo, iptl_hi) result(is_tagged)
        implicit none
        type(particle_type), intent(in) :: ptl
        integer, intent(out) :: iptl_lo, iptl_hi
        integer :: i1, i2, i3, i4, i5, i6
        integer :: nsplit
        logical :: is_tagged

        is_tagged = .false.
        iptl_lo = -1
        iptl_hi = -1
        nsplit = ptl%split_times
        if (nsplit <= split_times_max) then
            i1 = findloc(tags_tracking(1, :), ptl%origin, dim=1)
            if (i1 > 0) then
                i2 = findloc(tags_tracking(1, :), ptl%origin, dim=1, back=.true.)
                i3 = findloc(tags_tracking(2, i1:i2), abs(ptl%tag_injected), dim=1)
                if (i3 > 0) then
                    i4 = findloc(tags_tracking(2, i1:i2), abs(ptl%tag_injected), dim=1, back=.true.)
                    i3 = i3 + i1 - 1
                    i4 = i4 + i1 - 1
                    if (nsplit > 0) then
                        i5 = findloc(tags_tracking(nsplit+2, i3:i4), abs(ptl%tag_splitted), dim=1)
                        if (i5 > 0) then
                            i6 = findloc(tags_tracking(nsplit+2, i3:i4), abs(ptl%tag_splitted), &
                                dim=1, back=.true.)
                            iptl_lo = i5 + i3 - 1
                            iptl_hi = i6 + i3 - 1
                            is_tagged = .true.
                        endif
                    else
                        iptl_lo = i3
                        iptl_hi = i4
                        is_tagged = .true.
                    endif
                endif
            endif
        endif
    end function is_particle_selected

    !---------------------------------------------------------------------------
    !< Locate the particle in tags_tracking. Since it is only called by tracked
    !< particles, we don't need to check the arguments.
    !< Args:
    !<  ptl: particle information
    !---------------------------------------------------------------------------
    subroutine locate_particle(ptl, iptl_lo, iptl_hi)
        implicit none
        type(particle_type), intent(in) :: ptl
        integer, intent(out) :: iptl_lo, iptl_hi
        integer :: i1, i2, i3, i4, i5, i6
        integer :: nsplit

        i1 = findloc(tags_tracking(1, :), ptl%origin, dim=1)
        i2 = findloc(tags_tracking(1, :), ptl%origin, dim=1, back=.true.)
        i3 = findloc(tags_tracking(2, i1:i2), abs(ptl%tag_injected), dim=1)
        i4 = findloc(tags_tracking(2, i1:i2), abs(ptl%tag_injected), dim=1, back=.true.)
        i3 = i3 + i1 - 1
        i4 = i4 + i1 - 1
        nsplit = ptl%split_times
        if (nsplit > 0) then
            i5 = findloc(tags_tracking(nsplit+2, i3:i4), abs(ptl%tag_splitted), dim=1)
            i6 = findloc(tags_tracking(nsplit+2, i3:i4), abs(ptl%tag_splitted), dim=1, back=.true.)
            iptl_lo = i5 + i3 - 1
            iptl_hi = i6 + i3 - 1
        else
            iptl_lo = i3
            iptl_hi = i4
        endif
    end subroutine locate_particle

    !---------------------------------------------------------------------------
    !< Sync tracked particles across all MPI ranks
    !---------------------------------------------------------------------------
    subroutine sync_tracked_particles
        implicit none
        integer, allocatable, dimension(:) :: nsteps_tracked
        integer, allocatable, dimension(:) :: nsteps_tracked_offset
        integer :: iptl, i, j

        allocate(nsteps_tracked(nptl_tracking))
        allocate(nsteps_tracked_offset(nptl_tracking))

        nsteps_tracked = count(particles_tracked%tag_splitted < 0, dim=1)

        ! Shift all to the left
        do iptl = 1, nptl_tracking
            j = 1
            do i = 1, nsteps_tracking_max
                if (particles_tracked(i, iptl)%tag_splitted < 0) then
                    particles_tracked(j, iptl) = particles_tracked(i, iptl)
                    j = j + 1
                endif
            enddo
            particles_tracked(j:, iptl)%x = 0.0
            particles_tracked(j:, iptl)%y = 0.0
            particles_tracked(j:, iptl)%z = 0.0
            particles_tracked(j:, iptl)%p = 0.0
            particles_tracked(j:, iptl)%v = 0.0
            particles_tracked(j:, iptl)%mu = 0.0
            particles_tracked(j:, iptl)%weight = 0.0
            particles_tracked(j:, iptl)%t = 0.0
            particles_tracked(j:, iptl)%dt = 0.0
            particles_tracked(j:, iptl)%split_times = 0
            particles_tracked(j:, iptl)%count_flag = 0
            particles_tracked(j:, iptl)%origin = 0
            particles_tracked(j:, iptl)%nsteps_tracked = 0
            particles_tracked(j:, iptl)%nsteps_pushed = 0
            particles_tracked(j:, iptl)%tag_injected = 0
            particles_tracked(j:, iptl)%tag_splitted = 0
        enddo

        call MPI_SCAN(nsteps_tracked, nsteps_tracked_offset, nptl_tracking, &
            MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)
        nsteps_tracked_offset = nsteps_tracked_offset - nsteps_tracked

        ! Shift the tracked particles to match the offset
        particles_tracked = cshift(particles_tracked, &
            SHIFT=-nsteps_tracked_offset, dim=1)

        if (mpi_rank == master) then
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%x, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%y, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%z, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%p, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%v, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%mu, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%weight, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%t, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%dt, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%split_times, &
                nsteps_tracking_max*nptl_tracking, MPI_INTEGER1, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%count_flag, &
                nsteps_tracking_max*nptl_tracking, MPI_INTEGER1, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%origin, &
                nsteps_tracking_max*nptl_tracking, MPI_INTEGER, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%nsteps_tracked, &
                nsteps_tracking_max*nptl_tracking, MPI_INTEGER, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%nsteps_pushed, &
                nsteps_tracking_max*nptl_tracking, MPI_INTEGER, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%tag_injected, &
                nsteps_tracking_max*nptl_tracking, MPI_INTEGER, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(MPI_IN_PLACE, particles_tracked%tag_splitted, &
                nsteps_tracking_max*nptl_tracking, MPI_INTEGER, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
        else
            call MPI_REDUCE(particles_tracked%x, particles_tracked%x, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%y, particles_tracked%y, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%z, particles_tracked%z, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%p, particles_tracked%p, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%v, particles_tracked%v, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%mu, particles_tracked%mu, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%weight, particles_tracked%weight, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%t, particles_tracked%t, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%dt, particles_tracked%dt, &
                nsteps_tracking_max*nptl_tracking, MPI_DOUBLE_PRECISION, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%split_times, particles_tracked%split_times, &
                nsteps_tracking_max*nptl_tracking, MPI_INTEGER1, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%count_flag, particles_tracked%count_flag, &
                nsteps_tracking_max*nptl_tracking, MPI_INTEGER1, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%origin, particles_tracked%origin, &
                nsteps_tracking_max*nptl_tracking, MPI_INTEGER, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%nsteps_tracked, particles_tracked%nsteps_tracked, &
                nsteps_tracking_max*nptl_tracking, MPI_INTEGER, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%nsteps_pushed, particles_tracked%nsteps_pushed, &
                nsteps_tracking_max*nptl_tracking, MPI_INTEGER, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%tag_injected, particles_tracked%tag_injected, &
                nsteps_tracking_max*nptl_tracking, MPI_INTEGER, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(particles_tracked%tag_splitted, particles_tracked%tag_splitted, &
                nsteps_tracking_max*nptl_tracking, MPI_INTEGER, &
                MPI_SUM, master, MPI_COMM_WORLD, ierr)
            endif

        deallocate(nsteps_tracked)
        deallocate(nsteps_tracked_offset)
    end subroutine sync_tracked_particles

    !---------------------------------------------------------------------------
    !< Write one double element of the tracked particle data
    !---------------------------------------------------------------------------
    subroutine write_double_tracked_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: write_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(2), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        real(dp), dimension(:, :), intent(in) :: fdata
        integer(hid_t) :: dset_id, filespace
        integer :: error
        call h5screate_simple_f(2, dset_dims, filespace, error)
        call h5dcreate_f(file_id, trim(dset_name), H5T_NATIVE_DOUBLE, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .false., .false.)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)
    end subroutine write_double_tracked_element

    !---------------------------------------------------------------------------
    !< Write one char-length integer element of the tracked particle data
    !---------------------------------------------------------------------------
    subroutine write_integer1_tracked_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: write_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(2), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        integer(i1), dimension(:, :), intent(in) :: fdata
        integer(hid_t) :: dset_id, filespace
        integer :: error
        call h5screate_simple_f(2, dset_dims, filespace, error)
        call h5dcreate_f(file_id, trim(dset_name), H5T_STD_I8LE, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .false., .false.)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)
    end subroutine write_integer1_tracked_element

    !---------------------------------------------------------------------------
    !< Write one default-length integer element of the tracked particle data
    !---------------------------------------------------------------------------
    subroutine write_integer4_tracked_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: write_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(2), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        integer(i4), dimension(:, :), intent(in) :: fdata
        integer(hid_t) :: dset_id, filespace
        integer :: error
        call h5screate_simple_f(2, dset_dims, filespace, error)
        call h5dcreate_f(file_id, trim(dset_name), H5T_NATIVE_INTEGER, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .false., .false.)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)
    end subroutine write_integer4_tracked_element

    !---------------------------------------------------------------------------
    !< Write one long-length integer element of the tracked particle data
    !---------------------------------------------------------------------------
    subroutine write_integer8_tracked_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: write_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(2), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        integer(i8), dimension(:, :), intent(in) :: fdata
        integer(hid_t) :: dset_id, filespace
        integer :: error
        call h5screate_simple_f(2, dset_dims, filespace, error)
        call h5dcreate_f(file_id, trim(dset_name), H5T_STD_I64LE, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .false., .false.)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)
    end subroutine write_integer8_tracked_element

    !---------------------------------------------------------------------------
    !< Sync tracked particles across all MPI ranks
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine dump_tracked_particles(iframe, file_path)
        use hdf5_io, only: create_file_h5, open_file_h5, close_file_h5, write_data_h5
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        integer(hsize_t), dimension(2) :: dcount, doffset, dset_dims
        character(len=256) :: fname
        character(len=4) :: ctime
        integer(hid_t) :: file_id
        integer :: error
        logical :: dir_e

        call sync_tracked_particles

        ! Write the local number of particles
        write (ctime,'(i4.4)') iframe
        fname = trim(file_path)//'particles_tracked_'//ctime//'.h5'
        if (mpi_rank == master) then
            call h5open_f(error)
            call create_file_h5(fname, H5F_ACC_TRUNC_F, file_id, .false., MPI_COMM_WORLD)
            dcount(1) = nsteps_tracking_max
            dcount(2) = nptl_tracking
            doffset(1) = 0
            doffset(2) = 0
            dset_dims(1) = nsteps_tracking_max
            dset_dims(2) = nptl_tracking
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "x", particles_tracked%x)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "y", particles_tracked%y)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "z", particles_tracked%z)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "p", particles_tracked%p)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "v", particles_tracked%v)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "mu", particles_tracked%mu)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "weight", particles_tracked%weight)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "t", particles_tracked%t)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "dt", particles_tracked%dt)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "split_times", particles_tracked%split_times)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "count_flag", particles_tracked%count_flag)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "origin", particles_tracked%origin)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "nsteps_tracked", particles_tracked%nsteps_tracked)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "nsteps_pushed", particles_tracked%nsteps_pushed)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "tag_injected", particles_tracked%tag_injected)
            call write_tracked_ptl_element(file_id, dcount, doffset, dset_dims, &
                "tag_splitted", particles_tracked%tag_splitted)
            call close_file_h5(file_id)
            call h5close_f(error)
        endif

        ! after dumping the data, reset it
        call reset_tracked_particles
    end subroutine dump_tracked_particles
end module particle_module

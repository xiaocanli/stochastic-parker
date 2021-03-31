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
        init_particle_distributions, clean_particle_distributions, &
        free_particle_distributions, distributions_diagnostics, quick_check, &
        set_particle_datatype_mpi, free_particle_datatype_mpi, &
        select_particles_tracking, init_particle_tracking, &
        free_particle_tracking, init_tracked_particle_points, &
        free_tracked_particle_points, negative_particle_tags, &
        save_tracked_particle_points, inject_particles_at_shock, &
        set_mpi_io_data_sizes, init_local_particle_distributions, &
        free_local_particle_distributions, inject_particles_at_large_jz, &
        set_dpp_params, set_flags_params, set_drift_parameters, &
        get_pmax_global, set_flag_check_drift_2d, dump_particles, &
        init_escaped_particles, free_escaped_particles, &
        reset_escaped_particles, dump_escaped_particles, &
        record_tracked_particle_init

    type particle_type
        real(dp) :: x, y, z, p      !< Position and momentum
        real(dp) :: weight, t, dt   !< Particle weight, time and time step
        integer  :: split_times     !< Particle splitting times
        integer  :: count_flag      !< Only count particle when it is 1
        integer  :: tag             !< Particle tag
        integer  :: nsteps_tracking !< Total particle tracking steps
    end type particle_type

    integer, parameter :: COUNT_FLAG_INBOX  = 1  !< For in-box particles
    integer, parameter :: COUNT_FLAG_ESCAPE = -1 !< For escaped particles
    integer, parameter :: COUNT_FLAG_OTHERS = 0  !< For other particles

    integer :: particle_datatype_mpi
    type(particle_type), allocatable, dimension(:) :: ptls
    type(particle_type), allocatable, dimension(:, :) :: senders
    type(particle_type), allocatable, dimension(:, :) :: recvers
    integer, allocatable, dimension(:) :: nsenders, nrecvers
    !dir$ attributes align:64 :: ptls

    ! Particle tracking
    integer, allocatable, dimension(:) :: nsteps_tracked_ptls, noffsets_tracked_ptls
    integer, allocatable, dimension(:) :: tags_selected_ptls
    type(particle_type), allocatable, dimension(:) :: ptl_traj_points
    integer :: nsteps_tracked_tot   !< Total tracking steps for all tracked particles

    integer :: nptl_current         !< Number of particles currently in the box
    integer :: nptl_old             !< Number of particles without receivers
    integer :: nptl_max             !< Maximum number of particles allowed
    integer :: nptl_new             !< Number of particles from splitting
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

    !< Parameters for particle distributions
    real(dp) :: pmin  !< Minimum particle momentum
    real(dp) :: pmax  !< Maximum particle momentum
    real(dp) :: dx_diag, dy_diag, dz_diag
    real(dp) :: pmin_log, pmax_log, dp_log, dp_bands_log
    real(dp) :: pmax_local, pmax_global
    integer :: nx, ny, nz, npp, nreduce
    integer :: nx_mhd_reduced, ny_mhd_reduced, nz_mhd_reduced

    ! Particle distributions
    integer, parameter :: nbands = 12  ! Divide pmin to pmax by nbands
    real(dp), allocatable, dimension(:, :, :, :) :: fbands, fbands_sum
    real(dp), allocatable, dimension(:, :) :: fdpdt, fdpdt_sum
    real(dp), allocatable, dimension(:) :: fp_global, fp_global_sum
    real(dp), allocatable, dimension(:) :: parray, parray_bands
    ! Particle distribution in each cell
    real(dp), allocatable, dimension(:, :, :, :) :: fp_local, fp_local_sum

    !< MPI/IO data sizes
    integer, dimension(4) :: sizes_fxy, subsizes_fxy, starts_fxy
    integer, dimension(4) :: sizes_fp_local, subsizes_fp_local, starts_fp_local

    !< Momentum diffusion
    logical :: dpp_wave_flag, dpp_shear_flag
    !< Whether particle scattering is weak (tau*Omega >> 1)
    logical :: weak_scattering
    !< The scattering time for initial particles (kpara = v^2\tau/3)
    real(dp) :: tau0

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
    !dir$ attributes align:64 :: escaped_ptls

    interface write_ptl_element
        module procedure &
            write_integer_element, write_double_element
    end interface write_ptl_element

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
        ptls%weight = 0.0
        ptls%t = 0.0
        ptls%dt = 0.0
        ptls%split_times = 0
        ptls%count_flag = COUNT_FLAG_OTHERS
        ptls%tag = 0
        ptls%nsteps_tracking = 0
        nptl_current = 0     ! No particle initially
        nptl_new = 0
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
    !< Reset escaped particle data array
    !---------------------------------------------------------------------------
    subroutine reset_escaped_particles
        implicit none
        escaped_ptls%x = 0.0
        escaped_ptls%y = 0.0
        escaped_ptls%z = 0.0
        escaped_ptls%p = 0.0
        escaped_ptls%weight = 0.0
        escaped_ptls%t = 0.0
        escaped_ptls%dt = 0.0
        escaped_ptls%split_times = 0
        escaped_ptls%count_flag = COUNT_FLAG_OTHERS
        escaped_ptls%tag = 0
        escaped_ptls%nsteps_tracking = 0
        nptl_escaped = 0
    end subroutine reset_escaped_particles

    !---------------------------------------------------------------------------
    !< Initialize escaped particle data
    !---------------------------------------------------------------------------
    subroutine init_escaped_particles
        use simulation_setup_module, only: particle_can_escape
        implicit none
        if (particle_can_escape) then
            nptl_escaped_max = nptl_max
        else
            nptl_escaped_max = 1
        endif
        allocate(escaped_ptls(nptl_escaped_max))
        call reset_escaped_particles
    end subroutine init_escaped_particles

    !---------------------------------------------------------------------------
    !< Free escaped particle data
    !---------------------------------------------------------------------------
    subroutine free_escaped_particles
        implicit none
        deallocate(escaped_ptls)
    end subroutine free_escaped_particles

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
    !< Set parameters for momentum diffusion
    !< Args:
    !<  dpp_wave_int(integer): flag for momentum diffusion due to wave scattering
    !<  dpp_shear_int(integer): flag for momentum diffusion due to flow shear
    !<  weak_scattering_flag(integer): flag for weak-scattering regime
    !<  tau0_scattering: the scattering time for initial particles
    !---------------------------------------------------------------------------
    subroutine set_dpp_params(dpp_wave_int, dpp_shear_int, weak_scattering_flag, tau0_scattering)
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
    !< Set other flags and parameters
    !< Args:
    !<  deltab_flag_int(integer): flag for magnetic fluctuation
    !<  correlation_flag_int(integer): flag for turbulence correlation length
    !---------------------------------------------------------------------------
    subroutine set_flags_params(deltab_flag_int, correlation_flag_int)
        implicit none
        integer, intent(in) :: deltab_flag_int, correlation_flag_int
        deltab_flag = .false.
        correlation_flag = .false.
        if (deltab_flag_int == 1) deltab_flag = .true.
        if (correlation_flag_int == 1) correlation_flag = .true.
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
        ! Setup description of the 7 MPI_DOUBLE fields.
        offsets(0) = 0
        oldtypes(0) = MPI_DOUBLE_PRECISION
        blockcounts(0) = 7
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
    !< Inject particles which are spatially uniform
    !< Args:
    !<  nptl: number of particles to be injected
    !<  dt: the time interval
    !<  dist_flag: momentum distribution flag. 0 for Maxwellian, 1 for delta.
    !<  ct_mhd: MHD simulation time frame
    !<  part_box: box to inject particles
    !---------------------------------------------------------------------------
    subroutine inject_particles_spatial_uniform(nptl, dt, dist_flag, ct_mhd, part_box)
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: get_ncells_large_jz
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: nptl, dist_flag, ct_mhd
        real(dp), intent(in) :: dt
        real(dp), intent(in), dimension(6) :: part_box
        integer :: i, imod2
        real(dp) :: xmin, ymin, zmin, xmax, ymax, zmax
        real(dp) :: xmin_box, ymin_box, zmin_box
        real(dp) :: xmax_box, ymax_box, zmax_box
        real(dp) :: rands(2)
        real(dp) :: jz_min, xtmp, ytmp, ztmp
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
            ptls(nptl_current)%x = xtmp
            ptls(nptl_current)%y = ytmp
            ptls(nptl_current)%z = ztmp
            if (dist_flag == 0) then
                imod2 = mod(i, 2)
                if (imod2 == 1) rands = two_normals(0)
                ptls(nptl_current)%p = abs(rands(imod2+1)) * p0
            else
                ptls(nptl_current)%p = p0
            endif
            ptls(nptl_current)%weight = 1.0
            ptls(nptl_current)%t = ct_mhd * mhd_config%dt_out
            ptls(nptl_current)%dt = dt
            ptls(nptl_current)%split_times = 0
            ptls(nptl_current)%count_flag = COUNT_FLAG_INBOX
            tag_max = tag_max + 1
            ptls(nptl_current)%tag = tag_max
            ptls(nptl_current)%nsteps_tracking = 0
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
    !<  dist_flag: momentum distribution flag. 0 for Maxwellian, 1 for delta.
    !<  ct_mhd: MHD simulation time frame
    !---------------------------------------------------------------------------
    subroutine inject_particles_at_shock(nptl, dt, dist_flag, ct_mhd)
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: interp_shock_location
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: nptl, dist_flag, ct_mhd
        real(dp), intent(in) :: dt
        integer :: i, iy, iz, imod2
        real(dp) :: xmin, ymin, xmax, ymax, zmin, zmax
        real(dp) :: ry, rz, dpy, dpz, shock_xpos
        real(dp), dimension(2) :: rands
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
                imod2 = mod(i, 2)
                if (imod2 == 1) rands = two_normals(0)
                ptls(nptl_current)%p = abs(rands(imod2+1)) * p0
            else
                ptls(nptl_current)%p = p0
            endif
            ptls(nptl_current)%weight = 1.0
            ptls(nptl_current)%t = ct_mhd * mhd_config%dt_out
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
        real(dp) :: half_pi

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
        half_pi = pi * 0.5
        ! We assume here theta is from -pi/2 to pi/2
        ctheta1 = cos(ypos_local(pos(2)) + half_pi)
        ctheta2 = cos(ypos_local(pos(2)+1) + half_pi)
        ctheta = cos(y + half_pi)
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
    !<  dist_flag: momentum distribution flag. 0 for Maxwellian, 1 for delta.
    !<  ct_mhd: MHD simulation time frame
    !<  jz_min: the minimum jz
    !<  part_box: box to inject particles
    !---------------------------------------------------------------------------
    subroutine inject_particles_at_large_jz(nptl, dt, dist_flag, ct_mhd, jz_min, part_box)
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: interp_fields
        use mhd_data_parallel, only: get_ncells_large_jz
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: nptl, dist_flag, ct_mhd
        real(dp), intent(in) :: dt, jz_min
        real(dp), intent(in), dimension(6) :: part_box
        real(dp) :: xmin, ymin, zmin, xmax, ymax, zmax
        real(dp) :: xmin_box, ymin_box, zmin_box
        real(dp) :: xmax_box, ymax_box, zmax_box
        real(dp) :: xtmp, ytmp, ztmp, px, py, pz
        real(dp) :: dxm, dym, dzm, dby_dx, dbx_dy
        real(dp) :: rt, jz, by
        real(dp), dimension(2) :: rands
        real(dp), dimension(nfields) :: fields !< Fields at particle position
        real(dp), dimension(ngrads) :: gradf !< Field gradients at particle position
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights
        integer :: i, imod2
        integer :: ncells_large_jz, ncells_large_jz_g

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
        call MPI_ALLREDUCE(ncells_large_jz, ncells_large_jz_g, 1, &
            MPI_INTEGER, MPI_SUM, mpi_sub_comm, ierr)
        nptl_inject = int(nptl * mpi_sub_size * dble(ncells_large_jz) / dble(ncells_large_jz_g))

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
                    call interp_fields(pos, weights, rt, fields, gradf)
                    dbx_dy = gradf(14)
                    dby_dx = gradf(16)
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
            ptls(nptl_current)%x = xtmp
            ptls(nptl_current)%y = ytmp
            ptls(nptl_current)%z = ztmp
            if (dist_flag == 0) then
                imod2 = mod(i, 2)
                if (imod2 == 1) rands = two_normals(0)
                ptls(nptl_current)%p = abs(rands(imod2+1)) * p0
            else
                ptls(nptl_current)%p = p0
            endif
            ptls(nptl_current)%weight = 1.0
            ptls(nptl_current)%t = ct_mhd * mhd_config%dt_out
            ptls(nptl_current)%dt = dt
            ptls(nptl_current)%split_times = 0
            ptls(nptl_current)%count_flag = COUNT_FLAG_INBOX
            tag_max = tag_max + 1
            ptls(nptl_current)%tag = tag_max
            ptls(nptl_current)%nsteps_tracking = 0
        enddo
        if (mpi_rank == master) then
            write(*, "(A)") "Finished injecting particles where jz is large"
        endif
    end subroutine inject_particles_at_large_jz

    !---------------------------------------------------------------------------
    !< Particle mover in one cycle
    !< Args:
    !<  t0: the starting time
    !<  track_particle_flag: whether to track particles, 0 for no, 1 for yes
    !<  nptl_selected: number of selected particles
    !<  nsteps_interval: save particle points every nsteps_interval
    !<  num_fine_steps: number of fine time steps
    !---------------------------------------------------------------------------
    subroutine particle_mover_one_cycle(t0, track_particle_flag, nptl_selected, &
            nsteps_interval, num_fine_steps)
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: interp_fields, interp_magnetic_fluctuation, &
            interp_correlation_length
        implicit none
        real(dp), intent(in) :: t0
        integer, intent(in) :: track_particle_flag, nptl_selected, nsteps_interval
        integer, intent(in) :: num_fine_steps
        real(dp) :: dtf, dxm, dym, dzm, xmin, xmax, ymin, ymax, zmin, zmax
        real(dp) :: xmin1, xmax1, ymin1, ymax1, zmin1, zmax1
        real(dp) :: deltax, deltay, deltaz, deltap
        real(dp) :: dt_target, dt_fine
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights
        real(dp) :: px, py, pz, rt
        integer :: i, tracking_step, offset, step, thread_id
        type(particle_type) :: ptl
        type(kappa_type) :: kappa
        real(dp), dimension(nfields) :: fields !< Fields at particle position
        real(dp), dimension(ngrads) :: gradf !< Field gradients at particle position
        real(dp) :: db2, lc !< turbulence variance and correlation length
        real(dp), dimension(3) :: grad_db2, grad_lc !< Gradients of db2 and lc
        !dir$ attributes align:64 :: fields
        !dir$ attributes align:64 :: gradf

        dtf = mhd_config%dt_out
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
        !$OMP& fields, gradf, db2, lc, grad_db2, grad_lc, &
        !$OMP& deltax, deltay, deltaz, deltap, &
        !$OMP& dt_target, pos, weights, px, py, pz, rt, &
        !$OMP& tracking_step, offset, step, thread_id)
        thread_id = OMP_GET_THREAD_NUM()
        !$OMP DO
        do i = nptl_old + 1, nptl_current
            ptl = ptls(i)
            deltax = 0.0
            deltay = 0.0
            deltaz = 0.0
            deltap = 0.0

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
                    call interp_fields(pos, weights, rt, fields, gradf)
                    if (deltab_flag) then
                        call interp_magnetic_fluctuation(pos, weights, rt, db2, grad_db2)
                    endif
                    if (correlation_flag) then
                        call interp_correlation_length(pos, weights, rt, lc, grad_lc)
                    endif
                    call calc_spatial_diffusion_coefficients(ptl, fields, gradf, &
                        db2, grad_db2, lc, grad_lc, kappa)
                    call set_time_step(t0, dt_target, fields, gradf, kappa, ptl)
                    if (ndim_field == 1) then
                        call push_particle_1d(thread_id, rt, ptl, fields, gradf, db2, &
                            grad_db2, lc, grad_lc, kappa, deltax, deltap)
                    else if (ndim_field == 2) then
                        call push_particle_2d(thread_id, rt, ptl, fields, gradf, db2, &
                            grad_db2, lc, grad_lc, kappa, deltax, deltay, deltap)
                    else
                        call push_particle_3d(thread_id, rt, ptl, fields, gradf, db2, &
                            grad_db2, lc, grad_lc, kappa, deltax, deltay, deltaz, deltap)
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
                        call interp_fields(pos, weights, rt, fields, gradf)
                        if (deltab_flag) then
                            call interp_magnetic_fluctuation(pos, weights, rt, db2, grad_db2)
                        endif
                        if (correlation_flag) then
                            call interp_correlation_length(pos, weights, rt, lc, grad_lc)
                        endif
                        call calc_spatial_diffusion_coefficients(ptl, fields, gradf, &
                            db2, grad_db2, lc, grad_lc, kappa)
                        if (ndim_field == 1) then
                            call push_particle_1d(thread_id, rt, ptl, fields, gradf, db2, &
                                grad_db2, lc, grad_lc, kappa, deltax, deltap)
                        else if (ndim_field == 2) then
                            call push_particle_2d(thread_id, rt, ptl, fields, gradf, db2, &
                                grad_db2, lc, grad_lc, kappa, deltax, deltay, deltap)
                        else
                            call push_particle_3d(thread_id, rt, ptl, fields, gradf, db2, &
                                grad_db2, lc, grad_lc, kappa, deltax, deltay, deltaz, deltap)
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

                if (ptl%count_flag == COUNT_FLAG_INBOX) then
                    call energization_dist(ptl, t0)
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
    !<  track_particle_flag: whether to track particles, 0 for no, 1 for yes
    !<  nptl_selected: number of selected particles
    !<  nsteps_interval: save particle points every nsteps_interval
    !<  mhd_tframe: MHD time frame
    !<  num_fine_steps: number of fine time steps
    !---------------------------------------------------------------------------
    subroutine particle_mover(track_particle_flag, nptl_selected, &
            nsteps_interval, mhd_tframe, num_fine_steps)
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config
        implicit none
        integer, intent(in) :: track_particle_flag, nptl_selected, nsteps_interval
        integer, intent(in) :: mhd_tframe, num_fine_steps
        integer :: i, local_flag, global_flag, ncycle
        logical :: all_particles_in_box
        real(dp) :: t0, xmin, xmax, ymin, ymax, zmin, zmax
        type(particle_type) :: ptl
        all_particles_in_box = .false.
        nptl_old = 0
        nptl_new = 0

        t0 = mhd_config%dt_out * mhd_tframe

        ncycle = 0
        local_flag = 0
        global_flag = 0

        do while (.not. all_particles_in_box)
            ncycle = ncycle + 1
            nsenders = 0
            nrecvers = 0
            if (nptl_old < nptl_current) then
                call particle_mover_one_cycle(t0, track_particle_flag, &
                    nptl_selected, nsteps_interval, num_fine_steps)
            endif
            call remove_particles
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
        call remove_particles
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

        if (ndim_field == 3) then
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
    !<  fields: fields at particle position
    !<  gradf: gradients of the fields at particle position
    !<  db2: turbulence variance at particle position
    !<  grad_db2: gradients of db2
    !<  lc: turbulence correlation length at particle position
    !<  grad_lc: gradients of lc
    !<  kappa: kappa and related variables
    !---------------------------------------------------------------------------
    subroutine calc_spatial_diffusion_coefficients(ptl, fields, gradf, db2, &
            grad_db2, lc, grad_lc, kappa)
        implicit none
        type(particle_type), intent(in) :: ptl
        real(dp), dimension(*), intent(in) :: fields
        real(dp), dimension(*), intent(in) :: gradf
        real(dp), dimension(*), intent(in) :: grad_db2
        real(dp), dimension(*), intent(in) :: grad_lc
        real(dp), intent(in) :: db2, lc
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
            kappa%knorm0 = kappa%knorm0 / db2
        endif

        ! Turbulence correlation length
        ! Make sure that lc is non-zero in the data file!!!
        if (correlation_flag) then
            kappa%knorm0 = kappa%knorm0 * lc**(2./3.)
        endif

        ! Momentum-dependent kappa
        if (momentum_dependency == 1) then
            knorm = kappa%knorm0 * (ptl%p / p0)**pindex
        else
            knorm = kappa%knorm0
        endif

        kappa%kpara = kpara0 * knorm
        kappa%kperp = kappa%kpara * kret

        if (deltab_flag) then
            kappa%kperp = kappa%kperp * db2
        endif

        kappa%skpara = dsqrt(2.0 * kappa%kpara)
        kappa%skperp = dsqrt(2.0 * kappa%kperp)
        kappa%skpara_perp = dsqrt(2.0 * (kappa%kpara - kappa%kperp))

        if (ndim_field == 1) then
            dkdx = 0.0_dp
            if (mag_dependency == 1) then
                dkdx = -db_dx * ib1 / 3
            endif
            if (deltab_flag) then
                dkdx = dkdx - grad_db2(1) / db2
            endif
            if (correlation_flag) then
                dkdx = dkdx + 2 * grad_lc(1) / (3 * lc)
            endif
            kappa%dkxx_dx = kappa%kpara * dkdx
            if (spherical_coord_flag) then
                kappa%kxx = kappa%kpara
            endif
        else if (ndim_field == 2) then
            dbx_dx = gradf(13)
            dbx_dy = gradf(14)
            dby_dx = gradf(16)
            dby_dy = gradf(17)
            db_dx = gradf(22)
            db_dy = gradf(23)
            dkdx = 0.0_dp
            dkdy = 0.0_dp
            if (mag_dependency == 1) then
                dkdx = -db_dx * ib1 / 3
                dkdy = -db_dy * ib1 / 3
            endif
            if (deltab_flag) then
                dkdx = dkdx - grad_db2(1) / db2
                dkdy = dkdy - grad_db2(2) / db2
            endif
            if (correlation_flag) then
                dkdx = dkdx + 2 * grad_lc(1) / (3 * lc)
                dkdy = dkdy + 2 * grad_lc(2) / (3 * lc)
            endif
            kpp = kappa%kpara - kappa%kperp
            kappa%dkxx_dx = kappa%kperp*dkdx + kpp*dkdx*bx**2*ib2 + &
                2.0*kpp*bx*(dbx_dx*b-bx*db_dx)*ib3
            kappa%dkyy_dy = kappa%kperp*dkdy + kpp*dkdy*by**2*ib2 + &
                2.0*kpp*by*(dby_dy*b-by*db_dy)*ib3
            kappa%dkxy_dx = kpp*dkdx*bx*by*ib2 + kpp * &
                ((dbx_dx*by+bx*dby_dx)*ib2 - 2.0*bx*by*db_dx*ib3)
            kappa%dkxy_dy = kpp*dkdy*bx*by*ib2 + kpp * &
                ((dbx_dy*by+bx*dby_dy)*ib2 - 2.0*bx*by*db_dy*ib3)
            if (spherical_coord_flag) then
                kappa%kxx = kappa%kperp + kpp * bx * bx * ib2
                kappa%kyy = kappa%kperp + kpp * by * by * ib2
                kappa%kxy = kpp * bx * by * ib2
            endif
        else
            dbx_dx = gradf(13)
            dbx_dy = gradf(14)
            dbx_dz = gradf(15)
            dby_dx = gradf(16)
            dby_dy = gradf(17)
            dby_dz = gradf(18)
            dbz_dx = gradf(19)
            dbz_dy = gradf(20)
            dbz_dz = gradf(21)
            db_dx = gradf(22)
            db_dy = gradf(23)
            db_dz = gradf(24)
            dkdx = 0.0_dp
            dkdy = 0.0_dp
            dkdz = 0.0_dp
            if (mag_dependency == 1) then
                dkdx = -db_dx * ib1 / 3
                dkdy = -db_dy * ib1 / 3
                dkdz = -db_dz * ib1 / 3
            endif
            if (deltab_flag) then
                dkdx = dkdx - grad_db2(1) / db2
                dkdy = dkdy - grad_db2(2) / db2
                dkdz = dkdz - grad_db2(3) / db2
            endif
            if (correlation_flag) then
                dkdx = dkdx + 2 * grad_lc(1) / (3 * lc)
                dkdy = dkdy + 2 * grad_lc(2) / (3 * lc)
                dkdz = dkdz + 2 * grad_lc(3) / (3 * lc)
            endif
            kpp = kappa%kpara - kappa%kperp
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
            if (spherical_coord_flag) then
                kappa%kxx = kappa%kperp + kpp * bx * bx * ib2
                kappa%kyy = kappa%kperp + kpp * by * by * ib2
                kappa%kzz = kappa%kperp + kpp * bz * bz * ib2
                kappa%kxy = kpp * bx * by * ib2
                kappa%kxz = kpp * bx * bz * ib2
                kappa%kyz = kpp * by * bz * ib2
            endif
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
        use simulation_setup_module, only: mpi_sizex, mpi_sizey
        use mhd_config_module, only: mhd_config
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
            dt_min = get_variable(fh, 'dt_min', '=')
            dt_min_rel = get_variable(fh, 'dt_min_rel', '=')
            dt_max_rel = get_variable(fh, 'dt_max_rel', '=')

            dt_min = max(dt_min, dt_min_rel * mhd_config%dt_out)
            dt_max = dt_max_rel * mhd_config%dt_out

            temp = get_variable(fh, 'nreduce', '=')
            nreduce = int(temp)
            nx = (fconfig%nx + nreduce - 1) / nreduce
            ny = (fconfig%ny + nreduce - 1) / nreduce
            nz = (fconfig%nz + nreduce - 1) / nreduce
            temp = get_variable(fh, 'npp', '=')
            npp = int(temp)
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
            write(*, "(A,E13.6E2)") " Minimum time step = ", dt_min
            write(*, "(A,E13.6E2)") " maximum time step = ", dt_max
            write(*, "(A,I0,A,I0,A,I0)") " Dimensions of spatial distributions = ", &
                nx, " ", ny, " ", nz
            write(*, "(A,I0)") " Dimensions of momentum distributions = ", npp
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
        call MPI_BCAST(dt_min, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(dt_max, 1, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nreduce, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nx, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ny, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nz, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(npp, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(acc_region_flag, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)

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

        if (ndim_field == 1) then
            if (nx * nreduce /= fconfig%nx) then
                if (mpi_rank == master) then
                    write(*, "(A)") "Wrong factor 'nreduce' for particle distribution"
                endif
                call MPI_FINALIZE(ierr)
                stop
            endif
        else if (ndim_field == 2) then
            if (nx * nreduce /= fconfig%nx .or. &
                ny * nreduce /= fconfig%ny) then
                if (mpi_rank == master) then
                    write(*, "(A)") "Wrong factor 'nreduce' for particle distribution"
                endif
                call MPI_FINALIZE(ierr)
                stop
            endif
        else
            if (nx * nreduce /= fconfig%nx .or. &
                ny * nreduce /= fconfig%ny .or. &
                nz * nreduce /= fconfig%nz) then
                if (mpi_rank == master) then
                    write(*, "(A)") "Wrong factor 'nreduce' for particle distribution"
                endif
                call MPI_FINALIZE(ierr)
                stop
            endif
        endif
    end subroutine read_particle_params

    !---------------------------------------------------------------------------
    !< Determine the time step. X, Y, Z are for r, theta, phi in
    !< spherical coordinates
    !< Args;
    !<  t0: the initial time for current particle
    !<  dtf: the time interval between fine diagnostics
    !<  fields: fields at particle position
    !<  gradf: gradients of the fields at particle position
    !<  kappa: kappa and related variables
    !<  ptl: particle structure
    !---------------------------------------------------------------------------
    subroutine set_time_step(t0, dtf, fields, gradf, kappa, ptl)
        use constants, only: pi
        use mhd_config_module, only: mhd_config
        implicit none
        real(dp), intent(in) :: t0, dtf
        real(dp), dimension(*), intent(in) :: fields
        real(dp), dimension(*), intent(in) :: gradf
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
                    ptl%dt = min(dxm/(80.0*tmp40), (tmp30/tmp40)**2) * 0.5d0
                else
                    ptl%dt = dxm/(80.0*tmp40)
                endif
            else
                ptl%dt = dt_min
            endif
        else if (ndim_field == 2) then
            vx = fields(1)
            vy = fields(2)
            bx = fields(5)
            by = fields(6)
            b = dsqrt(bx**2 + by**2)
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
                ctheta = cos(ptl%y + 0.5*pi)
                istheta = 1.0 / sin(ptl%y + 0.5*pi)
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
                    (2.0*kappa%kxx + kappa%dkxy_dy*ir + kappa%kxy*ctheta*istheta)*ir)
            else
                tmp30 = kappa%skperp + kappa%skpara_perp * abs(bxn)
                tmp40 = abs(vx + kappa%dkxx_dx + kappa%dkxy_dy)
            endif
            if (tmp40 .ne. 0.0d0) then
                if (tmp30 > 0) then
                    dt1 = min(dxm/(80.0*tmp40), (tmp30/tmp40)**2) * 0.5d0
                else
                    dt1 = dxm/(80.0*tmp40)
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
                    dt2 = min(dym/(80.0*tmp40), (tmp30/tmp40)**2) * 0.5d0
                else
                    dt2 = dym/(80.0*tmp40)
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
                ctheta = cos(ptl%y + 0.5*pi)
                istheta = 1.0 / sin(ptl%y + 0.5*pi)
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
            dbx_dy = gradf(14)
            dbx_dz = gradf(15)
            dby_dx = gradf(16)
            dby_dz = gradf(18)
            dbz_dx = gradf(19)
            dbz_dy = gradf(20)
            db_dx = gradf(22)
            db_dy = gradf(23)
            db_dz = gradf(24)
            ib2 = ib * ib
            ib3 = ib * ib2
            vdp = 1.0 / (3 * pcharge) / dsqrt((drift1*p0/ptl%p)**2 + (drift2*p0**2/ptl%p**2)**2)
            if (spherical_coord_flag) then
                gbr = db_dx
                gbt = db_dy * ir
                gbp = db_dz * ir * istheta
                vdx = vdp * ((dbz_dy + bz*ctheta*istheta - dby_dz*istheta)*ir*ib2 - &
                             (gbt*bz - gbp*by*istheta)*2.0*ir*ib3)
                vdy = vdp * ((dbx_dz*istheta*ir - dbz_dx - bz*ir)*ib2 - &
                             (gbp*bx*istheta*ir - gbr*bz)*2.0*ib3)
                vdz = vdp * ((dby_dx + by*ir - dbx_dy*ir)*ib2 - &
                             (gbr*by - gbt*bx*ir)*2.0*ib3)
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
                    dt1 = min(dxm/(80.0*tmp40), (tmp30/tmp40)**2) * 0.5d0
                else
                    dt1 = dxm/(80.0*tmp40)
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
                    dt2 = min(dym/(80.0*tmp40), (tmp30/tmp40)**2) * 0.5d0
                else
                    dt2 = dym/(80.0*tmp40)
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
                    dt3 = min(dym/(80.0*tmp40), (tmp30/tmp40)**2) * 0.5d0
                else
                    dt3 = dym/(80.0*tmp40)
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
    !< Push particle for a single step for a 1D simulation
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
        inx = (xnorm >= acc_region(1)) .and. ((1.0d0-xnorm) <= acc_region(2))
        if (ndim_field > 1) then
            ynorm = (ptl%y - mhd_config%ymin) / mhd_config%ly
            iny = (ynorm >= acc_region(3)) .and. ((1.0d0-ynorm) <= acc_region(4))
        endif
        if (ndim_field == 3) then
            znorm = (ptl%z - mhd_config%zmin) / mhd_config%lz
            inz = (znorm >= acc_region(5)) .and. ((1.0d0-znorm) <= acc_region(6))
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
    !<  fields: fields at particle position
    !<  gradf: gradients of the fields at particle position
    !<  db2: turbulence variance at particle position
    !<  grad_db2: gradients of db2
    !<  lc: turbulence correlation length at particle position
    !<  grad_lc: gradients of lc
    !<  kappa: kappa and related variables
    !<  deltax, deltap: the change of x and p in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_1d(thread_id, rt, ptl, fields, gradf, db2, &
            grad_db2, lc, grad_lc, kappa, deltax, deltap)
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
        real(dp), dimension(*), intent(inout) :: gradf
        real(dp), dimension(*), intent(inout) :: grad_db2
        real(dp), dimension(*), intent(inout) :: grad_lc
        real(dp), intent(inout) :: db2, lc
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
        dvx_dx = gradf(1)
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
            call interp_fields(pos, weights, rt, fields, gradf)
            if (deltab_flag) then
                call interp_magnetic_fluctuation(pos, weights, rt, db2, grad_db2)
            endif
            if (correlation_flag) then
                call interp_correlation_length(pos, weights, rt, lc, grad_lc)
            endif

            skpara = kappa%skpara

            call calc_spatial_diffusion_coefficients(ptl, fields, gradf, &
                db2, grad_db2, lc, grad_lc, kappa)

            skpara1 = kappa%skpara ! Diffusion coefficient at predicted position

            deltax = deltax + ran1*skpara*sdt + &
                (skpara1-skpara)*(ran1*ran1-1.0)*sdt/2.0
        endif

        ptl%x = ptl%x + deltax
        ptl%t = ptl%t + ptl%dt

        ! Momentum
        if (spherical_coord_flag) then
            divv = dvx_dx
        else
            divv = dvx_dx + 2.0*vx/ptl%x
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
                ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
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
    end subroutine push_particle_1d

    !---------------------------------------------------------------------------
    !< Push particle for a single step for a 2D simulation
    !< Args:
    !<  thread_id: thread ID staring from 0
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  ptl: particle structure
    !<  fields: fields at particle position
    !<  gradf: gradients of the fields at particle position
    !<  db2: turbulence variance at particle position
    !<  grad_db2: gradients of db2
    !<  lc: turbulence correlation length at particle position
    !<  grad_lc: gradients of lc
    !<  kappa: kappa and related variables
    !<  deltax, deltay, deltap: the change of x, y and p in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_2d(thread_id, rt, ptl, fields, gradf, db2, &
            grad_db2, lc, grad_lc, kappa, deltax, deltay, deltap)
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
        real(dp), dimension(*), intent(inout) :: gradf
        real(dp), dimension(*), intent(inout) :: grad_db2
        real(dp), dimension(*), intent(inout) :: grad_lc
        real(dp), intent(inout) :: db2, lc
        type(kappa_type), intent(inout) :: kappa
        real(dp), intent(out) :: deltax, deltay, deltap
        real(dp) :: xtmp, ytmp
        real(dp) :: sdt, dvx_dx, dvx_dy, dvy_dx, dvy_dy, divv, gshear
        real(dp) :: bx, by, bz, b, vx, vy, px, py, pz, rt1, ib
        real(dp) :: bx1, by1, btot1, ib1
        real(dp) :: dbx_dy, dby_dx, db_dx, db_dy, vdz, vdp, ib2, ib3, gbr, gbt
        real(dp) :: xmin, ymin, xmax, ymax, dxm, dym
        reaL(dp) :: xmin1, ymin1, xmax1, ymax1, dxmh, dymh
        real(dp) :: skperp, skpara_perp, skperp1, skpara_perp1
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rho, va ! Plasma density and Alfven speed
        real(dp) :: rands(2)
        real(dp) :: a1, b1, c1, Qpp, Qpm, Qmp, Qmm, ctheta, istheta, ir, ir2
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
        dvx_dx = gradf(1)
        dvy_dy = gradf(5)
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
            ctheta = cos(ptl%y + 0.5*pi)
            istheta = 1.0 / sin(ptl%y + 0.5*pi)
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

        !< Check particle drift along the out-of-plane direction
        if (check_drift_2d) then
            dbx_dy = gradf(14)
            dby_dx = gradf(16)
            db_dx = gradf(22)
            db_dy = gradf(23)
            ib2 = ib * ib
            ib3 = ib * ib2
            vdp = 1.0 / (3 * pcharge) / dsqrt((drift1*p0/ptl%p)**2 + (drift2*p0**2/ptl%p**2)**2)
            if (spherical_coord_flag) then
                gbr = db_dx
                gbt = db_dy * ir
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
            deltax = (vx + kappa%dkxx_dx + &
                (2.0*kappa%kxx + kappa%dkxy_dy*ir + kappa%kxy*ctheta*istheta)*ir)*ptl%dt
            deltay = ((vy + kappa%dkxy_dx)*ir + &
                (kappa%kxy + kappa%dkyy_dy + kappa%kyy*ctheta*istheta)*ir2)*ptl%dt
            xtmp = ptl%x + deltax + (-Qmm*qtmp1 + Qpm*qtmp2)*sdt
            ytmp = ptl%y + deltay + (2*b1*qtmp1 + 2*b1*qtmp2)*sdt
        else
            deltax = (vx + kappa%dkxx_dx + kappa%dkxy_dy)*ptl%dt
            deltay = (vy + kappa%dkxy_dx + kappa%dkyy_dy)*ptl%dt
            xtmp = ptl%x + deltax + (kappa%skperp + kappa%skpara_perp*bx*ib)*sdt
            ytmp = ptl%y + deltay + (kappa%skperp + kappa%skpara_perp*by*ib)*sdt
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

        !< We originally tried to decrease the time step when xtmp or ytmp are out-of-bound,
        !< but decreasing the time step does not necessarily make the moving distance smaller.
        !< Therefore, we switch between first-order and second-order method.
        if (xtmp < xmin1 .or. xtmp > xmax1 .or. ytmp < ymin1 .or. ytmp > ymax1) then
            !< First-order method
            if (spherical_coord_flag) then
                deltax = deltax + (-Qmm*qtmp1*ran1 + Qpm*qtmp2*ran2)*sdt
                deltay = deltay + (2*b1*qtmp1*ran1 + 2*b1*qtmp2*ran2)*sdt
            else
                deltax = deltax + ran1*kappa%skperp*sdt + ran3*kappa%skpara_perp*sdt*bx*ib
                deltay = deltay + ran2*kappa%skperp*sdt + ran3*kappa%skpara_perp*sdt*by*ib
            endif
        else
            !< Second-order method. It requires xtmp and ytmp are in the local domain.
            px = (xtmp - xmin) / dxm
            py = (ytmp - ymin) / dym
            if (spherical_coord_flag) then
                call get_interp_paramters_spherical(ptl%x, ptl%y, 0.0_dp, pos, weights)
            else
                call get_interp_paramters(px, py, 0.0_dp, pos, weights)
            endif
            call interp_fields(pos, weights, rt, fields, gradf)
            if (deltab_flag) then
                call interp_magnetic_fluctuation(pos, weights, rt, db2, grad_db2)
            endif
            if (correlation_flag) then
                call interp_correlation_length(pos, weights, rt, lc, grad_lc)
            endif

            skperp = kappa%skperp
            skpara_perp = kappa%skpara_perp

            call calc_spatial_diffusion_coefficients(ptl, fields, gradf, &
                db2, grad_db2, lc, grad_lc, kappa)

            if (spherical_coord_flag) then
                ir = 1.0 / xtmp
                ir2 = ir * ir
                a1_1 = kappa%kxx
                b1_1 = kappa%kxy * ir
                c1_1 = kappa%kyy * ir2
                atmp = dsqrt((a1_1-c1_1)**2 + 4*b1_1**2)
                Qpp_1 = atmp + (a1_1 + c1_1)
                Qmp_1 = atmp - (a1_1 + c1_1)
                Qpm_1 = atmp + (a1_1 - c1_1)
                Qmm_1 = atmp - (a1_1 - c1_1)
                qtmp1_1 = dsqrt(-Qmp_1/(Qmm_1**2+4*b1_1**2))
                qtmp2_1 = dsqrt(Qpp_1/(Qpm_1**2+4*b1_1**2))
            endif

            !< Magnetic field at the predicted position
            bx1 = fields(5)
            by1 = fields(6)
            btot1 = fields(8)

            skperp1 = kappa%skperp
            skpara_perp1 = kappa%skpara_perp

            if (btot1 == 0) then
                ib1 = 0.0
            else
                ib1 = 1.0 / btot1
            endif

            if (spherical_coord_flag) then
                deltax = deltax - Qmm*qtmp1*ran1*sdt + Qpm*qtmp2*ran2*sdt + &
                         (-Qmm_1*qtmp1_1 + Qmm*qtmp1)*(ran1*ran1-1.0)*sdt/2.0 + &
                         (Qpm_1*qtmp2_1 - Qpm*qtmp2)*(ran2*ran2-1.0)*sdt/2.0
                deltay = deltay + 2*b1*qtmp1*ran1*sdt + 2*b1*qtmp2*ran2*sdt + &
                         (2*b1_1*qtmp1_1 - 2*b1_1*qtmp1_1)*(ran1*ran1-1.0)*sdt/2.0 + &
                         (2*b1_1*qtmp2_1 - 2*b1_1*qtmp2_1)*(ran2*ran2-1.0)*sdt/2.0
            else
                deltax = deltax + ran1*kappa%skperp*sdt + ran3*skpara_perp*sdt*bx*ib + &
                         (skperp1-skperp)*(ran1*ran1-1.0)*sdt/2.0 + &
                         (skpara_perp1*bx1*ib1-skpara_perp*bx*ib)*(ran3*ran3-1.0)*sdt/2.0
                deltay = deltay + ran2*kappa%skperp*sdt + ran3*skpara_perp*sdt*by*ib + &
                         (skperp1-skperp)*(ran2*ran2-1.0)*sdt/2.0 + &
                         (skpara_perp1*by1*ib1-skpara_perp*by*ib)*(ran3*ran3-1.0)*sdt/2.0
            endif
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
            if (momentum_dependency == 1) then
                deltap = deltap + (8*ptl%p / (27*kappa%kpara)) * va**2 * ptl%dt
            else
                deltap = deltap + (4*ptl%p / (9*kappa%kpara)) * va**2 * ptl%dt
            endif
            ! ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
            rands = two_normals(thread_id)
            ran1 = rands(1)
            deltap = deltap + ran1 * va * ptl%p * dsqrt(2/(9*kappa%kpara)) * sdt
        endif

        ! Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            dvx_dy = gradf(2)
            dvy_dx = gradf(4)
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
                ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
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
    end subroutine push_particle_2d

    !---------------------------------------------------------------------------
    !< Push particle for a single step for a 3D simulation
    !< Args:
    !<  thread_id: thread ID staring from 0
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  ptl: particle structure
    !<  fields: fields at particle position
    !<  gradf: gradients of the fields at particle position
    !<  db2: turbulence variance at particle position
    !<  grad_db2: gradients of db2
    !<  lc: turbulence correlation length at particle position
    !<  grad_lc: gradients of lc
    !<  kappa: kappa and related variables
    !<  deltax, deltay, deltaz, deltap: the change of x, y, z and p in this step
    !---------------------------------------------------------------------------
    subroutine push_particle_3d(thread_id, rt, ptl, fields, gradf, db2, &
            grad_db2, lc, grad_lc, kappa, deltax, deltay, deltaz, deltap)
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
        real(dp), dimension(*), intent(inout) :: gradf
        real(dp), dimension(*), intent(inout) :: grad_db2
        real(dp), dimension(*), intent(inout) :: grad_lc
        real(dp), intent(inout) :: db2, lc
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
        dvx_dx = gradf(1)
        dvy_dy = gradf(5)
        dvz_dz = gradf(9)
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
            ctheta = cos(ptl%y + 0.5*pi)
            istheta = 1.0 / sin(ptl%y + 0.5*pi)
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
        dbx_dy = gradf(14)
        dbx_dz = gradf(15)
        dby_dx = gradf(16)
        dby_dz = gradf(18)
        dbz_dx = gradf(19)
        dbz_dy = gradf(20)
        db_dx = gradf(22)
        db_dy = gradf(23)
        db_dz = gradf(24)
        ib2 = ib * ib
        ib3 = ib * ib2
        vdp = 1.0 / (3 * pcharge) / dsqrt((drift1*p0/ptl%p)**2 + (drift2*p0**2/ptl%p**2)**2)
        if (spherical_coord_flag) then
            gbr = db_dx
            gbt = db_dy * ir
            gbp = db_dz * ir * istheta
            vdx = vdp * ((dbz_dy + bz*ctheta*istheta - dby_dz*istheta)*ir*ib2 - &
                         (gbt*bz - gbp*by*istheta)*2.0*ir*ib3)
            vdy = vdp * ((dbx_dz*istheta*ir - dbz_dx - bz*ir)*ib2 - &
                         (gbp*bx*istheta*ir - gbr*bz)*2.0*ib3)
            vdz = vdp * ((dby_dx + by*ir - dbx_dy*ir)*ib2 - &
                         (gbr*by - gbt*bx*ir)*2.0*ib3)
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
            dvx_dy = gradf(2)
            dvx_dz = gradf(3)
            dvy_dx = gradf(4)
            dvy_dz = gradf(6)
            dvz_dx = gradf(7)
            dvz_dy = gradf(8)
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
            ran1 = (2.0_dp*unif_01(thread_id) - 1.0_dp) * sqrt3
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
    end subroutine push_particle_3d

    !---------------------------------------------------------------------------
    !< Remove particles from simulation if their count_flags are 0.
    !---------------------------------------------------------------------------
    subroutine remove_particles
        implicit none
        integer :: i, nremoved, ntail
        type(particle_type) :: ptl1

        ! Resize the escaped particle data array if necessary
        if ((nptl_escaped + nptl_current) > nptl_escaped_max) then
            call resize_escaped_particles
        endif

        if (nptl_current > 0) then
            nremoved = 0
            ntail = nptl_current - 1
            i = 1
            do while (i < nptl_current)
                if (ntail == (nremoved - 1)) then
                    exit
                endif
                if (ptls(i)%count_flag == COUNT_FLAG_INBOX) then
                    ntail = ntail - 1
                    i = i + 1
                else
                    ! Copy escaped particles
                    if (ptls(i)%count_flag == COUNT_FLAG_ESCAPE) then
                        nptl_escaped = nptl_escaped + 1
                        escaped_ptls(nptl_escaped) = ptls(i)
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
    !---------------------------------------------------------------------------
    subroutine split_particle
        implicit none
        integer :: i, nptl
        real(dp) :: p_threshold
        type(particle_type) :: ptl, ptl_new
        nptl = nptl_current
        do i = 1, nptl
            ptl = ptls(i)
            p_threshold = (1+2.72**ptl%split_times)*p0
            if (ptl%p > p_threshold .and. ptl%p <= pmax) then
                nptl_current = nptl_current + 1
                if (nptl_current > nptl_max) then
                    nptl_current = nptl_max
                    return
                endif
                nptl_new = nptl_new + 1
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
    !< Quick diagnostics of the particle number information
    !< Args:
    !<  iframe: the time frame
    !<  if_create_file: whether to create a file
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine quick_check(iframe, if_create_file, file_path)
        implicit none
        integer, intent(in) :: iframe
        logical, intent(in) :: if_create_file
        character(*), intent(in) :: file_path
        integer, parameter :: nvar = 6
        real(dp) :: pdt_min, pdt_max, pdt_min_g, pdt_max_g
        real(dp), dimension(nvar) :: var_local, var_global
        integer :: i
        logical :: dir_e

        inquire(file='./data/.', exist=dir_e)
        if (.not. dir_e) then
            call system('mkdir -p ./data')
        endif
        var_local = 0.0_dp
        var_global = 0.0_dp
        pdt_min = 1.0_dp
        pdt_max = 0.0_dp
        pdt_min_g = 0.0_dp
        pdt_max_g = 0.0_dp
        var_local(1) = nptl_current
        var_local(2) = nptl_new
        var_local(4) = leak
        var_local(5) = leak_negp
        do i = 1, nptl_current
            var_local(3) = var_local(3) + ptls(i)%weight
            if (ptls(i)%dt < pdt_min) pdt_min = ptls(i)%dt
            if (ptls(i)%dt > pdt_max) pdt_max = ptls(i)%dt
            var_local(6) = var_local(6) + ptls(i)%dt
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(var_local, var_global, nvar, MPI_DOUBLE_PRECISION, MPI_SUM, &
            master, MPI_COMM_WORLD, ierr)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(pdt_min, pdt_min_g, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
            master, MPI_COMM_WORLD, ierr)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(pdt_max, pdt_max_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
            master, MPI_COMM_WORLD, ierr)
        if (mpi_rank == master) then
            if (var_global(1) > 0) then
                var_global(6) = var_global(6) / var_global(1)
            else
                var_global(6) = 0.0_dp
            endif
            if (if_create_file) then
                open (17, file=trim(file_path)//'quick.dat', status='unknown')
                write(17, "(A6,8A13)") "iframe", "nptl_current", "nptl_new", &
                    "ntot", "leak", "leak_negp", "pdt_min", "pdt_max", "pdt_avg"
            else
                open (17, file=trim(file_path)//'quick.dat', status="old", &
                    position="append", action="write")
            endif
            write(17, "(I6.6,8E13.6)") iframe, var_global(1:5), &
                pdt_min_g, pdt_max_g, var_global(6)
            close(17)
        endif
    end subroutine quick_check

    !---------------------------------------------------------------------------
    !< Initialize the particle distributions
    !---------------------------------------------------------------------------
    subroutine init_particle_distributions
        use mhd_config_module, only: mhd_config
        implicit none
        integer :: i
        allocate(fbands(nx, ny, nz, nbands))
        allocate(fp_global(npp))
        allocate(fdpdt(npp, 2))
        if (mpi_cross_rank == master) then
            allocate(fbands_sum(nx, ny, nz, nbands))
            allocate(fp_global_sum(npp))
            allocate(fdpdt_sum(npp, 2))
        endif
        call clean_particle_distributions

        !< Intervals for distributions
        nx_mhd_reduced = (mhd_config%nx + nreduce - 1) / nreduce
        ny_mhd_reduced = (mhd_config%ny + nreduce - 1) / nreduce
        nz_mhd_reduced = (mhd_config%nz + nreduce - 1) / nreduce
        dx_diag = mhd_config%lx / nx_mhd_reduced
        dy_diag = mhd_config%ly / ny_mhd_reduced
        dz_diag = mhd_config%lz / nz_mhd_reduced
        pmin_log = log10(pmin)
        pmax_log = log10(pmax)
        dp_log = (pmax_log - pmin_log) / (npp - 1)

        !< Momentum bins for 1D distributions
        allocate(parray(npp))
        parray = 0.0_dp
        do i = 1, npp
            parray(i) = 10**(pmin_log + (i-1) * dp_log)
        enddo

        !< Momentum band edges for distributions in different energy band
        allocate(parray_bands(nbands+1))
        parray_bands = 0.0_dp
        dp_bands_log = (pmax_log - pmin_log) / nbands
        do i = 1, nbands+1
            parray_bands(i) = 10**(pmin_log + (i-1) * dp_bands_log)
        enddo

        if (mpi_rank == master) then
            write(*, "(A)") "Finished Initializing particle distributions."
        endif
    end subroutine init_particle_distributions

    !---------------------------------------------------------------------------
    !< Initialize local particle distributions
    !---------------------------------------------------------------------------
    subroutine init_local_particle_distributions
        implicit none
        allocate(fp_local(npp, nx, ny, nz))
        if (mpi_cross_rank == master) then
            allocate(fp_local_sum(npp, nx, ny, nz))
        endif
        call clean_local_particle_distribution
    end subroutine init_local_particle_distributions

    !---------------------------------------------------------------------------
    !< Set particle distributions to be zero
    !---------------------------------------------------------------------------
    subroutine clean_particle_distributions
        implicit none
        fbands = 0.0_dp
        fp_global = 0.0_dp
        fdpdt = 0.0_dp
        if (mpi_cross_rank == master) then
            fbands_sum = 0.0_dp
            fp_global_sum = 0.0_dp
            fdpdt_sum = 0.0_dp
        endif
    end subroutine clean_particle_distributions

    !---------------------------------------------------------------------------
    !< Set local particle distributions to be zero
    !---------------------------------------------------------------------------
    subroutine clean_local_particle_distribution
        implicit none
        fp_local = 0.0_dp
        if (mpi_cross_rank == master) then
            fp_local_sum = 0.0_dp
        endif
    end subroutine clean_local_particle_distribution

    !---------------------------------------------------------------------------
    !< Free particle distributions
    !---------------------------------------------------------------------------
    subroutine free_particle_distributions
        implicit none
        deallocate(fbands, parray_bands)
        deallocate(fp_global, parray)
        deallocate(fdpdt)
        if (mpi_cross_rank == master) then
            deallocate(fbands_sum)
            deallocate(fp_global_sum)
            deallocate(fdpdt_sum)
        endif
    end subroutine free_particle_distributions

    !---------------------------------------------------------------------------
    !< Free local particle distributions
    !---------------------------------------------------------------------------
    subroutine free_local_particle_distributions
        implicit none
        deallocate(fp_local)
        if (mpi_cross_rank == master) then
            deallocate(fp_local_sum)
        endif
    end subroutine free_local_particle_distributions

    !---------------------------------------------------------------------------
    !< Accumulate particle energization distributions
    !< Args:
    !<  ptl: particle structure
    !<  t0: initial time of current time interval
    !---------------------------------------------------------------------------
    subroutine energization_dist(ptl, t0)
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: interp_fields
        implicit none
        type(particle_type), intent(in) :: ptl
        real(dp), intent(in) :: t0
        real(dp) :: weight, px, py, pz, rt
        real(dp) :: dvx_dx, dvy_dy, dvz_dz
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights
        real(dp), dimension(nfields) :: fields !< Fields at particle position
        real(dp), dimension(ngrads) :: gradf !< Field gradients at particle position
        integer :: ip

        if (ptl%p > pmin .and. ptl%p <= pmax) then
            ip = ceiling((log10(ptl%p)-pmin_log) / dp_log)
            px = (ptl%x-fconfig%xmin) / mhd_config%dx
            py = (ptl%y-fconfig%ymin) / mhd_config%dy
            pz = (ptl%z-fconfig%zmin) / mhd_config%dz
            if (spherical_coord_flag) then
                call get_interp_paramters_spherical(ptl%x, ptl%y, ptl%z, pos, weights)
            else
                call get_interp_paramters(px, py, pz, pos, weights)
            endif
            rt = (ptl%t - t0) / mhd_config%dt_out
            call interp_fields(pos, weights, rt, fields, gradf)
            dvx_dx = gradf(1)
            dvy_dy = gradf(5)
            dvz_dz = gradf(9)

            fdpdt(ip, 1) = fdpdt(ip, 1) + ptl%weight
            fdpdt(ip, 2) = fdpdt(ip, 2) - &
                ptl%p * (dvx_dx + dvy_dy + dvz_dz) / 3.0d0 * ptl%weight
        endif
    end subroutine energization_dist

    !---------------------------------------------------------------------------
    !< Accumulate particle distributions
    !< Args:
    !<  local_dist: whether to accumulate local particle distribution
    !---------------------------------------------------------------------------
    subroutine calc_particle_distributions(local_dist)
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config
        implicit none
        integer, intent(in) :: local_dist
        integer :: i, ix, iy, iz, ip
        real(dp) :: weight, p, xmin, xmax, ymin, ymax, zmin, zmax
        real(dp) :: px, py, pz, rx, ry, rz, rt
        real(dp) :: dxm, dym, dzm
        type(particle_type) :: ptl

        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax
        zmin = fconfig%zmin
        zmax = fconfig%zmax
        dxm = mhd_config%dx
        dym = mhd_config%dy
        dzm = mhd_config%dz

        do i = 1, nptl_current
            ptl = ptls(i)
            ix = floor((ptl%x - xmin)/dx_diag) + 1
            iy = floor((ptl%y - ymin)/dy_diag) + 1
            iz = floor((ptl%z - zmin)/dz_diag) + 1
            p = ptl%p
            weight = ptl%weight

            !< Different momentum band
            if (ix >= 1 .and. ix <= nx .and. &
                iy >= 1 .and. iy <= ny .and. &
                iz >= 1 .and. iz <= nz) then
                ip = floor((log10(ptl%p)-pmin_log) / dp_bands_log)
                if (ip > 0 .and. ip < nbands) then
                    fbands(ix, iy, iz, ip+1) = fbands(ix, iy, iz, ip+1) + weight
                endif
            endif

            if (p > pmin .and. p <= pmax) then
                ip = ceiling((log10(ptl%p)-pmin_log) / dp_log)
                fp_global(ip) = fp_global(ip) + weight
                if (local_dist == 1) then
                    if (ix >= 1 .and. ix <= nx .and. &
                        iy >= 1 .and. iy <= ny .and. &
                        iz >= 1 .and. iz <= nz) then
                        fp_local(ip, ix, iy, iz) = fp_local(ip, ix, iy, iz) + weight
                    endif
                endif
            endif
        enddo

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fp_global, fp_global_sum, npp, MPI_DOUBLE_PRECISION, &
            MPI_SUM, master, MPI_COMM_WORLD, ierr)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fdpdt, fdpdt_sum, npp*2, MPI_DOUBLE_PRECISION, &
            MPI_SUM, master, MPI_COMM_WORLD, ierr)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fbands, fbands_sum, nx*ny*nz*nbands, MPI_DOUBLE_PRECISION, &
            MPI_SUM, master, mpi_cross_comm, ierr)
        if (local_dist == 1) then
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(fp_local, fp_local_sum, npp*nx*ny*nz, &
                MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
        endif
    end subroutine calc_particle_distributions

    !---------------------------------------------------------------------------
    !< Diagnostics of the particle distributions
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !<  local_dist: whether to accumulate local particle distribution
    !---------------------------------------------------------------------------
    subroutine distributions_diagnostics(iframe, file_path, local_dist)
        use mpi_io_module, only: set_mpi_datatype_double, set_mpi_info, fileinfo, &
            open_data_mpi_io, write_data_mpi_io
        use mhd_config_module, only: mhd_config
        use constants, only: fp, dp
        implicit none
        integer, intent(in) :: iframe, local_dist
        character(*), intent(in) :: file_path
        integer :: fh, pos1
        character(len=4) :: ctime
        character(len=128) :: fname
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
        integer :: mpi_datatype
        logical :: dir_e

        inquire(file='./data/.', exist=dir_e)
        if (.not. dir_e) then
            call system('mkdir -p ./data')
        endif

        call calc_particle_distributions(local_dist)

        write (ctime,'(i4.4)') iframe
        if (mpi_rank .eq. master) then
            fh = 15
            fname = trim(file_path)//'fp-'//ctime//'_sum.dat'
            open(fh, file=trim(fname), access='stream', status='unknown', &
                 form='unformatted', action='write')
            write(fh, pos=1) parray
            pos1 = npp * sizeof(1.0_dp) + 1
            write(fh, pos=pos1) fp_global_sum
            close(fh)

            fname = trim(file_path)//'fdpdt-'//ctime//'_sum.dat'
            open(fh, file=trim(fname), access='stream', status='unknown', &
                 form='unformatted', action='write')
            write(fh, pos=1) parray
            pos1 = npp * sizeof(1.0_dp) + 1
            write(fh, pos=pos1) fdpdt_sum
            close(fh)
        endif

        if (mpi_cross_rank == master) then
            call set_mpi_info
            ! fxy for different energy band
            fh = 18
            fname = trim(file_path)//'fxy-'//ctime//'_sum.dat'
            if (mpi_sub_rank == master) then
                fh = 18
                open(fh, file=trim(fname), access='stream', status='unknown', &
                     form='unformatted', action='write')
                write(fh, pos=1) (nbands + 0.0_dp)
                pos1 = sizeof(1.0_dp) + 1
                write(fh, pos=pos1) parray_bands
                close(fh)
            endif
            disp = (nbands + 2) * sizeof(1.0_dp)
            offset = 0
            mpi_datatype = set_mpi_datatype_double(sizes_fxy, subsizes_fxy, starts_fxy)
            call open_data_mpi_io(trim(fname), MPI_MODE_APPEND+MPI_MODE_WRONLY, &
                fileinfo, mpi_sub_comm, fh)
            call write_data_mpi_io(fh, mpi_datatype, subsizes_fxy, disp, offset, fbands_sum)
            call MPI_FILE_CLOSE(fh, ierror)

            if (local_dist == 1) then
                fh = 19
                fname = trim(file_path)//'fp_local_'//ctime//'_sum.dat'
                disp = 0
                offset = 0
                mpi_datatype = set_mpi_datatype_double(sizes_fp_local, subsizes_fp_local, starts_fp_local)
                call open_data_mpi_io(trim(fname), MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                    fileinfo, mpi_sub_comm, fh)
                call write_data_mpi_io(fh, mpi_datatype, subsizes_fp_local, disp, offset, fp_local_sum)
                call MPI_FILE_CLOSE(fh, ierror)
            endif
        endif

        call clean_particle_distributions
        if (local_dist == 1) then
            call clean_local_particle_distribution
        endif
    end subroutine distributions_diagnostics


    !---------------------------------------------------------------------------
    !< Get maximum particle momentum and write to file
    !< Args:
    !<  iframe: the time frame
    !<  if_create_file: whether to create a file
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine get_pmax_global(iframe, if_create_file, file_path)
        implicit none
        integer, intent(in) :: iframe
        logical, intent(in) :: if_create_file
        character(*), intent(in) :: file_path
        integer :: i
        logical :: dir_e

        inquire(file='./data/.', exist=dir_e)
        if (.not. dir_e) then
            call system('mkdir -p ./data')
        endif
        pmax_local = 0.0_dp
        do i = 1, nptl_current
            if (ptls(i)%p > pmax_local) then
                pmax_local = ptls(i)%p
            endif
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(pmax_local, pmax_global, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, master, MPI_COMM_WORLD, ierr)
        if (mpi_rank == master) then
            if (if_create_file) then
                open (17, file=trim(file_path)//'pmax_global.dat', status='unknown')
            else
                open (17, file=trim(file_path)//'pmax_global.dat', status="old", &
                    position="append", action="write")
            endif
            write(17, "(E13.6)") pmax_global
            close(17)
        endif
    end subroutine get_pmax_global

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
    !< Set MPI/IO data sizes for distribution diagnostics
    !---------------------------------------------------------------------------
    subroutine set_mpi_io_data_sizes
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: mpi_ix, mpi_iy, mpi_iz
        implicit none
        sizes_fxy = (/ nx_mhd_reduced, ny_mhd_reduced, nz_mhd_reduced, nbands /)
        subsizes_fxy = (/ nx, ny, nz, nbands /)
        starts_fxy = (/ nx * mpi_ix, ny * mpi_iy, nz * mpi_iz, 0 /)

        sizes_fp_local = (/ npp, nx_mhd_reduced, ny_mhd_reduced, nz_mhd_reduced /)
        subsizes_fp_local = (/ npp, nx, ny, nz /)
        starts_fp_local = (/ 0, nx * mpi_ix, ny * mpi_iy, nz * mpi_iz /)
    end subroutine set_mpi_io_data_sizes

    !---------------------------------------------------------------------------
    !< Write one double element of the particle data
    !---------------------------------------------------------------------------
    subroutine write_double_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: write_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        real(dp), dimension(:), intent(in) :: fdata
        integer(hid_t) :: dset_id, filespace
        integer :: error
        call h5screate_simple_f(1, dset_dims, filespace, error)
        call h5dcreate_f(file_id, trim(dset_name), H5T_NATIVE_DOUBLE, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .true., .true.)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)
    end subroutine write_double_element

    !---------------------------------------------------------------------------
    !< Write one integer element of the particle data
    !---------------------------------------------------------------------------
    subroutine write_integer_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: write_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        integer, dimension(:), intent(in) :: fdata
        integer(hid_t) :: dset_id, filespace
        integer :: error
        call h5screate_simple_f(1, dset_dims, filespace, error)
        call h5dcreate_f(file_id, trim(dset_name), H5T_NATIVE_INTEGER, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .true., .true.)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)
    end subroutine write_integer_element

    !---------------------------------------------------------------------------
    !< dump all the particles
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine dump_particles(iframe, file_path)
        use hdf5_io, only: create_file_h5, close_file_h5
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        integer(hsize_t), dimension(1) :: dcount, doffset, dset_dims
        real(dp) :: nptl_local, nptl_global, nptl_offset
        character(len=128) :: fname
        character(len=4) :: ctime
        integer(hid_t) :: file_id
        integer :: error
        logical :: dir_e

        inquire(file='./data/.', exist=dir_e)
        if (.not. dir_e) then
            call system('mkdir -p ./data')
        endif

        CALL h5open_f(error)

        write (ctime,'(i4.4)') iframe
        fname = trim(file_path)//'particles_'//ctime//'.h5'
        call create_file_h5(fname, H5F_ACC_TRUNC_F, file_id, .true., MPI_COMM_WORLD)

        nptl_local = nptl_current
        call MPI_ALLREDUCE(nptl_local, nptl_global, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_SCAN(nptl_local, nptl_offset, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
        nptl_offset = nptl_offset - nptl_local

        dcount(1) = nptl_local
        doffset(1) = nptl_offset
        dset_dims(1) = nptl_global

        call write_ptl_element(file_id, dcount, doffset, dset_dims, "x", ptls%x)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "y", ptls%y)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "z", ptls%z)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "p", ptls%p)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "weight", ptls%weight)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "t", ptls%t)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "dt", ptls%dt)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, &
            "split_times", ptls%split_times)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, &
            "count_flag", ptls%count_flag)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, &
            "tag", ptls%tag)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, &
            "nsteps_tracking", ptls%nsteps_tracking)

        call close_file_h5(file_id)
        CALL h5close_f(error)
    end subroutine dump_particles

    !---------------------------------------------------------------------------
    !< dump escaped particles
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine dump_escaped_particles(iframe, file_path)
        use hdf5_io, only: create_file_h5, close_file_h5
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        integer(hsize_t), dimension(1) :: dcount, doffset, dset_dims
        real(dp) :: nptl_local, nptl_global, nptl_offset
        character(len=128) :: fname
        character(len=4) :: ctime
        integer(hid_t) :: file_id
        integer :: error
        logical :: dir_e

        inquire(file='./data/.', exist=dir_e)
        if (.not. dir_e) then
            call system('mkdir -p ./data')
        endif

        nptl_local = nptl_escaped
        call MPI_ALLREDUCE(nptl_local, nptl_global, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_SCAN(nptl_local, nptl_offset, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
        nptl_offset = nptl_offset - nptl_local

        if (nptl_global > 0) then
            CALL h5open_f(error)

            write (ctime,'(i4.4)') iframe
            fname = trim(file_path)//'escaped_particles_'//ctime//'.h5'
            call create_file_h5(fname, H5F_ACC_TRUNC_F, file_id, .true., MPI_COMM_WORLD)

            dcount(1) = nptl_local
            doffset(1) = nptl_offset
            dset_dims(1) = nptl_global

            call write_ptl_element(file_id, dcount, doffset, dset_dims, "x", escaped_ptls%x)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "y", escaped_ptls%y)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "z", escaped_ptls%z)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "p", escaped_ptls%p)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "weight", escaped_ptls%weight)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "t", escaped_ptls%t)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "dt", escaped_ptls%dt)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, &
                "split_times", escaped_ptls%split_times)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, &
                "count_flag", escaped_ptls%count_flag)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, &
                "tag", escaped_ptls%tag)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, &
                "nsteps_tracking", escaped_ptls%nsteps_tracking)

            call close_file_h5(file_id)
            CALL h5close_f(error)
        else
            if (mpi_rank == master) then
                write(*, "(A)") "No escaped particles during this time interval"
            endif
        endif
    end subroutine dump_escaped_particles

end module particle_module

!*******************************************************************************
!< Module of particle data and methods to inject, remove and push particles
!*******************************************************************************
module particle_module
    use constants, only: fp, dp
    implicit none
    private
    save
    public init_particles, free_particles, inject_particles_spatial_uniform, &
        read_particle_params, particle_mover, remove_particles, split_particle, &
        init_particle_distributions, clean_particle_distributions, &
        free_particle_distributions, distributions_diagnostics, quick_check, &
        set_particle_datatype_mpi, free_particle_datatype_mpi, &
        select_particles_tracking, init_particle_tracking, &
        free_particle_tracking, init_tracked_particle_points, &
        free_tracked_particle_points, negative_particle_tags, &
        save_tracked_particle_points, inject_particles_at_shock, &
        set_mpi_io_data_sizes, init_local_particle_distributions, &
        free_local_particle_distributions, inject_particles_at_large_jz, &
        set_dpp_params, set_flags_params

    type particle_type
        real(dp) :: x, y, p         !< Position and momentum
        real(dp) :: weight, t, dt   !< Particle weight, time and time step
        integer  :: split_times     !< Particle splitting times
        integer  :: count_flag      !< Only count particle when it is 1
        integer  :: tag             !< Particle tag
        integer  :: nsteps_tracking !< Total particle tracking steps
    end type particle_type

    integer :: particle_datatype_mpi
    type(particle_type), allocatable, dimension(:) :: ptls
    type(particle_type), allocatable, dimension(:, :) :: senders
    type(particle_type), allocatable, dimension(:, :) :: recvers
    integer, dimension(4) :: nsenders, nrecvers
    !dir$ attributes align:64 :: ptls
    type(particle_type) :: ptl, ptl1

    ! Particle tracking
    integer, allocatable, dimension(:) :: nsteps_tracked_ptls, noffsets_tracked_ptls
    integer, allocatable, dimension(:) :: tags_selected_ptls
    type(particle_type), allocatable, dimension(:) :: ptl_traj_points
    integer :: nsteps_tracked_tot   !< Total tracking steps for all tracked particles

    integer :: nptl_current         !< Number of particles currently in the box
    integer :: nptl_old             !< Number of particles without receivers
    integer :: nptl_max             !< Maximum number of particles allowed
    integer :: nptl_new             !< Number of particles from splitting
    real(dp) :: leak                !< Leaking particles from boundary considering weight
    real(dp) :: leak_negp           !< Leaking particles with negative momentum

    real(dp) :: kpara0              !< Normalization for kappa parallel
    real(dp) :: kret                !< The ratio of kpara to kperp
    integer :: momentum_dependency  !< kappa dependency on particle momentum
    integer :: mag_dependency       !< kappa dependency on magnetic field
    real(dp) :: pindex              !< power index for the momentum dependency
    real(dp) :: p0    !< the standard deviation of the Gaussian distribution of momentum
    real(dp) :: b0    !< Initial magnetic field strength
    real(dp) :: kpara, kperp, dkxx_dx, dkyy_dy, dkxy_dx, dkxy_dy
    real(dp) :: skperp, skpara_perp

    real(dp) :: dt_min      !< Minimum time step
    real(dp) :: dt_max      !< Maximum time step
    real(dp) :: dt_min_rel  !< Minimum time step w.r.t. one field time interval
    real(dp) :: dt_max_rel  !< Maximum time step w.r.t. one field time interval

    !< Parameters for particle distributions
    real(dp) :: pmin  !< Minimum particle momentum
    real(dp) :: pmax  !< Maximum particle momentum
    real(dp) :: dx_diag, dy_diag, pmin_log, pmax_log, dp_log
    integer :: nx, ny, npp, nreduce
    integer :: nx_mhd_reduced, ny_mhd_reduced

    real(dp), allocatable, dimension(:, :) :: f0, f1, f2, f3, f4, f5
    real(dp), allocatable, dimension(:, :) :: f0_sum, f1_sum, f2_sum, f3_sum, f4_sum, f5_sum
    real(dp), allocatable, dimension(:) :: fp0, fp0_sum
    real(dp), allocatable, dimension(:) :: fdpdt, fdpdt_sum
    real(dp), allocatable, dimension(:, :) :: fp1, fp1_sum
    real(dp), allocatable, dimension(:, :) :: fp2, fp2_sum
    real(dp), allocatable, dimension(:) :: parray
    real(dp), allocatable, dimension(:) :: fnptl, fnptl_sum
    ! Particle distribution in each cell
    real(dp), allocatable, dimension(:, :, :) :: fp_local, fp_local_sum

    !< MPI/IO data sizes
    integer, dimension(2) :: sizes_fpx, subsizes_fpx, starts_fpx
    integer, dimension(2) :: sizes_fpy, subsizes_fpy, starts_fpy
    integer, dimension(2) :: sizes_fxy, subsizes_fxy, starts_fxy
    integer, dimension(3) :: sizes_fp_local, subsizes_fp_local, starts_fp_local

    !< Momentum diffusion
    logical :: dpp_wave_flag, dpp_shear_flag
    real(dp) :: dpp0_wave, dpp0_shear

    ! Other flags and parameters
    logical :: deltab_flag, correlation_flag
    real(dp) :: lc0 ! Normalization for turbulence correlation length

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
        ptls%p = 0.0
        ptls%weight = 0.0
        ptls%t = 0.0
        ptls%dt = 0.0
        ptls%split_times = 0
        ptls%count_flag = 0
        ptls%tag = 0
        ptls%nsteps_tracking = 0
        nptl_current = 0     ! No particle initially

        !< Particles crossing domain boundaries
        allocate(senders(nptl_max / 10, 4))
        allocate(recvers(nptl_max / 10, 4))
        senders%x = 0.0
        senders%y = 0.0
        senders%p = 0.0
        senders%weight = 0.0
        senders%t = 0.0
        senders%dt = 0.0
        senders%split_times = 0
        senders%count_flag = 0
        senders%tag = 0
        senders%nsteps_tracking = 0

        recvers%x = 0.0
        recvers%y = 0.0
        recvers%p = 0.0
        recvers%weight = 0.0
        recvers%t = 0.0
        recvers%dt = 0.0
        recvers%split_times = 0
        recvers%count_flag = 0
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
    end subroutine free_particles

    !---------------------------------------------------------------------------
    !< Set parameters for momentum diffusion
    !< Args:
    !<  dpp_wave_int(integer): flag for momentum diffusion due to wave scattering
    !<  dpp_shear_int(integer): flag for momentum diffusion due to flow shear
    !<  dpp1: normalization for Dpp due to wave scattering (v_A^2/kappa_0)
    !<  dpp2: normalization for Dpp due to flow shear (\tau_0p_0^2)
    !---------------------------------------------------------------------------
    subroutine set_dpp_params(dpp_wave_int, dpp_shear_int, dpp1, dpp2)
        implicit none
        integer, intent(in) :: dpp_wave_int, dpp_shear_int
        real(dp), intent(in) :: dpp1, dpp2
        dpp_wave_flag = .false.
        dpp_shear_flag = .false.
        if (dpp_wave_int) dpp_wave_flag = .true.
        if (dpp_shear_int) dpp_shear_flag = .true.
        dpp0_wave = dpp1
        dpp0_shear = dpp2
    end subroutine set_dpp_params

    !---------------------------------------------------------------------------
    !< Set other flags and parameters
    !< Args:
    !<  deltab_flag_int(integer): flag for magnetic fluctuation
    !<  correlation_flag_int(integer): flag for turbulence correlation length
    !<  lc0_in(real8): normalization for turbulence correlation length
    !---------------------------------------------------------------------------
    subroutine set_flags_params(deltab_flag_int, correlation_flag_int, lc0_in)
        implicit none
        integer, intent(in) :: deltab_flag_int, correlation_flag_int
        real(dp), intent(in) :: lc0_in
        deltab_flag = .false.
        correlation_flag = .false.
        if (deltab_flag_int) deltab_flag = .true.
        if (correlation_flag_int) correlation_flag = .true.
        lc0 = lc0_in
    end subroutine set_flags_params

    !---------------------------------------------------------------------------
    !< Set MPI datatype for particle type
    !---------------------------------------------------------------------------
    subroutine set_particle_datatype_mpi
        use mpi_module
        implicit none
        integer :: oldtypes(0:1), blockcounts(0:1)
        integer :: offsets(0:1), extent
        ! Setup description of the 8 MPI_DOUBLE fields.
        offsets(0) = 0
        oldtypes(0) = MPI_DOUBLE_PRECISION
        blockcounts(0) = 6
        ! Setup description of the 7 MPI_INTEGER fields.
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
        use mpi_module
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
    !---------------------------------------------------------------------------
    subroutine inject_particles_spatial_uniform(nptl, dt, dist_flag, ct_mhd)
        use simulation_setup_module, only: fconfig
        use mpi_module, only: mpi_rank, master
        use mhd_config_module, only: mhd_config
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: nptl, dist_flag, ct_mhd
        real(dp), intent(in) :: dt
        integer :: i, imod2
        real(dp) :: xmin, ymin, xmax, ymax
        real(dp) :: rands(2)

        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax

        do i = 1, nptl
            nptl_current = nptl_current + 1
            if (nptl_current > nptl_max) nptl_current = nptl_max
            ptls(nptl_current)%x = unif_01()*(xmax-xmin) + xmin
            ptls(nptl_current)%y = unif_01()*(ymax-ymin) + ymin
            if (dist_flag == 0) then
                imod2 = mod(i, 2)
                if (imod2 == 1) rands = two_normals()
                ptls(nptl_current)%p = abs(rands(imod2+1)) * p0
            else
                ptls(nptl_current)%p = p0
            endif
            ptls(nptl_current)%weight = 1.0
            ptls(nptl_current)%t = ct_mhd * mhd_config%dt_out
            ptls(nptl_current)%dt = dt
            ptls(nptl_current)%split_times = 0
            ptls(nptl_current)%count_flag = 1
            ptls(nptl_current)%tag = nptl_current
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
        use mpi_module, only: mpi_rank, master
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: interp_shock_location
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: nptl, dist_flag, ct_mhd
        real(dp), intent(in) :: dt
        integer :: i, imod2, iy, ix
        real(dp) :: xmin, ymin, xmax, ymax, dpy, shock_xpos
        real(dp), dimension(2) :: rands

        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax

        do i = 1, nptl
            nptl_current = nptl_current + 1
            if (nptl_current > nptl_max) nptl_current = nptl_max
            ptls(nptl_current)%y = unif_01()*(ymax-ymin) + ymin
            dpy = ptls(nptl_current)%y / mhd_config%dy
            iy = floor(dpy)
            !< We assume that particles are inject at the earlier shock location
            shock_xpos = interp_shock_location(iy, dpy - iy, 0.0d0) + 2  ! Two ghost cells
            ptls(nptl_current)%x = shock_xpos * (xmax - xmin) / fconfig%nxg
            if (dist_flag == 0) then
                imod2 = mod(i, 2)
                if (imod2 == 1) rands = two_normals()
                ptls(nptl_current)%p = abs(rands(imod2+1)) * p0
            else
                ptls(nptl_current)%p = p0
            endif
            ptls(nptl_current)%weight = 1.0
            ptls(nptl_current)%t = ct_mhd * mhd_config%dt_out
            ptls(nptl_current)%dt = dt
            ptls(nptl_current)%split_times = 0
            ptls(nptl_current)%count_flag = 1
            ptls(nptl_current)%tag = nptl_current
        enddo
        if (mpi_rank == master) then
            write(*, "(A)") "Finished injecting particles at the shock"
        endif
    end subroutine inject_particles_at_shock

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
        use mpi_module, only: mpi_rank, master
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: gradf, interp_fields
        use random_number_generator, only: unif_01, two_normals
        implicit none
        integer, intent(in) :: nptl, dist_flag, ct_mhd
        real(dp), intent(in) :: dt, jz_min
        real(dp), intent(in), dimension(4) :: part_box
        integer :: i, imod2, iy, ix
        real(dp) :: xmin, ymin, xmax, ymax, dpy, shock_xpos
        real(dp) :: xmin_box, ymin_box
        real(dp) :: xtmp, ytmp, px, py, rx, ry, rt, jz
        real(dp) :: dxm, dym, dby_dx, dbx_dy
        real(dp), dimension(2) :: rands

        xmin = part_box(1)
        xmax = part_box(3)
        ymin = part_box(2)
        ymax = part_box(4)
        xmin_box = fconfig%xmin
        ymin_box = fconfig%ymin
        dxm = mhd_config%dx
        dym = mhd_config%dy

        do i = 1, nptl
            nptl_current = nptl_current + 1
            if (nptl_current > nptl_max) nptl_current = nptl_max
            jz = 0.0_dp
            do while (jz < jz_min)
                xtmp = unif_01()*(xmax-xmin) + xmin
                ytmp = unif_01()*(ymax-ymin) + ymin
                ! Field interpolation parameters
                px = (xtmp-xmin_box) / dxm
                py = (ytmp-ymin_box) / dym
                ix = floor(px) + 1
                iy = floor(py) + 1
                rx = px + 1 - ix
                ry = py + 1 - iy
                rt = 0.0_dp
                call interp_fields(ix, iy, rx, ry, rt)
                dbx_dy = gradf(11)
                dby_dx = gradf(13)
                jz = abs(dby_dx - dbx_dy)
            enddo
            ptls(nptl_current)%x = xtmp
            ptls(nptl_current)%y = ytmp
            if (dist_flag == 0) then
                imod2 = mod(i, 2)
                if (imod2 == 1) rands = two_normals()
                ptls(nptl_current)%p = abs(rands(imod2+1)) * p0
            else
                ptls(nptl_current)%p = p0
            endif
            ptls(nptl_current)%weight = 1.0
            ptls(nptl_current)%t = ct_mhd * mhd_config%dt_out
            ptls(nptl_current)%dt = dt
            ptls(nptl_current)%split_times = 0
            ptls(nptl_current)%count_flag = 1
            ptls(nptl_current)%tag = nptl_current
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
        use mpi_module
        implicit none
        real(dp), intent(in) :: t0
        integer, intent(in) :: track_particle_flag, nptl_selected, nsteps_interval
        integer, intent(in) :: num_fine_steps
        real(dp) :: dtf, dxm, dym, xmin, xmax, ymin, ymax
        real(dp) :: xmin1, xmax1, ymin1, ymax1
        real(dp) :: px, py, rx, ry, rt, rt1
        real(dp) :: deltax, deltay, deltap
        real(dp) :: dt_target, dt_fine
        integer :: i, ix, iy, tracking_step, offset, tfine, step

        dtf = mhd_config%dt_out
        dt_fine = dtf / num_fine_steps
        dxm = mhd_config%dx
        dym = mhd_config%dy
        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax
        xmin1 = xmin - dxm * 0.5
        xmax1 = xmax + dxm * 0.5
        ymin1 = ymin - dym * 0.5
        ymax1 = ymax + dym * 0.5

        do i = nptl_old + 1, nptl_current
            ptl = ptls(i)
            deltax = 0.0
            deltay = 0.0
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
            if (ptl%p < 0.0 .and. ptl%count_flag /= 0) then
                ptl%count_flag = 0
                leak_negp = leak_negp + ptl%weight
            else
                call particle_boundary_condition(ptl, xmin1, xmax1, ymin1, ymax1)
            endif
            if (ptl%count_flag == 0) then
                ptls(i) = ptl
                cycle
            endif

            ! Loop over fine time steps between each MHD time interval
            ! The loop might be terminated early when dt_target is reached.
            do while (dt_target < (dtf + dt_fine*0.1)) ! 0.1 is for safe comparison
                if (ptl%count_flag == 0) then
                    exit
                endif
                do while ((ptl%t - t0) < dt_target .and. ptl%count_flag /= 0)
                    ! Check if particles are leaked or out of the local domain
                    if (ptl%p < 0.0) then
                        ptl%count_flag = 0
                        leak_negp = leak_negp + ptl%weight
                    else
                        call particle_boundary_condition(ptl, xmin1, xmax1, ymin1, ymax1)
                    endif
                    if (ptl%count_flag == 0) then
                        exit
                    endif

                    ! Field interpolation parameters
                    px = (ptl%x-xmin) / dxm
                    py = (ptl%y-ymin) / dym
                    ix = floor(px) + 1
                    iy = floor(py) + 1
                    rx = px + 1 - ix
                    ry = py + 1 - iy
                    rt = (ptl%t - t0) / dtf
                    rt1 = 1.0_dp - rt

                    call interp_fields(ix, iy, rx, ry, rt)
                    if (deltab_flag) then
                        call interp_magnetic_fluctuation(ix, iy, rx, ry, rt)
                    endif
                    if (correlation_flag) then
                        call interp_correlation_length(ix, iy, rx, ry, rt)
                    endif
                    call calc_spatial_diffusion_coefficients
                    call set_time_step(t0, dt_target)
                    call push_particle(rt, deltax, deltay, deltap)

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
                if ((ptl%t - t0) > dt_target .and. ptl%count_flag /= 0) then
                    ptl%x = ptl%x - deltax
                    ptl%y = ptl%y - deltay
                    ptl%p = ptl%p - deltap
                    ptl%t = ptl%t - ptl%dt
                    ptl%dt = t0 + dt_target - ptl%t
                    if (ptl%dt > 0) then
                        ptl%nsteps_tracking = ptl%nsteps_tracking - 1

                        px = (ptl%x-xmin) / dxm
                        py = (ptl%y-ymin) / dym
                        ix = floor(px) + 1
                        iy = floor(py) + 1

                        rx = px + 1 - ix
                        ry = py + 1 - iy
                        rt = (ptl%t - t0) / dtf
                        rt1 = 1.0_dp - rt

                        call interp_fields(ix, iy, rx, ry, rt)
                        if (deltab_flag) then
                            call interp_magnetic_fluctuation(ix, iy, rx, ry, rt)
                        endif
                        if (correlation_flag) then
                            call interp_correlation_length(ix, iy, rx, ry, rt)
                        endif
                        call calc_spatial_diffusion_coefficients
                        call push_particle(rt, deltax, deltay, deltap)
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
                        ptl%count_flag = 0
                        leak_negp = leak_negp + ptl%weight
                    else
                        call particle_boundary_condition(ptl, xmin1, xmax1, ymin1, ymax1)
                    endif
                endif

                if (ptl%count_flag /= 0) then
                    call energization_dist(t0)
                endif
                dt_target = dt_target + dt_fine

            enddo ! Fine time step loop
            ptls(i) = ptl
        enddo ! Loop over particles
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
        use mpi_module
        implicit none
        integer, intent(in) :: track_particle_flag, nptl_selected, nsteps_interval
        integer, intent(in) :: mhd_tframe, num_fine_steps
        integer :: i, local_flag, global_flag, ncycle
        logical :: all_particles_in_box
        real(dp) :: t0, xmin, xmax, ymin, ymax
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
                MPI_SUM, MPI_COMM_WORLD, ierr)
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
        nsenders = 0
        nrecvers = 0
        do i = 1, nptl_current
            ptl = ptls(i)
            if (ptl%p < 0.0) then
                ptl%count_flag = 0
                leak_negp = leak_negp + ptl%weight
            else
                call particle_boundary_condition(ptl, xmin, xmax, ymin, ymax)
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
    !---------------------------------------------------------------------------
    subroutine particle_boundary_condition(ptl, xmin, xmax, ymin, ymax)
        use simulation_setup_module, only: neighbors, mpi_ix, mpi_iy, &
            mpi_sizex, mpi_sizey
        use mhd_config_module, only: mhd_config
        use mpi_module
        implicit none
        real(dp), intent(in) :: xmin, xmax, ymin, ymax
        type(particle_type), intent(inout) :: ptl

        if (ptl%x < xmin .and. ptl%count_flag /= 0) then
            if (neighbors(1) < 0) then
                leak = leak + ptl%weight
                ptl%count_flag = 0
            else if (neighbors(1) == mpi_rank) then
                ptl%x = ptl%x - xmin + xmax
            else if (neighbors(1) == mpi_rank + mpi_sizex - 1) then
                ptl%x = ptl%x - mhd_config%xmin + mhd_config%xmax
                nsenders(1) = nsenders(1) + 1
                senders(nsenders(1), 1) = ptl
                ptl%count_flag = 0
            else
                nsenders(1) = nsenders(1) + 1
                senders(nsenders(1), 1) = ptl
                ptl%count_flag = 0
            endif
        else if (ptl%x > xmax .and. ptl%count_flag /= 0) then
            if (neighbors(2) < 0) then
                leak = leak + ptl%weight
                ptl%count_flag = 0 !< remove particle
            else if (neighbors(2) == mpi_rank) then
                ptl%x = ptl%x - xmax + xmin
            else if (neighbors(2) == mpi_iy * mpi_sizex) then !< simulation boundary
                ptl%x = ptl%x - mhd_config%xmax + mhd_config%xmin
                nsenders(2) = nsenders(2) + 1
                senders(nsenders(2), 2) = ptl
                ptl%count_flag = 0
            else
                nsenders(2) = nsenders(2) + 1
                senders(nsenders(2), 2) = ptl
                ptl%count_flag = 0
            endif
        endif

        !< We need to make sure the count_flag is not set to 0.
        !< Otherwise, we met send the particles to two different neighbors.
        if (ptl%y < ymin .and. ptl%count_flag /= 0) then
            if (neighbors(3) < 0) then
                leak = leak + ptl%weight
                ptl%count_flag = 0
            else if (neighbors(3) == mpi_rank) then
                ptl%y = ptl%y - ymin + ymax
            else if (neighbors(3) == mpi_rank + (mpi_sizey - 1) * mpi_sizex) then
                ptl%y = ptl%y - mhd_config%ymin + mhd_config%ymax
                nsenders(3) = nsenders(3) + 1
                senders(nsenders(3), 3) = ptl
                ptl%count_flag = 0
            else
                nsenders(3) = nsenders(3) + 1
                senders(nsenders(3), 3) = ptl
                ptl%count_flag = 0
            endif
        else if (ptl%y > ymax .and. ptl%count_flag /= 0) then
            if (neighbors(4) < 0) then
                leak = leak + ptl%weight
                ptl%count_flag = 0
            else if (neighbors(4) == mpi_rank) then
                ptl%y = ptl%y - ymax + ymin
            else if (neighbors(4) == mpi_ix) then
                ptl%y = ptl%y - mhd_config%ymax + mhd_config%ymin
                nsenders(4) = nsenders(4) + 1
                senders(nsenders(4), 4) = ptl
                ptl%count_flag = 0
            else
                nsenders(4) = nsenders(4) + 1
                senders(nsenders(4), 4) = ptl
                ptl%count_flag = 0
            endif
        endif
    end subroutine particle_boundary_condition

    !---------------------------------------------------------------------------
    !< Send and receiver particles from one neighbor
    !< Args:
    !<  send_id, recv_id: sender and receiver ID
    !<  mpi_direc: MPI rank along one direction
    !---------------------------------------------------------------------------
    subroutine send_recv_particle_one_neighbor(send_id, recv_id, mpi_direc)
        use mpi_module
        use simulation_setup_module, only: neighbors
        implicit none
        integer, intent(in) :: send_id, recv_id, mpi_direc
        integer :: nsend, nrecv
        nrecv = 0
        nsend = nsenders(send_id)
        if (neighbors(send_id) /= mpi_rank .and. neighbors(send_id) >= 0) then
            call MPI_SEND(nsend, 1, MPI_INTEGER, neighbors(send_id), &
                mpi_rank, MPI_COMM_WORLD, ierr)
        endif
        if (neighbors(recv_id) /= mpi_rank .and. neighbors(recv_id) >= 0) then
            call MPI_RECV(nrecv, 1, MPI_INTEGER, &
                neighbors(recv_id), neighbors(recv_id), MPI_COMM_WORLD, status, ierr)
        endif
        nrecvers(recv_id) = nrecv
        !< This assumes MPI size along this direction is even and is larger
        !< than or equal to 4
        if (mpi_direc / 2 == 0) then
            if (nsend > 0) then
                call MPI_SEND(senders(1:nsend, send_id), nsend, &
                    particle_datatype_mpi, neighbors(send_id), mpi_rank, &
                    MPI_COMM_WORLD, ierr)
            endif
            if (nrecv > 0) then
                call MPI_RECV(recvers(1:nrecv, recv_id), nrecv, &
                    particle_datatype_mpi, neighbors(recv_id), &
                    neighbors(recv_id), MPI_COMM_WORLD, status, ierr)
            endif
        else
            if (nrecv > 0) then
                call MPI_RECV(recvers(1:nrecv, recv_id), nrecv, &
                    particle_datatype_mpi, neighbors(recv_id), &
                    neighbors(recv_id), MPI_COMM_WORLD, status, ierr)
            endif
            if (nsend > 0) then
                call MPI_SEND(senders(1:nsend, send_id), nsend, &
                    particle_datatype_mpi, neighbors(send_id), mpi_rank, &
                    MPI_COMM_WORLD, ierr)
            endif
        endif
    end subroutine send_recv_particle_one_neighbor

    !---------------------------------------------------------------------------
    !< Send particles to neighbors
    !---------------------------------------------------------------------------
    subroutine send_recv_particles
        use mpi_module
        use simulation_setup_module, only: neighbors, mpi_ix, mpi_iy
        implicit none
        call send_recv_particle_one_neighbor(1, 2, mpi_ix) !< Right -> Left
        call send_recv_particle_one_neighbor(2, 1, mpi_ix) !< Left -> Right
        call send_recv_particle_one_neighbor(3, 4, mpi_iy) !< Top -> Bottom
        call send_recv_particle_one_neighbor(4, 3, mpi_iy) !< Bottom -> Top
    end subroutine send_recv_particles

    !---------------------------------------------------------------------------
    !< Calculate the spatial diffusion coefficients
    !---------------------------------------------------------------------------
    subroutine calc_spatial_diffusion_coefficients
        use mhd_data_parallel, only: fields, gradf, db2, lc, grad_db, grad_lc
        implicit none
        real(dp) :: pnorm
        real(dp) :: bx, by, b, ib1, ib2, ib3, ib4
        real(dp) :: dbx_dx, dby_dx, dbx_dy, dby_dy, db_dx, db_dy
        real(dp) :: dkdx, dkdy

        bx = fields(5)
        by = fields(6)
        b = fields(8)
        dbx_dx = gradf(10)
        dbx_dy = gradf(11)
        dby_dx = gradf(13)
        dby_dy = gradf(14)
        db_dx = gradf(19)
        db_dy = gradf(20)
        if (b == 0) then
            ib1 = 1.0
        else
            ib1 = 1.0_dp / b
        endif
        ib2 = ib1 * ib1
        ib3 = ib1 * ib2
        ib4 = ib2 * ib2

        pnorm = 1.0_dp
        if (mag_dependency == 1) then
            pnorm = pnorm * ib1**(1./3.)
        endif
        if (momentum_dependency == 1) then
            pnorm = pnorm * (ptl%p / p0)**pindex
        endif

        ! Magnetic fluctuation dB^2/B^2
        if (deltab_flag) then
            pnorm = pnorm / db2
        endif

        ! Turbulence correlation length
        if (correlation_flag) then
            pnorm = pnorm * (lc/lc0)**(2./3.)
        endif

        kpara = kpara0 * pnorm
        kperp = kpara * kret

        ! if (mag_dependency == 1) then
        !     dkxx_dx = -kperp*db_dx*ib1 + (kperp-kpara)*db_dx*bx*bx*ib3 + &
        !         2.0*(kpara-kperp)*bx*(dbx_dx*b-bx*db_dx)*ib3
        !     dkyy_dy = -kperp*db_dy*ib1 + (kperp-kpara)*db_dy*by*by*ib3 + &
        !         2.0*(kpara-kperp)*by*(dby_dy*b-by*db_dy)*ib3
        !     dkxy_dx = (kperp-kpara)*db_dx*bx*by*ib3 + (kpara-kperp) * &
        !         ((dbx_dx*by+bx*dby_dx)*ib2 - 2.0*bx*by*db_dx*ib3)
        !     dkxy_dy = (kperp-kpara)*db_dy*bx*by*ib3 + (kpara-kperp) * &
        !         ((dbx_dy*by+bx*dby_dy)*ib2 - 2.0*bx*by*db_dy*ib3)
        ! else
        !     dkxx_dx = 2.0*(kpara-kperp)*bx*(dbx_dx*b-bx*db_dx)*ib3
        !     dkyy_dy = 2.0*(kpara-kperp)*by*(dby_dy*b-by*db_dy)*ib3
        !     dkxy_dx = (kpara-kperp) * ((dbx_dx*by+bx*dby_dx)*ib2 - 2.0*bx*by*db_dx*ib3)
        !     dkxy_dy = (kpara-kperp) * ((dbx_dy*by+bx*dby_dy)*ib2 - 2.0*bx*by*db_dy*ib3)
        ! endif
        dkdx = 2 * grad_lc(1) / (3 * lc) - grad_db(1) / db2
        dkdy = 2 * grad_lc(2) / (3 * lc) - grad_db(2) / db2
        if (mag_dependency == 1) then
            dkdx = dkdx - db_dx * ib1 / 3
            dkdy = dkdy - db_dx * ib1 / 3
        endif
        dkxx_dx = kperp*dkdx + (kpara+kperp)*dkdx*bx**2*ib2 + &
            2.0*(kpara-kperp)*bx*(dbx_dx*b-bx*db_dx)*ib3
        dkyy_dy = kperp*dkdy + (kpara-kperp)*dkdy*by**2*ib2 + &
            2.0*(kpara-kperp)*by*(dby_dy*b-by*db_dy)*ib3
        dkxy_dx = (kpara-kperp)*dkdx*bx*by*ib2 + (kpara-kperp) * &
            ((dbx_dx*by+bx*dby_dx)*ib2 - 2.0*bx*by*db_dx*ib3)
        dkxy_dy = (kpara-kperp)*dkdy*bx*by*ib2 + (kpara-kperp) * &
            ((dbx_dy*by+bx*dby_dy)*ib2 - 2.0*bx*by*db_dy*ib3)
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
        use mpi_module
        implicit none
        character(*), intent(in) :: conf_file
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
            nx = fconfig%nx / nreduce
            ny = fconfig%ny / nreduce
            temp = get_variable(fh, 'npp', '=')
            npp = int(temp)
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
            write(*, "(A,I0,A,I0)") " Dimensions of spatial distributions = ", &
                nx, " ", ny
            write(*, "(A,I0)") " Dimensions of momentum distributions = ", npp
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
        call MPI_BCAST(npp, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)

        if (nx * nreduce /= fconfig%nx .or. ny * nreduce /= fconfig%ny) then
            if (mpi_rank == master) then
                write(*, "(A)") "Wrong factor 'nreduce' for particle distribution"
            endif
            call MPI_FINALIZE(ierr)
            stop
        endif
    end subroutine read_particle_params

    !---------------------------------------------------------------------------
    !< Determine the time step.
    !< Args;
    !<  t0: the initial time for current particle
    !<  dtf: the time interval between fine diagnostics
    !---------------------------------------------------------------------------
    subroutine set_time_step(t0, dtf)
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: fields
        implicit none
        real(dp), intent(in) :: t0, dtf
        real(dp) :: tmp30, tmp40, bx, by, b, vx, vy, dxm, dym, dt1, dt2
        skperp = dsqrt(2.0*kperp)
        skpara_perp = dsqrt(2.0*(kpara-kperp))

        vx = fields(1)
        vy = fields(2)
        bx = fields(5)
        by = fields(6)
        b = fields(8)
        dxm = mhd_config%dx
        dym = mhd_config%dy

        if (b == 0.0d0) then
            tmp30 = 0.0
        else
            tmp30 = skperp + skpara_perp * abs(bx/b)
        endif
        tmp40 = abs(vx + dkxx_dx + dkxy_dy)
        if (tmp40 .ne. 0.0d0) then
            if (tmp30 > 0) then
                dt1 = min(dxm/(80.0*tmp40), (tmp30/tmp40)**2) * 0.5d0
                ! dt1 = min(dxm/tmp40, (dxm/tmp30)**2) * 0.2
                ! dt1 = min(dt1, (tmp30/tmp40)**2) * 0.05
            else
                dt1 = dxm/(80.0*tmp40)
            endif
        else
            dt1 = dt_min
        endif
        if (b == 0.0d0) then
            tmp30 = 0.0
        else
            tmp30 = skperp + skpara_perp * abs(by/b)
        endif
        tmp40 = abs(vy + dkxy_dx + dkyy_dy)
        if (tmp40 .ne. 0.0d0) then
            if (tmp30 > 0) then
                dt2 = min(dym/(80.0*tmp40), (tmp30/tmp40)**2) * 0.5d0
                ! dt2 = min(dym/tmp40, (dym/tmp30)**2) * 0.2
                ! dt2 = min(dt2, (tmp30/tmp40)**2) * 0.05
            else
                dt2 = dym/(80.0*tmp40)
            endif
        else
            dt2 = dt_min
        endif
        ptl%dt = min(dt1, dt2)

        !< Make sure the time step is not too large. Adding dt_min to make
        !< sure to exit the while where this routine is called
        if ((ptl%t + ptl%dt - t0) > dtf) then
            ptl%dt = t0 + dtf - ptl%t + dt_min
        endif

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
    !< Push particle for a single step
    !< Args:
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  deltax, deltay, deltap: the change of x, y and p in this step
    !---------------------------------------------------------------------------
    subroutine push_particle(rt, deltax, deltay, deltap)
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: fconfig
        use mhd_data_parallel, only: fields, gradf, interp_fields, &
            interp_magnetic_fluctuation, interp_correlation_length
        use random_number_generator, only: unif_01, two_normals
        implicit none
        real(dp), intent(in) :: rt
        real(dp), intent(out) :: deltax, deltay, deltap
        real(dp) :: xtmp, ytmp
        real(dp) :: sdt, dvx_dx, dvx_dy, dvy_dx, dvy_dy, divv, gshear
        real(dp) :: bx, by, b, vx, vy, px, py, rx, ry, rt1, ib
        real(dp) :: bx1, by1, b1, ib1
        real(dp) :: xmin, ymin, xmax, ymax, dxm, dym, skperp1, skpara_perp1
        reaL(dp) :: xmin1, ymin1, xmax1, ymax1, dxmh, dymh
        real(dp) :: ran1, ran2, ran3, sqrt3
        real(dp) :: rho, va ! Plasma density and Alfven speed
        real(dp) :: rands(2)
        integer :: ix, iy

        vx = fields(1)
        vy = fields(2)
        bx = fields(5)
        by = fields(6)
        b = fields(8)
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
        xtmp = ptl%x + (vx+dkxx_dx+dkxy_dy)*ptl%dt + (skperp+skpara_perp*bx*ib)*sdt
        ytmp = ptl%y + (vy+dkxy_dx+dkyy_dy)*ptl%dt + (skperp+skpara_perp*by*ib)*sdt
        deltax = (vx+dkxx_dx+dkxy_dy)*ptl%dt
        deltay = (vy+dkxy_dx+dkyy_dy)*ptl%dt
        sqrt3 = dsqrt(3.0_dp)
        ran1 = (2.0_dp*unif_01() - 1.0_dp) * sqrt3
        ran2 = (2.0_dp*unif_01() - 1.0_dp) * sqrt3
        ran3 = (2.0_dp*unif_01() - 1.0_dp) * sqrt3

        ! rands = two_normals()
        ! ran1 = rands(1)
        ! ran2 = rands(2)
        ! rands = two_normals()
        ! ran3 = rands(1)

        !< We originally tried to decrease the time step when xtmp or ytmp are out-of-bound,
        !< but ecreasing the time step does not necessarily make the moving distance smaller.
        !< Therefore, we switch between first-order and second-order method.
        if (xtmp < xmin1 .or. xtmp > xmax1 .or. ytmp < ymin1 .or. ytmp > ymax1) then
            !< First-order method
            deltax = deltax + ran1*skperp*sdt + ran3*skpara_perp*sdt*bx*ib
            deltay = deltay + ran2*skperp*sdt + ran3*skpara_perp*sdt*by*ib
        else
            !< Second-order method. It requires xtmp and ytmp are in the local domain.
            px = (xtmp - xmin) / dxm
            py = (ytmp - ymin) / dym
            ix = floor(px) + 1
            iy = floor(py) + 1
            rx = px + 1 - ix
            ry = py + 1 - iy
            rt1 = 1.0_dp - rt
            call interp_fields(ix, iy, rx, ry, rt)
            if (deltab_flag) then
                call interp_magnetic_fluctuation(ix, iy, rx, ry, rt)
            endif
            if (correlation_flag) then
                call interp_correlation_length(ix, iy, rx, ry, rt)
            endif

            call calc_spatial_diffusion_coefficients

            !< Magnetic field at the predicted position
            bx1 = fields(5)
            by1 = fields(6)
            b1 = fields(8)

            skperp1 = dsqrt(2.0*kperp)
            skpara_perp1 = dsqrt(2.0*(kpara-kperp))

            if (b1 == 0) then
                ib1 = 0.0
            else
                ib1 = 1.0 / b1
            endif

            deltax = deltax + ran1*skperp*sdt + ran3*skpara_perp*sdt*bx*ib + &
                     (skperp1-skperp)*(ran1*ran1-1.0)*sdt/2.0 + &
                     (skpara_perp1*bx1*ib1-skpara_perp*bx*ib)*(ran3*ran3-1.0)*sdt/2.0
            deltay = deltay + ran2*skperp*sdt + ran3*skpara_perp*sdt*by*ib + &
                     (skperp1-skperp)*(ran2*ran2-1.0)*sdt/2.0 + &
                     (skpara_perp1*by1*ib1-skpara_perp*by*ib)*(ran3*ran3-1.0)*sdt/2.0
        endif

        ptl%x = ptl%x + deltax
        ptl%y = ptl%y + deltay
        ptl%t = ptl%t + ptl%dt

        ! Momentum
        divv = dvx_dx + dvy_dy
        deltap = -ptl%p * divv * ptl%dt / 3.0d0
        ! Momentum diffusion due to wave scattering
        if (dpp_wave_flag) then
            rho = fields(4)
            va = b / dsqrt(rho)
            if (momentum_dependency) then
                deltap = deltap + (8*ptl%p / (27*kpara)) * dpp0_wave * va**2 * ptl%dt
            else
                deltap = deltap + (4*ptl%p / (9*kpara)) * dpp0_wave * va**2 * ptl%dt
            endif
            ran1 = (2.0_dp*unif_01() - 1.0_dp) * sqrt3
            deltap = deltap + ran1 * va * ptl%p * dsqrt(2*dpp0_wave/(9*kpara)) * sdt
        endif

        ! Momentum diffusion due to flow shear
        if (dpp_shear_flag) then
            dvx_dy = gradf(2)
            dvy_dx = gradf(4)
            gshear = (2*(dvx_dy+dvy_dx)**2 + 4*(dvx_dx**2+dvy_dy**2))/30 - 2*divv**2/45
            deltap = deltap + 2 * gshear * dpp0_shear * ptl%dt / ptl%p
            ran1 = (2.0_dp*unif_01() - 1.0_dp) * sqrt3
            deltap = deltap + dsqrt(2 * gshear * dpp0_shear) * ran1 * sdt
        endif

        ptl%p = ptl%p + deltap
    end subroutine push_particle

    !---------------------------------------------------------------------------
    !< Remove particles from simulation if their count_flags are 0.
    !---------------------------------------------------------------------------
    subroutine remove_particles
        implicit none
        integer :: i
        i = 1
        do while (i <= nptl_current)
            if (ptls(i)%count_flag == 0) then
                !< Switch current with the last particle
                ptl1 = ptls(nptl_current)
                ptls(nptl_current) = ptls(i)
                ptls(i) = ptl1
                nptl_current = nptl_current - 1
            else
                i = i + 1
            endif
        enddo
    end subroutine remove_particles

    !---------------------------------------------------------------------------
    !< Add particles from neighbors
    !---------------------------------------------------------------------------
    subroutine add_neighbor_particles
        implicit none
        integer :: i, j, nrecv
        nptl_old = nptl_current
        do i = 1, 4
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
        type(particle_type) :: ptl_new
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
                ptls(nptl_current)%tag = nptl_current
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
        use mpi_module
        implicit none
        integer, intent(in) :: iframe
        logical, intent(in) :: if_create_file
        character(*), intent(in) :: file_path
        real(dp) :: pdt_min, pdt_max, pdt_avg, ntot
        real(dp) :: ntot_g, leak_g, leak_negp_g
        integer :: i, nptl_current_g
        logical :: dir_e

        inquire(file='./data/.', exist=dir_e)
        if (.not. dir_e) then
            call system('mkdir -p ./data')
        endif
        ntot = 0.0_dp
        pdt_min = 1.0_dp
        pdt_max = 0.0_dp
        pdt_avg = 0.0_dp
        ntot_g = 0.0_dp
        leak_g = 0.0_dp
        leak_negp_g = 0.0_dp
        nptl_current_g = 0
        do i = 1, nptl_current
            ntot = ntot + ptls(i)%weight
            if (ptls(i)%dt < pdt_min) pdt_min = ptls(i)%dt
            if (ptls(i)%dt > pdt_max) pdt_max = ptls(i)%dt
            pdt_avg = pdt_avg + ptls(i)%dt
        enddo
        if (nptl_current > 0) then
            pdt_avg = pdt_avg / nptl_current
        else
            pdt_avg = 0.0
        endif
        call MPI_REDUCE(ntot, ntot_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            master, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(leak, leak_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            master, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(leak_negp, leak_negp_g, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
            master, MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(nptl_current, nptl_current_g, 1, MPI_INTEGER, MPI_SUM, &
            master, MPI_COMM_WORLD, ierr)
        if (mpi_rank == master) then
            if (if_create_file) then
                open (17, file=trim(file_path)//'quick.dat', status='unknown')
                write(17, "(A,A)") "iframe, nptl_current, nptl_new, ntot, ", &
                    "leak, leak_negp, pdt_min, pdt_max, pdt_avg"
            else
                open (17, file=trim(file_path)//'quick.dat', status="old", &
                    position="append", action="write")
            endif
            write(17, "(I4.4,A,I9.9,A,I6.6,A,6E13.6E2)") &
                iframe, ' ', nptl_current_g, ' ', nptl_new, ' ', ntot_g, &
                leak_g, leak_negp_g, pdt_min, pdt_max, pdt_avg
            close(17)
        endif
    end subroutine quick_check

    !---------------------------------------------------------------------------
    !< Initialize the particle distributions
    !---------------------------------------------------------------------------
    subroutine init_particle_distributions
        use mhd_config_module, only: mhd_config
        use mpi_module, only: mpi_rank, master
        implicit none
        integer :: i
        allocate(f0(nx, ny))
        allocate(f1(nx, ny))
        allocate(f2(nx, ny))
        allocate(f3(nx, ny))
        allocate(f4(nx, ny))
        allocate(f5(nx, ny))
        allocate(fp0(npp))
        allocate(fp1(npp, nx))
        allocate(fp2(npp, ny))
        allocate(fdpdt(npp))
        allocate(fnptl(npp))
        allocate(f0_sum(nx, ny))
        allocate(f1_sum(nx, ny))
        allocate(f2_sum(nx, ny))
        allocate(f3_sum(nx, ny))
        allocate(f4_sum(nx, ny))
        allocate(f5_sum(nx, ny))
        allocate(fp0_sum(npp))
        allocate(fp1_sum(npp, nx))
        allocate(fp2_sum(npp, ny))
        allocate(fdpdt_sum(npp))
        allocate(fnptl_sum(npp))
        call clean_particle_distributions

        !< Intervals for distributions
        nx_mhd_reduced = mhd_config%nx / nreduce
        ny_mhd_reduced = mhd_config%ny / nreduce
        dx_diag = mhd_config%lx / nx_mhd_reduced
        dy_diag = mhd_config%ly / ny_mhd_reduced
        pmin_log = log10(pmin)
        pmax_log = log10(pmax)
        dp_log = (pmax_log - pmin_log) / (npp - 1)

        ! if (mpi_rank == master) then
            allocate(parray(npp))
            parray = 0.0_dp
            do i = 1, npp
                parray(i) = 10**(pmin_log + (i-1) * dp_log)
            enddo
        ! endif
        if (mpi_rank == master) then
            write(*, "(A)") "Finished Initializing particle distributions."
        endif
    end subroutine init_particle_distributions

    !---------------------------------------------------------------------------
    !< Initialize local particle distributions
    !---------------------------------------------------------------------------
    subroutine init_local_particle_distributions
        use mpi_module, only: mpi_rank, master
        implicit none
        allocate(fp_local(npp, nx, ny))
        if (mpi_rank == master) then
            allocate(fp_local_sum(npp, nx, ny))
        endif
        call clean_local_particle_distribution
    end subroutine init_local_particle_distributions

    !---------------------------------------------------------------------------
    !< Set particle distributions to be zero
    !---------------------------------------------------------------------------
    subroutine clean_particle_distributions
        implicit none
        f0 = 0.0_dp
        f1 = 0.0_dp
        f2 = 0.0_dp
        f3 = 0.0_dp
        f4 = 0.0_dp
        f5 = 0.0_dp
        fp0 = 0.0_dp
        fp1 = 0.0_dp
        fp2 = 0.0_dp
        fdpdt = 0.0_dp
        fnptl = 0.0_dp
        f0_sum = 0.0_dp
        f1_sum = 0.0_dp
        f2_sum = 0.0_dp
        f3_sum = 0.0_dp
        f4_sum = 0.0_dp
        f5_sum = 0.0_dp
        fp0_sum = 0.0_dp
        fp1_sum = 0.0_dp
        fp2_sum = 0.0_dp
        fdpdt_sum = 0.0_dp
        fnptl_sum = 0.0_dp
    end subroutine clean_particle_distributions

    !---------------------------------------------------------------------------
    !< Set local particle distributions to be zero
    !---------------------------------------------------------------------------
    subroutine clean_local_particle_distribution
        use mpi_module, only: mpi_rank, master
        implicit none
        fp_local = 0.0_dp
        if (mpi_rank == master) then
            fp_local_sum = 0.0_dp
        endif
    end subroutine clean_local_particle_distribution

    !---------------------------------------------------------------------------
    !< Free particle distributions
    !---------------------------------------------------------------------------
    subroutine free_particle_distributions
        use mpi_module, only: mpi_rank, master
        implicit none
        deallocate(f0, f1, f2, f3, f4, f5, fp0, fp1, fp2)
        deallocate(f0_sum, f1_sum, f2_sum, f3_sum, f4_sum, f5_sum)
        deallocate(fp0_sum, fp1_sum, fp2_sum)
        deallocate(fdpdt, fdpdt_sum)
        deallocate(fnptl, fnptl_sum)
        ! if (mpi_rank == master) then
            deallocate(parray)
        ! endif
    end subroutine free_particle_distributions

    !---------------------------------------------------------------------------
    !< Free local particle distributions
    !---------------------------------------------------------------------------
    subroutine free_local_particle_distributions
        use mpi_module, only: mpi_rank, master
        implicit none
        deallocate(fp_local)
        if (mpi_rank == master) then
            deallocate(fp_local_sum)
        endif
    end subroutine free_local_particle_distributions

    !---------------------------------------------------------------------------
    !< Accumulate particle energization distributions
    !---------------------------------------------------------------------------
    subroutine energization_dist(t0)
        use mpi_module
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config
        use mhd_data_parallel, only: gradf, interp_fields
        implicit none
        real(dp), intent(in) :: t0
        integer :: ix, iy, ip
        real(dp) :: weight, px, py, rx, ry, rt
        real(dp) :: dvx_dx, dvy_dy

        if (ptl%p > pmin .and. ptl%p <= pmax) then
            ip = ceiling((log10(ptl%p)-pmin_log) / dp_log)
            px = (ptl%x-fconfig%xmin) / mhd_config%dx
            py = (ptl%y-fconfig%ymin) / mhd_config%dy
            rt = (ptl%t - t0) / mhd_config%dt_out

            ix = ceiling((ptl%x - fconfig%xmin)/dx_diag)
            iy = ceiling((ptl%y - fconfig%ymin)/dy_diag)
            if (ix < 1) ix = 1
            if (ix > nx) ix = nx
            if (iy < 1) iy = 1
            if (iy > ny) iy = ny
            ix = floor(px) + 1
            iy = floor(py) + 1
            rx = px + 1 - ix
            ry = py + 1 - iy
            call interp_fields(ix, iy, rx, ry, rt)
            dvx_dx = gradf(1)
            dvy_dy = gradf(5)

            fnptl(ip) = fnptl(ip) + ptl%weight
            fdpdt(ip) = fdpdt(ip) - ptl%p * (dvx_dx + dvy_dy) / 3.0d0 * ptl%weight
        endif
    end subroutine energization_dist

    !---------------------------------------------------------------------------
    !< Accumulate particle distributions
    !< Args:
    !<  whole_mhd_data: whether each MPI process holds the whole MHD data
    !<  local_dist: whether to accumulate local particle distribution
    !---------------------------------------------------------------------------
    subroutine calc_particle_distributions(whole_mhd_data, local_dist)
        use mpi_module
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config
        implicit none
        integer, intent(in) :: whole_mhd_data, local_dist
        integer :: i, ix, iy, ip
        real(dp) :: weight, p, xmin, xmax, ymin, ymax
        real(dp) :: px, py, rx, ry, rt
        real(dp) :: dxm, dym, dvx_dx, dvy_dy

        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax
        dxm = mhd_config%dx
        dym = mhd_config%dy

        do i = 1, nptl_current
            ptl = ptls(i)
            ix = ceiling((ptl%x - xmin)/dx_diag)
            iy = ceiling((ptl%y - ymin)/dy_diag)
            p = ptl%p
            weight = ptl%weight

            !< Different momentum band
            if (ix >= 1 .and. ix <= nx .and. iy >= 1 .and. iy <= ny) then
                if (p <= 0.75*p0) then
                    f0(ix,iy) = f0(ix,iy) + weight
                else if (p > 0.75*p0 .and. p <= 1.5*p0) then
                    f1(ix,iy) = f1(ix,iy) + weight
                else if (p > 1.5*p0 .and. p <= 3.0*p0) then
                    f2(ix,iy) = f2(ix,iy) + weight
                else if (p > 3.0*p0 .and. p <= 6.0*p0) then
                    f3(ix,iy) = f3(ix,iy) + weight
                else if (p > 6.0*p0 .and. p <= 12.0*p0) then
                    f4(ix,iy) = f4(ix,iy) + weight
                else
                    f5(ix,iy) = f5(ix,iy) + weight
                endif
            endif

            if (p > pmin .and. p <= pmax) then
                ip = ceiling((log10(ptl%p)-pmin_log) / dp_log)
                fp0(ip) = fp0(ip) + weight
                if (ix >= 1 .and. ix <= nx) then
                    fp1(ip, ix) = fp1(ip, ix) + weight
                endif
                if (iy >= 1 .and. iy <= ny) then
                    fp2(ip, iy) = fp2(ip, iy) + weight
                endif
                if (local_dist) then
                    if (ix >= 1 .and. ix <= nx .and. iy >= 1 .and. iy <= ny) then
                        fp_local(ip, ix, iy) = fp_local(ip, ix, iy) + weight
                    endif
                endif
            endif
        enddo

        call MPI_REDUCE(fp0, fp0_sum, npp, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fnptl, fnptl_sum, npp, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fdpdt, fdpdt_sum, npp, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
            MPI_COMM_WORLD, ierr)
        if (whole_mhd_data == 1) then
            call MPI_REDUCE(f0, f0_sum, nx*ny, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
                MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(f1, f1_sum, nx*ny, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
                MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(f2, f2_sum, nx*ny, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
                MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(f3, f3_sum, nx*ny, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
                MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(f4, f4_sum, nx*ny, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
                MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(f5, f5_sum, nx*ny, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
                MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(fp1, fp1_sum, npp*nx, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
                MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(fp2, fp2_sum, npp*ny, MPI_DOUBLE_PRECISION, MPI_SUM, master, &
                MPI_COMM_WORLD, ierr)
            if (local_dist) then
                call MPI_REDUCE(fp_local, fp_local_sum, npp*nx*ny, MPI_DOUBLE_PRECISION, &
                    MPI_SUM, master, MPI_COMM_WORLD, ierr)
            endif
        else
            if (nx == mhd_config%nx / nreduce) then
                call MPI_REDUCE(fp1, fp1_sum, npp*nx, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    master, MPI_COMM_WORLD, ierr)
            endif
            if (ny == mhd_config%ny / nreduce) then
                call MPI_REDUCE(fp2, fp2_sum, npp*ny, MPI_DOUBLE_PRECISION, MPI_SUM, &
                    master, MPI_COMM_WORLD, ierr)
            endif
        endif
    end subroutine calc_particle_distributions

    !---------------------------------------------------------------------------
    !< Diagnostics of the particle distributions
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !<  whole_mhd_data: whether each MPI process holds the whole MHD data
    !<  local_dist: whether to accumulate local particle distribution
    !---------------------------------------------------------------------------
    subroutine distributions_diagnostics(iframe, file_path, whole_mhd_data, local_dist)
        use mpi_io_module, only: set_mpi_datatype_double, set_mpi_info, fileinfo, &
            open_data_mpi_io, write_data_mpi_io
        use mhd_config_module, only: mhd_config
        use constants, only: fp, dp
        use mpi_module
        implicit none
        integer, intent(in) :: iframe, whole_mhd_data, local_dist
        character(*), intent(in) :: file_path
        integer :: fh, pos1
        character(len=4) :: ctime, mrank
        character(len=128) :: fname
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
        integer :: mpi_datatype
        logical :: dir_e

        inquire(file='./data/.', exist=dir_e)
        if (.not. dir_e) then
            call system('mkdir -p ./data')
        endif

        call calc_particle_distributions(whole_mhd_data, local_dist)

        write (ctime,'(i4.4)') iframe
        write (mrank,'(i4.4)') mpi_rank
        if (mpi_rank .eq. 0) then
            fh = 15
            fname = trim(file_path)//'fp-'//ctime//'_sum.dat'
            open(fh, file=trim(fname), access='stream', status='unknown', &
                 form='unformatted', action='write')
            write(fh, pos=1) parray
            pos1 = npp * sizeof(1.0_dp) + 1
            write(fh, pos=pos1) fp0_sum
            close(fh)

            fname = trim(file_path)//'fdpdt-'//ctime//'_sum.dat'
            open(fh, file=trim(fname), access='stream', status='unknown', &
                 form='unformatted', action='write')
            write(fh, pos=1) parray
            pos1 = npp * sizeof(1.0_dp) + 1
            write(fh, pos=pos1) fnptl_sum
            pos1 = pos1 + npp * sizeof(1.0_dp)
            write(fh, pos=pos1) fdpdt_sum
            close(fh)
        endif

        if (whole_mhd_data == 1) then
            if (mpi_rank == master) then
                fh = 16
                fname = trim(file_path)//'fpx-'//ctime//'_sum.dat'
                open(fh, file=trim(fname), access='stream', status='unknown', &
                     form='unformatted', action='write')
                write(fh, pos=1) fp1_sum
                close(fh)

                fh = 17
                fname = trim(file_path)//'fpy-'//ctime//'_sum.dat'
                open(fh, file=trim(fname), access='stream', status='unknown', &
                     form='unformatted', action='write')
                write(fh, pos=1) fp2_sum
                close(fh)

                fh = 18
                fname = trim(file_path)//'fxy-'//ctime//'_sum.dat'
                open(fh, file=trim(fname), access='stream', status='unknown', &
                     form='unformatted', action='write')
                write(fh, pos=1) f0_sum
                pos1 = nx * ny * sizeof(1.0_dp) + 1
                write(fh, pos=pos1) f1_sum
                pos1 = pos1 + nx * ny * sizeof(1.0_dp)
                write(fh, pos=pos1) f2_sum
                pos1 = pos1 + nx * ny * sizeof(1.0_dp)
                write(fh, pos=pos1) f3_sum
                pos1 = pos1 + nx * ny * sizeof(1.0_dp)
                write(fh, pos=pos1) f4_sum
                pos1 = pos1 + nx * ny * sizeof(1.0_dp)
                write(fh, pos=pos1) f5_sum
                close(fh)

                if (local_dist) then
                    fh = 19
                    fname = trim(file_path)//'fp_local_'//ctime//'_sum.dat'
                    open(fh, file=trim(fname), access='stream', status='unknown', &
                         form='unformatted', action='write')
                    write(fh, pos=1) fp_local_sum
                    close(fh)
                endif
            endif
        else
            call set_mpi_info
            ! fpx
            fh = 16
            fname = trim(file_path)//'fpx-'//ctime//'_sum.dat'
            if (nx == mhd_config%nx) then
                open(fh, file=trim(fname), access='stream', status='unknown', &
                     form='unformatted', action='write')
                write(fh, pos=1) fp1_sum
                close(fh)
            else
                mpi_datatype = set_mpi_datatype_double(sizes_fpx, subsizes_fpx, starts_fpx)
                call open_data_mpi_io(trim(fname), MPI_MODE_CREATE+MPI_MODE_WRONLY, fileinfo, fh)
                disp = 0
                offset = 0
                call write_data_mpi_io(fh, mpi_datatype, subsizes_fpx, disp, offset, fp1)
                call MPI_FILE_CLOSE(fh, ierror)
            endif

            ! fpy
            fh = 17
            fname = trim(file_path)//'fpy-'//ctime//'_sum.dat'
            if (ny == mhd_config%ny) then
                open(fh, file=trim(fname), access='stream', status='unknown', &
                     form='unformatted', action='write')
                write(fh, pos=1) fp2_sum
                close(fh)
            else
                mpi_datatype = set_mpi_datatype_double(sizes_fpy, subsizes_fpy, starts_fpy)
                call open_data_mpi_io(trim(fname), MPI_MODE_CREATE+MPI_MODE_WRONLY, fileinfo, fh)
                disp = 0
                offset = 0
                call write_data_mpi_io(fh, mpi_datatype, subsizes_fpy, disp, offset, fp2)
                call MPI_FILE_CLOSE(fh, ierror)
            endif

            ! fxy for different energy band
            fh = 18
            fname = trim(file_path)//'fxy-'//ctime//'_sum.dat'
            disp = 0
            offset = 0
            mpi_datatype = set_mpi_datatype_double(sizes_fxy, subsizes_fxy, starts_fxy)
            call open_data_mpi_io(trim(fname), MPI_MODE_CREATE+MPI_MODE_WRONLY, fileinfo, fh)
            call write_data_mpi_io(fh, mpi_datatype, subsizes_fxy, disp, offset, f0)
            disp = disp + nx_mhd_reduced * ny_mhd_reduced * 8
            call write_data_mpi_io(fh, mpi_datatype, subsizes_fxy, disp, offset, f1)
            disp = disp + nx_mhd_reduced * ny_mhd_reduced * 8
            call write_data_mpi_io(fh, mpi_datatype, subsizes_fxy, disp, offset, f2)
            disp = disp + nx_mhd_reduced * ny_mhd_reduced * 8
            call write_data_mpi_io(fh, mpi_datatype, subsizes_fxy, disp, offset, f3)
            disp = disp + nx_mhd_reduced * ny_mhd_reduced * 8
            call write_data_mpi_io(fh, mpi_datatype, subsizes_fxy, disp, offset, f4)
            disp = disp + nx_mhd_reduced * ny_mhd_reduced * 8
            call write_data_mpi_io(fh, mpi_datatype, subsizes_fxy, disp, offset, f5)
            call MPI_FILE_CLOSE(fh, ierror)

            if (local_dist) then
                fh = 19
                fname = trim(file_path)//'fp_local_'//ctime//'_sum.dat'
                disp = 0
                offset = 0
                mpi_datatype = set_mpi_datatype_double(sizes_fp_local, subsizes_fp_local, starts_fp_local)
                call open_data_mpi_io(trim(fname), MPI_MODE_CREATE+MPI_MODE_WRONLY, fileinfo, fh)
                call write_data_mpi_io(fh, mpi_datatype, subsizes_fp_local, disp, offset, fp_local)
                call MPI_FILE_CLOSE(fh, ierror)
            endif
        endif

        call clean_particle_distributions
        if (local_dist) then
            call clean_local_particle_distribution
        endif
    end subroutine distributions_diagnostics

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

        !< Set particles to their initial values, but switch the tags for
        !< for high-energy particles
        ptls%x = 0.0
        ptls%y = 0.0
        ptls%p = 0.0
        ptls%weight = 0.0
        ptls%t = 0.0
        ptls%dt = 0.0
        ptls%split_times = 0
        ptls%count_flag = 0
        ptls%tag = 0
        ptls%nsteps_tracking = 0
    end subroutine select_particles_tracking

    !---------------------------------------------------------------------------
    !< Make the flags of tracked particles negative, so they can be easily
    !< identified.
    !< Args:
    !<  nptl_selected: number of selected particles
    !---------------------------------------------------------------------------
    subroutine negative_particle_tags(nptl_selected)
        implicit none
        integer, intent(in) :: nptl_selected
        integer :: i
        do i = 1, nptl_selected
            ptls(tags_selected_ptls(i))%tag = -i
        enddo
    end subroutine negative_particle_tags

    !---------------------------------------------------------------------------
    !< Initialize tracked particle points
    !< Args:
    !<  nptl_selected: number of selected particles
    !---------------------------------------------------------------------------
    subroutine init_tracked_particle_points(nptl_selected)
        implicit none
        integer, intent(in) :: nptl_selected
        integer :: i, offset, tag
        allocate(ptl_traj_points(nsteps_tracked_tot))
        ptl_traj_points%x = 0.0
        ptl_traj_points%y = 0.0
        ptl_traj_points%p = 0.0
        ptl_traj_points%weight = 0.0
        ptl_traj_points%t = 0.0
        ptl_traj_points%dt = 0.0
        ptl_traj_points%split_times = 0
        ptl_traj_points%count_flag = 0
        ptl_traj_points%tag = 0
        ptl_traj_points%nsteps_tracking = 0

        !< Save the initial information of tracked particles
        do i = 1, nptl_selected
            offset = noffsets_tracked_ptls(i)
            tag = tags_selected_ptls(i)
            ptl_traj_points(offset + 1) = ptls(tag)
        enddo
    end subroutine init_tracked_particle_points

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
        use mpi_module
        implicit none
        integer, intent(in) :: nptl_selected
        character(*), intent(in) :: file_path
        character(len=4) :: mrank
        character(len=128) :: fname
        integer :: fh, offset
        write (mrank,'(i4.4)') mpi_rank
        fh = 41
        fname = trim(file_path)//'tracked_particle_points_'//mrank//'.dat'
        open(unit=fh, file=fname, access='stream', status='unknown', &
            form='unformatted', action='write')
        write(fh, pos=1) nptl_selected
        write(fh, pos=5) nsteps_tracked_ptls
        offset = (nptl_selected + 1) * sizeof(fp)
        write(fh, pos=offset+1) ptl_traj_points
        close(fh)
    end subroutine save_tracked_particle_points

    !---------------------------------------------------------------------------
    !< Set MPI/IO data sizes for distribution diagnostics
    !---------------------------------------------------------------------------
    subroutine set_mpi_io_data_sizes
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: mpi_ix, mpi_iy
        implicit none
        sizes_fpx = (/ npp, nx_mhd_reduced /)
        subsizes_fpx = (/ npp, nx /)
        starts_fpx = (/ 0, nx * mpi_ix /)

        sizes_fpy = (/ npp, ny_mhd_reduced /)
        subsizes_fpy = (/ npp, ny /)
        starts_fpy = (/ 0, ny * mpi_iy /)

        sizes_fxy = (/ nx_mhd_reduced, ny_mhd_reduced /)
        subsizes_fxy = (/ nx, ny /)
        starts_fxy = (/ nx * mpi_ix, ny * mpi_iy /)

        sizes_fp_local = (/ npp, nx_mhd_reduced, ny_mhd_reduced /)
        subsizes_fp_local = (/ npp, nx, ny /)
        starts_fp_local = (/ 0, nx * mpi_ix, ny * mpi_iy /)
    end subroutine set_mpi_io_data_sizes
end module particle_module

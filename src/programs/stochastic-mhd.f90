!*******************************************************************************
!< Program for stochastic integration
!*******************************************************************************

program stochastic
    use constants, only: dp
    use mpi_module
    use omp_lib
    use mhd_config_module, only: load_mhd_config, mhd_config, &
        echo_mhd_config, set_mhd_grid_type, init_tstamps_mhd, &
        free_tstamps_mhd, load_tstamps_mhd, calc_tstamps_mhd
    use particle_module, only: init_particles, free_particles, &
        inject_particles_spatial_uniform, read_particle_params, &
        particle_mover, split_particle, set_particle_datatype_mpi, &
        free_particle_datatype_mpi, init_particle_tracking, &
        free_particle_tracking, dump_tracked_particles, &
        inject_particles_at_shock, inject_particles_at_large_jz, &
        set_dpp_params, set_duu_params, set_flags_params, &
        set_drift_parameters, set_flag_check_drift_2d, &
        inject_particles_at_large_db2, &
        inject_particles_at_large_divv, &
        inject_particles_at_large_rho, &
        inject_particles_at_large_absj, &
        read_particles, &
        save_particle_module_state, read_particle_module_state
    use diagnostics, only: distributions_diagnostics, quick_check, &
        init_particle_distributions, free_particle_distributions, &
        get_pmax_global, dump_particles, &
        init_escaped_particles, free_escaped_particles, &
        reset_escaped_particles, dump_escaped_particles, &
        read_diagnostics_params
    use random_number_generator, only: init_prng, delete_prng, &
        save_prng
    use mhd_data_parallel, only: init_field_data, free_field_data, &
        read_field_data_parallel, calc_fields_gradients, &
        calc_fields_gradients_nonuniform, &
        copy_fields, init_shock_xpos, free_shock_xpos, &
        locate_shock_xpos, init_magnetic_fluctuation, free_magnetic_fluctuation, &
        read_magnetic_fluctuation, copy_magnetic_fluctuation, &
        init_correlation_length, free_correlation_length, &
        read_correlation_length, copy_correlation_length, &
        calc_grad_deltab2, calc_grad_deltab2_nonuniform, &
        calc_grad_correl_length, calc_grad_correl_length_nonuniform, &
        init_grid_positions, free_grid_positions, &
        set_local_grid_positions
    use simulation_setup_module, only: read_simuation_mpi_topology, &
        set_field_configuration, fconfig, read_particle_boundary_conditions, &
        set_neighbors, check_particle_can_escape
    use acc_region_surface, only: init_acc_surface, free_acc_surface, &
        read_acc_surface, copy_acc_surface
    use flap, only : command_line_interface !< FLAP package
    use penf
    implicit none
    character(len=256) :: dir_mhd_data, filename
    character(len=128) :: diagnostics_directory
    character(len=128) :: particle_tags_file
    character(len=64) :: conf_file, mhd_config_filename
    character(len=64) :: surface_filename1, surface_filename2
    character(len=2) :: surface_norm1, surface_norm2 ! The 1st character is the orientation
    integer :: nptl_max, nptl
    real(dp) :: start, finish, uptime, step1, step2, dt
    real(dp) :: jz_min, absj_min, db2_min, divv_min, rho_min
    real(dp) :: ptl_xmin, ptl_xmax, ptl_ymin, ptl_ymax, ptl_zmin, ptl_zmax
    real(dp) :: tau0_scattering ! Scattering time for initial particles
    real(dp) :: drift_param1, drift_param2 ! Drift parameter for 3D simulation
    real(dp) :: power_index ! Power-law spectrum index for initial distribution
    real(dp) :: split_ratio ! Momentum increase ratio for particle splitting
    real(dp) :: pmin_split  ! The minimum momentum (in terms in p0) to start splitting
    real(dp) :: particle_v0 ! Initial particle velocity in the normalized velocity
    real(dp) :: duu0        ! Pitch-angle diffusion for particles with a momentum of p0
    real(dp) :: quota_hour  ! The maximum wall time in hours for the simulation
    real(dp), dimension(6) :: part_box
    integer :: ncells_large_jz_norm   ! Normalization for the number of cells with large jz
    integer :: ncells_large_absj_norm ! Normalization for the number of cells with large absj
    integer :: ncells_large_db2_norm  ! Normalization for the number of cells with large db2
    integer :: ncells_large_divv_norm ! Normalization for the number of cells with large divv
    integer :: ncells_large_rho_norm  ! Normalization for the number of cells with large rho
    integer :: nthreads, color
    integer :: t_start, t_end, tf, split_flag
    integer :: tmin  ! tmin = t_start typically. When restarting, tmin > t_start
    integer :: tmax_mhd           ! Maximum time frame for the MHD fields
    integer :: time_interp_flag   ! Whether to interpolate between two MHD frames
    integer :: nsteps_interval    ! Steps interval to track particles
    integer :: dist_flag          ! 0 for Maxwellian. 1 for delta function. 2 for power-law
    integer :: num_fine_steps     ! Number of the fine steps for diagnostics
    integer :: tmax_to_inject     ! Maximum time frame to inject particles
    integer :: dpp_wave           ! Whether to include momentum diffusion due to wave scattering
    integer :: dpp_shear          ! Whether to include momentum diffusion due to flow shear
    integer :: weak_scattering    ! Whether particle scattering is weak (tau*Omega >> 1)
    integer :: deltab_flag        ! Whether to include spatially dependent magnetic fluctuation
    integer :: correlation_flag   ! Whether to include turbulence correlation length
    integer :: ndim_field         ! Number of the dimension of the MHD fields
    integer :: charge             ! Particle change in the unit of e
    integer :: spherical_coord    ! Whether MHD simulations are in spherical coordinates
    integer :: uniform_grid       ! Whether the grid is uniform
    integer :: check_drift_2d     ! Whether to check particle drift in 2D simulations
    integer :: particle_data_dump ! Whether to dump particle data
    integer :: size_mpi_sub       ! Size of a MPI sub-communicator
    integer :: single_time_frame  ! Whether to use a single time frame of fields
    integer :: include_3rd_dim    ! Whether to include transport along the 3rd-dim in 2D runs
    integer :: acc_by_surface     ! Whether the acceleration region is separated by surfaces
    logical :: inject_at_shock    ! Inject particles at shock location
    logical :: inject_new_ptl     ! Whether to inject new particles every MHD step
    logical :: inject_large_jz    ! Whether to inject where jz is large
    logical :: inject_large_absj  ! Whether to inject where absj is large
    logical :: inject_large_db2   ! Whether to inject where db2 is large
    logical :: inject_large_divv  ! Whether to inject where divv is negatively large
    logical :: inject_large_rho   ! Whether to inject where rho is large
    logical :: inject_same_nptl   ! Whether to inject same of number particles every step
    logical :: inject_part_box    ! Whether to inject in part of the box
    logical :: restart_flag       ! Whether to restart a previous simulation
    logical :: local_dist         ! Whether to dump local particle distributions
    logical :: surface2_existed   ! Whether surface 2 is existed
    logical :: is_intersection    ! Intersection or Union of the two surfaces
    logical :: varying_dt_mhd     ! Whether the time interval for MHD fields is varying
    logical :: focused_transport  ! Whether to the Focused Transport equation
    logical :: dump_escaped_dist  ! Whether to dump the distributions of the escaped particles
    logical :: dump_escaped       ! Whether to dump escaped particles
    logical :: reached_quota      ! Whether the simulation reach the quota time
    logical :: track_particle_flag! Whether to track particles

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)

    start = MPI_Wtime()

    call get_cmd_args

    if (single_time_frame == 1) then
        time_interp_flag = 0  ! no time interpolation for a single time frame
    endif

    ! Split the communicator based on the color and use the original rank for ordering
    color = mpi_rank / size_mpi_sub
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, mpi_rank, mpi_sub_comm, ierr);
    call MPI_COMM_RANK(mpi_sub_comm, mpi_sub_rank, ierr)
    call MPI_COMM_SIZE(mpi_sub_comm, mpi_sub_size, ierr)
    if (mpi_rank == master) then
        print '(A,I10.3)', 'mpi_sub_size: ', mpi_sub_size
    endif

    ! Group the same mpi_sub_rank in each mpi_sub_comm into another communicator
    color = mpi_sub_rank
    call MPI_COMM_SPLIT(MPI_COMM_WORLD, color, mpi_rank, mpi_cross_comm, ierr);
    call MPI_COMM_RANK(mpi_cross_comm, mpi_cross_rank, ierr)
    call MPI_COMM_SIZE(mpi_cross_comm, mpi_cross_size, ierr)
    if (mpi_rank == master) then
        print '(A,I10.3)', 'mpi_cross_size: ', mpi_cross_size
    endif

    nthreads = 1
    !$OMP PARALLEL
    !$OMP MASTER
#if (defined USE_OPENMP)
    nthreads = OMP_GET_NUM_THREADS()
#endif
    !$OMP END MASTER
    !$OMP END PARALLEL
    call init_prng(mpi_rank, nthreads, restart_flag, diagnostics_directory)

    !< Configurations
    filename = trim(dir_mhd_data)//trim(mhd_config_filename)
    call load_mhd_config(filename)
    if (mpi_rank == master) then
        call echo_mhd_config
    endif
    ! The time stamps for the MHD fields
    call init_tstamps_mhd(t_start, t_end)
    if (varying_dt_mhd) then
        call load_tstamps_mhd(t_start, t_end, tmax_mhd, dir_mhd_data)
    else
        call calc_tstamps_mhd(t_start, t_end)
    endif
    call set_mhd_grid_type(uniform_grid, spherical_coord)
    call read_simuation_mpi_topology(conf_file)
    call read_particle_boundary_conditions(conf_file)
    call set_field_configuration(ndim=ndim_field)
    call set_neighbors
    call check_particle_can_escape
    call init_grid_positions
    call set_local_grid_positions(dir_mhd_data)

    !< Initialization
    call init_field_data(time_interp_flag)
    if (deltab_flag == 1) then
        call init_magnetic_fluctuation
    endif
    if (correlation_flag == 1) then
        call init_correlation_length
    endif
    if (acc_by_surface == 1) then
        if (surface2_existed) then
            call init_acc_surface(time_interp_flag, surface_norm1, surface_norm2, &
                is_intersection)
        else
            call init_acc_surface(time_interp_flag, surface_norm1)
        endif
    endif
    call read_particle_params(conf_file)
    call read_diagnostics_params(conf_file, focused_transport, t_end-t_start)
    call set_dpp_params(dpp_wave, dpp_shear, weak_scattering, tau0_scattering)
    call set_duu_params(duu0) ! duu0 is not used when solving Parker's transport
    call set_flags_params(deltab_flag, correlation_flag, include_3rd_dim, &
        acc_by_surface)
    call set_drift_parameters(drift_param1, drift_param2, charge)
    if (ndim_field == 2 .and. check_drift_2d == 1) then
        call set_flag_check_drift_2d(check_drift_2d)
    endif

    call init_particles(nptl_max)
    if (dump_escaped_dist) then
        call init_escaped_particles
    endif
    call init_particle_distributions(local_dist, dump_escaped_dist)
    call set_particle_datatype_mpi

    reached_quota = .false.  ! might be updated in solve_transport_equation

    if (restart_flag) then
        ! Read the time frame before previous restart
        if (mpi_rank == master) then
            open(unit=20,file=trim(diagnostics_directory)//"restart/latest_restart",&
                access='stream',status='unknown',form='unformatted',action='read')
            read(20, pos=1) tmin
            close(20)
        endif
        call MPI_BCAST(tmin, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)

        !< Read particles
        call read_particles(tmin, trim(diagnostics_directory)//"restart/")
        call read_particle_module_state(tmin, trim(diagnostics_directory)//"restart/")
    else
        tmin = t_start
    endif

    if (track_particle_flag) then
        filename = trim(diagnostics_directory)//"particle_tracking/"//trim(particle_tags_file)
        call init_particle_tracking(trim(filename), nsteps_interval)
    endif

    call solve_transport_equation(track_particle_flag)

    if (track_particle_flag) then
        call free_particle_tracking
    endif

    ! dump files for restart if necessary
    if (.not. reached_quota) then
        tf =  tf - 1  ! when finished time loop, tf is one over t_end
        if (mpi_rank == master) then
            print '(A)', 'Dumping restart files at the end of the simulation'
        endif
    else
        if (mpi_rank == master) then
            print '(A)', 'Reached quota time. Dumping restart files.'
        endif
    endif
    call save_prng(mpi_rank, nthreads, diagnostics_directory)
    call dump_particles(t_end, trim(diagnostics_directory)//"restart/")
    call save_particle_module_state(t_end, trim(diagnostics_directory)//"restart/")
    if (mpi_rank == master) then
        ! Write the finished time frame
        open(unit=20,file=trim(diagnostics_directory)//"restart/latest_restart",&
            access='stream',status='unknown',form='unformatted',action='write')
        write(20, pos=1) tf
        close(20)
    endif

    if (inject_at_shock) then
        call free_shock_xpos
    endif

    call free_particle_datatype_mpi
    call free_particle_distributions(local_dist, dump_escaped_dist)
    call free_particles
    if (dump_escaped_dist) then
        call free_escaped_particles
    endif
    if (acc_by_surface == 1) then
        call free_acc_surface
    endif
    call free_field_data
    if (deltab_flag == 1) then
        call free_magnetic_fluctuation
    endif
    if (correlation_flag == 1) then
        call free_correlation_length
    endif

    call free_tstamps_mhd
    call free_grid_positions

    finish = MPI_Wtime()
    if (mpi_rank == master) then
        print '("Time = ",f14.7," seconds.")',finish-start
    endif

    call delete_prng

    call MPI_FINALIZE(ierr)

    contains

    !---------------------------------------------------------------------------
    !< Loop over all time steps and solve the transport equation
    !---------------------------------------------------------------------------
    subroutine solve_transport_equation(track_particles)
        implicit none
        logical, intent(in) :: track_particles
        character(len=256) :: mhd_fname1, mhd_fname2
        character(len=256) :: acc_surface_fname11, acc_surface_fname12
        character(len=256) :: acc_surface_fname21, acc_surface_fname22

        ! Read the MHD fields
        write(mhd_fname1, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', tmin
        call read_field_data_parallel(mhd_fname1, 0)

        ! Read the surface data to separate the regions
        if (acc_by_surface == 1) then
            write(acc_surface_fname11, "(A,I4.4,A)") &
                trim(dir_mhd_data)//trim(surface_filename1)//'_', tmin, ".dat"
            if (surface2_existed) then
                write(acc_surface_fname21, "(A,I4.4,A)") &
                    trim(dir_mhd_data)//trim(surface_filename2)//'_', tmin, ".dat"
                call read_acc_surface(0, acc_surface_fname11, acc_surface_fname21)
            else
                call read_acc_surface(0, acc_surface_fname11)
            endif
        endif

        ! Read the turbulent fluctuations if needed
        if (deltab_flag == 1) then
            write(mhd_fname1, "(A,I4.4)") trim(dir_mhd_data)//'deltab_', tmin
            call read_magnetic_fluctuation(mhd_fname1, 0)
        endif

        ! Read the turbulent correlation length if needed
        if (correlation_flag == 1) then
            write(mhd_fname1, "(A,I4.4)") trim(dir_mhd_data)//'lc_', tmin
            call read_correlation_length(mhd_fname1, 0)
        endif

        ! Calculate the gradients
        if (uniform_grid == 1) then
            call calc_fields_gradients(0)
        else
            call calc_fields_gradients_nonuniform(0)
        endif
        if (deltab_flag == 1) then
            if (uniform_grid == 1) then
                call calc_grad_deltab2(0)
            else
                call calc_grad_deltab2_nonuniform(0)
            endif
        endif
        if (correlation_flag == 1) then
            if (uniform_grid == 1) then
                call calc_grad_correl_length(0)
            else
                call calc_grad_correl_length_nonuniform(0)
            endif
        endif

        ! Prepare for injecting particles
        if (inject_at_shock) then ! Whether it is a shock problem
            call init_shock_xpos
        else
            if (inject_part_box) then  ! only inject in part of the box
                part_box(1) = ptl_xmin
                part_box(2) = ptl_ymin
                part_box(3) = ptl_zmin
                part_box(4) = ptl_xmax
                part_box(5) = ptl_ymax
                part_box(6) = ptl_zmax
            else
                part_box(1) = fconfig%xmin
                part_box(2) = fconfig%ymin
                part_box(3) = fconfig%zmin
                part_box(4) = fconfig%xmax
                part_box(5) = fconfig%ymax
                part_box(6) = fconfig%zmax
            endif
        endif

        step1 = MPI_Wtime()

        !< Time loop
        do tf = tmin+1, t_end
            if (mpi_rank == master) then
                write(*, "(A,I0)") " Starting step ", tf
            endif
            if (single_time_frame == 0 .and. tf <= tmax_mhd) then
                ! Read the second frame of fields data for time interpolation
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                write(mhd_fname2, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', tf
                call read_field_data_parallel(mhd_fname2, var_flag=time_interp_flag)
                if (acc_by_surface == 1) then
                    write(acc_surface_fname12, "(A,I4.4,A)") &
                        trim(dir_mhd_data)//trim(surface_filename1)//'_', tf, ".dat"
                    if (surface2_existed) then
                        write(acc_surface_fname22, "(A,I4.4,A)") &
                            trim(dir_mhd_data)//trim(surface_filename2)//'_', tf, ".dat"
                        call read_acc_surface(time_interp_flag, acc_surface_fname12, acc_surface_fname22)
                    else
                        call read_acc_surface(time_interp_flag, acc_surface_fname12)
                    endif
                endif
                if (deltab_flag == 1) then
                    write(mhd_fname2, "(A,I4.4)") trim(dir_mhd_data)//'deltab_', tf
                    call read_magnetic_fluctuation(mhd_fname2, var_flag=time_interp_flag)
                endif
                if (correlation_flag == 1) then
                    write(mhd_fname2, "(A,I4.4)") trim(dir_mhd_data)//'lc_', tf
                    call read_correlation_length(mhd_fname2, var_flag=time_interp_flag)
                endif
                if (uniform_grid == 1) then
                    call calc_fields_gradients(var_flag=time_interp_flag)
                else
                    call calc_fields_gradients_nonuniform(var_flag=time_interp_flag)
                endif
                if (deltab_flag == 1) then
                    if (uniform_grid == 1) then
                        call calc_grad_deltab2(var_flag=time_interp_flag)
                    else
                        call calc_grad_deltab2_nonuniform(var_flag=time_interp_flag)
                    endif
                endif
                if (correlation_flag == 1) then
                    if (uniform_grid == 1) then
                        call calc_grad_correl_length(var_flag=time_interp_flag)
                    else
                        call calc_grad_correl_length_nonuniform(var_flag=time_interp_flag)
                    endif
                endif
            endif  ! if (single_time_frame == 0)

            ! Inject particles at the starting point of the this time interval
            if (inject_at_shock) then
                call locate_shock_xpos
                call inject_particles_at_shock(nptl, dt, dist_flag, &
                    particle_v0, tf-t_start, power_index)
            else
                if (tf == t_start+1 .or. inject_new_ptl) then
                    if (tf <= tmax_to_inject) then
                        if (inject_large_jz) then
                            call inject_particles_at_large_jz(nptl, dt, dist_flag, &
                                particle_v0, tf-t_start, inject_same_nptl, jz_min, &
                                ncells_large_jz_norm, part_box, power_index)
                        else if (inject_large_absj) then
                            call inject_particles_at_large_absj(nptl, dt, dist_flag, &
                                particle_v0, tf-t_start, inject_same_nptl, absj_min, &
                                ncells_large_absj_norm, part_box, power_index)
                        else if (inject_large_db2) then
                            call inject_particles_at_large_db2(nptl, dt, dist_flag, &
                                particle_v0, tf-t_start, inject_same_nptl, db2_min, &
                                ncells_large_db2_norm, part_box, power_index)
                        else if (inject_large_divv) then
                            call inject_particles_at_large_divv(nptl, dt, dist_flag, &
                                particle_v0, tf-t_start, inject_same_nptl, divv_min, &
                                ncells_large_divv_norm, part_box, power_index)
                        else if (inject_large_rho) then
                            call inject_particles_at_large_rho(nptl, dt, dist_flag, &
                                particle_v0, tf-t_start, inject_same_nptl, rho_min, &
                                ncells_large_rho_norm, part_box, power_index)
                        else
                            call inject_particles_spatial_uniform(nptl, dt, dist_flag, &
                                particle_v0, tf-t_start, part_box, power_index)
                        endif
                    endif
                endif
            endif

            ! Initial diagnostics if we are not tracking particles
            if (tf == t_start + 1 .and. (.not. track_particles)) then
                call distributions_diagnostics(t_start, diagnostics_directory, &
                    local_dist, dump_escaped_dist)
                if (particle_data_dump == 1) then
                    call dump_particles(t_start, diagnostics_directory)
                endif
                call quick_check(t_start, .true., diagnostics_directory)
                call get_pmax_global(t_start, .true., diagnostics_directory)
            endif

            ! Move particles
            if (track_particles) then
                call particle_mover(focused_transport, nsteps_interval, &
                    tf-t_start, 1, dump_escaped_dist)
            else
                call particle_mover(focused_transport, nsteps_interval, &
                    tf-t_start, num_fine_steps, dump_escaped_dist)
            endif

            if (mpi_rank == master) then
                write(*, "(A)") " Finishing moving particles "
            endif

            if (track_particles) then
                call dump_tracked_particles(tf, trim(diagnostics_directory)//"particle_tracking/")
            endif

            if (split_flag == 1) then
                call split_particle(split_ratio, pmin_split, nsteps_interval)
            endif
            if (.not. track_particles) then
                call quick_check(tf, .false., diagnostics_directory)
                call get_pmax_global(tf, .false., diagnostics_directory)
                call distributions_diagnostics(tf, diagnostics_directory, &
                    local_dist, dump_escaped_dist)
                if (mpi_rank == master) then
                    write(*, "(A)") " Finishing distribution diagnostics "
                endif
                if (particle_data_dump == 1) then
                    call dump_particles(tf, diagnostics_directory)
                endif
                if (dump_escaped_dist) then
                    if (dump_escaped) then
                        call dump_escaped_particles(tf, diagnostics_directory)
                    endif
                    call reset_escaped_particles
                endif
            endif

            ! Copy fields for next step if necessary
            if (time_interp_flag == 1) then
                call copy_fields
                if (acc_by_surface == 1) then
                    call copy_acc_surface
                endif
                if (deltab_flag == 1) then
                    call copy_magnetic_fluctuation
                endif
                if (correlation_flag == 1) then
                    call copy_correlation_length
                endif
                if (mpi_rank == master) then
                    write(*, "(A)") " Finishing copying fields "
                endif
            endif

            step2 = MPI_Wtime()
            if (mpi_rank == master) then
                print '("Step ", I0, " takes ", f9.4, " seconds.")', tf, step2 - step1
            endif
            step1 = step2

            ! Check whether to dump restart files
            uptime = step2 - start
            if (uptime > (quota_hour - 0.5) * 3600) then
                !< 30 mins before the quota time limit
                reached_quota = .true.
                return
            endif
        enddo  ! Time loop
    end subroutine solve_transport_equation

    !---------------------------------------------------------------------------
    !< Get commandline arguments
    !---------------------------------------------------------------------------
    subroutine get_cmd_args
        implicit none
        type(command_line_interface) :: cli     !< Command Line Interface (CLI).
        integer(I4P)                 :: error   !< Error trapping flag.
        call cli%init(&
            progname = 'stochastic-mhd', &
            authors     = 'Xiaocan Li', &
            help        = 'Usage: ', &
            description = "Solving Parker's transport equation "// &
                          "using stochastic differential equation", &
            examples = ['stochastic-mhd.exec '// &
                        '-qh quota_hour -rf restart_flag '//&
                        '-ft focused_transport -pv particle_v0 '//&
                        '-sm size_mpi_sub '//&
                        '-dm dir_mhd_data -mc mhd_config_filename '//&
                        '-nm nptl_max -np nptl '//&
                        '-ti time_interp_flag -dt dt -ts t_start -te t_end '//&
                        '-tm tmax_mhd '//&
                        '-st single_time_frame -df dist_flag -pi power_index '//&
                        '-sf split_flag -sr split_ratio -ps pmin_split '//&
                        '-tf track_particle_flag '//&
                        '-ptf particle_tags_file -ni nsteps_interval '//&
                        '-dd diagnostics_directory -is inject_at_shock '//&
                        '-in inject_new_ptl -ij inject_large_jz '//&
                        '-sn inject_same_nptl '//&
                        '-tti tmax_to_inject '//&
                        '-ip inject_part_box -jz jz_min '//&
                        '-nn ncells_large_jz_norm -ib inject_large_db2 '//&
                        '-db2 db2_min -nb ncells_large_db2_norm '//&
                        '-iv inject_large_divv -dv divv_min '//&
                        '-nv ncells_large_divv_norm '//&
                        '-ir inject_large_rho -rm rho_min '//&
                        '-nr ncells_large_rho_norm '//&
                        '-iaj inject_large_absj -ajm absj_min '//&
                        '-naj ncells_large_absj_norm '//&
                        '-xs ptl_xmin -xe ptl_xmax '//&
                        '-ys ptl_ymin -ye ptl_ymax -zs ptl_zmin -ze ptl_zmax '//&
                        '-cf conf_file -nf num_fine_steps '//&
                        '-ld local_dist -dw dpp_wave '//&
                        '-ds dpp_shear -ws weak_scattering -t0 tau0_scattering '//&
                        '-db deltab_flag -co correlation_flag -nd ndim_field '//&
                        '-dp1 drift_param1 -dp2 drift_param2 -ch charge '//&
                        '-sc spherical_coord -ug uniform_grid '//&
                        '-cd check_drift_2d -pd particle_data_dump '//&
                        '-i3 include_3rd_dim -as acc_by_surface '//&
                        '-s2e surface2_existed -ii is_intersection '//&
                        '-sf1 surface_filename1 -sn1 surface_norm1 '//&
                        '-sf2 surface_filename2 -sn2 surface_norm2 '//&
                        '-vdt varying_dt_mhd -du duu0 '//&
                        '-ded dump_escaped_dist -de dump_escaped'])
        call cli%add(switch='--quota_hour', switch_ab='-qh', &
            help='The maximum wall time for the simulation', &
            required=.false., act='store', def='48.0', error=error)
        if (error/=0) stop
        call cli%add(switch='--restart_flag', switch_ab='-rf', &
            help='whether to restart a previous simulation', &
            required=.false., act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--focused_transport', switch_ab='-ft', &
            help='whether to solve the Focused Transport Equation', &
            required=.false., act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--particle_v0', switch_ab='-pv', &
            help='Initial particle velocity in the normalized velocity', &
            required=.false., act='store', def='5.00', error=error)
        if (error/=0) stop
        call cli%add(switch='--size_mpi_sub', switch_ab='-sm', &
            help='Size of the a MPI sub-communicator', required=.false., &
            act='store', def='1', error=error)
        if (error/=0) stop
        call cli%add(switch='--dir_mhd_data', switch_ab='-dm', &
            help='MHD simulation data file directory', required=.true., &
            act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--mhd_config_filename', switch_ab='-mc', &
            help='MHD configuration filename', required=.false., &
            act='store', def='mhd_config.dat', error=error)
        if (error/=0) stop
        call cli%add(switch='--nptl_max', switch_ab='-nm', &
            help='Maximum number of particles', required=.false., &
            act='store', def='1E7', error=error)
        if (error/=0) stop
        call cli%add(switch='--nptl', switch_ab='-np', &
            help='Number of particles', required=.false., &
            act='store', def='1E4', error=error)
        if (error/=0) stop
        call cli%add(switch='--time_interp_flag', switch_ab='-ti', &
            help='Flag for time interpolation', required=.false., &
            act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--tinterval', switch_ab='-dt', &
            help='Time interval to push particles', required=.false., &
            act='store', def='1E-7', error=error)
        if (error/=0) stop
        call cli%add(switch='--tstart', switch_ab='-ts', &
            help='Starting time frame', required=.false., &
            act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--tend', switch_ab='-te', &
            help='The last time frame', required=.false., &
            act='store', def='200', error=error)
        if (error/=0) stop
        call cli%add(switch='--tmax_mhd', switch_ab='-tm', &
            help='The maximum time frame for the MHD fields', required=.false., &
            act='store', def='100000', error=error)
        if (error/=0) stop
        call cli%add(switch='--single_time_frame', switch_ab='-st', &
            help='Whether to use a single time frame of fields', required=.false., &
            act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--dist_flag', switch_ab='-df', &
            help='Flag to indicate momentum distribution', required=.false., &
            act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--power_index', switch_ab='-pi', &
            help='Power-law spectrum index for initial distribution', &
            required=.false., act='store', def='7.0', error=error)
        if (error/=0) stop
        call cli%add(switch='--split_flag', switch_ab='-sf', &
            help='Flag to split the particles', required=.false., &
            act='store', def='1', error=error)
        if (error/=0) stop
        call cli%add(switch='--split_ratio', switch_ab='-sr', &
            help='(> 1.0) momentum increase ratio for particle splitting', &
            required=.false., act='store', def='2.72', error=error)
        if (error/=0) stop
        call cli%add(switch='--pmin_split', switch_ab='-ps', &
            help='(> 1.0) the minimum momentum (in terms in p0) to start splitting particles', &
            required=.false., act='store', def='2.0', error=error)
        if (error/=0) stop
        call cli%add(switch='--track_particle_flag', switch_ab='-tf', &
            help='Whether to track some particles', required=.false., &
            act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--particle_tags_file', switch_ab='-ptf', &
            help='File containing the selected particle tags', required=.false., &
            act='store', def='tags_selected_01.h5', error=error)
        if (error/=0) stop
        call cli%add(switch='--nsteps_interval', switch_ab='-ni', &
            help='Steps interval to track particles', required=.false., &
            act='store', def='10', error=error)
        if (error/=0) stop
        call cli%add(switch='--diagnostics_directory', switch_ab='-dd', &
            help='Diagnostics directory', required=.false., &
            act='store', def='data/', error=error)
        if (error/=0) stop
        call cli%add(switch='--inject_at_shock', switch_ab='-is', &
            help='whether to inject particles at shock', required=.false., &
            act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--inject_new_ptl', switch_ab='-in', &
            help='whether to inject new particles every MHD step', required=.false., &
            act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--inject_large_jz', switch_ab='-ij', &
            help='whether to inject where jz is large', required=.false., &
            act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--inject_same_nptl', switch_ab='-sn', &
            help='whether to inject the same number of particles every step', &
            required=.false., act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--tmax_to_inject', switch_ab='-tti', &
            help='The maximum time frame to inject particles', required=.false., &
            act='store', def='100000', error=error)
        if (error/=0) stop
        call cli%add(switch='--inject_part_box', switch_ab='-ip', &
            help='whether to inject in a part of the box', required=.false., &
            act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--jz_min', switch_ab='-jz', &
            help='Minimum jz in regions to inject new particles', &
            required=.false., def='100.0', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--ncells_large_jz_norm', switch_ab='-nn', &
            help='Normalization for the number of cells with large jz', &
            required=.false., def='800', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--inject_large_db2', switch_ab='-ib', &
            help='whether to inject where db2 is large', required=.false., &
            act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--db2_min', switch_ab='-db2', &
            help='Minimum db2 in regions to inject new particles', &
            required=.false., def='0.03', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--ncells_large_db2_norm', switch_ab='-nb', &
            help='Normalization for the number of cells with large db2', &
            required=.false., def='800', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--inject_large_divv', switch_ab='-iv', &
            help='whether to inject where divv is large', required=.false., &
            act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--divv_min', switch_ab='-dv', &
            help='Minimum divv in regions to inject new particles', &
            required=.false., def='10.0', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--ncells_large_divv_norm', switch_ab='-nv', &
            help='Normalization for the number of cells with large divv', &
            required=.false., def='800', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--inject_large_rho', switch_ab='-ir', &
            help='whether to inject where rho is large', required=.false., &
            act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--rho_min', switch_ab='-rm', &
            help='Minimum rho in regions to inject new particles', &
            required=.false., def='5.0', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--ncells_large_rho_norm', switch_ab='-nr', &
            help='Normalization for the number of cells with large rho', &
            required=.false., def='800', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--inject_large_absj', switch_ab='-iaj', &
            help='whether to inject where absj is large', required=.false., &
            act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--absj_min', switch_ab='-ajm', &
            help='Minimum absj in regions to inject new particles', &
            required=.false., def='5.0', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--ncells_large_absj_norm', switch_ab='-naj', &
            help='Normalization for the number of cells with large absj', &
            required=.false., def='800', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--ptl_xmin', switch_ab='-xs', &
            help='Minimum x of the region to inject new particles', &
            required=.false., def='0.0', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--ptl_xmax', switch_ab='-xe', &
            help='Maximum x of the region to inject new particles', &
            required=.false., def='1.0', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--ptl_ymin', switch_ab='-ys', &
            help='Minimum y of the region to inject new particles', &
            required=.false., def='0.0', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--ptl_ymax', switch_ab='-ye', &
            help='Maximum y of the region to inject new particles', &
            required=.false., def='1.0', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--ptl_zmin', switch_ab='-zs', &
            help='Minimum z of the region to inject new particles', &
            required=.false., def='0.0', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--ptl_zmax', switch_ab='-ze', &
            help='Maximum z of the region to inject new particles', &
            required=.false., def='1.0', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--conf_file', switch_ab='-cf', &
            help='Configuration file name', required=.false., &
            act='store', def='conf.dat', error=error)
        if (error/=0) stop
        call cli%add(switch='--num_fine_steps', switch_ab='-nf', &
            help='number of fine steps', required=.false., &
            act='store', def='1', error=error)
        if (error/=0) stop
        call cli%add(switch='--local_dist', switch_ab='-ld', &
            help='whether to diagnose local distribution', &
            required=.false., act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--dpp_wave', switch_ab='-dw', &
            help='whether to include momentum diffusion due to waves', &
            required=.false., act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--dpp_shear', switch_ab='-ds', &
            help='whether to include momentum diffusion due to flow shear', &
            required=.false., act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--weak_scattering', switch_ab='-ws', &
            help='whether particle scattering is weak (tau*Omega >> 1)', &
            required=.false., act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--tau0_scattering', switch_ab='-t0', &
            help='The scattering time for initial particles', &
            required=.false., def='1.0', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--deltab_flag', switch_ab='-db', &
            help='whether to include spatially dependent magnetic fluctuation', &
            required=.false., act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--correlation', switch_ab='-co', &
            help='whether to include turbulence correlation length', &
            required=.false., act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--ndim_field', switch_ab='-nd', &
            help='the number of dimensions of the MHD field', &
            required=.false., act='store', def='2', error=error)
        if (error/=0) stop
        call cli%add(switch='--drift_param1', switch_ab='-dp1', &
            help='Drift parameter for 3D simulation (ev_ABL_0/pc)', &
            required=.false., def='4E7', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--drift_param2', switch_ab='-dp2', &
            help='Drift parameter for 3D simulation (mev_ABL_0/p^2)', &
            required=.false., def='2E8', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--charge', switch_ab='-ch', &
            help='Particle charge in the unit of e', &
            required=.false., act='store', def='-1', error=error)
        if (error/=0) stop
        call cli%add(switch='--spherical_coord', switch_ab='-sc', &
            help='whether MHD simulations are spherical coordinates', &
            required=.false., act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--uniform_grid', switch_ab='-ug', &
            help='whether the grid is uniform', &
            required=.false., act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--check_drift_2d', switch_ab='-cd', &
            help='whether to check drift in 2D simulations', &
            required=.false., act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--particle_data_dump', switch_ab='-pd', &
            help='whether to dump particle data', &
            required=.false., act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--include_3rd_dim', switch_ab='-i3', &
            help='whether to include 3rd-dim transport in 2D runs', &
            required=.false., act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--acc_by_surface', switch_ab='-as', &
            help='whether the acceleration region is separated by a surface', &
            required=.false., act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--surface_filename1', switch_ab='-sf1', &
            help='the filename for the surface 1 data', required=.false., &
            act='store', def="sname_default", error=error)
        if (error/=0) stop
        call cli%add(switch='--surface_norm1', switch_ab='-sn1', &
            help='the norm direction of surface 1', required=.false., &
            act='store', def='+y', error=error)
        if (error/=0) stop
        call cli%add(switch='--surface2_existed', switch_ab='-s2e', &
            help='whether there is a surface 2', required=.false., &
            act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--is_intersection', switch_ab='-ii', &
            help='whether it is the intersection of the two regions', &
            required=.false., act='store', def='.true.', error=error)
        if (error/=0) stop
        call cli%add(switch='--surface_filename2', switch_ab='-sf2', &
            help='the filename for the surface 2 data', required=.false., &
            act='store', def="sname_default", error=error)
        if (error/=0) stop
        call cli%add(switch='--surface_norm2', switch_ab='-sn2', &
            help='the norm direction of surface 2', required=.false., &
            act='store', def='-y', error=error)
        if (error/=0) stop
        call cli%add(switch='--varying_dt_mhd', switch_ab='-vdt', &
            help='whether the time interval for MHD fields is varying', &
            required=.false., act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--duu_init', switch_ab='-du', &
            help='pitch-angle diffusion for particles with p0', &
            required=.false., def='1E1', act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--dump_escaped_dist', switch_ab='-ded', &
            help='whether to dump the distributions of the escaped particles', &
            required=.false., act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%add(switch='--dump_escaped', switch_ab='-de', &
            help='whether to dump escaped particles', &
            required=.false., act='store', def='.false.', error=error)
        if (error/=0) stop
        call cli%get(switch='-qh', val=quota_hour, error=error)
        if (error/=0) stop
        call cli%get(switch='-rf', val=restart_flag, error=error)
        if (error/=0) stop
        call cli%get(switch='-ft', val=focused_transport, error=error)
        if (error/=0) stop
        call cli%get(switch='-pv', val=particle_v0, error=error)
        if (error/=0) stop
        call cli%get(switch='-sm', val=size_mpi_sub, error=error)
        if (error/=0) stop
        call cli%get(switch='-dm', val=dir_mhd_data, error=error)
        if (error/=0) stop
        call cli%get(switch='-mc', val=mhd_config_filename, error=error)
        if (error/=0) stop
        call cli%get(switch='-nm', val=nptl_max, error=error)
        if (error/=0) stop
        call cli%get(switch='-np', val=nptl, error=error)
        if (error/=0) stop
        call cli%get(switch='-ti', val=time_interp_flag, error=error)
        if (error/=0) stop
        call cli%get(switch='-dt', val=dt, error=error)
        if (error/=0) stop
        call cli%get(switch='-ts', val=t_start, error=error)
        if (error/=0) stop
        call cli%get(switch='-te', val=t_end, error=error)
        if (error/=0) stop
        call cli%get(switch='-tm', val=tmax_mhd, error=error)
        if (error/=0) stop
        call cli%get(switch='-st', val=single_time_frame, error=error)
        if (error/=0) stop
        call cli%get(switch='-df', val=dist_flag, error=error)
        if (error/=0) stop
        call cli%get(switch='-pi', val=power_index, error=error)
        if (error/=0) stop
        call cli%get(switch='-sf', val=split_flag, error=error)
        if (error/=0) stop
        call cli%get(switch='-sr', val=split_ratio, error=error)
        if (error/=0) stop
        call cli%get(switch='-ps', val=pmin_split, error=error)
        if (error/=0) stop
        call cli%get(switch='-tf', val=track_particle_flag, error=error)
        if (error/=0) stop
        call cli%get(switch='-ptf', val=particle_tags_file, error=error)
        if (error/=0) stop
        call cli%get(switch='-ni', val=nsteps_interval, error=error)
        if (error/=0) stop
        call cli%get(switch='-dd', val=diagnostics_directory, error=error)
        if (error/=0) stop
        call cli%get(switch='-is', val=inject_at_shock, error=error)
        if (error/=0) stop
        call cli%get(switch='-in', val=inject_new_ptl, error=error)
        if (error/=0) stop
        call cli%get(switch='-ij', val=inject_large_jz, error=error)
        if (error/=0) stop
        call cli%get(switch='-sn', val=inject_same_nptl, error=error)
        if (error/=0) stop
        call cli%get(switch='-tti', val=tmax_to_inject, error=error)
        if (error/=0) stop
        call cli%get(switch='-ip', val=inject_part_box, error=error)
        if (error/=0) stop
        call cli%get(switch='-jz', val=jz_min, error=error)
        if (error/=0) stop
        call cli%get(switch='-nn', val=ncells_large_jz_norm, error=error)
        if (error/=0) stop
        call cli%get(switch='-ib', val=inject_large_db2, error=error)
        if (error/=0) stop
        call cli%get(switch='-db2', val=db2_min, error=error)
        if (error/=0) stop
        call cli%get(switch='-nb', val=ncells_large_db2_norm, error=error)
        if (error/=0) stop
        call cli%get(switch='-iv', val=inject_large_divv, error=error)
        if (error/=0) stop
        call cli%get(switch='-dv', val=divv_min, error=error)
        if (error/=0) stop
        call cli%get(switch='-nv', val=ncells_large_divv_norm, error=error)
        if (error/=0) stop
        call cli%get(switch='-ir', val=inject_large_rho, error=error)
        if (error/=0) stop
        call cli%get(switch='-rm', val=rho_min, error=error)
        if (error/=0) stop
        call cli%get(switch='-nr', val=ncells_large_rho_norm, error=error)
        if (error/=0) stop
        call cli%get(switch='-iaj', val=inject_large_absj, error=error)
        if (error/=0) stop
        call cli%get(switch='-ajm', val=absj_min, error=error)
        if (error/=0) stop
        call cli%get(switch='-naj', val=ncells_large_absj_norm, error=error)
        if (error/=0) stop
        call cli%get(switch='-xs', val=ptl_xmin, error=error)
        if (error/=0) stop
        call cli%get(switch='-xe', val=ptl_xmax, error=error)
        if (error/=0) stop
        call cli%get(switch='-ys', val=ptl_ymin, error=error)
        if (error/=0) stop
        call cli%get(switch='-ye', val=ptl_ymax, error=error)
        if (error/=0) stop
        call cli%get(switch='-zs', val=ptl_zmin, error=error)
        if (error/=0) stop
        call cli%get(switch='-ze', val=ptl_zmax, error=error)
        if (error/=0) stop
        call cli%get(switch='-cf', val=conf_file, error=error)
        if (error/=0) stop
        call cli%get(switch='-nf', val=num_fine_steps, error=error)
        if (error/=0) stop
        call cli%get(switch='-ld', val=local_dist, error=error)
        if (error/=0) stop
        call cli%get(switch='-dw', val=dpp_wave, error=error)
        if (error/=0) stop
        call cli%get(switch='-ds', val=dpp_shear, error=error)
        if (error/=0) stop
        call cli%get(switch='-ws', val=weak_scattering, error=error)
        if (error/=0) stop
        call cli%get(switch='-t0', val=tau0_scattering, error=error)
        if (error/=0) stop
        call cli%get(switch='-db', val=deltab_flag, error=error)
        if (error/=0) stop
        call cli%get(switch='-co', val=correlation_flag, error=error)
        if (error/=0) stop
        call cli%get(switch='-nd', val=ndim_field, error=error)
        if (error/=0) stop
        call cli%get(switch='-dp1', val=drift_param1, error=error)
        if (error/=0) stop
        call cli%get(switch='-dp2', val=drift_param2, error=error)
        if (error/=0) stop
        call cli%get(switch='-ch', val=charge, error=error)
        if (error/=0) stop
        call cli%get(switch='-sc', val=spherical_coord, error=error)
        if (error/=0) stop
        call cli%get(switch='-ug', val=uniform_grid, error=error)
        if (error/=0) stop
        call cli%get(switch='-cd', val=check_drift_2d, error=error)
        if (error/=0) stop
        call cli%get(switch='-pd', val=particle_data_dump, error=error)
        if (error/=0) stop
        call cli%get(switch='-i3', val=include_3rd_dim, error=error)
        if (error/=0) stop
        call cli%get(switch='-as', val=acc_by_surface, error=error)
        if (error/=0) stop
        call cli%get(switch='-sf1', val=surface_filename1, error=error)
        if (error/=0) stop
        call cli%get(switch='-sn1', val=surface_norm1, error=error)
        if (error/=0) stop
        call cli%get(switch='-s2e', val=surface2_existed, error=error)
        if (error/=0) stop
        call cli%get(switch='-ii', val=is_intersection, error=error)
        if (error/=0) stop
        call cli%get(switch='-sf2', val=surface_filename2, error=error)
        if (error/=0) stop
        call cli%get(switch='-sn2', val=surface_norm2, error=error)
        if (error/=0) stop
        call cli%get(switch='-vdt', val=varying_dt_mhd, error=error)
        if (error/=0) stop
        call cli%get(switch='-du', val=duu0, error=error)
        if (error/=0) stop
        call cli%get(switch='-ded', val=dump_escaped_dist, error=error)
        if (error/=0) stop
        call cli%get(switch='-de', val=dump_escaped, error=error)
        if (error/=0) stop

        if (mpi_rank == master) then
            print '(A)', '---------------Commandline Arguments:---------------'
            if (restart_flag) then
                print '(A)', 'This is a restart of a previous simulation'
            endif
            print '(A,E14.7)', ' The maximum wall time: ', quota_hour
            if (focused_transport) then
                print '(A)', 'Solving the Focused Transport Equation'
                print '(A,E14.7)', ' Initial particle velocity in the normalized velocity: ', &
                    particle_v0
                print '(A,E14.7)', ' Pitch-angle diffusion for particles with a momentum of p0: ', &
                    duu0
            else
                print '(A)', "Solving Parker's Transport Equation"
            endif
            print '(A,I10.3)', 'Size of the a MPI sub-communicator: ', size_mpi_sub
            print '(A,A)', 'Direcotry of MHD data files: ', trim(dir_mhd_data)
            print '(A,A)', 'MHD configuration filename: ', trim(mhd_config_filename)
            print '(A,I15.3)', 'Maximum number of particles: ', nptl_max
            print '(A,I15.3)', 'Initial number of particles: ', nptl
            if (dist_flag == 0) then
                print '(A)', 'Initial distribution is Maxwellian'
            else if (dist_flag == 1) then
                print '(A)', 'All particles with the same momentum initially'
            else
                print '(A,F14.7)', 'Initial distribution is a power-law with index: ', &
                    power_index
            endif
            print '(A,E14.7)', 'Time interval to push particles: ', dt
            print '(A,I10.3)', 'Starting time frame: ', t_start
            print '(A,I10.3)', 'The last time frame: ', t_end
            print '(A,I10.3)', 'The maximum time frame for the MHD fields: ', tmax_mhd
            if (single_time_frame == 1) then
                print '(A)', 'Use only one time of frame of fields data'
                print '(A,I10.3)', 'We will use the data from starting time frame: ', t_start
            else
                print '(A)', 'Use multiple time frames of fields data'
                if (time_interp_flag == 1) then
                    print '(A)', 'Interpolate the fields data inbetween time frames'
                else
                    print '(A)', 'Do not interpolate the fields data inbetween time frames'
                endif
                if (tmax_mhd < t_end) then
                    print '(A)', 'Since tmax_mhd < t_end, the simulation will use the last MHD frame after t_end.'
                endif
            endif

            if (split_flag == 1) then
                print '(A)', 'Particles will be splitted when reaching certain energies'
                print '(A,E14.7)', ' Momentum increase ratio for particle splitting: ', &
                    split_ratio
                print '(A,E14.7)', &
                    ' The minimum momentum (in terms in p0) to start splitting particles: ', &
                    pmin_split
            endif
            if (track_particle_flag) then
                print '(A)', 'The program will tracking high-energy particles'
                print '(A,A)', ' * File containing the selected particle tags: ', particle_tags_file
                print '(A,I10.3)', ' * Steps interval to track particles: ', nsteps_interval
            endif
            print '(A,A)', 'Diagnostic file directory is: ', trim(diagnostics_directory)
            print '(A,I10.3)', 'Number of fine steps: ', num_fine_steps
            if (inject_at_shock) then
                print '(A)', 'Inject particles at shock location'
            endif
            if (inject_new_ptl) then
                print '(A)', 'Inject new particles at every MHD time step'
                if (tmax_to_inject < t_end) then
                    print '(A, I0)', 'The simulation will stop injecting particles after frame', &
                        tmax_to_inject
                endif
            endif
            if (inject_part_box) then
                print '(A)', 'Inject new particles in part of the simulation box'
                print '(A)', 'The local region to inject new particles (in code unit):'
                print '(A,E14.7E2,E14.7E2)', ' Min and Max x: ', ptl_xmin, ptl_xmax
                print '(A,E14.7E2,E14.7E2)', ' Min and Max y: ', ptl_ymin, ptl_ymax
                print '(A,E14.7E2,E14.7E2)', ' Min and Max z: ', ptl_zmin, ptl_zmax
            endif
            if (inject_large_jz) then
                print '(A)', 'Inject new particles where jz is large'
                if (inject_same_nptl) then
                    print '(A)', 'Inject the same number of particles every step'
                else
                    print '(A, A)', 'Inject different number of particles every step, ', &
                        'depending on the number of cells with |jz| > jz_min'
                    print '(A,I10.3)', 'Normalization for the number of cells with large jz: ', &
                        ncells_large_jz_norm
                endif
                print '(A,E14.7E2)', 'Minumum jz in regions to inject new particles', jz_min
            endif
            if (inject_large_db2) then
                print '(A)', 'Inject new particles where db2 is large'
                if (inject_same_nptl) then
                    print '(A)', 'Inject the same number of particles every step'
                else
                    print '(A, A)', 'Inject different number of particles every step, ', &
                        'depending on the number of cells with db2 > db2_min'
                    print '(A,I10.3)', 'Normalization for the number of cells with large db2: ', &
                        ncells_large_db2_norm
                endif
                print '(A,E14.7E2)', 'Minumum db2 in regions to inject new particles', db2_min
            endif
            if (inject_large_divv) then
                print '(A)', 'Inject new particles where divv is negatively large'
                if (inject_same_nptl) then
                    print '(A)', 'Inject the same number of particles every step'
                else
                    print '(A, A)', 'Inject different number of particles every step, ', &
                        'depending on the number of cells with -divv > divv_min'
                    print '(A,I10.3)', 'Normalization for the number of cells with large divv: ', &
                        ncells_large_divv_norm
                endif
                print '(A,E14.7E2)', 'Minumum divv in regions to inject new particles', divv_min
            endif
            if (inject_large_rho) then
                print '(A)', 'Inject new particles where rho large'
                if (inject_same_nptl) then
                    print '(A)', 'Inject the same number of particles every step'
                else
                    print '(A, A)', 'Inject different number of particles every step, ', &
                        'depending on the number of cells with rho > rho_min'
                    print '(A,I10.3)', 'Normalization for the number of cells with large rho: ', &
                        ncells_large_rho_norm
                endif
                print '(A,E14.7E2)', 'Minumum rho in regions to inject new particles', rho_min
            endif
            if (inject_large_absj) then
                print '(A)', 'Inject new particles where absj large'
                if (inject_same_nptl) then
                    print '(A)', 'Inject the same number of particles every step'
                else
                    print '(A, A)', 'Inject different number of particles every step, ', &
                        'depending on the number of cells with absj > absj_min'
                    print '(A,I10.3)', 'Normalization for the number of cells with large absj: ', &
                        ncells_large_absj_norm
                endif
                print '(A,E14.7E2)', 'Minumum absj in regions to inject new particles', absj_min
            endif
            if (dpp_wave == 1) then
                print '(A)', 'Include momentum diffusion due to wave scattering'
            endif
            if (dpp_shear == 1) then
                print '(A)', 'Include momentum diffusion due to flow shear'
                print '(A,E14.7E2)', 'The scattering time for initial particles: ', tau0_scattering
            endif
            if (weak_scattering == 1) then
                print '(A)', 'Particle transport is in the weak-scattering regime (tau\Omega >> 1)'
            else
                print '(A)', 'Particle transport is in the strong-scattering regime (tau\Omega << 1)'
            endif
            print '(A,A)', 'Configuration file name: ', trim(conf_file)
            if (local_dist) then
                print '(A)', 'Simulation will diagnose local particle distribution'
            endif
            if (dump_escaped_dist) then
                print '(A)', 'Simulation will dump the distribution of the escaped particles for open BC'
                if (dump_escaped) then
                    print '(A)', 'Simulation will dump escaped particles for open BC'
                endif
            endif
            if (deltab_flag == 1) then
                print '(A)', 'Including spatially dependent magnetic fluctuation'
            endif
            if (correlation_flag == 1) then
                print '(A)', 'Include spatially dependent turbulence correlation length'
            endif
            print '(A,I10.3)', 'MHD field dimension:', ndim_field
            if (ndim_field == 3 .or. (ndim_field == 2 .and. check_drift_2d == 1)) then
                print '(A,E14.7E2,E14.7E2)', 'Parameters for particle drift: ', &
                    drift_param1, drift_param2
                print '(A,I10.3)', 'Particle change in the unit of e: ', charge
            endif
            if (spherical_coord == 1) then
                print '(A)', 'MHD simulations are in spherical coordinates'
            endif
            if (uniform_grid == 1) then
                print '(A)', 'The grid is uniform'
            else
                print '(A)', 'The grid is non-uniform'
            endif
            if (ndim_field == 2 .and. check_drift_2d == 1) then
                print '(A)', 'Check particle drift in 2D simulations'
            endif
            if (particle_data_dump == 1) then
                print '(A)', 'Particle data will be dumped'
            endif
            if (ndim_field == 2 .and. include_3rd_dim == 1) then
                print '(A)', 'Include transport along the 3rd-dim in 2D simulations'
            endif
            if (acc_by_surface == 1) then
                if (surface2_existed) then
                    print '(A)', 'The acceleration region is separated by two surfaces'
                    if (is_intersection) then
                        print '(A)', 'The acceleration is the intersection of the two regions'
                    else
                        print '(A)', 'The acceleration is the union of the two regions'
                    endif
                    print '(A,A)', 'The filename of the surface 1 data is ', surface_filename1
                    print '(A,A)', 'The norm of surface 1 is along ', surface_norm1
                    print '(A,A)', 'The filename of the surface 2 data is ', surface_filename2
                    print '(A,A)', 'The norm of surface 2 is along ', surface_norm2
                else
                    print '(A)', 'The acceleration region is separated by a surface'
                    print '(A,A)', 'The filename of the surface data is ', surface_filename1
                    print '(A,A)', 'The surface norm is along ', surface_norm1
                endif
            endif
            if (varying_dt_mhd) then
                print '(A)', 'The time interval for MHD fields is varying'
            else
                print '(A)', 'The time interval for MHD fields is fixed'
            endif
            print '(A)', '----------------------------------------------------'
        endif
    end subroutine get_cmd_args
end program stochastic

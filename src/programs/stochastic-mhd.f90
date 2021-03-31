!*******************************************************************************
!< Program for stochastic integration
!*******************************************************************************

program stochastic
    use constants, only: fp, dp
    use mpi_module
    use omp_lib
    use mhd_config_module, only: load_mhd_config, mhd_config, &
        echo_mhd_config, set_mhd_grid_type
    use particle_module, only: init_particles, free_particles, &
        inject_particles_spatial_uniform, read_particle_params, &
        particle_mover, split_particle, &
        init_particle_distributions, free_particle_distributions, &
        distributions_diagnostics, quick_check, set_particle_datatype_mpi, &
        free_particle_datatype_mpi, select_particles_tracking, &
        init_particle_tracking, free_particle_tracking, &
        init_tracked_particle_points, free_tracked_particle_points, &
        negative_particle_tags, save_tracked_particle_points, &
        inject_particles_at_shock, set_mpi_io_data_sizes, &
        init_local_particle_distributions, free_local_particle_distributions, &
        inject_particles_at_large_jz, set_dpp_params, set_flags_params, &
        set_drift_parameters, get_pmax_global, set_flag_check_drift_2d, &
        dump_particles, init_escaped_particles, free_escaped_particles, &
        reset_escaped_particles, dump_escaped_particles, &
        record_tracked_particle_init
    use random_number_generator, only: init_prng, delete_prng
    use mhd_data_parallel, only: init_field_data, free_field_data, &
        read_field_data_parallel, init_fields_gradients, free_fields_gradients, &
        calc_fields_gradients, calc_fields_gradients_nonuniform, &
        copy_fields, init_shock_xpos, free_shock_xpos, &
        locate_shock_xpos, init_magnetic_fluctuation, free_magnetic_fluctuation, &
        read_magnetic_fluctuation, copy_magnetic_fluctuation, &
        init_correlation_length, free_correlation_length, &
        read_correlation_length, copy_correlation_length, &
        init_grad_deltab2, free_grad_deltab2, &
        calc_grad_deltab2, calc_grad_deltab2_nonuniform, &
        init_grad_correl_length, free_grad_correl_length, &
        calc_grad_correl_length, calc_grad_correl_length_nonuniform, &
        init_grid_positions, free_grid_positions, &
        set_local_grid_positions
    use simulation_setup_module, only: read_simuation_mpi_topology, &
        set_field_configuration, fconfig, read_particle_boundary_conditions, &
        set_neighbors, check_particle_can_escape
    use flap, only : command_line_interface !< FLAP package
    use penf
    implicit none
    character(len=256) :: dir_mhd_data
    character(len=256) :: fname1, fname2
    character(len=128) :: diagnostics_directory
    character(len=32) :: conf_file
    integer :: nptl_max, nptl
    real(dp) :: start, finish, step1, step2, dt, jz_min
    real(dp) :: ptl_xmin, ptl_xmax, ptl_ymin, ptl_ymax, ptl_zmin, ptl_zmax
    real(dp) :: tau0_scattering ! Scattering time for initial particles
    real(dp) :: drift_param1, drift_param2 ! Drift parameter for 3D simulation
    real(dp), dimension(6) :: part_box
    integer :: nthreads, color
    integer :: t_start, t_end, tf, dist_flag, split_flag
    integer :: interp_flag, local_dist
    integer :: track_particle_flag, nptl_selected, nsteps_interval
    integer :: inject_at_shock    ! Inject particles at shock location
    integer :: num_fine_steps     ! Number of the fine steps for diagnostics
    integer :: inject_new_ptl     ! Whether to inject new particles every MHD step
    integer :: inject_large_jz    ! Whether to inject where jz is large
    integer :: inject_part_box    ! Whether to inject in part of the box
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

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)

    start = MPI_Wtime()

    call get_cmd_args

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
    call init_prng(mpi_rank, nthreads)

    !< Configurations
    write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'mhd_config.dat'
    call load_mhd_config(fname2)
    if (mpi_rank == master) then
        call echo_mhd_config
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
    interp_flag = 1 ! Two time are needed for interpolation
    call init_field_data(interp_flag)
    call init_fields_gradients(interp_flag)
    if (deltab_flag == 1) then
        call init_magnetic_fluctuation(interp_flag)
        call init_grad_deltab2(interp_flag)
    endif
    if (correlation_flag == 1) then
        call init_correlation_length(interp_flag)
        call init_grad_correl_length(interp_flag)
    endif
    call read_particle_params(conf_file)
    call set_dpp_params(dpp_wave, dpp_shear, weak_scattering, tau0_scattering)
    call set_flags_params(deltab_flag, correlation_flag)
    if (ndim_field == 3) then
        call set_drift_parameters(drift_param1, drift_param2, charge)
    else if (ndim_field == 2 .and. check_drift_2d == 1) then
        call set_drift_parameters(drift_param1, drift_param2, charge)
        call set_flag_check_drift_2d(check_drift_2d)
    endif

    call init_particles(nptl_max)
    call init_escaped_particles
    call init_particle_distributions
    if (local_dist == 1) then
        call init_local_particle_distributions
    endif
    call set_particle_datatype_mpi

    call solve_transport_equation(.false.)

    if (track_particle_flag == 1) then
        !< We need to reset the random number generator
        call delete_prng
        call init_prng(mpi_rank, nthreads)

        call init_particle_tracking(nptl_selected)
        call select_particles_tracking(nptl, nptl_selected, nsteps_interval)
        call init_tracked_particle_points(nptl_selected)

        call solve_transport_equation(.true.)

        call save_tracked_particle_points(nptl_selected, diagnostics_directory)
        call free_tracked_particle_points
        call free_particle_tracking
        call delete_prng
    endif

    if (inject_at_shock == 1) then
        call free_shock_xpos(interp_flag)
    endif

    call free_particle_datatype_mpi
    call free_particle_distributions
    if (local_dist == 1) then
        call free_local_particle_distributions
    endif
    call free_particles
    call free_escaped_particles
    call free_fields_gradients(interp_flag)
    if (deltab_flag == 1) then
        call free_grad_deltab2(interp_flag)
    endif
    if (correlation_flag == 1) then
        call free_grad_correl_length(interp_flag)
    endif
    call free_field_data(interp_flag)
    if (deltab_flag == 1) then
        call free_magnetic_fluctuation(interp_flag)
    endif
    if (correlation_flag == 1) then
        call free_correlation_length(interp_flag)
    endif

    call free_grid_positions

    finish = MPI_Wtime()
    if (mpi_rank == master) then
        print '("Time = ",f9.4," seconds.")',finish-start
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
        write(fname1, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', t_start
        write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', t_start + 1
        call read_field_data_parallel(fname1, var_flag=0)
        call read_field_data_parallel(fname2, var_flag=1)
        if (deltab_flag == 1) then
            write(fname1, "(A,I4.4)") trim(dir_mhd_data)//'deltab_', t_start
            write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'deltab_', t_start + 1
            call read_magnetic_fluctuation(fname1, var_flag=0)
            call read_magnetic_fluctuation(fname2, var_flag=1)
        endif
        if (correlation_flag == 1) then
            write(fname1, "(A,I4.4)") trim(dir_mhd_data)//'lc_', t_start
            write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'lc_', t_start + 1
            call read_correlation_length(fname1, var_flag=0)
            call read_correlation_length(fname2, var_flag=1)
        endif

        if (uniform_grid == 1) then
            call calc_fields_gradients(var_flag=0)
            call calc_fields_gradients(var_flag=1)
        else
            call calc_fields_gradients_nonuniform(var_flag=0)
            call calc_fields_gradients_nonuniform(var_flag=1)
        endif
        if (deltab_flag == 1) then
            if (uniform_grid == 1) then
                call calc_grad_deltab2(var_flag=0)
                call calc_grad_deltab2(var_flag=1)
            else
                call calc_grad_deltab2_nonuniform(var_flag=0)
                call calc_grad_deltab2_nonuniform(var_flag=1)
            endif
        endif
        if (correlation_flag == 1) then
            if (uniform_grid == 1) then
                call calc_grad_correl_length(var_flag=0)
                call calc_grad_correl_length(var_flag=1)
            else
                call calc_grad_correl_length_nonuniform(var_flag=0)
                call calc_grad_correl_length_nonuniform(var_flag=1)
            endif
        endif

        if (inject_at_shock == 1) then ! Whether it is a shock problem
            call init_shock_xpos(interp_flag)
            call locate_shock_xpos(interp_flag)
            call inject_particles_at_shock(nptl, dt, dist_flag, t_start)
        else
            if (inject_part_box == 1) then
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
            if (inject_large_jz == 1) then
                call inject_particles_at_large_jz(nptl, dt, dist_flag, t_start, &
                    jz_min, part_box)
            else
                call inject_particles_spatial_uniform(nptl, dt, dist_flag, t_start, part_box)
            endif
        endif

        if (track_particles) then
            call negative_particle_tags(nptl_selected)
            call record_tracked_particle_init(nptl_selected)
        endif

        step1 = MPI_Wtime()

        if (.not. track_particles) then
            call set_mpi_io_data_sizes
            call distributions_diagnostics(t_start, diagnostics_directory, local_dist)
            if (particle_data_dump == 1) then
                call dump_particles(t_start, diagnostics_directory)
            endif
            call quick_check(t_start, .true., diagnostics_directory)
            call get_pmax_global(t_start, .true., diagnostics_directory)
        endif

        !< Time loop
        do tf = t_start, t_end
            if (mpi_rank == master) then
                write(*, "(A,I0)") " Starting step ", tf
            endif
            if (track_particles) then
                call particle_mover(1, nptl_selected, nsteps_interval, tf, 1)
            else
                call particle_mover(0, nptl_selected, nsteps_interval, tf, num_fine_steps)
            endif
            if (mpi_rank == master) then
                write(*, "(A)") " Finishing moving particles "
            endif
            if (split_flag == 1) then
                call split_particle
            endif
            if (.not. track_particles) then
                call quick_check(tf+1, .false., diagnostics_directory)
                call get_pmax_global(tf+1, .false., diagnostics_directory)
                call distributions_diagnostics(tf+1, diagnostics_directory, local_dist)
                if (mpi_rank == master) then
                    write(*, "(A)") " Finishing distribution diagnostics "
                endif
                if (particle_data_dump == 1) then
                    call dump_particles(tf+1, diagnostics_directory)
                endif
                call dump_escaped_particles(tf+1, diagnostics_directory)
                call reset_escaped_particles
            endif
            call copy_fields
            if (deltab_flag == 1) then
                call copy_magnetic_fluctuation
            endif
            if (correlation_flag == 1) then
                call copy_correlation_length
            endif
            if (mpi_rank == master) then
                write(*, "(A)") " Finishing copying fields "
            endif
            if (tf < t_end) then
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', tf + 1
                call read_field_data_parallel(fname2, var_flag=1)
                if (deltab_flag == 1) then
                    write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'deltab_', tf + 1
                    call read_magnetic_fluctuation(fname2, var_flag=1)
                endif
                if (correlation_flag == 1) then
                    write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'lc_', tf + 1
                    call read_correlation_length(fname2, var_flag=1)
                endif
                if (uniform_grid == 1) then
                    call calc_fields_gradients(var_flag=1)
                else
                    call calc_fields_gradients_nonuniform(var_flag=1)
                endif
                if (deltab_flag == 1) then
                    if (uniform_grid == 1) then
                        call calc_grad_deltab2(var_flag=1)
                    else
                        call calc_grad_deltab2_nonuniform(var_flag=1)
                    endif
                endif
                if (correlation_flag == 1) then
                    if (uniform_grid == 1) then
                        call calc_grad_correl_length(var_flag=1)
                    else
                        call calc_grad_correl_length_nonuniform(var_flag=1)
                    endif
                endif
            endif
            if (inject_at_shock == 1) then
                call locate_shock_xpos(interp_flag)
                call inject_particles_at_shock(nptl, dt, dist_flag, tf+1)
            else
                if (inject_new_ptl == 1) then
                    if (inject_large_jz == 1) then
                        call inject_particles_at_large_jz(nptl, dt, dist_flag, tf+1, &
                            jz_min, part_box)
                    else
                        call inject_particles_spatial_uniform(nptl, dt, dist_flag, tf+1, part_box)
                    endif
                endif
            endif
            if (track_particles) then
                call negative_particle_tags(nptl_selected)
                call record_tracked_particle_init(nptl_selected)
            endif
            step2 = MPI_Wtime()
            if (mpi_rank == master) then
                print '("Step ", I0, " takes ", f9.4, " seconds.")', tf, step2 - step1
            endif
            step1 = step2
        enddo  ! Time loop
    end subroutine solve_transport_equation

    !---------------------------------------------------------------------------
    !< Get commandline arguments
    !---------------------------------------------------------------------------
    subroutine get_cmd_args
        implicit none
        type(command_line_interface) :: cli     !< Command Line Interface (CLI).
        integer(I4P)                 :: error   !< Error trapping flag.
        integer :: tstart
        call cli%init(&
            progname = 'stochastic-mhd', &
            authors     = 'Xiaocan Li', &
            help        = 'Usage: ', &
            description = "Solving Parker's transport equation "// &
                          "using stochastic differential equation", &
            examples = ['stochastic-mhd.exec -sm size_mpi_sub -dm dir_mhd_data '//&
                        '-nm nptl_max -np nptl -dt dt -ts t_start -te t_end '//&
                        '-df dist_flag -tf track_particle_flag '//&
                        '-ns nptl_selected -ni nsteps_interval '//&
                        '-dd diagnostics_directory -is inject_at_shock '//&
                        '-in inject_new_ptl -ij inject_at_shock '//&
                        '-ip inject_part_box -xs ptl_xmin -xe ptl_xmax '//&
                        '-ys ptl_ymin -ye ptl_ymax -zs ptl_zmin -ze ptl_zmax '//&
                        '-cf conf_file -nf num_fine_steps -dw dpp_wave '//&
                        '-ds dpp_shear -t0 tau0_scattering '//&
                        '-db deltab_flag -co correlation_flag '//&
                        '-dp1 drift_param1 -dp2 drift_param2 -ch charge '//&
                        '-sc spherical_coord -ug uniform_grid '//&
                        '-cd check_drift_2d -pd particle_data_dump'])
        call cli%add(switch='--size_mpi_sub', switch_ab='-sm', &
            help='Size of the a MPI sub-communicator', required=.false., &
            act='store', def='1', error=error)
        if (error/=0) stop
        call cli%add(switch='--dir_mhd_data', switch_ab='-dm', &
            help='MHD simulation data file directory', required=.true., &
            act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--nptl_max', switch_ab='-nm', &
            help='Maximum number of particles', required=.false., &
            act='store', def='1E7', error=error)
        if (error/=0) stop
        call cli%add(switch='--nptl', switch_ab='-np', &
            help='Number of particles', required=.false., &
            act='store', def='1E4', error=error)
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
        call cli%add(switch='--dist_flag', switch_ab='-df', &
            help='Flag to indicate momentum distribution', required=.false., &
            act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--split_flag', switch_ab='-sf', &
            help='Flag to split the particles', required=.false., &
            act='store', def='1', error=error)
        if (error/=0) stop
        call cli%add(switch='--track_particle_flag', switch_ab='-tf', &
            help='Flag to track some particles', required=.false., &
            act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--nptl_selected', switch_ab='-ns', &
            help='Number of selected particles to track', required=.false., &
            act='store', def='100', error=error)
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
            act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--inject_new_ptl', switch_ab='-in', &
            help='whether to inject new particles every MHD step', required=.false., &
            act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--inject_large_jz', switch_ab='-ij', &
            help='whether to inject where jz is large', required=.false., &
            act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--inject_part_box', switch_ab='-ip', &
            help='whether to inject in a part of the box', required=.false., &
            act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--jz_min', switch_ab='-jz', &
            help='Minimum jz in regions to inject new particles', &
            required=.false., def='100.0', act='store', error=error)
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
            help='whether to diagnose local distribution', required=.false., &
            act='store', def='0', error=error)
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
        call cli%get(switch='-sm', val=size_mpi_sub, error=error)
        if (error/=0) stop
        call cli%get(switch='-dm', val=dir_mhd_data, error=error)
        if (error/=0) stop
        call cli%get(switch='-nm', val=nptl_max, error=error)
        if (error/=0) stop
        call cli%get(switch='-np', val=nptl, error=error)
        if (error/=0) stop
        call cli%get(switch='-dt', val=dt, error=error)
        if (error/=0) stop
        call cli%get(switch='-ts', val=t_start, error=error)
        if (error/=0) stop
        call cli%get(switch='-te', val=t_end, error=error)
        if (error/=0) stop
        call cli%get(switch='-df', val=dist_flag, error=error)
        if (error/=0) stop
        call cli%get(switch='-sf', val=split_flag, error=error)
        if (error/=0) stop
        call cli%get(switch='-tf', val=track_particle_flag, error=error)
        if (error/=0) stop
        call cli%get(switch='-ns', val=nptl_selected, error=error)
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
        call cli%get(switch='-ip', val=inject_part_box, error=error)
        if (error/=0) stop
        call cli%get(switch='-jz', val=jz_min, error=error)
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

        if (mpi_rank == master) then
            print '(A,I10.3)', 'Size of the a MPI sub-communicator: ', size_mpi_sub
            print '(A,A)', 'Direcotry of MHD data files: ', trim(dir_mhd_data)
            print '(A,I10.3)', 'Maximum number of particles: ', nptl_max
            print '(A,I10.3)', 'Initial number of particles: ', nptl
            if (dist_flag == 0) then
                print '(A)', 'Initial distribution in Maxwellian'
            else
                print '(A)', 'All particles with the same momentum initially'
            endif
            print '(A,E14.7E2)', 'Time interval to push particles', dt
            print '(A,I10.3)', 'Starting time frame: ', t_start
            print '(A,I10.3)', 'The last time frame: ', t_end
            if (split_flag == 1) then
                print '(A)', 'Particles will be splitted when reaching certain energies'
            endif
            if (track_particle_flag == 1) then
                print '(A,I10.3)', 'Number of particles to track: ', nptl_selected
                print '(A,I10.3)', 'Steps interval to track particles: ', nsteps_interval
            endif
            print '(A,A)', 'Diagnostic file directory is: ', trim(diagnostics_directory)
            print '(A,I10.3)', 'Number of fine steps: ', num_fine_steps
            if (inject_at_shock == 1) then
                print '(A)', 'Inject particles at shock location'
            endif
            if (inject_new_ptl == 1) then
                print '(A)', 'Inject new particles at every MHD time step'
            endif
            if (inject_part_box == 1) then
                print '(A)', 'Inject new particles in part of the simulation box'
                print '(A,E14.7E2,E14.7E2)', 'Minimum and maximum x of the region to inject new particles', &
                    ptl_xmin, ptl_xmax
                print '(A,E14.7E2,E14.7E2)', 'Minimum and maximum y of the region to inject new particles', &
                    ptl_ymin, ptl_ymax
                print '(A,E14.7E2,E14.7E2)', 'Minimum and maximum z of the region to inject new particles', &
                    ptl_zmin, ptl_zmax
            endif
            if (inject_large_jz == 1) then
                print '(A)', 'Inject new particles where jz is large'
                print '(A,E14.7E2)', 'Minumum jz in regions to inject new particles', jz_min
            endif
            if (dpp_wave == 1) then
                print '(A)', 'Include momentum diffusion due to wave scattering'
            endif
            if (dpp_shear == 1) then
                print '(A)', 'Include momentum diffusion due to flow shear'
                print '(A,E14.7E2)', 'The scattering time for initial particles', tau0_scattering
            endif
            if (weak_scattering == 1) then
                print '(A)', 'Particle transport is in the weak-scattering regime (tau\Omega >> 1)'
            else
                print '(A)', 'Particle transport is in the strong-scattering regime (tau\Omega << 1)'
            endif
            print '(A,A)', 'Configuration file name: ', trim(conf_file)
            if (local_dist == 1) then
                print '(A)', 'Simulation will diagnose local particle distribution'
            endif
            if (deltab_flag == 1) then
                print '(A)', 'Including spatially dependent magnetic fluctuation'
            endif
            if (correlation_flag == 1) then
                print '(A)', 'Include spatially dependent turbulence correlation length'
            endif
            print '(A,I10.3)', 'MHD field dimension:', ndim_field
            if (ndim_field == 3 .or. (ndim_field == 2 .and. check_drift_2d == 1)) then
                print '(A,E14.7E2,E14.7E2)', 'Parameters for particle drift', &
                    drift_param1, drift_param2
                print '(A,I10.3)', 'Particle change in the unit of e:', charge
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
        endif
    end subroutine get_cmd_args
end program stochastic

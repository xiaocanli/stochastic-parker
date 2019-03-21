!*******************************************************************************
!< Program for stochastic integration
!*******************************************************************************

program stochastic
    use constants, only: fp, dp
    use mpi_module
    use mhd_config_module, only: load_mhd_config, mhd_config, echo_mhd_config
    use particle_module, only: init_particles, free_particles, &
        inject_particles_spatial_uniform, read_particle_params, &
        particle_mover, remove_particles, split_particle, &
        init_particle_distributions, free_particle_distributions, &
        distributions_diagnostics, quick_check, set_particle_datatype_mpi, &
        free_particle_datatype_mpi, select_particles_tracking, &
        init_particle_tracking, free_particle_tracking, &
        init_tracked_particle_points, free_tracked_particle_points, &
        negative_particle_tags, save_tracked_particle_points, &
        inject_particles_at_shock, set_mpi_io_data_sizes, &
        init_local_particle_distributions, free_local_particle_distributions
    use random_number_generator, only: init_prng, delete_prng
    use mhd_data_parallel, only: init_field_data, free_field_data, &
        read_field_data_parallel, init_fields_gradients, free_fields_gradients, &
        calc_fields_gradients, copy_fields, init_shock_xpos, free_shock_xpos, &
        locate_shock_xpos
    use simulation_setup_module, only: read_simuation_mpi_topology, &
        set_field_configuration, fconfig, read_particle_boundary_conditions, &
        set_neighbors
    use flap, only : command_line_interface !< FLAP package
    use penf
    implicit none
    character(len=256) :: dir_mhd_data
    character(len=256) :: fname1, fname2
    character(len=128) :: diagnostics_directory
    character(len=32) :: conf_file
    integer :: nptl_max, nptl
    real(dp) :: start, finish, step1, step2, dt
    integer :: t_start, t_end, tf, dist_flag, split_flag, whole_mhd_data
    integer :: interp_flag, nx, ny, nz, local_dist
    integer :: track_particle_flag, nptl_selected, nsteps_interval
    integer :: inject_at_shock  ! Inject particles at shock location
    integer :: num_fine_steps   ! Number of the fine steps for diagnostics

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)

    call cpu_time(start)

    call get_cmd_args
    call init_prng

    !< Configurations
    write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'mhd_config.dat'
    call load_mhd_config(fname2)
    if (mpi_rank == master) then
        call echo_mhd_config
    endif
    call read_simuation_mpi_topology(conf_file)
    call read_particle_boundary_conditions(conf_file)
    call set_field_configuration(whole_data_flag=whole_mhd_data, ndim=2)
    call set_neighbors(whole_data_flag=whole_mhd_data)

    !< Initialization
    interp_flag = 1 ! Two time are needed for interpolation
    nx = fconfig%nx
    ny = fconfig%ny
    nz = fconfig%nz
    call init_field_data(interp_flag, nx, ny, nz, ndim=2)
    call init_fields_gradients(interp_flag, nx, ny, nz, ndim=2)
    call read_particle_params(conf_file)

    call init_particles(nptl_max)
    call init_particle_distributions
    if (local_dist) then
        call init_local_particle_distributions
    endif
    call set_particle_datatype_mpi

    write(fname1, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', t_start
    write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', t_start + 1
    call read_field_data_parallel(fname1, var_flag=0)
    call read_field_data_parallel(fname2, var_flag=1)

    call calc_fields_gradients(var_flag=0)
    call calc_fields_gradients(var_flag=1)

    ! Whether it is a shock problem
    if (inject_at_shock == 1) then
        call init_shock_xpos(interp_flag, nx, ny, nz, ndim=2)
        call locate_shock_xpos(interp_flag, nx, ny, nz, ndim=2)
        call inject_particles_at_shock(nptl, dt, dist_flag, t_start)
    else
        call inject_particles_spatial_uniform(nptl, dt, dist_flag, t_start)
    endif

    call cpu_time(step1)

    call set_mpi_io_data_sizes
    call distributions_diagnostics(t_start, diagnostics_directory, whole_mhd_data, local_dist)

    !< Time loop
    do tf = t_start, t_end
        if (mpi_rank == master) then
            write(*, "(A,I0)") " Starting step ", tf
        endif
        call particle_mover(0, nptl_selected, nsteps_interval, num_fine_steps)
        if (mpi_rank == master) then
            write(*, "(A)") " Finishing moving particles "
        endif
        if (split_flag == 1) then
            call split_particle
        endif
        if (tf == t_start) then
            call quick_check(tf, .true., diagnostics_directory)
        else
            call quick_check(tf, .false., diagnostics_directory)
        endif
        call distributions_diagnostics(tf + 1, diagnostics_directory, whole_mhd_data, local_dist)
        if (mpi_rank == master) then
            write(*, "(A)") " Finishing distribution diagnostics "
        endif
        call copy_fields
        if (mpi_rank == master) then
            write(*, "(A)") " Finishing copying fields "
        endif
        if (tf < t_end) then
            write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', tf + 1
            call read_field_data_parallel(fname2, var_flag=1)
            call calc_fields_gradients(var_flag=1)
        endif
        if (inject_at_shock == 1) then
            call locate_shock_xpos(interp_flag, nx, ny, nz, ndim=2)
            call inject_particles_at_shock(nptl, dt, dist_flag, tf)
        endif
        call cpu_time(step2)
        if (mpi_rank == master) then
            print '("Step ", I0, " takes ", f9.4, " seconds.")', tf, step2 - step1
        endif
        step1 = step2
    enddo

    if (track_particle_flag == 1) then
        !< We need to reset the random number generator
        call delete_prng
        call init_prng

        call init_particle_tracking(nptl_selected)
        call select_particles_tracking(nptl, nptl_selected, nsteps_interval)

        ! It uses the initial part of the random number series
        if (inject_at_shock == 1) then
            call locate_shock_xpos(interp_flag, nx, ny, nz, ndim=2)
            call inject_particles_at_shock(nptl, dt, dist_flag, t_start)
        else
            call inject_particles_spatial_uniform(nptl, dt, dist_flag, t_start)
        endif

        call init_tracked_particle_points(nptl_selected)
        call negative_particle_tags(nptl_selected)

        write(fname1, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', t_start
        write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', t_start + 1
        call read_field_data_parallel(fname1, var_flag=0)
        call read_field_data_parallel(fname2, var_flag=1)

        call calc_fields_gradients(var_flag=0)
        call calc_fields_gradients(var_flag=1)
        !< Time loop
        do tf = t_start, t_end
            call particle_mover(1, nptl_selected, nsteps_interval, 1)
            if (split_flag == 1) then
                call split_particle
            endif
            call copy_fields
            if (tf < t_end) then
                write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', tf + 1
                call read_field_data_parallel(fname2, var_flag=1)
                call calc_fields_gradients(var_flag=1)
            endif
            if (inject_at_shock == 1) then
                call locate_shock_xpos(interp_flag, nx, ny, nz, ndim=2)
                call inject_particles_at_shock(nptl, dt, dist_flag, tf)
            endif
            call cpu_time(step2)
            if (mpi_rank == master) then
                print '("Step ", I0, " takes ", f9.4, " seconds.")', tf, step2 - step1
            endif
            step1 = step2
        enddo
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
    if (local_dist) then
        call free_local_particle_distributions
    endif
    call free_particles
    call free_fields_gradients(interp_flag)
    call free_field_data(interp_flag)

    call cpu_time(finish)
    if (mpi_rank == master) then
        print '("Time = ",f9.4," seconds.")',finish-start
    endif

    call delete_prng

    call MPI_FINALIZE(ierr)

    contains

    !---------------------------------------------------------------------------
    !< Get commandline arguments
    !---------------------------------------------------------------------------
    subroutine get_cmd_args
        implicit none
        type(command_line_interface) :: cli     !< Command Line Interface (CLI).
        integer(I4P)                 :: error   !< Error trapping flag.
        integer :: tstart
        call cli%init(&
            progname = 'stochastic-2dmhd', &
            authors     = 'Xiaocan Li', &
            help        = 'Usage: ', &
            description = "Solving Parker's transport equation "// &
                          "using stochastic differential equation", &
            examples = ['stochastic-2dmhd.exec -dm dir_mhd_data -nm nptl_max '//&
                        '-np nptl -dt dt -ts t_start -te t_end -df dist_flag '//&
                        '-wm whole_mhd_data -tf track_particle_flag '//&
                        '-ns nptl_selected -ni nsteps_interval '//&
                        '-dd diagnostics_directory -is inject_at_shock '//&
                        '-cf conf_file -nf num_fine_steps'])
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
            act='store', def='1E-6', error=error)
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
        call cli%add(switch='--whole_mhd_data', switch_ab='-wm', &
            help='whether to read the whole MHD data', required=.false., &
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
        call cli%get(switch='-wm', val=whole_mhd_data, error=error)
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
        call cli%get(switch='-cf', val=conf_file, error=error)
        if (error/=0) stop
        call cli%get(switch='-nf', val=num_fine_steps, error=error)
        if (error/=0) stop
        call cli%get(switch='-ld', val=local_dist, error=error)
        if (error/=0) stop

        if (mpi_rank == master) then
            print '(A,A)', 'Direcotry of MHD data files: ', trim(dir_mhd_data)
            print '(A,I0)', 'Maximum number of particles: ', nptl_max
            print '(A,I0)', 'Initial number of particles: ', nptl
            if (dist_flag == 0) then
                print '(A)', 'Initial distribution in Maxwellian'
            else
                print '(A)', 'All particles with the same momentum initially'
            endif
            print '(A,E13.6E2)', 'Time interval to push particles', dt
            print '(A,I0)', 'Starting time frame: ', t_start
            print '(A,I0)', 'The last time frame: ', t_end
            if (split_flag == 1) then
                print '(A)', 'Particles will be splitted when reaching certain energies'
            endif
            if (whole_mhd_data) then
                print '(A)', 'Each process reads the whole MHD simulation data'
            else
                print '(A)', 'Each process reads only part of the MHD simulation data'
            endif
            if (track_particle_flag) then
                print '(A,I0)', 'Number of particles to track: ', nptl_selected
                print '(A,I0)', 'Steps interval to track particles: ', nsteps_interval
            endif
            print '(A,A)', 'Diagnostic file directory is: ', diagnostics_directory
            print '(A,I0)', 'Number of fine steps: ', num_fine_steps
            if (inject_at_shock) then
                print '(A)', 'Inject particles at shock location'
            endif
            print '(A,A)', 'Configuration file name: ', trim(conf_file)
            if (local_dist) then
                print '(A)', 'Simulation will diagnose local particle distribution'
            endif
        endif
    end subroutine get_cmd_args
end program stochastic

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
        distributions_diagnostics, quick_check
    use random_number_generator, only: init_prng, delete_prng
    use mhd_data_parallel, only: init_field_data, free_field_data, &
        read_field_data_parallel, init_fields_gradients, free_fields_gradients, &
        calc_fields_gradients, copy_fields
    use simulation_setup_module, only: read_simuation_mpi_topology, &
        set_field_configuration, fconfig, read_particle_boundary_conditions, &
        set_neighbors
    implicit none
    character(len=256) :: dir_mhd_data
    character(len=256) :: fname1, fname2
    integer :: nptl_max, nptl
    real(dp) :: start, finish, step1, step2, dt
    integer :: t_start, t_end, tf, dist_flag, split_flag
    integer :: interp_flag, nx, ny, nz

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
    call read_simuation_mpi_topology
    call read_particle_boundary_conditions
    call set_field_configuration(whole_data_flag=1, ndim=2)
    call set_neighbors(whole_data_flag=1)

    !< Initialization
    interp_flag = 1 ! Two time are needed for interpolation
    nx = mhd_config%nx
    ny = mhd_config%ny
    nz = mhd_config%nz
    call init_field_data(interp_flag, nx, ny, nz, ndim=2)
    call init_fields_gradients(interp_flag, nx, ny, nz, ndim=2)
    call read_particle_params

    call init_particles(nptl_max)
    call inject_particles_spatial_uniform(nptl, dt, dist_flag)
    call init_particle_distributions

    write(fname1, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', t_start
    write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', t_start + 1
    call read_field_data_parallel(fname1, var_flag=0)
    call read_field_data_parallel(fname2, var_flag=1)

    call calc_fields_gradients(var_flag=0)
    call calc_fields_gradients(var_flag=1)

    call cpu_time(step1)

    call distributions_diagnostics(t_start)

    !< Time loop
    do tf = t_start + 1, t_end
        call particle_mover
        call remove_particles
        if (split_flag == 1) then
            call split_particle
        endif
        if (tf == t_start + 1) then
            call quick_check(tf, .true.)
        else
            call quick_check(tf, .false.)
        endif
        call distributions_diagnostics(tf)
        call copy_fields
        if (tf < t_end) then
            write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', tf + 1
            call read_field_data_parallel(fname2, var_flag=1)
            call calc_fields_gradients(var_flag=1)
        endif
        call cpu_time(step2)
        if (mpi_rank == master) then
            print '("Step ", I0, " takes ", f9.4, " seconds.")', tf, step2 - step1
        endif
        step1 = step2
    enddo

    call free_particle_distributions
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
        use flap                                !< FLAP package
        use penf
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
                        '-np nptl -dt dt -ts t_start -te t_end -df dist_flag'])
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
            print '(A,I0)', 'Starting time frame', t_start
            print '(A,I0)', 'The last time frame', t_end
            if (split_flag == 1) then
                print '(A)', 'Particles will be splitted when reaching certain energies'
            endif
        endif
    end subroutine get_cmd_args
end program stochastic

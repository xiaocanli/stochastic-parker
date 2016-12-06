!*******************************************************************************
!< Program for stochastic integration
!*******************************************************************************

program stochastic
    use constants, only: fp, dp
    use mpi_module
    use mhd_data_sli, only: read_mhd_config, read_mhd_config_from_outfile, &
        broadcast_mhd_config, init_mhd_data, free_mhd_data, read_mhd_data, &
        init_fields_gradients, free_fields_gradients, calc_fields_gradients
    use particle_module, only: init_particles, free_particles, &
        inject_particles_spatial_uniform, read_particle_params, &
        particle_mover, remove_particles
    use random_number_generator, only: init_prng, delete_prng
    implicit none
    character(len=256) :: dir_mhd_data
    character(len=256) :: fname, fname1
    integer :: nptl_max
    real(dp) :: start, finish

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, ierr)

    call cpu_time(start)

    call get_cmd_args
    call init_prng

    fname = trim(dir_mhd_data)//'bin_out0000'
    fname1 = trim(dir_mhd_data)//'bin_out0001'
    if (mpi_rank == master) then
        ! call read_mhd_config
        call read_mhd_config_from_outfile(fname, fname1)
    endif
    call broadcast_mhd_config
    call init_mhd_data
    call init_fields_gradients
    call read_particle_params

    call init_particles(nptl_max)
    call inject_particles_spatial_uniform(100, 0.001_dp)

    call read_mhd_data(fname, var_flag=0)
    call read_mhd_data(fname1, var_flag=1)
    call calc_fields_gradients(var_flag=0)
    call calc_fields_gradients(var_flag=1)

    call particle_mover
    call remove_particles

    call free_particles
    call free_fields_gradients
    call free_mhd_data

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
            examples = ['stochastic-2dmhd.exec -dm dir_mhd_data -nm nptl_max'])
        call cli%add(switch='--dir_mhd_data', switch_ab='-dm', &
            help='MHD simulation data file directory', required=.true., &
            act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--nptl_max', switch_ab='-nm', &
            help='Maximum number of particles', required=.false., &
            act='store', def='1E7', error=error)
        if (error/=0) stop
        call cli%get(switch='-dm', val=dir_mhd_data, error=error)
        if (error/=0) stop
        call cli%get(switch='-nm', val=nptl_max, error=error)
        if (error/=0) stop

        if (mpi_rank == master) then
            print '(A,A)', 'Direcotry of MHD data files: ', trim(dir_mhd_data)
            print '(A,I0)', 'Maximum number of particles: ', nptl_max
        endif
    end subroutine get_cmd_args
end program stochastic

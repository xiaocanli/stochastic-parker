!*******************************************************************************
!< Program for organizing MHD simulation data
!*******************************************************************************
program organize_mhd_data
    use constants, only: fp, dp
    use mhd_data_sli, only: read_mhd_config_from_outfile, init_mhd_data, &
        free_mhd_data, read_mhd_data, save_organized_mhd_data, &
        adjust_mhd_data_boundary
    use mhd_config_module, only: save_mhd_config, echo_mhd_config
    implicit none
    character(len=256) :: dir_mhd_data
    character(len=256) :: fname1, fname2
    real(dp) :: start, finish, step1, step2
    integer :: t_start, t_end, tf

    call cpu_time(start)

    call get_cmd_args

    !< Read and save MHD configuration
    write(fname1, "(A,I4.4)") trim(dir_mhd_data)//'bin_out', t_start
    write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'bin_out', t_start + 1
    call read_mhd_config_from_outfile(fname1, fname2)
    call echo_mhd_config

    write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'mhd_config.dat'
    call save_mhd_config(fname2)

    call init_mhd_data

    call cpu_time(step1)

    do tf = t_start, t_end
        write(fname1, "(A,I4.4)") trim(dir_mhd_data)//'bin_out', tf
        call read_mhd_data(fname1)
        call adjust_mhd_data_boundary
        write(fname2, "(A,I4.4)") trim(dir_mhd_data)//'mhd_data_', tf
        call save_organized_mhd_data(fname2)
        call cpu_time(step2)
        print '("Step ", I0, " takes ", f9.4, " seconds.")', tf, step2 - step1
        step1 = step2
    enddo

    call free_mhd_data

    call cpu_time(finish)
    print '("Time = ",f9.4," seconds.")',finish - start

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
            progname = 'organize_mhd_data', &
            authors     = 'Xiaocan Li', &
            help        = 'Usage: ', &
            description = "Organize the MHD simulation data ", &
            examples = ['organize_mhd_data.exec -dm dir_mhd_data '//&
                        '-ts t_start -te t_end'])
        call cli%add(switch='--dir_mhd_data', switch_ab='-dm', &
            help='MHD simulation data file directory', required=.true., &
            act='store', error=error)
        if (error/=0) stop
        call cli%add(switch='--tstart', switch_ab='-ts', &
            help='Starting time frame', required=.false., &
            act='store', def='0', error=error)
        if (error/=0) stop
        call cli%add(switch='--tend', switch_ab='-te', &
            help='The last time frame', required=.false., &
            act='store', def='200', error=error)
        if (error/=0) stop
        call cli%get(switch='-dm', val=dir_mhd_data, error=error)
        if (error/=0) stop
        call cli%get(switch='-ts', val=t_start, error=error)
        if (error/=0) stop
        call cli%get(switch='-te', val=t_end, error=error)
        if (error/=0) stop

        print '(A,A)', 'Direcotry of MHD data files: ', trim(dir_mhd_data)
        print '(A,I0)', 'Starting time frame ', t_start
        print '(A,I0)', 'The last time frame ', t_end
    end subroutine get_cmd_args
end program organize_mhd_data

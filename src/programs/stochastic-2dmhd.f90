!*******************************************************************************
!< Program for stochastic integration
!*******************************************************************************

program stochastic
    use constants, only: fp, dp
    use mpi_module
    use mhd_data, only: read_mhd_config

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, numprocs, ierr)


    call read_mhd_config

    ! CALL init
    ! CALL read_mhd_data(sframe,3,1)
    ! lframe = 1000-dnMHD*2
    ! DO i = sframe, nframes-dnMHD-2, dnMHD
    !     iframe = i
    !     vxMHD = vxMHD1; vyMHD = vyMHD1
    !     bxMHD = bxMHD1; byMHD = byMHD1
    !     DivVxMHD = DivVxMHD1
    !     DivVyMHD = DivVyMHD1
    !     DivBxxMHD = DivBxxMHD1
    !     DivBxyMHD = DivBxyMHD1
    !     DivByxMHD = DivByxMHD1
    !     DivByyMHD = DivByyMHD1
    !     iMinDivVx = iMinDivVx1
    !     CALL read_mhd_data(iframe+dnMHD,3,1)
    !     CALL inject
    !     CALL pmov(3,1)
    !     CALL reduce_part
    !     CALL split
    !     ! IF(MOD(i,1000) .EQ. lframe) THEN
    !         CALL output
    !     ! ENDIF
    !     IF (mype .EQ. 0) THEN
    !         total = 0
    !         DO j = 1,npart
    !             total = total + weight(i)
    !         ENDDO
    !         OPEN (17,FILE='quick.dat',STATUS='unknown')
    !         WRITE (17,*) iframe,npart,ndaughter,total,leak
    !         CLOSE(17)
    !     ENDIF
    ! ENDDO
    call MPI_FINALIZE(ierr)

    contains

    !---------------------------------------------------------------------------
    !< Get commandline arguments
    !---------------------------------------------------------------------------
    subroutine get_cmd_args
        ! use flap                                !< FLAP package
        ! use penf
        implicit none
        ! type(command_line_interface) :: cli     !< Command Line Interface (CLI).
        ! integer(I4P)                 :: error   !< Error trapping flag.
        ! integer :: tstart
        ! call cli%init(&
        !     progname = 'stochastic-2dmhd', &
        !     authors     = 'Xiaocan Li', &
        !     help        = 'Usage: ', &
        !     description = "Solving Parker's transport equation "// &
        !                   "using stochastic differential equation", &
        !     examples = ['haha'])
        ! call cli%add(switch='--tstart', switch_ab='-ts', &
        !     help='Starting time frame', required=.false., act='store', &
        !     def='0', error=error)
        ! if (error/=0) stop
        ! call cli%get(switch='-ts', val=tstart, error=error)
        ! if (error/=0) stop

        ! if (myid == 0) then
        !     print '(A,L1)', 'Whether using translated fields file: ', is_translated_file
        !     print '(A,L1)', 'Whether reading one-fluid velocity: ', &
        !         if_read_one_fluid_velocity
        !     print '(A,A)', 'Tracer directory: ', dir_tracer_hdf5
        !     print '(A,I0,A,I0,A,I0)', 'Min, max and interval: ', &
        !         tstart, ' ', tend, ' ', tinterval
        !     if (species == 'e') then
        !         print '(A,A)', 'Particle: electron'
        !     else if (species == 'h' .or. species == 'i') then
        !         print '(A,A)', 'Particle: ion'
        !     endif
        !     print '(A,A)', 'Tracer filename: ', trim(fname_tracer)
        !     print '(A,A)', 'Metadata filename: ', trim(fname_metadata)
        !     print '(A,A)', 'Root directory: ', trim(rpath)
        !     print '(A,A)', 'EMF data directory: ', trim(dir_emf)
        !     print '(A,A)', 'Hydro data directory: ', trim(dir_hydro)
        ! endif
    end subroutine get_cmd_args
end program stochastic

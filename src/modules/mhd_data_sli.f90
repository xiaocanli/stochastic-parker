!*******************************************************************************
!< Module of MHD data for Shengtai Li's MHD code
!*******************************************************************************
module mhd_data_sli
    use constants, only: fp, dp
    implicit none
    private
    public mhd_config
    public read_mhd_config, read_mhd_config_from_outfile, init_mhd_data, &
           free_mhd_data, read_mhd_data, broadcast_mhd_config

    type mhd_configuration
        real(dp) :: dx, dy, xmin, xmax, ymin, ymax, lx, ly  ! Grid sizes
        real(dp) :: dt_out      ! Time interval for MHD data output
        integer :: nx, ny       ! Grid dimensions
        integer :: nxs, nys     ! Grid dimensions for a single MPI process
        integer :: nvar         ! Number of output variables
        integer :: topox, topoy ! MHD simulation topology
    end type mhd_configuration

    type file_header
        integer, dimension(4) :: nx4
        real(fp) :: time
        real(fp), dimension(4) :: bbox
    end type file_header

    type(mhd_configuration) :: mhd_config
    type(file_header) :: fheader

    real(fp), allocatable, dimension(:, :) :: rho, pres, vx, vy, vz, bx, by, bz
    real(fp), allocatable, dimension(:, :, :) :: mhd_data_single_core

    contains

    !---------------------------------------------------------------------------
    !< Read MHD simulation configuration from the configuration file
    !---------------------------------------------------------------------------
    subroutine read_mhd_config
        use read_config, only: get_variable
        implicit none
        real(fp) :: temp
        integer :: fh
        fh = 10
        open(unit=fh, file='config/mhd_config.dat', status='old')
        mhd_config%xmin = get_variable(fh, 'xmin', '=')
        mhd_config%xmax = get_variable(fh, 'xmax', '=')
        mhd_config%ymin = get_variable(fh, 'ymin', '=')
        mhd_config%ymax = get_variable(fh, 'ymax', '=')
        mhd_config%lx = mhd_config%xmax - mhd_config%xmin
        mhd_config%ly = mhd_config%ymax - mhd_config%ymin
        temp = get_variable(fh, 'nx', '=')
        mhd_config%nx = int(temp)
        temp = get_variable(fh, 'ny', '=')
        mhd_config%ny = int(temp)
        mhd_config%dx = mhd_config%lx / mhd_config%nx
        mhd_config%dy = mhd_config%ly / mhd_config%ny
        temp = get_variable(fh, 'nvar', '=')
        mhd_config%nvar = int(temp)
        temp = get_variable(fh, 'topology_x', '=')
        mhd_config%topox = int(temp)
        temp = get_variable(fh, 'topology_y', '=')
        mhd_config%topoy = int(temp)
        mhd_config%dt_out = get_variable(fh, 'dt_out', '=')

        mhd_config%nxs = mhd_config%nx / mhd_config%topox
        mhd_config%nys = mhd_config%ny / mhd_config%topoy
        close(fh)

        call echo_mhd_config_info
    end subroutine read_mhd_config

    !---------------------------------------------------------------------------
    !< Echo MHD configuration information
    !---------------------------------------------------------------------------
    subroutine echo_mhd_config_info
        implicit none
        print *, "---------------------------------------------------"
        write(*, "(A)") " MHD simulation information."
        write(*, "(A,F7.2,A,F7.2)") " lx, ly = ", &
            mhd_config%lx, ',', mhd_config%ly
        write(*, "(A,I0,A,I0)") " nx, ny = ", &
            mhd_config%nx, ',', mhd_config%ny
        write(*, "(A,F9.6,A,F9.6)") " dx, dy = ", &
            mhd_config%dx, ',', mhd_config%dy
        write(*, "(A,I0)") " Number of output variables: ", mhd_config%nvar
        write(*, "(A,I0,A,I0)") " MHD topology: ", mhd_config%topox, &
            " * ", mhd_config%topoy
        write(*, "(A,I0,A,I0)") " Grid dimensions at each MPI rank: ", &
            mhd_config%nxs, ",", mhd_config%nys
        write(*, "(A,F7.3)") " Time interval for MHD data output: ", &
            mhd_config%dt_out
        print *, "---------------------------------------------------"
    end subroutine echo_mhd_config_info

    !---------------------------------------------------------------------------
    !< Read MHD simulation configuration from the simulation output file
    !< Args:
    !<  filename: the name of one output file
    !<  filename_next_step: the name of the output file of the next time step
    !---------------------------------------------------------------------------
    subroutine read_mhd_config_from_outfile(filename, filename_next_step)
        implicit none
        character(*), intent(in) :: filename, filename_next_step
        integer, dimension(4) :: n4
        real(dp) :: t1
        integer :: fh
        fh = 20
        open(unit=fh, file=filename, access='stream', status='unknown', &
             form='unformatted', action='read')
        read(fh) fheader
        read(fh) n4
        mhd_config%nvar = fheader%nx4(1)
        mhd_config%nx = fheader%nx4(2)
        mhd_config%ny = fheader%nx4(3)
        mhd_config%nxs = n4(3)
        mhd_config%nys = n4(4)
        mhd_config%topox = mhd_config%nx / mhd_config%nxs
        mhd_config%topoy = mhd_config%ny / mhd_config%nys
        mhd_config%xmin = fheader%bbox(1)
        mhd_config%xmax = fheader%bbox(2)
        mhd_config%ymin = fheader%bbox(3)
        mhd_config%ymax = fheader%bbox(4)
        mhd_config%lx = mhd_config%xmax - mhd_config%xmin
        mhd_config%ly = mhd_config%ymax - mhd_config%ymin
        mhd_config%dx = mhd_config%lx / mhd_config%nx
        mhd_config%dy = mhd_config%ly / mhd_config%ny
        close(fh)

        t1 = fheader%time
        open(unit=fh, file=filename_next_step, access='stream', &
             status='unknown', form='unformatted', action='read')
        read(fh) fheader
        close(fh)
        mhd_config%dt_out = fheader%time - t1

        call echo_mhd_config_info
    end subroutine read_mhd_config_from_outfile

    !---------------------------------------------------------------------------
    !< Broadcast MHD configuration
    !---------------------------------------------------------------------------
    subroutine broadcast_mhd_config
        use mpi_module
        implicit none
        integer :: mhd_config_type, oldtypes(0:1), blockcounts(0:1)
        integer :: offsets(0:1), extent
        ! Setup description of the 8 MPI_DOUBLE fields.
        offsets(0) = 0
        oldtypes(0) = MPI_DOUBLE_PRECISION
        blockcounts(0) = 9
        ! Setup description of the 7 MPI_INTEGER fields.
        call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent, ierr)
        offsets(1) = blockcounts(0) * extent
        oldtypes(1) = MPI_INTEGER
        blockcounts(1) = 7
        ! Define structured type and commit it. 
        call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, &
            mhd_config_type, ierr)
        call MPI_TYPE_COMMIT(mhd_config_type, ierr)
        call MPI_BCAST(mhd_config, 1, mhd_config_type, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_TYPE_FREE(mhd_config_type, ierr)
    end subroutine broadcast_mhd_config

    !---------------------------------------------------------------------------
    !< Initialize MHD data arrays. Note that we do not read Az.
    !---------------------------------------------------------------------------
    subroutine init_mhd_data
        implicit none
        integer :: nx, ny, nvar, nxs, nys
        nx = mhd_config%nx
        ny = mhd_config%ny
        nvar = mhd_config%nvar
        nxs = mhd_config%nxs
        nys = mhd_config%nys
        allocate(rho(nx, ny))
        allocate(pres(nx, ny))
        allocate(vx(nx, ny))
        allocate(vy(nx, ny))
        allocate(bx(nx, ny))
        allocate(by(nx, ny))
        rho = 0.0
        pres = 0.0
        vx = 0.0
        vy = 0.0
        bx = 0.0
        by = 0.0
        if (mhd_config%nvar .gt. 7) then !< rho, p, vx, vy, vz, bx, by, bz, Az
            allocate(vz(nx, ny))
            allocate(bz(nx, ny))
            vz = 0.0
            bz = 0.0
        endif
        allocate(mhd_data_single_core(nxs, nys, nvar))
        mhd_data_single_core = 0.0
    end subroutine init_mhd_data

    !---------------------------------------------------------------------------
    !< Free MHD data arrays
    !---------------------------------------------------------------------------
    subroutine free_mhd_data
        implicit none
        if (mhd_config%nvar .eq. 7) then
            deallocate(rho, pres, vx, vy, bx, by)
        else
            deallocate(rho, pres, vx, vy, vz, bx, by, bz)
        endif
        deallocate(mhd_data_single_core)
    end subroutine free_mhd_data

    !---------------------------------------------------------------------------
    !< Read MHD data arrays from a file
    !---------------------------------------------------------------------------
    subroutine read_mhd_data(filename)
        implicit none
        character(*), intent(in) :: filename
        integer, dimension(4) :: n4
        integer :: fh, rank, ix, iy, nx1, ny1, ncpus
        fh = 30
        open(unit=fh, file=filename, access='stream', status='unknown', &
             form='unformatted', action='read')
        read(fh) fheader
        ncpus = mhd_config%topox * mhd_config%topoy
        do rank = 1, ncpus
            read(fh) n4
            ix = n4(1)
            iy = n4(2)
            nx1 = n4(3)
            ny1 = n4(4)
            read(fh) mhd_data_single_core
            rho(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 1)
            pres(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 2)
            if (mhd_config%nvar .eq. 7) then
                vx(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 3)
                vy(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 4)
                bx(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 5)
                by(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 6)
            else
                vx(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 3)
                vy(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 4)
                vz(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 5)
                bx(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 6)
                by(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 7)
                bz(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 8)
            endif
        enddo
        close(fh)
    end subroutine read_mhd_data
end module mhd_data_sli

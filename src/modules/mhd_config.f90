!*******************************************************************************
!< Module of MHD configuration
!*******************************************************************************
module mhd_config_module
    use constants, only: dp
    implicit none
    private
    save
    public mhd_config
    public broadcast_mhd_config, save_mhd_config, load_mhd_config, &
           set_mhd_config, echo_mhd_config, set_mhd_grid_type
    public uniform_grid_flag, spherical_coord_flag
    public tstamps_mhd, init_tstamps_mhd, free_tstamps_mhd, load_tstamps_mhd, &
        calc_tstamps_mhd, dt_mhd_max, tstart_mhd, tend_mhd

    type mhd_configuration
        real(dp) :: dx = 1.0_dp, dy = 1.0_dp, dz = 1.0_dp
        real(dp) :: xmin = 0.0_dp, ymin = 0.0_dp, zmin = 0.0_dp
        real(dp) :: xmax = 1.0_dp, ymax = 1.0_dp, zmax = 1.0_dp
        real(dp) :: lx = 1.0_dp, ly = 1.0_dp, lz = 1.0_dp
        real(dp) :: dt_out = 0.25_dp               ! Time interval for MHD data output
        integer :: nx = 1, ny = 1, nz = 1          ! Grid dimensions
        integer :: nxs = 1, nys = 1, nzs = 1       ! Grid dimensions for a single MPI process
        integer :: topox = 1, topoy = 1, topoz = 1 ! MHD simulation topology
        integer :: nvar = 9                        ! Number of output variables
        integer :: bcx = 0, bcy = 0, bcz = 0       ! 0 for periodic, 1 for others
    end type mhd_configuration

    type(mhd_configuration) :: mhd_config
    logical :: uniform_grid_flag     ! whether the grid is uniform
    logical :: spherical_coord_flag  ! whether the grid is in spherical coordinates

    real(dp) :: dt_mhd_max  !< Maximum MHD fields interval
    real(dp), allocatable, dimension(:) :: tstamps_mhd ! Time stamps for MHD fields

    integer :: tstart_mhd, tend_mhd

    contains

    !---------------------------------------------------------------------------
    !< Set MHD configuration information
    !---------------------------------------------------------------------------
    subroutine set_mhd_config(dx, dy, dz, xmin, ymin, zmin, xmax, ymax, zmax, &
            lx, ly, lz, dt_out, nx, ny, nz, nxs, nys, nzs, topox, topoy, topoz, &
            nvar, bcx, bcy, bcz)
        implicit none
        real(dp), intent(in) :: dx, dy, dz, xmin, ymin, zmin, xmax, ymax, zmax
        real(dp), intent(in) :: lx, ly, lz, dt_out
        integer, intent(in) :: nx, ny, nz, nxs, nys, nzs
        integer, intent(in) :: topox, topoy, topoz, nvar
        integer, intent(in) :: bcx, bcy, bcz
        mhd_config%dx = dx
        mhd_config%dy = dy
        mhd_config%dz = dz
        mhd_config%xmin = xmin
        mhd_config%ymin = ymin
        mhd_config%zmin = zmin
        mhd_config%xmax = xmax
        mhd_config%ymax = ymax
        mhd_config%zmax = zmax
        mhd_config%lx = lx
        mhd_config%ly = ly
        mhd_config%lz = lz
        mhd_config%dt_out = dt_out
        mhd_config%nx = nx
        mhd_config%ny = ny
        mhd_config%nz = nz
        mhd_config%nxs = nxs
        mhd_config%nys = nys
        mhd_config%nzs = nzs
        mhd_config%topox = topox
        mhd_config%topoy = topoy
        mhd_config%topoz = topoz
        mhd_config%nvar = nvar
        mhd_config%bcx = bcx
        mhd_config%bcy = bcy
        mhd_config%bcz = bcz
    end subroutine set_mhd_config

    !---------------------------------------------------------------------------
    !< Echo MHD configuration information
    !---------------------------------------------------------------------------
    subroutine echo_mhd_config
        implicit none
        print *, "---------------------------------------------------"
        write(*, "(A)") " MHD simulation information."
        write(*, "(A,F7.2,A,F7.2,A,F7.2)") " lx, ly, lz = ", &
            mhd_config%lx, ',', mhd_config%ly, ',', mhd_config%lz
        write(*, "(A,F9.6,A,F9.6,A,F9.6)") " dx, dy, dz = ", &
            mhd_config%dx, ',', mhd_config%dy, ',', mhd_config%dz
        write(*, "(A,F12.6,A,F12.6,A,F12.6)") " xmin, ymin, zmin = ", &
            mhd_config%xmin, ',', mhd_config%ymin, ',', mhd_config%zmin
        write(*, "(A,F12.6,A,F12.6,A,F12.6)") " xmax, ymax, zmax = ", &
            mhd_config%xmax, ',', mhd_config%ymax, ',', mhd_config%zmax
        write(*, "(A,I0,A,I0,A,I0)") " nx, ny, nz = ", &
            mhd_config%nx, ',', mhd_config%ny, ',', mhd_config%nz
        write(*, "(A,I0)") " Number of output variables: ", mhd_config%nvar
        write(*, "(A,I0,A,I0,A,I0)") " MHD topology: ", mhd_config%topox, &
            " * ", mhd_config%topoy, " * ", mhd_config%topoz
        write(*, "(A,I0,A,I0,A,I0)") " Grid dimensions at each MPI rank: ", &
            mhd_config%nxs, ",", mhd_config%nys, ",", mhd_config%nzs
        write(*, "(A,E13.6)") " Time interval for MHD data output: ", &
            mhd_config%dt_out
        if (mhd_config%bcx == 0) then
            write(*, "(A)") " MHD simulation is periodic along x"
        else
            write(*, "(A)") " MHD simulation is not periodic along x"
        endif
        if (mhd_config%bcy == 0) then
            write(*, "(A)") " MHD simulation is periodic along y"
        else
            write(*, "(A)") " MHD simulation is not periodic along y"
        endif
        if (mhd_config%bcz == 0) then
            write(*, "(A)") " MHD simulation is periodic along z"
        else
            write(*, "(A)") " MHD simulation is not periodic along z"
        endif
        print *, "---------------------------------------------------"
    end subroutine echo_mhd_config

    !---------------------------------------------------------------------------
    !< Save MHD configuration to file
    !---------------------------------------------------------------------------
    subroutine save_mhd_config(filename)
        implicit none
        character(*), intent(in) :: filename
        integer :: fh
        fh = 25
        open(unit=fh, file=filename, access='stream', status='unknown', &
             form='unformatted', action='write')
        write(fh, pos=1) mhd_config
        close(fh)
    end subroutine save_mhd_config

    !---------------------------------------------------------------------------
    !< Load MHD configuration
    !---------------------------------------------------------------------------
    subroutine load_mhd_config(filename)
        implicit none
        character(*), intent(in) :: filename
        integer :: fh
        fh = 25
        open(unit=fh, file=filename, access='stream', status='unknown', &
             form='unformatted', action='read')
        read(fh, pos=1) mhd_config
        close(fh)
    end subroutine load_mhd_config

    !---------------------------------------------------------------------------
    !< Broadcast MHD configuration
    !---------------------------------------------------------------------------
    subroutine broadcast_mhd_config
        use mpi_module
        implicit none
        integer :: mhd_config_type, oldtypes(0:1), blockcounts(0:1)
        integer :: offsets(0:1), extent
        ! Setup description of the MPI_DOUBLE fields.
        offsets(0) = 0
        oldtypes(0) = MPI_DOUBLE_PRECISION
        blockcounts(0) = 13
        ! Setup description of the MPI_INTEGER fields.
        call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, extent, ierr)
        offsets(1) = blockcounts(0) * extent
        oldtypes(1) = MPI_INTEGER
        blockcounts(1) = 13
        ! Define structured type and commit it.
        call MPI_TYPE_STRUCT(2, blockcounts, offsets, oldtypes, &
            mhd_config_type, ierr)
        call MPI_TYPE_COMMIT(mhd_config_type, ierr)
        call MPI_BCAST(mhd_config, 1, mhd_config_type, master, &
            MPI_COMM_WORLD, ierr)
        call MPI_TYPE_FREE(mhd_config_type, ierr)
    end subroutine broadcast_mhd_config

    !---------------------------------------------------------------------------
    !< Set MHD grid type
    !---------------------------------------------------------------------------
    subroutine set_mhd_grid_type(uniform_grid, spherical_coord)
        implicit none
        integer, intent(in) :: uniform_grid, spherical_coord
        uniform_grid_flag = .false.
        spherical_coord_flag = .false.
        if (uniform_grid == 1) uniform_grid_flag = .true.
        if (spherical_coord == 1) spherical_coord_flag = .true.
    end subroutine set_mhd_grid_type

    !---------------------------------------------------------------------------
    !< Initialize the time stamps of MHD fields
    !< Args:
    !<  t_start: starting time frame
    !<  t_end: ending time frame
    !---------------------------------------------------------------------------
    subroutine init_tstamps_mhd(t_start, t_end)
        implicit none
        integer, intent(in) :: t_start, t_end
        tstart_mhd = t_start
        tend_mhd = t_end
        allocate(tstamps_mhd(t_end-t_start+1))
    end subroutine init_tstamps_mhd

    !---------------------------------------------------------------------------
    !< Free the time stamps of MHD fields
    !---------------------------------------------------------------------------
    subroutine free_tstamps_mhd
        implicit none
        deallocate(tstamps_mhd)
    end subroutine free_tstamps_mhd

    !---------------------------------------------------------------------------
    !< Load the time stamps of MHD fields
    !< The default filename for the time stamps is "time_stamps.dat"
    !< The first two data points are the starting and ending frame of the MHD
    !< simulation output.
    !< Args:
    !<  t_start: starting time frame for the current simulation
    !<  dir_mhd_data: directory for the MHD data
    !---------------------------------------------------------------------------
    subroutine load_tstamps_mhd(t_start, dir_mhd_data)
        use mpi_module
        implicit none
        integer, intent(in) :: t_start
        character(*), intent(in) :: dir_mhd_data
        character(256) :: filename
        integer :: fh, pos, ts_mhd, i
        real(dp) :: dt_mhd
        if (mpi_rank == master) then
            filename = trim(dir_mhd_data)//"time_stamps.dat"
            fh = 25
            open(unit=fh, file=filename, access='stream', status='unknown', &
                 form='unformatted', action='read')
            pos = 1
            read(fh, pos=pos) ts_mhd
            pos = pos + (t_start - ts_mhd + 2) * 8
            read(fh, pos=pos) tstamps_mhd
            close(fh)
        endif
        call MPI_BCAST(tstamps_mhd, size(tstamps_mhd), &
            MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        dt_mhd_max = 0.0_dp
        do i = 2, size(tstamps_mhd)
            dt_mhd = tstamps_mhd(i) - tstamps_mhd(i-1)
            if (dt_mhd > dt_mhd_max) then
                dt_mhd_max = dt_mhd
            endif
        enddo
    end subroutine load_tstamps_mhd

    !---------------------------------------------------------------------------
    !< Calculate the time stamps of MHD fields when the time interval is uniform
    !< Args:
    !<  t_start: starting time frame for the current simulation
    !<  t_end: ending time frame for the current simulation
    !---------------------------------------------------------------------------
    subroutine calc_tstamps_mhd(t_start, t_end)
        implicit none
        integer, intent(in) :: t_start, t_end
        integer :: i
        do i = t_start, t_end
            tstamps_mhd(i-t_start+1) = i * mhd_config%dt_out
        enddo
        dt_mhd_max = mhd_config%dt_out
    end subroutine calc_tstamps_mhd
end module mhd_config_module

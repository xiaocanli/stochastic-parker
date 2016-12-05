!*******************************************************************************
!< Module of MHD data for Shengtai Li's MHD code
!*******************************************************************************
module mhd_data_sli
    use constants, only: fp, dp
    implicit none
    private
    public mhd_config
    public read_mhd_config, read_mhd_config_from_outfile, init_mhd_data, &
           free_mhd_data, read_mhd_data, broadcast_mhd_config, &
           init_fields_gradients, free_fields_gradients, calc_fields_gradients

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

    !< Two sets of data for two different time frames
    real(fp), allocatable, dimension(:, :) :: rho, pres, rho1, pres1
    real(fp), allocatable, dimension(:, :) :: vx, vy, vz, vx1, vy1, vz1
    real(fp), allocatable, dimension(:, :) :: bx, by, bz, btot
    real(fp), allocatable, dimension(:, :) :: bx1, by1, bz1, btot1
    real(fp), allocatable, dimension(:, :, :) :: mhd_data_single_core

    !< Gradients of the fields
    real(fp), allocatable, dimension(:, :) :: dvx_dx, dvy_dy, dvx1_dx, dvy1_dy
    real(fp), allocatable, dimension(:, :) :: dbx_dx, dbx_dy, dbx1_dx, dbx1_dy
    real(fp), allocatable, dimension(:, :) :: dby_dx, dby_dy, dby1_dx, dby1_dy
    real(fp), allocatable, dimension(:, :) :: dbtot_dx, dbtot_dy
    real(fp), allocatable, dimension(:, :) :: dbtot1_dx, dbtot1_dy

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
        allocate(btot(nx, ny))
        allocate(rho1(nx, ny))
        allocate(pres1(nx, ny))
        allocate(vx1(nx, ny))
        allocate(vy1(nx, ny))
        allocate(bx1(nx, ny))
        allocate(by1(nx, ny))
        allocate(btot1(nx, ny))

        rho = 0.0
        pres = 0.0
        vx = 0.0
        vy = 0.0
        bx = 0.0
        by = 0.0
        btot = 0.0
        rho1 = 0.0
        pres1 = 0.0
        vx1 = 0.0
        vy1 = 0.0
        bx1 = 0.0
        by1 = 0.0
        btot1 = 0.0

        if (mhd_config%nvar .gt. 7) then !< rho, p, vx, vy, vz, bx, by, bz, Az
            allocate(vz(nx, ny))
            allocate(bz(nx, ny))
            allocate(vz1(nx, ny))
            allocate(bz1(nx, ny))
            vz = 0.0
            bz = 0.0
            vz1 = 0.0
            bz1 = 0.0
        endif
        allocate(mhd_data_single_core(nxs, nys, nvar))
        mhd_data_single_core = 0.0
    end subroutine init_mhd_data

    !---------------------------------------------------------------------------
    !< Initialize the gradients of the MHD data arrays.
    !---------------------------------------------------------------------------
    subroutine init_fields_gradients
        implicit none
        integer :: nx, ny, nvar, nxs, nys
        nx = mhd_config%nx
        ny = mhd_config%ny
        allocate(dvx_dx(nx, ny))
        allocate(dvy_dy(nx, ny))
        allocate(dvx1_dx(nx, ny))
        allocate(dvy1_dy(nx, ny))
        allocate(dbx_dx(nx, ny))
        allocate(dbx_dy(nx, ny))
        allocate(dby_dx(nx, ny))
        allocate(dby_dy(nx, ny))
        allocate(dbtot_dx(nx, ny))
        allocate(dbtot_dy(nx, ny))
        allocate(dbx1_dx(nx, ny))
        allocate(dbx1_dy(nx, ny))
        allocate(dby1_dx(nx, ny))
        allocate(dby1_dy(nx, ny))
        allocate(dbtot1_dx(nx, ny))
        allocate(dbtot1_dy(nx, ny))
        dvx_dx = 0.0
        dvy_dy = 0.0
        dvx1_dx = 0.0
        dvy1_dy = 0.0
        dbx_dx = 0.0
        dbx_dy = 0.0
        dby_dx = 0.0
        dby_dy = 0.0
        dbtot_dx = 0.0
        dbtot_dy = 0.0
        dbx1_dx = 0.0
        dbx1_dy = 0.0
        dby1_dx = 0.0
        dby1_dy = 0.0
        dbtot1_dx = 0.0
        dbtot1_dy = 0.0
    end subroutine init_fields_gradients

    !---------------------------------------------------------------------------
    !< Free MHD data arrays
    !---------------------------------------------------------------------------
    subroutine free_mhd_data
        implicit none
        if (mhd_config%nvar .eq. 7) then
            deallocate(rho, pres, vx, vy, bx, by, btot)
            deallocate(rho1, pres1, vx1, vy1, bx1, by1, btot1)
        else
            deallocate(rho, pres, vx, vy, vz, bx, by, bz, btot)
            deallocate(rho1, pres1, vx1, vy1, vz1, bx1, by1, bz1, btot1)
        endif
        deallocate(mhd_data_single_core)
    end subroutine free_mhd_data

    !---------------------------------------------------------------------------
    !< Free the gradients of the MHD data arrays.
    !---------------------------------------------------------------------------
    subroutine free_fields_gradients
        implicit none
        deallocate(dvx_dx, dvy_dy, dvx1_dx, dvy1_dy)
        deallocate(dbx_dx, dbx_dy, dbx1_dx, dbx1_dy)
        deallocate(dby_dx, dby_dy, dby1_dx, dby1_dy)
        deallocate(dbtot_dx, dbtot_dy, dbtot1_dx, dbtot1_dy)
    end subroutine free_fields_gradients

    !---------------------------------------------------------------------------
    !< Read MHD data arrays from a file
    !< Args:
    !<  filename: file name to get the data
    !<  var_flag: indicating which set of variables to save the data. 0 for
    !<            rho, pres etc. and other numbers for rho1, pres1 etc.
    !---------------------------------------------------------------------------
    subroutine read_mhd_data(filename, var_flag)
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: var_flag
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
            if (var_flag == 0) then
                rho(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 1)
                pres(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 2)
            else
                rho1(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 1)
                pres1(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 2)
            endif
            if (mhd_config%nvar .eq. 7) then
                if (var_flag == 0) then
                    vx(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 3)
                    vy(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 4)
                    bx(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 5)
                    by(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 6)
                else
                    vx1(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 3)
                    vy1(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 4)
                    bx1(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 5)
                    by1(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 6)
                endif
            else
                if (var_flag == 0) then
                    vx(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 3)
                    vy(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 4)
                    vz(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 5)
                    bx(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 6)
                    by(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 7)
                    bz(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 8)
                else
                    vx1(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 3)
                    vy1(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 4)
                    vz1(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 5)
                    bx1(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 6)
                    by1(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 7)
                    bz1(ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 8)
                endif
            endif
        enddo

        !< Calculate the magnitude of the magnetic field
        if (mhd_config%nvar .eq. 7) then
            if (var_flag == 0) then
                btot = sqrt(bx**2 + by**2)
            else
                btot1 = sqrt(bx1**2 + by1**2)
            endif
        else
            if (var_flag == 0) then
                btot = sqrt(bx**2 + by**2 + bz**2)
            else
                btot1 = sqrt(bx1**2 + by1**2 + bz1**2)
            endif
        endif
        close(fh)
    end subroutine read_mhd_data

    !---------------------------------------------------------------------------
    !< Calculate the gradients of the MHD data arrays.
    !< Args:
    !<  var_flag: indicating which set of variables.
    !<            0 for dvx_dx etc. and other numbers for dvx1_dx etc.
    !---------------------------------------------------------------------------
    subroutine calc_fields_gradients(var_flag)
        implicit none
        integer, intent(in) :: var_flag
        real(dp) :: idxh, idyh
        integer :: nx, ny, nx1, nx2, ny1, ny2
        idxh = 0.5_dp / mhd_config%dx
        idyh = 0.5_dp / mhd_config%dy
        nx = mhd_config%nx
        ny = mhd_config%ny
        nx1 = nx - 1
        nx2 = nx - 2
        ny1 = ny - 1
        ny2 = ny - 2
        if (var_flag == 0) then
            dvx_dx(2:nx1, :) = (vx(3:nx, :) - vx(1:nx2, :)) * idxh
            dvx_dx(1, :) = (-3.0*vx(1, :) + 4.0*vx(2, :) - vx(3, :)) * idxh
            dvx_dx(nx, :) = (3.0*vx(nx, :) - 4.0*vx(nx1, :) + vx(nx2, :)) * idxh
            dvy_dy(2:ny1, :) = (vy(:, 3:ny)-vy(:, 1:ny2)) * idyh
            dvy_dy(1, :) = (-3.0*vy(:, 1) + 4.0*vy(:, 2) - vy(:, 3)) * idyh
            dvy_dy(ny, :) = (3.0*vy(:, ny) - 4.0*vy(:, ny1) + vy(:, ny2)) * idyh

            dbx_dx(2:nx1, :) = (bx(3:nx, :) - bx(1:nx2, :)) * idxh
            dbx_dx(1, :) = (-3.0*bx(1, :) + 4.0*bx(2, :) - bx(3, :)) * idxh
            dbx_dx(nx, :) = (3.0*bx(nx, :) - 4.0*bx(nx1, :) + bx(nx2, :)) * idxh
            dbx_dy(2:ny1, :) = (bx(:, 3:ny)-bx(:, 1:ny2)) * idyh
            dbx_dy(1, :) = (-3.0*bx(:, 1) + 4.0*bx(:, 2) - bx(:, 3)) * idyh
            dbx_dy(ny, :) = (3.0*bx(:, ny) - 4.0*bx(:, ny1) + bx(:, ny2)) * idyh

            dby_dx(2:nx1, :) = (by(3:nx, :) - by(1:nx2, :)) * idxh
            dby_dx(1, :) = (-3.0*by(1, :) + 4.0*by(2, :) - by(3, :)) * idxh
            dby_dx(nx, :) = (3.0*by(nx, :) - 4.0*by(nx1, :) + by(nx2, :)) * idxh
            dby_dy(2:ny1, :) = (by(:, 3:ny)-by(:, 1:ny2)) * idyh
            dby_dy(1, :) = (-3.0*by(:, 1) + 4.0*by(:, 2) - by(:, 3)) * idyh
            dby_dy(ny, :) = (3.0*by(:, ny) - 4.0*by(:, ny1) + by(:, ny2)) * idyh

            dbtot_dx(2:nx1, :) = (btot(3:nx, :) - btot(1:nx2, :)) * idxh
            dbtot_dx(1, :) = (-3.0*btot(1, :) + 4.0*btot(2, :) - btot(3, :)) * idxh
            dbtot_dx(nx, :) = (3.0*btot(nx, :) - 4.0*btot(nx1, :) + btot(nx2, :)) * idxh
            dbtot_dy(2:ny1, :) = (btot(:, 3:ny)-btot(:, 1:ny2)) * idyh
            dbtot_dy(1, :) = (-3.0*btot(:, 1) + 4.0*btot(:, 2) - btot(:, 3)) * idyh
            dbtot_dy(ny, :) = (3.0*btot(:, ny) - 4.0*btot(:, ny1) + btot(:, ny2)) * idyh
        else
            dvx1_dx(2:nx1, :) = (vx1(3:nx, :) - vx1(1:nx2, :)) * idxh
            dvx1_dx(1, :) = (-3.0*vx1(1, :) + 4.0*vx1(2, :) - vx1(3, :)) * idxh
            dvx1_dx(nx, :) = (3.0*vx1(nx, :) - 4.0*vx1(nx1, :) + vx1(nx2, :)) * idxh
            dvy1_dy(2:ny1, :) = (vy1(:, 3:ny)-vy1(:, 1:ny2)) * idyh
            dvy1_dy(1, :) = (-3.0*vy1(:, 1) + 4.0*vy1(:, 2) - vy1(:, 3)) * idyh
            dvy1_dy(ny, :) = (3.0*vy1(:, ny) - 4.0*vy1(:, ny1) + vy1(:, ny2)) * idyh

            dbx1_dx(2:nx1, :) = (bx1(3:nx, :) - bx1(1:nx2, :)) * idxh
            dbx1_dx(1, :) = (-3.0*bx1(1, :) + 4.0*bx1(2, :) - bx1(3, :)) * idxh
            dbx1_dx(nx, :) = (3.0*bx1(nx, :) - 4.0*bx1(nx1, :) + bx1(nx2, :)) * idxh
            dbx1_dy(2:ny1, :) = (bx1(:, 3:ny)-bx1(:, 1:ny2)) * idyh
            dbx1_dy(1, :) = (-3.0*bx1(:, 1) + 4.0*bx1(:, 2) - bx1(:, 3)) * idyh
            dbx1_dy(ny, :) = (3.0*bx1(:, ny) - 4.0*bx1(:, ny1) + bx1(:, ny2)) * idyh

            dby1_dx(2:nx1, :) = (by1(3:nx, :) - by1(1:nx2, :)) * idxh
            dby1_dx(1, :) = (-3.0*by1(1, :) + 4.0*by1(2, :) - by1(3, :)) * idxh
            dby1_dx(nx, :) = (3.0*by1(nx, :) - 4.0*by1(nx1, :) + by1(nx2, :)) * idxh
            dby1_dy(2:ny1, :) = (by1(:, 3:ny)-by1(:, 1:ny2)) * idyh
            dby1_dy(1, :) = (-3.0*by1(:, 1) + 4.0*by1(:, 2) - by1(:, 3)) * idyh
            dby1_dy(ny, :) = (3.0*by1(:, ny) - 4.0*by1(:, ny1) + by1(:, ny2)) * idyh

            dbtot1_dx(2:nx1, :) = (btot1(3:nx, :) - btot1(1:nx2, :)) * idxh
            dbtot1_dx(1, :) = (-3.0*btot1(1, :) + 4.0*btot1(2, :) - btot1(3, :)) * idxh
            dbtot1_dx(nx, :) = (3.0*btot1(nx, :) - 4.0*btot1(nx1, :) + btot1(nx2, :)) * idxh
            dbtot1_dy(2:ny1, :) = (btot1(:, 3:ny)-btot1(:, 1:ny2)) * idyh
            dbtot1_dy(1, :) = (-3.0*btot1(:, 1) + 4.0*btot1(:, 2) - btot1(:, 3)) * idyh
            dbtot1_dy(ny, :) = (3.0*btot1(:, ny) - 4.0*btot1(:, ny1) + btot1(:, ny2)) * idyh
        endif
    end subroutine calc_fields_gradients
end module mhd_data_sli

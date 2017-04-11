!*******************************************************************************
!< Module of MHD data for Shengtai Li's MHD code
!*******************************************************************************
module mhd_data_sli
    use constants, only: fp, dp
    implicit none
    private
    public mhd_config, fields, gradf
    public read_mhd_config, read_mhd_config_from_outfile, init_mhd_data, &
           free_mhd_data, read_mhd_data, broadcast_mhd_config, &
           init_fields_gradients, free_fields_gradients, calc_fields_gradients, &
           interp_fields, copy_fields

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

    type mhd_fields
        real(dp) :: vx, vy, bx, by, btot
    end type mhd_fields

    type fields_gradients
        real(dp) :: dvx_dx, dvy_dy
        real(dp) :: dbx_dx, dbx_dy
        real(dp) :: dby_dx, dby_dy
        real(dp) :: dbtot_dx, dbtot_dy
    end type fields_gradients

    type(mhd_configuration) :: mhd_config
    type(file_header) :: fheader
    type(mhd_fields) :: fields, fields1
    type(fields_gradients) :: gradf, gradf1

    type fields_type
        real(fp) :: vx, vy, vz
        real(fp) :: bx, by, bz, btot
        real(fp) :: vx1, vy1, vz1
        real(fp) :: bx1, by1, bz1, btot1
        real(fp), dimension(2) :: pad1  ! For 64-byte alignment and future expansion
    end type fields_type

    type fields_gradients_type
        real(fp) :: dvx_dx, dvy_dy
        real(fp) :: dbx_dx, dbx_dy
        real(fp) :: dby_dx, dby_dy
        real(fp) :: dbtot_dx, dbtot_dy
        real(fp) :: dvx1_dx, dvy1_dy
        real(fp) :: dbx1_dx, dbx1_dy
        real(fp) :: dby1_dx, dby1_dy
        real(fp) :: dbtot1_dx, dbtot1_dy
    end type fields_gradients_type

    type(fields_type), allocatable, dimension(:, :) :: f_array
    type(fields_gradients_type), allocatable, dimension(:, :) :: fgrad_array
    !dir$ attributes align:64 :: f_array
    !dir$ attributes align:64 :: fgrad_array

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
        mhd_config%dx = mhd_config%lx / (mhd_config%nx - 1)
        mhd_config%dy = mhd_config%ly / (mhd_config%ny - 1)
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
        mhd_config%dx = mhd_config%lx / (mhd_config%nx - 1)
        mhd_config%dy = mhd_config%ly / (mhd_config%ny - 1)
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

        allocate(f_array(nx, ny))

        f_array%vx = 0.0
        f_array%vy = 0.0
        f_array%vz = 0.0
        f_array%bx = 0.0
        f_array%by = 0.0
        f_array%bz = 0.0
        f_array%btot = 0.0
        f_array%vx1 = 0.0
        f_array%vy1 = 0.0
        f_array%vz1 = 0.0
        f_array%bx1 = 0.0
        f_array%by1 = 0.0
        f_array%bz1 = 0.0
        f_array%btot1 = 0.0

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

        allocate(fgrad_array(nx, ny))

        fgrad_array%dvx_dx = 0.0
        fgrad_array%dvy_dy = 0.0
        fgrad_array%dbx_dx = 0.0
        fgrad_array%dbx_dy = 0.0
        fgrad_array%dby_dx = 0.0
        fgrad_array%dby_dy = 0.0
        fgrad_array%dbtot_dx = 0.0
        fgrad_array%dbtot_dy = 0.0
        fgrad_array%dvx1_dx = 0.0
        fgrad_array%dvy1_dy = 0.0
        fgrad_array%dbx1_dx = 0.0
        fgrad_array%dbx1_dy = 0.0
        fgrad_array%dby1_dx = 0.0
        fgrad_array%dby1_dy = 0.0
        fgrad_array%dbtot1_dx = 0.0
        fgrad_array%dbtot1_dy = 0.0
    end subroutine init_fields_gradients

    !---------------------------------------------------------------------------
    !< Free MHD data arrays
    !---------------------------------------------------------------------------
    subroutine free_mhd_data
        implicit none
        deallocate(f_array)
        deallocate(mhd_data_single_core)
    end subroutine free_mhd_data

    !---------------------------------------------------------------------------
    !< Free the gradients of the MHD data arrays.
    !---------------------------------------------------------------------------
    subroutine free_fields_gradients
        implicit none
        deallocate(fgrad_array)
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
            if (mhd_config%nvar .eq. 7) then
                if (var_flag == 0) then
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%vx = mhd_data_single_core(:, :, 3)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%vy = mhd_data_single_core(:, :, 4)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%bx = mhd_data_single_core(:, :, 5)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%by = mhd_data_single_core(:, :, 6)
                else
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%vx1 = mhd_data_single_core(:, :, 3)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%vy1 = mhd_data_single_core(:, :, 4)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%bx1 = mhd_data_single_core(:, :, 5)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%by1 = mhd_data_single_core(:, :, 6)
                endif
            else
                if (var_flag == 0) then
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%vx = mhd_data_single_core(:, :, 3)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%vy = mhd_data_single_core(:, :, 4)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%vz = mhd_data_single_core(:, :, 5)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%bx = mhd_data_single_core(:, :, 6)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%by = mhd_data_single_core(:, :, 7)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%bz = mhd_data_single_core(:, :, 8)
                else
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%vx1 = mhd_data_single_core(:, :, 3)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%vy1 = mhd_data_single_core(:, :, 4)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%vz1 = mhd_data_single_core(:, :, 5)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%bx1 = mhd_data_single_core(:, :, 6)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%by1 = mhd_data_single_core(:, :, 7)
                    f_array(ix+1:ix+nx1, iy+1:iy+ny1)%bz1 = mhd_data_single_core(:, :, 8)
                endif
            endif
        enddo

        !< Calculate the magnitude of the magnetic field
        if (mhd_config%nvar .eq. 7) then
            if (var_flag == 0) then
                f_array%btot = sqrt(f_array%bx**2 + f_array%by**2)
            else
                f_array%btot1 = sqrt(f_array%bx1**2 + f_array%by1**2)
            endif
        else
            if (var_flag == 0) then
                f_array%btot = sqrt(f_array%bx**2 + f_array%by**2 + f_array%bz**2)
            else
                f_array%btot1 = sqrt(f_array%bx1**2 + f_array%by1**2 + f_array%bz1**2)
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
            fgrad_array(2:nx1, :)%dvx_dx = (f_array(3:nx, :)%vx - f_array(1:nx2, :)%vx) * idxh
            fgrad_array(1, :)%dvx_dx = (-3.0*f_array(1, :)%vx + 4.0*f_array(2, :)%vx - f_array(3, :)%vx) * idxh
            fgrad_array(nx, :)%dvx_dx = (3.0*f_array(nx, :)%vx - 4.0*f_array(nx1, :)%vx + f_array(nx2, :)%vx) * idxh
            fgrad_array(:, 2:ny1)%dvy_dy = (f_array(:, 3:ny)%vy - f_array(:, 1:ny2)%vy) * idyh
            fgrad_array(:, 1)%dvy_dy = (-3.0*f_array(:, 1)%vy + 4.0*f_array(:, 2)%vy - f_array(:, 3)%vy) * idyh
            fgrad_array(:, ny)%dvy_dy = (3.0*f_array(:, ny)%vy - 4.0*f_array(:, ny1)%vy + f_array(:, ny2)%vy) * idyh

            fgrad_array(2:nx1, :)%dbx_dx = (f_array(3:nx, :)%bx - f_array(1:nx2, :)%bx) * idxh
            fgrad_array(1, :)%dbx_dx = (-3.0*f_array(1, :)%bx + 4.0*f_array(2, :)%bx - f_array(3, :)%bx) * idxh
            fgrad_array(nx, :)%dbx_dx = (3.0*f_array(nx, :)%bx - 4.0*f_array(nx1, :)%bx + f_array(nx2, :)%bx) * idxh
            fgrad_array(:, 2:ny1)%dbx_dy = (f_array(:, 3:ny)%bx-f_array(:, 1:ny2)%bx) * idyh
            fgrad_array(:, 1)%dbx_dy = (-3.0*f_array(:, 1)%bx + 4.0*f_array(:, 2)%bx - f_array(:, 3)%bx) * idyh
            fgrad_array(:, ny)%dbx_dy = (3.0*f_array(:, ny)%bx - 4.0*f_array(:, ny1)%bx + f_array(:, ny2)%bx) * idyh

            fgrad_array(2:nx1, :)%dby_dx = (f_array(3:nx, :)%by - f_array(1:nx2, :)%by) * idxh
            fgrad_array(1, :)%dby_dx = (-3.0*f_array(1, :)%by + 4.0*f_array(2, :)%by - f_array(3, :)%by) * idxh
            fgrad_array(nx, :)%dby_dx = (3.0*f_array(nx, :)%by - 4.0*f_array(nx1, :)%by + f_array(nx2, :)%by) * idxh
            fgrad_array(:, 2:ny1)%dby_dy = (f_array(:, 3:ny)%by - f_array(:, 1:ny2)%by) * idyh
            fgrad_array(:, 1)%dby_dy = (-3.0*f_array(:, 1)%by + 4.0*f_array(:, 2)%by - f_array(:, 3)%by) * idyh
            fgrad_array(:, ny)%dby_dy = (3.0*f_array(:, ny)%by - 4.0*f_array(:, ny1)%by + f_array(:, ny2)%by) * idyh

            fgrad_array(2:nx1, :)%dbtot_dx = (f_array(3:nx, :)%btot - f_array(1:nx2, :)%btot) * idxh
            fgrad_array(1, :)%dbtot_dx = (-3.0*f_array(1, :)%btot + 4.0*f_array(2, :)%btot - f_array(3, :)%btot) * idxh
            fgrad_array(nx, :)%dbtot_dx = (3.0*f_array(nx, :)%btot - 4.0*f_array(nx1, :)%btot + f_array(nx2, :)%btot) * idxh
            fgrad_array(:, 2:ny1)%dbtot_dy = (f_array(:, 3:ny)%btot-f_array(:, 1:ny2)%btot) * idyh
            fgrad_array(:, 1)%dbtot_dy = (-3.0*f_array(:, 1)%btot + 4.0*f_array(:, 2)%btot - f_array(:, 3)%btot) * idyh
            fgrad_array(:, ny)%dbtot_dy = (3.0*f_array(:, ny)%btot - 4.0*f_array(:, ny1)%btot + f_array(:, ny2)%btot) * idyh
        else
            fgrad_array(2:nx1, :)%dvx1_dx = (f_array(3:nx, :)%vx1 - f_array(1:nx2, :)%vx1) * idxh
            fgrad_array(1, :)%dvx1_dx = (-3.0*f_array(1, :)%vx1 + 4.0*f_array(2, :)%vx1 - f_array(3, :)%vx1) * idxh
            fgrad_array(nx, :)%dvx1_dx = (3.0*f_array(nx, :)%vx1 - 4.0*f_array(nx1, :)%vx1 + f_array(nx2, :)%vx1) * idxh
            fgrad_array(:, 2:ny1)%dvy1_dy = (f_array(:, 3:ny)%vy1-f_array(:, 1:ny2)%vy1) * idyh
            fgrad_array(:, 1)%dvy1_dy = (-3.0*f_array(:, 1)%vy1 + 4.0*f_array(:, 2)%vy1 - f_array(:, 3)%vy1) * idyh
            fgrad_array(:, ny)%dvy1_dy = (3.0*f_array(:, ny)%vy1 - 4.0*f_array(:, ny1)%vy1 + f_array(:, ny2)%vy1) * idyh

            fgrad_array(2:nx1, :)%dbx1_dx = (f_array(3:nx, :)%bx1 - f_array(1:nx2, :)%bx1) * idxh
            fgrad_array(1, :)%dbx1_dx = (-3.0*f_array(1, :)%bx1 + 4.0*f_array(2, :)%bx1 - f_array(3, :)%bx1) * idxh
            fgrad_array(nx, :)%dbx1_dx = (3.0*f_array(nx, :)%bx1 - 4.0*f_array(nx1, :)%bx1 + f_array(nx2, :)%bx1) * idxh
            fgrad_array(:, 2:ny1)%dbx1_dy = (f_array(:, 3:ny)%bx1-f_array(:, 1:ny2)%bx1) * idyh
            fgrad_array(:, 1)%dbx1_dy = (-3.0*f_array(:, 1)%bx1 + 4.0*f_array(:, 2)%bx1 - f_array(:, 3)%bx1) * idyh
            fgrad_array(:, ny)%dbx1_dy = (3.0*f_array(:, ny)%bx1 - 4.0*f_array(:, ny1)%bx1 + f_array(:, ny2)%bx1) * idyh

            fgrad_array(2:nx1, :)%dby1_dx = (f_array(3:nx, :)%by1 - f_array(1:nx2, :)%by1) * idxh
            fgrad_array(1, :)%dby1_dx = (-3.0*f_array(1, :)%by1 + 4.0*f_array(2, :)%by1 - f_array(3, :)%by1) * idxh
            fgrad_array(nx, :)%dby1_dx = (3.0*f_array(nx, :)%by1 - 4.0*f_array(nx1, :)%by1 + f_array(nx2, :)%by1) * idxh
            fgrad_array(:, 2:ny1)%dby1_dy = (f_array(:, 3:ny)%by1-f_array(:, 1:ny2)%by1) * idyh
            fgrad_array(:, 1)%dby1_dy = (-3.0*f_array(:, 1)%by1 + 4.0*f_array(:, 2)%by1 - f_array(:, 3)%by1) * idyh
            fgrad_array(:, ny)%dby1_dy = (3.0*f_array(:, ny)%by1 - 4.0*f_array(:, ny1)%by1 + f_array(:, ny2)%by1) * idyh

            fgrad_array(2:nx1, :)%dbtot1_dx = (f_array(3:nx, :)%btot1 - f_array(1:nx2, :)%btot1) * idxh
            fgrad_array(1, :)%dbtot1_dx = (-3.0*f_array(1, :)%btot1 + 4.0*f_array(2, :)%btot1 - f_array(3, :)%btot1) * idxh
            fgrad_array(nx, :)%dbtot1_dx = (3.0*f_array(nx, :)%btot1 - 4.0*f_array(nx1, :)%btot1 + f_array(nx2, :)%btot1) * idxh
            fgrad_array(:, 2:ny1)%dbtot1_dy = (f_array(:, 3:ny)%btot1-f_array(:, 1:ny2)%btot1) * idyh
            fgrad_array(:, 1)%dbtot1_dy = (-3.0*f_array(:, 1)%btot1 + 4.0*f_array(:, 2)%btot1 - f_array(:, 3)%btot1) * idyh
            fgrad_array(:, ny)%dbtot1_dy = (3.0*f_array(:, ny)%btot1 - 4.0*f_array(:, ny1)%btot1 + f_array(:, ny2)%btot1) * idyh
        endif
    end subroutine calc_fields_gradients

    !---------------------------------------------------------------------------
    !< Interpolate the MHD fields and their gradients on one position
    !< Args:
    !<  ix, iy: the lower-left corner of the grid.
    !<  rx, ry: the offset to the lower-left corner of the grid where the
    !<          the position is. They are normalized to the grid sizes.
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !---------------------------------------------------------------------------
    subroutine interp_fields(ix, iy, rx, ry, rt)
        implicit none
        real(dp), intent(in) :: rx, ry, rt
        integer, intent(in) :: ix, iy
        real(dp) :: rx1, ry1, rt1, w1, w2, w3, w4
        integer :: ix1, iy1
        rx1 = 1.0_dp - rx
        ry1 = 1.0_dp - ry
        rt1 = 1.0_dp - rt
        ix1 = ix + 1
        iy1 = iy + 1
        w1 = rx1 * ry1
        w2 = rx * ry1
        w3 = rx1 * ry
        w4 = rx * ry

        fields%vx    = f_array(ix, iy)%vx    * w1
        fields%vy    = f_array(ix, iy)%vy    * w1
        fields%bx    = f_array(ix, iy)%bx    * w1
        fields%by    = f_array(ix, iy)%by    * w1
        fields%btot  = f_array(ix, iy)%btot  * w1
        fields1%vx   = f_array(ix, iy)%vx1   * w1
        fields1%vy   = f_array(ix, iy)%vy1   * w1
        fields1%bx   = f_array(ix, iy)%bx1   * w1
        fields1%by   = f_array(ix, iy)%by1   * w1
        fields1%btot = f_array(ix, iy)%btot1 * w1

        fields%vx    = fields%vx     + f_array(ix1, iy)%vx    * w2
        fields%vy    = fields%vy     + f_array(ix1, iy)%vy    * w2
        fields%bx    = fields%bx     + f_array(ix1, iy)%bx    * w2
        fields%by    = fields%by     + f_array(ix1, iy)%by    * w2
        fields%btot  = fields%btot   + f_array(ix1, iy)%btot  * w2
        fields1%vx   = fields1%vx    + f_array(ix1, iy)%vx1   * w2
        fields1%vy   = fields1%vy    + f_array(ix1, iy)%vy1   * w2
        fields1%bx   = fields1%bx    + f_array(ix1, iy)%bx1   * w2
        fields1%by   = fields1%by    + f_array(ix1, iy)%by1   * w2
        fields1%btot = fields1%btot  + f_array(ix1, iy)%btot1 * w2

        fields%vx    = fields%vx     + f_array(ix, iy1)%vx    * w3
        fields%vy    = fields%vy     + f_array(ix, iy1)%vy    * w3
        fields%bx    = fields%bx     + f_array(ix, iy1)%bx    * w3
        fields%by    = fields%by     + f_array(ix, iy1)%by    * w3
        fields%btot  = fields%btot   + f_array(ix, iy1)%btot  * w3
        fields1%vx   = fields1%vx    + f_array(ix, iy1)%vx1   * w3
        fields1%vy   = fields1%vy    + f_array(ix, iy1)%vy1   * w3
        fields1%bx   = fields1%bx    + f_array(ix, iy1)%bx1   * w3
        fields1%by   = fields1%by    + f_array(ix, iy1)%by1   * w3
        fields1%btot = fields1%btot  + f_array(ix, iy1)%btot1 * w3

        fields%vx    = fields%vx     + f_array(ix1, iy1)%vx    * w4
        fields%vy    = fields%vy     + f_array(ix1, iy1)%vy    * w4
        fields%bx    = fields%bx     + f_array(ix1, iy1)%bx    * w4
        fields%by    = fields%by     + f_array(ix1, iy1)%by    * w4
        fields%btot  = fields%btot   + f_array(ix1, iy1)%btot  * w4
        fields1%vx   = fields1%vx    + f_array(ix1, iy1)%vx1   * w4
        fields1%vy   = fields1%vy    + f_array(ix1, iy1)%vy1   * w4
        fields1%bx   = fields1%bx    + f_array(ix1, iy1)%bx1   * w4
        fields1%by   = fields1%by    + f_array(ix1, iy1)%by1   * w4
        fields1%btot = fields1%btot  + f_array(ix1, iy1)%btot1 * w4

        !< Time interpolation
        fields%vx   = fields1%vx * rt1 + fields%vx * rt
        fields%vy   = fields1%vy * rt1 + fields%vy * rt
        fields%bx   = fields1%bx * rt1 + fields%bx * rt
        fields%by   = fields1%by * rt1 + fields%by * rt
        fields%btot = fields1%btot * rt1 + fields%btot * rt

        gradf%dvx_dx    = fgrad_array(ix, iy)%dvx_dx    * w1
        gradf%dvy_dy    = fgrad_array(ix, iy)%dvy_dy    * w1
        gradf%dbx_dx    = fgrad_array(ix, iy)%dbx_dx    * w1
        gradf%dbx_dy    = fgrad_array(ix, iy)%dbx_dy    * w1
        gradf%dby_dx    = fgrad_array(ix, iy)%dby_dx    * w1
        gradf%dby_dy    = fgrad_array(ix, iy)%dby_dy    * w1
        gradf%dbtot_dx  = fgrad_array(ix, iy)%dbtot_dx  * w1
        gradf%dbtot_dy  = fgrad_array(ix, iy)%dbtot_dy  * w1
        gradf1%dvx_dx   = fgrad_array(ix, iy)%dvx1_dx   * w1
        gradf1%dvy_dy   = fgrad_array(ix, iy)%dvy1_dy   * w1
        gradf1%dbx_dx   = fgrad_array(ix, iy)%dbx1_dx   * w1
        gradf1%dbx_dy   = fgrad_array(ix, iy)%dbx1_dy   * w1
        gradf1%dby_dx   = fgrad_array(ix, iy)%dby1_dx   * w1
        gradf1%dby_dy   = fgrad_array(ix, iy)%dby1_dy   * w1
        gradf1%dbtot_dx = fgrad_array(ix, iy)%dbtot1_dx * w1
        gradf1%dbtot_dy = fgrad_array(ix, iy)%dbtot1_dy * w1

        gradf%dvx_dx    = gradf%dvx_dx    + fgrad_array(ix1, iy)%dvx_dx    * w2
        gradf%dvy_dy    = gradf%dvy_dy    + fgrad_array(ix1, iy)%dvy_dy    * w2
        gradf%dbx_dx    = gradf%dbx_dx    + fgrad_array(ix1, iy)%dbx_dx    * w2
        gradf%dbx_dy    = gradf%dbx_dy    + fgrad_array(ix1, iy)%dbx_dy    * w2
        gradf%dby_dx    = gradf%dby_dx    + fgrad_array(ix1, iy)%dby_dx    * w2
        gradf%dby_dy    = gradf%dby_dy    + fgrad_array(ix1, iy)%dby_dy    * w2
        gradf%dbtot_dx  = gradf%dbtot_dx  + fgrad_array(ix1, iy)%dbtot_dx  * w2
        gradf%dbtot_dy  = gradf%dbtot_dy  + fgrad_array(ix1, iy)%dbtot_dy  * w2
        gradf1%dvx_dx   = gradf1%dvx_dx   + fgrad_array(ix1, iy)%dvx1_dx   * w2
        gradf1%dvy_dy   = gradf1%dvy_dy   + fgrad_array(ix1, iy)%dvy1_dy   * w2
        gradf1%dbx_dx   = gradf1%dbx_dx   + fgrad_array(ix1, iy)%dbx1_dx   * w2
        gradf1%dbx_dy   = gradf1%dbx_dy   + fgrad_array(ix1, iy)%dbx1_dy   * w2
        gradf1%dby_dx   = gradf1%dby_dx   + fgrad_array(ix1, iy)%dby1_dx   * w2
        gradf1%dby_dy   = gradf1%dby_dy   + fgrad_array(ix1, iy)%dby1_dy   * w2
        gradf1%dbtot_dx = gradf1%dbtot_dx + fgrad_array(ix1, iy)%dbtot1_dx * w2
        gradf1%dbtot_dy = gradf1%dbtot_dy + fgrad_array(ix1, iy)%dbtot1_dy * w2

        gradf%dvx_dx    = gradf%dvx_dx    + fgrad_array(ix, iy1)%dvx_dx    * w3
        gradf%dvy_dy    = gradf%dvy_dy    + fgrad_array(ix, iy1)%dvy_dy    * w3
        gradf%dbx_dx    = gradf%dbx_dx    + fgrad_array(ix, iy1)%dbx_dx    * w3
        gradf%dbx_dy    = gradf%dbx_dy    + fgrad_array(ix, iy1)%dbx_dy    * w3
        gradf%dby_dx    = gradf%dby_dx    + fgrad_array(ix, iy1)%dby_dx    * w3
        gradf%dby_dy    = gradf%dby_dy    + fgrad_array(ix, iy1)%dby_dy    * w3
        gradf%dbtot_dx  = gradf%dbtot_dx  + fgrad_array(ix, iy1)%dbtot_dx  * w3
        gradf%dbtot_dy  = gradf%dbtot_dy  + fgrad_array(ix, iy1)%dbtot_dy  * w3
        gradf1%dvx_dx   = gradf1%dvx_dx   + fgrad_array(ix, iy1)%dvx1_dx   * w3
        gradf1%dvy_dy   = gradf1%dvy_dy   + fgrad_array(ix, iy1)%dvy1_dy   * w3
        gradf1%dbx_dx   = gradf1%dbx_dx   + fgrad_array(ix, iy1)%dbx1_dx   * w3
        gradf1%dbx_dy   = gradf1%dbx_dy   + fgrad_array(ix, iy1)%dbx1_dy   * w3
        gradf1%dby_dx   = gradf1%dby_dx   + fgrad_array(ix, iy1)%dby1_dx   * w3
        gradf1%dby_dy   = gradf1%dby_dy   + fgrad_array(ix, iy1)%dby1_dy   * w3
        gradf1%dbtot_dx = gradf1%dbtot_dx + fgrad_array(ix, iy1)%dbtot1_dx * w3
        gradf1%dbtot_dy = gradf1%dbtot_dy + fgrad_array(ix, iy1)%dbtot1_dy * w3

        gradf%dvx_dx    = gradf%dvx_dx    + fgrad_array(ix1, iy1)%dvx_dx    * w4
        gradf%dvy_dy    = gradf%dvy_dy    + fgrad_array(ix1, iy1)%dvy_dy    * w4
        gradf%dbx_dx    = gradf%dbx_dx    + fgrad_array(ix1, iy1)%dbx_dx    * w4
        gradf%dbx_dy    = gradf%dbx_dy    + fgrad_array(ix1, iy1)%dbx_dy    * w4
        gradf%dby_dx    = gradf%dby_dx    + fgrad_array(ix1, iy1)%dby_dx    * w4
        gradf%dby_dy    = gradf%dby_dy    + fgrad_array(ix1, iy1)%dby_dy    * w4
        gradf%dbtot_dx  = gradf%dbtot_dx  + fgrad_array(ix1, iy1)%dbtot_dx  * w4
        gradf%dbtot_dy  = gradf%dbtot_dy  + fgrad_array(ix1, iy1)%dbtot_dy  * w4
        gradf1%dvx_dx   = gradf1%dvx_dx   + fgrad_array(ix1, iy1)%dvx1_dx   * w4
        gradf1%dvy_dy   = gradf1%dvy_dy   + fgrad_array(ix1, iy1)%dvy1_dy   * w4
        gradf1%dbx_dx   = gradf1%dbx_dx   + fgrad_array(ix1, iy1)%dbx1_dx   * w4
        gradf1%dbx_dy   = gradf1%dbx_dy   + fgrad_array(ix1, iy1)%dbx1_dy   * w4
        gradf1%dby_dx   = gradf1%dby_dx   + fgrad_array(ix1, iy1)%dby1_dx   * w4
        gradf1%dby_dy   = gradf1%dby_dy   + fgrad_array(ix1, iy1)%dby1_dy   * w4
        gradf1%dbtot_dx = gradf1%dbtot_dx + fgrad_array(ix1, iy1)%dbtot1_dx * w4
        gradf1%dbtot_dy = gradf1%dbtot_dy + fgrad_array(ix1, iy1)%dbtot1_dy * w4

        !< Time interpolation
        gradf%dvx_dx   = gradf1%dvx_dx * rt1 + gradf%dvx_dx * rt
        gradf%dvy_dy   = gradf1%dvy_dy * rt1 + gradf%dvy_dy * rt
        gradf%dbx_dx   = gradf1%dbx_dx * rt1 + gradf%dbx_dx * rt
        gradf%dbx_dy   = gradf1%dbx_dy * rt1 + gradf%dbx_dy * rt
        gradf%dby_dx   = gradf1%dby_dx * rt1 + gradf%dby_dx * rt
        gradf%dby_dy   = gradf1%dby_dy * rt1 + gradf%dby_dy * rt
        gradf%dbtot_dx = gradf1%dbtot_dx * rt1 + gradf%dbtot_dx * rt
        gradf%dbtot_dy = gradf1%dbtot_dy * rt1 + gradf%dbtot_dy * rt
    end subroutine interp_fields

    !---------------------------------------------------------------------------
    ! Copy fields for usage in the next time interval
    !---------------------------------------------------------------------------
    subroutine copy_fields
        implicit none
        f_array%vx1 = f_array%vx
        f_array%vy1 = f_array%vy
        f_array%bx1 = f_array%bx
        f_array%by1 = f_array%by
        f_array%btot1 = f_array%btot
        if (mhd_config%nvar .gt. 7) then
            f_array%vz1 = f_array%vz
            f_array%bz1 = f_array%bz
        endif

        fgrad_array%dvx1_dx = fgrad_array%dvx_dx
        fgrad_array%dvy1_dy = fgrad_array%dvy_dy
        fgrad_array%dbx1_dx = fgrad_array%dbx_dx
        fgrad_array%dby1_dx = fgrad_array%dby_dx
        fgrad_array%dbx1_dy = fgrad_array%dbx_dy
        fgrad_array%dby1_dy = fgrad_array%dby_dy
        fgrad_array%dbtot1_dx = fgrad_array%dbtot_dx
        fgrad_array%dbtot1_dy = fgrad_array%dbtot_dy
    end subroutine copy_fields
end module mhd_data_sli

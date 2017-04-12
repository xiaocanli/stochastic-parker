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
           interp_fields, copy_fields, save_mhd_config

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
        real(dp) :: pad1, pad2, pad3
    end type mhd_fields

    type fields_gradients
        real(dp) :: dvx_dx, dvy_dy
        real(dp) :: dbx_dx, dbx_dy
        real(dp) :: dby_dx, dby_dy
        real(dp) :: dbtot_dx, dbtot_dy
    end type fields_gradients

    type(mhd_configuration) :: mhd_config
    type(file_header) :: fheader
    ! type(mhd_fields) :: fields, fields1
    ! type(fields_gradients) :: gradf, gradf1
    real(dp), dimension(8) :: fields, fields1
    real(dp), dimension(8) :: gradf, gradf1

    type fields_type
        real(fp) :: vx, vy, vz
        real(fp) :: bx, by, bz, btot
        real(fp) :: pad1
        real(fp) :: vx1, vy1, vz1
        real(fp) :: bx1, by1, bz1, btot1
        real(fp) :: pad2  ! For 64-byte alignment
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

    ! type(fields_type), allocatable, dimension(:, :) :: f_array
    ! type(fields_gradients_type), allocatable, dimension(:, :) :: fgrad_array
    real(fp), allocatable, dimension(:, :, :) :: f_array
    real(fp), allocatable, dimension(:, :, :) :: fgrad_array
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

        ! vx, vy, vz, pad1, bx, by, bz, btot, vx1, vy1, vz1, pad2, bx1, by1, bz1, btot1
        allocate(f_array(16, nx, ny))
        f_array = 0.0

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

        ! dvx_dx, dvy_dy dbx_dx, dbx_dy dby_dx, dby_dy dbtot_dx, dbtot_dy
        ! dvx1_dx, dvy1_dy dbx1_dx, dbx1_dy dby1_dx, dby1_dy dbtot1_dx, dbtot1_dy
        allocate(fgrad_array(16, nx, ny))
        fgrad_array = 0.0
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
                    f_array(1, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 3)
                    f_array(2, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 4)
                    f_array(5, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 5)
                    f_array(6, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 6)
                else
                    f_array(9,  ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 3)
                    f_array(10, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 4)
                    f_array(13, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 5)
                    f_array(14, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 6)
                endif
            else
                if (var_flag == 0) then
                    f_array(1, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 3)
                    f_array(2, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 4)
                    f_array(3, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 5)
                    f_array(5, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 6)
                    f_array(6, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 7)
                    f_array(7, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 8)
                else
                    f_array(9,  ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 3)
                    f_array(10, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 4)
                    f_array(11, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 5)
                    f_array(13, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 6)
                    f_array(14, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 7)
                    f_array(15, ix+1:ix+nx1, iy+1:iy+ny1) = mhd_data_single_core(:, :, 8)
                endif
            endif
        enddo

        !< Calculate the magnitude of the magnetic field
        if (mhd_config%nvar .eq. 7) then
            if (var_flag == 0) then
                f_array(8, :, :) = sqrt(f_array(5, :, :)**2 + f_array(6, :, :)**2)
            else
                f_array(16, :, :) = sqrt(f_array(13, :, :)**2 + f_array(14, :, :)**2)
            endif
        else
            if (var_flag == 0) then
                f_array(8, :, :) = sqrt(f_array(5, :, :)**2 + &
                                        f_array(6, :, :)**2 + &
                                        f_array(7, :, :)**2)
            else
                f_array(16, :, :) = sqrt(f_array(13, :, :)**2 + &
                                         f_array(14, :, :)**2 + &
                                         f_array(15, :, :)**2)
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
            fgrad_array(1, 2:nx1, :) =  (f_array(1, 3:nx, :) - &
                                         f_array(1, 1:nx2, :)) * idxh
            fgrad_array(1, 1, :) = (-3.0*f_array(1, 1, :) + &
                                     4.0*f_array(1, 2, :) - &
                                         f_array(1, 3, :)) * idxh
            fgrad_array(1, nx, :) = (3.0*f_array(1, nx, :) - &
                                     4.0*f_array(1, nx1, :) + &
                                         f_array(1, nx2, :)) * idxh

            fgrad_array(2, :, 2:ny1) =  (f_array(2, :, 3:ny) - &
                                         f_array(2, :, 1:ny2)) * idyh
            fgrad_array(2, :, 1) = (-3.0*f_array(2, :, 1) + &
                                     4.0*f_array(2, :, 2) - &
                                         f_array(2, :, 3)) * idyh
            fgrad_array(2, :, ny) = (3.0*f_array(2, :, ny) - &
                                     4.0*f_array(2, :, ny1) + &
                                         f_array(2, :, ny2)) * idyh

            fgrad_array(3, 2:nx1, :) =  (f_array(5, 3:nx, :) - &
                                         f_array(5, 1:nx2, :)) * idxh
            fgrad_array(3, 1, :) = (-3.0*f_array(5, 1, :) + &
                                     4.0*f_array(5, 2, :) - &
                                         f_array(5, 3, :)) * idxh
            fgrad_array(3, nx, :) = (3.0*f_array(5, nx, :) - &
                                     4.0*f_array(5, nx1, :) + &
                                         f_array(5, nx2, :)) * idxh

            fgrad_array(4, :, 2:ny1) =  (f_array(5, :, 3:ny) - &
                                         f_array(5, :, 1:ny2)) * idyh
            fgrad_array(4, :, 1) = (-3.0*f_array(5, :, 1) + &
                                     4.0*f_array(5, :, 2) - &
                                         f_array(5, :, 3)) * idyh
            fgrad_array(4, :, ny) = (3.0*f_array(5, :, ny) - &
                                     4.0*f_array(5, :, ny1) + &
                                         f_array(5, :, ny2)) * idyh

            fgrad_array(5, 2:nx1, :) =  (f_array(6, 3:nx, :) - &
                                         f_array(6, 1:nx2, :)) * idxh
            fgrad_array(5, 1, :) = (-3.0*f_array(6, 1, :) + &
                                     4.0*f_array(6, 2, :) - &
                                         f_array(6, 3, :)) * idxh
            fgrad_array(5, nx, :) = (3.0*f_array(6, nx, :) - &
                                     4.0*f_array(6, nx1, :) + &
                                         f_array(6, nx2, :)) * idxh

            fgrad_array(6, :, 2:ny1) =  (f_array(6, :, 3:ny) - &
                                         f_array(6, :, 1:ny2)) * idyh
            fgrad_array(6, :, 1) = (-3.0*f_array(6, :, 1) + &
                                     4.0*f_array(6, :, 2) - &
                                         f_array(6, :, 3)) * idyh
            fgrad_array(6, :, ny) = (3.0*f_array(6, :, ny) - &
                                     4.0*f_array(6, :, ny1) + &
                                         f_array(6, :, ny2)) * idyh

            fgrad_array(7, 2:nx1, :) =  (f_array(8, 3:nx, :) - &
                                         f_array(8, 1:nx2, :)) * idxh
            fgrad_array(7, 1, :) = (-3.0*f_array(8, 1, :) + &
                                     4.0*f_array(8, 2, :) - &
                                         f_array(8, 3, :)) * idxh
            fgrad_array(7, nx, :) = (3.0*f_array(8, nx, :) - &
                                     4.0*f_array(8, nx1, :) + &
                                         f_array(8, nx2, :)) * idxh

            fgrad_array(8, :, 2:ny1) =  (f_array(8, :, 3:ny) - &
                                         f_array(8, :, 1:ny2)) * idyh
            fgrad_array(8, :, 1) = (-3.0*f_array(8, :, 1) + &
                                     4.0*f_array(8, :, 2) - &
                                         f_array(8, :, 3)) * idyh
            fgrad_array(8, :, ny) = (3.0*f_array(8, :, ny) - &
                                     4.0*f_array(8, :, ny1) + &
                                         f_array(8, :, ny2)) * idyh
        else
            fgrad_array(9, 2:nx1, :) =  (f_array(9, 3:nx, :) - &
                                         f_array(9, 1:nx2, :)) * idxh
            fgrad_array(9, 1, :) = (-3.0*f_array(9, 1, :) + &
                                     4.0*f_array(9, 2, :) - &
                                         f_array(9, 3, :)) * idxh
            fgrad_array(9, nx, :) = (3.0*f_array(9, nx, :) - &
                                     4.0*f_array(9, nx1, :) + &
                                         f_array(9, nx2, :)) * idxh

            fgrad_array(10, :, 2:ny1) =  (f_array(10, :, 3:ny) - &
                                          f_array(10, :, 1:ny2)) * idyh
            fgrad_array(10, :, 1) = (-3.0*f_array(10, :, 1) + &
                                      4.0*f_array(10, :, 2) - &
                                          f_array(10, :, 3)) * idyh
            fgrad_array(10, :, ny) = (3.0*f_array(10, :, ny) - &
                                      4.0*f_array(10, :, ny1) + &
                                          f_array(10, :, ny2)) * idyh

            fgrad_array(11, 2:nx1, :) =  (f_array(13, 3:nx, :) - &
                                          f_array(13, 1:nx2, :)) * idxh
            fgrad_array(11, 1, :) = (-3.0*f_array(13, 1, :) + &
                                      4.0*f_array(13, 2, :) - &
                                          f_array(13, 3, :)) * idxh
            fgrad_array(11, nx, :) = (3.0*f_array(13, nx, :) - &
                                      4.0*f_array(13, nx1, :) + &
                                          f_array(13, nx2, :)) * idxh

            fgrad_array(12, :, 2:ny1) =  (f_array(13, :, 3:ny) - &
                                          f_array(13, :, 1:ny2)) * idyh
            fgrad_array(12, :, 1) = (-3.0*f_array(13, :, 1) + &
                                      4.0*f_array(13, :, 2) - &
                                          f_array(13, :, 3)) * idyh
            fgrad_array(12, :, ny) = (3.0*f_array(13, :, ny) - &
                                      4.0*f_array(13, :, ny1) + &
                                          f_array(13, :, ny2)) * idyh

            fgrad_array(13, 2:nx1, :) =  (f_array(14, 3:nx, :) - &
                                          f_array(14, 1:nx2, :)) * idxh
            fgrad_array(13, 1, :) = (-3.0*f_array(14, 1, :) + &
                                      4.0*f_array(14, 2, :) - &
                                          f_array(14, 3, :)) * idxh
            fgrad_array(13, nx, :) = (3.0*f_array(14, nx, :) - &
                                      4.0*f_array(14, nx1, :) + &
                                          f_array(14, nx2, :)) * idxh

            fgrad_array(14, :, 2:ny1) =  (f_array(14, :, 3:ny) - &
                                          f_array(14, :, 1:ny2)) * idyh
            fgrad_array(14, :, 1) = (-3.0*f_array(14, :, 1) + &
                                      4.0*f_array(14, :, 2) - &
                                          f_array(14, :, 3)) * idyh
            fgrad_array(14, :, ny) = (3.0*f_array(14, :, ny) - &
                                      4.0*f_array(14, :, ny1) + &
                                          f_array(14, :, ny2)) * idyh

            fgrad_array(15, 2:nx1, :) =  (f_array(16, 3:nx, :) - &
                                          f_array(16, 1:nx2, :)) * idxh
            fgrad_array(15, 1, :) = (-3.0*f_array(16, 1, :) + &
                                      4.0*f_array(16, 2, :) - &
                                          f_array(16, 3, :)) * idxh
            fgrad_array(15, nx, :) = (3.0*f_array(16, nx, :) - &
                                      4.0*f_array(16, nx1, :) + &
                                          f_array(16, nx2, :)) * idxh

            fgrad_array(16, :, 2:ny1) =  (f_array(16, :, 3:ny) - &
                                          f_array(16, :, 1:ny2)) * idyh
            fgrad_array(16, :, 1) = (-3.0*f_array(16, :, 1) + &
                                      4.0*f_array(16, :, 2) - &
                                          f_array(16, :, 3)) * idyh
            fgrad_array(16, :, ny) = (3.0*f_array(16, :, ny) - &
                                      4.0*f_array(16, :, ny1) + &
                                          f_array(16, :, ny2)) * idyh
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
        integer :: ix1, iy1, i
        rx1 = 1.0_dp - rx
        ry1 = 1.0_dp - ry
        rt1 = 1.0_dp - rt
        ix1 = ix + 1
        iy1 = iy + 1
        w1 = rx1 * ry1
        w2 = rx * ry1
        w3 = rx1 * ry
        w4 = rx * ry

        fields  = f_array(1:8,  ix, iy) * w1
        fields1 = f_array(9:16, ix, iy) * w1

        fields  = fields  + f_array(1:8,  ix1, iy) * w2
        fields1 = fields1 + f_array(9:16, ix1, iy) * w2

        fields  = fields  + f_array(1:8,  ix, iy1) * w3
        fields1 = fields1 + f_array(9:16, ix, iy1) * w3

        fields  = fields  + f_array(1:8,  ix1, iy1) * w4
        fields1 = fields1 + f_array(9:16, ix1, iy1) * w4

        !< Time interpolation
        fields = fields * rt1 + fields1 * rt

        gradf  = fgrad_array(1:8,  ix, iy) * w1
        gradf1 = fgrad_array(9:16, ix, iy) * w1

        gradf  = gradf  + fgrad_array(1:8,  ix1, iy) * w2
        gradf1 = gradf1 + fgrad_array(9:16, ix1, iy) * w2

        gradf  = gradf  + fgrad_array(1:8,  ix, iy1) * w3
        gradf1 = gradf1 + fgrad_array(9:16, ix, iy1) * w3

        gradf  = gradf  + fgrad_array(1:8,  ix1, iy1) * w4
        gradf1 = gradf1 + fgrad_array(9:16, ix1, iy1) * w4

        !< Time interpolation
        gradf = gradf * rt1 + gradf1 * rt
    end subroutine interp_fields

    !---------------------------------------------------------------------------
    ! Copy fields for usage in the next time interval
    !---------------------------------------------------------------------------
    subroutine copy_fields
        implicit none
        f_array(9,  :, :) = f_array(1, :, :)
        f_array(10, :, :) = f_array(2, :, :)
        if (mhd_config%nvar .gt. 7) then
            f_array(11, :, :) = f_array(3, :, :)
        endif
        f_array(12, :, :) = f_array(4, :, :)
        f_array(13, :, :) = f_array(5, :, :)
        if (mhd_config%nvar .gt. 7) then
            f_array(14, :, :) = f_array(6, :, :)
        endif
        f_array(15, :, :) = f_array(7, :, :)

        fgrad_array(9:16,  :, :) = fgrad_array(1:8, :, :)
    end subroutine copy_fields
end module mhd_data_sli

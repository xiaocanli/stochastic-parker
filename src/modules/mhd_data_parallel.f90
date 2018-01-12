!*******************************************************************************
!< Module for dealing with MHD data in parallel
!*******************************************************************************
module mhd_data_parallel
    use constants, only: fp, dp
    implicit none
    private
    save

    public init_field_data, free_field_data, read_field_data_parallel, &
           init_fields_gradients, free_fields_gradients, calc_fields_gradients, & 
           interp_fields, copy_fields, init_shock_xpos, free_shock_xpos, &
           locate_shock_xpos, interp_shock_location
    public fields, gradf

    real(dp), dimension(8) :: fields, fields1
    real(dp), dimension(16) :: gradf, gradf1
    !dir$ attributes align:64 :: fields
    !dir$ attributes align:64 :: fields1
    !dir$ attributes align:128 :: gradf
    !dir$ attributes align:128 :: gradf1

    real(fp), allocatable, dimension(:, :, :, :) :: f_array1, f_array2 ! Current,next
    real(fp), allocatable, dimension(:, :, :, :) :: fgrad_array1, fgrad_array2
    !dir$ attributes align:32 :: f_array1
    !dir$ attributes align:32 :: f_array2
    !dir$ attributes align:64 :: fgrad_array1
    !dir$ attributes align:64 :: fgrad_array2

    integer, allocatable, dimension(:, :) :: shock_xpos1, shock_xpos2  ! Shock x-position indices

    contains

    !---------------------------------------------------------------------------
    !< Initialize MHD field data arrays
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !<  nx, ny, nz: the dimensions of the data
    !<  ndim: number of actual dimension of the data. 1, 2 or 3
    !---------------------------------------------------------------------------
    subroutine init_field_data(interp_flag, nx, ny, nz, ndim)
        implicit none
        integer, intent(in) :: interp_flag, nx, ny, nz, ndim

        !< vx, vy, vz, pad1, bx, by, bz, btot
        if (ndim == 1) then
            allocate(f_array1(8, -1:nx+2, 1, 1))
        else if (ndim == 2) then
            allocate(f_array1(8, -1:nx+2, -1:ny+2, 1))
        else
            allocate(f_array1(8, -1:nx+2, -1:ny+2, -1:nz+2))
        endif
        f_array1 = 0.0
        ! Next time step
        if (interp_flag == 1) then
            if (ndim == 1) then
                allocate(f_array2(8, -1:nx+2, 1, 1))
            else if (ndim == 2) then
                allocate(f_array2(8, -1:nx+2, -1:ny+2, 1))
            else
                allocate(f_array2(8, -1:nx+2, -1:ny+2, -1:nz+2))
            endif
            f_array2 = 0.0
        endif
    end subroutine init_field_data

    !---------------------------------------------------------------------------
    !< Initialize the gradients of the MHD data arrays.
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !<  nx, ny, nz: the dimensions of the data
    !<  ndim: number of actual dimension of the data. 1, 2 or 3
    !---------------------------------------------------------------------------
    subroutine init_fields_gradients(interp_flag, nx, ny, nz, ndim)
        implicit none
        integer, intent(in) :: interp_flag, nx, ny, nz, ndim

        !< dvx_dx, dvy_dy, dvz_dz, dbx_dx, dbx_dy, dbx_dz
        !< dby_dx, dby_dy, dby_dz, dbz_dx, dbz_dy, dbz_dz
        !< dbtot_dx, dbtot_dy, dbtot_dz, pad1
        if (ndim == 1) then
            allocate(fgrad_array1(16, -1:nx+2, 1, 1))
        else if (ndim == 2) then
            allocate(fgrad_array1(16, -1:nx+2, -1:ny+2, 1))
        else
            allocate(fgrad_array1(16, -1:nx+2, -1:ny+2, -1:nz+2))
        endif
        fgrad_array1 = 0.0
        ! Next time step
        if (interp_flag == 1) then
            if (ndim == 1) then
                allocate(fgrad_array2(16, -1:nx+2, 1, 1))
            else if (ndim == 2) then
                allocate(fgrad_array2(16, -1:nx+2, -1:ny+2, 1))
            else
                allocate(fgrad_array2(16, -1:nx+2, -1:ny+2, -1:nz+2))
            endif
            fgrad_array2 = 0.0
        endif
    end subroutine init_fields_gradients

    !---------------------------------------------------------------------------
    !< Free MHD field data arrays
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine free_field_data(interp_flag)
        implicit none
        integer, intent(in) :: interp_flag
        deallocate(f_array1)
        if (interp_flag == 1) then
            deallocate(f_array2)
        endif
    end subroutine free_field_data

    !---------------------------------------------------------------------------
    !< Free the gradients of the MHD data arrays.
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine free_fields_gradients(interp_flag)
        implicit none
        integer, intent(in) :: interp_flag
        deallocate(fgrad_array1)
        if (interp_flag == 1) then
            deallocate(fgrad_array2)
        endif
    end subroutine free_fields_gradients

    !---------------------------------------------------------------------------
    !< Read MHD field data in parallel
    !< Args:
    !<  filename: file name to get the data
    !<  var_flag: indicating which set of variables to save the data. 0 for
    !<            f_array1 and other numbers for f_array2.
    !---------------------------------------------------------------------------
    subroutine read_field_data_parallel(filename, var_flag)
        use mpi_io_module, only: set_mpi_datatype_real, set_mpi_info, fileinfo, &
            open_data_mpi_io, read_data_mpi_io
        use simulation_setup_module, only: fconfig
        use mpi_module
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: var_flag
        integer :: mpi_datatype, fh
        integer, dimension(4) :: sizes, subsizes, starts
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
        sizes(1) = 8 ! vx, vy, vz, pad1, bx, by, bz, btot
        sizes(2) = fconfig%nxg
        sizes(3) = fconfig%nyg
        sizes(4) = fconfig%nzg
        subsizes(1) = 8
        subsizes(2) = fconfig%nxf
        subsizes(3) = fconfig%nyf
        subsizes(4) = fconfig%nzf
        starts(1) = 0
        starts(2) = fconfig%ix_min - 1
        starts(3) = fconfig%iy_min - 1
        starts(4) = fconfig%iz_min - 1
        fh = 11
        if (all(sizes == subsizes)) then
            if (mpi_rank == master) then
                write(*, "(A)") " Each MPI rank reads the whole data set."
            endif
            open(unit=fh, file=filename, access='stream', status='unknown', &
                 form='unformatted', action='read')
            if (var_flag == 0) then
                if (mpi_rank == master) then
                    read(fh, pos=1) f_array1
                endif
                call MPI_BCAST(f_array1, product(sizes), MPI_REAL4, master, MPI_COMM_WORLD, ierr)
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            else
                if (mpi_rank == master) then
                    read(fh, pos=1) f_array2
                endif
                call MPI_BCAST(f_array2, product(sizes), MPI_REAL4, master, MPI_COMM_WORLD, ierr)
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            endif
            close(fh)
        else
            mpi_datatype = set_mpi_datatype_real(sizes, subsizes, starts)
            call set_mpi_info
            call open_data_mpi_io(filename, MPI_MODE_RDONLY, fileinfo, fh)
            disp = 0
            offset = 0
            if (var_flag == 0) then
                call read_data_mpi_io(fh, mpi_datatype, subsizes, disp, offset, f_array1)
            else
                call read_data_mpi_io(fh, mpi_datatype, subsizes, disp, offset, f_array2)
            endif
            call MPI_FILE_CLOSE(fh, ierror)
        endif
        if (mpi_rank == master) then
            write(*, "(A)") "Finished reading MHD data."
        endif
    end subroutine read_field_data_parallel

    !---------------------------------------------------------------------------
    !< Calculate the gradients of the MHD data arrays.
    !< Args:
    !<  var_flag: indicating which set of variables.
    !<            0 for fgrad_array1 and other numbers for fgrad_array2.
    !---------------------------------------------------------------------------
    subroutine calc_fields_gradients(var_flag)
        use mhd_config_module, only: mhd_config
        implicit none
        integer, intent(in) :: var_flag
        real(dp) :: idxh, idyh
        integer :: unx, uny, lnx, lny
        integer :: unx1, unx2, uny1, uny2
        integer :: lnx1, lnx2, lny1, lny2
        idxh = 0.5_dp / mhd_config%dx
        idyh = 0.5_dp / mhd_config%dy
        unx = ubound(f_array1, 2)
        uny = ubound(f_array1, 3)
        lnx = lbound(f_array1, 2)
        lny = lbound(f_array1, 3)
        unx1 = unx - 1
        unx2 = unx - 2
        uny1 = uny - 1
        uny2 = uny - 2
        lnx1 = lnx + 1
        lnx2 = lnx + 2
        lny1 = lny + 1
        lny2 = lny + 2
        if (var_flag == 0) then
            fgrad_array1(1, lnx1:unx1, :, :) = (f_array1(1, lnx2:unx, :, :) - &
                                                f_array1(1, lnx:unx2, :, :)) * idxh
            fgrad_array1(1, lnx, :, :) =  (-3.0*f_array1(1, lnx, :, :) + &
                                            4.0*f_array1(1, lnx1, :, :) - &
                                                f_array1(1, lnx2, :, :)) * idxh
            fgrad_array1(1, unx, :, :) =   (3.0*f_array1(1, unx, :, :) - &
                                            4.0*f_array1(1, unx1, :, :) + &
                                                f_array1(1, unx2, :, :)) * idxh

            fgrad_array1(2, :, lny1:uny1, :) = (f_array1(2, :, lny2:uny, :) - &
                                                f_array1(2, :, lny:uny2, :)) * idyh
            fgrad_array1(2, :, lny, :) =  (-3.0*f_array1(2, :, lny, :) + &
                                            4.0*f_array1(2, :, lny1, :) - &
                                                f_array1(2, :, lny2, :)) * idyh
            fgrad_array1(2, :, uny, :) =   (3.0*f_array1(2, :, uny, :) - &
                                            4.0*f_array1(2, :, uny1, :) + &
                                                f_array1(2, :, uny2, :)) * idyh

            fgrad_array1(4, lnx1:unx1, :, :) = (f_array1(5, lnx2:unx, :, :) - &
                                                f_array1(5, lnx:unx2, :, :)) * idxh
            fgrad_array1(4, lnx, :, :) =  (-3.0*f_array1(5, lnx, :, :) + &
                                            4.0*f_array1(5, lnx1, :, :) - &
                                                f_array1(5, lnx2, :, :)) * idxh
            fgrad_array1(4, unx, :, :) =   (3.0*f_array1(5, unx, :, :) - &
                                            4.0*f_array1(5, unx1, :, :) + &
                                                f_array1(5, unx2, :, :)) * idxh

            fgrad_array1(5, :, lny1:uny1, :) = (f_array1(5, :, lny2:uny, :) - &
                                                f_array1(5, :, lny:uny2, :)) * idyh
            fgrad_array1(5, :, lny, :) =  (-3.0*f_array1(5, :, lny, :) + &
                                            4.0*f_array1(5, :, lny1, :) - &
                                                f_array1(5, :, lny2, :)) * idyh
            fgrad_array1(5, :, uny, :) =   (3.0*f_array1(5, :, uny, :) - &
                                            4.0*f_array1(5, :, uny1, :) + &
                                                f_array1(5, :, uny2, :)) * idyh

            fgrad_array1(7, lnx1:unx1, :, :) = (f_array1(6, lnx2:unx, :, :) - &
                                                f_array1(6, lnx:unx2, :, :)) * idxh
            fgrad_array1(7, lnx, :, :) =  (-3.0*f_array1(6, lnx, :, :) + &
                                            4.0*f_array1(6, lnx1, :, :) - &
                                                f_array1(6, lnx2, :, :)) * idxh
            fgrad_array1(7, unx, :, :) =   (3.0*f_array1(6, unx, :, :) - &
                                            4.0*f_array1(6, unx1, :, :) + &
                                                f_array1(6, unx2, :, :)) * idxh

            fgrad_array1(8, :, lny1:uny1, :) = (f_array1(6, :, lny2:uny, :) - &
                                                f_array1(6, :, lny:uny2, :)) * idyh
            fgrad_array1(8, :, lny, :) =  (-3.0*f_array1(6, :, lny, :) + &
                                            4.0*f_array1(6, :, lny1, :) - &
                                                f_array1(6, :, lny2, :)) * idyh
            fgrad_array1(8, :, uny, :) =   (3.0*f_array1(6, :, uny, :) - &
                                            4.0*f_array1(6, :, uny1, :) + &
                                                f_array1(6, :, uny2, :)) * idyh

            fgrad_array1(13, lnx1:unx1, :, :) = (f_array1(8, lnx2:unx, :, :) - &
                                                 f_array1(8, lnx:unx2, :, :)) * idxh
            fgrad_array1(13, lnx, :, :) =  (-3.0*f_array1(8, lnx, :, :) + &
                                             4.0*f_array1(8, lnx1, :, :) - &
                                                 f_array1(8, lnx2, :, :)) * idxh
            fgrad_array1(13, unx, :, :) =   (3.0*f_array1(8, unx, :, :) - &
                                             4.0*f_array1(8, unx1, :, :) + &
                                                 f_array1(8, unx2, :, :)) * idxh

            fgrad_array1(14, :, lny1:uny1, :) = (f_array1(8, :, lny2:uny, :) - &
                                                 f_array1(8, :, lny:uny2, :)) * idyh
            fgrad_array1(14, :, lny, :) =  (-3.0*f_array1(8, :, lny, :) + &
                                             4.0*f_array1(8, :, lny1, :) - &
                                                 f_array1(8, :, lny2, :)) * idyh
            fgrad_array1(14, :, uny, :) =   (3.0*f_array1(8, :, uny, :) - &
                                             4.0*f_array1(8, :, uny1, :) + &
                                                 f_array1(8, :, uny2, :)) * idyh
        else
            fgrad_array2(1, lnx1:unx1, :, :) = (f_array2(1, lnx2:unx, :, :) - &
                                                f_array2(1, lnx:unx2, :, :)) * idxh
            fgrad_array2(1, lnx, :, :) =  (-3.0*f_array2(1, lnx, :, :) + &
                                            4.0*f_array2(1, lnx1, :, :) - &
                                                f_array2(1, lnx2, :, :)) * idxh
            fgrad_array2(1, unx, :, :) =   (3.0*f_array2(1, unx, :, :) - &
                                            4.0*f_array2(1, unx1, :, :) + &
                                                f_array2(1, unx2, :, :)) * idxh

            fgrad_array2(2, :, lny1:uny1, :) = (f_array2(2, :, lny2:uny, :) - &
                                                f_array2(2, :, lny:uny2, :)) * idyh
            fgrad_array2(2, :, lny, :) =  (-3.0*f_array2(2, :, lny, :) + &
                                            4.0*f_array2(2, :, lny1, :) - &
                                                f_array2(2, :, lny2, :)) * idyh
            fgrad_array2(2, :, uny, :) =   (3.0*f_array2(2, :, uny, :) - &
                                            4.0*f_array2(2, :, uny1, :) + &
                                                f_array2(2, :, uny2, :)) * idyh

            fgrad_array2(4, lnx1:unx1, :, :) = (f_array2(5, lnx2:unx, :, :) - &
                                                f_array2(5, lnx:unx2, :, :)) * idxh
            fgrad_array2(4, lnx, :, :) =  (-3.0*f_array2(5, lnx, :, :) + &
                                            4.0*f_array2(5, lnx1, :, :) - &
                                                f_array2(5, lnx2, :, :)) * idxh
            fgrad_array2(4, unx, :, :) =   (3.0*f_array2(5, unx, :, :) - &
                                            4.0*f_array2(5, unx1, :, :) + &
                                                f_array2(5, unx2, :, :)) * idxh

            fgrad_array2(5, :, lny1:uny1, :) = (f_array2(5, :, lny2:uny, :) - &
                                                f_array2(5, :, lny:uny2, :)) * idyh
            fgrad_array2(5, :, lny, :) =  (-3.0*f_array2(5, :, lny, :) + &
                                            4.0*f_array2(5, :, lny1, :) - &
                                                f_array2(5, :, lny2, :)) * idyh
            fgrad_array2(5, :, uny, :) =   (3.0*f_array2(5, :, uny, :) - &
                                            4.0*f_array2(5, :, uny1, :) + &
                                                f_array2(5, :, uny2, :)) * idyh

            fgrad_array2(7, lnx1:unx1, :, :) = (f_array2(6, lnx2:unx, :, :) - &
                                                f_array2(6, lnx:unx2, :, :)) * idxh
            fgrad_array2(7, lnx, :, :) =  (-3.0*f_array2(6, lnx, :, :) + &
                                            4.0*f_array2(6, lnx1, :, :) - &
                                                f_array2(6, lnx2, :, :)) * idxh
            fgrad_array2(7, unx, :, :) =   (3.0*f_array2(6, unx, :, :) - &
                                            4.0*f_array2(6, unx1, :, :) + &
                                                f_array2(6, unx2, :, :)) * idxh

            fgrad_array2(8, :, lny1:uny1, :) = (f_array2(6, :, lny2:uny, :) - &
                                                f_array2(6, :, lny:uny2, :)) * idyh
            fgrad_array2(8, :, lny, :) =  (-3.0*f_array2(6, :, lny, :) + &
                                            4.0*f_array2(6, :, lny1, :) - &
                                                f_array2(6, :, lny2, :)) * idyh
            fgrad_array2(8, :, uny, :) =   (3.0*f_array2(6, :, uny, :) - &
                                            4.0*f_array2(6, :, uny1, :) + &
                                                f_array2(6, :, uny2, :)) * idyh

            fgrad_array2(13, lnx1:unx1, :, :) = (f_array2(8, lnx2:unx, :, :) - &
                                                 f_array2(8, lnx:unx2, :, :)) * idxh
            fgrad_array2(13, lnx, :, :) =  (-3.0*f_array2(8, lnx, :, :) + &
                                             4.0*f_array2(8, lnx1, :, :) - &
                                                 f_array2(8, lnx2, :, :)) * idxh
            fgrad_array2(13, unx, :, :) =   (3.0*f_array2(8, unx, :, :) - &
                                             4.0*f_array2(8, unx1, :, :) + &
                                                 f_array2(8, unx2, :, :)) * idxh

            fgrad_array2(14, :, lny1:uny1, :) = (f_array2(8, :, lny2:uny, :) - &
                                                 f_array2(8, :, lny:uny2, :)) * idyh
            fgrad_array2(14, :, lny, :) =  (-3.0*f_array2(8, :, lny, :) + &
                                             4.0*f_array2(8, :, lny1, :) - &
                                                 f_array2(8, :, lny2, :)) * idyh
            fgrad_array2(14, :, uny, :) =   (3.0*f_array2(8, :, uny, :) - &
                                             4.0*f_array2(8, :, uny1, :) + &
                                                 f_array2(8, :, uny2, :)) * idyh

        endif
    end subroutine calc_fields_gradients

    !---------------------------------------------------------------------------
    !< Interpolate the MHD fields and their gradients on one position
    !< nz is assumed to be 1.
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
        ix1 = ix + 1
        iy1 = iy + 1
        rx1 = 1.0 - rx
        ry1 = 1.0 - ry
        rt1 = 1.0 - rt
        w1 = rx1 * ry1
        w2 = rx * ry1
        w3 = rx1 * ry
        w4 = rx * ry

        fields = f_array1(:, ix, iy, 1) * w1
        fields = fields + f_array1(:, ix1, iy, 1) * w2
        fields = fields + f_array1(:, ix, iy1, 1) * w3
        fields = fields + f_array1(:, ix1, iy1, 1) * w4

        fields1 = f_array2(:, ix, iy, 1) * w1
        fields1 = fields1 + f_array2(:, ix1, iy, 1) * w2
        fields1 = fields1 + f_array2(:, ix, iy1, 1) * w3
        fields1 = fields1 + f_array2(:, ix1, iy1, 1) * w4

        !< Time interpolation
        fields = fields * rt1 + fields1 * rt

        gradf = fgrad_array1(:, ix, iy, 1) * w1
        gradf = gradf + fgrad_array1(:, ix1, iy, 1) * w2
        gradf = gradf + fgrad_array1(:, ix, iy1, 1) * w3
        gradf = gradf + fgrad_array1(:, ix1, iy1, 1) * w4

        gradf1 = fgrad_array2(:, ix, iy, 1) * w1
        gradf1 = gradf1 + fgrad_array2(:, ix1, iy, 1) * w2
        gradf1 = gradf1 + fgrad_array2(:, ix, iy1, 1) * w3
        gradf1 = gradf1 + fgrad_array2(:, ix1, iy1, 1) * w4

        !< Time interpolation
        gradf = gradf * rt1 + gradf1 * rt
    end subroutine interp_fields

    !<--------------------------------------------------------------------------
    !< Copy fields for usage in the next time interval
    !<--------------------------------------------------------------------------
    subroutine copy_fields
        implicit none
        f_array1 = f_array2
        fgrad_array1 = fgrad_array2
    end subroutine copy_fields

    !---------------------------------------------------------------------------
    !< Initialize shock x-position indices
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !<  nx, ny, nz: the dimensions of the data
    !<  ndim: number of actual dimension of the data. 1, 2 or 3
    !---------------------------------------------------------------------------
    subroutine init_shock_xpos(interp_flag, nx, ny, nz, ndim)
        implicit none
        integer, intent(in) :: interp_flag, nx, ny, nz, ndim

        if (ndim == 1) then
            allocate(shock_xpos1(1, 1))
        else if (ndim == 2) then
            allocate(shock_xpos1(-1:ny+2, 1))
        else
            allocate(shock_xpos1(-1:ny+2, -1:nz+2))
        endif
        shock_xpos1 = 0
        ! Next time step
        if (interp_flag == 1) then
            if (ndim == 1) then
                allocate(shock_xpos2(1, 1))
            else if (ndim == 2) then
                allocate(shock_xpos2(-1:ny+2, 1))
            else
                allocate(shock_xpos2(-1:ny+2, -1:nz+2))
            endif
        endif
        shock_xpos2 = 0
    end subroutine init_shock_xpos

    !---------------------------------------------------------------------------
    !< Free shock x-position indices
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine free_shock_xpos(interp_flag)
        implicit none
        integer, intent(in) :: interp_flag
        deallocate(shock_xpos1)
        if (interp_flag == 1) then
            deallocate(shock_xpos2)
        endif
    end subroutine free_shock_xpos

    !---------------------------------------------------------------------------
    !< Locate shock x-position indices. We use Vx to locate the shock here.
    !< Args:
    !<  nx, ny, nz: the dimensions of the data
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine locate_shock_xpos(interp_flag, nx, ny, nz, ndim)
        use mpi_module
        implicit none
        integer, intent(in) :: interp_flag, nx, ny, nz, ndim
        if (ndim == 1) then
            shock_xpos1(1, 1) = maxloc(abs(fgrad_array1(1, :, 1, 1)), dim=1)
        else if (ndim == 2) then
            shock_xpos1(:, 1) = maxloc(abs(fgrad_array1(1, :, :, 1)), dim=1)
        else
            shock_xpos1(:, :) = maxloc(abs(fgrad_array1(1, :, :, :)), dim=1)
        endif
        if (interp_flag == 1) then
            if (ndim == 1) then
                shock_xpos2(1, 1) = maxloc(abs(fgrad_array2(1, :, 1, 1)), dim=1)
            else if (ndim == 2) then
                shock_xpos2(:, 1) = maxloc(abs(fgrad_array2(1, :, :, 1)), dim=1)
            else
                shock_xpos2(:, :) = maxloc(abs(fgrad_array2(1, :, :, :)), dim=1)
            endif
        endif
    end subroutine locate_shock_xpos

    !---------------------------------------------------------------------------
    !< Interpolate shock location, assuming a 2D shock
    !< nz is assumed to be 1. We assume shock is propagating along x-direction
    !< Args:
    !<  ix, iy: the lower-left corner of the grid.
    !<  rx, ry: the offset to the lower-left corner of the grid where the
    !<          the position is. They are normalized to the grid sizes.
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !---------------------------------------------------------------------------
    function interp_shock_location(iy, ry, rt) result(shock_xpos)
        implicit none
        real(dp), intent(in) :: ry, rt
        integer, intent(in) :: iy
        real(dp) :: ry1, rt1, sx1, sx2, shock_xpos
        integer :: iy1
        iy1 = iy + 1
        ry1 = 1.0 - ry
        rt1 = 1.0 - rt

        !< iy starts at 0, so we need to consider the ghost cells
        sx1 = shock_xpos1(iy-1, 1) * ry1 + shock_xpos1(iy1-1, 1) * ry
        sx2 = shock_xpos2(iy-1, 1) * ry1 + shock_xpos2(iy1-1, 1) * ry

        !< Time interpolation
        shock_xpos = sx2 * rt1 + sx1 * rt
    end function interp_shock_location

end module mhd_data_parallel

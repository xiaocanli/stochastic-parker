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
           interp_fields, copy_fields
    public fields, gradf

    real(dp), dimension(8) :: fields, fields1
    real(dp), dimension(16) :: gradf, gradf1

    real(fp), allocatable, dimension(:, :, :, :) :: f_array1, f_array2 ! Current,next
    real(fp), allocatable, dimension(:, :, :, :) :: fgrad_array1, fgrad_array2
    !dir$ attributes align:64 :: f_array1
    !dir$ attributes align:64 :: f_array2
    !dir$ attributes align:64 :: fgrad_array1
    !dir$ attributes align:64 :: fgrad_array2

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
        use mpi_io_module, only: set_mpi_datatype, set_mpi_info, fileinfo, &
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
            open(unit=fh, file=filename, access='stream', status='unknown', &
                 form='unformatted', action='read')
            if (var_flag == 0) then
                if (mpi_rank == master) then
                    read(fh, pos=1) f_array1
                endif
                call MPI_BCAST(f_array1, product(sizes), MPI_REAL4, master, MPI_COMM_WORLD, ierr)
            else
                if (mpi_rank == master) then
                    read(fh, pos=1) f_array2
                endif
                call MPI_BCAST(f_array2, product(sizes), MPI_REAL4, master, MPI_COMM_WORLD, ierr)
            endif
            close(fh)
        else
            mpi_datatype = set_mpi_datatype(sizes, subsizes, starts)
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
        integer :: nx, ny, nx1, nx2, ny1, ny2
        idxh = 0.5_dp / mhd_config%dx
        idyh = 0.5_dp / mhd_config%dy
        nx = ubound(f_array1, 2)
        ny = ubound(f_array1, 3)
        nx1 = nx - 1
        nx2 = nx - 2
        ny1 = ny - 1
        ny2 = ny - 2
        if (var_flag == 0) then
            fgrad_array1(1, 2:nx1, :, :) =  (f_array1(1, 3:nx, :, :) - &
                                             f_array1(1, 1:nx2, :, :)) * idxh
            fgrad_array1(1, 1, :, :) = (-3.0*f_array1(1, 1, :, :) + &
                                         4.0*f_array1(1, 2, :, :) - &
                                             f_array1(1, 3, :, :)) * idxh
            fgrad_array1(1, nx, :, :) = (3.0*f_array1(1, nx, :, :) - &
                                         4.0*f_array1(1, nx1, :, :) + &
                                             f_array1(1, nx2, :, :)) * idxh

            fgrad_array1(2, :, 2:ny1, :) =  (f_array1(2, :, 3:ny, :) - &
                                             f_array1(2, :, 1:ny2, :)) * idyh
            fgrad_array1(2, :, 1, :) = (-3.0*f_array1(2, :, 1, :) + &
                                         4.0*f_array1(2, :, 2, :) - &
                                             f_array1(2, :, 3, :)) * idyh
            fgrad_array1(2, :, ny, :) = (3.0*f_array1(2, :, ny, :) - &
                                         4.0*f_array1(2, :, ny1, :) + &
                                             f_array1(2, :, ny2, :)) * idyh

            fgrad_array1(4, 2:nx1, :, :) =  (f_array1(5, 3:nx, :, :) - &
                                             f_array1(5, 1:nx2, :, :)) * idxh
            fgrad_array1(4, 1, :, :) = (-3.0*f_array1(5, 1, :, :) + &
                                         4.0*f_array1(5, 2, :, :) - &
                                             f_array1(5, 3, :, :)) * idxh
            fgrad_array1(4, nx, :, :) = (3.0*f_array1(5, nx, :, :) - &
                                         4.0*f_array1(5, nx1, :, :) + &
                                             f_array1(5, nx2, :, :)) * idxh

            fgrad_array1(5, :, 2:ny1, :) =  (f_array1(5, :, 3:ny, :) - &
                                             f_array1(5, :, 1:ny2, :)) * idyh
            fgrad_array1(5, :, 1, :) = (-3.0*f_array1(5, :, 1, :) + &
                                         4.0*f_array1(5, :, 2, :) - &
                                             f_array1(5, :, 3, :)) * idyh
            fgrad_array1(5, :, ny, :) = (3.0*f_array1(5, :, ny, :) - &
                                         4.0*f_array1(5, :, ny1, :) + &
                                             f_array1(5, :, ny2, :)) * idyh

            fgrad_array1(7, 2:nx1, :, :) =  (f_array1(6, 3:nx, :, :) - &
                                             f_array1(6, 1:nx2, :, :)) * idxh
            fgrad_array1(7, 1, :, :) = (-3.0*f_array1(6, 1, :, :) + &
                                         4.0*f_array1(6, 2, :, :) - &
                                             f_array1(6, 3, :, :)) * idxh
            fgrad_array1(7, nx, :, :) = (3.0*f_array1(6, nx, :, :) - &
                                         4.0*f_array1(6, nx1, :, :) + &
                                             f_array1(6, nx2, :, :)) * idxh

            fgrad_array1(8, :, 2:ny1, :) =  (f_array1(6, :, 3:ny, :) - &
                                             f_array1(6, :, 1:ny2, :)) * idyh
            fgrad_array1(8, :, 1, :) = (-3.0*f_array1(6, :, 1, :) + &
                                         4.0*f_array1(6, :, 2, :) - &
                                             f_array1(6, :, 3, :)) * idyh
            fgrad_array1(8, :, ny, :) = (3.0*f_array1(6, :, ny, :) - &
                                         4.0*f_array1(6, :, ny1, :) + &
                                             f_array1(6, :, ny2, :)) * idyh

            fgrad_array1(13, 2:nx1, :, :) =  (f_array1(8, 3:nx, :, :) - &
                                              f_array1(8, 1:nx2, :, :)) * idxh
            fgrad_array1(13, 1, :, :) = (-3.0*f_array1(8, 1, :, :) + &
                                          4.0*f_array1(8, 2, :, :) - &
                                              f_array1(8, 3, :, :)) * idxh
            fgrad_array1(13, nx, :, :) = (3.0*f_array1(8, nx, :, :) - &
                                          4.0*f_array1(8, nx1, :, :) + &
                                              f_array1(8, nx2, :, :)) * idxh

            fgrad_array1(14, :, 2:ny1, :) =  (f_array1(8, :, 3:ny, :) - &
                                              f_array1(8, :, 1:ny2, :)) * idyh
            fgrad_array1(14, :, 1, :) = (-3.0*f_array1(8, :, 1, :) + &
                                          4.0*f_array1(8, :, 2, :) - &
                                              f_array1(8, :, 3, :)) * idyh
            fgrad_array1(14, :, ny, :) = (3.0*f_array1(8, :, ny, :) - &
                                          4.0*f_array1(8, :, ny1, :) + &
                                              f_array1(8, :, ny2, :)) * idyh
        else
            fgrad_array2(1, 2:nx1, :, :) =  (f_array2(1, 3:nx, :, :) - &
                                             f_array2(1, 1:nx2, :, :)) * idxh
            fgrad_array2(1, 1, :, :) = (-3.0*f_array2(1, 1, :, :) + &
                                         4.0*f_array2(1, 2, :, :) - &
                                             f_array2(1, 3, :, :)) * idxh
            fgrad_array2(1, nx, :, :) = (3.0*f_array2(1, nx, :, :) - &
                                         4.0*f_array2(1, nx1, :, :) + &
                                             f_array2(1, nx2, :, :)) * idxh

            fgrad_array2(2, :, 2:ny1, :) =  (f_array2(2, :, 3:ny, :) - &
                                             f_array2(2, :, 1:ny2, :)) * idyh
            fgrad_array2(2, :, 1, :) = (-3.0*f_array2(2, :, 1, :) + &
                                         4.0*f_array2(2, :, 2, :) - &
                                             f_array2(2, :, 3, :)) * idyh
            fgrad_array2(2, :, ny, :) = (3.0*f_array2(2, :, ny, :) - &
                                         4.0*f_array2(2, :, ny1, :) + &
                                             f_array2(2, :, ny2, :)) * idyh

            fgrad_array2(4, 2:nx1, :, :) =  (f_array2(5, 3:nx, :, :) - &
                                             f_array2(5, 1:nx2, :, :)) * idxh
            fgrad_array2(4, 1, :, :) = (-3.0*f_array2(5, 1, :, :) + &
                                         4.0*f_array2(5, 2, :, :) - &
                                             f_array2(5, 3, :, :)) * idxh
            fgrad_array2(4, nx, :, :) = (3.0*f_array2(5, nx, :, :) - &
                                         4.0*f_array2(5, nx1, :, :) + &
                                             f_array2(5, nx2, :, :)) * idxh

            fgrad_array2(5, :, 2:ny1, :) =  (f_array2(5, :, 3:ny, :) - &
                                             f_array2(5, :, 1:ny2, :)) * idyh
            fgrad_array2(5, :, 1, :) = (-3.0*f_array2(5, :, 1, :) + &
                                         4.0*f_array2(5, :, 2, :) - &
                                             f_array2(5, :, 3, :)) * idyh
            fgrad_array2(5, :, ny, :) = (3.0*f_array2(5, :, ny, :) - &
                                         4.0*f_array2(5, :, ny1, :) + &
                                             f_array2(5, :, ny2, :)) * idyh

            fgrad_array2(7, 2:nx1, :, :) =  (f_array2(6, 3:nx, :, :) - &
                                             f_array2(6, 1:nx2, :, :)) * idxh
            fgrad_array2(7, 1, :, :) = (-3.0*f_array2(6, 1, :, :) + &
                                         4.0*f_array2(6, 2, :, :) - &
                                             f_array2(6, 3, :, :)) * idxh
            fgrad_array2(7, nx, :, :) = (3.0*f_array2(6, nx, :, :) - &
                                         4.0*f_array2(6, nx1, :, :) + &
                                             f_array2(6, nx2, :, :)) * idxh

            fgrad_array2(8, :, 2:ny1, :) =  (f_array2(6, :, 3:ny, :) - &
                                             f_array2(6, :, 1:ny2, :)) * idyh
            fgrad_array2(8, :, 1, :) = (-3.0*f_array2(6, :, 1, :) + &
                                         4.0*f_array2(6, :, 2, :) - &
                                             f_array2(6, :, 3, :)) * idyh
            fgrad_array2(8, :, ny, :) = (3.0*f_array2(6, :, ny, :) - &
                                         4.0*f_array2(6, :, ny1, :) + &
                                             f_array2(6, :, ny2, :)) * idyh

            fgrad_array2(13, 2:nx1, :, :) =  (f_array2(8, 3:nx, :, :) - &
                                              f_array2(8, 1:nx2, :, :)) * idxh
            fgrad_array2(13, 1, :, :) = (-3.0*f_array2(8, 1, :, :) + &
                                          4.0*f_array2(8, 2, :, :) - &
                                              f_array2(8, 3, :, :)) * idxh
            fgrad_array2(13, nx, :, :) = (3.0*f_array2(8, nx, :, :) - &
                                          4.0*f_array2(8, nx1, :, :) + &
                                              f_array2(8, nx2, :, :)) * idxh

            fgrad_array2(14, :, 2:ny1, :) =  (f_array2(8, :, 3:ny, :) - &
                                              f_array2(8, :, 1:ny2, :)) * idyh
            fgrad_array2(14, :, 1, :) = (-3.0*f_array2(8, :, 1, :) + &
                                          4.0*f_array2(8, :, 2, :) - &
                                              f_array2(8, :, 3, :)) * idyh
            fgrad_array2(14, :, ny, :) = (3.0*f_array2(8, :, ny, :) - &
                                          4.0*f_array2(8, :, ny1, :) + &
                                              f_array2(8, :, ny2, :)) * idyh
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

        fields = f_array1(:,  ix, iy, 1) * w1
        fields = fields + f_array1(:,  ix1, iy, 1) * w2
        fields = fields + f_array1(:,  ix, iy1, 1) * w3
        fields = fields + f_array1(:,  ix1, iy1, 1) * w4

        fields1  = f_array2(:,  ix, iy, 1) * w1
        fields1  = fields1  + f_array2(:,  ix1, iy, 1) * w2
        fields1  = fields1  + f_array2(:,  ix, iy1, 1) * w3
        fields1  = fields1  + f_array2(:,  ix1, iy1, 1) * w4

        !< Time interpolation
        fields = fields * rt1 + fields1 * rt

        gradf  = fgrad_array1(:,  ix, iy, 1) * w1
        gradf  = gradf  + fgrad_array1(:,  ix1, iy, 1) * w2
        gradf  = gradf  + fgrad_array1(:,  ix, iy1, 1) * w3
        gradf  = gradf  + fgrad_array1(:,  ix1, iy1, 1) * w4

        gradf1  = fgrad_array2(:,  ix, iy, 1) * w1
        gradf1  = gradf1  + fgrad_array2(:,  ix1, iy, 1) * w2
        gradf1  = gradf1  + fgrad_array2(:,  ix, iy1, 1) * w3
        gradf1  = gradf1  + fgrad_array2(:,  ix1, iy1, 1) * w4

        !< Time interpolation
        gradf = gradf * rt1 + gradf1 * rt
    end subroutine interp_fields

    !---------------------------------------------------------------------------
    ! Copy fields for usage in the next time interval
    !---------------------------------------------------------------------------
    subroutine copy_fields
        implicit none
        f_array1 = f_array2
        fgrad_array1 = fgrad_array2
    end subroutine copy_fields
end module mhd_data_parallel

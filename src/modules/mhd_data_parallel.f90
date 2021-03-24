!*******************************************************************************
!< Module for dealing with MHD data in parallel
!*******************************************************************************
module mhd_data_parallel
    use constants, only: fp, dp
    use simulation_setup_module, only: fconfig, ndim_field
    implicit none
    private
    save

    public init_field_data, free_field_data, read_field_data_parallel, &
           init_fields_gradients, free_fields_gradients, calc_fields_gradients, &
           calc_fields_gradients_nonuniform, interp_fields, &
           copy_fields, init_shock_xpos, free_shock_xpos, &
           locate_shock_xpos, interp_shock_location, init_magnetic_fluctuation, &
           free_magnetic_fluctuation, read_magnetic_fluctuation, &
           interp_magnetic_fluctuation, copy_magnetic_fluctuation, &
           init_grad_deltab2, free_grad_deltab2, calc_grad_deltab2, &
           calc_grad_deltab2_nonuniform, &
           init_correlation_length, free_correlation_length, &
           read_correlation_length, interp_correlation_length, &
           copy_correlation_length, init_grad_correl_length, &
           free_grad_correl_length, calc_grad_correl_length, &
           calc_grad_correl_length_nonuniform, &
           init_grid_positions, free_grid_positions, &
           set_local_grid_positions, get_ncells_large_jz
    public fields, gradf, db2, lc, grad_db, grad_lc
    public xpos_local, ypos_local, zpos_local

    real(dp) :: db2, db2_1
    real(dp) :: lc0_1, lc0_2, lc
    integer, parameter :: nfields=8
    integer, parameter :: ngrads=24
    real(dp), dimension(nfields) :: fields, fields1
    real(dp), dimension(ngrads) :: gradf, gradf1
    real(dp), dimension(3) :: grad_db, grad_db1
    real(dp), dimension(3) :: grad_lc, grad_lc1
    !dir$ attributes align:64 :: fields
    !dir$ attributes align:64 :: fields1
    !dir$ attributes align:128 :: gradf
    !dir$ attributes align:128 :: gradf1

    !< 80 variables in total
    real(fp), allocatable, dimension(:, :, :, :) :: f_array1, f_array2 ! Current,next
    real(fp), allocatable, dimension(:, :, :, :) :: fgrad_array1, fgrad_array2
    real(fp), allocatable, dimension(:, :, :) :: deltab1, deltab2 ! Current,next
    real(fp), allocatable, dimension(:, :, :) :: lc1, lc2 ! Current,next
    real(fp), allocatable, dimension(:, :, :, :) :: grad_deltab1, grad_deltab2
    real(fp), allocatable, dimension(:, :, :, :) :: grad_correl1, grad_correl2
    !dir$ attributes align:32 :: f_array1
    !dir$ attributes align:32 :: f_array2
    !dir$ attributes align:32 :: deltab1
    !dir$ attributes align:32 :: deltab2
    !dir$ attributes align:32 :: lc1
    !dir$ attributes align:32 :: lc2
    !dir$ attributes align:64 :: fgrad_array1
    !dir$ attributes align:64 :: fgrad_array2
    !dir$ attributes align:64 :: grad_deltab1
    !dir$ attributes align:64 :: grad_deltab2
    !dir$ attributes align:64 :: grad_correl1
    !dir$ attributes align:64 :: grad_correl2

    integer, allocatable, dimension(:, :) :: shock_xpos1, shock_xpos2  ! Shock x-position indices

    ! Grid positions along each direction for non-uniform grid
    real(dp), allocatable, dimension(:) :: xpos_global, ypos_global, zpos_global
    real(dp), allocatable, dimension(:) :: xpos_local, ypos_local, zpos_local
    real(dp), allocatable, dimension(:) :: idxh_grid, idyh_grid, idzh_grid

    contains
    !---------------------------------------------------------------------------
    !< Initialize MHD field data arrays
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine init_field_data(interp_flag)
        implicit none
        integer, intent(in) :: interp_flag
        integer :: nx, ny, nz

        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz

        !< vx, vy, vz, rho, bx, by, bz, btot
        if (ndim_field == 1) then
            allocate(f_array1(nfields, -1:nx+2, 1, 1))
        else if (ndim_field == 2) then
            allocate(f_array1(nfields, -1:nx+2, -1:ny+2, 1))
        else
            allocate(f_array1(nfields, -1:nx+2, -1:ny+2, -1:nz+2))
        endif
        f_array1 = 0.0
        ! Next time step
        if (interp_flag == 1) then
            if (ndim_field == 1) then
                allocate(f_array2(nfields, -1:nx+2, 1, 1))
            else if (ndim_field == 2) then
                allocate(f_array2(nfields, -1:nx+2, -1:ny+2, 1))
            else
                allocate(f_array2(nfields, -1:nx+2, -1:ny+2, -1:nz+2))
            endif
            f_array2 = 0.0
        endif
    end subroutine init_field_data

    !---------------------------------------------------------------------------
    !< Initialize the gradients of the MHD data arrays.
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine init_fields_gradients(interp_flag)
        implicit none
        integer, intent(in) :: interp_flag
        integer :: nx, ny, nz

        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz

        !< dvx_dx, dvx_dy, dvx_dz, dvy_dx, dvy_dy, dvy_dz,
        !< dvz_dx, dvz_dy, dvz_dz, drho_dx, drho_d, drho_dz,
        !< dbx_dx, dbx_dy, dbx_dz, dby_dx, dby_dy, dby_dz,
        !< dbz_dx, dbz_dy, dbz_dz, dbtot_dx, dbtot_dy, dbtot_dz
        if (ndim_field == 1) then
            allocate(fgrad_array1(ngrads, -1:nx+2, 1, 1))
        else if (ndim_field == 2) then
            allocate(fgrad_array1(ngrads, -1:nx+2, -1:ny+2, 1))
        else
            allocate(fgrad_array1(ngrads, -1:nx+2, -1:ny+2, -1:nz+2))
        endif
        fgrad_array1 = 0.0
        ! Next time step
        if (interp_flag == 1) then
            if (ndim_field == 1) then
                allocate(fgrad_array2(ngrads, -1:nx+2, 1, 1))
            else if (ndim_field == 2) then
                allocate(fgrad_array2(ngrads, -1:nx+2, -1:ny+2, 1))
            else
                allocate(fgrad_array2(ngrads, -1:nx+2, -1:ny+2, -1:nz+2))
            endif
            fgrad_array2 = 0.0
        endif
    end subroutine init_fields_gradients

    !---------------------------------------------------------------------------
    !< Initialize the magnetic fluctuation, defined as deltaB**2 / B0**2
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine init_magnetic_fluctuation(interp_flag)
        implicit none
        integer, intent(in) :: interp_flag
        integer :: nx, ny, nz

        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz

        if (ndim_field == 1) then
            allocate(deltab1(-1:nx+2, 1, 1))
        else if (ndim_field == 2) then
            allocate(deltab1(-1:nx+2, -1:ny+2, 1))
        else
            allocate(deltab1(-1:nx+2, -1:ny+2, -1:nz+2))
        endif
        deltab1 = 1.0
        ! Next time step
        if (interp_flag == 1) then
            if (ndim_field == 1) then
                allocate(deltab2(-1:nx+2, 1, 1))
            else if (ndim_field == 2) then
                allocate(deltab2(-1:nx+2, -1:ny+2, 1))
            else
                allocate(deltab2(-1:nx+2, -1:ny+2, -1:nz+2))
            endif
            deltab2 = 1.0
        endif
    end subroutine init_magnetic_fluctuation

    !---------------------------------------------------------------------------
    !< Initialize the gradients of the magnetic fluctuation.
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine init_grad_deltab2(interp_flag)
        implicit none
        integer, intent(in) :: interp_flag
        integer :: nx, ny, nz

        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz

        !< 3 components
        if (ndim_field == 1) then
            allocate(grad_deltab1(3, -1:nx+2, 1, 1))
        else if (ndim_field == 2) then
            allocate(grad_deltab1(3, -1:nx+2, -1:ny+2, 1))
        else
            allocate(grad_deltab1(3, -1:nx+2, -1:ny+2, -1:nz+2))
        endif
        grad_deltab1 = 0.0
        ! Next time step
        if (interp_flag == 1) then
            if (ndim_field == 1) then
                allocate(grad_deltab2(3, -1:nx+2, 1, 1))
            else if (ndim_field == 2) then
                allocate(grad_deltab2(3, -1:nx+2, -1:ny+2, 1))
            else
                allocate(grad_deltab2(3, -1:nx+2, -1:ny+2, -1:nz+2))
            endif
            grad_deltab2 = 0.0
        endif
    end subroutine init_grad_deltab2

    !---------------------------------------------------------------------------
    !< Initialize the turbulence correlation length
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine init_correlation_length(interp_flag)
        implicit none
        integer, intent(in) :: interp_flag
        integer :: nx, ny, nz

        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz

        if (ndim_field == 1) then
            allocate(lc1(-1:nx+2, 1, 1))
        else if (ndim_field == 2) then
            allocate(lc1(-1:nx+2, -1:ny+2, 1))
        else
            allocate(lc1(-1:nx+2, -1:ny+2, -1:nz+2))
        endif
        lc1 = 1.0
        ! Next time step
        if (interp_flag == 1) then
            if (ndim_field == 1) then
                allocate(lc2(-1:nx+2, 1, 1))
            else if (ndim_field == 2) then
                allocate(lc2(-1:nx+2, -1:ny+2, 1))
            else
                allocate(lc2(-1:nx+2, -1:ny+2, -1:nz+2))
            endif
            lc2 = 1.0
        endif
    end subroutine init_correlation_length

    !---------------------------------------------------------------------------
    !< Initialize the gradients of turbulence correlation length.
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine init_grad_correl_length(interp_flag)
        implicit none
        integer, intent(in) :: interp_flag
        integer :: nx, ny, nz

        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz

        !< 3 components
        if (ndim_field == 1) then
            allocate(grad_correl1(3, -1:nx+2, 1, 1))
        else if (ndim_field == 2) then
            allocate(grad_correl1(3, -1:nx+2, -1:ny+2, 1))
        else
            allocate(grad_correl1(3, -1:nx+2, -1:ny+2, -1:nz+2))
        endif
        grad_correl1 = 0.0
        ! Next time step
        if (interp_flag == 1) then
            if (ndim_field == 1) then
                allocate(grad_correl2(3, -1:nx+2, 1, 1))
            else if (ndim_field == 2) then
                allocate(grad_correl2(3, -1:nx+2, -1:ny+2, 1))
            else
                allocate(grad_correl2(3, -1:nx+2, -1:ny+2, -1:nz+2))
            endif
            grad_correl2 = 0.0
        endif
    end subroutine init_grad_correl_length

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
    !< Free the data array for magnetic fluctuation
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine free_magnetic_fluctuation(interp_flag)
        implicit none
        integer, intent(in) :: interp_flag
        deallocate(deltab1)
        if (interp_flag == 1) then
            deallocate(deltab2)
        endif
    end subroutine free_magnetic_fluctuation

    !---------------------------------------------------------------------------
    !< Free the gradients of magnetic fluctuation.
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine free_grad_deltab2(interp_flag)
        implicit none
        integer, intent(in) :: interp_flag
        deallocate(grad_deltab1)
        if (interp_flag == 1) then
            deallocate(grad_deltab2)
        endif
    end subroutine free_grad_deltab2

    !---------------------------------------------------------------------------
    !< Free the data array for turbulence correlation length
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine free_correlation_length(interp_flag)
        implicit none
        integer, intent(in) :: interp_flag
        deallocate(lc1)
        if (interp_flag == 1) then
            deallocate(lc2)
        endif
    end subroutine free_correlation_length

    !---------------------------------------------------------------------------
    !< Free the gradients of turbulence correlation length.
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine free_grad_correl_length(interp_flag)
        implicit none
        integer, intent(in) :: interp_flag
        deallocate(grad_correl1)
        if (interp_flag == 1) then
            deallocate(grad_correl2)
        endif
    end subroutine free_grad_correl_length

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
        use mpi_module
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: var_flag
        integer :: mpi_datatype, fh
        integer, dimension(4) :: sizes, subsizes, starts
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
        sizes(1) = 8 ! vx, vy, vz, rho, bx, by, bz, btot
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
    !< Read magnetic fluctuation data
    !< Args:
    !<  filename: file name to get the data
    !<  var_flag: indicating which set of variables to save the data. 0 for
    !<            detalb1 and other numbers for detalb2.
    !---------------------------------------------------------------------------
    subroutine read_magnetic_fluctuation(filename, var_flag)
        use mpi_io_module, only: set_mpi_datatype_real, set_mpi_info, fileinfo, &
            open_data_mpi_io, read_data_mpi_io
        use mpi_module
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: var_flag
        integer :: mpi_datatype, fh
        integer, dimension(3) :: sizes, subsizes, starts
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
        sizes(1) = fconfig%nxg
        sizes(2) = fconfig%nyg
        sizes(3) = fconfig%nzg
        subsizes(1) = fconfig%nxf
        subsizes(2) = fconfig%nyf
        subsizes(3) = fconfig%nzf
        starts(1) = fconfig%ix_min - 1
        starts(2) = fconfig%iy_min - 1
        starts(3) = fconfig%iz_min - 1
        fh = 11
        if (all(sizes == subsizes)) then
            if (mpi_rank == master) then
                write(*, "(A)") " Each MPI rank reads the whole data set."
            endif
            open(unit=fh, file=filename, access='stream', status='unknown', &
                 form='unformatted', action='read')
            if (var_flag == 0) then
                if (mpi_rank == master) then
                    read(fh, pos=1) deltab1
                endif
                call MPI_BCAST(deltab1, product(sizes), MPI_REAL4, master, MPI_COMM_WORLD, ierr)
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            else
                if (mpi_rank == master) then
                    read(fh, pos=1) deltab2
                endif
                call MPI_BCAST(deltab2, product(sizes), MPI_REAL4, master, MPI_COMM_WORLD, ierr)
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
                call read_data_mpi_io(fh, mpi_datatype, subsizes, disp, offset, deltab1)
            else
                call read_data_mpi_io(fh, mpi_datatype, subsizes, disp, offset, deltab2)
            endif
            call MPI_FILE_CLOSE(fh, ierror)
        endif
        if (mpi_rank == master) then
            write(*, "(A)") "Finished reading magnetic fluctuation data."
        endif
    end subroutine read_magnetic_fluctuation

    !---------------------------------------------------------------------------
    !< Read turbulence correlation length
    !< Args:
    !<  filename: file name to get the data
    !<  var_flag: indicating which set of variables to save the data. 0 for
    !<            lc1 and other numbers for lc2.
    !---------------------------------------------------------------------------
    subroutine read_correlation_length(filename, var_flag)
        use mpi_io_module, only: set_mpi_datatype_real, set_mpi_info, fileinfo, &
            open_data_mpi_io, read_data_mpi_io
        use mpi_module
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: var_flag
        integer :: mpi_datatype, fh
        integer, dimension(3) :: sizes, subsizes, starts
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
        sizes(1) = fconfig%nxg
        sizes(2) = fconfig%nyg
        sizes(3) = fconfig%nzg
        subsizes(1) = fconfig%nxf
        subsizes(2) = fconfig%nyf
        subsizes(3) = fconfig%nzf
        starts(1) = fconfig%ix_min - 1
        starts(2) = fconfig%iy_min - 1
        starts(3) = fconfig%iz_min - 1
        fh = 11
        if (all(sizes == subsizes)) then
            if (mpi_rank == master) then
                write(*, "(A)") " Each MPI rank reads the whole data set."
            endif
            open(unit=fh, file=filename, access='stream', status='unknown', &
                 form='unformatted', action='read')
            if (var_flag == 0) then
                if (mpi_rank == master) then
                    read(fh, pos=1) lc1
                endif
                call MPI_BCAST(lc1, product(sizes), MPI_REAL4, master, MPI_COMM_WORLD, ierr)
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            else
                if (mpi_rank == master) then
                    read(fh, pos=1) lc2
                endif
                call MPI_BCAST(lc2, product(sizes), MPI_REAL4, master, MPI_COMM_WORLD, ierr)
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
                call read_data_mpi_io(fh, mpi_datatype, subsizes, disp, offset, lc1)
            else
                call read_data_mpi_io(fh, mpi_datatype, subsizes, disp, offset, lc2)
            endif
            call MPI_FILE_CLOSE(fh, ierror)
        endif
        if (mpi_rank == master) then
            write(*, "(A)") "Finished reading turbulence correlation length."
        endif
    end subroutine read_correlation_length

    !---------------------------------------------------------------------------
    !< Calculate the gradients of the MHD data arrays.
    !< Args:
    !<  var_flag: indicating which set of variables.
    !<            0 for fgrad_array1 and other numbers for fgrad_array2.
    !---------------------------------------------------------------------------
    subroutine calc_fields_gradients(var_flag)
        use mhd_config_module, only: mhd_config
        use mpi_module
        implicit none
        integer, intent(in) :: var_flag
        real(dp) :: idxh, idyh, idzh
        integer :: unx, uny, unz, lnx, lny, lnz
        integer :: unx1, unx2, uny1, uny2, unz1, unz2
        integer :: lnx1, lnx2, lny1, lny2, lnz1, lnz2
        idxh = 0.5_dp / mhd_config%dx
        idyh = 0.5_dp / mhd_config%dy
        idzh = 0.5_dp / mhd_config%dz
        unx = ubound(f_array1, 2)
        uny = ubound(f_array1, 3)
        unz = ubound(f_array1, 4)
        lnx = lbound(f_array1, 2)
        lny = lbound(f_array1, 3)
        lnz = lbound(f_array1, 4)
        unx1 = unx - 1
        unx2 = unx - 2
        uny1 = uny - 1
        uny2 = uny - 2
        unz1 = unz - 1
        unz2 = unz - 2
        lnx1 = lnx + 1
        lnx2 = lnx + 2
        lny1 = lny + 1
        lny2 = lny + 2
        lnz1 = lnz + 1
        lnz2 = lnz + 2
        if (var_flag == 0) then
            ! d/dx
            fgrad_array1(1::3, lnx1:unx1, :, :) = (f_array1(:, lnx2:unx, :, :) - &
                                                   f_array1(:, lnx:unx2, :, :)) * idxh
            fgrad_array1(1::3, lnx, :, :) =  (-3.0*f_array1(:, lnx, :, :) + &
                                               4.0*f_array1(:, lnx1, :, :) - &
                                                   f_array1(:, lnx2, :, :)) * idxh
            fgrad_array1(1::3, unx, :, :) =   (3.0*f_array1(:, unx, :, :) - &
                                               4.0*f_array1(:, unx1, :, :) + &
                                                   f_array1(:, unx2, :, :)) * idxh

            ! d/dy
            if (uny > lny) then
                fgrad_array1(2::3, :, lny1:uny1, :) = (f_array1(:, :, lny2:uny, :) - &
                                                       f_array1(:, :, lny:uny2, :)) * idyh
                fgrad_array1(2::3, :, lny, :) =  (-3.0*f_array1(:, :, lny, :) + &
                                                   4.0*f_array1(:, :, lny1, :) - &
                                                       f_array1(:, :, lny2, :)) * idyh
                fgrad_array1(2::3, :, uny, :) =   (3.0*f_array1(:, :, uny, :) - &
                                                   4.0*f_array1(:, :, uny1, :) + &
                                                       f_array1(:, :, uny2, :)) * idyh
            endif

            ! d/dz
            if (unz > lnz) then
                fgrad_array1(3::3, :, :, lnz1:unz1) = (f_array1(:, :, :, lnz2:unz) - &
                                                       f_array1(:, :, :, lnz:unz2)) * idzh
                fgrad_array1(3::3, :, :, lnz) =  (-3.0*f_array1(:, :, :, lnz) + &
                                                   4.0*f_array1(:, :, :, lnz1) - &
                                                       f_array1(:, :, :, lnz2)) * idzh
                fgrad_array1(3::3, :, :, unz) =   (3.0*f_array1(:, :, :, unz) - &
                                                   4.0*f_array1(:, :, :, unz1) + &
                                                       f_array1(:, :, :, unz2)) * idzh
            endif
        else
            ! d/dx
            fgrad_array2(1::3, lnx1:unx1, :, :) = (f_array2(:, lnx2:unx, :, :) - &
                                                   f_array2(:, lnx:unx2, :, :)) * idxh
            fgrad_array2(1::3, lnx, :, :) =  (-3.0*f_array2(:, lnx, :, :) + &
                                               4.0*f_array2(:, lnx1, :, :) - &
                                                   f_array2(:, lnx2, :, :)) * idxh
            fgrad_array2(1::3, unx, :, :) =   (3.0*f_array2(:, unx, :, :) - &
                                               4.0*f_array2(:, unx1, :, :) + &
                                                   f_array2(:, unx2, :, :)) * idxh

            ! d/dy
            if (uny > lny) then
                fgrad_array2(2::3, :, lny1:uny1, :) = (f_array2(:, :, lny2:uny, :) - &
                                                       f_array2(:, :, lny:uny2, :)) * idyh
                fgrad_array2(2::3, :, lny, :) =  (-3.0*f_array2(:, :, lny, :) + &
                                                   4.0*f_array2(:, :, lny1, :) - &
                                                       f_array2(:, :, lny2, :)) * idyh
                fgrad_array2(2::3, :, uny, :) =   (3.0*f_array2(:, :, uny, :) - &
                                                   4.0*f_array2(:, :, uny1, :) + &
                                                       f_array2(:, :, uny2, :)) * idyh
            endif

            ! d/dz
            if (unz > lnz) then
                fgrad_array2(3::3, :, :, lnz1:unz1) = (f_array2(:, :, :, lnz2:unz) - &
                                                       f_array2(:, :, :, lnz:unz2)) * idzh
                fgrad_array2(3::3, :, :, lnz) =  (-3.0*f_array2(:, :, :, lnz) + &
                                                   4.0*f_array2(:, :, :, lnz1) - &
                                                       f_array2(:, :, :, lnz2)) * idzh
                fgrad_array2(3::3, :, :, unz) =   (3.0*f_array2(:, :, :, unz) - &
                                                   4.0*f_array2(:, :, :, unz1) + &
                                                       f_array2(:, :, :, unz2)) * idzh
            endif
        endif
        if (mpi_rank == master) then
            write(*, "(A)") "Finished calculating fields gradients."
        endif
    end subroutine calc_fields_gradients

    !---------------------------------------------------------------------------
    !< Calculate the gradients of the MHD data arrays for non-uniform grid.
    !< Args:
    !<  var_flag: indicating which set of variables.
    !<            0 for fgrad_array1 and other numbers for fgrad_array2.
    !---------------------------------------------------------------------------
    subroutine calc_fields_gradients_nonuniform(var_flag)
        use mhd_config_module, only: mhd_config
        use mpi_module
        implicit none
        integer, intent(in) :: var_flag
        integer :: unx, uny, unz, lnx, lny, lnz
        integer :: unx1, unx2, uny1, uny2, unz1, unz2
        integer :: lnx1, lnx2, lny1, lny2, lnz1, lnz2
        integer :: iv, ix, iy, iz
        unx = ubound(f_array1, 2)
        uny = ubound(f_array1, 3)
        unz = ubound(f_array1, 4)
        lnx = lbound(f_array1, 2)
        lny = lbound(f_array1, 3)
        lnz = lbound(f_array1, 4)
        unx1 = unx - 1
        unx2 = unx - 2
        uny1 = uny - 1
        uny2 = uny - 2
        unz1 = unz - 1
        unz2 = unz - 2
        lnx1 = lnx + 1
        lnx2 = lnx + 2
        lny1 = lny + 1
        lny2 = lny + 2
        lnz1 = lnz + 1
        lnz2 = lnz + 2
        if (var_flag == 0) then
            ! d/dx
            do iz = lnz, unz
            do iy = lny, uny
            do iv = 1, nfields
                fgrad_array1(3*iv-2, lnx1:unx1, iy, iz) = &
                    (f_array1(iv, lnx2:unx, iy, iz) - &
                     f_array1(iv, lnx:unx2, iy, iz)) * idxh_grid(lnx1:unx1)
            enddo
            enddo
            enddo
            fgrad_array1(1::3, lnx, :, :) =  (-3.0*f_array1(:, lnx, :, :) + &
                                               4.0*f_array1(:, lnx1, :, :) - &
                                                   f_array1(:, lnx2, :, :)) * idxh_grid(lnx)
            fgrad_array1(1::3, unx, :, :) =   (3.0*f_array1(:, unx, :, :) - &
                                               4.0*f_array1(:, unx1, :, :) + &
                                                   f_array1(:, unx2, :, :)) * idxh_grid(unx)

            ! d/dy
            if (uny > lny) then
                do iz = lnz, unz
                do ix = lnx, unx
                do iv = 1, nfields
                    fgrad_array1(3*iv-1, ix, lny1:uny1, iz) = &
                        (f_array1(iv, ix, lny2:uny, iz) - &
                         f_array1(iv, ix, lny:uny2, iz)) * idyh_grid(lny1:uny1)
                enddo
                enddo
                enddo
                fgrad_array1(2::3, :, lny, :) =  (-3.0*f_array1(:, :, lny, :) + &
                                                   4.0*f_array1(:, :, lny1, :) - &
                                                       f_array1(:, :, lny2, :)) * idyh_grid(lny)
                fgrad_array1(2::3, :, uny, :) =   (3.0*f_array1(:, :, uny, :) - &
                                                   4.0*f_array1(:, :, uny1, :) + &
                                                       f_array1(:, :, uny2, :)) * idyh_grid(uny)
            endif

            ! d/dz
            if (unz > lnz) then
                do iy = lny, uny
                do ix = lnx, unx
                do iv = 1, nfields
                    fgrad_array1(3*iv, ix, iy, lnz1:unz1) = &
                        (f_array1(iv, ix, iy, lnz2:unz) - &
                         f_array1(iv, ix, iy, lnz:unz2)) * idzh_grid(lnz1:unz1)
                enddo
                enddo
                enddo
                fgrad_array1(3::3, :, :, lnz) =  (-3.0*f_array1(:, :, :, lnz) + &
                                                   4.0*f_array1(:, :, :, lnz1) - &
                                                       f_array1(:, :, :, lnz2)) * idzh_grid(lnz)
                fgrad_array1(3::3, :, :, unz) =   (3.0*f_array1(:, :, :, unz) - &
                                                   4.0*f_array1(:, :, :, unz1) + &
                                                       f_array1(:, :, :, unz2)) * idzh_grid(unz)
            endif
        else
            ! d/dx
            do iz = lnz, unz
            do iy = lny, uny
            do iv = 1, nfields
                fgrad_array2(3*iv-2, lnx1:unx1, iy, iz) = &
                    (f_array2(iv, lnx2:unx, iy, iz) - &
                     f_array2(iv, lnx:unx2, iy, iz)) * idxh_grid(lnx1:unx1)
            enddo
            enddo
            enddo
            fgrad_array2(1::3, lnx, :, :) =  (-3.0*f_array2(:, lnx, :, :) + &
                                               4.0*f_array2(:, lnx1, :, :) - &
                                                   f_array2(:, lnx2, :, :)) * idxh_grid(lnx)
            fgrad_array2(1::3, unx, :, :) =   (3.0*f_array2(:, unx, :, :) - &
                                               4.0*f_array2(:, unx1, :, :) + &
                                                   f_array2(:, unx2, :, :)) * idxh_grid(unx)

            ! d/dy
            if (uny > lny) then
                do iz = lnz, unz
                do ix = lnx, unx
                do iv = 1, nfields
                    fgrad_array2(3*iv-1, ix, lny1:uny1, iz) = &
                        (f_array2(iv, ix, lny2:uny, iz) - &
                         f_array2(iv, ix, lny:uny2, iz)) * idyh_grid(lny1:uny1)
                enddo
                enddo
                enddo
                fgrad_array2(2::3, :, lny, :) =  (-3.0*f_array2(:, :, lny, :) + &
                                                   4.0*f_array2(:, :, lny1, :) - &
                                                       f_array2(:, :, lny2, :)) * idyh_grid(lny)
                fgrad_array2(2::3, :, uny, :) =   (3.0*f_array2(:, :, uny, :) - &
                                                   4.0*f_array2(:, :, uny1, :) + &
                                                       f_array2(:, :, uny2, :)) * idyh_grid(uny)
            endif

            ! d/dz
            if (unz > lnz) then
                do iy = lny, uny
                do ix = lnx, unx
                do iv = 1, nfields
                    fgrad_array2(3*iv, ix, iy, lnz1:unz1) = &
                        (f_array2(iv, ix, iy, lnz2:unz) - &
                         f_array2(iv, ix, iy, lnz:unz2)) * idzh_grid(lnz1:unz1)
                enddo
                enddo
                enddo
                fgrad_array2(3::3, :, :, lnz) =  (-3.0*f_array2(:, :, :, lnz) + &
                                                   4.0*f_array2(:, :, :, lnz1) - &
                                                       f_array2(:, :, :, lnz2)) * idzh_grid(lnz)
                fgrad_array2(3::3, :, :, unz) =   (3.0*f_array2(:, :, :, unz) - &
                                                   4.0*f_array2(:, :, :, unz1) + &
                                                       f_array2(:, :, :, unz2)) * idzh_grid(unz)
            endif
        endif
        if (mpi_rank == master) then
            write(*, "(A)") "Finished calculating fields gradients."
        endif
    end subroutine calc_fields_gradients_nonuniform

    !---------------------------------------------------------------------------
    !< Calculate the gradients of magnetic fluctuation.
    !< Args:
    !<  var_flag: indicating which set of variables.
    !<            0 for grad_deltab1 and other numbers for grad_deltab2.
    !---------------------------------------------------------------------------
    subroutine calc_grad_deltab2(var_flag)
        use mhd_config_module, only: mhd_config
        use mpi_module
        implicit none
        integer, intent(in) :: var_flag
        real(dp) :: idxh, idyh, idzh
        integer :: unx, uny, unz, lnx, lny, lnz
        integer :: unx1, unx2, uny1, uny2, unz1, unz2
        integer :: lnx1, lnx2, lny1, lny2, lnz1, lnz2
        idxh = 0.5_dp / mhd_config%dx
        idyh = 0.5_dp / mhd_config%dy
        idzh = 0.5_dp / mhd_config%dz
        unx = ubound(deltab1, 1)
        uny = ubound(deltab1, 2)
        unz = ubound(deltab1, 3)
        lnx = lbound(deltab1, 1)
        lny = lbound(deltab1, 2)
        lnz = lbound(deltab1, 3)
        unx1 = unx - 1
        unx2 = unx - 2
        uny1 = uny - 1
        uny2 = uny - 2
        unz1 = unz - 1
        unz2 = unz - 2
        lnx1 = lnx + 1
        lnx2 = lnx + 2
        lny1 = lny + 1
        lny2 = lny + 2
        lnz1 = lnz + 1
        lnz2 = lnz + 2
        if (var_flag == 0) then
            ! d/dx
            grad_deltab1(1, lnx1:unx1, :, :) = (deltab1(lnx2:unx, :, :) - &
                                                deltab1(lnx:unx2, :, :)) * idxh
            grad_deltab1(1, lnx, :, :) =  (-3.0*deltab1(lnx, :, :) + &
                                            4.0*deltab1(lnx1, :, :) - &
                                                deltab1(lnx2, :, :)) * idxh
            grad_deltab1(1, unx, :, :) =   (3.0*deltab1(unx, :, :) - &
                                            4.0*deltab1(unx1, :, :) + &
                                                deltab1(unx2, :, :)) * idxh

            ! d/dy
            if (uny > lny) then
                grad_deltab1(2, :, lny1:uny1, :) = (deltab1(:, lny2:uny, :) - &
                                                    deltab1(:, lny:uny2, :)) * idyh
                grad_deltab1(2, :, lny, :) =  (-3.0*deltab1(:, lny, :) + &
                                                4.0*deltab1(:, lny1, :) - &
                                                    deltab1(:, lny2, :)) * idyh
                grad_deltab1(2, :, uny, :) =   (3.0*deltab1(:, uny, :) - &
                                                4.0*deltab1(:, uny1, :) + &
                                                    deltab1(:, uny2, :)) * idyh
            endif

            ! d/dz
            if (unz > lnz) then
                grad_deltab1(3, :, :, lnz1:unz1) = (deltab1(:, :, lnz2:unz) - &
                                                    deltab1(:, :, lnz:unz2)) * idzh
                grad_deltab1(3, :, :, lnz) =  (-3.0*deltab1(:, :, lnz) + &
                                                4.0*deltab1(:, :, lnz1) - &
                                                    deltab1(:, :, lnz2)) * idzh
                grad_deltab1(3, :, :, unz) =   (3.0*deltab1(:, :, unz) - &
                                                4.0*deltab1(:, :, unz1) + &
                                                    deltab1(:, :, unz2)) * idzh
            endif
        else
            ! d/dx
            grad_deltab2(1, lnx1:unx1, :, :) = (deltab2(lnx2:unx, :, :) - &
                                                deltab2(lnx:unx2, :, :)) * idxh
            grad_deltab2(1, lnx, :, :) =  (-3.0*deltab2(lnx, :, :) + &
                                            4.0*deltab2(lnx1, :, :) - &
                                                deltab2(lnx2, :, :)) * idxh
            grad_deltab2(1, unx, :, :) =   (3.0*deltab2(unx, :, :) - &
                                            4.0*deltab2(unx1, :, :) + &
                                                deltab2(unx2, :, :)) * idxh

            ! d/dy
            if (uny > lny) then
                grad_deltab2(2, :, lny1:uny1, :) = (deltab2(:, lny2:uny, :) - &
                                                    deltab2(:, lny:uny2, :)) * idyh
                grad_deltab2(2, :, lny, :) =  (-3.0*deltab2(:, lny, :) + &
                                                4.0*deltab2(:, lny1, :) - &
                                                    deltab2(:, lny2, :)) * idyh
                grad_deltab2(2, :, uny, :) =   (3.0*deltab2(:, uny, :) - &
                                                4.0*deltab2(:, uny1, :) + &
                                                    deltab2(:, uny2, :)) * idyh
            endif

            ! d/dz
            if (unz > lnz) then
                grad_deltab2(3, :, :, lnz1:unz1) = (deltab2(:, :, lnz2:unz) - &
                                                    deltab2(:, :, lnz:unz2)) * idzh
                grad_deltab2(3, :, :, lnz) =  (-3.0*deltab2(:, :, lnz) + &
                                                4.0*deltab2(:, :, lnz1) - &
                                                    deltab2(:, :, lnz2)) * idzh
                grad_deltab2(3, :, :, unz) =   (3.0*deltab2(:, :, unz) - &
                                                4.0*deltab2(:, :, unz1) + &
                                                    deltab2(:, :, unz2)) * idzh
            endif
        endif
        if (mpi_rank == master) then
            write(*, "(A)") "Finished calculating gradients of magnetic fluctuation."
        endif
    end subroutine calc_grad_deltab2

    !---------------------------------------------------------------------------
    !< Calculate the gradients of magnetic fluctuation for nonuniform grid
    !< Args:
    !<  var_flag: indicating which set of variables.
    !<            0 for grad_deltab1 and other numbers for grad_deltab2.
    !---------------------------------------------------------------------------
    subroutine calc_grad_deltab2_nonuniform(var_flag)
        use mhd_config_module, only: mhd_config
        use mpi_module
        implicit none
        integer, intent(in) :: var_flag
        integer :: unx, uny, unz, lnx, lny, lnz
        integer :: unx1, unx2, uny1, uny2, unz1, unz2
        integer :: lnx1, lnx2, lny1, lny2, lnz1, lnz2
        integer :: ix, iy, iz
        unx = ubound(deltab1, 1)
        uny = ubound(deltab1, 2)
        unz = ubound(deltab1, 3)
        lnx = lbound(deltab1, 1)
        lny = lbound(deltab1, 2)
        lnz = lbound(deltab1, 3)
        unx1 = unx - 1
        unx2 = unx - 2
        uny1 = uny - 1
        uny2 = uny - 2
        unz1 = unz - 1
        unz2 = unz - 2
        lnx1 = lnx + 1
        lnx2 = lnx + 2
        lny1 = lny + 1
        lny2 = lny + 2
        lnz1 = lnz + 1
        lnz2 = lnz + 2
        if (var_flag == 0) then
            ! d/dx
            do iz = lnz, unz
            do iy = lny, uny
                grad_deltab1(1, lnx1:unx1, iy, iz) = &
                    (deltab1(lnx2:unx, iy, iz) - &
                     deltab1(lnx:unx2, iy, iz)) * idxh_grid(lnx1:unx1)
            enddo
            enddo
            grad_deltab1(1, lnx, :, :) =  (-3.0*deltab1(lnx, :, :) + &
                                            4.0*deltab1(lnx1, :, :) - &
                                                deltab1(lnx2, :, :)) * idxh_grid(lnx)
            grad_deltab1(1, unx, :, :) =   (3.0*deltab1(unx, :, :) - &
                                            4.0*deltab1(unx1, :, :) + &
                                                deltab1(unx2, :, :)) * idxh_grid(unx)

            ! d/dy
            if (uny > lny) then
                do iz = lnz, unz
                do ix = lnx, unx
                    grad_deltab1(2, ix, lny1:uny1, iz) = &
                        (deltab1(ix, lny2:uny, iz) - &
                         deltab1(ix, lny:uny2, iz)) * idyh_grid(lny1:uny1)
                enddo
                enddo
                grad_deltab1(2, :, lny, :) =  (-3.0*deltab1(:, lny, :) + &
                                                4.0*deltab1(:, lny1, :) - &
                                                    deltab1(:, lny2, :)) * idyh_grid(lny)
                grad_deltab1(2, :, uny, :) =   (3.0*deltab1(:, uny, :) - &
                                                4.0*deltab1(:, uny1, :) + &
                                                    deltab1(:, uny2, :)) * idyh_grid(uny)
            endif

            ! d/dz
            if (unz > lnz) then
                do iy = lny, uny
                do ix = lnx, unx
                    grad_deltab1(3, ix, iy, lnz1:unz1) = &
                        (deltab1(ix, iy, lnz2:unz) - &
                         deltab1(ix, iy, lnz:unz2)) * idzh_grid(lnz1:unz1)
                enddo
                enddo
                grad_deltab1(3, :, :, lnz) =  (-3.0*deltab1(:, :, lnz) + &
                                                4.0*deltab1(:, :, lnz1) - &
                                                    deltab1(:, :, lnz2)) * idzh_grid(lnz)
                grad_deltab1(3, :, :, unz) =   (3.0*deltab1(:, :, unz) - &
                                                4.0*deltab1(:, :, unz1) + &
                                                    deltab1(:, :, unz2)) * idzh_grid(unz)
            endif
        else
            ! d/dx
            do iz = lnz, unz
            do iy = lny, uny
                grad_deltab2(1, lnx1:unx1, iy, iz) = &
                    (deltab2(lnx2:unx, iy, iz) - &
                     deltab2(lnx:unx2, iy, iz)) * idxh_grid(lnx1:unx1)
            enddo
            enddo
            grad_deltab2(1, lnx, :, :) =  (-3.0*deltab2(lnx, :, :) + &
                                            4.0*deltab2(lnx1, :, :) - &
                                                deltab2(lnx2, :, :)) * idxh_grid(lnx)
            grad_deltab2(1, unx, :, :) =   (3.0*deltab2(unx, :, :) - &
                                            4.0*deltab2(unx1, :, :) + &
                                                deltab2(unx2, :, :)) * idxh_grid(unx)

            ! d/dy
            if (uny > lny) then
                do iz = lnz, unz
                do ix = lnx, unx
                    grad_deltab2(2, ix, lny1:uny1, iz) = &
                        (deltab2(ix, lny2:uny, iz) - &
                         deltab2(ix, lny:uny2, iz)) * idyh_grid(lny1:uny1)
                enddo
                enddo
                grad_deltab2(2, :, lny, :) =  (-3.0*deltab2(:, lny, :) + &
                                                4.0*deltab2(:, lny1, :) - &
                                                    deltab2(:, lny2, :)) * idyh_grid(lny)
                grad_deltab2(2, :, uny, :) =   (3.0*deltab2(:, uny, :) - &
                                                4.0*deltab2(:, uny1, :) + &
                                                    deltab2(:, uny2, :)) * idyh_grid(uny)
            endif

            ! d/dz
            if (unz > lnz) then
                do iy = lny, uny
                do ix = lnx, unx
                    grad_deltab2(3, ix, iy, lnz1:unz1) = &
                        (deltab2(ix, iy, lnz2:unz) - &
                         deltab2(ix, iy, lnz:unz2)) * idzh_grid(lnz1:unz1)
                enddo
                enddo
                grad_deltab2(3, :, :, lnz) =  (-3.0*deltab2(:, :, lnz) + &
                                                4.0*deltab2(:, :, lnz1) - &
                                                    deltab2(:, :, lnz2)) * idzh_grid(lnz)
                grad_deltab2(3, :, :, unz) =   (3.0*deltab2(:, :, unz) - &
                                                4.0*deltab2(:, :, unz1) + &
                                                    deltab2(:, :, unz2)) * idzh_grid(unz)
            endif
        endif
        if (mpi_rank == master) then
            write(*, "(A)") "Finished calculating gradients of magnetic fluctuation."
        endif
    end subroutine calc_grad_deltab2_nonuniform

    !---------------------------------------------------------------------------
    !< Calculate the gradients of turbulence correlation length.
    !< Args:
    !<  var_flag: indicating which set of variables.
    !<            0 for grad_correl1 and other numbers for grad_correl2.
    !---------------------------------------------------------------------------
    subroutine calc_grad_correl_length(var_flag)
        use mhd_config_module, only: mhd_config
        use mpi_module
        implicit none
        integer, intent(in) :: var_flag
        real(dp) :: idxh, idyh, idzh
        integer :: unx, uny, unz, lnx, lny, lnz
        integer :: unx1, unx2, uny1, uny2, unz1, unz2
        integer :: lnx1, lnx2, lny1, lny2, lnz1, lnz2
        idxh = 0.5_dp / mhd_config%dx
        idyh = 0.5_dp / mhd_config%dy
        idzh = 0.5_dp / mhd_config%dz
        unx = ubound(lc1, 1)
        uny = ubound(lc1, 2)
        unz = ubound(lc1, 3)
        lnx = lbound(lc1, 1)
        lny = lbound(lc1, 2)
        lnz = lbound(lc1, 3)
        unx1 = unx - 1
        unx2 = unx - 2
        uny1 = uny - 1
        uny2 = uny - 2
        unz1 = unz - 1
        unz2 = unz - 2
        lnx1 = lnx + 1
        lnx2 = lnx + 2
        lny1 = lny + 1
        lny2 = lny + 2
        lnz1 = lnz + 1
        lnz2 = lnz + 2
        if (var_flag == 0) then
            ! d/dx
            grad_correl1(1, lnx1:unx1, :, :) = (lc1(lnx2:unx, :, :) - &
                                                lc1(lnx:unx2, :, :)) * idxh
            grad_correl1(1, lnx, :, :) =  (-3.0*lc1(lnx, :, :) + &
                                            4.0*lc1(lnx1, :, :) - &
                                                lc1(lnx2, :, :)) * idxh
            grad_correl1(1, unx, :, :) =   (3.0*lc1(unx, :, :) - &
                                            4.0*lc1(unx1, :, :) + &
                                                lc1(unx2, :, :)) * idxh

            ! d/dy
            if (uny > lny) then
                grad_correl1(2, :, lny1:uny1, :) = (lc1(:, lny2:uny, :) - &
                                                    lc1(:, lny:uny2, :)) * idyh
                grad_correl1(2, :, lny, :) =  (-3.0*lc1(:, lny, :) + &
                                                4.0*lc1(:, lny1, :) - &
                                                    lc1(:, lny2, :)) * idyh
                grad_correl1(2, :, uny, :) =   (3.0*lc1(:, uny, :) - &
                                                4.0*lc1(:, uny1, :) + &
                                                    lc1(:, uny2, :)) * idyh
            endif

            ! d/dz
            if (unz > lnz) then
                grad_correl1(3, :, :, lnz1:unz1) = (lc1(:, :, lnz2:unz) - &
                                                    lc1(:, :, lnz:unz2)) * idzh
                grad_correl1(3, :, :, lnz) =  (-3.0*lc1(:, :, lnz) + &
                                                4.0*lc1(:, :, lnz1) - &
                                                    lc1(:, :, lnz2)) * idzh
                grad_correl1(3, :, :, unz) =   (3.0*lc1(:, :, unz) - &
                                                4.0*lc1(:, :, unz1) + &
                                                    lc1(:, :, unz2)) * idzh
            endif
        else
            ! d/dx
            grad_correl2(1, lnx1:unx1, :, :) = (lc2(lnx2:unx, :, :) - &
                                                lc2(lnx:unx2, :, :)) * idxh
            grad_correl2(1, lnx, :, :) =  (-3.0*lc2(lnx, :, :) + &
                                            4.0*lc2(lnx1, :, :) - &
                                                lc2(lnx2, :, :)) * idxh
            grad_correl2(1, unx, :, :) =   (3.0*lc2(unx, :, :) - &
                                            4.0*lc2(unx1, :, :) + &
                                                lc2(unx2, :, :)) * idxh

            ! d/dy
            if (uny > lny) then
                grad_correl2(2, :, lny1:uny1, :) = (lc2(:, lny2:uny, :) - &
                                                    lc2(:, lny:uny2, :)) * idyh
                grad_correl2(2, :, lny, :) =  (-3.0*lc2(:, lny, :) + &
                                                4.0*lc2(:, lny1, :) - &
                                                    lc2(:, lny2, :)) * idyh
                grad_correl2(2, :, uny, :) =   (3.0*lc2(:, uny, :) - &
                                                4.0*lc2(:, uny1, :) + &
                                                    lc2(:, uny2, :)) * idyh
            endif

            ! d/dz
            if (unz > lnz) then
                grad_correl2(3, :, :, lnz1:unz1) = (lc2(:, :, lnz2:unz) - &
                                                    lc2(:, :, lnz:unz2)) * idzh
                grad_correl2(3, :, :, lnz) =  (-3.0*lc2(:, :, lnz) + &
                                                4.0*lc2(:, :, lnz1) - &
                                                    lc2(:, :, lnz2)) * idzh
                grad_correl2(3, :, :, unz) =   (3.0*lc2(:, :, unz) - &
                                                4.0*lc2(:, :, unz1) + &
                                                    lc2(:, :, unz2)) * idzh
            endif
        endif
        if (mpi_rank == master) then
            write(*, "(A)") "Finished calculating gradients of turbulence correlation length."
        endif
    end subroutine calc_grad_correl_length

    !---------------------------------------------------------------------------
    !< Calculate the gradients of turbulence correlation length for nonuniform
    !< grid
    !< Args:
    !<  var_flag: indicating which set of variables.
    !<            0 for grad_correl1 and other numbers for grad_correl2.
    !---------------------------------------------------------------------------
    subroutine calc_grad_correl_length_nonuniform(var_flag)
        use mhd_config_module, only: mhd_config
        use mpi_module
        implicit none
        integer, intent(in) :: var_flag
        integer :: unx, uny, unz, lnx, lny, lnz
        integer :: unx1, unx2, uny1, uny2, unz1, unz2
        integer :: lnx1, lnx2, lny1, lny2, lnz1, lnz2
        integer :: ix, iy, iz
        unx = ubound(lc1, 1)
        uny = ubound(lc1, 2)
        unz = ubound(lc1, 3)
        lnx = lbound(lc1, 1)
        lny = lbound(lc1, 2)
        lnz = lbound(lc1, 3)
        unx1 = unx - 1
        unx2 = unx - 2
        uny1 = uny - 1
        uny2 = uny - 2
        unz1 = unz - 1
        unz2 = unz - 2
        lnx1 = lnx + 1
        lnx2 = lnx + 2
        lny1 = lny + 1
        lny2 = lny + 2
        lnz1 = lnz + 1
        lnz2 = lnz + 2
        if (var_flag == 0) then
            ! d/dx
            do iz = lnz, unz
            do iy = lny, uny
                grad_correl1(1, lnx1:unx1, iy, iz) = &
                    (lc1(lnx2:unx, iy, iz) - &
                     lc1(lnx:unx2, iy, iz)) * idxh_grid(lnx1:unx1)
            enddo
            enddo
            grad_correl1(1, lnx, :, :) =  (-3.0*lc1(lnx, :, :) + &
                                            4.0*lc1(lnx1, :, :) - &
                                                lc1(lnx2, :, :)) * idxh_grid(lnx)
            grad_correl1(1, unx, :, :) =   (3.0*lc1(unx, :, :) - &
                                            4.0*lc1(unx1, :, :) + &
                                                lc1(unx2, :, :)) * idxh_grid(unx)

            ! d/dy
            if (uny > lny) then
                do iz = lnz, unz
                do ix = lnx, unx
                    grad_correl1(2, ix, lny1:uny1, iz) = &
                        (lc1(ix, lny2:uny, iz) - &
                         lc1(ix, lny:uny2, iz)) * idyh_grid(lny1:uny1)
                enddo
                enddo
                grad_correl1(2, :, lny, :) =  (-3.0*lc1(:, lny, :) + &
                                                4.0*lc1(:, lny1, :) - &
                                                    lc1(:, lny2, :)) * idyh_grid(lny)
                grad_correl1(2, :, uny, :) =   (3.0*lc1(:, uny, :) - &
                                                4.0*lc1(:, uny1, :) + &
                                                    lc1(:, uny2, :)) * idyh_grid(uny)
            endif

            ! d/dz
            if (unz > lnz) then
                do iy = lny, uny
                do ix = lnx, unx
                    grad_correl1(3, ix, iy, lnz1:unz1) = &
                        (lc1(ix, iy, lnz2:unz) - &
                         lc1(ix, iy, lnz:unz2)) * idzh_grid(lnz1:unz1)
                enddo
                enddo
                grad_correl1(3, :, :, lnz) =  (-3.0*lc1(:, :, lnz) + &
                                                4.0*lc1(:, :, lnz1) - &
                                                    lc1(:, :, lnz2)) * idzh_grid(lnz)
                grad_correl1(3, :, :, unz) =   (3.0*lc1(:, :, unz) - &
                                                4.0*lc1(:, :, unz1) + &
                                                    lc1(:, :, unz2)) * idzh_grid(unz)
            endif
        else
            ! d/dx
            do iz = lnz, unz
            do iy = lny, uny
                grad_correl2(1, lnx1:unx1, iy, iz) = &
                    (lc2(lnx2:unx, iy, iz) - &
                     lc2(lnx:unx2, iy, iz)) * idxh_grid(lnx1:unx1)
            enddo
            enddo
            grad_correl2(1, lnx, :, :) =  (-3.0*lc2(lnx, :, :) + &
                                            4.0*lc2(lnx1, :, :) - &
                                                lc2(lnx2, :, :)) * idxh_grid(lnx)
            grad_correl2(1, unx, :, :) =   (3.0*lc2(unx, :, :) - &
                                            4.0*lc2(unx1, :, :) + &
                                                lc2(unx2, :, :)) * idxh_grid(unx)

            ! d/dy
            if (uny > lny) then
                do iz = lnz, unz
                do ix = lnx, unx
                    grad_correl2(2, ix, lny1:uny1, iz) = &
                        (lc2(ix, lny2:uny, iz) - &
                         lc2(ix, lny:uny2, iz)) * idyh_grid(lny1:uny1)
                enddo
                enddo
                grad_correl2(2, :, lny, :) =  (-3.0*lc2(:, lny, :) + &
                                                4.0*lc2(:, lny1, :) - &
                                                    lc2(:, lny2, :)) * idyh_grid(lny)
                grad_correl2(2, :, uny, :) =   (3.0*lc2(:, uny, :) - &
                                                4.0*lc2(:, uny1, :) + &
                                                    lc2(:, uny2, :)) * idyh_grid(uny)
            endif

            ! d/dz
            if (unz > lnz) then
                do iy = lny, uny
                do ix = lnx, unx
                    grad_correl2(3, ix, iy, lnz1:unz1) = &
                        (lc2(ix, iy, lnz2:unz) - &
                         lc2(ix, iy, lnz:unz2)) * idzh_grid(lnz1:unz1)
                enddo
                enddo
                grad_correl2(3, :, :, lnz) =  (-3.0*lc2(:, :, lnz) + &
                                                4.0*lc2(:, :, lnz1) - &
                                                    lc2(:, :, lnz2)) * idzh_grid(lnz)
                grad_correl2(3, :, :, unz) =   (3.0*lc2(:, :, unz) - &
                                                4.0*lc2(:, :, unz1) + &
                                                    lc2(:, :, unz2)) * idzh_grid(unz)
            endif
        endif
        if (mpi_rank == master) then
            write(*, "(A)") "Finished calculating gradients of turbulence correlation length."
        endif
    end subroutine calc_grad_correl_length_nonuniform

    !---------------------------------------------------------------------------
    !< Interpolate the MHD fields and their gradients on one position
    !< nz is assumed to be 1.
    !< Args:
    !<  pos: the lower-left corner of the grid in grid indices.
    !<  weights: for linear interpolation
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !---------------------------------------------------------------------------
    subroutine interp_fields(pos, weights, rt)
        implicit none
        integer, dimension(3), intent(in) :: pos
        real(dp), dimension(8), intent(in) :: weights
        real(dp), intent(in) :: rt
        integer :: ix, iy, iz, i, j, k

        ix = pos(1)
        iy = pos(2)
        iz = pos(3)

        fields = 0.0_dp
        fields1 = 0.0_dp
        gradf = 0.0_dp
        gradf1 = 0.0_dp
        if (ndim_field .eq. 1) then
            do i = 0, 1
                fields = fields + f_array1(:, ix+i, 1, 1) * weights(i+1)
                fields1 = fields1 + f_array2(:, ix+i, 1, 1) * weights(i+1)
                gradf = gradf + fgrad_array1(:, ix+i, 1, 1) * weights(i+1)
                gradf1 = gradf1 + fgrad_array2(:, ix+i, 1, 1) * weights(i+1)
            enddo
        else if (ndim_field .eq. 2) then
            do j = 0, 1
            do i = 0, 1
                fields = fields + &
                    f_array1(:, ix+i, iy+j, 1) * weights(j*2+i+1)
                fields1 = fields1 + &
                    f_array2(:, ix+i, iy+j, 1) * weights(j*2+i+1)
                gradf = gradf + &
                    fgrad_array1(:, ix+i, iy+j, 1) * weights(j*2+i+1)
                gradf1 = gradf1 + &
                    fgrad_array2(:, ix+i, iy+j, 1) * weights(j*2+i+1)
            enddo
            enddo
        else
            do k = 0, 1
            do j = 0, 1
            do i = 0, 1
                fields = fields + &
                    f_array1(:, ix+i, iy+j, iz+k) * weights(k*4+j*2+i+1)
                fields1 = fields1 + &
                    f_array2(:, ix+i, iy+j, iz+k) * weights(k*4+j*2+i+1)
                gradf = gradf + &
                    fgrad_array1(:, ix+i, iy+j, iz+k) * weights(k*4+j*2+i+1)
                gradf1 = gradf1 + &
                    fgrad_array2(:, ix+i, iy+j, iz+k) * weights(k*4+j*2+i+1)
            enddo
            enddo
            enddo
        endif
        !< Time interpolation
        fields = fields * (1.0 - rt) + fields1 * rt
        gradf = gradf * (1.0 - rt) + gradf1 * rt
    end subroutine interp_fields

    !---------------------------------------------------------------------------
    !< Interpolate the magnetic fluctuation at one position
    !< nz is assumed to be 1.
    !< Args:
    !<  pos: the lower-left corner of the grid in grid indices.
    !<  weights: for linear interpolation
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !---------------------------------------------------------------------------
    subroutine interp_magnetic_fluctuation(pos, weights, rt)
        implicit none
        integer, dimension(3), intent(in) :: pos
        real(dp), dimension(8), intent(in) :: weights
        real(dp), intent(in) :: rt
        integer :: ix, iy, iz, i, j, k

        ix = pos(1)
        iy = pos(2)
        iz = pos(3)

        db2 = 0.0_dp
        db2_1 = 0.0_dp
        grad_db = 0.0_dp
        grad_db1 = 0.0_dp
        if (ndim_field .eq. 1) then
            do i = 0, 1
                db2 = db2 + deltab1(ix+i, 1, 1) * weights(i+1)
                db2_1 = db2_1 + deltab2(ix+i, 1, 1) * weights(i+1)
                grad_db = grad_db + grad_deltab1(:, ix+i, 1, 1) * weights(i+1)
                grad_db1 = grad_db1 + grad_deltab2(:, ix+i, 1, 1) * weights(i+1)
            enddo
        else if (ndim_field .eq. 2) then
            do j = 0, 1
            do i = 0, 1
                db2 = db2 + deltab1(ix+i, iy+j, 1) * weights(j*2+i+1)
                db2_1 = db2_1 + deltab2(ix+i, iy+j, 1) * weights(j*2+i+1)
                grad_db = grad_db + &
                    grad_deltab1(:, ix+i, iy+j, 1) * weights(j*2+i+1)
                grad_db1 = grad_db1 + &
                    grad_deltab2(:, ix+i, iy+j, 1) * weights(j*2+i+1)
            enddo
            enddo
        else
            do k = 0, 1
            do j = 0, 1
            do i = 0, 1
                db2 = db2 + deltab1(ix+i, iy+j, iz+k) * weights(k*4+j*2+i+1)
                db2_1 = db2_1 + deltab2(ix+i, iy+j, iz+k) * weights(k*4+j*2+i+1)
                grad_db = grad_db + &
                    grad_deltab1(:, ix+i, iy+j, iz+k) * weights(k*4+j*2+i+1)
                grad_db1 = grad_db1 + &
                    grad_deltab2(:, ix+i, iy+j, iz+k) * weights(k*4+j*2+i+1)
            enddo
            enddo
            enddo
        endif
        !< Time interpolation
        db2 = db2 * (1.0 - rt) + db2_1 * rt
        grad_db = grad_db * (1.0 - rt) + grad_db1 * rt
    end subroutine interp_magnetic_fluctuation

    !---------------------------------------------------------------------------
    !< Interpolate the turbulence correlation length at one position
    !< nz is assumed to be 1.
    !< Args:
    !<  pos: the lower-left corner of the grid in grid indices.
    !<  weights: for linear interpolation
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !---------------------------------------------------------------------------
    subroutine interp_correlation_length(pos, weights, rt)
        implicit none
        integer, dimension(3), intent(in) :: pos
        real(dp), dimension(8), intent(in) :: weights
        real(dp), intent(in) :: rt
        integer :: ix, iy, iz, i, j, k

        ix = pos(1)
        iy = pos(2)
        iz = pos(3)

        lc0_1 = 0.0_dp
        lc0_2 = 0.0_dp
        grad_lc = 0.0_dp
        grad_lc1 = 0.0_dp
        if (ndim_field .eq. 1) then
            do i = 0, 1
                lc0_1 = lc0_1 + lc1(ix+i, 1, 1) * weights(i+1)
                lc0_2 = lc0_2 + lc2(ix+i, 1, 1) * weights(i+1)
                grad_lc = grad_lc + grad_correl1(:, ix+i, 1, 1) * weights(i+1)
                grad_lc1 = grad_lc1 + grad_correl2(:, ix+i, 1, 1) * weights(i+1)
            enddo
        else if (ndim_field .eq. 2) then
            do j = 0, 1
            do i = 0, 1
                lc0_1 = lc0_1 + lc1(ix+i, iy+j, 1) * weights(j*2+i+1)
                lc0_2 = lc0_2 + lc2(ix+i, iy+j, 1) * weights(j*2+i+1)
                grad_lc = grad_lc + &
                    grad_correl1(:, ix+i, iy+j, 1) * weights(j*2+i+1)
                grad_lc1 = grad_lc1 + &
                    grad_correl2(:, ix+i, iy+j, 1) * weights(j*2+i+1)
            enddo
            enddo
        else
            do k = 0, 1
            do j = 0, 1
            do i = 0, 1
                lc0_1 = lc0_1 + lc1(ix+i, iy+j, iz+k) * weights(k*4+j*2+i+1)
                lc0_2 = lc0_2 + lc2(ix+i, iy+j, iz+k) * weights(k*4+j*2+i+1)
                grad_lc = grad_lc + &
                    grad_correl1(:, ix+i, iy+j, iz+k) * weights(k*4+j*2+i+1)
                grad_lc1 = grad_lc1 + &
                    grad_correl2(:, ix+i, iy+j, iz+k) * weights(k*4+j*2+i+1)
            enddo
            enddo
            enddo
        endif
        !< Time interpolation
        lc = lc0_1 * (1.0 - rt) + lc0_2 * rt
        grad_lc = grad_lc * (1.0 - rt) + grad_lc1 * rt
    end subroutine interp_correlation_length

    !<--------------------------------------------------------------------------
    !< Copy fields for usage in the next time interval
    !<--------------------------------------------------------------------------
    subroutine copy_fields
        implicit none
        f_array1 = f_array2
        fgrad_array1 = fgrad_array2
    end subroutine copy_fields

    !<--------------------------------------------------------------------------
    !< Copy magnetic fluctuation data for usage in the next time interval
    !<--------------------------------------------------------------------------
    subroutine copy_magnetic_fluctuation
        implicit none
        deltab1 = deltab2
        grad_deltab1 = grad_deltab2
    end subroutine copy_magnetic_fluctuation

    !<--------------------------------------------------------------------------
    !< Copy turbulence correlation length for usage in the next time interval
    !<--------------------------------------------------------------------------
    subroutine copy_correlation_length
        implicit none
        lc1 = lc2
        grad_correl1 = grad_correl2
    end subroutine copy_correlation_length

    !---------------------------------------------------------------------------
    !< Initialize shock x-position indices
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine init_shock_xpos(interp_flag)
        implicit none
        integer, intent(in) :: interp_flag
        integer :: ny, nz

        ny = fconfig%ny
        nz = fconfig%nz

        if (ndim_field == 1) then
            allocate(shock_xpos1(1, 1))
        else if (ndim_field == 2) then
            allocate(shock_xpos1(-1:ny+2, 1))
        else
            allocate(shock_xpos1(-1:ny+2, -1:nz+2))
        endif
        shock_xpos1 = 0
        ! Next time step
        if (interp_flag == 1) then
            if (ndim_field == 1) then
                allocate(shock_xpos2(1, 1))
            else if (ndim_field == 2) then
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
    !<  interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine locate_shock_xpos(interp_flag)
        use mpi_module
        implicit none
        integer, intent(in) :: interp_flag
        if (ndim_field == 1) then
            shock_xpos1(1, 1) = maxloc(abs(fgrad_array1(1, :, 1, 1)), dim=1)
        else if (ndim_field == 2) then
            shock_xpos1(:, 1) = maxloc(abs(fgrad_array1(1, :, :, 1)), dim=1)
        else
            shock_xpos1(:, :) = maxloc(abs(fgrad_array1(1, :, :, :)), dim=1)
        endif
        if (interp_flag == 1) then
            if (ndim_field == 1) then
                shock_xpos2(1, 1) = maxloc(abs(fgrad_array2(1, :, 1, 1)), dim=1)
            else if (ndim_field == 2) then
                shock_xpos2(:, 1) = maxloc(abs(fgrad_array2(1, :, :, 1)), dim=1)
            else
                shock_xpos2(:, :) = maxloc(abs(fgrad_array2(1, :, :, :)), dim=1)
            endif
        endif
    end subroutine locate_shock_xpos

    !---------------------------------------------------------------------------
    !< Interpolate shock location
    !< nz is assumed to be 1. We assume shock is propagating along x-direction
    !< Args:
    !<  pos: the lower-left corner of the grid in grid indices.
    !<  weights: for linear interpolation
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !---------------------------------------------------------------------------
    function interp_shock_location(pos, weights, rt) result(shock_xpos)
        implicit none
        integer, dimension(2), intent(in) :: pos
        real(dp), dimension(4), intent(in) :: weights
        real(dp), intent(in) :: rt
        real(dp) :: sx1, sx2, shock_xpos
        integer :: iy, iz, j, k

        iy = pos(1)
        iz = pos(2)

        if (ndim_field .eq. 1) then
            sx1 = shock_xpos1(1, 1)
            sx2 = shock_xpos2(1, 1)
        else if (ndim_field .eq. 2) then
            do j = 0, 1
                sx1 = sx1 + shock_xpos1(iy+j-1, 1) * weights(j+1)
                sx2 = sx2 + shock_xpos2(iy+j-1, 1) * weights(j+1)
            enddo
        else
            do k = 0, 1
            do j = 0, 1
                sx1 = sx1 + shock_xpos1(iy+j-1, iz+k-1) * weights(k*2+j+1)
                sx2 = sx2 + shock_xpos2(iy+j-1, iz+k-1) * weights(k*2+j+1)
            enddo
            enddo
        endif
        !< Time interpolation
        shock_xpos = sx2 * (1.0 - rt) + sx1 * rt
    end function interp_shock_location

    !---------------------------------------------------------------------------
    !< Initialize grid positions
    !---------------------------------------------------------------------------
    subroutine init_grid_positions
        implicit none
        integer :: nx, ny, nz, nxg, nyg, nzg

        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz
        nxg = fconfig%nxg
        nyg = fconfig%nyg
        nzg = fconfig%nzg

        allocate(xpos_global(nxg))
        allocate(ypos_global(nyg))
        allocate(zpos_global(nzg))
        allocate(xpos_local(-1:nx+2))
        allocate(ypos_local(-1:ny+2))
        allocate(zpos_local(-1:nz+2))
        allocate(idxh_grid(-1:nx+2))
        allocate(idyh_grid(-1:ny+2))
        allocate(idzh_grid(-1:nz+2))
    end subroutine init_grid_positions

    !---------------------------------------------------------------------------
    !< Free grid positions
    !---------------------------------------------------------------------------
    subroutine free_grid_positions
        implicit none
        deallocate(xpos_local, ypos_local, zpos_local)
        deallocate(xpos_global, ypos_global, zpos_global)
        deallocate(idxh_grid, idyh_grid, idzh_grid)
    end subroutine free_grid_positions

    !---------------------------------------------------------------------------
    !< Read grid positions for simulations with non-uniform grids.
    !---------------------------------------------------------------------------
    subroutine read_global_grid_positions(dir_mhd_data)
        use mpi_module
        implicit none
        character(*), intent(in) :: dir_mhd_data
        character(len=256) :: filename
        integer :: fh

        if (mpi_rank == master) then
            fh = 33
            write(filename, "(A,A)") trim(dir_mhd_data)//'xpos.dat'
            open(unit=fh, file=filename, access='stream', status='unknown', &
                 form='unformatted', action='read')
            read(fh, pos=1) xpos_global
            close(fh)
            if (ndim_field > 1) then
                write(filename, "(A,A)") trim(dir_mhd_data)//'ypos.dat'
                open(unit=fh, file=filename, access='stream', status='unknown', &
                     form='unformatted', action='read')
                read(fh, pos=1) ypos_global
                close(fh)
            endif
            if (ndim_field == 3) then
                write(filename, "(A,A)") trim(dir_mhd_data)//'zpos.dat'
                open(unit=fh, file=filename, access='stream', status='unknown', &
                     form='unformatted', action='read')
                read(fh, pos=1) zpos_global
                close(fh)
            endif
        endif
        call MPI_BCAST(xpos_global, fconfig%nxg, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        if (ndim_field > 1) then
            call MPI_BCAST(ypos_global, fconfig%nyg, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        endif
        if (ndim_field == 3) then
            call MPI_BCAST(zpos_global, fconfig%nzg, MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        endif
    end subroutine read_global_grid_positions

    !---------------------------------------------------------------------------
    !< Calculate grid positions for simulations with uniform grids.
    !---------------------------------------------------------------------------
    subroutine calc_global_grid_positions
        use mpi_module
        use mhd_config_module, only: mhd_config
        implicit none
        integer :: i
        do i = 1, fconfig%nxg
            xpos_global(i) = mhd_config%dx*(i-3) + mhd_config%xmin
        enddo
        if (ndim_field > 1) then
            do i = 1, fconfig%nyg
                ypos_global(i) = mhd_config%dy*(i-3) + mhd_config%ymin
            enddo
        endif
        if (ndim_field == 3) then
            do i = 1, fconfig%nzg
                zpos_global(i) = mhd_config%dz*(i-3) + mhd_config%zmin
            enddo
        endif
    end subroutine calc_global_grid_positions

    !---------------------------------------------------------------------------
    !< Get local grid positions and inverse interval for calculating gradients
    !---------------------------------------------------------------------------
    subroutine set_local_grid_positions(dir_mhd_data)
        use mpi_module
        use mhd_config_module, only: uniform_grid_flag
        implicit none
        character(*), intent(in) :: dir_mhd_data
        integer :: i, ix_min, iy_min, iz_min, ix_max, iy_max, iz_max

        if (uniform_grid_flag) then
            call calc_global_grid_positions
        else
            call read_global_grid_positions(dir_mhd_data)
        endif

        ix_min = fconfig%ix_min
        iy_min = fconfig%iy_min
        iz_min = fconfig%iz_min
        ix_max = fconfig%ix_max
        iy_max = fconfig%iy_max
        iz_max = fconfig%iz_max

        ! Local grid positions
        xpos_local = xpos_global(ix_min:ix_max)
        ypos_local = ypos_global(iy_min:iy_max)
        zpos_local = zpos_global(iz_min:iz_max)

        ! Inverse interval for calculating gradients
        do i = 2, fconfig%nxf-1
            idxh_grid(i-2) = 1.0_dp / (xpos_global(ix_min+i) - xpos_global(ix_min+i-2))
        enddo
        idxh_grid(-1) = &
            1.0_dp / (-3*xpos_global(ix_min) + 4*xpos_global(ix_min+1) - xpos_global(ix_min+2))
        idxh_grid(fconfig%nx+2) = &
            1.0_dp / (3*xpos_global(ix_max) - 4*xpos_global(ix_max-1) + xpos_global(ix_max-2))

        if (ndim_field > 1) then
            do i = 2, fconfig%nyf-1
                idyh_grid(i-2) = 1.0_dp / (ypos_global(iy_min+i) - ypos_global(iy_min+i-2))
            enddo
            idyh_grid(-1) = &
                1.0_dp / (-3*ypos_global(iy_min) + 4*ypos_global(iy_min+1) - ypos_global(iy_min+2))
            idyh_grid(fconfig%ny+2) = &
                1.0_dp / (3*ypos_global(iy_max) - 4*ypos_global(iy_max-1) + ypos_global(iy_max-2))
        endif

        if (ndim_field == 3) then
            do i = 2, fconfig%nzf-1
                idzh_grid(i-2) = 1.0_dp / (zpos_global(iz_min+i) - zpos_global(iz_min+i-2))
            enddo
            idzh_grid(-1) = &
                1.0_dp / (-3*zpos_global(iz_min) + 4*zpos_global(iz_min+1) - zpos_global(iz_min+2))
            idzh_grid(fconfig%nz+2) = &
                1.0_dp / (3*zpos_global(iz_max) - 4*zpos_global(iz_max-1) + zpos_global(iz_max-2))
        endif
    end subroutine set_local_grid_positions

    !---------------------------------------------------------------------------
    !< Get the number of cells with current density |jz| > jz_min
    !< Args:
    !<  jz_min: the minimum jz
    !<  part_box: box to inject particles
    !---------------------------------------------------------------------------
    function get_ncells_large_jz(jz_min, spherical_coord_flag, part_box) result (ncells_large_jz)
        implicit none
        real(dp), intent(in) :: jz_min
        logical, intent(in) :: spherical_coord_flag
        real(dp), intent(in), dimension(6) :: part_box
        integer :: nx, ny, nz, ix, iy, iz, ncells_large_jz
        real(dp), allocatable, dimension(:, :, :) :: abs_jz
        logical :: inbox_x, inbox_y, inbox_z, inbox

        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz

        ncells_large_jz = 0
        allocate(abs_jz(nx, ny, nz))
        if (spherical_coord_flag) then
            do iz = 1, nz
                do iy = 1, ny
                    abs_jz(:, iy, iz) = &
                        abs(fgrad_array1(16, 1:nx, iy, iz) + &
                            (f_array1(6, 1:nx, iy, iz) - &
                             fgrad_array1(14, 1:nx, iy, iz)) / xpos_local(1:nx))
                enddo
            enddo
        else
            abs_jz = abs(fgrad_array1(16, 1:nx, 1:ny, 1:nz) - &
                         fgrad_array1(14, 1:nx, 1:ny, 1:nz))
        endif
        do iz = 1, nz
            if (ndim_field > 2) then
                inbox_z = zpos_local(iz) > part_box(3) .and. zpos_local(iz) < part_box(6)
            else
                inbox_z = .true.
            endif
            do iy = 1, ny
                if (ndim_field > 1) then
                    inbox_y = ypos_local(iy) > part_box(2) .and. ypos_local(iy) < part_box(5)
                else
                    inbox_y = .true.
                endif
                do ix = 1, nx
                    inbox_x = xpos_local(ix) > part_box(1) .and. xpos_local(ix) < part_box(4)
                    inbox = inbox_z .and. inbox_y .and. inbox_x
                    if (inbox .and. abs_jz(ix, iy, iz) > jz_min) then
                        ncells_large_jz = ncells_large_jz + 1
                    endif
                enddo
            enddo
        enddo
        deallocate(abs_jz)
    end function get_ncells_large_jz
end module mhd_data_parallel

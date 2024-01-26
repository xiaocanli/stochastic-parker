!*******************************************************************************
!< Module for dealing with MHD data in parallel
!*******************************************************************************
module mhd_data_parallel
    use constants, only: sp, dp
    use simulation_setup_module, only: fconfig, ndim_field
    use mpi_module
    implicit none
    private
    save

    public init_field_data, free_field_data, read_field_data_parallel, &
           calc_fields_gradients, calc_fields_gradients_nonuniform, interp_fields, &
           copy_fields, init_shock_xpos, free_shock_xpos, &
           locate_shock_xpos, interp_shock_location, init_magnetic_fluctuation, &
           free_magnetic_fluctuation, read_magnetic_fluctuation, &
           interp_magnetic_fluctuation, copy_magnetic_fluctuation, &
           calc_grad_deltab2, calc_grad_deltab2_nonuniform, &
           init_correlation_length, free_correlation_length, &
           read_correlation_length, interp_correlation_length, &
           copy_correlation_length, calc_grad_correl_length, &
           calc_grad_correl_length_nonuniform, &
           init_grid_positions, free_grid_positions, &
           set_local_grid_positions, get_ncells_large_jz, &
           get_ncells_large_db2, get_ncells_large_divv
    public xpos_local, ypos_local, zpos_local
    public nfields, ngrads

    integer, parameter :: nfields=8
    integer, parameter :: ngrads=24

    !< 80 variables in total if there are two time frames
    real(sp), allocatable, dimension(:, :, :, :) :: farray1, farray2 ! Current,next
    real(sp), allocatable, dimension(:, :, :, :) :: deltab2_1, deltab2_2 ! Current,next
    real(sp), allocatable, dimension(:, :, :, :) :: lcorr1, lcorr2 ! Current,next
    !dir$ attributes align:128 :: farray1
    !dir$ attributes align:128 :: farray2
    !dir$ attributes align:16 :: deltab2_1
    !dir$ attributes align:16 :: deltab2_2
    !dir$ attributes align:16 :: lcorr1
    !dir$ attributes align:16 :: lcorr2

    logical :: time_interp

    integer, allocatable, dimension(:, :) :: shock_xpos1, shock_xpos2  ! Shock x-position indices

    ! Grid positions along each direction for non-uniform grid
    real(dp), allocatable, dimension(:) :: xpos_global, ypos_global, zpos_global
    real(dp), allocatable, dimension(:) :: xpos_local, ypos_local, zpos_local
    real(dp), allocatable, dimension(:) :: idxh_grid, idyh_grid, idzh_grid

    contains
    !---------------------------------------------------------------------------
    !< Initialize MHD field data arrays and their gradients
    !< Args:
    !<  time_interp_flag: whether two time steps are needed for interpolation
    !---------------------------------------------------------------------------
    subroutine init_field_data(time_interp_flag)
        implicit none
        integer, intent(in) :: time_interp_flag
        integer :: nx, ny, nz

        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz

        !< vx, vy, vz, rho, bx, by, bz, btot
        !< dvx_dx, dvx_dy, dvx_dz, dvy_dx, dvy_dy, dvy_dz,
        !< dvz_dx, dvz_dy, dvz_dz, drho_dx, drho_d, drho_dz,
        !< dbx_dx, dbx_dy, dbx_dz, dby_dx, dby_dy, dby_dz,
        !< dbz_dx, dbz_dy, dbz_dz, dbtot_dx, dbtot_dy, dbtot_dz
        if (ndim_field == 1) then
            allocate(farray1(nfields+ngrads, -1:nx+2, 1, 1))
        else if (ndim_field == 2) then
            allocate(farray1(nfields+ngrads, -1:nx+2, -1:ny+2, 1))
        else
            allocate(farray1(nfields+ngrads, -1:nx+2, -1:ny+2, -1:nz+2))
        endif
        farray1 = 0.0
        ! Next time step
        time_interp = .false.
        if (time_interp_flag == 1) then
            time_interp = .true.
            if (ndim_field == 1) then
                allocate(farray2(nfields+ngrads, -1:nx+2, 1, 1))
            else if (ndim_field == 2) then
                allocate(farray2(nfields+ngrads, -1:nx+2, -1:ny+2, 1))
            else
                allocate(farray2(nfields+ngrads, -1:nx+2, -1:ny+2, -1:nz+2))
            endif
            farray2 = 0.0
        endif
    end subroutine init_field_data

    !---------------------------------------------------------------------------
    !< Initialize the magnetic fluctuation, defined as deltaB**2 / B0**2, and
    !< its gradients.
    !---------------------------------------------------------------------------
    subroutine init_magnetic_fluctuation
        implicit none
        integer :: nx, ny, nz

        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz

        if (ndim_field == 1) then
            allocate(deltab2_1(4, -1:nx+2, 1, 1))
        else if (ndim_field == 2) then
            allocate(deltab2_1(4, -1:nx+2, -1:ny+2, 1))
        else
            allocate(deltab2_1(4, -1:nx+2, -1:ny+2, -1:nz+2))
        endif
        deltab2_1 = 1.0
        ! Next time step
        if (time_interp) then
            if (ndim_field == 1) then
                allocate(deltab2_2(4, -1:nx+2, 1, 1))
            else if (ndim_field == 2) then
                allocate(deltab2_2(4, -1:nx+2, -1:ny+2, 1))
            else
                allocate(deltab2_2(4, -1:nx+2, -1:ny+2, -1:nz+2))
            endif
            deltab2_2 = 1.0
        endif
    end subroutine init_magnetic_fluctuation

    !---------------------------------------------------------------------------
    !< Initialize the turbulence correlation length and its gradients
    !---------------------------------------------------------------------------
    subroutine init_correlation_length
        implicit none
        integer :: nx, ny, nz

        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz

        if (ndim_field == 1) then
            allocate(lcorr1(4, -1:nx+2, 1, 1))
        else if (ndim_field == 2) then
            allocate(lcorr1(4, -1:nx+2, -1:ny+2, 1))
        else
            allocate(lcorr1(4, -1:nx+2, -1:ny+2, -1:nz+2))
        endif
        lcorr1 = 1.0
        ! Next time step
        if (time_interp) then
            if (ndim_field == 1) then
                allocate(lcorr2(4, -1:nx+2, 1, 1))
            else if (ndim_field == 2) then
                allocate(lcorr2(4, -1:nx+2, -1:ny+2, 1))
            else
                allocate(lcorr2(4, -1:nx+2, -1:ny+2, -1:nz+2))
            endif
            lcorr2 = 1.0
        endif
    end subroutine init_correlation_length

    !---------------------------------------------------------------------------
    !< Free MHD field data arrays
    !---------------------------------------------------------------------------
    subroutine free_field_data
        implicit none
        deallocate(farray1)
        if (time_interp) then
            deallocate(farray2)
        endif
    end subroutine free_field_data

    !---------------------------------------------------------------------------
    !< Free the data array for magnetic fluctuation
    !---------------------------------------------------------------------------
    subroutine free_magnetic_fluctuation
        implicit none
        deallocate(deltab2_1)
        if (time_interp) then
            deallocate(deltab2_2)
        endif
    end subroutine free_magnetic_fluctuation

    !---------------------------------------------------------------------------
    !< Free the data array for turbulence correlation length
    !---------------------------------------------------------------------------
    subroutine free_correlation_length
        implicit none
        deallocate(lcorr1)
        if (time_interp) then
            deallocate(lcorr2)
        endif
    end subroutine free_correlation_length

    !---------------------------------------------------------------------------
    !< Read MHD field data in parallel
    !< Args:
    !<  filename: file name to get the data
    !<  var_flag: indicating which set of variables to save the data. 0 for
    !<            farray1 and other numbers for farray2.
    !---------------------------------------------------------------------------
    subroutine read_field_data_parallel(filename, var_flag)
        use mpi_io_module, only: set_mpi_datatype_real, set_mpi_info, fileinfo, &
            open_data_mpi_io, read_data_mpi_io
        implicit none
        character(*), intent(in) :: filename
        integer, intent(in) :: var_flag
        integer :: mpi_datatype, fh
        integer, dimension(4) :: sizes, subsizes, starts
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
        sizes(1) = nfields ! vx, vy, vz, rho, bx, by, bz, btot
        sizes(2) = fconfig%nxg
        sizes(3) = fconfig%nyg
        sizes(4) = fconfig%nzg
        subsizes(1) = nfields
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
                    read(fh, pos=1) farray1(:nfields, :, :, :)
                endif
                call MPI_BCAST(farray1(:nfields, :, :, :), product(sizes), &
                    MPI_REAL4, master, MPI_COMM_WORLD, ierr)
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            else
                if (mpi_rank == master) then
                    read(fh, pos=1) farray2(:nfields, :, :, :)
                endif
                call MPI_BCAST(farray2(:nfields, :, :, :), product(sizes), &
                    MPI_REAL4, master, MPI_COMM_WORLD, ierr)
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            endif
            close(fh)
        else
            if (mpi_cross_rank == master) then
                mpi_datatype = set_mpi_datatype_real(sizes, subsizes, starts)
                call set_mpi_info
                call open_data_mpi_io(filename, MPI_MODE_RDONLY, fileinfo, mpi_sub_comm, fh)
                disp = 0
                offset = 0
                if (var_flag == 0) then
                    call read_data_mpi_io(fh, mpi_datatype, subsizes, disp, &
                        offset, farray1(:nfields, :, :, :))
                else
                    call read_data_mpi_io(fh, mpi_datatype, subsizes, disp, &
                        offset, farray2(:nfields, :, :, :))
                endif
                call MPI_FILE_CLOSE(fh, ierror)
            endif
            if (var_flag == 0) then
                call MPI_BCAST(farray1(:nfields, :, :, :), product(subsizes), &
                    MPI_REAL4, master, mpi_cross_comm, ierr)
                call MPI_BARRIER(mpi_cross_comm, ierr)
            else
                call MPI_BCAST(farray2(:nfields, :, :, :), product(subsizes), &
                    MPI_REAL4, master, mpi_cross_comm, ierr)
                call MPI_BARRIER(mpi_cross_comm, ierr)
            endif
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
                    read(fh, pos=1) deltab2_1(1, :, :, :)
                endif
                call MPI_BCAST(deltab2_1(1, :, :, :), product(sizes), &
                    MPI_REAL4, master, MPI_COMM_WORLD, ierr)
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            else
                if (mpi_rank == master) then
                    read(fh, pos=1) deltab2_2(1, :, :, :)
                endif
                call MPI_BCAST(deltab2_2(1, :, :, :), product(sizes), &
                    MPI_REAL4, master, MPI_COMM_WORLD, ierr)
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            endif
            close(fh)
        else
            if (mpi_cross_rank == master) then
                mpi_datatype = set_mpi_datatype_real(sizes, subsizes, starts)
                call set_mpi_info
                call open_data_mpi_io(filename, MPI_MODE_RDONLY, fileinfo, mpi_sub_comm, fh)
                disp = 0
                offset = 0
                if (var_flag == 0) then
                    call read_data_mpi_io(fh, mpi_datatype, subsizes, &
                        disp, offset, deltab2_1(1, :, :, :))
                else
                    call read_data_mpi_io(fh, mpi_datatype, subsizes, &
                        disp, offset, deltab2_2(1, :, :, :))
                endif
                call MPI_FILE_CLOSE(fh, ierror)
            endif
            if (var_flag == 0) then
                call MPI_BCAST(deltab2_1(1, :, :, :), product(subsizes), &
                    MPI_REAL4, master, mpi_cross_comm, ierr)
                call MPI_BARRIER(mpi_cross_comm, ierr)
            else
                call MPI_BCAST(deltab2_2(1, :, :, :), product(subsizes), &
                    MPI_REAL4, master, mpi_cross_comm, ierr)
                call MPI_BARRIER(mpi_cross_comm, ierr)
            endif
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
                    read(fh, pos=1) lcorr1(1, :, :, :)
                endif
                call MPI_BCAST(lcorr1(1, :, :, :), product(sizes), &
                    MPI_REAL4, master, MPI_COMM_WORLD, ierr)
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            else
                if (mpi_rank == master) then
                    read(fh, pos=1) lcorr2(1, :, :, :)
                endif
                call MPI_BCAST(lcorr2(1, :, :, :), product(sizes), &
                    MPI_REAL4, master, MPI_COMM_WORLD, ierr)
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            endif
            close(fh)
        else
            if (mpi_cross_rank == master) then
                mpi_datatype = set_mpi_datatype_real(sizes, subsizes, starts)
                call set_mpi_info
                call open_data_mpi_io(filename, MPI_MODE_RDONLY, fileinfo, mpi_sub_comm, fh)
                disp = 0
                offset = 0
                if (var_flag == 0) then
                    call read_data_mpi_io(fh, mpi_datatype, subsizes, &
                        disp, offset, lcorr1(1, :, :, :))
                else
                    call read_data_mpi_io(fh, mpi_datatype, subsizes, &
                        disp, offset, lcorr2(1, :, :, :))
                endif
                call MPI_FILE_CLOSE(fh, ierror)
            endif
            if (var_flag == 0) then
                call MPI_BCAST(lcorr1(1, :, :, :), product(subsizes), &
                    MPI_REAL4, master, mpi_cross_comm, ierr)
                call MPI_BARRIER(mpi_cross_comm, ierr)
            else
                call MPI_BCAST(lcorr2(1, :, :, :), product(subsizes), &
                    MPI_REAL4, master, mpi_cross_comm, ierr)
                call MPI_BARRIER(mpi_cross_comm, ierr)
            endif
        endif
        if (mpi_rank == master) then
            write(*, "(A)") "Finished reading turbulence correlation length."
        endif
    end subroutine read_correlation_length

    !---------------------------------------------------------------------------
    !< Calculate the gradients of the MHD data arrays.
    !< Args:
    !<  var_flag: indicating which set of variables.
    !---------------------------------------------------------------------------
    subroutine calc_fields_gradients(var_flag)
        use mhd_config_module, only: mhd_config
        implicit none
        integer, intent(in) :: var_flag
        real(dp) :: idxh, idyh, idzh
        integer :: unx, uny, unz, lnx, lny, lnz
        integer :: unx1, unx2, uny1, uny2, unz1, unz2
        integer :: lnx1, lnx2, lny1, lny2, lnz1, lnz2
        idxh = 0.5_dp / mhd_config%dx
        idyh = 0.5_dp / mhd_config%dy
        idzh = 0.5_dp / mhd_config%dz
        unx = ubound(farray1, 2)
        uny = ubound(farray1, 3)
        unz = ubound(farray1, 4)
        lnx = lbound(farray1, 2)
        lny = lbound(farray1, 3)
        lnz = lbound(farray1, 4)
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
            farray1(nfields+1::3, lnx1:unx1, :, :) = (farray1(:nfields, lnx2:unx, :, :) - &
                                                      farray1(:nfields, lnx:unx2, :, :)) * idxh
            farray1(nfields+1::3, lnx, :, :) =  (-3.0*farray1(:nfields, lnx, :, :) + &
                                                  4.0*farray1(:nfields, lnx1, :, :) - &
                                                      farray1(:nfields, lnx2, :, :)) * idxh
            farray1(nfields+1::3, unx, :, :) =   (3.0*farray1(:nfields, unx, :, :) - &
                                                  4.0*farray1(:nfields, unx1, :, :) + &
                                                      farray1(:nfields, unx2, :, :)) * idxh

            ! d/dy
            if (uny > lny) then
                farray1(nfields+2::3, :, lny1:uny1, :) = (farray1(:nfields, :, lny2:uny, :) - &
                                                          farray1(:nfields, :, lny:uny2, :)) * idyh
                farray1(nfields+2::3, :, lny, :) =  (-3.0*farray1(:nfields, :, lny, :) + &
                                                      4.0*farray1(:nfields, :, lny1, :) - &
                                                          farray1(:nfields, :, lny2, :)) * idyh
                farray1(nfields+2::3, :, uny, :) =   (3.0*farray1(:nfields, :, uny, :) - &
                                                      4.0*farray1(:nfields, :, uny1, :) + &
                                                          farray1(:nfields, :, uny2, :)) * idyh
            endif

            ! d/dz
            if (unz > lnz) then
                farray1(nfields+3::3, :, :, lnz1:unz1) = (farray1(:nfields, :, :, lnz2:unz) - &
                                                          farray1(:nfields, :, :, lnz:unz2)) * idzh
                farray1(nfields+3::3, :, :, lnz) =  (-3.0*farray1(:nfields, :, :, lnz) + &
                                                      4.0*farray1(:nfields, :, :, lnz1) - &
                                                          farray1(:nfields, :, :, lnz2)) * idzh
                farray1(nfields+3::3, :, :, unz) =   (3.0*farray1(:nfields, :, :, unz) - &
                                                      4.0*farray1(:nfields, :, :, unz1) + &
                                                          farray1(:nfields, :, :, unz2)) * idzh
            endif
        else
            ! d/dx
            farray2(nfields+1::3, lnx1:unx1, :, :) = (farray2(:nfields, lnx2:unx, :, :) - &
                                                      farray2(:nfields, lnx:unx2, :, :)) * idxh
            farray2(nfields+1::3, lnx, :, :) =  (-3.0*farray2(:nfields, lnx, :, :) + &
                                                  4.0*farray2(:nfields, lnx1, :, :) - &
                                                      farray2(:nfields, lnx2, :, :)) * idxh
            farray2(nfields+1::3, unx, :, :) =   (3.0*farray2(:nfields, unx, :, :) - &
                                                  4.0*farray2(:nfields, unx1, :, :) + &
                                                      farray2(:nfields, unx2, :, :)) * idxh

            ! d/dy
            if (uny > lny) then
                farray2(nfields+2::3, :, lny1:uny1, :) = (farray2(:nfields, :, lny2:uny, :) - &
                                                          farray2(:nfields, :, lny:uny2, :)) * idyh
                farray2(nfields+2::3, :, lny, :) =  (-3.0*farray2(:nfields, :, lny, :) + &
                                                      4.0*farray2(:nfields, :, lny1, :) - &
                                                          farray2(:nfields, :, lny2, :)) * idyh
                farray2(nfields+2::3, :, uny, :) =   (3.0*farray2(:nfields, :, uny, :) - &
                                                      4.0*farray2(:nfields, :, uny1, :) + &
                                                          farray2(:nfields, :, uny2, :)) * idyh
            endif

            ! d/dz
            if (unz > lnz) then
                farray2(nfields+3::3, :, :, lnz1:unz1) = (farray2(:nfields, :, :, lnz2:unz) - &
                                                          farray2(:nfields, :, :, lnz:unz2)) * idzh
                farray2(nfields+3::3, :, :, lnz) =  (-3.0*farray2(:nfields, :, :, lnz) + &
                                                      4.0*farray2(:nfields, :, :, lnz1) - &
                                                          farray2(:nfields, :, :, lnz2)) * idzh
                farray2(nfields+3::3, :, :, unz) =   (3.0*farray2(:nfields, :, :, unz) - &
                                                      4.0*farray2(:nfields, :, :, unz1) + &
                                                          farray2(:nfields, :, :, unz2)) * idzh
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
    !---------------------------------------------------------------------------
    subroutine calc_fields_gradients_nonuniform(var_flag)
        use mhd_config_module, only: mhd_config
        implicit none
        integer, intent(in) :: var_flag
        integer :: unx, uny, unz, lnx, lny, lnz
        integer :: unx1, unx2, uny1, uny2, unz1, unz2
        integer :: lnx1, lnx2, lny1, lny2, lnz1, lnz2
        integer :: iv, ix, iy, iz
        unx = ubound(farray1, 2)
        uny = ubound(farray1, 3)
        unz = ubound(farray1, 4)
        lnx = lbound(farray1, 2)
        lny = lbound(farray1, 3)
        lnz = lbound(farray1, 4)
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
                farray1(nfields+3*iv-2, lnx1:unx1, iy, iz) = &
                    (farray1(iv, lnx2:unx, iy, iz) - &
                     farray1(iv, lnx:unx2, iy, iz)) * idxh_grid(lnx1:unx1)
            enddo
            enddo
            enddo
            farray1(nfields+1::3, lnx, :, :) =  (-3.0*farray1(:nfields, lnx, :, :) + &
                                                  4.0*farray1(:nfields, lnx1, :, :) - &
                                                      farray1(:nfields, lnx2, :, :)) * &
                                                      idxh_grid(lnx)
            farray1(nfields+1::3, unx, :, :) =   (3.0*farray1(:nfields, unx, :, :) - &
                                                  4.0*farray1(:nfields, unx1, :, :) + &
                                                      farray1(:nfields, unx2, :, :)) * &
                                                      idxh_grid(unx)

            ! d/dy
            if (uny > lny) then
                do iz = lnz, unz
                do ix = lnx, unx
                do iv = 1, nfields
                    farray1(nfields+3*iv-1, ix, lny1:uny1, iz) = &
                        (farray1(iv, ix, lny2:uny, iz) - &
                         farray1(iv, ix, lny:uny2, iz)) * idyh_grid(lny1:uny1)
                enddo
                enddo
                enddo
                farray1(nfields+2::3, :, lny, :) =  (-3.0*farray1(:nfields, :, lny, :) + &
                                                      4.0*farray1(:nfields, :, lny1, :) - &
                                                          farray1(:nfields, :, lny2, :)) * &
                                                          idyh_grid(lny)
                farray1(nfields+2::3, :, uny, :) =   (3.0*farray1(:nfields, :, uny, :) - &
                                                      4.0*farray1(:nfields, :, uny1, :) + &
                                                          farray1(:nfields, :, uny2, :)) * &
                                                          idyh_grid(uny)
            endif

            ! d/dz
            if (unz > lnz) then
                do iy = lny, uny
                do ix = lnx, unx
                do iv = 1, nfields
                    farray1(nfields+3*iv, ix, iy, lnz1:unz1) = &
                        (farray1(iv, ix, iy, lnz2:unz) - &
                         farray1(iv, ix, iy, lnz:unz2)) * idzh_grid(lnz1:unz1)
                enddo
                enddo
                enddo
                farray1(nfields+3::3, :, :, lnz) =  (-3.0*farray1(:nfields, :, :, lnz) + &
                                                      4.0*farray1(:nfields, :, :, lnz1) - &
                                                          farray1(:nfields, :, :, lnz2)) * &
                                                          idzh_grid(lnz)
                farray1(nfields+3::3, :, :, unz) =   (3.0*farray1(:nfields, :, :, unz) - &
                                                      4.0*farray1(:nfields, :, :, unz1) + &
                                                          farray1(:nfields, :, :, unz2)) * &
                                                          idzh_grid(unz)
            endif
        else
            ! d/dx
            do iz = lnz, unz
            do iy = lny, uny
            do iv = 1, nfields
                farray2(nfields+3*iv-2, lnx1:unx1, iy, iz) = &
                    (farray2(iv, lnx2:unx, iy, iz) - &
                     farray2(iv, lnx:unx2, iy, iz)) * idxh_grid(lnx1:unx1)
            enddo
            enddo
            enddo
            farray2(nfields+1::3, lnx, :, :) =  (-3.0*farray2(:nfields, lnx, :, :) + &
                                                  4.0*farray2(:nfields, lnx1, :, :) - &
                                                      farray2(:nfields, lnx2, :, :)) * &
                                                      idxh_grid(lnx)
            farray2(nfields+1::3, unx, :, :) =   (3.0*farray2(:nfields, unx, :, :) - &
                                                  4.0*farray2(:nfields, unx1, :, :) + &
                                                      farray2(:nfields, unx2, :, :)) * &
                                                      idxh_grid(unx)

            ! d/dy
            if (uny > lny) then
                do iz = lnz, unz
                do ix = lnx, unx
                do iv = 1, nfields
                    farray2(nfields+3*iv-1, ix, lny1:uny1, iz) = &
                        (farray2(iv, ix, lny2:uny, iz) - &
                         farray2(iv, ix, lny:uny2, iz)) * idyh_grid(lny1:uny1)
                enddo
                enddo
                enddo
                farray2(nfields+2::3, :, lny, :) =  (-3.0*farray2(:nfields, :, lny, :) + &
                                                      4.0*farray2(:nfields, :, lny1, :) - &
                                                          farray2(:nfields, :, lny2, :)) * &
                                                          idyh_grid(lny)
                farray2(nfields+2::3, :, uny, :) =   (3.0*farray2(:nfields, :, uny, :) - &
                                                      4.0*farray2(:nfields, :, uny1, :) + &
                                                          farray2(:nfields, :, uny2, :)) * &
                                                          idyh_grid(uny)
            endif

            ! d/dz
            if (unz > lnz) then
                do iy = lny, uny
                do ix = lnx, unx
                do iv = 1, nfields
                    farray2(nfields+3*iv, ix, iy, lnz1:unz1) = &
                        (farray2(iv, ix, iy, lnz2:unz) - &
                         farray2(iv, ix, iy, lnz:unz2)) * idzh_grid(lnz1:unz1)
                enddo
                enddo
                enddo
                farray2(nfields+3::3, :, :, lnz) =  (-3.0*farray2(:nfields, :, :, lnz) + &
                                                      4.0*farray2(:nfields, :, :, lnz1) - &
                                                          farray2(:nfields, :, :, lnz2)) * &
                                                          idzh_grid(lnz)
                farray2(nfields+3::3, :, :, unz) =   (3.0*farray2(:nfields, :, :, unz) - &
                                                      4.0*farray2(:nfields, :, :, unz1) + &
                                                          farray2(:nfields, :, :, unz2)) * &
                                                          idzh_grid(unz)
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
    !---------------------------------------------------------------------------
    subroutine calc_grad_deltab2(var_flag)
        use mhd_config_module, only: mhd_config
        implicit none
        integer, intent(in) :: var_flag
        real(dp) :: idxh, idyh, idzh
        integer :: unx, uny, unz, lnx, lny, lnz
        integer :: unx1, unx2, uny1, uny2, unz1, unz2
        integer :: lnx1, lnx2, lny1, lny2, lnz1, lnz2
        idxh = 0.5_dp / mhd_config%dx
        idyh = 0.5_dp / mhd_config%dy
        idzh = 0.5_dp / mhd_config%dz
        unx = ubound(deltab2_1, 2)
        uny = ubound(deltab2_1, 3)
        unz = ubound(deltab2_1, 4)
        lnx = lbound(deltab2_1, 2)
        lny = lbound(deltab2_1, 3)
        lnz = lbound(deltab2_1, 4)
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
            deltab2_1(2, lnx1:unx1, :, :) = (deltab2_1(1, lnx2:unx, :, :) - &
                                             deltab2_1(1, lnx:unx2, :, :)) * idxh
            deltab2_1(2, lnx, :, :) =  (-3.0*deltab2_1(1, lnx, :, :) + &
                                         4.0*deltab2_1(1, lnx1, :, :) - &
                                             deltab2_1(1, lnx2, :, :)) * idxh
            deltab2_1(2, unx, :, :) =   (3.0*deltab2_1(1, unx, :, :) - &
                                         4.0*deltab2_1(1, unx1, :, :) + &
                                             deltab2_1(1, unx2, :, :)) * idxh

            ! d/dy
            if (uny > lny) then
                deltab2_1(3, :, lny1:uny1, :) = (deltab2_1(1, :, lny2:uny, :) - &
                                                 deltab2_1(1, :, lny:uny2, :)) * idyh
                deltab2_1(3, :, lny, :) =  (-3.0*deltab2_1(1, :, lny, :) + &
                                             4.0*deltab2_1(1, :, lny1, :) - &
                                                 deltab2_1(1, :, lny2, :)) * idyh
                deltab2_1(3, :, uny, :) =   (3.0*deltab2_1(1, :, uny, :) - &
                                             4.0*deltab2_1(1, :, uny1, :) + &
                                                 deltab2_1(1, :, uny2, :)) * idyh
            endif

            ! d/dz
            if (unz > lnz) then
                deltab2_1(4, :, :, lnz1:unz1) = (deltab2_1(1, :, :, lnz2:unz) - &
                                                 deltab2_1(1, :, :, lnz:unz2)) * idzh
                deltab2_1(4, :, :, lnz) =  (-3.0*deltab2_1(1, :, :, lnz) + &
                                             4.0*deltab2_1(1, :, :, lnz1) - &
                                                 deltab2_1(1, :, :, lnz2)) * idzh
                deltab2_1(4, :, :, unz) =   (3.0*deltab2_1(1, :, :, unz) - &
                                             4.0*deltab2_1(1, :, :, unz1) + &
                                                 deltab2_1(1, :, :, unz2)) * idzh
            endif
        else
            ! d/dx
            deltab2_2(2, lnx1:unx1, :, :) = (deltab2_2(1, lnx2:unx, :, :) - &
                                             deltab2_2(1, lnx:unx2, :, :)) * idxh
            deltab2_2(2, lnx, :, :) =  (-3.0*deltab2_2(1, lnx, :, :) + &
                                         4.0*deltab2_2(1, lnx1, :, :) - &
                                             deltab2_2(1, lnx2, :, :)) * idxh
            deltab2_2(2, unx, :, :) =   (3.0*deltab2_2(1, unx, :, :) - &
                                         4.0*deltab2_2(1, unx1, :, :) + &
                                             deltab2_2(1, unx2, :, :)) * idxh

            ! d/dy
            if (uny > lny) then
                deltab2_2(3, :, lny1:uny1, :) = (deltab2_2(1, :, lny2:uny, :) - &
                                                 deltab2_2(1, :, lny:uny2, :)) * idyh
                deltab2_2(3, :, lny, :) =  (-3.0*deltab2_2(1, :, lny, :) + &
                                             4.0*deltab2_2(1, :, lny1, :) - &
                                                 deltab2_2(1, :, lny2, :)) * idyh
                deltab2_2(3, :, uny, :) =   (3.0*deltab2_2(1, :, uny, :) - &
                                             4.0*deltab2_2(1, :, uny1, :) + &
                                                 deltab2_2(1, :, uny2, :)) * idyh
            endif

            ! d/dz
            if (unz > lnz) then
                deltab2_2(4, :, :, lnz1:unz1) = (deltab2_2(1, :, :, lnz2:unz) - &
                                                 deltab2_2(1, :, :, lnz:unz2)) * idzh
                deltab2_2(4, :, :, lnz) =  (-3.0*deltab2_2(1, :, :, lnz) + &
                                             4.0*deltab2_2(1, :, :, lnz1) - &
                                                 deltab2_2(1, :, :, lnz2)) * idzh
                deltab2_2(4, :, :, unz) =   (3.0*deltab2_2(1, :, :, unz) - &
                                             4.0*deltab2_2(1, :, :, unz1) + &
                                                 deltab2_2(1, :, :, unz2)) * idzh
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
    !---------------------------------------------------------------------------
    subroutine calc_grad_deltab2_nonuniform(var_flag)
        use mhd_config_module, only: mhd_config
        implicit none
        integer, intent(in) :: var_flag
        integer :: unx, uny, unz, lnx, lny, lnz
        integer :: unx1, unx2, uny1, uny2, unz1, unz2
        integer :: lnx1, lnx2, lny1, lny2, lnz1, lnz2
        integer :: ix, iy, iz
        unx = ubound(deltab2_1, 2)
        uny = ubound(deltab2_1, 3)
        unz = ubound(deltab2_1, 4)
        lnx = lbound(deltab2_1, 2)
        lny = lbound(deltab2_1, 3)
        lnz = lbound(deltab2_1, 4)
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
                deltab2_1(2, lnx1:unx1, iy, iz) = &
                    (deltab2_1(1, lnx2:unx, iy, iz) - &
                     deltab2_1(1, lnx:unx2, iy, iz)) * idxh_grid(lnx1:unx1)
            enddo
            enddo
            deltab2_1(2, lnx, :, :) =  (-3.0*deltab2_1(1, lnx, :, :) + &
                                         4.0*deltab2_1(1, lnx1, :, :) - &
                                             deltab2_1(1, lnx2, :, :)) * idxh_grid(lnx)
            deltab2_1(2, unx, :, :) =   (3.0*deltab2_1(1, unx, :, :) - &
                                         4.0*deltab2_1(1, unx1, :, :) + &
                                             deltab2_1(1, unx2, :, :)) * idxh_grid(unx)

            ! d/dy
            if (uny > lny) then
                do iz = lnz, unz
                do ix = lnx, unx
                    deltab2_1(3, ix, lny1:uny1, iz) = &
                        (deltab2_1(1, ix, lny2:uny, iz) - &
                         deltab2_1(1, ix, lny:uny2, iz)) * idyh_grid(lny1:uny1)
                enddo
                enddo
                deltab2_1(3, :, lny, :) =  (-3.0*deltab2_1(1, :, lny, :) + &
                                             4.0*deltab2_1(1, :, lny1, :) - &
                                                 deltab2_1(1, :, lny2, :)) * idyh_grid(lny)
                deltab2_1(3, :, uny, :) =   (3.0*deltab2_1(1, :, uny, :) - &
                                             4.0*deltab2_1(1, :, uny1, :) + &
                                                 deltab2_1(1, :, uny2, :)) * idyh_grid(uny)
            endif

            ! d/dz
            if (unz > lnz) then
                do iy = lny, uny
                do ix = lnx, unx
                    deltab2_1(4, ix, iy, lnz1:unz1) = &
                        (deltab2_1(1, ix, iy, lnz2:unz) - &
                         deltab2_1(1, ix, iy, lnz:unz2)) * idzh_grid(lnz1:unz1)
                enddo
                enddo
                deltab2_1(4, :, :, lnz) =  (-3.0*deltab2_1(1, :, :, lnz) + &
                                             4.0*deltab2_1(1, :, :, lnz1) - &
                                                 deltab2_1(1, :, :, lnz2)) * idzh_grid(lnz)
                deltab2_1(4, :, :, unz) =   (3.0*deltab2_1(1, :, :, unz) - &
                                             4.0*deltab2_1(1, :, :, unz1) + &
                                                 deltab2_1(1, :, :, unz2)) * idzh_grid(unz)
            endif
        else
            ! d/dx
            do iz = lnz, unz
            do iy = lny, uny
                deltab2_2(2, lnx1:unx1, iy, iz) = &
                    (deltab2_2(1, lnx2:unx, iy, iz) - &
                     deltab2_2(1, lnx:unx2, iy, iz)) * idxh_grid(lnx1:unx1)
            enddo
            enddo
            deltab2_2(2, lnx, :, :) =  (-3.0*deltab2_2(1, lnx, :, :) + &
                                         4.0*deltab2_2(1, lnx1, :, :) - &
                                             deltab2_2(1, lnx2, :, :)) * idxh_grid(lnx)
            deltab2_2(2, unx, :, :) =   (3.0*deltab2_2(1, unx, :, :) - &
                                         4.0*deltab2_2(1, unx1, :, :) + &
                                             deltab2_2(1, unx2, :, :)) * idxh_grid(unx)

            ! d/dy
            if (uny > lny) then
                do iz = lnz, unz
                do ix = lnx, unx
                    deltab2_2(3, ix, lny1:uny1, iz) = &
                        (deltab2_2(1, ix, lny2:uny, iz) - &
                         deltab2_2(1, ix, lny:uny2, iz)) * idyh_grid(lny1:uny1)
                enddo
                enddo
                deltab2_2(3, :, lny, :) =  (-3.0*deltab2_2(1, :, lny, :) + &
                                             4.0*deltab2_2(1, :, lny1, :) - &
                                                 deltab2_2(1, :, lny2, :)) * idyh_grid(lny)
                deltab2_2(3, :, uny, :) =   (3.0*deltab2_2(1, :, uny, :) - &
                                             4.0*deltab2_2(1, :, uny1, :) + &
                                                 deltab2_2(1, :, uny2, :)) * idyh_grid(uny)
            endif

            ! d/dz
            if (unz > lnz) then
                do iy = lny, uny
                do ix = lnx, unx
                    deltab2_2(4, ix, iy, lnz1:unz1) = &
                        (deltab2_2(1, ix, iy, lnz2:unz) - &
                         deltab2_2(1, ix, iy, lnz:unz2)) * idzh_grid(lnz1:unz1)
                enddo
                enddo
                deltab2_2(4, :, :, lnz) =  (-3.0*deltab2_2(1, :, :, lnz) + &
                                             4.0*deltab2_2(1, :, :, lnz1) - &
                                                 deltab2_2(1, :, :, lnz2)) * idzh_grid(lnz)
                deltab2_2(4, :, :, unz) =   (3.0*deltab2_2(1, :, :, unz) - &
                                             4.0*deltab2_2(1, :, :, unz1) + &
                                                 deltab2_2(1, :, :, unz2)) * idzh_grid(unz)
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
    !---------------------------------------------------------------------------
    subroutine calc_grad_correl_length(var_flag)
        use mhd_config_module, only: mhd_config
        implicit none
        integer, intent(in) :: var_flag
        real(dp) :: idxh, idyh, idzh
        integer :: unx, uny, unz, lnx, lny, lnz
        integer :: unx1, unx2, uny1, uny2, unz1, unz2
        integer :: lnx1, lnx2, lny1, lny2, lnz1, lnz2
        idxh = 0.5_dp / mhd_config%dx
        idyh = 0.5_dp / mhd_config%dy
        idzh = 0.5_dp / mhd_config%dz
        unx = ubound(lcorr1, 2)
        uny = ubound(lcorr1, 3)
        unz = ubound(lcorr1, 4)
        lnx = lbound(lcorr1, 2)
        lny = lbound(lcorr1, 3)
        lnz = lbound(lcorr1, 4)
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
            lcorr1(2, lnx1:unx1, :, :) = (lcorr1(1, lnx2:unx, :, :) - &
                                          lcorr1(1, lnx:unx2, :, :)) * idxh
            lcorr1(2, lnx, :, :) =  (-3.0*lcorr1(1, lnx, :, :) + &
                                      4.0*lcorr1(1, lnx1, :, :) - &
                                          lcorr1(1, lnx2, :, :)) * idxh
            lcorr1(2, unx, :, :) =   (3.0*lcorr1(1, unx, :, :) - &
                                      4.0*lcorr1(1, unx1, :, :) + &
                                          lcorr1(1, unx2, :, :)) * idxh

            ! d/dy
            if (uny > lny) then
                lcorr1(3, :, lny1:uny1, :) = (lcorr1(1, :, lny2:uny, :) - &
                                              lcorr1(1, :, lny:uny2, :)) * idyh
                lcorr1(3, :, lny, :) =  (-3.0*lcorr1(1, :, lny, :) + &
                                          4.0*lcorr1(1, :, lny1, :) - &
                                              lcorr1(1, :, lny2, :)) * idyh
                lcorr1(3, :, uny, :) =   (3.0*lcorr1(1, :, uny, :) - &
                                          4.0*lcorr1(1, :, uny1, :) + &
                                              lcorr1(1, :, uny2, :)) * idyh
            endif

            ! d/dz
            if (unz > lnz) then
                lcorr1(4, :, :, lnz1:unz1) = (lcorr1(1, :, :, lnz2:unz) - &
                                              lcorr1(1, :, :, lnz:unz2)) * idzh
                lcorr1(4, :, :, lnz) =  (-3.0*lcorr1(1, :, :, lnz) + &
                                          4.0*lcorr1(1, :, :, lnz1) - &
                                              lcorr1(1, :, :, lnz2)) * idzh
                lcorr1(4, :, :, unz) =   (3.0*lcorr1(1, :, :, unz) - &
                                          4.0*lcorr1(1, :, :, unz1) + &
                                              lcorr1(1, :, :, unz2)) * idzh
            endif
        else
            ! d/dx
            lcorr2(2, lnx1:unx1, :, :) = (lcorr2(1, lnx2:unx, :, :) - &
                                          lcorr2(1, lnx:unx2, :, :)) * idxh
            lcorr2(2, lnx, :, :) =  (-3.0*lcorr2(1, lnx, :, :) + &
                                      4.0*lcorr2(1, lnx1, :, :) - &
                                          lcorr2(1, lnx2, :, :)) * idxh
            lcorr2(2, unx, :, :) =   (3.0*lcorr2(1, unx, :, :) - &
                                      4.0*lcorr2(1, unx1, :, :) + &
                                          lcorr2(1, unx2, :, :)) * idxh

            ! d/dy
            if (uny > lny) then
                lcorr2(3, :, lny1:uny1, :) = (lcorr2(1, :, lny2:uny, :) - &
                                              lcorr2(1, :, lny:uny2, :)) * idyh
                lcorr2(3, :, lny, :) =  (-3.0*lcorr2(1, :, lny, :) + &
                                          4.0*lcorr2(1, :, lny1, :) - &
                                              lcorr2(1, :, lny2, :)) * idyh
                lcorr2(3, :, uny, :) =   (3.0*lcorr2(1, :, uny, :) - &
                                          4.0*lcorr2(1, :, uny1, :) + &
                                              lcorr2(1, :, uny2, :)) * idyh
            endif

            ! d/dz
            if (unz > lnz) then
                lcorr2(4, :, :, lnz1:unz1) = (lcorr2(1, :, :, lnz2:unz) - &
                                              lcorr2(1, :, :, lnz:unz2)) * idzh
                lcorr2(4, :, :, lnz) =  (-3.0*lcorr2(1, :, :, lnz) + &
                                          4.0*lcorr2(1, :, :, lnz1) - &
                                              lcorr2(1, :, :, lnz2)) * idzh
                lcorr2(4, :, :, unz) =   (3.0*lcorr2(1, :, :, unz) - &
                                          4.0*lcorr2(1, :, :, unz1) + &
                                              lcorr2(1, :, :, unz2)) * idzh
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
    !---------------------------------------------------------------------------
    subroutine calc_grad_correl_length_nonuniform(var_flag)
        use mhd_config_module, only: mhd_config
        implicit none
        integer, intent(in) :: var_flag
        integer :: unx, uny, unz, lnx, lny, lnz
        integer :: unx1, unx2, uny1, uny2, unz1, unz2
        integer :: lnx1, lnx2, lny1, lny2, lnz1, lnz2
        integer :: ix, iy, iz
        unx = ubound(lcorr1, 2)
        uny = ubound(lcorr1, 3)
        unz = ubound(lcorr1, 4)
        lnx = lbound(lcorr1, 2)
        lny = lbound(lcorr1, 3)
        lnz = lbound(lcorr1, 4)
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
                lcorr1(2, lnx1:unx1, iy, iz) = &
                    (lcorr1(1, lnx2:unx, iy, iz) - &
                     lcorr1(1, lnx:unx2, iy, iz)) * idxh_grid(lnx1:unx1)
            enddo
            enddo
            lcorr1(2, lnx, :, :) =  (-3.0*lcorr1(1, lnx, :, :) + &
                                      4.0*lcorr1(1, lnx1, :, :) - &
                                          lcorr1(1, lnx2, :, :)) * idxh_grid(lnx)
            lcorr1(2, unx, :, :) =   (3.0*lcorr1(1, unx, :, :) - &
                                      4.0*lcorr1(1, unx1, :, :) + &
                                          lcorr1(1, unx2, :, :)) * idxh_grid(unx)

            ! d/dy
            if (uny > lny) then
                do iz = lnz, unz
                do ix = lnx, unx
                    lcorr1(3, ix, lny1:uny1, iz) = &
                        (lcorr1(1, ix, lny2:uny, iz) - &
                         lcorr1(1, ix, lny:uny2, iz)) * idyh_grid(lny1:uny1)
                enddo
                enddo
                lcorr1(3, :, lny, :) =  (-3.0*lcorr1(1, :, lny, :) + &
                                          4.0*lcorr1(1, :, lny1, :) - &
                                              lcorr1(1, :, lny2, :)) * idyh_grid(lny)
                lcorr1(3, :, uny, :) =   (3.0*lcorr1(1, :, uny, :) - &
                                          4.0*lcorr1(1, :, uny1, :) + &
                                              lcorr1(1, :, uny2, :)) * idyh_grid(uny)
            endif

            ! d/dz
            if (unz > lnz) then
                do iy = lny, uny
                do ix = lnx, unx
                    lcorr1(4, ix, iy, lnz1:unz1) = &
                        (lcorr1(1, ix, iy, lnz2:unz) - &
                         lcorr1(1, ix, iy, lnz:unz2)) * idzh_grid(lnz1:unz1)
                enddo
                enddo
                lcorr1(4, :, :, lnz) =  (-3.0*lcorr1(1, :, :, lnz) + &
                                          4.0*lcorr1(1, :, :, lnz1) - &
                                              lcorr1(1, :, :, lnz2)) * idzh_grid(lnz)
                lcorr1(4, :, :, unz) =   (3.0*lcorr1(1, :, :, unz) - &
                                          4.0*lcorr1(1, :, :, unz1) + &
                                              lcorr1(1, :, :, unz2)) * idzh_grid(unz)
            endif
        else
            ! d/dx
            do iz = lnz, unz
            do iy = lny, uny
                lcorr2(2, lnx1:unx1, iy, iz) = &
                    (lcorr2(1, lnx2:unx, iy, iz) - &
                     lcorr2(1, lnx:unx2, iy, iz)) * idxh_grid(lnx1:unx1)
            enddo
            enddo
            lcorr2(2, lnx, :, :) =  (-3.0*lcorr2(1, lnx, :, :) + &
                                      4.0*lcorr2(1, lnx1, :, :) - &
                                          lcorr2(1, lnx2, :, :)) * idxh_grid(lnx)
            lcorr2(2, unx, :, :) =   (3.0*lcorr2(1, unx, :, :) - &
                                      4.0*lcorr2(1, unx1, :, :) + &
                                          lcorr2(1, unx2, :, :)) * idxh_grid(unx)

            ! d/dy
            if (uny > lny) then
                do iz = lnz, unz
                do ix = lnx, unx
                    lcorr2(3, ix, lny1:uny1, iz) = &
                        (lcorr2(1, ix, lny2:uny, iz) - &
                         lcorr2(1, ix, lny:uny2, iz)) * idyh_grid(lny1:uny1)
                enddo
                enddo
                lcorr2(3, :, lny, :) =  (-3.0*lcorr2(1, :, lny, :) + &
                                          4.0*lcorr2(1, :, lny1, :) - &
                                              lcorr2(1, :, lny2, :)) * idyh_grid(lny)
                lcorr2(3, :, uny, :) =   (3.0*lcorr2(1, :, uny, :) - &
                                          4.0*lcorr2(1, :, uny1, :) + &
                                              lcorr2(1, :, uny2, :)) * idyh_grid(uny)
            endif

            ! d/dz
            if (unz > lnz) then
                do iy = lny, uny
                do ix = lnx, unx
                    lcorr2(4, ix, iy, lnz1:unz1) = &
                        (lcorr2(1, ix, iy, lnz2:unz) - &
                         lcorr2(1, ix, iy, lnz:unz2)) * idzh_grid(lnz1:unz1)
                enddo
                enddo
                lcorr2(4, :, :, lnz) =  (-3.0*lcorr2(1, :, :, lnz) + &
                                          4.0*lcorr2(1, :, :, lnz1) - &
                                              lcorr2(1, :, :, lnz2)) * idzh_grid(lnz)
                lcorr2(4, :, :, unz) =   (3.0*lcorr2(1, :, :, unz) - &
                                          4.0*lcorr2(1, :, :, unz1) + &
                                              lcorr2(1, :, :, unz2)) * idzh_grid(unz)
            endif
        endif
        if (mpi_rank == master) then
            write(*, "(A)") "Finished calculating gradients of turbulence correlation length."
        endif
    end subroutine calc_grad_correl_length_nonuniform

    !---------------------------------------------------------------------------
    !< Interpolate the MHD fields and their gradients on one position
    !< Args:
    !<  pos: the lower-left corner of the grid in grid indices.
    !<  weights: for linear interpolation
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  fields(output): fields and their gradients at particle position
    !---------------------------------------------------------------------------
    subroutine interp_fields(pos, weights, rt, fields)
        implicit none
        integer, dimension(3), intent(in) :: pos
        real(dp), dimension(8), intent(in) :: weights
        real(dp), intent(in) :: rt
        real(dp), dimension(nfields+ngrads), intent(out) :: fields
        integer :: ix, iy, iz, i, j, k, ye, ze, index_1d
        real(dp), dimension(nfields+ngrads) :: fields2
        !dir$ attributes align:256 :: fields2

        ix = pos(1)
        iy = pos(2)
        iz = pos(3)
        if (ndim_field > 1) then
            ye = 1
        else
            ye = 0
        endif

        if (ndim_field > 2) then
            ze = 1
        else
            ze = 0
        endif

        fields = 0.0_dp
        fields2 = 0.0_dp
        do k = 0, ze
        do j = 0, ye
        do i = 0, 1
            index_1d = k*4 + j*2 + i + 1
            fields = fields + farray1(:, ix+i, iy+j, iz+k) * weights(index_1d)
            if (time_interp) then
                fields2 = fields2 + farray2(:, ix+i, iy+j, iz+k) * weights(index_1d)
            endif
        enddo
        enddo
        enddo
        !< Time interpolation
        if (time_interp) then
            fields = fields * (1.0 - rt) + fields2 * rt
        endif
    end subroutine interp_fields

    !---------------------------------------------------------------------------
    !< Interpolate the magnetic fluctuation at one position
    !< nz is assumed to be 1.
    !< Args:
    !<  pos: the lower-left corner of the grid in grid indices.
    !<  weights: for linear interpolation
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  db2(output): turbulence variance and its gradients at particle position
    !---------------------------------------------------------------------------
    subroutine interp_magnetic_fluctuation(pos, weights, rt, db2)
        implicit none
        integer, dimension(3), intent(in) :: pos
        real(dp), dimension(8), intent(in) :: weights
        real(dp), intent(in) :: rt
        real(dp), dimension(4), intent(out) :: db2
        integer :: ix, iy, iz, i, j, k, ye, ze, index_1d
        real(dp), dimension(4) :: db2_2
        !dir$ attributes align:32 :: db2_2

        ix = pos(1)
        iy = pos(2)
        iz = pos(3)
        if (ndim_field > 1) then
            ye = 1
        else
            ye = 0
        endif

        if (ndim_field > 2) then
            ze = 1
        else
            ze = 0
        endif

        db2 = 0.0_dp
        db2_2 = 0.0_dp
        do k = 0, ze
        do j = 0, ye
        do i = 0, 1
            index_1d = k*4 + j*2 + i + 1
            db2 = db2 + deltab2_1(:, ix+i, iy+j, iz+k) * weights(index_1d)
            if (time_interp) then
                db2_2 = db2_2 + deltab2_2(:, ix+i, iy+j, iz+k) * weights(index_1d)
            endif
        enddo
        enddo
        enddo
        !< Time interpolation
        if (time_interp) then
            db2 = db2 * (1.0 - rt) + db2_2 * rt
        endif
    end subroutine interp_magnetic_fluctuation

    !---------------------------------------------------------------------------
    !< Interpolate the turbulence correlation length at one position
    !< nz is assumed to be 1.
    !< Args:
    !<  pos: the lower-left corner of the grid in grid indices.
    !<  weights: for linear interpolation
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  lc(output): turbulence correlation length and its graidents at particle position
    !---------------------------------------------------------------------------
    subroutine interp_correlation_length(pos, weights, rt, lc)
        implicit none
        integer, dimension(3), intent(in) :: pos
        real(dp), dimension(8), intent(in) :: weights
        real(dp), intent(in) :: rt
        real(dp), dimension(4), intent(out) :: lc
        integer :: ix, iy, iz, i, j, k, ye, ze, index_1d
        real(dp), dimension(4) :: lc2
        !dir$ attributes align:32 :: lc2

        ix = pos(1)
        iy = pos(2)
        iz = pos(3)
        if (ndim_field > 1) then
            ye = 1
        else
            ye = 0
        endif

        if (ndim_field > 2) then
            ze = 1
        else
            ze = 0
        endif

        lc = 0.0_dp
        lc2 = 0.0_dp
        do k = 0, ze
        do j = 0, ye
        do i = 0, 1
            index_1d = k*4 + j*2 + i + 1
            lc = lc + lcorr1(:, ix+i, iy+j, iz+k) * weights(index_1d)
            if (time_interp) then
                lc2 = lc2 + lcorr2(:, ix+i, iy+j, iz+k) * weights(index_1d)
            endif
        enddo
        enddo
        enddo
        !< Time interpolation
        if (time_interp) then
            lc = lc * (1.0 - rt) + lc2 * rt
        endif
    end subroutine interp_correlation_length

    !<--------------------------------------------------------------------------
    !< Copy fields for usage in the next time interval
    !<--------------------------------------------------------------------------
    subroutine copy_fields
        implicit none
        farray1 = farray2
    end subroutine copy_fields

    !<--------------------------------------------------------------------------
    !< Copy magnetic fluctuation data for usage in the next time interval
    !<--------------------------------------------------------------------------
    subroutine copy_magnetic_fluctuation
        implicit none
        deltab2_1 = deltab2_2
    end subroutine copy_magnetic_fluctuation

    !<--------------------------------------------------------------------------
    !< Copy turbulence correlation length for usage in the next time interval
    !<--------------------------------------------------------------------------
    subroutine copy_correlation_length
        implicit none
        lcorr1 = lcorr2
    end subroutine copy_correlation_length

    !---------------------------------------------------------------------------
    !< Initialize shock x-position indices
    !---------------------------------------------------------------------------
    subroutine init_shock_xpos
        implicit none
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
        if (time_interp) then
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
    !---------------------------------------------------------------------------
    subroutine free_shock_xpos
        implicit none
        deallocate(shock_xpos1)
        if (time_interp) then
            deallocate(shock_xpos2)
        endif
    end subroutine free_shock_xpos

    !---------------------------------------------------------------------------
    !< Locate shock x-position indices. We use Vx to locate the shock here.
    !---------------------------------------------------------------------------
    subroutine locate_shock_xpos
        implicit none
        if (ndim_field == 1) then
            shock_xpos1(1, 1) = maxloc(abs(farray1(nfields+1, :, 1, 1)), dim=1)
        else if (ndim_field == 2) then
            shock_xpos1(:, 1) = maxloc(abs(farray1(nfields+1, :, :, 1)), dim=1)
        else
            shock_xpos1(:, :) = maxloc(abs(farray1(nfields+1, :, :, :)), dim=1)
        endif
        if (time_interp) then
            if (ndim_field == 1) then
                shock_xpos2(1, 1) = maxloc(abs(farray2(nfields+1, :, 1, 1)), dim=1)
            else if (ndim_field == 2) then
                shock_xpos2(:, 1) = maxloc(abs(farray2(nfields+1, :, :, 1)), dim=1)
            else
                shock_xpos2(:, :) = maxloc(abs(farray2(nfields+1, :, :, :)), dim=1)
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
                        abs(farray1(nfields+16, 1:nx, iy, iz) + &
                            (farray1(6, 1:nx, iy, iz) - &
                             farray1(nfields+14, 1:nx, iy, iz)) / xpos_local(1:nx))
                enddo
            enddo
        else
            abs_jz = abs(farray1(nfields+16, 1:nx, 1:ny, 1:nz) - &
                         farray1(nfields+14, 1:nx, 1:ny, 1:nz))
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

    !---------------------------------------------------------------------------
    !< Get the number of cells with current density db2 > db2_min
    !< Args:
    !<  db2_min: the minimum db2
    !<  part_box: box to inject particles
    !---------------------------------------------------------------------------
    function get_ncells_large_db2(db2_min, spherical_coord_flag, part_box) result (ncells_large_db2)
        implicit none
        real(dp), intent(in) :: db2_min
        logical, intent(in) :: spherical_coord_flag
        real(dp), intent(in), dimension(6) :: part_box
        integer :: nx, ny, nz, ix, iy, iz, ncells_large_db2
        logical :: inbox_x, inbox_y, inbox_z, inbox

        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz

        ncells_large_db2 = 0
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
                    if (inbox .and. deltab2_1(1, ix, iy, iz) > db2_min) then
                        ncells_large_db2 = ncells_large_db2 + 1
                    endif
                enddo
            enddo
        enddo
    end function get_ncells_large_db2

    !---------------------------------------------------------------------------
    !< Get the number of cells with compression -divv < divv_min
    !< Args:
    !<  divv_min: the minimum divv
    !<  part_box: box to inject particles
    !---------------------------------------------------------------------------
    function get_ncells_large_divv(divv_min, spherical_coord_flag, part_box) result (ncells_large_divv)
        use constants, only: pi
        implicit none
        real(dp), intent(in) :: divv_min
        logical, intent(in) :: spherical_coord_flag
        real(dp), intent(in), dimension(6) :: part_box
        integer :: nx, ny, nz, ix, iy, iz, ncells_large_divv
        logical :: inbox_x, inbox_y, inbox_z, inbox
        real(dp), allocatable, dimension(:, :, :) :: divv
        real(dp), allocatable, dimension(:) :: ctheta, istheta

        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz

        ncells_large_divv = 0
        allocate(divv(nx, ny, nz))
        allocate(ctheta(ny))
        allocate(istheta(ny))
        if (spherical_coord_flag) then
            ctheta = cos(ypos_local)
            istheta = 1.0 / sin(ypos_local)
            do iz = 1, nz
                do iy = 1, ny
                    divv(:, iy, iz) = farray1(nfields+1, 1:nx, iy, iz) + &
                        2.0 * farray1(1, 1:nx, iy, iz) / xpos_local(1:nx)
                    if (ndim_field > 1) then
                        divv(:, iy, iz) = divv(:, iy, iz) + &
                            (farray1(nfields+5, 1:nx, iy, iz) + &
                             farray1(2, 1:nx, iy, iz)*ctheta(iy)*istheta(iy)) / xpos_local(1:nx)
                    endif
                    if (ndim_field > 1) then
                        divv(:, iy, iz) = divv(:, iy, iz) + &
                            farray1(nfields+9, 1:nx, iy, iz)*istheta(iy) / xpos_local(1:nx)
                    endif
                enddo
            enddo
        else
            divv = farray1(nfields+1, :, :, :)
            if (ndim_field > 1) then
                divv = divv + farray1(nfields+5, :, :, :)
                if (ndim_field > 2) then
                    divv = divv + farray1(nfields+9, :, :, :)
                endif
            endif
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
                    if (inbox .and. -divv(ix, iy, iz) > divv_min) then
                        ncells_large_divv = ncells_large_divv + 1
                    endif
                enddo
            enddo
        enddo
        deallocate(divv)
        deallocate(ctheta)
        deallocate(istheta)
    end function get_ncells_large_divv
end module mhd_data_parallel

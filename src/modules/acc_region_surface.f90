!*******************************************************************************
!< Module for the surface to separate the acceleration region
!< This module is for 3D simulations only.
!*******************************************************************************
module acc_region_surface
    use constants, only: fp, dp
    use simulation_setup_module, only: fconfig
    use mpi_module
    implicit none
    private
    save
    public init_acc_surface, free_acc_surface, read_acc_surface, &
        interp_acc_surface, check_above_acc_surface, copy_acc_surface

    real(dp), allocatable, dimension(:, :) :: acc_surface11, acc_surface12 ! Current,next
    real(dp), allocatable, dimension(:, :) :: acc_surface21, acc_surface22 ! second surface
    character(len=2) :: surface_norm1, surface_norm2
    logical :: surface2_existed, time_interp, is_intersection

    contains
    !---------------------------------------------------------------------------
    !< Initialize the acceleration surface data
    !< Args:
    !<  interp_flag: whether two time steps are needed for interpolation
    !<  surface_norm_flag1: the norm direction of the surface 1
    !<  surface_norm_flag2(optional): the norm direction of the surface 2
    !<  is_intersection_flag(optional): whether it is the intersection of union of the two regions
    !---------------------------------------------------------------------------
    subroutine init_acc_surface(interp_flag, surface_norm_flag1, surface_norm_flag2, &
            is_intersection_flag)
        implicit none
        integer, intent(in) :: interp_flag
        character(*), intent(in) :: surface_norm_flag1
        character(*), intent(in), optional :: surface_norm_flag2
        logical, intent(in), optional :: is_intersection_flag
        integer :: nx, ny, nz

        surface_norm1 = surface_norm_flag1

        nx = fconfig%nx
        ny = fconfig%ny
        nz = fconfig%nz

        if (surface_norm1(2:2) == "x") then
            allocate(acc_surface11(-1:ny+2, -1:nz+2))
        else if (surface_norm1(2:2) == "y") then
            allocate(acc_surface11(-1:nx+2, -1:nz+2))
        else
            allocate(acc_surface11(-1:nx+2, -1:ny+2))
        endif
        acc_surface11 = 0.0_dp
        time_interp = .false.
        ! Next time step
        if (interp_flag == 1) then
            time_interp = .true.
            if (surface_norm1(2:2) == "x") then
                allocate(acc_surface12(-1:ny+2, -1:nz+2))
            else if (surface_norm1(2:2) == "y") then
                allocate(acc_surface12(-1:nx+2, -1:nz+2))
            else
                allocate(acc_surface12(-1:nx+2, -1:ny+2))
            endif
            acc_surface12 = 0.0_dp
        endif

        surface2_existed = .false.
        is_intersection = .false.
        if (present(surface_norm_flag2)) then
            surface_norm2 = surface_norm_flag2
            surface2_existed = .true.
            if (surface_norm2(2:2) == "x") then
                allocate(acc_surface21(-1:ny+2, -1:nz+2))
            else if (surface_norm2(2:2) == "y") then
                allocate(acc_surface21(-1:nx+2, -1:nz+2))
            else
                allocate(acc_surface21(-1:nx+2, -1:ny+2))
            endif
            acc_surface21 = 0.0_dp
            ! Next time step
            if (interp_flag == 1) then
                if (surface_norm2(2:2) == "x") then
                    allocate(acc_surface22(-1:ny+2, -1:nz+2))
                else if (surface_norm2(2:2) == "y") then
                    allocate(acc_surface22(-1:nx+2, -1:nz+2))
                else
                    allocate(acc_surface22(-1:nx+2, -1:ny+2))
                endif
                acc_surface22 = 0.0_dp
            endif
            if (present(is_intersection_flag)) then
                is_intersection = is_intersection_flag
            endif
        endif
    end subroutine init_acc_surface

    !---------------------------------------------------------------------------
    !< Free the acceleration surface data
    !---------------------------------------------------------------------------
    subroutine free_acc_surface
        implicit none
        deallocate(acc_surface11)
        if (time_interp) then
            deallocate(acc_surface12)
        endif
        if (surface2_existed) then
            deallocate(acc_surface21)
            if (time_interp) then
                deallocate(acc_surface22)
            endif
        endif
    end subroutine free_acc_surface

    !---------------------------------------------------------------------------
    !< Read the acceleration region surface data
    !< Args:
    !<  var_flag: indicating which set of variables to save the data. 0 for
    !<            acc_surface*1 and other numbers for acc_surface*2.
    !<  filename1: file name for surface 1
    !<  filename2(optional): file name for surface 2
    !---------------------------------------------------------------------------
    subroutine read_acc_surface(var_flag, filename1, filename2)
        use mpi_io_module, only: set_mpi_datatype_real, set_mpi_info, fileinfo, &
            open_data_mpi_io, read_data_mpi_io
        implicit none
        integer, intent(in) :: var_flag
        character(*), intent(in) :: filename1
        character(*), intent(in), optional :: filename2
        integer :: mpi_datatype, fh
        integer, dimension(2) :: sizes, subsizes, starts
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
        if (surface_norm1(2:2) == "x") then
            sizes(1) = fconfig%nyg
            sizes(2) = fconfig%nzg
            subsizes(1) = fconfig%nyf
            subsizes(2) = fconfig%nzf
            starts(1) = fconfig%iy_min - 1
            starts(2) = fconfig%iz_min - 1
        else if (surface_norm1(2:2) == "y") then
            sizes(1) = fconfig%nxg
            sizes(2) = fconfig%nzg
            subsizes(1) = fconfig%nxf
            subsizes(2) = fconfig%nzf
            starts(1) = fconfig%ix_min - 1
            starts(2) = fconfig%iz_min - 1
        else
            sizes(1) = fconfig%nxg
            sizes(2) = fconfig%nyg
            subsizes(1) = fconfig%nxf
            subsizes(2) = fconfig%nyf
            starts(1) = fconfig%ix_min - 1
            starts(2) = fconfig%iy_min - 1
        endif
        fh = 11
        if (all(sizes == subsizes)) then
            open(unit=fh, file=filename1, access='stream', status='unknown', &
                 form='unformatted', action='read')
            if (var_flag == 0) then
                if (mpi_rank == master) then
                    read(fh, pos=1) acc_surface11
                endif
                call MPI_BCAST(acc_surface11, product(sizes), &
                    MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            else
                if (mpi_rank == master) then
                    read(fh, pos=1) acc_surface12
                endif
                call MPI_BCAST(acc_surface12, product(sizes), &
                    MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
                call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            endif
            close(fh)
            if (present(filename2)) then
                open(unit=fh, file=filename2, access='stream', status='unknown', &
                     form='unformatted', action='read')
                if (var_flag == 0) then
                    if (mpi_rank == master) then
                        read(fh, pos=1) acc_surface21
                    endif
                    call MPI_BCAST(acc_surface21, product(sizes), &
                        MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
                    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                else
                    if (mpi_rank == master) then
                        read(fh, pos=1) acc_surface22
                    endif
                    call MPI_BCAST(acc_surface22, product(sizes), &
                        MPI_DOUBLE_PRECISION, master, MPI_COMM_WORLD, ierr)
                    call MPI_BARRIER(MPI_COMM_WORLD, ierr)
                endif
                close(fh)
            endif
        else
            if (mpi_cross_rank == master) then
                mpi_datatype = set_mpi_datatype_real(sizes, subsizes, starts)
                call open_data_mpi_io(filename1, MPI_MODE_RDONLY, fileinfo, mpi_sub_comm, fh)
                disp = 0
                offset = 0
                if (var_flag == 0) then
                    call read_data_mpi_io(fh, mpi_datatype, subsizes, disp, &
                        offset, acc_surface11)
                else
                    call read_data_mpi_io(fh, mpi_datatype, subsizes, disp, &
                        offset, acc_surface12)
                endif
                call MPI_FILE_CLOSE(fh, ierror)
            endif
            if (var_flag == 0) then
                call MPI_BCAST(acc_surface11, product(subsizes), &
                    MPI_DOUBLE_PRECISION, master, mpi_cross_comm, ierr)
                call MPI_BARRIER(mpi_cross_comm, ierr)
            else
                call MPI_BCAST(acc_surface12, product(subsizes), &
                    MPI_DOUBLE_PRECISION, master, mpi_cross_comm, ierr)
                call MPI_BARRIER(mpi_cross_comm, ierr)
            endif
            if (present(filename2)) then
                if (mpi_cross_rank == master) then
                    call open_data_mpi_io(filename2, MPI_MODE_RDONLY, fileinfo, mpi_sub_comm, fh)
                    if (var_flag == 0) then
                        call read_data_mpi_io(fh, mpi_datatype, subsizes, disp, &
                            offset, acc_surface21)
                    else
                        call read_data_mpi_io(fh, mpi_datatype, subsizes, disp, &
                            offset, acc_surface22)
                    endif
                    call MPI_FILE_CLOSE(fh, ierror)
                endif
                if (var_flag == 0) then
                    call MPI_BCAST(acc_surface21, product(subsizes), &
                        MPI_DOUBLE_PRECISION, master, mpi_cross_comm, ierr)
                    call MPI_BARRIER(mpi_cross_comm, ierr)
                else
                    call MPI_BCAST(acc_surface22, product(subsizes), &
                        MPI_DOUBLE_PRECISION, master, mpi_cross_comm, ierr)
                    call MPI_BARRIER(mpi_cross_comm, ierr)
                endif
            endif
        endif
        if (mpi_rank == master) then
            write(*, "(A)") "Finished reading the acceleration region surface."
        endif
    end subroutine read_acc_surface

    !---------------------------------------------------------------------------
    !< Interpolate the surface data on one position
    !< Args:
    !<  pos: the lower-left corner of the grid in grid indices.
    !<  weights: for linear interpolation
    !<  rt: the offset to the earlier time point of the MHD data. It is
    !<      normalized to the time interval of the MHD data output.
    !<  surface_height1: the height on surface 1
    !<  surface_height2: the height on surface 2
    !---------------------------------------------------------------------------
    subroutine interp_acc_surface(pos, weights, rt, surface_height1, surface_height2)
        implicit none
        integer, dimension(3), intent(in) :: pos
        real(dp), dimension(8), intent(in) :: weights
        real(dp), intent(in) :: rt
        real(dp), intent(out) :: surface_height1, surface_height2
        real(dp) :: surface_height11, surface_height12
        real(dp) :: surface_height21, surface_height22
        real(dp), dimension(2, 2) :: weights_2d
        integer :: i1, j1

        if (surface_norm1(2:2) == "x") then
            i1 = pos(2)
            j1 = pos(3)
            weights_2d(1, 1) = weights(1) + weights(2)
            weights_2d(2, 1) = weights(3) + weights(4)
            weights_2d(1, 2) = weights(5) + weights(6)
            weights_2d(2, 2) = weights(7) + weights(8)
        else if (surface_norm1(2:2) == "y") then
            i1 = pos(1)
            j1 = pos(3)
            weights_2d(1, 1) = weights(1) + weights(3)
            weights_2d(2, 1) = weights(2) + weights(4)
            weights_2d(1, 2) = weights(5) + weights(7)
            weights_2d(2, 2) = weights(6) + weights(8)
        else
            i1 = pos(1)
            j1 = pos(2)
            weights_2d(1, 1) = weights(1) + weights(5)
            weights_2d(2, 1) = weights(2) + weights(6)
            weights_2d(1, 2) = weights(3) + weights(7)
            weights_2d(2, 2) = weights(4) + weights(8)
        endif

        surface_height11 = sum(acc_surface11(i1:i1+1, j1:j1+1) * weights_2d)
        if (time_interp) then
            surface_height12 = sum(acc_surface12(i1:i1+1, j1:j1+1) * weights_2d)
            surface_height1 = surface_height11 * (1.0 - rt) + surface_height12 * rt
        else
            surface_height1 = surface_height11
        endif

        ! if there is a second surface
        if (surface2_existed) then
            if (surface_norm2(2:2) /= surface_norm1(2:2)) then
                ! the two surfaces are not along the same direction
                if (surface_norm2(2:2) == "x") then
                    i1 = pos(2)
                    j1 = pos(3)
                    weights_2d(1, 1) = weights(1) + weights(2)
                    weights_2d(2, 1) = weights(3) + weights(4)
                    weights_2d(1, 2) = weights(5) + weights(6)
                    weights_2d(2, 2) = weights(7) + weights(8)
                else if (surface_norm2(2:2) == "y") then
                    i1 = pos(1)
                    j1 = pos(3)
                    weights_2d(1, 1) = weights(1) + weights(3)
                    weights_2d(2, 1) = weights(2) + weights(4)
                    weights_2d(1, 2) = weights(5) + weights(7)
                    weights_2d(2, 2) = weights(6) + weights(8)
                else
                    i1 = pos(1)
                    j1 = pos(2)
                    weights_2d(1, 1) = weights(1) + weights(5)
                    weights_2d(2, 1) = weights(2) + weights(6)
                    weights_2d(1, 2) = weights(3) + weights(7)
                    weights_2d(2, 2) = weights(4) + weights(8)
                endif
            endif
            surface_height21 = sum(acc_surface21(i1:i1+1, j1:j1+1) * weights_2d)
            if (time_interp) then
                surface_height22 = sum(acc_surface22(i1:i1+1, j1:j1+1) * weights_2d)
                surface_height2 = surface_height21 * (1.0 - rt) + surface_height22 * rt
            else
                surface_height2 = surface_height21
            endif
        else
            surface_height2 = 0.0_dp
        endif
    end subroutine interp_acc_surface

    !---------------------------------------------------------------------------
    !< Check whether the particle is above the acceleration region surface
    !< Args:
    !<  x, y, z: location of the particle
    !<  surface_height1, surface_height2: the heights on the two surfaces
    !---------------------------------------------------------------------------
    function check_above_acc_surface(x, y, z, surface_height1, surface_height2) result (in_acc_region)
        implicit none
        real(dp), intent(in) :: x, y, z, surface_height1, surface_height2
        real(dp) :: ptl_height
        logical :: in_acc_region

        in_acc_region = .true.

        if (surface_norm1(2:2) == "x") then
            ptl_height = x
        else if (surface_norm1(2:2) == "y") then
            ptl_height = y
        else
            ptl_height = z
        endif
        if (surface_norm1(1:1) == "+") then
            in_acc_region = ptl_height > surface_height1
        else
            in_acc_region = ptl_height < surface_height1
        endif

        if (surface2_existed) then
            if (surface_norm2(2:2) == "x") then
                ptl_height = x
            else if (surface_norm2(2:2) == "y") then
                ptl_height = y
            else
                ptl_height = z
            endif
            if (is_intersection) then ! interaction of the two
                if (surface_norm2(1:1) == "+") then
                    in_acc_region = in_acc_region .and. (ptl_height > surface_height2)
                else
                    in_acc_region = in_acc_region .and. (ptl_height < surface_height2)
                endif
            else ! union of the two
                if (surface_norm2(1:1) == "+") then
                    in_acc_region = in_acc_region .or. (ptl_height > surface_height2)
                else
                    in_acc_region = in_acc_region .or. (ptl_height < surface_height2)
                endif
            endif
        endif
    end function check_above_acc_surface

    !<--------------------------------------------------------------------------
    !< Copy surface data for usage in the next time interval
    !<--------------------------------------------------------------------------
    subroutine copy_acc_surface
        implicit none
        acc_surface11 = acc_surface12
        if (surface2_existed) then
            acc_surface21 = acc_surface22
        endif
    end subroutine copy_acc_surface
end module acc_region_surface

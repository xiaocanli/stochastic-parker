!*******************************************************************************
! This contains one subroutine to open field data file using MPI I/O, one
! subroutine to read data using MPI I/O.
!*******************************************************************************
module mpi_io_module
    use mpi_module
    implicit none
    private
    public set_mpi_datatype
    public open_data_mpi_io, read_data_mpi_io, write_data_mpi_io
    public fileinfo, set_mpi_info

    interface read_data_mpi_io
        module procedure &
            read_data_mpi_io_real_4d, read_data_mpi_io_real_3d, &
            read_data_mpi_io_real_2d, read_data_mpi_io_real_1d
    end interface read_data_mpi_io

    interface write_data_mpi_io
        module procedure &
            write_data_mpi_io_real_3d, write_data_mpi_io_real_2d, &
            write_data_mpi_io_real_1d
    end interface write_data_mpi_io

    integer :: fileinfo     ! MPI_INFO object

    contains

    !---------------------------------------------------------------------------
    ! Create a MPI data type and commit it.    
    !---------------------------------------------------------------------------
    function set_mpi_datatype(sizes, subsizes, starts) result(datatype)
        implicit none
        integer, dimension(:), intent(in) :: sizes, subsizes, starts
        integer :: datatype
        integer :: sz

        sz = size(sizes)

        call MPI_TYPE_CREATE_SUBARRAY(sz, sizes, subsizes, starts, &
                MPI_ORDER_FORTRAN, MPI_REAL, datatype, ierror)
        call MPI_TYPE_COMMIT(datatype, ierror)
    end function set_mpi_datatype

    !---------------------------------------------------------------------------
    ! Create a MPI_INFO object and have proper settings for ROMIO's data-sieving
    ! and collective buffering.
    !
    ! Outputs:
    !   fileinfo: the MPI_INFO.
    !---------------------------------------------------------------------------
    subroutine set_mpi_info
        implicit none
        call MPI_INFO_CREATE(fileinfo, ierror)
        !! Disable ROMIO's data-sieving
        call MPI_INFO_SET(fileinfo, "romio_ds_read", "automatic", ierror)
        call MPI_INFO_SET(fileinfo, "romio_ds_write", "automatic", ierror)
        !! Enable ROMIO's collective buffering
        call MPI_INFO_SET(fileinfo, "romio_cb_read", "enable", ierror)
        call MPI_INFO_SET(fileinfo, "romio_cb_write", "enable", ierror)

        ! ! For panfs parallel file system.
        ! call MPI_INFO_SET(fileinfo, "panfs_concurrent_write", "1", ierror)
    end subroutine set_mpi_info

    !---------------------------------------------------------------------------
    ! Open one data file using MPI/IO.
    ! Input:
    !   fname: file name.
    !   amode: file access mode.
    !   fileinfo: MPI_INFO
    ! Output:
    !   fh: file handler.
    !---------------------------------------------------------------------------
    subroutine open_data_mpi_io(fname, amode, fileinfo, fh)
        implicit none
        character(*), intent(in) :: fname
        integer, intent(in) :: amode, fileinfo
        integer, intent(out) :: fh
        call MPI_FILE_OPEN(MPI_COMM_WORLD, fname, amode, &
            fileinfo, fh, ierror)
        if (ierror /= 0) then
            call MPI_ERROR_STRING(ierror, err_msg, err_length, ierror2)
            print*, "Error in MPI_FILE_OPEN: ", trim(err_msg)
        endif
    end subroutine open_data_mpi_io

    !---------------------------------------------------------------------------
    ! Set view for the file to read.
    ! Inputs:
    !   fh: file handler.
    !   datatype: MPI data type.
    !   disp: displacement form the beginning of the file (in bytes).
    !---------------------------------------------------------------------------
    subroutine set_file_view_real(fh, datatype, disp)
        implicit none
        integer, intent(in) :: fh, datatype
        integer(kind=MPI_OFFSET_KIND), intent(in) :: disp
        call MPI_FILE_SET_VIEW(fh, disp, MPI_REAL, datatype, 'native', &
            MPI_INFO_NULL, ierror)
        if (ierror /= 0) then
            call MPI_ERROR_STRING(ierror, err_msg, err_length, ierror2)
            print*, "Error in MPI_FILE_SET_VIEW: ", trim(err_msg)
        endif
    end subroutine set_file_view_real

    !---------------------------------------------------------------------------
    ! Handle MPI_FILE_READ error.
    !---------------------------------------------------------------------------
    subroutine handle_read_error
        implicit none
        if (ierror /= 0) then
            call MPI_ERROR_STRING(ierror, err_msg, err_length, ierror2)
            print*, "Error in MPI_FILE_READ: ", trim(err_msg)
        endif
    end subroutine handle_read_error

    !---------------------------------------------------------------------------
    ! Handle MPI_FILE_WRITE error.
    !---------------------------------------------------------------------------
    subroutine handle_write_error
        implicit none
        if (ierror /= 0) then
            call MPI_ERROR_STRING(ierror, err_msg, err_length, ierror2)
            print*, "Error in MPI_FILE_WRITE: ", trim(err_msg)
        endif
    end subroutine handle_write_error

    !---------------------------------------------------------------------------
    ! Read data from files using MPI/IO for 4D REAL data.
    ! Inputs:
    !   fh: file handler.
    !   datatype: MPI data type.
    !   subsizes: the sub-sizes of the data in current MPI process.
    !   disp: displacement form the beginning of the file (in bytes).
    !   offset: offset from current file view (in data etypes (e.g. int, real)).
    ! Output:
    !   rdata: the data read from the file.
    !---------------------------------------------------------------------------
    subroutine read_data_mpi_io_real_4d(fh, datatype, subsizes, disp, offset, rdata)
        use constants, only: fp
        implicit none
        integer, intent(in) :: fh, datatype
        integer, dimension(4), intent(in) :: subsizes
        integer(kind=MPI_OFFSET_KIND), intent(in) :: disp, offset
        real(fp), dimension(:, :, :, :), intent(out) :: rdata

        call set_file_view_real(fh, datatype, disp)
        call MPI_FILE_READ_AT_ALL(fh, offset, rdata, &
                product(subsizes), MPI_REAL, status, ierror)
        call handle_read_error
    end subroutine read_data_mpi_io_real_4d

    !---------------------------------------------------------------------------
    ! Read data from files using MPI/IO for 3D REAL data.
    ! Inputs:
    !   fh: file handler.
    !   datatype: MPI data type.
    !   subsizes: the sub-sizes of the data in current MPI process.
    !   disp: displacement form the beginning of the file (in bytes).
    !   offset: offset from current file view (in data etypes (e.g. int, real)).
    ! Output:
    !   rdata: the data read from the file.
    !---------------------------------------------------------------------------
    subroutine read_data_mpi_io_real_3d(fh, datatype, subsizes, disp, offset, rdata)
        use constants, only: fp
        implicit none
        integer, intent(in) :: fh, datatype
        integer, dimension(3), intent(in) :: subsizes
        integer(kind=MPI_OFFSET_KIND), intent(in) :: disp, offset
        real(fp), dimension(:, :, :), intent(out) :: rdata

        call set_file_view_real(fh, datatype, disp)
        call MPI_FILE_READ_AT_ALL(fh, offset, rdata, &
                product(subsizes), MPI_REAL, status, ierror)
        call handle_read_error
    end subroutine read_data_mpi_io_real_3d

    !---------------------------------------------------------------------------
    ! Read data from files using MPI/IO for 2D REAL data.
    !---------------------------------------------------------------------------
    subroutine read_data_mpi_io_real_2d(fh, datatype, subsizes, disp, offset, rdata)
        use constants, only: fp
        implicit none
        integer, intent(in) :: fh, datatype
        integer, dimension(2), intent(in) :: subsizes
        integer(kind=MPI_OFFSET_KIND), intent(in) :: disp, offset
        real(fp), dimension(:, :), intent(out) :: rdata

        call set_file_view_real(fh, datatype, disp)
        call MPI_FILE_READ_AT_ALL(fh, offset, rdata, &
                product(subsizes), MPI_REAL, status, ierror)
        call handle_read_error
    end subroutine read_data_mpi_io_real_2d

    !---------------------------------------------------------------------------
    ! Read data from files using MPI/IO for 1D REAL data.
    !---------------------------------------------------------------------------
    subroutine read_data_mpi_io_real_1d(fh, datatype, subsizes, disp, offset, rdata)
        use constants, only: fp
        implicit none
        integer, intent(in) :: fh, datatype
        integer, dimension(1), intent(in) :: subsizes
        integer(kind=MPI_OFFSET_KIND), intent(in) :: disp, offset
        real(fp), dimension(:), intent(out) :: rdata

        call set_file_view_real(fh, datatype, disp)
        call MPI_FILE_READ_AT_ALL(fh, offset, rdata, &
                subsizes(1), MPI_REAL, status, ierror)
        call handle_read_error
    end subroutine read_data_mpi_io_real_1d

    !---------------------------------------------------------------------------
    ! Write data to files using MPI/IO for 3D REAL data.
    ! Input:
    !   fh: file handler.
    !   datatype: MPI data type.
    !   subsizes: the sub-sizes of the data in current MPI process.
    !   disp: displacement form the beginning of the file (in bytes).
    !   offset: offset from current file view (in data etypes (e.g. int, real)).
    ! Output:
    !   wdata: the data to write to file.
    !---------------------------------------------------------------------------
    subroutine write_data_mpi_io_real_3d(fh, datatype, subsizes, disp, offset, wdata)
        use constants, only: fp
        implicit none
        integer, intent(in) :: fh, datatype
        integer, dimension(3), intent(in) :: subsizes
        integer(kind=MPI_OFFSET_KIND), intent(in) :: disp, offset
        real(fp), dimension(:,:,:), intent(in) :: wdata

        call set_file_view_real(fh, datatype, disp)

        call MPI_FILE_WRITE_AT_ALL(fh, offset, wdata, &
                product(subsizes), MPI_REAL, status, ierror)
        call handle_write_error
    end subroutine write_data_mpi_io_real_3d

    !---------------------------------------------------------------------------
    ! Write data to files using MPI/IO for 2D REAL data.
    !---------------------------------------------------------------------------
    subroutine write_data_mpi_io_real_2d(fh, datatype, subsizes, disp, offset, wdata)
        use constants, only: fp
        implicit none
        integer, intent(in) :: fh, datatype
        integer, dimension(2), intent(in) :: subsizes
        integer(kind=MPI_OFFSET_KIND), intent(in) :: disp, offset
        real(fp), dimension(:,:), intent(in) :: wdata

        call set_file_view_real(fh, datatype, disp)

        call MPI_FILE_WRITE_AT_ALL(fh, offset, wdata, &
                product(subsizes), MPI_REAL, status, ierror)
        call handle_write_error
    end subroutine write_data_mpi_io_real_2d

    !---------------------------------------------------------------------------
    ! Write data to files using MPI/IO for 1D REAL data.
    !---------------------------------------------------------------------------
    subroutine write_data_mpi_io_real_1d(fh, datatype, subsizes, disp, offset, wdata)
        use constants, only: fp
        implicit none
        integer, intent(in) :: fh, datatype
        integer, dimension(1), intent(in) :: subsizes
        integer(kind=MPI_OFFSET_KIND), intent(in) :: disp, offset
        real(fp), dimension(:), intent(in) :: wdata

        call set_file_view_real(fh, datatype, disp)

        call MPI_FILE_WRITE_AT_ALL(fh, offset, wdata, &
                subsizes(1), MPI_REAL, status, ierror)
        call handle_write_error
    end subroutine write_data_mpi_io_real_1d

end module mpi_io_module

!<******************************************************************************
!< This contains subroutines to create, open, close, read, and write HDF5 files
!<******************************************************************************
module hdf5_io
    use constants, only: fp, dp
    use hdf5
    implicit none
    private
    public create_file_h5, open_file_h5, close_file_h5, read_data_h5, write_data_h5

    interface read_data_h5
        module procedure &
            read_integer_h5_1d, read_integer_h5_2d, read_integer_h5_3d, &
            read_integer_h5_4d, read_real_h5_1d, read_real_h5_2d, &
            read_real_h5_3d, read_real_h5_4d, read_double_h5_1d, &
            read_double_h5_2d, read_double_h5_3d, read_double_h5_4d
    end interface read_data_h5

    interface write_data_h5
        module procedure &
            write_integer_h5_1d, write_integer_h5_2d, write_integer_h5_3d, &
            write_integer_h5_4d, write_real_h5_1d, write_real_h5_2d, &
            write_real_h5_3d, write_real_h5_4d, write_double_h5_1d, &
            write_double_h5_2d, write_double_h5_3d, write_double_h5_4d
    end interface write_data_h5

    contains

    !<--------------------------------------------------------------------------
    !< Create one data file using serial or parallel HDF5
    !< Input:
    !<   filename: file name
    !<   access_flag: file access flag
    !<   parallel_hdf5(optional): whether to use parallel HDF5
    !<   mpi_comm(optional): MPI communicator
    !< Output:
    !<   file_id: file handler
    !<--------------------------------------------------------------------------
    subroutine create_file_h5(filename, access_flag, file_id, parallel_hdf5, mpi_comm)
        implicit none
        character(*), intent(in) :: filename
        logical, intent(in), optional :: parallel_hdf5
        integer, intent(in), optional :: mpi_comm
        integer, intent(in) :: access_flag
        integer(hid_t), intent(out) :: file_id
        logical :: use_parallel_hdf5
        integer :: fileinfo, ierror
        integer(hid_t) :: plist_id
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif

        if (use_parallel_hdf5) then
            call MPI_INFO_CREATE(fileinfo, ierror)
            call MPI_INFO_SET(fileinfo, "romio_ds_read", "disable", ierror)
            call MPI_INFO_SET(fileinfo, "romio_cb_read", "enable", ierror)
            call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierror)
            call h5pset_fapl_mpio_f(plist_id, mpi_comm, fileinfo, ierror)
            call MPI_INFO_FREE(fileinfo, ierror)
            call h5fcreate_f(trim(filename), access_flag, file_id, ierror, access_prp=plist_id)
            call h5pclose_f(plist_id, ierror)
        else
            call h5fcreate_f(trim(filename), access_flag, file_id, ierror, &
                access_prp=h5p_default_f)
        endif
    end subroutine create_file_h5

    !<--------------------------------------------------------------------------
    !< Open one data file using serial or parallel HDF5
    !< Input:
    !<   filename: file name
    !<   access_flag: file access flag
    !<   parallel_hdf5(optional): whether to use parallel HDF5
    !<   mpi_comm(optional): MPI communicator
    !< Output:
    !<   file_id: file handler
    !<--------------------------------------------------------------------------
    subroutine open_file_h5(filename, access_flag, file_id, parallel_hdf5, mpi_comm)
        implicit none
        character(*), intent(in) :: filename
        logical, intent(in), optional :: parallel_hdf5
        integer, intent(in), optional :: mpi_comm
        integer, intent(in) :: access_flag
        integer(hid_t), intent(out) :: file_id
        logical :: use_parallel_hdf5
        integer :: fileinfo, ierror
        integer(hid_t) :: plist_id
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif

        if (use_parallel_hdf5) then
            call MPI_INFO_CREATE(fileinfo, ierror)
            call MPI_INFO_SET(fileinfo, "romio_ds_read", "disable", ierror)
            call MPI_INFO_SET(fileinfo, "romio_cb_read", "enable", ierror)
            call h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, ierror)
            call h5pset_fapl_mpio_f(plist_id, mpi_comm, fileinfo, ierror)
            call MPI_INFO_FREE(fileinfo, ierror)
            call h5fopen_f(trim(filename), access_flag, file_id, ierror, access_prp=plist_id)
            call h5pclose_f(plist_id, ierror)
        else
            call h5fopen_f(trim(filename), access_flag, file_id, ierror, &
                access_prp=h5p_default_f)
        endif
    end subroutine open_file_h5

    !<--------------------------------------------------------------------------
    !< Close a HDF5 file
    !<--------------------------------------------------------------------------
    subroutine close_file_h5(file_id)
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer :: ierror
        call h5fclose_f(file_id, ierror)
    end subroutine close_file_h5

    !<--------------------------------------------------------------------------
    !< Initial setup for reading HDF5 file in serial
    !< Input:
    !<   dset_id: dataset ID
    !<   dcount: 1D array specifying the size of each dimension of the dataset
    !<   doffset: offset of start of hyperslab
    !< Output:
    !<   filespace: dataspace identifier for the whole dataset
    !<   memspace: dataspace identifier for current hyperslab
    !<--------------------------------------------------------------------------
    subroutine init_read_serial_h5(dset_id, dcount, doffset, filespace, memspace)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(:), intent(in) :: dcount, doffset
        integer(hid_t), intent(out) :: filespace, memspace
        integer :: rank, error
        rank = size(dcount)
        call h5screate_simple_f(rank, dcount, memspace, error)
        call h5dget_space_f(dset_id, filespace, error)
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, &
            dcount, error)
    end subroutine init_read_serial_h5

    !<--------------------------------------------------------------------------
    !< Initial setup for writing HDF5 file in serial
    !< Input:
    !<   dset_id: dataset ID
    !<   dcount: 1D array specifying the size of each dimension of the dataset
    !<   doffset: offset of start of hyperslab
    !<   dset_dims: offset of start of hyperslab
    !< Output:
    !<   filespace: dataspace identifier for the whole dataset
    !<   memspace: dataspace identifier for current hyperslab
    !<--------------------------------------------------------------------------
    subroutine init_write_serial_h5(dset_id, dcount, doffset, dset_dims, &
            filespace, memspace)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(:), intent(in) :: dcount, doffset, dset_dims
        integer(hid_t), intent(out) :: filespace, memspace
        integer :: rank, error
        rank = size(dcount)
        call h5screate_simple_f(rank, dset_dims, filespace, error)
        call h5screate_simple_f(rank, dcount, memspace, error)
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, &
            dcount, error)
    end subroutine init_write_serial_h5

    !<--------------------------------------------------------------------------
    !< Initial setup for reading HDF5 file in parallel
    !< Input:
    !<   dset_id: dataset ID
    !<   dcount: 1D array specifying the size of each dimension of the dataset
    !<   doffset: offset of start of hyperslab
    !<   collective_io(optional): whether to use collective IO
    !< Output:
    !<   filespace: dataspace identifier for the whole dataset
    !<   memspace: dataspace identifier for current hyperslab
    !<   plist_id: property list for dataset read/write
    !<--------------------------------------------------------------------------
    subroutine init_read_parallel_h5(dset_id, dcount, doffset, filespace, &
            memspace, plist_id, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(:), intent(in) :: dcount, doffset
        integer(hid_t), intent(out) :: filespace, memspace, plist_id
        logical, intent(in), optional :: collective_io
        logical :: use_collective_io
        integer :: rank, error

        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        rank = size(dcount)
        ! Create property list for collective dataset write
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        if (use_collective_io) then
            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        else
            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
        endif

        call h5screate_simple_f(rank, dcount, memspace, error)
        call h5dget_space_f(dset_id, filespace, error)
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, &
            dcount, error)
    end subroutine init_read_parallel_h5

    !<--------------------------------------------------------------------------
    !< Initial setup for writing HDF5 file in parallel
    !< Input:
    !<   dset_id: dataset ID
    !<   dcount: 1D array specifying the size of each dimension of the dataset
    !<   doffset: offset of start of hyperslab
    !<   dset_dims: offset of start of hyperslab
    !<   collective_io(optional): whether to use collective IO
    !< Output:
    !<   filespace: dataspace identifier for the whole dataset
    !<   memspace: dataspace identifier for current hyperslab
    !<   plist_id: property list for dataset read/write
    !<--------------------------------------------------------------------------
    subroutine init_write_parallel_h5(dset_id, dcount, doffset, dset_dims, &
            filespace, memspace, plist_id, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(:), intent(in) :: dcount, doffset, dset_dims
        integer(hid_t), intent(out) :: filespace, memspace, plist_id
        logical, intent(in), optional :: collective_io
        logical :: use_collective_io
        integer :: rank, error

        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        rank = size(dcount)
        ! Create property list for collective dataset write
        call h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        if (use_collective_io) then
            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
        else
            call h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
        endif

        call h5screate_simple_f(rank, dset_dims, filespace, error)
        call h5screate_simple_f(rank, dcount, memspace, error)
        call h5sselect_hyperslab_f(filespace, H5S_SELECT_SET_F, doffset, &
            dcount, error)
    end subroutine init_write_parallel_h5

    !<--------------------------------------------------------------------------
    !< Finalize reading or writing a HDF5 file
    !<--------------------------------------------------------------------------
    subroutine final_rw_h5(filespace, memspace, plist_id)
        implicit none
        integer(hid_t), intent(in) :: filespace, memspace
        integer(hid_t), intent(in), optional :: plist_id
        integer :: error
        call h5sclose_f(filespace, error)
        call h5sclose_f(memspace, error)
        if (present(plist_id)) then
            call h5pclose_f(plist_id, error)
        endif
    end subroutine final_rw_h5

    !<--------------------------------------------------------------------------
    !< Read HDF5 dataset for 1D integer data
    !<--------------------------------------------------------------------------
    subroutine read_integer_h5_1d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        integer, dimension(:), intent(out) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_read_parallel_h5(dset_id, dcount, doffset, filespace, &
                memspace, plist_id, use_collective_io)
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_read_serial_h5(dset_id, dcount, doffset, filespace, memspace)
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine read_integer_h5_1d

    !<--------------------------------------------------------------------------
    !< Read HDF5 dataset for 2D integer data
    !<--------------------------------------------------------------------------
    subroutine read_integer_h5_2d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(2), intent(in) :: dcount, doffset, dset_dims
        integer, dimension(:, :), intent(out) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_read_parallel_h5(dset_id, dcount, doffset, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_read_serial_h5(dset_id, dcount, doffset, filespace, memspace)
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine read_integer_h5_2d

    !<--------------------------------------------------------------------------
    !< Read HDF5 dataset for 3D integer data
    !<--------------------------------------------------------------------------
    subroutine read_integer_h5_3d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(3), intent(in) :: dcount, doffset, dset_dims
        integer, dimension(:, :, :), intent(out) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_read_parallel_h5(dset_id, dcount, doffset, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_read_serial_h5(dset_id, dcount, doffset, filespace, memspace)
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine read_integer_h5_3d

    !<--------------------------------------------------------------------------
    !< Read HDF5 dataset for 4D integer data
    !<--------------------------------------------------------------------------
    subroutine read_integer_h5_4d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(4), intent(in) :: dcount, doffset, dset_dims
        integer, dimension(:, :, :, :), intent(out) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_read_parallel_h5(dset_id, dcount, doffset, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_read_serial_h5(dset_id, dcount, doffset, filespace, memspace)
            call h5dread_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine read_integer_h5_4d

    !<--------------------------------------------------------------------------
    !< Read HDF5 dataset for 1D real data
    !<--------------------------------------------------------------------------
    subroutine read_real_h5_1d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        real(fp), dimension(:), intent(out) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_read_parallel_h5(dset_id, dcount, doffset, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dread_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_read_serial_h5(dset_id, dcount, doffset, filespace, memspace)
            call h5dread_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine read_real_h5_1d

    !<--------------------------------------------------------------------------
    !< Read HDF5 dataset for 2D real data
    !<--------------------------------------------------------------------------
    subroutine read_real_h5_2d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(2), intent(in) :: dcount, doffset, dset_dims
        real(fp), dimension(:, :), intent(out) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_read_parallel_h5(dset_id, dcount, doffset, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dread_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_read_serial_h5(dset_id, dcount, doffset, filespace, memspace)
            call h5dread_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine read_real_h5_2d

    !<--------------------------------------------------------------------------
    !< Read HDF5 dataset for 3D real data
    !<--------------------------------------------------------------------------
    subroutine read_real_h5_3d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(3), intent(in) :: dcount, doffset, dset_dims
        real(fp), dimension(:, :, :), intent(out) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_read_parallel_h5(dset_id, dcount, doffset, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dread_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_read_serial_h5(dset_id, dcount, doffset, filespace, memspace)
            call h5dread_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine read_real_h5_3d

    !<--------------------------------------------------------------------------
    !< Read HDF5 dataset for 4D real data
    !<--------------------------------------------------------------------------
    subroutine read_real_h5_4d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(4), intent(in) :: dcount, doffset, dset_dims
        real(fp), dimension(:, :, :, :), intent(out) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_read_parallel_h5(dset_id, dcount, doffset, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dread_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_read_serial_h5(dset_id, dcount, doffset, filespace, memspace)
            call h5dread_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine read_real_h5_4d

    !<--------------------------------------------------------------------------
    !< Read HDF5 dataset for 1D double data
    !<--------------------------------------------------------------------------
    subroutine read_double_h5_1d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        real(dp), dimension(:), intent(out) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_read_parallel_h5(dset_id, dcount, doffset, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_read_serial_h5(dset_id, dcount, doffset, filespace, memspace)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine read_double_h5_1d

    !<--------------------------------------------------------------------------
    !< Read HDF5 dataset for 2D double data
    !<--------------------------------------------------------------------------
    subroutine read_double_h5_2d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(2), intent(in) :: dcount, doffset, dset_dims
        real(dp), dimension(:, :), intent(out) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_read_parallel_h5(dset_id, dcount, doffset, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_read_serial_h5(dset_id, dcount, doffset, filespace, memspace)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine read_double_h5_2d

    !<--------------------------------------------------------------------------
    !< Read HDF5 dataset for 3D double data
    !<--------------------------------------------------------------------------
    subroutine read_double_h5_3d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(3), intent(in) :: dcount, doffset, dset_dims
        real(dp), dimension(:, :, :), intent(out) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_read_parallel_h5(dset_id, dcount, doffset, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_read_serial_h5(dset_id, dcount, doffset, filespace, memspace)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine read_double_h5_3d

    !<--------------------------------------------------------------------------
    !< Read HDF5 dataset for 4D double data
    !<--------------------------------------------------------------------------
    subroutine read_double_h5_4d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(4), intent(in) :: dcount, doffset, dset_dims
        real(dp), dimension(:, :, :, :), intent(out) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_read_parallel_h5(dset_id, dcount, doffset, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_read_serial_h5(dset_id, dcount, doffset, filespace, memspace)
            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine read_double_h5_4d

    !<--------------------------------------------------------------------------
    !< Write HDF5 dataset for 1D integer data
    !<--------------------------------------------------------------------------
    subroutine write_integer_h5_1d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        integer, dimension(:), intent(in) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_write_parallel_h5(dset_id, dcount, doffset, dset_dims, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_write_serial_h5(dset_id, dcount, doffset, dset_dims, filespace, memspace)
            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine write_integer_h5_1d

    !<--------------------------------------------------------------------------
    !< Write HDF5 dataset for 2D integer data
    !<--------------------------------------------------------------------------
    subroutine write_integer_h5_2d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(2), intent(in) :: dcount, doffset, dset_dims
        integer, dimension(:, :), intent(in) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_write_parallel_h5(dset_id, dcount, doffset, dset_dims, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_write_serial_h5(dset_id, dcount, doffset, dset_dims, filespace, memspace)
            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine write_integer_h5_2d

    !<--------------------------------------------------------------------------
    !< Write HDF5 dataset for 3D integer data
    !<--------------------------------------------------------------------------
    subroutine write_integer_h5_3d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(3), intent(in) :: dcount, doffset, dset_dims
        integer, dimension(:, :, :), intent(in) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_write_parallel_h5(dset_id, dcount, doffset, dset_dims, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_write_serial_h5(dset_id, dcount, doffset, dset_dims, filespace, memspace)
            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine write_integer_h5_3d

    !<--------------------------------------------------------------------------
    !< Write HDF5 dataset for 4D integer data
    !<--------------------------------------------------------------------------
    subroutine write_integer_h5_4d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(4), intent(in) :: dcount, doffset, dset_dims
        integer, dimension(:, :, :, :), intent(in) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_write_parallel_h5(dset_id, dcount, doffset, dset_dims, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_write_serial_h5(dset_id, dcount, doffset, dset_dims, filespace, memspace)
            call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine write_integer_h5_4d

    !<--------------------------------------------------------------------------
    !< Write HDF5 dataset for 1D real data
    !<--------------------------------------------------------------------------
    subroutine write_real_h5_1d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        real(fp), dimension(:), intent(in) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_write_parallel_h5(dset_id, dcount, doffset, dset_dims, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_write_serial_h5(dset_id, dcount, doffset, dset_dims, filespace, memspace)
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine write_real_h5_1d

    !<--------------------------------------------------------------------------
    !< Write HDF5 dataset for 2D real data
    !<--------------------------------------------------------------------------
    subroutine write_real_h5_2d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(2), intent(in) :: dcount, doffset, dset_dims
        real(fp), dimension(:, :), intent(in) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_write_parallel_h5(dset_id, dcount, doffset, dset_dims, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_write_serial_h5(dset_id, dcount, doffset, dset_dims, filespace, memspace)
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine write_real_h5_2d

    !<--------------------------------------------------------------------------
    !< Write HDF5 dataset for 3D real data
    !<--------------------------------------------------------------------------
    subroutine write_real_h5_3d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(3), intent(in) :: dcount, doffset, dset_dims
        real(fp), dimension(:, :, :), intent(in) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_write_parallel_h5(dset_id, dcount, doffset, dset_dims, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_write_serial_h5(dset_id, dcount, doffset, dset_dims, filespace, memspace)
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine write_real_h5_3d

    !<--------------------------------------------------------------------------
    !< Write HDF5 dataset for 4D real data
    !<--------------------------------------------------------------------------
    subroutine write_real_h5_4d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(4), intent(in) :: dcount, doffset, dset_dims
        real(fp), dimension(:, :, :, :), intent(in) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_write_parallel_h5(dset_id, dcount, doffset, dset_dims, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_write_serial_h5(dset_id, dcount, doffset, dset_dims, filespace, memspace)
            call h5dwrite_f(dset_id, H5T_NATIVE_REAL, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine write_real_h5_4d

    !<--------------------------------------------------------------------------
    !< Write HDF5 dataset for 1D double data
    !<--------------------------------------------------------------------------
    subroutine write_double_h5_1d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        real(dp), dimension(:), intent(in) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_write_parallel_h5(dset_id, dcount, doffset, dset_dims, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_write_serial_h5(dset_id, dcount, doffset, dset_dims, filespace, memspace)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine write_double_h5_1d

    !<--------------------------------------------------------------------------
    !< Write HDF5 dataset for 2D double data
    !<--------------------------------------------------------------------------
    subroutine write_double_h5_2d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(2), intent(in) :: dcount, doffset, dset_dims
        real(dp), dimension(:, :), intent(in) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_write_parallel_h5(dset_id, dcount, doffset, dset_dims, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_write_serial_h5(dset_id, dcount, doffset, dset_dims, filespace, memspace)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine write_double_h5_2d

    !<--------------------------------------------------------------------------
    !< Write HDF5 dataset for 3D double data
    !<--------------------------------------------------------------------------
    subroutine write_double_h5_3d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(3), intent(in) :: dcount, doffset, dset_dims
        real(dp), dimension(:, :, :), intent(in) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_write_parallel_h5(dset_id, dcount, doffset, dset_dims, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_write_serial_h5(dset_id, dcount, doffset, dset_dims, filespace, memspace)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine write_double_h5_3d

    !<--------------------------------------------------------------------------
    !< Write HDF5 dataset for 4D double data
    !<--------------------------------------------------------------------------
    subroutine write_double_h5_4d(dset_id, dcount, doffset, dset_dims, fdata, &
            parallel_hdf5, collective_io)
        implicit none
        integer(hid_t), intent(in) :: dset_id
        integer(hsize_t), dimension(4), intent(in) :: dcount, doffset, dset_dims
        real(dp), dimension(:, :, :, :), intent(in) :: fdata
        logical, intent(in), optional :: parallel_hdf5, collective_io
        integer(hid_t) :: filespace, memspace, plist_id
        logical :: use_parallel_hdf5, use_collective_io
        integer :: error
        if (.not. present(parallel_hdf5)) then
            use_parallel_hdf5 = .false.
        else
            use_parallel_hdf5 = parallel_hdf5
        endif
        if (.not. present(collective_io)) then
            use_collective_io = .false.
        else
            use_collective_io = collective_io
        endif

        if (use_parallel_hdf5) then
            call init_write_parallel_h5(dset_id, dcount, doffset, dset_dims, &
                filespace, memspace, plist_id, use_collective_io)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace, xfer_prp=plist_id)
            call final_rw_h5(filespace, memspace, plist_id)
        else
            call init_write_serial_h5(dset_id, dcount, doffset, dset_dims, filespace, memspace)
            call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, fdata, dset_dims, error, &
                file_space_id=filespace, mem_space_id=memspace)
            call final_rw_h5(filespace, memspace)
        endif
    end subroutine write_double_h5_4d
end module hdf5_io

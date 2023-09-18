!*******************************************************************************
!< Module of diagnostics
!*******************************************************************************
module diagnostics
    use constants, only: fp, dp
    use simulation_setup_module, only: ndim_field
    use mhd_data_parallel, only: nfields, ngrads
    use particle_module, only: particle_type, ptls, escaped_ptls, &
        nptl_current, nptl_escaped, nptl_escaped_max, nptl_max, &
        spherical_coord_flag, leak, leak_negp, nptl_split, &
        get_interp_paramters, get_interp_paramters_spherical, &
        COUNT_FLAG_INBOX, COUNT_FLAG_OTHERS
    use mpi_module
    use hdf5
    implicit none
    private
    save
    public distributions_diagnostics, quick_check, &
        init_particle_distributions, clean_particle_distributions, &
        set_mpi_io_data_sizes, init_local_particle_distributions, &
        free_particle_distributions, free_local_particle_distributions, &
        get_pmax_global, dump_particles, &
        init_escaped_particles, free_escaped_particles, &
        reset_escaped_particles, dump_escaped_particles, &
        read_diagnostics_params

    !< Parameters for particle distributions
    real(dp) :: dx_diag, dy_diag, dz_diag
    real(dp) :: pmin_log, pmax_log, dp_log, dp_bands_log
    real(dp) :: pmax_local, pmax_global
    integer :: nx, ny, nz, npp, nreduce
    integer :: nx_mhd_reduced, ny_mhd_reduced, nz_mhd_reduced

    ! Particle distributions
    integer, parameter :: nbands = 12  ! Divide pmin to pmax by nbands
    real(dp), allocatable, dimension(:, :, :, :) :: fbands, fbands_sum
    real(dp), allocatable, dimension(:, :) :: fdpdt, fdpdt_sum
    real(dp), allocatable, dimension(:) :: fp_global, fp_global_sum
    real(dp), allocatable, dimension(:) :: parray, parray_bands
    ! Particle distribution in each cell
    real(dp), allocatable, dimension(:, :, :, :) :: fp_local, fp_local_sum

    real(dp) :: pmin  !< Minimum particle momentum
    real(dp) :: pmax  !< Maximum particle momentum

    !< MPI/IO data sizes
    integer, dimension(4) :: sizes_fxy, subsizes_fxy, starts_fxy
    integer, dimension(4) :: sizes_fp_local, subsizes_fp_local, starts_fp_local

    interface write_ptl_element
        module procedure &
            write_integer_element, write_double_element
    end interface write_ptl_element

    contains

    !---------------------------------------------------------------------------
    !< Quick diagnostics of the particle number information
    !< Args:
    !<  iframe: the time frame
    !<  if_create_file: whether to create a file
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine quick_check(iframe, if_create_file, file_path)
        implicit none
        integer, intent(in) :: iframe
        logical, intent(in) :: if_create_file
        character(*), intent(in) :: file_path
        integer, parameter :: nvar = 6
        real(dp) :: pdt_min, pdt_max, pdt_min_g, pdt_max_g
        real(dp), dimension(nvar) :: var_local, var_global
        integer :: i
        logical :: dir_e

        var_local = 0.0_dp
        var_global = 0.0_dp
        pdt_min = 1.0_dp
        pdt_max = 0.0_dp
        pdt_min_g = 0.0_dp
        pdt_max_g = 0.0_dp
        var_local(1) = nptl_current
        var_local(2) = nptl_split
        var_local(4) = leak
        var_local(5) = leak_negp
        do i = 1, nptl_current
            var_local(3) = var_local(3) + ptls(i)%weight
            if (ptls(i)%dt < pdt_min) pdt_min = ptls(i)%dt
            if (ptls(i)%dt > pdt_max) pdt_max = ptls(i)%dt
            var_local(6) = var_local(6) + ptls(i)%dt
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(var_local, var_global, nvar, MPI_DOUBLE_PRECISION, MPI_SUM, &
            master, MPI_COMM_WORLD, ierr)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(pdt_min, pdt_min_g, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
            master, MPI_COMM_WORLD, ierr)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(pdt_max, pdt_max_g, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
            master, MPI_COMM_WORLD, ierr)
        if (mpi_rank == master) then
            if (var_global(1) > 0) then
                var_global(6) = var_global(6) / var_global(1)
            else
                var_global(6) = 0.0_dp
            endif
            if (if_create_file) then
                open (17, file=trim(file_path)//'quick.dat', status='unknown')
                write(17, "(A6,8A13)") "iframe", "nptl_current", "nptl_split", &
                    "ntot", "leak", "leak_negp", "pdt_min", "pdt_max", "pdt_avg"
            else
                open (17, file=trim(file_path)//'quick.dat', status="old", &
                    position="append", action="write")
            endif
            write(17, "(I6.6,8E13.6)") iframe, var_global(1:5), &
                pdt_min_g, pdt_max_g, var_global(6)
            close(17)
        endif
    end subroutine quick_check

    !---------------------------------------------------------------------------
    !< Initialize the particle distributions
    !---------------------------------------------------------------------------
    subroutine init_particle_distributions
        use mhd_config_module, only: mhd_config
        implicit none
        integer :: i
        allocate(fbands(nx, ny, nz, nbands))
        allocate(fp_global(npp))
        allocate(fdpdt(npp, 2))
        if (mpi_cross_rank == master) then
            allocate(fbands_sum(nx, ny, nz, nbands))
            allocate(fp_global_sum(npp))
            allocate(fdpdt_sum(npp, 2))
        endif
        call clean_particle_distributions

        !< Intervals for distributions
        nx_mhd_reduced = (mhd_config%nx + nreduce - 1) / nreduce
        ny_mhd_reduced = (mhd_config%ny + nreduce - 1) / nreduce
        nz_mhd_reduced = (mhd_config%nz + nreduce - 1) / nreduce
        dx_diag = mhd_config%lx / nx_mhd_reduced
        dy_diag = mhd_config%ly / ny_mhd_reduced
        dz_diag = mhd_config%lz / nz_mhd_reduced
        pmin_log = log10(pmin)
        pmax_log = log10(pmax)
        dp_log = (pmax_log - pmin_log) / (npp - 1)

        !< Momentum bins for 1D distributions
        allocate(parray(npp))
        parray = 0.0_dp
        do i = 1, npp
            parray(i) = 10**(pmin_log + (i-1) * dp_log)
        enddo

        !< Momentum band edges for distributions in different energy band
        allocate(parray_bands(nbands+1))
        parray_bands = 0.0_dp
        dp_bands_log = (pmax_log - pmin_log) / nbands
        do i = 1, nbands+1
            parray_bands(i) = 10**(pmin_log + (i-1) * dp_bands_log)
        enddo

        if (mpi_rank == master) then
            write(*, "(A)") "Finished Initializing particle distributions."
        endif
    end subroutine init_particle_distributions

    !---------------------------------------------------------------------------
    !< Initialize local particle distributions
    !---------------------------------------------------------------------------
    subroutine init_local_particle_distributions
        implicit none
        allocate(fp_local(npp, nx, ny, nz))
        if (mpi_cross_rank == master) then
            allocate(fp_local_sum(npp, nx, ny, nz))
        endif
        call clean_local_particle_distribution
    end subroutine init_local_particle_distributions

    !---------------------------------------------------------------------------
    !< Set particle distributions to be zero
    !---------------------------------------------------------------------------
    subroutine clean_particle_distributions
        implicit none
        fbands = 0.0_dp
        fp_global = 0.0_dp
        fdpdt = 0.0_dp
        if (mpi_cross_rank == master) then
            fbands_sum = 0.0_dp
            fp_global_sum = 0.0_dp
            fdpdt_sum = 0.0_dp
        endif
    end subroutine clean_particle_distributions

    !---------------------------------------------------------------------------
    !< Set local particle distributions to be zero
    !---------------------------------------------------------------------------
    subroutine clean_local_particle_distribution
        implicit none
        fp_local = 0.0_dp
        if (mpi_cross_rank == master) then
            fp_local_sum = 0.0_dp
        endif
    end subroutine clean_local_particle_distribution

    !---------------------------------------------------------------------------
    !< Free particle distributions
    !---------------------------------------------------------------------------
    subroutine free_particle_distributions
        implicit none
        deallocate(fbands, parray_bands)
        deallocate(fp_global, parray)
        deallocate(fdpdt)
        if (mpi_cross_rank == master) then
            deallocate(fbands_sum)
            deallocate(fp_global_sum)
            deallocate(fdpdt_sum)
        endif
    end subroutine free_particle_distributions

    !---------------------------------------------------------------------------
    !< Free local particle distributions
    !---------------------------------------------------------------------------
    subroutine free_local_particle_distributions
        implicit none
        deallocate(fp_local)
        if (mpi_cross_rank == master) then
            deallocate(fp_local_sum)
        endif
    end subroutine free_local_particle_distributions

    !---------------------------------------------------------------------------
    !< Accumulate particle energization distributions
    !< Args:
    !<  iframe: time frame index (can be larger 1)
    !---------------------------------------------------------------------------
    subroutine energization_dist(iframe)
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config, tstamps_mhd, tstart_mhd
        use mhd_data_parallel, only: interp_fields
        implicit none
        integer, intent(in) :: iframe
        real(dp) :: weight, px, py, pz, rt
        real(dp) :: dvx_dx, dvy_dy, dvz_dz
        real(dp) :: t0, dtf
        type(particle_type) :: ptl
        integer, dimension(3) :: pos
        real(dp), dimension(8) :: weights
        real(dp), dimension(nfields+ngrads) :: fields !< Fields at particle position
        integer :: i, ip, mhd_tframe

        mhd_tframe = iframe - tstart_mhd
        t0 = tstamps_mhd(mhd_tframe)
        dtf = tstamps_mhd(mhd_tframe+1) - t0

        do i = 1, nptl_current
            ptl = ptls(i)
            if (ptl%count_flag == COUNT_FLAG_INBOX .and. &
                ptl%p > pmin .and. ptl%p <= pmax) then
                ip = ceiling((log10(ptl%p)-pmin_log) / dp_log)
                px = (ptl%x-fconfig%xmin) / mhd_config%dx
                py = (ptl%y-fconfig%ymin) / mhd_config%dy
                pz = (ptl%z-fconfig%zmin) / mhd_config%dz
                if (spherical_coord_flag) then
                    call get_interp_paramters_spherical(ptl%x, ptl%y, ptl%z, pos, weights)
                else
                    call get_interp_paramters(px, py, pz, pos, weights)
                endif
                rt = (ptl%t - t0) / dtf
                call interp_fields(pos, weights, rt, fields)
                dvx_dx = fields(nfields+1)
                dvy_dy = fields(nfields+5)
                dvz_dz = fields(nfields+9)

                fdpdt(ip, 1) = fdpdt(ip, 1) + ptl%weight
                fdpdt(ip, 2) = fdpdt(ip, 2) - &
                    ptl%p * (dvx_dx + dvy_dy + dvz_dz) / 3.0d0 * ptl%weight
            endif
        enddo
    end subroutine energization_dist

    !---------------------------------------------------------------------------
    !< Accumulate particle distributions
    !< Args:
    !<  local_dist: whether to accumulate local particle distribution
    !---------------------------------------------------------------------------
    subroutine calc_particle_distributions(local_dist)
        use simulation_setup_module, only: fconfig
        use mhd_config_module, only: mhd_config
        implicit none
        integer, intent(in) :: local_dist
        integer :: i, ix, iy, iz, ip
        real(dp) :: weight, p, xmin, xmax, ymin, ymax, zmin, zmax
        real(dp) :: px, py, pz, rx, ry, rz, rt
        real(dp) :: dxm, dym, dzm
        type(particle_type) :: ptl

        xmin = fconfig%xmin
        xmax = fconfig%xmax
        ymin = fconfig%ymin
        ymax = fconfig%ymax
        zmin = fconfig%zmin
        zmax = fconfig%zmax
        dxm = mhd_config%dx
        dym = mhd_config%dy
        dzm = mhd_config%dz

        do i = 1, nptl_current
            ptl = ptls(i)
            ix = floor((ptl%x - xmin)/dx_diag) + 1
            iy = floor((ptl%y - ymin)/dy_diag) + 1
            iz = floor((ptl%z - zmin)/dz_diag) + 1
            p = ptl%p
            weight = ptl%weight
            if (spherical_coord_flag) then
                ! For spherical coordinates, we solve for F=f*p^2*r^2 for 1D and
                ! f*p^2*r^2*sin(theta) for 2D and 3D
                weight = weight / (ptl%x)**2
                if (ndim_field >= 2) then
                    weight = weight / sin(ptl%y)
                endif
            endif

            !< Different momentum band
            if (ix >= 1 .and. ix <= nx .and. &
                iy >= 1 .and. iy <= ny .and. &
                iz >= 1 .and. iz <= nz) then
                ip = floor((log10(ptl%p)-pmin_log) / dp_bands_log)
                if (ip > 0 .and. ip < nbands) then
                    fbands(ix, iy, iz, ip+1) = fbands(ix, iy, iz, ip+1) + weight
                endif
            endif

            if (p > pmin .and. p <= pmax) then
                ip = ceiling((log10(ptl%p)-pmin_log) / dp_log)
                fp_global(ip) = fp_global(ip) + weight
                if (local_dist == 1) then
                    if (ix >= 1 .and. ix <= nx .and. &
                        iy >= 1 .and. iy <= ny .and. &
                        iz >= 1 .and. iz <= nz) then
                        fp_local(ip, ix, iy, iz) = fp_local(ip, ix, iy, iz) + weight
                    endif
                endif
            endif
        enddo

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fp_global, fp_global_sum, npp, MPI_DOUBLE_PRECISION, &
            MPI_SUM, master, MPI_COMM_WORLD, ierr)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fdpdt, fdpdt_sum, npp*2, MPI_DOUBLE_PRECISION, &
            MPI_SUM, master, MPI_COMM_WORLD, ierr)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(fbands, fbands_sum, nx*ny*nz*nbands, MPI_DOUBLE_PRECISION, &
            MPI_SUM, master, mpi_cross_comm, ierr)
        if (local_dist == 1) then
            call MPI_BARRIER(MPI_COMM_WORLD, ierr)
            call MPI_REDUCE(fp_local, fp_local_sum, npp*nx*ny*nz, &
                MPI_DOUBLE_PRECISION, MPI_SUM, master, mpi_cross_comm, ierr)
        endif
    end subroutine calc_particle_distributions

    !---------------------------------------------------------------------------
    !< Diagnostics of the particle distributions
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !<  local_dist: whether to accumulate local particle distribution
    !---------------------------------------------------------------------------
    subroutine distributions_diagnostics(iframe, file_path, local_dist)
        use mpi_io_module, only: set_mpi_datatype_double, set_mpi_info, fileinfo, &
            open_data_mpi_io, write_data_mpi_io
        use mhd_config_module, only: mhd_config
        use constants, only: fp, dp
        implicit none
        integer, intent(in) :: iframe, local_dist
        character(*), intent(in) :: file_path
        integer :: fh, pos1
        character(len=4) :: ctime
        character(len=128) :: fname
        integer(kind=MPI_OFFSET_KIND) :: disp, offset
        integer :: mpi_datatype
        logical :: dir_e

        call calc_particle_distributions(local_dist)
        call energization_dist(iframe)

        write (ctime,'(i4.4)') iframe
        if (mpi_rank .eq. master) then
            fh = 15
            fname = trim(file_path)//'fp-'//ctime//'_sum.dat'
            open(fh, file=trim(fname), access='stream', status='unknown', &
                 form='unformatted', action='write')
            write(fh, pos=1) parray
            pos1 = npp * sizeof(1.0_dp) + 1
            write(fh, pos=pos1) fp_global_sum
            close(fh)

            fname = trim(file_path)//'fdpdt-'//ctime//'_sum.dat'
            open(fh, file=trim(fname), access='stream', status='unknown', &
                 form='unformatted', action='write')
            write(fh, pos=1) parray
            pos1 = npp * sizeof(1.0_dp) + 1
            write(fh, pos=pos1) fdpdt_sum
            close(fh)
        endif

        if (mpi_cross_rank == master) then
            call set_mpi_info
            ! fxy for different energy band
            fh = 18
            fname = trim(file_path)//'fxy-'//ctime//'_sum.dat'
            if (mpi_sub_rank == master) then
                fh = 18
                open(fh, file=trim(fname), access='stream', status='unknown', &
                     form='unformatted', action='write')
                write(fh, pos=1) (nbands + 0.0_dp)
                pos1 = sizeof(1.0_dp) + 1
                write(fh, pos=pos1) parray_bands
                close(fh)
            endif
            disp = (nbands + 2) * sizeof(1.0_dp)
            offset = 0
            mpi_datatype = set_mpi_datatype_double(sizes_fxy, subsizes_fxy, starts_fxy)
            call open_data_mpi_io(trim(fname), MPI_MODE_APPEND+MPI_MODE_WRONLY, &
                fileinfo, mpi_sub_comm, fh)
            call write_data_mpi_io(fh, mpi_datatype, subsizes_fxy, disp, offset, fbands_sum)
            call MPI_FILE_CLOSE(fh, ierror)

            if (local_dist == 1) then
                fh = 19
                fname = trim(file_path)//'fp_local_'//ctime//'_sum.dat'
                disp = 0
                offset = 0
                mpi_datatype = set_mpi_datatype_double(sizes_fp_local, subsizes_fp_local, starts_fp_local)
                call open_data_mpi_io(trim(fname), MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                    fileinfo, mpi_sub_comm, fh)
                call write_data_mpi_io(fh, mpi_datatype, subsizes_fp_local, disp, offset, fp_local_sum)
                call MPI_FILE_CLOSE(fh, ierror)
            endif
        endif

        call clean_particle_distributions
        if (local_dist == 1) then
            call clean_local_particle_distribution
        endif
    end subroutine distributions_diagnostics


    !---------------------------------------------------------------------------
    !< Get maximum particle momentum and write to file
    !< Args:
    !<  iframe: the time frame
    !<  if_create_file: whether to create a file
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine get_pmax_global(iframe, if_create_file, file_path)
        implicit none
        integer, intent(in) :: iframe
        logical, intent(in) :: if_create_file
        character(*), intent(in) :: file_path
        integer :: i
        logical :: dir_e

        pmax_local = 0.0_dp
        do i = 1, nptl_current
            if (ptls(i)%p > pmax_local) then
                pmax_local = ptls(i)%p
            endif
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_REDUCE(pmax_local, pmax_global, 1, MPI_DOUBLE_PRECISION, &
            MPI_MAX, master, MPI_COMM_WORLD, ierr)
        if (mpi_rank == master) then
            if (if_create_file) then
                open (17, file=trim(file_path)//'pmax_global.dat', status='unknown')
            else
                open (17, file=trim(file_path)//'pmax_global.dat', status="old", &
                    position="append", action="write")
            endif
            write(17, "(E13.6)") pmax_global
            close(17)
        endif
    end subroutine get_pmax_global

    !---------------------------------------------------------------------------
    !< Set MPI/IO data sizes for distribution diagnostics
    !---------------------------------------------------------------------------
    subroutine set_mpi_io_data_sizes
        use mhd_config_module, only: mhd_config
        use simulation_setup_module, only: mpi_ix, mpi_iy, mpi_iz
        implicit none
        sizes_fxy = (/ nx_mhd_reduced, ny_mhd_reduced, nz_mhd_reduced, nbands /)
        subsizes_fxy = (/ nx, ny, nz, nbands /)
        starts_fxy = (/ nx * mpi_ix, ny * mpi_iy, nz * mpi_iz, 0 /)

        sizes_fp_local = (/ npp, nx_mhd_reduced, ny_mhd_reduced, nz_mhd_reduced /)
        subsizes_fp_local = (/ npp, nx, ny, nz /)
        starts_fp_local = (/ 0, nx * mpi_ix, ny * mpi_iy, nz * mpi_iz /)
    end subroutine set_mpi_io_data_sizes

    !---------------------------------------------------------------------------
    !< Write one double element of the particle data
    !---------------------------------------------------------------------------
    subroutine write_double_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: write_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        real(dp), dimension(:), intent(in) :: fdata
        integer(hid_t) :: dset_id, filespace
        integer :: error
        call h5screate_simple_f(1, dset_dims, filespace, error)
        call h5dcreate_f(file_id, trim(dset_name), H5T_NATIVE_DOUBLE, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .true., .true.)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)
    end subroutine write_double_element

    !---------------------------------------------------------------------------
    !< Write one integer element of the particle data
    !---------------------------------------------------------------------------
    subroutine write_integer_element(file_id, dcount, doffset, dset_dims, &
            dset_name, fdata)
        use hdf5_io, only: write_data_h5
        implicit none
        integer(hid_t), intent(in) :: file_id
        integer(hsize_t), dimension(1), intent(in) :: dcount, doffset, dset_dims
        character(*), intent(in) :: dset_name
        integer, dimension(:), intent(in) :: fdata
        integer(hid_t) :: dset_id, filespace
        integer :: error
        call h5screate_simple_f(1, dset_dims, filespace, error)
        call h5dcreate_f(file_id, trim(dset_name), H5T_NATIVE_INTEGER, &
            filespace, dset_id, error)
        call write_data_h5(dset_id, dcount, doffset, dset_dims, fdata, .true., .true.)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(filespace, error)
    end subroutine write_integer_element

    !---------------------------------------------------------------------------
    !< dump all the particles
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine dump_particles(iframe, file_path)
        use hdf5_io, only: create_file_h5, close_file_h5
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        integer(hsize_t), dimension(1) :: dcount, doffset, dset_dims
        real(dp) :: nptl_local, nptl_global, nptl_offset
        character(len=128) :: fname
        character(len=4) :: ctime
        integer(hid_t) :: file_id
        integer :: error
        logical :: dir_e

        CALL h5open_f(error)

        write (ctime,'(i4.4)') iframe
        fname = trim(file_path)//'particles_'//ctime//'.h5'
        call create_file_h5(fname, H5F_ACC_TRUNC_F, file_id, .true., MPI_COMM_WORLD)

        nptl_local = nptl_current
        call MPI_ALLREDUCE(nptl_local, nptl_global, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_SCAN(nptl_local, nptl_offset, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
        nptl_offset = nptl_offset - nptl_local

        dcount(1) = nptl_local
        doffset(1) = nptl_offset
        dset_dims(1) = nptl_global

        call write_ptl_element(file_id, dcount, doffset, dset_dims, "x", ptls%x)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "y", ptls%y)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "z", ptls%z)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "p", ptls%p)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "weight", ptls%weight)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "t", ptls%t)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, "dt", ptls%dt)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, &
            "split_times", ptls%split_times)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, &
            "count_flag", ptls%count_flag)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, &
            "tag", ptls%tag)
        call write_ptl_element(file_id, dcount, doffset, dset_dims, &
            "nsteps_tracking", ptls%nsteps_tracking)

        call close_file_h5(file_id)
        CALL h5close_f(error)
    end subroutine dump_particles

    !---------------------------------------------------------------------------
    !< dump escaped particles
    !< Args:
    !<  iframe: time frame index
    !<  file_path: save data files to this path
    !---------------------------------------------------------------------------
    subroutine dump_escaped_particles(iframe, file_path)
        use hdf5_io, only: create_file_h5, close_file_h5
        implicit none
        integer, intent(in) :: iframe
        character(*), intent(in) :: file_path
        integer(hsize_t), dimension(1) :: dcount, doffset, dset_dims
        real(dp) :: nptl_local, nptl_global, nptl_offset
        character(len=128) :: fname
        character(len=4) :: ctime
        integer(hid_t) :: file_id
        integer :: error
        logical :: dir_e

        nptl_local = nptl_escaped
        call MPI_ALLREDUCE(nptl_local, nptl_global, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
        call MPI_SCAN(nptl_local, nptl_offset, 1, MPI_DOUBLE_PRECISION, &
            MPI_SUM, MPI_COMM_WORLD, ierr)
        nptl_offset = nptl_offset - nptl_local

        if (nptl_global > 0) then
            CALL h5open_f(error)

            write (ctime,'(i4.4)') iframe
            fname = trim(file_path)//'escaped_particles_'//ctime//'.h5'
            call create_file_h5(fname, H5F_ACC_TRUNC_F, file_id, .true., MPI_COMM_WORLD)

            dcount(1) = nptl_local
            doffset(1) = nptl_offset
            dset_dims(1) = nptl_global

            call write_ptl_element(file_id, dcount, doffset, dset_dims, "x", escaped_ptls%x)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "y", escaped_ptls%y)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "z", escaped_ptls%z)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "p", escaped_ptls%p)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "weight", escaped_ptls%weight)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "t", escaped_ptls%t)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, "dt", escaped_ptls%dt)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, &
                "split_times", escaped_ptls%split_times)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, &
                "count_flag", escaped_ptls%count_flag)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, &
                "tag", escaped_ptls%tag)
            call write_ptl_element(file_id, dcount, doffset, dset_dims, &
                "nsteps_tracking", escaped_ptls%nsteps_tracking)

            call close_file_h5(file_id)
            CALL h5close_f(error)
        else
            if (mpi_rank == master) then
                write(*, "(A)") "No escaped particles during this time interval"
            endif
        endif
    end subroutine dump_escaped_particles

    !---------------------------------------------------------------------------
    !< Read parameters for diagnostics
    !< Args:
    !<  conf_file: configuration file name
    !---------------------------------------------------------------------------
    subroutine read_diagnostics_params(conf_file)
        use read_config, only: get_variable
        use simulation_setup_module, only: fconfig
        implicit none
        character(*), intent(in) :: conf_file
        logical :: condx, condy, condz
        real(fp) :: temp
        integer :: fh

        if (mpi_rank == master) then
            fh = 10
            open(unit=fh, file='config/'//trim(conf_file), status='old')
            pmin = get_variable(fh, 'pmin', '=')
            pmax = get_variable(fh, 'pmax', '=')
            temp = get_variable(fh, 'nreduce', '=')
            nreduce = int(temp)
            nx = (fconfig%nx + nreduce - 1) / nreduce
            ny = (fconfig%ny + nreduce - 1) / nreduce
            nz = (fconfig%nz + nreduce - 1) / nreduce
            temp = get_variable(fh, 'npp', '=')
            npp = int(temp)
            close(fh)
            !< echo the information
            print *, "---------------------------------------------------"
            write(*, "(A,I0,A,I0,A,I0)") " Dimensions of spatial distributions = ", &
                nx, " ", ny, " ", nz
            write(*, "(A,I0)") " Dimensions of momentum distributions = ", npp
            print *, "---------------------------------------------------"
        endif
        call MPI_BCAST(nreduce, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nx, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ny, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nz, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(npp, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr)

        if (ndim_field == 1) then
            if (nx * nreduce /= fconfig%nx) then
                if (mpi_rank == master) then
                    write(*, "(A)") "Wrong factor 'nreduce' for particle distribution"
                endif
                call MPI_FINALIZE(ierr)
                stop
            endif
        else if (ndim_field == 2) then
            if (nx * nreduce /= fconfig%nx .or. &
                ny * nreduce /= fconfig%ny) then
                if (mpi_rank == master) then
                    write(*, "(A)") "Wrong factor 'nreduce' for particle distribution"
                endif
                call MPI_FINALIZE(ierr)
                stop
            endif
        else
            if (nx * nreduce /= fconfig%nx .or. &
                ny * nreduce /= fconfig%ny .or. &
                nz * nreduce /= fconfig%nz) then
                if (mpi_rank == master) then
                    write(*, "(A)") "Wrong factor 'nreduce' for particle distribution"
                endif
                call MPI_FINALIZE(ierr)
                stop
            endif
        endif
    end subroutine read_diagnostics_params

    !---------------------------------------------------------------------------
    !< Reset escaped particle data array
    !---------------------------------------------------------------------------
    subroutine reset_escaped_particles
        implicit none
        escaped_ptls%x = 0.0
        escaped_ptls%y = 0.0
        escaped_ptls%z = 0.0
        escaped_ptls%p = 0.0
        escaped_ptls%weight = 0.0
        escaped_ptls%t = 0.0
        escaped_ptls%dt = 0.0
        escaped_ptls%split_times = 0
        escaped_ptls%count_flag = COUNT_FLAG_OTHERS
        escaped_ptls%tag = 0
        escaped_ptls%nsteps_tracking = 0
        nptl_escaped = 0
    end subroutine reset_escaped_particles

    !---------------------------------------------------------------------------
    !< Initialize escaped particle data
    !---------------------------------------------------------------------------
    subroutine init_escaped_particles
        use simulation_setup_module, only: particle_can_escape
        implicit none
        if (particle_can_escape) then
            nptl_escaped_max = nptl_max
        else
            nptl_escaped_max = 1
        endif
        allocate(escaped_ptls(nptl_escaped_max))
        call reset_escaped_particles
    end subroutine init_escaped_particles

    !---------------------------------------------------------------------------
    !< Free escaped particle data
    !---------------------------------------------------------------------------
    subroutine free_escaped_particles
        implicit none
        deallocate(escaped_ptls)
    end subroutine free_escaped_particles

end module diagnostics
